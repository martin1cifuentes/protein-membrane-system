using System.Text.Json;
using System.Text.Json.Serialization;
using System.Threading.Channels;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;

namespace ProteinInMembrane.Host;

// The browser presents the root's account. It never establishes scientific
// support, completion, or exportability from a successful HTTP exchange.
public static class BrowserExchange
{
    private const long MaximumCommandBytes = 1_000_000;
    private const long MaximumCoordinateBytes = 100_000_000;
    private static readonly JsonSerializerOptions WireJson = CreateWireJson();
    private static readonly HashSet<string> CoordinateExtensions =
        new(StringComparer.OrdinalIgnoreCase) { ".pdb", ".cif", ".mmcif" };

    public static void Map(WebApplication app, ProductRoot system, Uri origin)
    {
        app.MapGet("/api/state", () => Results.Json(system.Snapshot(), WireJson));

        app.MapPost("/api/commands", async (HttpContext context, CancellationToken cancellationToken) =>
        {
            if (!HasSameOrigin(context.Request, origin)) return Refused("The command must come from this local workspace.", 403);
            ActorCommand? command;
            string? suppliedKind;
            try
            {
                var bytes = await ReadBoundedAsync(context.Request.Body, MaximumCommandBytes, cancellationToken);
                using var document = JsonDocument.Parse(bytes);
                if (document.RootElement.ValueKind != JsonValueKind.Object ||
                    !document.RootElement.TryGetProperty("kind", out var kindValue) ||
                    kindValue.ValueKind != JsonValueKind.String)
                    return Refused("A recognized actor action is required.", 400);
                suppliedKind = kindValue.GetString();
                command = JsonSerializer.Deserialize<ActorCommand>(bytes, WireJson);
            }
            catch (JsonException)
            {
                return Refused("The command is not valid JSON.", 400);
            }
            catch (RequestTooLargeException)
            {
                return Refused("The command exceeds the local request limit.", 413);
            }

            if (command is null || !Enum.IsDefined(command.Kind) ||
                command.Kind == ActorActionKind.DeclinePreparationChange ||
                JsonSerializer.Serialize(command.Kind, WireJson) != JsonSerializer.Serialize(suppliedKind) ||
                command.Data.ValueKind != JsonValueKind.Object)
                return Refused("A recognized actor action, object data and account revision are required.", 400);

            var outcome = await system.ExecuteAsync(command, cancellationToken);
            return outcome.Value is { } account
                ? Results.Json(account, WireJson)
                : Refused(outcome.Reason, 422);
        });

        app.MapPost("/api/uploads", async (HttpContext context, CancellationToken cancellationToken) =>
        {
            if (!HasSameOrigin(context.Request, origin)) return Refused("The upload must come from this local workspace.", 403);
            if (!context.Request.HasFormContentType) return Refused("Supply one coordinate file as multipart form data.", 400);

            IFormCollection form;
            try
            {
                form = await context.Request.ReadFormAsync(cancellationToken);
            }
            catch (InvalidDataException)
            {
                return Refused("The uploaded coordinate file could not be read.", 400);
            }
            catch (BadHttpRequestException)
            {
                return Refused("The uploaded coordinate file exceeds the local request limit.", 413);
            }

            if (form.Files.Count != 1 || form.Files[0].Name != "file")
                return Refused("Supply exactly one coordinate file.", 400);
            if (!form.TryGetValue("provenance", out var provenanceValues) || provenanceValues.Count != 1 ||
                ParseUploadOrigin(provenanceValues[0]) is null ||
                form.TryGetValue("provenanceNote", out var noteValues) && noteValues.Count > 1)
                return Refused("Declare the uploaded structure's provenance as predicted, experimental, or unknown.", 400);
            var file = form.Files[0];
            if (file.Length is <= 0 or > MaximumCoordinateBytes ||
                !CoordinateExtensions.Contains(Path.GetExtension(file.FileName)))
                return Refused("Supply a nonempty PDB or mmCIF file within the 100 MB limit.", 413);

            try
            {
                await using var content = file.OpenReadStream();
                var token = await system.UploadAsync(content, file.FileName, ParseUploadOrigin(provenanceValues[0])!.Value,
                    form.TryGetValue("provenanceNote", out var note) && note.Count == 1 ? note[0] : null,
                    cancellationToken);
                return Results.Json(new { uploadToken = token }, WireJson);
            }
            catch (InvalidDataException exception)
            {
                return Refused(exception.Message, 400);
            }
        });

        app.MapGet("/api/structures/{artifactId}", (string artifactId) =>
        {
            if (!IsOpaqueSegment(artifactId)) return Results.NotFound();
            var path = system.StructurePath(artifactId);
            return path is not null && File.Exists(path)
                ? Results.File(path, Path.GetExtension(path).ToLowerInvariant() switch
                {
                    ".cif" or ".mmcif" => "chemical/x-mmcif",
                    ".pdb" or ".ent" => "chemical/x-pdb",
                    _ => "application/octet-stream"
                }, enableRangeProcessing: true)
                : Results.NotFound();
        });

        app.MapGet("/api/export/{stageId}", (string stageId) =>
        {
            if (!IsOpaqueSegment(stageId)) return Results.NotFound();
            var path = system.ExportPath(stageId);
            return path is not null && File.Exists(path)
                ? Results.File(path, "application/zip", $"protein-membrane-{stageId}.zip", enableRangeProcessing: true)
                : Results.NotFound();
        });

        app.MapGet("/api/events", async (HttpContext context) =>
        {
            context.Response.Headers.CacheControl = "no-store";
            context.Response.Headers["X-Accel-Buffering"] = "no";
            context.Response.ContentType = "text/event-stream";
            var updates = Channel.CreateBounded<byte>(new BoundedChannelOptions(1)
            {
                SingleReader = true,
                SingleWriter = false,
                FullMode = BoundedChannelFullMode.DropOldest
            });
            void Changed() => updates.Writer.TryWrite(1);
            system.Changed += Changed;
            try
            {
                // The first notice closes the gap between the initial account
                // read and subscription; every later notice asks for a fresh account.
                await context.Response.WriteAsync("data: refresh\n\n", context.RequestAborted);
                await context.Response.Body.FlushAsync(context.RequestAborted);
                await foreach (var _ in updates.Reader.ReadAllAsync(context.RequestAborted))
                {
                    await context.Response.WriteAsync("data: refresh\n\n", context.RequestAborted);
                    await context.Response.Body.FlushAsync(context.RequestAborted);
                }
            }
            catch (OperationCanceledException) when (context.RequestAborted.IsCancellationRequested)
            {
                // Disconnecting a view does not stop the application-owned run.
            }
            finally
            {
                system.Changed -= Changed;
                updates.Writer.TryComplete();
            }
        });
    }

    private static bool HasSameOrigin(HttpRequest request, Uri origin)
    {
        if (!Uri.TryCreate(request.Headers["Origin"].ToString(), UriKind.Absolute, out var submitted)) return false;
        return submitted.Scheme == origin.Scheme &&
               submitted.Host.Equals(origin.Host, StringComparison.OrdinalIgnoreCase) &&
               submitted.Port == origin.Port &&
               request.Host.Host.Equals(origin.Host, StringComparison.OrdinalIgnoreCase) &&
               request.Host.Port == origin.Port;
    }

    private static bool IsOpaqueSegment(string value) =>
        value.Length is > 0 and <= 160 &&
        char.IsAsciiLetterOrDigit(value[0]) &&
        value.All(character => char.IsAsciiLetterOrDigit(character) || character is '-' or '_' or '.');

    private static IResult Refused(string reason, int status) =>
        Results.Json(new { reason }, WireJson, statusCode: status);

    private static JsonSerializerOptions CreateWireJson()
    {
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        return options;
    }

    private static UploadOriginKind? ParseUploadOrigin(string? value) => value switch
    {
        "predicted" => UploadOriginKind.Predicted,
        "experimental" => UploadOriginKind.Experimental,
        "unknown" => UploadOriginKind.Unknown,
        _ => null
    };

    private static async Task<byte[]> ReadBoundedAsync(Stream input, long maximum, CancellationToken cancellationToken)
    {
        using var buffer = new MemoryStream();
        var chunk = new byte[16_384];
        int count;
        while ((count = await input.ReadAsync(chunk, cancellationToken)) != 0)
        {
            if (buffer.Length + count > maximum) throw new RequestTooLargeException();
            await buffer.WriteAsync(chunk.AsMemory(0, count), cancellationToken);
        }
        return buffer.ToArray();
    }

    private sealed class RequestTooLargeException : Exception { }
}
