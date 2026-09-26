using System.Net;
using System.Diagnostics;
using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;
using Microsoft.Extensions.FileProviders;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;

var originText = Environment.GetEnvironmentVariable("ASPNETCORE_URLS") ?? "http://127.0.0.1:4185";
if (!Uri.TryCreate(originText, UriKind.Absolute, out var origin) ||
    origin.Scheme != Uri.UriSchemeHttp || origin.Host != IPAddress.Loopback.ToString() ||
    origin.Port is < 1024 or > 65535 || origin.AbsolutePath != "/" ||
    origin.Query.Length != 0 || origin.Fragment.Length != 0)
    throw new InvalidOperationException("The local host must use one HTTP address on 127.0.0.1 and an unprivileged port.");

var python = RequiredEnvironmentPath("PIM_WORKER_PYTHON");
var ownerSourceRoot = RequiredEnvironmentPath("PIM_OWNER_SOURCE_ROOT");
var workspaceRoot = RequiredEnvironmentPath("PIM_WORKSPACE_ROOT");
var browserRoot = Path.Combine(AppContext.BaseDirectory, "wwwroot");
if (!File.Exists(Path.Combine(browserRoot, "index.html")))
    throw new InvalidOperationException("The local browser bundle is absent. Build the product before starting it.");

var builder = WebApplication.CreateBuilder(new WebApplicationOptions
{
    Args = args,
    WebRootPath = browserRoot
});
builder.WebHost.ConfigureKestrel(options =>
{
    options.Listen(IPAddress.Loopback, origin.Port);
    options.Limits.MaxRequestBodySize = 101_000_000;
});

// Source URLs are selected by this product. Do not follow a provider redirect
// to an unreviewed host or local address as part of structure retrieval.
using var remoteClient = new HttpClient(new HttpClientHandler { AllowAutoRedirect = false });
var sources = new ExternalSourceExchange(remoteClient);
var worker = new ScientificWorkerExchange(python, ownerSourceRoot);
var system = new ProductRoot(worker, sources, workspaceRoot,
    Environment.GetEnvironmentVariable("PIM_POLICY_CATALOGUE") ?? string.Empty,
    Environment.GetEnvironmentVariable("PIM_PPM_EXECUTABLE") ?? string.Empty,
    () => ProbeConstructionProvider(python));

var app = builder.Build();
app.UseExceptionHandler(error => error.Run(async context =>
{
    context.Response.StatusCode = StatusCodes.Status500InternalServerError;
    await context.Response.WriteAsJsonAsync(new { reason = "The local product encountered an unexpected error." });
}));
app.Use(async (context, next) =>
{
    if (!context.Request.Host.Host.Equals(IPAddress.Loopback.ToString(), StringComparison.OrdinalIgnoreCase) ||
        context.Request.Host.Port != origin.Port ||
        !IPAddress.IsLoopback(context.Connection.RemoteIpAddress ?? IPAddress.None))
    {
        context.Response.StatusCode = StatusCodes.Status403Forbidden;
        return;
    }
    context.Response.Headers["X-Content-Type-Options"] = "nosniff";
    context.Response.Headers["X-Frame-Options"] = "DENY";
    context.Response.Headers["Referrer-Policy"] = "no-referrer";
    context.Response.Headers.CacheControl = "no-store";
    await next();
});

BrowserExchange.Map(app, system, origin);
var browserFiles = new PhysicalFileProvider(browserRoot);
app.UseDefaultFiles(new DefaultFilesOptions { FileProvider = browserFiles });
app.UseStaticFiles(new StaticFileOptions { FileProvider = browserFiles });
app.MapFallback((HttpContext context) =>
    context.Request.Path.StartsWithSegments("/api", StringComparison.OrdinalIgnoreCase)
        ? Results.NotFound()
        : Results.File(Path.Combine(browserRoot, "index.html"), "text/html"));

await app.RunAsync();

static string RequiredEnvironmentPath(string name)
{
    var path = Environment.GetEnvironmentVariable(name);
    if (string.IsNullOrWhiteSpace(path))
        throw new InvalidOperationException($"The local prerequisite {name} is not configured.");
    return Path.GetFullPath(path);
}

static ConstructionProviderInstallation? ProbeConstructionProvider(string python)
{
    const string script = "import json, importlib.metadata, openmm, openmm.app; " +
        "from pathlib import Path; " +
        "print(json.dumps({'version': openmm.version.full_version, " +
        "'packageVersion': importlib.metadata.version('openmm'), " +
        "'dataPath': str((Path(openmm.app.__file__).resolve().parent / 'data').resolve())}))";
    using var process = new Process();
    process.StartInfo = new ProcessStartInfo
    {
        FileName = python,
        UseShellExecute = false,
        RedirectStandardOutput = true,
        RedirectStandardError = true,
        CreateNoWindow = true
    };
    process.StartInfo.ArgumentList.Add("-c");
    process.StartInfo.ArgumentList.Add(script);
    using var deadline = new CancellationTokenSource(TimeSpan.FromSeconds(10));
    try
    {
        process.Start();
        var output = process.StandardOutput.ReadToEndAsync(deadline.Token);
        var errors = process.StandardError.ReadToEndAsync(deadline.Token);
        process.WaitForExitAsync(deadline.Token).GetAwaiter().GetResult();
        var response = output.GetAwaiter().GetResult();
        _ = errors.GetAwaiter().GetResult();
        if (process.ExitCode != 0) return null;
        using var document = JsonDocument.Parse(response);
        var root = document.RootElement;
        if (root.GetProperty("packageVersion").GetString() != "8.6.0") return null;
        var version = root.GetProperty("version").GetString();
        var dataPath = root.GetProperty("dataPath").GetString();
        if (string.IsNullOrWhiteSpace(version) || string.IsNullOrWhiteSpace(dataPath) ||
            !Path.IsPathFullyQualified(dataPath)) return null;
        string PatchPath(string species) => Path.Combine(dataPath, species + ".pdb");
        string? PatchHash(string path)
        {
            if (!File.Exists(path)) return null;
            using var patch = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(patch)).ToLowerInvariant();
        }
        var dmpcPath = PatchPath("DMPC");
        var popcPath = PatchPath("POPC");
        var popcHash = PatchHash(popcPath);
        return new ConstructionProviderInstallation(version, dmpcPath,
            PatchHash(dmpcPath) ?? string.Empty,
            popcHash is null ? ImmutableArray<NativePatchInstallation>.Empty :
                ImmutableArray.Create(new NativePatchInstallation("POPC", popcPath, popcHash)));
    }
    catch (Exception exception) when (exception is IOException or UnauthorizedAccessException or
        System.ComponentModel.Win32Exception or InvalidOperationException or JsonException or
        KeyNotFoundException or OperationCanceledException)
    {
        try { if (!process.HasExited) process.Kill(entireProcessTree: true); }
        catch (InvalidOperationException) { }
        return null;
    }
}
