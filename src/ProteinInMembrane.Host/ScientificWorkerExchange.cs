using System.Collections.Immutable;
using System.Diagnostics;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Text.Json.Serialization;
using ProteinInMembrane.Host.ProteinInMembraneSystem;

namespace ProteinInMembrane.Host;

// This adapter transports one bounded request to one local worker process. It
// verifies the physical result; scientific conclusions remain with C# children.
public sealed class ScientificWorkerExchange : IScientificWorkerExchange
{
    private static readonly JsonSerializerOptions WireJson = new(JsonSerializerDefaults.Web)
    {
        Converters = { new JsonStringEnumConverter(JsonNamingPolicy.CamelCase) }
    };

    private readonly string _pythonExecutable;
    private readonly string _ownerSourceRoot;

    public ScientificWorkerExchange(string pythonExecutable, string ownerSourceRoot)
    {
        _pythonExecutable = Path.GetFullPath(pythonExecutable);
        _ownerSourceRoot = Path.GetFullPath(ownerSourceRoot);
    }

    public event Action<string, JsonElement>? ProgressObserved;

    public Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
        ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<SourceInspectionPayload, SourceInspectionObservations>("inspect_source", request, cancellationToken);

    public Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
        ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<PredictionRegionSummaryPayload, PredictionRegionSummaryObservations>(
            "summarize_prediction_evidence", request, cancellationToken);

    public Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
        ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<PreparationChangeInspectionPayload, PreparationChangeObservations>("inspect_preparation_changes", request, cancellationToken);

    public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
        ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<MembraneAssessmentPayload, MembraneAssessmentObservations>("assess_membrane", request, cancellationToken);

    public Task<WorkerResult<ConstructionInputObservations>> MeasureConstructionInputsAsync(
        ScientificWorkRequest<ConstructionInputPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<ConstructionInputPayload, ConstructionInputObservations>("measure_construction_inputs", request, cancellationToken);

    public Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
        ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<ProteinPreparationPayload, ProteinPreparationObservations>("prepare_protein", request, cancellationToken);

    public Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
        ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<PlacementPayload, PlacementObservations>("place_ppm", request, cancellationToken);

    public Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
        ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<PlacementAdjustmentPayload, PlacementAdjustmentObservations>("adjust_placement", request, cancellationToken);

    public Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
        ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<PlacementMeasurementPayload, PlacementMeasurementObservations>("measure_placement", request, cancellationToken);

    public Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
        ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<ConstructionPayload, ConstructionObservations>("construct_system", request, cancellationToken);

    public Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
        ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<MinimizationPayload, MinimizationObservations>("minimize", request, cancellationToken);

    public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
        ScientificWorkRequest<EquilibrationPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<EquilibrationPayload, EquilibrationObservations>("equilibrate", request, cancellationToken);

    public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<StageObservationPayload, StageObservationObservations>("observe_stage", request, cancellationToken);

    public Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
        ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken)
        => InvokeAsync<ExportVerificationPayload, ExportVerificationObservations>("verify_export", request, cancellationToken);

    private async Task<WorkerResult<TObservation>> InvokeAsync<TPayload, TObservation>(
        string operation,
        ScientificWorkRequest<TPayload> request,
        CancellationToken cancellationToken)
        where TPayload : class
        where TObservation : class
    {
        if (!File.Exists(_pythonExecutable) || !Directory.Exists(_ownerSourceRoot))
            return Unobserved<TObservation>(request.RequestId, "worker-unavailable", "The selected local Python worker is not installed.");

        var workingDirectory = Path.GetFullPath(request.WorkingDirectory);
        if (!Directory.Exists(workingDirectory))
            return Unobserved<TObservation>(request.RequestId, "working-directory-unavailable", "The identified attempt directory is unavailable.");

        var start = new ProcessStartInfo(_pythonExecutable)
        {
            WorkingDirectory = _ownerSourceRoot,
            RedirectStandardInput = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            UseShellExecute = false,
            CreateNoWindow = true
        };
        start.ArgumentList.Add("-m");
        start.ArgumentList.Add("ProteinInMembraneSystem.worker");
        start.Environment["PYTHONPATH"] = _ownerSourceRoot;

        using var process = new Process { StartInfo = start };
        var started = false;
        try
        {
            if (!process.Start())
                return Unobserved<TObservation>(request.RequestId, "worker-not-started", "The local scientific worker did not start.");
            started = true;

            var stderrTask = ReadBoundedErrorAsync(process.StandardError, cancellationToken);
            var stagedPayload = await StageInputsAsync(request.Payload, workingDirectory, cancellationToken);
            var envelope = new { requestId = request.RequestId, operation, workingDirectory, payload = stagedPayload };
            await process.StandardInput.WriteLineAsync(JsonSerializer.Serialize(envelope, WireJson));
            await process.StandardInput.FlushAsync(cancellationToken);
            process.StandardInput.Close();

            WorkerResult<TObservation>? terminal = null;
            while (await process.StandardOutput.ReadLineAsync(cancellationToken) is { } line)
            {
                if (line.Length > 8_000_000)
                    return Unobserved<TObservation>(request.RequestId, "oversized-worker-message", "A worker message exceeded the accepted exchange size.");

                using var message = JsonDocument.Parse(line);
                var root = message.RootElement;
                if (!root.TryGetProperty("requestId", out var echoed) || echoed.GetString() != request.RequestId ||
                    !root.TryGetProperty("kind", out var kindElement) || !root.TryGetProperty("payload", out var payload))
                    return Unobserved<TObservation>(request.RequestId, "uncorrelated-worker-message", "A worker message did not match its request.");

                var kind = kindElement.GetString();
                if (kind == "progress")
                {
                    ProgressObserved?.Invoke(request.RequestId, payload.Clone());
                    continue;
                }
                if (terminal is not null)
                    return Unobserved<TObservation>(request.RequestId, "duplicate-worker-terminal", "The worker reported more than one terminal outcome.");

                terminal = kind switch
                {
                    "result" => payload.Deserialize<WorkerResult<TObservation>>(WireJson),
                    "error" => new WorkerResult<TObservation>(
                        request.RequestId,
                        OptionalString(payload, "studyRevisionId"),
                        OptionalString(payload, "attemptId"),
                        OptionalString(payload, "stageId"),
                        WorkerResultStanding.Failed,
                        ImmutableArray<WorkerArtifact>.Empty,
                        null,
                        null,
                        OptionalString(payload, "failureCode"),
                        OptionalString(payload, "failureMessage")),
                    _ => null
                };
                if (terminal is null)
                    return Unobserved<TObservation>(request.RequestId, "unknown-worker-message", "The worker used an unrecognized terminal message.");
            }

            await process.WaitForExitAsync(cancellationToken);
            var stderr = await stderrTask;
            if (process.ExitCode != 0 && terminal?.Standing == WorkerResultStanding.Observed)
                return Unobserved<TObservation>(request.RequestId, "worker-exit-disagreement", $"Worker exit status contradicted its result. {stderr}");
            if (terminal is null || terminal.RequestId != request.RequestId)
                return Unobserved<TObservation>(request.RequestId, "worker-result-unobserved", $"No corresponding terminal worker result was observed. {stderr}");
            if (terminal.Standing != WorkerResultStanding.Observed)
                return terminal;
            if (terminal.Observations is null || terminal.Artifacts.IsDefault)
                return Unobserved<TObservation>(request.RequestId, "incomplete-worker-result", "A successful worker result lacked observations or an artifact account.");

            foreach (var artifact in terminal.Artifacts)
            {
                var path = Path.GetFullPath(Path.IsPathRooted(artifact.Path) ? artifact.Path : Path.Combine(workingDirectory, artifact.Path));
                var relative = Path.GetRelativePath(workingDirectory, path);
                if (Path.IsPathRooted(relative) || relative == ".." || relative.StartsWith($"..{Path.DirectorySeparatorChar}", StringComparison.Ordinal) || !File.Exists(path))
                    return Unobserved<TObservation>(request.RequestId, "invalid-worker-artifact", "A reported artifact is missing or outside its attempt directory.");
                await using var stream = File.OpenRead(path);
                var actualHash = Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken));
                if (!actualHash.Equals(artifact.Sha256, StringComparison.OrdinalIgnoreCase))
                    return Unobserved<TObservation>(request.RequestId, "worker-artifact-mismatch", "A reported artifact did not match its content hash.");
            }
            return terminal;
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            if (started && !process.HasExited)
                process.Kill(entireProcessTree: true);
            return new WorkerResult<TObservation>(request.RequestId, null, null, null,
                WorkerResultStanding.Stopped, ImmutableArray<WorkerArtifact>.Empty, null, null,
                "worker-stopped", "The local work request was stopped before an observed result.");
        }
        catch (Exception exception) when (exception is IOException or InvalidDataException or UnauthorizedAccessException or
            ArgumentException or NotSupportedException or JsonException or InvalidOperationException or System.ComponentModel.Win32Exception)
        {
            if (started && !process.HasExited)
                process.Kill(entireProcessTree: true);
            return Unobserved<TObservation>(request.RequestId, "worker-exchange-failed", exception.Message);
        }
        finally
        {
            if (started && !process.HasExited)
                process.Kill(entireProcessTree: true);
        }
    }

    private static string? OptionalString(JsonElement payload, string property)
        => payload.TryGetProperty(property, out var value) && value.ValueKind == JsonValueKind.String ? value.GetString() : null;

    // Every instruction reads an immutable request-local copy of its molecular
    // inputs. A previous stage's artifact may be reused, but the worker never
    // opens a mutable external path while carrying out this request. Explicit
    // executable paths and OpenMM's built-in force-field selectors are not
    // molecular file inputs and retain their declared meanings.
    private static async Task<JsonNode> StageInputsAsync<TPayload>(
        TPayload payload, string workingDirectory, CancellationToken cancellationToken)
        where TPayload : class
    {
        var node = JsonSerializer.SerializeToNode(payload, WireJson)
            ?? throw new InvalidDataException("The scientific payload was not serializable.");
        var inputDirectory = Path.Combine(workingDirectory, "inputs");
        Directory.CreateDirectory(inputDirectory);
        await StageNodeAsync(node, inputDirectory, cancellationToken);
        return node;
    }

    private static async Task StageNodeAsync(JsonNode node, string inputDirectory, CancellationToken cancellationToken)
    {
        if (node is JsonArray array)
        {
            foreach (var child in array)
                if (child is not null) await StageNodeAsync(child, inputDirectory, cancellationToken);
            return;
        }
        if (node is not JsonObject objectNode) return;

        foreach (var entry in objectNode.ToArray())
        {
            if (entry.Value is null) continue;
            if (entry.Key == "forceFieldFiles" && entry.Value is JsonArray forceFields)
            {
                foreach (var asset in forceFields.OfType<JsonObject>())
                {
                    var path = asset["path"]?.GetValue<string>()
                        ?? throw new InvalidDataException("A force-field asset lacks its identified path.");
                    var hash = asset["sha256"]?.GetValue<string>()
                        ?? throw new InvalidDataException("A force-field asset lacks its content hash.");
                    asset["path"] = await StageFileAsync(path, hash, inputDirectory, cancellationToken);
                }
                continue;
            }
            if (IsMolecularFileProperty(entry.Key) && entry.Value is JsonValue pathValue &&
                pathValue.TryGetValue<string>(out var originalPath) && !string.IsNullOrWhiteSpace(originalPath))
            {
                var hashName = entry.Key switch
                {
                    "sourcePath" => "sourceSha256",
                    "preparedPdbPath" => objectNode.ContainsKey("preparedPdbSha256")
                        ? "preparedPdbSha256" : "preparedSha256",
                    "preparedBondGraphPath" => "preparedBondGraphSha256",
                    "templatePath" => "templateSha256",
                    "coordinateTemplatePath" => "coordinateTemplateSha256",
                    "templateCoordinatePath" => "templateCoordinateSha256",
                    "paePath" => "paeSha256",
                    "paeMappingPath" => "paeMappingSha256",
                    "correspondencePath" => "correspondenceSha256",
                    "coordinatePath" => "coordinateSha256",
                    "topologyCifPath" => "topologyCifSha256",
                    "topologyJsonPath" => "topologyJsonSha256",
                    "systemXmlPath" => "systemXmlSha256",
                    "stateXmlPath" => "stateXmlSha256",
                    "minimizedStateXmlPath" => "minimizedStateXmlSha256",
                    _ => null
                };
                var expectedHash = hashName is not null && objectNode[hashName] is JsonValue hashValue &&
                    hashValue.TryGetValue<string>(out var hash) ? hash : null;
                if ((entry.Key is "topologyCifPath" or "topologyJsonPath" or "systemXmlPath" or
                    "stateXmlPath" or "minimizedStateXmlPath") && string.IsNullOrWhiteSpace(expectedHash))
                    throw new InvalidDataException("An identified stage input lacks its expected content hash.");
                objectNode[entry.Key] = await StageFileAsync(originalPath, expectedHash, inputDirectory, cancellationToken);
            }
            else if (entry.Value is JsonObject or JsonArray)
                await StageNodeAsync(entry.Value, inputDirectory, cancellationToken);
        }
    }

    private static bool IsMolecularFileProperty(string name)
        => name is "sourcePath" or "preparedPdbPath" or "preparedBondGraphPath" or "orientedPdbPath" or
            "topologyPdbPath" or "systemXmlPath" or "stateXmlPath" or
            "topologyCifPath" or "topologyJsonPath" or "minimizedStateXmlPath" or
            "templatePdbPath" or "templatePath" or "templateCoordinatePath" or
            "coordinateTemplatePath" or "coordinatePath" or "paePath" or
            "paeMappingPath" or "correspondencePath";

    private static async Task<string> StageFileAsync(
        string sourcePath, string? expectedHash, string inputDirectory, CancellationToken cancellationToken)
    {
        var fullPath = Path.GetFullPath(sourcePath);
        if (!File.Exists(fullPath)) throw new InvalidDataException($"The required molecular input is absent: {fullPath}");
        await using var input = File.OpenRead(fullPath);
        var hash = Convert.ToHexString(await SHA256.HashDataAsync(input, cancellationToken));
        if (expectedHash is not null && !hash.Equals(expectedHash, StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("A molecular input no longer matches its identified content hash.");
        var extension = Path.GetExtension(fullPath);
        var target = Path.Combine(inputDirectory, hash + extension);
        if (Path.GetFullPath(target) == fullPath) return target;
        if (File.Exists(target))
        {
            await using var existing = File.OpenRead(target);
            var existingHash = Convert.ToHexString(await SHA256.HashDataAsync(existing, cancellationToken));
            if (!existingHash.Equals(hash, StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException("A request-local molecular input does not match its content address.");
        }
        else
        {
            input.Position = 0;
            await using var output = new FileStream(target, FileMode.CreateNew, FileAccess.Write, FileShare.None);
            await input.CopyToAsync(output, cancellationToken);
            await output.FlushAsync(cancellationToken);
        }
        await using var staged = File.OpenRead(target);
        var stagedHash = Convert.ToHexString(await SHA256.HashDataAsync(staged, cancellationToken));
        if (!stagedHash.Equals(hash, StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("A staged molecular input does not match its content address.");
        return target;
    }

    private static async Task<string> ReadBoundedErrorAsync(StreamReader reader, CancellationToken cancellationToken)
    {
        var tail = string.Empty;
        while (await reader.ReadLineAsync(cancellationToken) is { } line)
        {
            tail += line + Environment.NewLine;
            if (tail.Length > 4096)
                tail = tail[^4096..];
        }
        return tail;
    }

    private static WorkerResult<TObservation> Unobserved<TObservation>(string requestId, string code, string message)
        where TObservation : class
        => new(requestId, null, null, null, WorkerResultStanding.Unobserved,
            ImmutableArray<WorkerArtifact>.Empty, null, null, code, message);
}
