using System.Diagnostics;
using System.Collections.Immutable;
using System.Runtime.Versioning;
using System.Security.Cryptography;
using System.Text.Json;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

[SupportedOSPlatform("linux")]
public sealed class ScientificWorkerCorrelationTests
{
    [Fact]
    public async Task Custom_popc_patch_is_staged_by_digest_while_installed_source_remains_canonical()
    {
        using var directory = new TemporaryDirectory();
        using var fixture = new ConstructionFixture(mixed: true);
        var executable = Path.Combine(directory.Path, "patch-staging-worker");
        File.WriteAllText(executable, """
            #!/usr/bin/env python3
            import json
            import os
            import sys
            request = json.loads(sys.stdin.readline())
            with open(os.path.join(request["workingDirectory"], "observed-request.json"), "w") as stream:
                json.dump(request, stream)
            print(json.dumps({"requestId": request["requestId"], "kind": "error",
                              "payload": {"requestId": request["requestId"],
                                          "failureCode": "controlled-observation",
                                          "failureMessage": "staged paths observed"}}), flush=True)
            """);
        File.SetUnixFileMode(executable, UnixFileMode.UserRead | UnixFileMode.UserWrite |
            UnixFileMode.UserExecute);
        var installed = Path.Combine(directory.Path, "installed-POPC.pdb");
        var derived = Path.Combine(directory.Path, "derived-POPC.pdb");
        File.WriteAllText(installed, "controlled installed provider resource");
        File.WriteAllText(derived, "controlled deletion-only patch");
        var policy = fixture.Policy;
        var protein = fixture.Protein.Molecule;
        var oriented = fixture.Placement.Proposal.OrientedProtein;
        var lipid = fixture.Membrane.SpeciesRepresentations.Single(item => item.SpeciesId == "POPC");
        MolecularRepresentation Ready(MolecularRepresentation item) => item with
        { StereoChecks = item.StereoChecks.IsDefault ? ImmutableArray<MolecularStereoCheck>.Empty : item.StereoChecks };
        ConstructionPayload Payload(string patch, string digest, string? mode,
            string? source, string? sourceDigest) => new(
            fixture.Revision.Id, fixture.Attempt.Id,
            oriented.CoordinatePath, oriented.CoordinateSha256,
            protein.CoordinatePath, protein.CoordinateSha256,
            fixture.Protein.Correspondence, oriented.TopologyPath!, oriented.TopologySha256!,
            Ready(lipid), Ready(policy.Water), Ready(policy.Sodium), Ready(policy.Chloride),
            patch, digest, policy.Construction.ProviderName, policy.Construction.ProviderVersion,
            "POPC", "Na+", "Cl-", 0, 1, 0.15,
            policy.ForceFieldFiles, policy.SystemSettings, policy.LocalStateObservation,
            100000, 180, mode, source, sourceDigest,
            mode is null ? null : ImmutableArray.Create("62", "109"));
        var exchange = new ScientificWorkerExchange(executable, directory.Path);

        var customWork = Path.Combine(directory.Path, "custom-attempt");
        Directory.CreateDirectory(customWork);
        var customPayload = Payload(derived, ConstructionFixture.Hash(derived),
            "popc-62-109-deletion", installed, ConstructionFixture.Hash(installed));
        var customResult = await exchange.ConstructSystemAsync(
            new ScientificWorkRequest<ConstructionPayload>("custom-patch-request", customWork,
                customPayload), TestContext.Current.CancellationToken);
        Assert.True(customResult.FailureCode == "controlled-observation",
            $"{customResult.FailureCode}: {customResult.FailureMessage}");
        using (var observed = JsonDocument.Parse(File.ReadAllText(
                   Path.Combine(customWork, "observed-request.json"))))
        {
            var staged = observed.RootElement.GetProperty("payload");
            var stagedPatch = staged.GetProperty("nativePatchPath").GetString()!;
            Assert.StartsWith(Path.Combine(customWork, "inputs") + Path.DirectorySeparatorChar,
                stagedPatch, StringComparison.Ordinal);
            Assert.Equal(customPayload.NativePatchSha256, ConstructionFixture.Hash(stagedPatch));
            Assert.Equal(installed, staged.GetProperty("nativeSourcePatchPath").GetString());
            Assert.NotEqual(derived, stagedPatch);
            File.WriteAllText(derived, "changed after request completion");
            Assert.Equal(customPayload.NativePatchSha256, ConstructionFixture.Hash(stagedPatch));
        }

        var installedWork = Path.Combine(directory.Path, "installed-attempt");
        Directory.CreateDirectory(installedWork);
        var installedResult = await exchange.ConstructSystemAsync(
            new ScientificWorkRequest<ConstructionPayload>("installed-patch-request", installedWork,
                Payload(installed, ConstructionFixture.Hash(installed), null, null, null)),
            TestContext.Current.CancellationToken);
        Assert.Equal("controlled-observation", installedResult.FailureCode);
        using (var observed = JsonDocument.Parse(File.ReadAllText(
                   Path.Combine(installedWork, "observed-request.json"))))
            Assert.Equal(installed, observed.RootElement.GetProperty("payload")
                .GetProperty("nativePatchPath").GetString());

        var refusedWork = Path.Combine(directory.Path, "refused-attempt");
        Directory.CreateDirectory(refusedWork);
        var refused = await exchange.ConstructSystemAsync(
            new ScientificWorkRequest<ConstructionPayload>("wrong-derived-digest", refusedWork,
                Payload(derived, new string('0', 64), "popc-62-109-deletion",
                    installed, ConstructionFixture.Hash(installed))),
            TestContext.Current.CancellationToken);
        Assert.Equal(WorkerResultStanding.Unobserved, refused.Standing);
        Assert.Equal("worker-exchange-failed", refused.FailureCode);
        Assert.False(File.Exists(Path.Combine(refusedWork, "observed-request.json")));

        var aliasedWork = Path.Combine(directory.Path, "aliased-attempt");
        Directory.CreateDirectory(aliasedWork);
        var aliased = await exchange.ConstructSystemAsync(
            new ScientificWorkRequest<ConstructionPayload>("aliased-source-derived", aliasedWork,
                Payload(installed, ConstructionFixture.Hash(installed), "popc-62-109-deletion",
                    installed, ConstructionFixture.Hash(installed))),
            TestContext.Current.CancellationToken);
        Assert.Equal(WorkerResultStanding.Unobserved, aliased.Standing);
        Assert.Equal("worker-exchange-failed", aliased.FailureCode);
        Assert.False(File.Exists(Path.Combine(aliasedWork, "observed-request.json")));
    }

    [Fact]
    public async Task Optional_step_progress_is_request_bound_and_cannot_follow_a_terminal_message()
    {
        using var directory = new TemporaryDirectory();
        using var fixture = new OptionalEquilibrationFixture();
        var executable = Path.Combine(directory.Path, "optional-progress-worker");
        File.WriteAllText(executable, """
            #!/usr/bin/env python3
            import json
            import sys
            request = json.loads(sys.stdin.readline())
            source = request["payload"]
            detail = {"window": "unrestrained", "completedSteps": 1, "requestedSteps": 3}
            base = {"operation": "equilibrate", "stage": "equilibrationProgress",
                    "studyRevisionId": source["studyRevisionId"],
                    "attemptId": source["attemptId"], "stageId": source["stageId"], "detail": detail}
            print(json.dumps({"requestId": request["requestId"], "kind": "progress",
                              "payload": {**base, "attemptId": "foreign-attempt"}}), flush=True)
            print(json.dumps({"requestId": request["requestId"], "kind": "progress",
                              "payload": base}), flush=True)
            print(json.dumps({"requestId": request["requestId"], "kind": "error",
                              "payload": {"requestId": request["requestId"],
                                          "failureCode": "controlled-failure",
                                          "failureMessage": "controlled terminal"}}), flush=True)
            print(json.dumps({"requestId": request["requestId"], "kind": "progress",
                              "payload": {**base, "detail": {**detail, "completedSteps": 2}}}), flush=True)
            """);
        File.SetUnixFileMode(executable, UnixFileMode.UserRead | UnixFileMode.UserWrite |
            UnixFileMode.UserExecute);
        var exchange = new ScientificWorkerExchange(executable, directory.Path);
        var progress = new OptionalProgressCollector();
        var molecule = fixture.Minimized.Molecule;
        var request = new ScientificWorkRequest<EquilibrationPayload>("exact-optional-request",
            fixture.Directory, new EquilibrationPayload(fixture.Minimized.Attempt.StudyRevisionId,
                fixture.Minimized.Attempt.Id, "optional-stage", fixture.Minimized.Id,
                molecule.CoordinatePath, molecule.CoordinateSha256,
                molecule.TopologyPath!, molecule.TopologySha256!,
                molecule.SystemXmlPath!, molecule.SystemXmlSha256!,
                molecule.StateXmlPath!, molecule.StateXmlSha256!, fixture.Protocol.Id,
                ImmutableArray<int>.Empty, ImmutableArray<int>.Empty, ImmutableArray<int>.Empty,
                ImmutableArray<ResolvedEquilibrationObservable>.Empty, fixture.Protocol));

        var result = await exchange.EquilibrateAsync(request, progress,
            TestContext.Current.CancellationToken);

        Assert.Equal(WorkerResultStanding.Unobserved, result.Standing);
        Assert.Equal("late-worker-progress", result.FailureCode);
        var observed = Assert.Single(progress.Items);
        Assert.Equal("optional-stage", observed.StageId);
        Assert.Equal(fixture.Minimized.Attempt.Id, observed.AttemptId);
        Assert.Equal(1, observed.CompletedSteps);
    }

    [Fact]
    public async Task Transport_rejects_a_terminal_result_that_does_not_echo_the_request_id()
    {
        using var directory = new TemporaryDirectory();
        var executable = Path.Combine(directory.Path, "controlled-worker");
        File.WriteAllText(executable, """
            #!/usr/bin/env python3
            import json
            import sys
            request = json.loads(sys.stdin.readline())
            print(json.dumps({"requestId": request["requestId"], "kind": "progress",
                              "payload": {"stage": "started"}}), flush=True)
            print(json.dumps({"requestId": "another-request", "kind": "result",
                              "payload": {"requestId": "another-request", "standing": "observed",
                                          "artifacts": [], "observations": {}}}), flush=True)
            """);
        File.SetUnixFileMode(executable, UnixFileMode.UserRead | UnixFileMode.UserWrite |
            UnixFileMode.UserExecute);
        var work = Path.Combine(directory.Path, "attempt");
        Directory.CreateDirectory(work);
        var source = Path.Combine(work, "source.pdb");
        File.WriteAllText(source, "ATOM\n");
        var digest = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(source))).ToLowerInvariant();
        var exchange = new ScientificWorkerExchange(executable, directory.Path);
        var progressIds = new List<string>();
        exchange.ProgressObserved += (id, _) => progressIds.Add(id);
        var request = new ScientificWorkRequest<SourceInspectionPayload>("expected-request", work,
            new SourceInspectionPayload(source, digest, 100, SourceRouteKind.Upload, null));

        var result = await exchange.InspectSourceAsync(request, TestContext.Current.CancellationToken);

        Assert.Equal(WorkerResultStanding.Unobserved, result.Standing);
        Assert.Equal("expected-request", result.RequestId);
        Assert.Equal("uncorrelated-worker-message", result.FailureCode);
        Assert.Equal(["expected-request"], progressIds);
        Assert.Null(result.Observations);
        Assert.Empty(result.Artifacts);
    }

    [Fact]
    public async Task Transport_rejects_an_artifact_outside_the_identified_attempt_directory()
    {
        using var directory = new TemporaryDirectory();
        var executable = Path.Combine(directory.Path, "outside-artifact-worker");
        File.WriteAllText(executable, """
            #!/usr/bin/env python3
            import hashlib
            import json
            import os
            import sys
            request = json.loads(sys.stdin.readline())
            outside = os.path.join(os.getcwd(), "outside.pdb")
            with open(outside, "rb") as stream:
                digest = hashlib.sha256(stream.read()).hexdigest()
            print(json.dumps({"requestId": request["requestId"], "kind": "result",
                              "payload": {"requestId": request["requestId"],
                                          "standing": "observed",
                                          "artifacts": [{"role": "preparedPdb", "path": outside,
                                                         "sha256": digest}],
                                          "observations": {"sourceFormat": "pdb", "models": []}}}),
                  flush=True)
            """);
        File.SetUnixFileMode(executable, UnixFileMode.UserRead | UnixFileMode.UserWrite |
            UnixFileMode.UserExecute);
        File.WriteAllText(Path.Combine(directory.Path, "outside.pdb"), "outside attempt");
        var work = Path.Combine(directory.Path, "attempt");
        Directory.CreateDirectory(work);
        var source = Path.Combine(work, "source.pdb");
        File.WriteAllText(source, "ATOM\n");
        var digest = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(source))).ToLowerInvariant();
        var exchange = new ScientificWorkerExchange(executable, directory.Path);
        var request = new ScientificWorkRequest<SourceInspectionPayload>("identified-attempt", work,
            new SourceInspectionPayload(source, digest, 100, SourceRouteKind.Upload, null));

        var result = await exchange.InspectSourceAsync(request, TestContext.Current.CancellationToken);

        Assert.Equal(WorkerResultStanding.Unobserved, result.Standing);
        Assert.Equal("identified-attempt", result.RequestId);
        Assert.Equal("invalid-worker-artifact", result.FailureCode);
        Assert.Empty(result.Artifacts);
        Assert.Null(result.Observations);
    }

    [Fact]
    public async Task Cancelling_a_blocked_request_kills_its_worker_process_without_a_result()
    {
        using var directory = new TemporaryDirectory();
        var executable = Path.Combine(directory.Path, "blocked-worker");
        File.WriteAllText(executable, """
            #!/usr/bin/env python3
            import json
            import os
            import sys
            import time
            request = json.loads(sys.stdin.readline())
            print(json.dumps({"requestId": request["requestId"], "kind": "progress",
                              "payload": {"stage": "blocked", "pid": os.getpid()}}), flush=True)
            time.sleep(60)
            """);
        File.SetUnixFileMode(executable, UnixFileMode.UserRead | UnixFileMode.UserWrite |
            UnixFileMode.UserExecute);
        var work = Path.Combine(directory.Path, "attempt");
        Directory.CreateDirectory(work);
        var source = Path.Combine(work, "source.pdb");
        File.WriteAllText(source, "ATOM\n");
        var digest = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(source))).ToLowerInvariant();
        var exchange = new ScientificWorkerExchange(executable, directory.Path);
        var started = new TaskCompletionSource<int>(TaskCreationOptions.RunContinuationsAsynchronously);
        exchange.ProgressObserved += (requestId, payload) =>
        {
            if (requestId == "blocked-request" && payload.GetProperty("stage").GetString() == "blocked")
                started.TrySetResult(payload.GetProperty("pid").GetInt32());
        };
        var request = new ScientificWorkRequest<SourceInspectionPayload>("blocked-request", work,
            new SourceInspectionPayload(source, digest, 100, SourceRouteKind.Upload, null));
        using var cancellation = CancellationTokenSource.CreateLinkedTokenSource(
            TestContext.Current.CancellationToken);
        cancellation.CancelAfter(TimeSpan.FromSeconds(10));
        var pending = exchange.InspectSourceAsync(request, cancellation.Token);
        int pid;
        try
        {
            pid = await started.Task.WaitAsync(TimeSpan.FromSeconds(3),
                TestContext.Current.CancellationToken);
        }
        catch
        {
            cancellation.Cancel();
            await pending.WaitAsync(TimeSpan.FromSeconds(3), TestContext.Current.CancellationToken);
            throw;
        }
        cancellation.Cancel();
        var result = await pending.WaitAsync(TimeSpan.FromSeconds(5),
            TestContext.Current.CancellationToken);

        Assert.Equal("blocked-request", result.RequestId);
        Assert.Equal(WorkerResultStanding.Stopped, result.Standing);
        Assert.Equal("worker-stopped", result.FailureCode);
        Assert.Null(result.Observations);
        Assert.Empty(result.Artifacts);
        var deadline = DateTime.UtcNow + TimeSpan.FromSeconds(3);
        while (WorkerStillAlive(pid) && DateTime.UtcNow < deadline)
            await Task.Delay(20, TestContext.Current.CancellationToken);
        Assert.False(WorkerStillAlive(pid), $"Cancelled worker PID {pid} remains running.");

        static bool WorkerStillAlive(int pid)
        {
            try { using var process = Process.GetProcessById(pid); return !process.HasExited; }
            catch (ArgumentException) { return false; }
        }
    }
}

internal sealed class OptionalProgressCollector : IProgress<EquilibrationWorkProgress>
{
    public List<EquilibrationWorkProgress> Items { get; } = [];
    public void Report(EquilibrationWorkProgress value) => Items.Add(value);
}
