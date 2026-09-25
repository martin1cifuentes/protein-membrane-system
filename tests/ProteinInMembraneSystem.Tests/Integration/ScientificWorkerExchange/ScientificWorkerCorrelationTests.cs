using System.Runtime.Versioning;
using System.Security.Cryptography;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

[SupportedOSPlatform("linux")]
public sealed class ScientificWorkerCorrelationTests
{
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
}
