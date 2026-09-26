using System.Collections.Immutable;
using System.Reflection;
using System.Text.Json;
using System.Text.Json.Serialization;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class InterruptedHostRestartTests
{
    [Fact]
    public async Task Fresh_root_does_not_resume_or_complete_an_interrupted_attempt_from_workspace_files()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var workspacePath = Path.Combine(fixture.Directory, "interrupted-workspace");
        var cataloguePath = WriteCatalogue(fixture);
        var firstWorker = new AttemptRouteWorker(fixture) { BlockConstruction = true };
        var first = NewRoot(fixture, firstWorker, http, workspacePath, cataloguePath);
        EstablishStudy(first, fixture);

        var started = await Command(first, ActorActionKind.StartPreparation, new { });
        Assert.True(started.Established, started.Reason);
        await firstWorker.ConstructionEntered.Task.WaitAsync(TimeSpan.FromSeconds(4),
            TestContext.Current.CancellationToken);
        var running = first.Snapshot();
        var attemptId = Assert.IsType<AttemptAccount>(running.Attempt).AttemptId;
        Assert.Equal("running", running.Attempt.Status);
        Assert.Empty(running.Stages);

        var attemptPath = Path.Combine(workspacePath, "attempts", attemptId);
        Assert.True(File.Exists(Path.Combine(attemptPath, "constructed.cif")));
        Assert.True(File.Exists(Path.Combine(attemptPath, "state.xml")));
        var stopped = await Command(first, ActorActionKind.StopAttempt, new { attemptId });
        Assert.True(stopped.Established, stopped.Reason);
        await Until(() => first.Snapshot().Attempt?.Status == "stopped");
        Assert.Empty(first.Snapshot().Stages);

        // The interrupted worker left plausible candidate files. Even a partial
        // final-state file and a completion-shaped sidecar are not owner-retained
        // stage and assessment evidence after the process is gone.
        File.WriteAllText(Path.Combine(attemptPath, "minimized.cif"),
            "partial final coordinates from interrupted work");
        File.WriteAllText(Path.Combine(attemptPath, "minimized-state.xml"),
            "partial final state from interrupted work");
        File.WriteAllText(Path.Combine(attemptPath, "stage-observation.json"),
            JsonSerializer.Serialize(new { attemptId, stageId = "plausible-stage",
                standing = "completed" }));

        var freshWorker = new AttemptRouteWorker(fixture);
        var reopened = NewRoot(fixture, freshWorker, http, workspacePath, cataloguePath);
        var state = reopened.Snapshot();
        Assert.Equal(0, state.Revision);
        Assert.Null(state.Attempt);
        Assert.Empty(state.Stages);
        Assert.Null(state.Inspection);
        Assert.False(state.Actions.Single(item =>
            item.Kind == ActorActionKind.ContinueMinimization).Enabled);
        var retained = Workspace(reopened).Snapshot();
        Assert.Null(retained.CurrentAttempt);
        Assert.Null(retained.CurrentExecution);
        Assert.Empty(retained.CompletedStages);
        Assert.Empty(retained.Assessments);
        Assert.Empty(freshWorker.ConstructionRequests);
        Assert.Empty(freshWorker.MinimizationRequests);
    }

    private static string WriteCatalogue(ConstructionFixture fixture)
    {
        var policy = fixture.Policy with
        {
            Water = fixture.Policy.Water with
                { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Sodium = fixture.Policy.Sodium with
                { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Chloride = fixture.Policy.Chloride with
                { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty }
        };
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        var path = Path.Combine(fixture.Directory, "interrupted-workspace-catalogue.json");
        File.WriteAllText(path, JsonSerializer.Serialize(new
        {
            version = "controlled-restart-proof",
            evidenceReferences = new[] { "controlled interrupted attempt proof" },
            lipids = Array.Empty<object>(), proteinChemicalStates = Array.Empty<object>(),
            proteinStructuralPolicies = Array.Empty<object>(), membranePolicies = Array.Empty<object>(),
            placementPolicies = Array.Empty<object>(), placementWitnesses = Array.Empty<object>(),
            preparationPolicies = new[] { policy }, equilibrationQualifications = Array.Empty<object>(),
            ppmVersion = "", ppmExecutableSha256 = "", maximumSourceAtoms = 100000
        }, options));
        return path;
    }

    private static ProductRoot NewRoot(ConstructionFixture fixture, AttemptRouteWorker worker,
        HttpClient http, string workspacePath, string cataloguePath) =>
        new(worker, new ExternalSourceExchange(http), workspacePath, cataloguePath, "",
            () => new ConstructionProviderInstallation(
                fixture.Policy.Construction.ProviderVersion,
                fixture.NativePatchPath, fixture.NativePatchSha));

    private static void EstablishStudy(ProductRoot root, ConstructionFixture fixture)
    {
        var revision = fixture.Revision with { Number = root.Snapshot().Study!.Number + 1 };
        Set(root, "_study", revision);
        Set(root, "_protein", fixture.Protein);
        Set(root, "_membrane", fixture.Membrane);
        Set(root, "_placement", fixture.Placement);
        var revisions = (Dictionary<string, StudyRevision>)typeof(ProductRoot)
            .GetField("_revisions", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        revisions.Add(revision.Id, revision);
        Workspace(root).RetainStudy(revision);
        Assert.True(root.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
    }

    private static LocalRunWorkspace Workspace(ProductRoot root) =>
        (LocalRunWorkspace)typeof(ProductRoot)
            .GetField("_workspace", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;

    private static void Set(ProductRoot root, string name, object value) =>
        typeof(ProductRoot).GetField(name, BindingFlags.Instance | BindingFlags.NonPublic)!
            .SetValue(root, value);

    private static Task<BoundaryOutcome<WorkspaceState>> Command(ProductRoot root,
        ActorActionKind kind, object data) => root.ExecuteAsync(new ActorCommand(kind,
            JsonSerializer.SerializeToElement(data), root.Snapshot().Revision),
            TestContext.Current.CancellationToken);

    private static async Task Until(Func<bool> ready)
    {
        var deadline = DateTimeOffset.UtcNow + TimeSpan.FromSeconds(4);
        while (!ready() && DateTimeOffset.UtcNow < deadline)
            await Task.Delay(10, TestContext.Current.CancellationToken);
        Assert.True(ready(), "The interrupted attempt did not report its actual stopped standing.");
    }
}
