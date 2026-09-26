using System.Collections.Immutable;
using System.Reflection;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class RetainedWorkspaceRouteTests
{
    [Fact]
    public async Task Completed_optional_stage_account_identifies_the_exact_earlier_minimized_stage()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);
        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var ready = root.Snapshot().Attempt!;
        Assert.True((await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId,
                constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Stages.Length == 1);
        await Until(() => root.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);

        var stages = (Dictionary<string, CompletedStage>)typeof(ProductRoot)
            .GetField("_stages", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        var source = Assert.Single(stages.Values);
        var sourceAssessment = Assert.Single(root.Snapshot().Stages).Assessment!;
        Assert.False(root.Snapshot().Actions.Single(action =>
            action.Kind == ActorActionKind.RequestEquilibration &&
            action.SubjectId == source.Id).Enabled);
        var unqualified = await Command(root, ActorActionKind.RequestEquilibration,
            new { stageId = source.Id });
        Assert.False(unqualified.Established);
        Assert.Contains("validated, applicable", unqualified.Reason ?? "", StringComparison.Ordinal);
        Assert.Single(root.Snapshot().Stages);
        var bundleBytes = System.Text.Encoding.UTF8.GetBytes("controlled earlier-stage bundle");
        var bundleSha = Convert.ToHexString(SHA256.HashData(bundleBytes)).ToLowerInvariant();
        var workspaceRoot = (string)typeof(ProductRoot)
            .GetField("_workspaceRoot", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        var bundlePath = Path.Combine(workspaceRoot, "exports", "earlier-stage.zip");
        Directory.CreateDirectory(Path.GetDirectoryName(bundlePath)!);
        File.WriteAllBytes(bundlePath, bundleBytes);
        var bundles = (Dictionary<string, CompletedStageBundle>)typeof(ProductRoot)
            .GetField("_bundles", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        bundles[source.Id] = new CompletedStageBundle(source.Id, source.Attempt.StudyRevisionId,
            source.Attempt.Id, sourceAssessment.Id, bundlePath, bundleSha, bundleBytes.Length,
            DateTimeOffset.UtcNow);
        var exportAccounts = (Dictionary<string, StageExportAccount>)typeof(ProductRoot)
            .GetField("_exportAccounts", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        exportAccounts[source.Id] = new StageExportAccount(source.Id, sourceAssessment.Id,
            "verified", null, bundleSha, bundleBytes.Length);
        Assert.Equal(ExportDeliveryStanding.Available,
            root.VerifiedExportContent(source.Id, bundleSha).Standing);
        var optionalId = "optional-" + Guid.NewGuid().ToString("N");
        var optional = source with
        {
            Id = optionalId,
            Kind = StageKind.Equilibration,
            Molecule = source.Molecule with { Id = optionalId },
            Observation = source.Observation with { StageId = optionalId,
                Kind = StageKind.Equilibration, Termination = StageTermination.Completed },
            Correspondence = source.Correspondence with { ResultId = optionalId },
            SourceStageId = source.Id,
            Findings = ImmutableArray.Create(new ScientificFinding("optional-only-warning", optionalId,
                "optional-evidence", "The later stage has its own material warning.",
                "Review this later stage.", FindingDisposition.Challenges, true,
                DateTimeOffset.UtcNow)),
            CompletedAt = DateTimeOffset.UtcNow
        };
        Workspace(root).RetainExecution(new StageExecutionState(source.Attempt.Id, optionalId,
            StageKind.Equilibration, StageExecutionStanding.Running, "Controlled optional stage",
            null, DateTimeOffset.UtcNow));
        var constructed = (Dictionary<string, ConstructedExplicitSystem>)typeof(ProductRoot)
            .GetField("_constructedByAttempt", BindingFlags.Instance | BindingFlags.NonPublic)!
            .GetValue(root)!;
        var policies = (Dictionary<string, ApplicablePreparationPolicy>)typeof(ProductRoot)
            .GetField("_policyByAttempt", BindingFlags.Instance | BindingFlags.NonPublic)!
            .GetValue(root)!;
        typeof(ProductRoot).GetMethod("RetainCompleted", BindingFlags.Instance |
            BindingFlags.NonPublic)!.Invoke(root,
            [optional, constructed[source.Attempt.Id], fixture.Protein, fixture.Membrane,
                fixture.Placement, policies[source.Attempt.Id]]);

        var accounts = root.Snapshot().Stages;
        Assert.Null(accounts.Single(item => item.StageId == source.Id).SourceStageId);
        var later = accounts.Single(item => item.StageId == optional.Id);
        Assert.Equal(StageKind.Equilibration, later.Kind);
        Assert.Equal(source.Id, later.SourceStageId);
        Assert.Equal(source.Attempt.Id, later.AttemptId);
        Assert.Equal(sourceAssessment.Id, accounts.Single(item => item.StageId == source.Id).Assessment?.Id);
        Assert.Equal("verified", accounts.Single(item => item.StageId == source.Id).Export?.Status);
        Assert.Equal(bundleBytes, root.VerifiedExportContent(source.Id, bundleSha).Bytes);
    }

    [Fact]
    public async Task Completion_is_exposed_only_with_its_retained_stage_and_assessment()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);
        var witnessed = new List<WorkspaceState>();
        root.Changed += () => witnessed.Add(root.Snapshot());

        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var ready = root.Snapshot();
        var retainedReady = Workspace(root).Snapshot();
        Assert.Equal(ready.Attempt!.AttemptId, retainedReady.CurrentAttempt?.Id);
        Assert.Equal(StageExecutionStanding.ReadyForMinimization, retainedReady.CurrentExecution?.Standing);
        Assert.Empty(retainedReady.CompletedStages);

        Assert.True((await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.Attempt.AttemptId,
                constructedSubjectId = ready.Attempt.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        var running = root.Snapshot();
        Assert.Equal("running", running.Attempt?.Status);
        Assert.Equal(StageExecutionStanding.Running, Workspace(root).Snapshot().CurrentExecution?.Standing);
        var receive = typeof(ProductRoot).GetMethod("SetExecution",
            BindingFlags.Instance | BindingFlags.NonPublic)!;
        receive.Invoke(root, [new StageExecutionState("foreign-attempt", running.Attempt!.CurrentStageId,
            StageKind.Minimization, StageExecutionStanding.Running, "foreign progress", 0.5,
            DateTimeOffset.UtcNow)]);
        receive.Invoke(root, [new StageExecutionState(running.Attempt.AttemptId, "foreign-stage",
            StageKind.Minimization, StageExecutionStanding.Running, "late progress", 0.5,
            DateTimeOffset.UtcNow)]);
        Assert.Equal(running.Revision, root.Snapshot().Revision);
        Assert.Equal(running.Attempt.CurrentStageId, root.Snapshot().Attempt?.CurrentStageId);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Attempt?.Status == "completed");

        var completed = root.Snapshot();
        var stage = Assert.Single(completed.Stages);
        var retained = Workspace(root).Snapshot();
        Assert.Equal(completed.Attempt!.AttemptId, retained.CurrentAttempt?.Id);
        Assert.Equal(stage.StageId, retained.CurrentExecution?.StageId);
        Assert.Equal(StageExecutionStanding.Completed, retained.CurrentExecution?.Standing);
        Assert.Equal(stage.StageId, Assert.Single(retained.CompletedStages).Id);
        Assert.Equal(stage.Assessment?.Id, Assert.Single(retained.Assessments).Id);
        Assert.All(witnessed.Where(item => item.Attempt?.Status == "completed"), item =>
        {
            var actual = Assert.Single(item.Stages);
            Assert.NotNull(actual.Assessment);
            Assert.Equal(actual.StageId, item.Attempt?.CurrentStageId);
        });
        Assert.Single(worker.MinimizationRequests);

        var rootStages = (Dictionary<string, CompletedStage>)typeof(ProductRoot)
            .GetField("_stages", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        var rootConstructed = (Dictionary<string, ConstructedExplicitSystem>)typeof(ProductRoot)
            .GetField("_constructedByAttempt", BindingFlags.Instance | BindingFlags.NonPublic)!
            .GetValue(root)!;
        var real = rootStages[stage.StageId];
        var falseStageId = Guid.NewGuid().ToString("N");
        var unobserved = real with
        {
            Id = falseStageId,
            Molecule = real.Molecule with { Id = falseStageId },
            Observation = real.Observation with { StageId = falseStageId },
            Correspondence = real.Correspondence with { ResultId = falseStageId }
        };
        var rejected = Assert.Throws<TargetInvocationException>(() =>
            typeof(ProductRoot).GetMethod("RetainCompleted", BindingFlags.Instance |
                BindingFlags.NonPublic)!.Invoke(root,
                [unobserved, rootConstructed[real.Attempt.Id], fixture.Protein,
                    fixture.Membrane, fixture.Placement, fixture.Policy]));
        Assert.Contains("Incomplete or mismatched", rejected.InnerException?.Message ?? "",
            StringComparison.Ordinal);
        Assert.Equal(stage.StageId, Assert.Single(root.Snapshot().Stages).StageId);
        Assert.Equal(stage.Assessment?.Id, Assert.Single(root.Snapshot().Stages).Assessment?.Id);
        Assert.Equal("completed", root.Snapshot().Attempt?.Status);
        Assert.Equal(stage.StageId, Assert.Single(Workspace(root).Snapshot().CompletedStages).Id);

        worker.BlockConstruction = true;
        await Until(() => root.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt is { Status: "running" } next &&
            next.AttemptId != completed.Attempt.AttemptId);
        var later = root.Snapshot();
        Assert.Equal(later.Attempt?.AttemptId, Workspace(root).Snapshot().CurrentAttempt?.Id);
        Assert.Equal(stage.StageId, Assert.Single(later.Stages).StageId);
        Assert.Equal(stage.Assessment?.Id, Assert.Single(later.Stages).Assessment?.Id);
        Assert.True((await Command(root, ActorActionKind.StopAttempt,
            new { attemptId = later.Attempt!.AttemptId })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "stopped");
        Assert.Equal(stage.StageId, Assert.Single(root.Snapshot().Stages).StageId);
    }

    [Fact]
    public async Task Historical_stage_keeps_its_origin_and_assessment_after_a_new_study_choice()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);

        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var ready = root.Snapshot().Attempt!;
        Assert.True((await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId,
                constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Stages.Length == 1);
        var original = Assert.Single(root.Snapshot().Stages);
        await Until(() => root.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);

        var proposed = await Command(root, ActorActionKind.ProposeMembrane, new
        {
            upper = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            lower = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            scientificPurpose = "A revised membrane premise"
        });
        Assert.True(proposed.Established, proposed.Reason);
        var changed = await Command(root, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value!.Membrane!.ModelId });
        Assert.True(changed.Established, changed.Reason);
        var current = root.Snapshot();
        Assert.NotEqual(original.StudyRevisionId, current.Study!.Id);
        var historical = Assert.Single(current.Stages);
        Assert.Equal(original.StageId, historical.StageId);
        Assert.Equal(original.AttemptId, historical.AttemptId);
        Assert.Equal(original.StudyRevisionId, historical.StudyRevisionId);
        Assert.Equal(original.Assessment?.Id, historical.Assessment?.Id);
        Assert.Contains("Historical", historical.Summary, StringComparison.Ordinal);
        Assert.Equal(historical.StageId, Assert.Single(Workspace(root).Snapshot().CompletedStages).Id);

        var selected = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = historical.StageId });
        Assert.True(selected.Established, selected.Reason);
        Assert.Equal(historical.StageId, selected.Value!.Inspection?.SubjectId);
        Assert.Equal(historical.StudyRevisionId, selected.Value.Inspection?.StudyRevisionId);
        Assert.Equal(historical.Assessment?.Id, selected.Value.Inspection?.Assessment?.Id);
        Assert.Equal(current.Study.Id, root.Snapshot().Study?.Id);
        Assert.Single(worker.MinimizationRequests);
    }

    [Fact]
    public async Task A_surviving_stage_finishes_under_its_origin_after_the_current_study_changes()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);
        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var ready = root.Snapshot().Attempt!;
        Assert.True((await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId,
                constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        var running = root.Snapshot();
        Assert.Equal("running", running.Attempt?.Status);

        var proposed = await Command(root, ActorActionKind.ProposeMembrane, new
        {
            upper = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            lower = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            scientificPurpose = "A changed premise during surviving work"
        });
        Assert.True(proposed.Established, proposed.Reason);
        Assert.True((await Command(root, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value!.Membrane!.ModelId })).Established);
        Assert.NotEqual(running.Study!.Id, root.Snapshot().Study!.Id);
        Assert.Equal(running.Attempt!.AttemptId, root.Snapshot().Attempt?.AttemptId);
        Assert.Equal("running", root.Snapshot().Attempt?.Status);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Stages.Length == 1);

        var completed = root.Snapshot();
        var stage = Assert.Single(completed.Stages);
        Assert.Equal(running.Attempt.AttemptId, stage.AttemptId);
        Assert.Equal(running.Study.Id, stage.StudyRevisionId);
        Assert.NotEqual(completed.Study!.Id, stage.StudyRevisionId);
        Assert.Contains("Historical", stage.Summary, StringComparison.Ordinal);
        Assert.Equal(stage.StageId, Workspace(root).Snapshot().CurrentExecution?.StageId);
        Assert.Equal(stage.Assessment?.Id, Assert.Single(Workspace(root).Snapshot().Assessments).Id);
        Assert.Single(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Reselecting_a_completed_stage_preserves_derived_evidence_identity_and_values()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);
        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var ready = root.Snapshot().Attempt!;
        Assert.True((await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId,
                constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(root.Snapshot().Stages);

        var first = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = stage.StageId });
        var second = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = stage.StageId });
        Assert.True(first.Established, first.Reason);
        Assert.True(second.Established, second.Reason);
        var original = Assert.IsType<InspectionAccount>(first.Value!.Inspection);
        var repeated = Assert.IsType<InspectionAccount>(second.Value!.Inspection);
        AssertStableDerivedAccount(original, repeated);

        var proposed = await Command(root, ActorActionKind.ProposeMembrane, new
        {
            upper = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            lower = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            scientificPurpose = "A later study premise"
        });
        Assert.True(proposed.Established, proposed.Reason);
        Assert.True((await Command(root, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value!.Membrane!.ModelId })).Established);
        var historical = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = stage.StageId });
        Assert.True(historical.Established, historical.Reason);
        var later = Assert.IsType<InspectionAccount>(historical.Value!.Inspection);
        AssertStableDerivedAccount(original, later);
        Assert.Equal(stage.StudyRevisionId, later.StudyRevisionId);
        Assert.NotEqual(historical.Value.Study!.Id, later.StudyRevisionId);
    }

    private static void AssertStableDerivedAccount(InspectionAccount original,
        InspectionAccount later)
    {
        var originalEvidence = original.Evidence.Where(item =>
            item.Id.StartsWith("inspection-derived-", StringComparison.Ordinal)).ToArray();
        var laterEvidence = later.Evidence.Where(item =>
            item.Id.StartsWith("inspection-derived-", StringComparison.Ordinal)).ToArray();
        Assert.True(originalEvidence.Length >= 2);
        Assert.Equal(originalEvidence.Length, originalEvidence.Select(item => item.Id)
            .Distinct(StringComparer.Ordinal).Count());
        Assert.Equal(originalEvidence, laterEvidence);
        var ids = originalEvidence.Select(item => item.Id).ToHashSet(StringComparer.Ordinal);
        Assert.All(originalEvidence, item => Assert.Contains(original.Metrics,
            metric => metric.EvidenceId == item.Id));
        Assert.Equal(original.Metrics.Where(item => item.EvidenceId is not null &&
                ids.Contains(item.EvidenceId)).ToArray(),
            later.Metrics.Where(item => item.EvidenceId is not null &&
                ids.Contains(item.EvidenceId)).ToArray());
        Assert.Equal(original.Assessment?.Id, later.Assessment?.Id);
        Assert.Equal(original.Assessment?.Qualification, later.Assessment?.Qualification);
        Assert.Equal(original.Assessment?.Reason, later.Assessment?.Reason);
    }

    [Fact]
    public void Issued_structure_url_serves_only_the_selected_bytes()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var root = Product(fixture, new AttemptRouteWorker(fixture), http);
        var path = Path.Combine(fixture.Directory, "root-workspace", "selected.cif");
        File.WriteAllText(path, "selected coordinates");
        var url = (string)typeof(ProductRoot).GetMethod("StructureUrlLocked",
            BindingFlags.Instance | BindingFlags.NonPublic)!.Invoke(root, [path])!;
        var token = url.Split('/')[3].Split('?')[0];
        var served = root.VerifiedStructureContent(token);
        Assert.NotNull(served);
        Assert.Equal(".cif", served.Value.Extension);
        Assert.Equal("selected coordinates", System.Text.Encoding.UTF8.GetString(served.Value.Bytes));

        File.WriteAllText(path, "mutated coordinates");
        Assert.Null(root.VerifiedStructureContent(token));
        File.Delete(path);
        Assert.Null(root.VerifiedStructureContent(token));
    }

    [Fact]
    public async Task Picked_full_system_rows_resolve_only_for_the_selected_unchanged_subject()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var root = Product(fixture, worker, http);
        Assert.True((await Command(root, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => root.Snapshot().Attempt?.Status == "readyForMinimization");
        var constructed = root.Snapshot().Attempt!.Constructed!;
        var selected = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = constructed.SubjectId });
        Assert.True(selected.Established, selected.Reason);
        var token = Token(selected.Value!.Inspection!.StructureUrl!);

        var protein = root.ResolveInspectionAtom(constructed.SubjectId, token, 0);
        Assert.Equal(AtomOriginKind.Source, protein?.Atom.Role);
        Assert.Equal(MoleculeRoleKind.Protein, protein?.Atom.MoleculeRole);
        Assert.NotNull(protein?.Atom.SourceAtomId);
        var upper = root.ResolveInspectionAtom(constructed.SubjectId, token, 3);
        Assert.Equal(AtomOriginKind.Generated, upper?.Atom.Role);
        Assert.Equal(MoleculeRoleKind.Lipid, upper?.Atom.MoleculeRole);
        Assert.Equal(LeafletSide.Upper, upper?.Atom.PhysicalSide);
        Assert.Equal("DMPC", upper?.Atom.GeneratedSpeciesId);
        var lower = root.ResolveInspectionAtom(constructed.SubjectId, token, 9);
        Assert.Equal(LeafletSide.Lower, lower?.Atom.PhysicalSide);
        Assert.Null(root.ResolveInspectionAtom("another-subject", token, 3));
        Assert.Null(root.ResolveInspectionAtom(constructed.SubjectId, "another-token", 3));
        Assert.Null(root.ResolveInspectionAtom(constructed.SubjectId, token, constructed.AtomCount));

        var stageReady = await Command(root, ActorActionKind.ContinueMinimization,
            new { attemptId = root.Snapshot().Attempt!.AttemptId,
                constructedSubjectId = constructed.SubjectId });
        Assert.True(stageReady.Established, stageReady.Reason);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        worker.ReleaseMinimization();
        await Until(() => root.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(root.Snapshot().Stages);
        var stageSelected = await Command(root, ActorActionKind.SelectInspectionSubject,
            new { subjectId = stage.StageId });
        Assert.True(stageSelected.Established, stageSelected.Reason);
        var stageToken = Token(stageSelected.Value!.Inspection!.StructureUrl!);
        Assert.Equal(stage.StageId, root.ResolveInspectionAtom(stage.StageId, stageToken, 3)?.SubjectId);
        Assert.Equal(LeafletSide.Upper,
            root.ResolveInspectionAtom(stage.StageId, stageToken, 3)?.Atom.PhysicalSide);
        Assert.Null(root.ResolveInspectionAtom(constructed.SubjectId, token, 3));
        Assert.Null(root.ResolveInspectionAtom(stage.StageId, token, 3));
        var retainedStages = (Dictionary<string, CompletedStage>)typeof(ProductRoot)
            .GetField("_stages", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(root)!;
        File.AppendAllText(retainedStages[stage.StageId].Molecule.CoordinatePath,
            "changed after inspection selection");
        Assert.Null(root.ResolveInspectionAtom(stage.StageId, stageToken, 3));
    }

    private static string Token(string url) => url.Split('/')[3].Split('?')[0];

    private static ProductRoot Product(ConstructionFixture fixture, AttemptRouteWorker worker,
        HttpClient http)
    {
        var policy = fixture.Policy with
        {
            Water = fixture.Policy.Water with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Sodium = fixture.Policy.Sodium with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Chloride = fixture.Policy.Chloride with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty }
        };
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        var catalogue = Path.Combine(fixture.Directory, "retained-workspace-catalogue.json");
        File.WriteAllText(catalogue, JsonSerializer.Serialize(new
        {
            version = "controlled-slice5",
            evidenceReferences = new[] { "controlled root continuity proof" },
            lipids = Array.Empty<object>(), proteinChemicalStates = Array.Empty<object>(),
            proteinStructuralPolicies = Array.Empty<object>(), membranePolicies = Array.Empty<object>(),
            placementPolicies = Array.Empty<object>(), placementWitnesses = Array.Empty<object>(),
            preparationPolicies = new[] { policy }, equilibrationQualifications = Array.Empty<object>(),
            ppmVersion = "", ppmExecutableSha256 = "", maximumSourceAtoms = 100000
        }, options));
        var root = new ProductRoot(worker, new ExternalSourceExchange(http),
            Path.Combine(fixture.Directory, "root-workspace"), catalogue, "",
            () => new ConstructionProviderInstallation(fixture.Policy.Construction.ProviderVersion,
                fixture.NativePatchPath, fixture.NativePatchSha));
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
        return root;
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
        Assert.True(ready(), "The controlled root did not reach the required observed standing.");
    }
}
