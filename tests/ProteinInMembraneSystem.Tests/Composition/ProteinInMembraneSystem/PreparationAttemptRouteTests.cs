using System.Collections.Immutable;
using System.Diagnostics;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using ExportOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport.CompletedStageExport;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class PreparationAttemptRouteTests
{
    [Fact]
    public async Task Wrong_installed_OpenMM_full_version_withholds_start_before_attempt_or_worker()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var calls = 0;
        var product = Product(fixture, worker, http, () =>
        {
            calls++;
            return new ConstructionProviderInstallation(
                fixture.Policy.Construction.ProviderVersion + ".other",
                fixture.NativePatchPath, fixture.NativePatchSha);
        }, expectStartAvailable: false);

        Assert.Equal(1, calls);
        Assert.False(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.Equal(1, calls);
        var refused = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.False(refused.Established);
        Assert.Contains("OpenMM", refused.Reason ?? "", StringComparison.Ordinal);
        Assert.Equal(1, calls);
        Assert.Null(product.Snapshot().Attempt);
        Assert.Empty(worker.ConstructionRequests);
    }

    [Fact]
    public async Task OpenMM_full_version_drift_is_refused_at_admission_before_worker()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var calls = 0;
        var installed = new ConstructionProviderInstallation(
            fixture.Policy.Construction.ProviderVersion, fixture.NativePatchPath,
            fixture.NativePatchSha);
        var product = Product(fixture, worker, http, () =>
        {
            calls++;
            return installed;
        });
        Assert.Equal(1, calls);
        Assert.True(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.Equal(1, calls);

        installed = installed with { FullVersion = installed.FullVersion + ".changed" };
        var refused = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.False(refused.Established);
        Assert.Contains("changed", refused.Reason ?? "", StringComparison.Ordinal);
        Assert.Equal(2, calls);
        Assert.False(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.Null(product.Snapshot().Attempt);
        Assert.Empty(worker.ConstructionRequests);
    }

    [Fact]
    public async Task Start_requires_the_exact_current_adopted_supported_placement_identity()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        var stale = fixture.Revision with { AdoptedPlacementProposalId = "another-proposal" };
        Set(product, "_study", stale);
        Assert.False(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        var refused = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.False(refused.Established);
        Assert.Empty(worker.ConstructionRequests);
        Assert.Empty(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Root_exposes_native_candidate_before_exact_continue_and_reports_completed_stage_separately()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);

        var started = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.True(started.Established, started.Reason);
        await WaitForReady(product);
        var ready = product.Snapshot();
        var attempt = Assert.IsType<AttemptAccount>(ready.Attempt);
        Assert.Equal("readyForMinimization", attempt.Status);
        Assert.Equal(fixture.Revision.Id, attempt.StudyRevisionId);
        Assert.Equal(fixture.Policy.Id, attempt.PolicyId);
        Assert.Equal(fixture.Policy.Version, attempt.PolicyVersion);
        var accepted = (PreparationAttempt)typeof(ProductRoot)
            .GetField("_currentAttempt", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        Assert.Equal(fixture.Revision.Id, accepted.StudyRevisionId);
        Assert.Equal(fixture.Protein.Id, accepted.ProteinId);
        Assert.Equal(fixture.Membrane.Id, accepted.MembraneId);
        Assert.Equal(fixture.Placement.Id, accepted.PlacementId);
        Assert.Equal(PreparationPolicyFingerprint.Compute(fixture.Policy),
            accepted.PolicyFingerprintSha256);
        Assert.Equal(fixture.Policy.ForceFieldFiles, accepted.ForceFieldFiles);
        Assert.Equal(fixture.Policy.Construction.ProviderVersion, accepted.ConstructionProviderVersion);
        Assert.Equal(fixture.NativePatchSha, accepted.NativePatchSha256);
        Assert.NotNull(attempt.Derivation);
        var constructed = Assert.IsType<ConstructedSystemAccount>(attempt.Constructed);
        Assert.Equal(attempt.AttemptId, constructed.AttemptId);
        Assert.Equal(attempt.Derivation.LipidCounts, constructed.AchievedComposition);
        Assert.Equal(attempt.Derivation.CellAngstrom, constructed.CellAngstrom);
        Assert.Equal(80.0, constructed.CellAngstrom[0]);
        Assert.Equal(81.0, constructed.CellAngstrom[1]);
        Assert.Equal(100.0, constructed.CellAngstrom[2]);
        Assert.Equal(367, constructed.WaterCount);
        Assert.Equal(2, constructed.SodiumCount);
        Assert.Equal(1, constructed.ChlorideCount);
        Assert.Empty(ready.Stages);
        var continueAction = ready.Actions.Single(item => item.Kind == ActorActionKind.ContinueMinimization);
        Assert.True(continueAction.Enabled);
        Assert.Equal(attempt.AttemptId, continueAction.SubjectId);
        Assert.False(ready.Actions.Single(item => item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.False((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        Assert.Single(worker.ConstructionRequests);
        Assert.Equal(fixture.NativePatchSha, worker.ConstructionRequests[0].Payload.NativePatchSha256);
        Assert.Equal(fixture.Policy.Construction.ProviderVersion,
            worker.ConstructionRequests[0].Payload.ProviderVersion);
        Assert.Equal(fixture.Policy.Construction.MinimumPaddingNanometers,
            worker.ConstructionRequests[0].Payload.MinimumPaddingNanometers);
        Assert.Empty(worker.MinimizationRequests);
        Assert.Empty(worker.ObservationRequests);
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = constructed.SubjectId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal("constructedSystem", inspected.Value!.Inspection?.RepresentationKind);
        Assert.NotNull(inspected.Value.Inspection?.StructureUrl);

        var continued = await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = attempt.AttemptId, constructedSubjectId = constructed.SubjectId });
        Assert.True(continued.Established, continued.Reason);
        await WaitForMinimization(product, worker);
        Assert.Equal("running", product.Snapshot().Attempt?.Status);
        Assert.Single(worker.MinimizationRequests);

        var boundFingerprint = accepted.PolicyFingerprintSha256;
        File.WriteAllText(fixture.Policy.ForceFieldFiles[0].Path,
            "changed after the complete system was constructed");
        Assert.Equal(boundFingerprint, ((PreparationAttempt)typeof(ProductRoot)
            .GetField("_currentAttempt", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!).PolicyFingerprintSha256);

        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var completed = product.Snapshot();
        var stage = Assert.Single(completed.Stages);
        Assert.Equal("completed", stage.Status);
        Assert.Equal(attempt.AttemptId, stage.AttemptId);
        Assert.Equal(fixture.Revision.Id, stage.StudyRevisionId);
        Assert.Equal(StageKind.Minimization, stage.Kind);
        Assert.Equal(StageTermination.Converged, stage.Observation?.Termination);
        Assert.Equal(PreparationQualification.Indeterminate, stage.Assessment?.Qualification);
        Assert.True(stage.Assessment?.CurrentlyApplicable);
        Assert.Equal(constructed.SubjectId, stage.Constructed?.SubjectId);
        Assert.Equal(constructed.AchievedComposition, stage.Constructed?.AchievedComposition);
        Assert.Equal("completed", completed.Attempt?.Status);
        Assert.Single(worker.ObservationRequests);
        Assert.False(completed.Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
    }

    [Fact]
    public async Task Completed_minimization_state_serializes_for_the_browser()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);

        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        using var document = JsonDocument.Parse(JsonSerializer.Serialize(product.Snapshot(), options));
        var stage = document.RootElement.GetProperty("stages")[0];
        Assert.Equal("completed", stage.GetProperty("status").GetString());
        var observation = stage.GetProperty("observation");
        Assert.Equal(JsonValueKind.Array, observation.GetProperty("equilibrationAssessments").ValueKind);
        Assert.Equal(JsonValueKind.Array, observation.GetProperty("equilibrationSamples").ValueKind);
        Assert.Equal(0, observation.GetProperty("equilibrationAssessments").GetArrayLength());
        Assert.Equal(0, observation.GetProperty("equilibrationSamples").GetArrayLength());
    }

    [Fact]
    public async Task Completed_assessed_stage_offers_export_but_unknown_stage_is_refused()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);

        var state = product.Snapshot();
        var stage = Assert.Single(state.Stages);
        Assert.True(stage.Assessment?.CurrentlyApplicable);
        var export = Assert.Single(state.Actions.Where(item =>
            item.Kind == ActorActionKind.ExportStage && item.SubjectId == stage.StageId));
        Assert.True(export.Enabled);

        var direct = await Command(product, ActorActionKind.ExportStage,
            new { stageId = "unknown-completed-stage" });
        Assert.False(direct.Established);
        Assert.Contains("identified completed", direct.Reason ?? "", StringComparison.Ordinal);
        Assert.Equal(stage.StageId, Assert.Single(product.Snapshot().Stages).StageId);
    }

    [Fact]
    public async Task Constructed_running_and_stopped_partials_never_offer_or_deliver_stage_export()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);

        Assert.Empty(product.Snapshot().Stages);
        Assert.Empty(product.Snapshot().Actions.Where(item => item.Kind == ActorActionKind.ExportStage));
        Assert.False((await Command(product, ActorActionKind.ExportStage,
            new { stageId = "not-yet-constructed" })).Established);

        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        var candidateId = product.Snapshot().Attempt!.Constructed!.SubjectId;
        await RefusePartial(candidateId, "readyForMinimization");

        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        var unfinishedStageId = Assert.Single(worker.MinimizationRequests).Payload.StageId;
        await RefusePartial(unfinishedStageId, "running");

        Assert.True((await Command(product, ActorActionKind.StopAttempt,
            new { attemptId = product.Snapshot().Attempt!.AttemptId })).Established);
        await Until(() => product.Snapshot().Attempt?.Status == "stopped");
        await RefusePartial(unfinishedStageId, "stopped");
        Assert.Empty(worker.ExportRequests);

        async Task RefusePartial(string selectedId, string expectedStanding)
        {
            var before = product.Snapshot();
            Assert.Equal(expectedStanding, before.Attempt?.Status);
            Assert.Empty(before.Stages);
            Assert.Empty(before.Actions.Where(item => item.Kind == ActorActionKind.ExportStage));
            var refused = await Command(product, ActorActionKind.ExportStage,
                new { stageId = selectedId });
            Assert.False(refused.Established);
            Assert.Empty(worker.ExportRequests);
            Assert.Empty(product.Snapshot().Stages);
        }
    }

    [Fact]
    public async Task Unobserved_minimization_partial_cannot_be_exported_as_a_completed_stage()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture)
        { MinimizationStanding = WorkerResultStanding.Unobserved };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        var unfinishedStageId = Assert.Single(worker.MinimizationRequests).Payload.StageId;
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Attempt?.Status == "unobserved");

        Assert.Empty(product.Snapshot().Stages);
        Assert.Empty(product.Snapshot().Actions.Where(item => item.Kind == ActorActionKind.ExportStage));
        Assert.False((await Command(product, ActorActionKind.ExportStage,
            new { stageId = unfinishedStageId })).Established);
        Assert.Empty(worker.ExportRequests);
    }

    [Fact]
    public async Task Stale_actor_revision_and_absent_current_assessment_each_refuse_exact_export()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        Assert.True(stage.Assessment?.CurrentlyApplicable);

        var staleRevision = product.Snapshot().Revision - 1;
        var stale = await product.ExecuteAsync(new ActorCommand(ActorActionKind.ExportStage,
            JsonSerializer.SerializeToElement(new { stageId = stage.StageId }), staleRevision),
            TestContext.Current.CancellationToken);
        Assert.False(stale.Established);
        Assert.Empty(worker.ExportRequests);
        Assert.Null(Assert.Single(product.Snapshot().Stages).Export);

        var assessments = (Dictionary<string, PreparationAssessmentResult>)typeof(ProductRoot)
            .GetField("_stageAssessments", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        Assert.True(assessments.Remove(stage.StageId));
        var withoutAssessment = product.Snapshot();
        Assert.Equal(stage.Assessment?.Id, Assert.Single(withoutAssessment.Stages).Assessment?.Id);
        // The retained historical assessment can remain visible, but it no
        // longer authorizes a current export when the root has none to bind.
        Assert.False(withoutAssessment.Actions.Single(item =>
            item.Kind == ActorActionKind.ExportStage && item.SubjectId == stage.StageId).Enabled);
        var refused = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.False(refused.Established);
        Assert.Empty(worker.ExportRequests);
        Assert.Null(Assert.Single(product.Snapshot().Stages).Export);
    }

    [Fact]
    public async Task Export_keeps_the_exact_indeterminate_stage_and_repairs_a_tampered_published_bundle()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var assessment = Assert.IsType<PreparationAssessmentResult>(stage.Assessment);
        Assert.Equal(PreparationQualification.Indeterminate, assessment.Qualification);

        var exportDirectory = Path.Combine(fixture.Directory, "root-workspace", "exports",
            stage.StageId + "-" + assessment.Id);
        Directory.CreateDirectory(exportDirectory);
        var orphan = Path.Combine(exportDirectory,
            $"completed-stage-{stage.StageId}-{assessment.Id}.zip");
        File.WriteAllText(orphan, "orphan from an interrupted prior delivery");

        var delivered = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(delivered.Established, delivered.Reason);
        var export = Assert.IsType<StageExportAccount>(Assert.Single(delivered.Value!.Stages).Export);
        Assert.Equal("verified", export.Status);
        Assert.Equal(stage.StageId, export.StageId);
        Assert.Equal(assessment.Id, export.AssessmentId);
        var request = Assert.Single(worker.ExportRequests);
        Assert.Equal(stage.StageId, request.Payload.StageId);
        Assert.Equal(stage.AttemptId, request.Payload.AttemptId);
        Assert.Equal(stage.StudyRevisionId, request.Payload.StudyRevisionId);
        Assert.Equal(fixture.Policy.ExportCoordinateReadBackToleranceAngstrom,
            request.Payload.CoordinateReadBackToleranceAngstrom);
        var content = product.VerifiedExportContent(stage.StageId, export.Sha256);
        Assert.Equal(ExportDeliveryStanding.Available, content.Standing);
        Assert.Equal(export.Sha256, Convert.ToHexString(SHA256.HashData(content.Bytes!)).ToLowerInvariant());
        var repeat = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(repeat.Established, repeat.Reason);
        Assert.Equal(export.Sha256, Assert.Single(repeat.Value!.Stages).Export?.Sha256);
        Assert.Single(worker.ExportRequests);
        Assert.Equal(content.Bytes,
            product.VerifiedExportContent(stage.StageId, export.Sha256).Bytes);
        using (var archive = new ZipArchive(new MemoryStream(content.Bytes!), ZipArchiveMode.Read))
        {
            Assert.Equal(5, archive.Entries.Count);
            Assert.NotNull(archive.GetEntry("structure.cif"));
            Assert.NotNull(archive.GetEntry("topology.json"));
            Assert.NotNull(archive.GetEntry("system.xml"));
            Assert.NotNull(archive.GetEntry("state.xml"));
            using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
            var root = manifest.RootElement;
            Assert.Equal(stage.StageId, root.GetProperty("stage").GetProperty("id").GetString());
            Assert.Equal(stage.AttemptId, root.GetProperty("attempt").GetProperty("id").GetString());
            Assert.Equal(assessment.Id, root.GetProperty("assessment").GetProperty("id").GetString());
            Assert.Equal("Indeterminate", root.GetProperty("assessment").GetProperty("qualification").GetString());
            Assert.Equal(fixture.Protein.Intended.Source.Sha256,
                root.GetProperty("lineage").GetProperty("sourceCoordinateSha256").GetString());
            Assert.Equal(fixture.Placement.Proposal.OrientedProtein.CoordinateSha256,
                root.GetProperty("lineage").GetProperty("orientedCoordinateSha256").GetString());
            Assert.Equal(fixture.Membrane.Id, root.GetProperty("membrane").GetProperty("id").GetString());
            Assert.Equal(fixture.Placement.Id, root.GetProperty("placement").GetProperty("id").GetString());
            Assert.Equal(3, root.GetProperty("forceFieldAssets").GetArrayLength());
            Assert.True(root.GetProperty("attribution").GetProperty("incorporatedData").GetArrayLength() >= 5);
        }
        Assert.Equal(ExportDeliveryStanding.IdentityChanged,
            product.VerifiedExportContent(stage.StageId, new string('0', 64)).Standing);
        var published = (Dictionary<string, CompletedStageBundle>)typeof(ProductRoot)
            .GetField("_bundles", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        File.AppendAllText(published[stage.StageId].BundlePath, "tampered");
        Assert.Equal(ExportDeliveryStanding.Corrupt,
            product.VerifiedExportContent(stage.StageId, export.Sha256).Standing);
        var failed = Assert.Single(product.Snapshot().Stages);
        Assert.Equal("failed", failed.Export?.Status);
        Assert.Equal(assessment.Id, failed.Assessment?.Id);
        Assert.Equal("completed", failed.Status);

        var retry = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(retry.Established, retry.Reason);
        var renewed = Assert.Single(retry.Value!.Stages);
        Assert.Equal("verified", renewed.Export?.Status);
        Assert.Equal(assessment.Id, renewed.Assessment?.Id);
        Assert.Equal(stage.AttemptId, renewed.AttemptId);
        Assert.Equal(2, worker.ExportRequests.Count);
        Assert.Equal(ExportDeliveryStanding.Available,
            product.VerifiedExportContent(stage.StageId, renewed.Export?.Sha256).Standing);
    }

    [Fact]
    public async Task Export_failure_preserves_stage_and_refuses_stale_assessment_before_worker()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture) { FailNextExport = true };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var initialAssessment = stage.Assessment!.Id;
        var failed = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.False(failed.Established);
        Assert.Equal("failed", Assert.Single(product.Snapshot().Stages).Export?.Status);
        Assert.Equal(initialAssessment, Assert.Single(product.Snapshot().Stages).Assessment?.Id);
        Assert.Equal(ExportDeliveryStanding.Absent,
            product.VerifiedExportContent(stage.StageId, new string('0', 64)).Standing);

        var assessments = (Dictionary<string, PreparationAssessmentResult>)typeof(ProductRoot)
            .GetField("_stageAssessments", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        assessments[stage.StageId] = assessments[stage.StageId] with { CurrentlyApplicable = false };
        Assert.False(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.ExportStage && item.SubjectId == stage.StageId).Enabled);
        var stale = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.False(stale.Established);
        Assert.Single(worker.ExportRequests);
    }

    [Fact]
    public async Task Two_completed_stages_export_under_their_own_attempts_and_assessments()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var first = Assert.Single(product.Snapshot().Stages);
        await Until(() => product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);

        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await Until(() => product.Snapshot().Stages.Length == 2);
        var stages = product.Snapshot().Stages;
        var second = stages.Single(item => item.StageId != first.StageId);
        Assert.NotEqual(first.AttemptId, second.AttemptId);

        foreach (var stage in new[] { first, second })
        {
            var delivered = await Command(product, ActorActionKind.ExportStage,
                new { stageId = stage.StageId });
            Assert.True(delivered.Established, delivered.Reason);
            var export = delivered.Value!.Stages.Single(item => item.StageId == stage.StageId).Export!;
            var content = product.VerifiedExportContent(stage.StageId, export.Sha256);
            Assert.Equal(ExportDeliveryStanding.Available, content.Standing);
            using var archive = new ZipArchive(new MemoryStream(content.Bytes!), ZipArchiveMode.Read);
            using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
            Assert.Equal(stage.StageId,
                manifest.RootElement.GetProperty("stage").GetProperty("id").GetString());
            Assert.Equal(stage.AttemptId,
                manifest.RootElement.GetProperty("attempt").GetProperty("id").GetString());
            Assert.Equal(export.AssessmentId,
                manifest.RootElement.GetProperty("assessment").GetProperty("id").GetString());
        }
        Assert.Equal(2, worker.ExportRequests.Count);
        Assert.Equal(2, worker.ExportRequests.Select(item => item.Payload.StageId).Distinct().Count());
    }

    [Fact]
    public async Task Changed_export_worker_artifact_is_a_delivery_fault_and_retry_preserves_science()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture) { ChangeNextExportArtifactAfterHash = true };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var initialAssessment = stage.Assessment!.Id;
        var failed = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.False(failed.Established);
        Assert.Equal("failed", Assert.Single(product.Snapshot().Stages).Export?.Status);
        Assert.Equal(ExportDeliveryStanding.Absent,
            product.VerifiedExportContent(stage.StageId, new string('0', 64)).Standing);
        Assert.Empty(Directory.GetFiles(Path.Combine(fixture.Directory, "root-workspace", "exports"), "*.zip",
            SearchOption.AllDirectories));
        var retry = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(retry.Established, retry.Reason);
        Assert.Equal(initialAssessment, retry.Value!.Stages.Single().Assessment?.Id);
        Assert.Equal("verified", retry.Value.Stages.Single().Export?.Status);
        Assert.Equal(2, worker.ExportRequests.Count);
    }

    [Fact]
    public async Task Final_publish_conflict_leaves_no_deliverable_and_retry_succeeds_for_same_stage()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var assessmentId = stage.Assessment!.Id;
        var destination = Path.Combine(fixture.Directory, "root-workspace", "exports",
            stage.StageId + "-" + assessmentId,
            $"completed-stage-{stage.StageId}-{assessmentId}.zip");
        Directory.CreateDirectory(destination);
        var failed = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.False(failed.Established);
        Assert.Equal("failed", Assert.Single(product.Snapshot().Stages).Export?.Status);
        Assert.Equal(assessmentId, Assert.Single(product.Snapshot().Stages).Assessment?.Id);
        Assert.Empty(Directory.GetFiles(Path.GetDirectoryName(destination)!, "*.tmp"));
        Assert.Equal(ExportDeliveryStanding.Absent,
            product.VerifiedExportContent(stage.StageId, new string('0', 64)).Standing);
        Directory.Delete(destination);
        var retry = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(retry.Established, retry.Reason);
        Assert.Equal(assessmentId, Assert.Single(retry.Value!.Stages).Export?.AssessmentId);
        Assert.Equal(2, worker.ExportRequests.Count);
    }

    [Fact]
    public async Task Current_not_qualified_assessment_exports_its_reason_and_limitations_as_given()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var assessed = stage.Assessment! with
        {
            Id = "controlled-not-qualified-assessment",
            Qualification = PreparationQualification.NotQualified,
            Reason = "Controlled disqualifying scientific assessment premise",
            Limitations = ImmutableArray.Create("Controlled limitation remains disclosed")
        };
        var assessments = (Dictionary<string, PreparationAssessmentResult>)typeof(ProductRoot)
            .GetField("_stageAssessments", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        assessments[stage.StageId] = assessed;
        var workspace = (ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace)
            typeof(ProductRoot).GetField("_workspace", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        workspace.RetainAssessment(assessed);
        Assert.True(product.Snapshot().Actions.Single(action => action.Kind == ActorActionKind.ExportStage &&
            action.SubjectId == stage.StageId).Enabled);
        var delivered = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(delivered.Established, delivered.Reason);
        var resultingStage = Assert.Single(delivered.Value!.Stages);
        Assert.Equal(PreparationQualification.NotQualified, resultingStage.Assessment?.Qualification);
        Assert.Equal(assessed.Reason, resultingStage.Assessment?.Reason);
        var content = product.VerifiedExportContent(stage.StageId, resultingStage.Export?.Sha256);
        Assert.Equal(ExportDeliveryStanding.Available, content.Standing);
        using var archive = new ZipArchive(new MemoryStream(content.Bytes!), ZipArchiveMode.Read);
        using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
        var declared = manifest.RootElement.GetProperty("assessment");
        Assert.Equal("NotQualified", declared.GetProperty("qualification").GetString());
        Assert.Equal(assessed.Reason, declared.GetProperty("reason").GetString());
        Assert.Equal(assessed.Limitations[0],
            manifest.RootElement.GetProperty("limitations")[0].GetString());
    }

    [Fact]
    public async Task Reassessment_requires_a_new_export_identity_and_does_not_reuse_old_status()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var first = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(first.Established, first.Reason);
        var oldExport = Assert.Single(first.Value!.Stages).Export!;
        var oldBundle = product.VerifiedExportContent(stage.StageId, oldExport.Sha256);
        Assert.Equal(ExportDeliveryStanding.Available, oldBundle.Standing);

        var finding = new ScientificFinding("later-material-finding", stage.StageId,
            "later-evidence", "A later material observation challenges the stage.",
            "Reassess the selected completed stage.", FindingDisposition.Challenges,
            true, DateTimeOffset.UtcNow);
        var renewed = stage.Assessment! with
        {
            Id = "renewed-assessment",
            Qualification = PreparationQualification.NotQualified,
            Reason = "Current material finding disqualifies this stage.",
            Findings = stage.Assessment.Findings.Add(finding),
            Limitations = stage.Assessment.Limitations.Add("Later material finding applies.")
        };
        var assessments = (Dictionary<string, PreparationAssessmentResult>)typeof(ProductRoot)
            .GetField("_stageAssessments", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        assessments[stage.StageId] = renewed;
        var workspace = (ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace)
            typeof(ProductRoot).GetField("_workspace", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        workspace.RetainAssessment(renewed);
        Assert.Equal(ExportDeliveryStanding.Absent,
            product.VerifiedExportContent(stage.StageId, oldExport.Sha256).Standing);
        Assert.Null(Assert.Single(product.Snapshot().Stages).Export);

        var redelivered = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(redelivered.Established, redelivered.Reason);
        var current = Assert.Single(redelivered.Value!.Stages);
        Assert.Equal("renewed-assessment", current.Export?.AssessmentId);
        Assert.Equal(PreparationQualification.NotQualified, current.Assessment?.Qualification);
        Assert.NotEqual(oldExport.Sha256, current.Export?.Sha256);
        var content = product.VerifiedExportContent(stage.StageId, current.Export?.Sha256);
        using var archive = new ZipArchive(new MemoryStream(content.Bytes!), ZipArchiveMode.Read);
        using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
        Assert.Equal("renewed-assessment",
            manifest.RootElement.GetProperty("assessment").GetProperty("id").GetString());
        Assert.Contains("later-material-finding",
            manifest.RootElement.GetProperty("findings").EnumerateArray()
                .Select(item => item.GetProperty("id").GetString()));
        Assert.Equal(2, worker.ExportRequests.Count);
    }

    [Fact]
    public async Task Export_owner_refuses_each_identity_and_correspondence_violation_before_worker()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stageId = Assert.Single(product.Snapshot().Stages).StageId;
        var rootType = typeof(ProductRoot);
        T Field<T>(string name) => (T)rootType.GetField(name,
            System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic)!
            .GetValue(product)!;
        var stage = Field<Dictionary<string, CompletedStage>>("_stages")[stageId];
        var assessment = Field<Dictionary<string, PreparationAssessmentResult>>("_stageAssessments")[stageId];
        var revision = Field<Dictionary<string, StudyRevision>>("_revisions")[stage.Attempt.StudyRevisionId];
        var policy = Field<Dictionary<string, ApplicablePreparationPolicy>>("_policyByAttempt")[stage.Attempt.Id];
        var constructed = Field<Dictionary<string, ConstructedExplicitSystem>>("_constructedByAttempt")[stage.Attempt.Id];
        var decisions = Field<Dictionary<string, ImmutableArray<ResearcherDecision>>>("_decisionsByAttempt")
            [stage.Attempt.Id];
        var owner = new ExportOwner(worker);
        async Task Reject(string caseName, CompletedStage? selected = null,
            PreparationAssessmentResult? current = null, StudyRevision? origin = null,
            AssessedPreparedProtein? protein = null, ApplicablePreparationPolicy? selectedPolicy = null)
        {
            var before = worker.ExportRequests.Count;
            var outcome = await owner.ExportAsync(selected ?? stage, current ?? assessment,
                origin ?? revision, protein ?? fixture.Protein, fixture.Membrane,
                fixture.Placement, constructed, constructed.Derivation, selectedPolicy ?? policy,
                decisions, null, null, null,
                Path.Combine(fixture.Directory, "owner-refusal-" + caseName),
                TestContext.Current.CancellationToken);
            Assert.False(outcome.Established);
            Assert.Equal(before, worker.ExportRequests.Count);
        }

        await Reject("stale", current: assessment with { CurrentlyApplicable = false });
        await Reject("wrong-assessment-stage", current: assessment with { StageId = "other-stage" });
        await Reject("wrong-revision", origin: revision with { Id = "other-revision" });
        await Reject("wrong-attempt", selected: stage with
            { Attempt = stage.Attempt with { Id = "other-attempt" } });
        await Reject("wrong-protein", protein: fixture.Protein with { Id = "other-protein" });
        await Reject("wrong-policy", selectedPolicy: policy with { Id = "other-policy" });
        await Reject("incomplete", selected: stage with
            { Correspondence = stage.Correspondence with { Complete = false } });
        await Reject("duplicate-atom-id", selected: stage with
            { Correspondence = stage.Correspondence with { Atoms = stage.Correspondence.Atoms.SetItem(1,
                stage.Correspondence.Atoms[1] with
                { ResultAtomId = stage.Correspondence.Atoms[0].ResultAtomId }) } });
        await Reject("misordered-index", selected: stage with
            { Correspondence = stage.Correspondence with { Atoms = stage.Correspondence.Atoms.SetItem(1,
                stage.Correspondence.Atoms[1] with { ResultAtomIndex = 0 }) } });
        await Reject("unapproved-change", selected: stage with
            { Correspondence = stage.Correspondence with { Atoms = stage.Correspondence.Atoms.SetItem(0,
                stage.Correspondence.Atoms[0] with { ApprovedChangeId = "unapproved-change" }) } });
        await Reject("missing-topology", selected: stage with
            { Molecule = stage.Molecule with { TopologyPath = "not-present.json" } });
        await Reject("missing-system", selected: stage with
            { Molecule = stage.Molecule with { SystemXmlPath = "not-present.xml" } });
        await Reject("missing-state", selected: stage with
            { Molecule = stage.Molecule with { StateXmlPath = "not-present.xml" } });
        await Reject("bad-tolerance", selectedPolicy: policy with
            { ExportCoordinateReadBackToleranceAngstrom = -1 });
        Assert.Empty(worker.ExportRequests);
        using var canceled = new CancellationTokenSource();
        canceled.Cancel();
        await Assert.ThrowsAnyAsync<OperationCanceledException>(() => owner.ExportAsync(stage,
            assessment, revision, fixture.Protein, fixture.Membrane, fixture.Placement,
            constructed, constructed.Derivation, policy, decisions, null, null, null,
            Path.Combine(fixture.Directory, "owner-canceled"), canceled.Token));
        Assert.Empty(worker.ExportRequests);
    }

    [Fact]
    public async Task Export_worker_mismatch_classes_each_refuse_delivery_without_changing_assessment()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        var assessmentId = stage.Assessment!.Id;
        var changes = new (string Name, Func<WorkerResult<ExportVerificationObservations>,
            WorkerResult<ExportVerificationObservations>> Change)[]
        {
            ("request", result => result with { RequestId = "wrong-request" }),
            ("revision", result => result with { StudyRevisionId = "wrong-revision" }),
            ("attempt", result => result with { AttemptId = "wrong-attempt" }),
            ("stage", result => result with { StageId = "wrong-stage" }),
            ("unobserved", result => result with { Standing = WorkerResultStanding.Unobserved }),
            ("stopped", result => result with { Standing = WorkerResultStanding.Stopped }),
            ("timeout", result => result with { Standing = WorkerResultStanding.Failed,
                FailureCode = "timeout", FailureMessage = "Bounded worker timeout." }),
            ("missing-mmcif", result => result with { Artifacts = ImmutableArray<WorkerArtifact>.Empty }),
            ("count", result => result with { Observations = result.Observations! with
                { ExportedAtomCount = result.Observations.SourceAtomCount - 1 } }),
            ("order", result => result with { Observations = result.Observations! with
                { AtomOrderMatched = false } }),
            ("bonds", result => result with { Observations = result.Observations! with
                { BondsMatched = false } }),
            ("cell", result => result with { Observations = result.Observations! with
                { CellMatched = false } }),
            ("coordinate", result => result with { Observations = result.Observations! with
                { CoordinateMaxDeviationAngstrom = fixture.Policy.ExportCoordinateReadBackToleranceAngstrom + 1 } }),
            ("warning", result => result with { Observations = result.Observations! with
                { Warnings = ImmutableArray.Create("unresolved read-back") } })
        };
        foreach (var (name, change) in changes)
        {
            worker.ChangeNextExportResult = change;
            var refused = await Command(product, ActorActionKind.ExportStage,
                new { stageId = stage.StageId });
            Assert.False(refused.Established, name);
            var after = Assert.Single(product.Snapshot().Stages);
            Assert.Equal("failed", after.Export?.Status);
            Assert.Equal("completed", after.Status);
            Assert.Equal(assessmentId, after.Assessment?.Id);
            Assert.Equal(ExportDeliveryStanding.Absent,
                product.VerifiedExportContent(stage.StageId, new string('0', 64)).Standing);
        }
        Assert.Equal(changes.Length, worker.ExportRequests.Count);
        var retry = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(retry.Established, retry.Reason);
        Assert.Equal("verified", Assert.Single(retry.Value!.Stages).Export?.Status);
    }

    [Fact]
    public async Task Root_stop_matches_only_the_active_attempt_and_keeps_constructed_without_a_stage()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        var running = product.Snapshot();
        var id = running.Attempt!.AttemptId;
        var second = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.False(second.Established);
        var wrong = await Command(product, ActorActionKind.StopAttempt,
            new { attemptId = "other" });
        Assert.False(wrong.Established);
        Assert.Equal("running", product.Snapshot().Attempt?.Status);

        var stop = await Command(product, ActorActionKind.StopAttempt, new { attemptId = id });
        Assert.True(stop.Established, stop.Reason);
        await Until(() => product.Snapshot().Attempt?.Status == "stopped");
        var stopped = product.Snapshot();
        Assert.Equal(id, stopped.Attempt?.AttemptId);
        Assert.NotNull(stopped.Attempt?.Constructed);
        Assert.Empty(stopped.Stages);
        Assert.Single(worker.MinimizationRequests);
        Assert.Empty(worker.ObservationRequests);
    }

    [Fact]
    public async Task Ready_candidate_can_be_declined_without_starting_minimization()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        var ready = product.Snapshot();
        var candidateId = ready.Attempt!.Constructed!.SubjectId;
        Assert.True(ready.Actions.Single(item => item.Kind == ActorActionKind.StopAttempt).Enabled);

        var declined = await Command(product, ActorActionKind.StopAttempt,
            new { attemptId = ready.Attempt.AttemptId });
        Assert.True(declined.Established, declined.Reason);
        var stopped = product.Snapshot();
        Assert.Equal("stopped", stopped.Attempt?.Status);
        Assert.Equal(candidateId, stopped.Attempt?.Constructed?.SubjectId);
        Assert.Empty(stopped.Stages);
        Assert.Empty(worker.MinimizationRequests);
        Assert.False(stopped.Actions.Single(item => item.Kind == ActorActionKind.ContinueMinimization).Enabled);
        Assert.True(stopped.Actions.Single(item => item.Kind == ActorActionKind.StartPreparation).Enabled);
        Assert.False((await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.Attempt.AttemptId, constructedSubjectId = candidateId })).Established);
        Assert.True((await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = candidateId })).Established);
    }

    [Fact]
    public async Task Stop_during_active_construction_cancels_only_the_exact_attempt()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture) { BlockConstruction = true };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await worker.ConstructionEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        var active = product.Snapshot();
        var id = active.Attempt!.AttemptId;
        Assert.Equal("running", active.Attempt.Status);
        Assert.True(active.Actions.Single(item => item.Kind == ActorActionKind.StopAttempt).Enabled);
        Assert.False((await Command(product, ActorActionKind.StopAttempt,
            new { attemptId = "another-attempt" })).Established);
        Assert.True((await Command(product, ActorActionKind.StopAttempt,
            new { attemptId = id })).Established);
        await Until(() => product.Snapshot().Attempt?.Status == "stopped");
        Assert.Null(product.Snapshot().Attempt?.Constructed);
        Assert.Empty(product.Snapshot().Stages);
        Assert.Empty(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Declared_construction_timeout_cancels_work_without_a_candidate_or_stage()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture) { BlockConstruction = true };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http, maximumConstructionSeconds: 1);

        var started = await Command(product, ActorActionKind.StartPreparation, new { });
        Assert.True(started.Established, started.Reason);
        await worker.ConstructionEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
            TestContext.Current.CancellationToken);
        await Until(() => product.Snapshot().Attempt?.Status == "failed");

        var terminal = product.Snapshot();
        Assert.True(worker.ConstructionCancellationObserved);
        Assert.Single(worker.ConstructionRequests);
        Assert.Equal("failed", terminal.Attempt?.Status);
        Assert.Contains("exceeded its declared execution bound", terminal.Attempt?.Message ?? "",
            StringComparison.Ordinal);
        Assert.Null(terminal.Attempt?.Constructed);
        Assert.Empty(terminal.Stages);
        Assert.Empty(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Continue_requires_exact_current_candidate_and_refuses_duplicate_or_stale_binding()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        var ready = product.Snapshot().Attempt!;
        Assert.False((await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = "other", constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        Assert.False((await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId, constructedSubjectId = "other" })).Established);
        Assert.Empty(worker.MinimizationRequests);

        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        Assert.False((await Continue(product)).Established);
        Assert.Single(worker.MinimizationRequests);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1 &&
            product.Snapshot().Actions.Single(item => item.Kind == ActorActionKind.StartPreparation).Enabled);
    }

    [Theory]
    [InlineData("coordinates")]
    [InlineData("topology")]
    [InlineData("system")]
    [InlineData("state")]
    [InlineData("correspondence")]
    public async Task Continue_refuses_changed_constructed_artifact_bytes(string artifact)
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        var ready = product.Snapshot().Attempt!;
        var candidates = (Dictionary<string, ConstructedExplicitSystem>)typeof(ProductRoot)
            .GetField("_constructedByAttempt", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        var molecule = candidates[ready.AttemptId].Molecule;
        var path = artifact switch
        {
            "coordinates" => molecule.CoordinatePath,
            "topology" => molecule.TopologyPath!,
            "system" => molecule.SystemXmlPath!,
            "state" => molecule.StateXmlPath!,
            "correspondence" => molecule.CorrespondencePath!,
            _ => throw new ArgumentOutOfRangeException(nameof(artifact))
        };
        File.AppendAllText(path, "changed after candidate review");

        var refused = await Continue(product);
        Assert.False(refused.Established);
        Assert.Contains("changed", refused.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Equal("readyForMinimization", product.Snapshot().Attempt?.Status);
        Assert.Empty(product.Snapshot().Stages);
        Assert.Empty(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Continue_refuses_a_candidate_after_its_study_revision_is_replaced()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        var ready = product.Snapshot().Attempt!;
        Set(product, "_study", fixture.Revision with { Id = "later-revision", Number = 3 });
        Assert.False(product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.ContinueMinimization).Enabled);
        Assert.False((await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = ready.AttemptId, constructedSubjectId = ready.Constructed!.SubjectId })).Established);
        Assert.Empty(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Failed_retry_retains_the_earlier_completed_stage_under_its_original_attempt()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        await Until(() => product.Snapshot().Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        var original = Assert.Single(product.Snapshot().Stages);
        var originalAssessmentId = original.Assessment!.Id;
        worker.FailNextConstruction = true;

        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => product.Snapshot().Attempt?.Status == "failed");
        var after = product.Snapshot();
        var retained = Assert.Single(after.Stages);
        Assert.NotEqual(original.AttemptId, after.Attempt?.AttemptId);
        Assert.Null(after.Attempt?.Constructed);
        Assert.Equal(original.StageId, retained.StageId);
        Assert.Equal(original.AttemptId, retained.AttemptId);
        Assert.Equal(originalAssessmentId, retained.Assessment?.Id);
        Assert.Equal(original.Constructed?.SubjectId, retained.Constructed?.SubjectId);
        Assert.Single(worker.MinimizationRequests);
    }

    [Fact]
    public async Task Later_actor_membrane_revision_retains_prior_stage_and_assessment_only_under_original_identity()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var before = product.Snapshot();
        var stage = Assert.Single(before.Stages);
        var attemptId = stage.AttemptId;
        var assessmentId = stage.Assessment!.Id;
        var originalExport = await Command(product, ActorActionKind.ExportStage,
            new { stageId = stage.StageId });
        Assert.True(originalExport.Established, originalExport.Reason);
        var bundleSha = Assert.Single(originalExport.Value!.Stages).Export!.Sha256;
        var originalBytes = product.VerifiedExportContent(stage.StageId, bundleSha).Bytes;
        Assert.NotNull(originalBytes);

        var proposal = await Command(product, ActorActionKind.ProposeMembrane, new
        {
            upper = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            lower = new[] { new { speciesId = "DMPC", fraction = 1.0 } },
            scientificPurpose = "A later membrane decision"
        });
        Assert.True(proposal.Established, proposal.Reason);
        var proposedId = proposal.Value!.Membrane!.ModelId;
        var changed = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = proposedId });
        Assert.True(changed.Established, changed.Reason);
        var after = product.Snapshot();
        Assert.NotEqual(before.Study!.Id, after.Study!.Id);
        var retained = Assert.Single(after.Stages);
        Assert.Equal(stage.StageId, retained.StageId);
        Assert.Equal(attemptId, retained.AttemptId);
        Assert.Equal(fixture.Revision.Id, retained.StudyRevisionId);
        Assert.Equal(assessmentId, retained.Assessment?.Id);
        Assert.Equal(PreparationQualification.Indeterminate, retained.Assessment?.Qualification);
        Assert.Equal(bundleSha, retained.Export?.Sha256);
        Assert.Equal(originalBytes, product.VerifiedExportContent(stage.StageId, bundleSha).Bytes);
        Assert.True(after.Actions.Single(action => action.Kind == ActorActionKind.ExportStage &&
            action.SubjectId == stage.StageId).Enabled);
        Assert.False(after.Actions.Single(action =>
            action.Kind == ActorActionKind.StartPreparation).Enabled);
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = stage.StageId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal(stage.StageId, inspected.Value!.Inspection?.SubjectId);
    }

    [Fact]
    public async Task Completed_stage_warning_reaches_preparation_assessment_without_erasing_the_stage()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture) { ContactWarning = "observed unresolved contact" };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var stage = Assert.Single(product.Snapshot().Stages);
        Assert.Equal("completed", stage.Status);
        Assert.Equal(PreparationQualification.Indeterminate, stage.Assessment?.Qualification);
        Assert.Contains(stage.Assessment!.Findings, finding =>
            finding.SubjectId == stage.StageId && finding.Material &&
            finding.Disposition == FindingDisposition.Challenges);
    }

    [Fact]
    public async Task An_explicitly_proposal_attributed_stage_finding_reassesses_placement_and_preserves_prior_account()
    {
        using var fixture = new ConstructionFixture();
        var worker = new AttemptRouteWorker(fixture);
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);
        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1);
        var originalAccount = Assert.Single(product.Snapshot().Stages);
        var originalAssessmentId = originalAccount.Assessment!.Id;
        Set(product, "_placementProposal", fixture.Placement.Proposal);
        Assert.Equal("supported", product.Snapshot().Placement?.Status);

        Set(product, "_placementMeasurement", new PlacementMeasurementReport(
            fixture.Placement.Proposal.Id, ImmutableArray<MeasuredValue>.Empty,
            ImmutableArray<PlacementResidueObservation>.Empty,
            ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<string>.Empty));
        var stages = (Dictionary<string, CompletedStage>)typeof(ProductRoot)
            .GetField("_stages", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        var constructed = (Dictionary<string, ConstructedExplicitSystem>)typeof(ProductRoot)
            .GetField("_constructedByAttempt", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        var finding = new ScientificFinding(Guid.NewGuid().ToString("N"),
            fixture.Placement.Proposal.Id, "controlled-exact-measurement",
            "An independently measured premise contradicts the original proposal.",
            "Reassess the placement and dependent preparation qualification.",
            FindingDisposition.Disqualifies, true, DateTimeOffset.UtcNow);
        var source = stages[originalAccount.StageId];
        var downstreamStageId = Guid.NewGuid().ToString("N");
        // This controlled handoff supplies a distinct downstream stage with an
        // independently attributed finding. Warning strings alone do not do so.
        var completed = source with
        {
            Id = downstreamStageId,
            Kind = StageKind.Equilibration,
            Molecule = source.Molecule with { Id = downstreamStageId },
            Observation = source.Observation with { StageId = downstreamStageId,
                Kind = StageKind.Equilibration, Termination = StageTermination.Completed },
            Correspondence = source.Correspondence with { ResultId = downstreamStageId },
            SourceStageId = source.Id,
            Findings = ImmutableArray.Create(finding)
        };
        var workspace = (ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace)
            typeof(ProductRoot).GetField("_workspace", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;
        var running = new StageExecutionState(source.Attempt.Id, downstreamStageId,
            StageKind.Equilibration, StageExecutionStanding.Running,
            "Controlled downstream stage is running.", 0.5, DateTimeOffset.UtcNow);
        workspace.RetainExecution(running);
        Set(product, "_execution", running);
        typeof(ProductRoot).GetMethod("RetainCompleted", System.Reflection.BindingFlags.Instance |
            System.Reflection.BindingFlags.NonPublic)!.Invoke(product,
            [completed, constructed[originalAccount.AttemptId], fixture.Protein,
                fixture.Membrane, fixture.Placement, fixture.Policy]);

        var revised = product.Snapshot();
        Assert.Equal("unsupported", revised.Placement?.Status);
        Assert.Equal(2, revised.Stages.Length);
        var originalNow = revised.Stages.Single(item => item.StageId == originalAccount.StageId);
        var downstream = revised.Stages.Single(item => item.StageId == downstreamStageId);
        Assert.Equal("completed", originalNow.Status);
        Assert.Equal("completed", downstream.Status);
        Assert.Equal(fixture.Revision.Id, downstream.StudyRevisionId);
        Assert.NotEqual(PreparationQualification.QualifiedPrepared, originalNow.Assessment?.Qualification);
        Assert.Contains(finding, downstream.Assessment!.Findings);
        Assert.NotEqual(originalAssessmentId, originalNow.Assessment?.Id);
        Assert.Equal(originalAssessmentId, originalAccount.Assessment.Id);
        Assert.DoesNotContain(finding, originalAccount.Assessment.Findings);
        Assert.False(revised.Actions.Single(action =>
            action.Kind == ActorActionKind.StartPreparation).Enabled);
    }

    [Fact]
    public async Task Failed_or_unobserved_minimization_retains_only_the_verified_construction()
    {
        foreach (var (workerStanding, expectedStatus) in new[]
                 {
                     (WorkerResultStanding.Failed, "failed"),
                     (WorkerResultStanding.Unobserved, "unobserved")
                 })
        {
            using var fixture = new ConstructionFixture();
            var worker = new AttemptRouteWorker(fixture)
            { MinimizationStanding = workerStanding };
            using var http = new HttpClient();
            var product = Product(fixture, worker, http);
            Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
            await WaitForReady(product);
            Assert.True((await Continue(product)).Established);
            await WaitForMinimization(product, worker);
            worker.ReleaseMinimization();
            await Until(() => product.Snapshot().Attempt?.Status == expectedStatus);
            var state = product.Snapshot();
            Assert.Equal(expectedStatus, state.Attempt?.Status);
            Assert.NotNull(state.Attempt?.Constructed);
            Assert.Empty(state.Stages);
            Assert.Empty(worker.ObservationRequests);
        }
    }

    [Fact]
    public async Task Correlated_worker_process_loss_after_construction_reports_unobserved_without_a_stage()
    {
        using var fixture = new ConstructionFixture();
        var packageRoot = Path.Combine(fixture.Directory, "process-loss-worker");
        var package = Path.Combine(packageRoot, "ProteinInMembraneSystem");
        Directory.CreateDirectory(package);
        File.WriteAllText(Path.Combine(package, "__init__.py"), "");
        File.WriteAllText(Path.Combine(package, "worker.py"),
            "import json, os, sys, time\n" +
            "request = json.loads(sys.stdin.readline())\n" +
            "print(json.dumps({'requestId': request['requestId'], 'kind': 'progress', " +
            "'payload': {'attemptId': request['payload']['attemptId'], " +
            "'stageId': request['payload']['stageId'], 'pid': os.getpid()}}), flush=True)\n" +
            "time.sleep(60)\n");
        var python = Environment.GetEnvironmentVariable("PATH")!.Split(Path.PathSeparator)
            .Select(directory => Path.Combine(directory, "python3"))
            .First(File.Exists);
        var exchange = new ScientificWorkerExchange(python, packageRoot);
        var progress = new TaskCompletionSource<(string RequestId, string AttemptId, string StageId, int Pid)>(
            TaskCreationOptions.RunContinuationsAsynchronously);
        exchange.ProgressObserved += (requestId, payload) =>
        {
            var observed = (requestId,
                payload.GetProperty("attemptId").GetString()!,
                payload.GetProperty("stageId").GetString()!,
                payload.GetProperty("pid").GetInt32());
            progress.TrySetResult(observed);
            Process.GetProcessById(observed.Item4).Kill();
        };
        var worker = new AttemptRouteWorker(fixture) { MinimizationExchange = exchange };
        using var http = new HttpClient();
        var product = Product(fixture, worker, http);

        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await WaitForReady(product);
        Assert.True((await Continue(product)).Established);
        await WaitForMinimization(product, worker);
        var admitted = product.Snapshot();
        var attemptId = admitted.Attempt!.AttemptId;
        var constructedId = admitted.Attempt.Constructed!.SubjectId;
        Assert.Empty(admitted.Stages);

        var eventIdentity = await progress.Task.WaitAsync(TimeSpan.FromSeconds(10),
            TestContext.Current.CancellationToken);
        Assert.Equal(attemptId, eventIdentity.AttemptId);
        Assert.Equal(worker.MinimizationRequests.Single().RequestId, eventIdentity.RequestId);
        Assert.Equal(worker.MinimizationRequests.Single().Payload.StageId, eventIdentity.StageId);
        await Until(() => product.Snapshot().Attempt?.Status == "unobserved");
        var terminal = product.Snapshot();
        Assert.Equal(attemptId, terminal.Attempt?.AttemptId);
        Assert.Equal(constructedId, terminal.Attempt?.Constructed?.SubjectId);
        Assert.Contains("worker", terminal.Attempt?.Message ?? "", StringComparison.OrdinalIgnoreCase);
        Assert.Empty(terminal.Stages);
        Assert.Empty(worker.ObservationRequests);
    }

    private static ProductRoot Product(ConstructionFixture fixture, AttemptRouteWorker worker,
        HttpClient http, Func<ConstructionProviderInstallation?>? providerProbe = null,
        bool expectStartAvailable = true, int? maximumConstructionSeconds = null)
    {
        var policy = fixture.Policy with
        {
            Water = fixture.Policy.Water with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Sodium = fixture.Policy.Sodium with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Chloride = fixture.Policy.Chloride with { StereoChecks = ImmutableArray<MolecularStereoCheck>.Empty },
            Construction = fixture.Policy.Construction with
            { MaximumConstructionSeconds = maximumConstructionSeconds ??
                fixture.Policy.Construction.MaximumConstructionSeconds }
        };
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        var catalogue = Path.Combine(fixture.Directory, "controlled-catalogue.json");
        File.WriteAllText(catalogue, JsonSerializer.Serialize(new
        {
            version = "controlled-slice4",
            evidenceReferences = new[] { "controlled C# composition proof" },
            lipids = Array.Empty<object>(), proteinChemicalStates = Array.Empty<object>(),
            proteinStructuralPolicies = Array.Empty<object>(), membranePolicies = Array.Empty<object>(),
            placementPolicies = Array.Empty<object>(), placementWitnesses = Array.Empty<object>(),
            preparationPolicies = new[] { policy }, equilibrationQualifications = Array.Empty<object>(),
            ppmVersion = "", ppmExecutableSha256 = "", maximumSourceAtoms = 100000
        }, options));
        providerProbe ??= () => new ConstructionProviderInstallation(
            fixture.Policy.Construction.ProviderVersion, fixture.NativePatchPath, fixture.NativePatchSha);
        var root = new ProductRoot(worker, new ExternalSourceExchange(http),
            Path.Combine(fixture.Directory, "root-workspace"), catalogue, "", providerProbe);
        // The root has already retained its initial revision 1. This controlled
        // revision is a successor with the same scientific identity as the fixture.
        var revision = fixture.Revision with { Number = 2 };
        Set(root, "_study", revision);
        Set(root, "_protein", fixture.Protein);
        Set(root, "_membrane", fixture.Membrane);
        Set(root, "_placement", fixture.Placement);
        var revisions = (Dictionary<string, StudyRevision>)typeof(ProductRoot)
            .GetField("_revisions", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(root)!;
        revisions.Add(revision.Id, revision);
        var workspace = (ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace)
            typeof(ProductRoot).GetField("_workspace", System.Reflection.BindingFlags.Instance |
                System.Reflection.BindingFlags.NonPublic)!.GetValue(root)!;
        workspace.RetainStudy(revision);
        var snapshot = root.Snapshot();
        Assert.Equal(expectStartAvailable, snapshot.Actions.Single(item =>
            item.Kind == ActorActionKind.StartPreparation).Enabled);
        return root;
    }

    private static void Set(ProductRoot product, string name, object value) =>
        typeof(ProductRoot).GetField(name, System.Reflection.BindingFlags.Instance |
            System.Reflection.BindingFlags.NonPublic)!.SetValue(product, value);

    private static Task<BoundaryOutcome<WorkspaceState>> Command(ProductRoot product,
        ActorActionKind kind, object data) => product.ExecuteAsync(new ActorCommand(kind,
            JsonSerializer.SerializeToElement(data), product.Snapshot().Revision),
            TestContext.Current.CancellationToken);

    private static Task<BoundaryOutcome<WorkspaceState>> Continue(ProductRoot product)
    {
        var attempt = Assert.IsType<AttemptAccount>(product.Snapshot().Attempt);
        var candidate = Assert.IsType<ConstructedSystemAccount>(attempt.Constructed);
        return Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = attempt.AttemptId, constructedSubjectId = candidate.SubjectId });
    }

    private static async Task Until(Func<bool> done)
    {
        var deadline = DateTime.UtcNow + TimeSpan.FromSeconds(10);
        while (!done() && DateTime.UtcNow < deadline)
            await Task.Delay(20, TestContext.Current.CancellationToken);
        Assert.True(done(), "The corresponding root state did not reach its expected terminal account.");
    }

    private static async Task WaitForMinimization(ProductRoot product, AttemptRouteWorker worker)
    {
        try
        {
            await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(3),
                TestContext.Current.CancellationToken);
        }
        catch (TimeoutException)
        {
            var state = product.Snapshot();
            throw new TimeoutException($"Minimization was not entered: status={state.Attempt?.Status}, " +
                $"message={state.Attempt?.Message}, constructionRequests={worker.ConstructionRequests.Count}");
        }
    }

    private static Task WaitForReady(ProductRoot product) => Until(() =>
        product.Snapshot().Attempt is { Status: "readyForMinimization", Constructed: not null } &&
        product.Snapshot().Actions.Single(item => item.Kind == ActorActionKind.ContinueMinimization).Enabled);
}

internal sealed class AttemptRouteWorker(ConstructionFixture fixture) : IScientificWorkerExchange
{
    private readonly ConstructionWorker _construction = new(fixture);
    private readonly TaskCompletionSource _releaseConstruction = new(TaskCreationOptions.RunContinuationsAsynchronously);
    private readonly TaskCompletionSource _releaseMinimization = new(TaskCreationOptions.RunContinuationsAsynchronously);
    public TaskCompletionSource ConstructionEntered { get; } =
        new(TaskCreationOptions.RunContinuationsAsynchronously);
    public TaskCompletionSource MinimizationEntered { get; } =
        new(TaskCreationOptions.RunContinuationsAsynchronously);
    public List<ScientificWorkRequest<ConstructionPayload>> ConstructionRequests =>
        _construction.ConstructionRequests;
    public List<ScientificWorkRequest<MinimizationPayload>> MinimizationRequests { get; } = [];
    public List<ScientificWorkRequest<StageObservationPayload>> ObservationRequests { get; } = [];
    public List<ScientificWorkRequest<ExportVerificationPayload>> ExportRequests { get; } = [];
    public bool FailNextExport { get; set; }
    public bool ChangeNextExportArtifactAfterHash { get; set; }
    public Func<WorkerResult<ExportVerificationObservations>,
        WorkerResult<ExportVerificationObservations>>? ChangeNextExportResult { get; set; }
    public bool FailNextConstruction { get; set; }
    public bool BlockConstruction { get; set; }
    public bool ConstructionCancellationObserved { get; private set; }
    public int? ConstructedAtomCount { get; private set; }
    public string? ContactWarning { get; init; }
    public WorkerResultStanding MinimizationStanding { get; init; } = WorkerResultStanding.Observed;
    public ScientificWorkerExchange? MinimizationExchange { get; init; }
    public void ReleaseConstruction() => _releaseConstruction.TrySetResult();
    public void ReleaseMinimization() => _releaseMinimization.TrySetResult();

    public async Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
        ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken)
    {
        WorkerResult<ConstructionObservations> result;
        if (!FailNextConstruction)
            result = await _construction.ConstructSystemAsync(request, cancellationToken);
        else
        {
            FailNextConstruction = false;
            result = new WorkerResult<ConstructionObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, null,
            WorkerResultStanding.Failed, ImmutableArray<WorkerArtifact>.Empty, null,
            new ProviderIdentity("controlled native OpenMM", "8.6"), "constructionFailed",
            "The controlled retry failed before construction.");
        }
        ConstructedAtomCount = result.Observations?.AtomCount;
        ConstructionEntered.TrySetResult();
        if (BlockConstruction)
        {
            try { await _releaseConstruction.Task.WaitAsync(cancellationToken); }
            catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
            {
                ConstructionCancellationObserved = true;
                throw;
            }
        }
        return result;
    }

    public async Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
        ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken)
    {
        MinimizationRequests.Add(request);
        MinimizationEntered.TrySetResult();
        if (MinimizationExchange is not null)
            return await MinimizationExchange.MinimizeAsync(request, cancellationToken);
        await _releaseMinimization.Task.WaitAsync(cancellationToken);
        if (MinimizationStanding != WorkerResultStanding.Observed)
            return new WorkerResult<MinimizationObservations>(request.RequestId,
                request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
                MinimizationStanding, ImmutableArray<WorkerArtifact>.Empty, null,
                new ProviderIdentity("controlled OpenMM", "8.6"), "noFinalState",
                "No completed minimized state was observed.");
        var coordinate = Path.Combine(request.WorkingDirectory, "minimized.cif");
        var state = Path.Combine(request.WorkingDirectory, "minimized-state.xml");
        File.WriteAllText(coordinate, "controlled minimized coordinates");
        File.WriteAllText(state, "controlled minimized state");
        var atomCount = ConstructedAtomCount ?? throw new InvalidOperationException("No native candidate atom count was observed.");
        return new WorkerResult<MinimizationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed, ImmutableArray.Create(
                Artifact(coordinate, "minimizedCif"), Artifact(state, "minimizedStateXml")),
            new MinimizationObservations(-100, -110, 8, 19, StageTermination.Converged,
                true, atomCount, ImmutableArray<string>.Empty,
                FinalRawRmsForceKjMolNm: 500,
                MaximumRelativeConstraintError: 0.0000001,
                AppliedConstraintTolerance: 0.00001,
                FinalRmsForceMethod: "constraint-tangent-per-particle-v1"),
            new ProviderIdentity("controlled OpenMM", "8.6"), null, null);
    }

    public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken)
    {
        ObservationRequests.Add(request);
        var atomCount = ConstructedAtomCount ?? throw new InvalidOperationException("No native candidate atom count was observed.");
        var local = new LocalStateObservations(ObservationStanding.Unavailable,
            "No positive qualification observation in this controlled composition proof",
            ImmutableArray<MeasuredValue>.Empty, ImmutableArray<LocalContactObservation>.Empty,
            ImmutableArray<LocalContactRolePair>.Empty, ImmutableArray<string>.Empty,
            ImmutableArray<LocalRolePairMeasurement>.Empty);
        var geometry = new ProteinGeometryObservations(ObservationStanding.Unavailable,
            ImmutableArray<ProteinGeometryKindObservation>.Empty,
            ImmutableArray<GeometryDistanceObservation>.Empty, ImmutableArray<string>.Empty);
        return Task.FromResult(new WorkerResult<StageObservationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed, ImmutableArray<WorkerArtifact>.Empty,
            new StageObservationObservations(atomCount, true, true,
                ImmutableArray<MeasuredValue>.Empty, ImmutableArray<string>.Empty,
                ContactWarning is null ? ImmutableArray<string>.Empty :
                    ImmutableArray.Create(ContactWarning), ImmutableArray<string>.Empty, local, geometry),
            new ProviderIdentity("controlled observation", "1"), null, null));
    }

    private static WorkerArtifact Artifact(string path, string role) => new(role, path,
        Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant());

    public Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
        ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Source selection is pre-established in this composition fixture.");
    public Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
        ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Source selection is pre-established in this composition fixture.");
    public Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
        ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Protein is pre-established in this composition fixture.");
    public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
        ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Membrane is pre-established in this composition fixture.");
    public Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
        ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Placement is pre-established in this composition fixture.");
    public Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
        ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Placement is pre-established in this composition fixture.");
    public Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
        ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Placement is pre-established in this composition fixture.");
    public Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
        ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Placement is pre-established in this composition fixture.");
    public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
        ScientificWorkRequest<EquilibrationPayload> request,
        IProgress<EquilibrationWorkProgress>? progress, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Optional equilibration is outside Slice 4.");
    public Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
        ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken) =>
        Task.FromResult(VerifyExport(request));

    private WorkerResult<ExportVerificationObservations> VerifyExport(
        ScientificWorkRequest<ExportVerificationPayload> request)
    {
        ExportRequests.Add(request);
        if (FailNextExport)
        {
            FailNextExport = false;
            return new WorkerResult<ExportVerificationObservations>(request.RequestId,
                request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
                WorkerResultStanding.Failed, ImmutableArray<WorkerArtifact>.Empty, null,
                new ProviderIdentity("controlled export read-back", "1"), "controlledFailure",
                "The controlled export read-back failed.");
        }
        Directory.CreateDirectory(request.WorkingDirectory);
        var path = Path.Combine(request.WorkingDirectory, "prepared-stage.cif");
        File.Copy(request.Payload.TopologyCifPath, path, overwrite: true);
        var artifact = Artifact(path, "stageMmcif");
        if (ChangeNextExportArtifactAfterHash)
        {
            ChangeNextExportArtifactAfterHash = false;
            File.AppendAllText(path, "changed after worker hashing");
        }
        var atoms = ConstructedAtomCount ?? 0;
        var observed = new WorkerResult<ExportVerificationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed, ImmutableArray.Create(artifact),
            new ExportVerificationObservations(atoms, atoms, atoms,
                true, true, true, true, 0, ImmutableArray<string>.Empty),
            new ProviderIdentity("controlled export read-back", "1"), null, null);
        var change = ChangeNextExportResult;
        ChangeNextExportResult = null;
        return change is null ? observed : change(observed);
    }
}
