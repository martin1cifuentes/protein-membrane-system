using System.Collections.Immutable;
using System.Text.Json;
using System.Text.Json.Serialization;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class UnavailableOptionalRouteTests
{
    [Fact]
    public async Task No_applicable_optional_policy_refuses_request_and_keeps_minimized_assessment_and_export()
    {
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        var worker = new AttemptRouteWorker(fixture);
        var product = Product(fixture, worker, http);
        Assert.Null(fixture.Policy.OptionalEquilibration);

        Assert.True((await Command(product, ActorActionKind.StartPreparation, new { })).Established);
        await Until(() => product.Snapshot().Attempt?.Status == "readyForMinimization");
        var candidate = product.Snapshot().Attempt!;
        Assert.True((await Command(product, ActorActionKind.ContinueMinimization,
            new { attemptId = candidate.AttemptId,
                constructedSubjectId = candidate.Constructed!.SubjectId })).Established);
        await worker.MinimizationEntered.Task.WaitAsync(TimeSpan.FromSeconds(10),
            TestContext.Current.CancellationToken);
        worker.ReleaseMinimization();
        await Until(() => product.Snapshot().Stages.Length == 1 &&
            product.Snapshot().Actions.Single(action =>
                action.Kind == ActorActionKind.StartPreparation).Enabled);

        var before = Assert.Single(product.Snapshot().Stages);
        var assessment = Assert.IsType<PreparationAssessmentResult>(before.Assessment);
        Assert.Equal(StageKind.Minimization, before.Kind);
        Assert.Equal("completed", before.Status);
        Assert.True(assessment.CurrentlyApplicable);
        Assert.Equal(PreparationQualification.Indeterminate, assessment.Qualification);
        var optionalAction = product.Snapshot().Actions.Single(action =>
            action.Kind == ActorActionKind.RequestEquilibration &&
            action.SubjectId == before.StageId);
        Assert.False(optionalAction.Enabled);
        Assert.Contains("No validated applicable optional procedure", optionalAction.Reason ?? "",
            StringComparison.Ordinal);

        var refused = await Command(product, ActorActionKind.RequestEquilibration,
            new { stageId = before.StageId });
        Assert.False(refused.Established);
        Assert.Contains("No validated, applicable optional equilibration procedure",
            refused.Reason ?? "", StringComparison.Ordinal);
        var afterRefusal = Assert.Single(product.Snapshot().Stages);
        Assert.Equal(before.StageId, afterRefusal.StageId);
        Assert.Equal(before.AttemptId, afterRefusal.AttemptId);
        Assert.Equal("completed", afterRefusal.Status);
        Assert.Equal(assessment, afterRefusal.Assessment);
        Assert.Empty(worker.ExportRequests);

        var selected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = before.StageId });
        Assert.True(selected.Established, selected.Reason);
        Assert.Equal(before.StageId, selected.Value!.Inspection?.SubjectId);
        var delivered = await Command(product, ActorActionKind.ExportStage,
            new { stageId = before.StageId });
        Assert.True(delivered.Established, delivered.Reason);
        var exported = Assert.Single(delivered.Value!.Stages);
        var export = Assert.IsType<StageExportAccount>(exported.Export);
        Assert.Equal("verified", export.Status);
        Assert.Equal(before.StageId, export.StageId);
        Assert.Equal(assessment.Id, export.AssessmentId);
        Assert.Equal(assessment, exported.Assessment);
        Assert.Single(worker.ExportRequests);
        Assert.Equal(before.StageId, worker.ExportRequests[0].Payload.StageId);
        var content = product.VerifiedExportContent(before.StageId, export.Sha256);
        Assert.Equal(ExportDeliveryStanding.Available, content.Standing);
        Assert.NotEmpty(content.Bytes!);

        var refusedAgain = await Command(product, ActorActionKind.RequestEquilibration,
            new { stageId = before.StageId });
        Assert.False(refusedAgain.Established);
        var retained = Assert.Single(product.Snapshot().Stages);
        Assert.Equal(before.StageId, retained.StageId);
        Assert.Equal(assessment, retained.Assessment);
        Assert.Equal("verified", retained.Export?.Status);
        Assert.Equal(export.Sha256, retained.Export?.Sha256);
        Assert.Equal(content.Bytes, product.VerifiedExportContent(before.StageId, export.Sha256).Bytes);
        Assert.Single(worker.ExportRequests);
    }

    private static ProductRoot Product(ConstructionFixture fixture, AttemptRouteWorker worker,
        HttpClient http)
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
        var catalogue = Path.Combine(fixture.Directory, "controlled-unavailable-optional-catalogue.json");
        File.WriteAllText(catalogue, JsonSerializer.Serialize(new
        {
            version = "controlled-slice7-unavailable",
            evidenceReferences = new[] { "controlled C# composition proof" },
            lipids = Array.Empty<object>(), proteinChemicalStates = Array.Empty<object>(),
            proteinStructuralPolicies = Array.Empty<object>(), membranePolicies = Array.Empty<object>(),
            placementPolicies = Array.Empty<object>(), placementWitnesses = Array.Empty<object>(),
            preparationPolicies = new[] { policy }, equilibrationQualifications = Array.Empty<object>(),
            ppmVersion = "", ppmExecutableSha256 = "", maximumSourceAtoms = 100000
        }, options));
        var provider = new ConstructionProviderInstallation(
            policy.Construction.ProviderVersion, fixture.NativePatchPath, fixture.NativePatchSha);
        var product = new ProductRoot(worker, new ExternalSourceExchange(http),
            Path.Combine(fixture.Directory, "root-workspace"), catalogue, "", () => provider);
        var revision = fixture.Revision with { Number = 2 };
        Set(product, "_study", revision);
        Set(product, "_protein", fixture.Protein);
        Set(product, "_membrane", fixture.Membrane);
        Set(product, "_placement", fixture.Placement);
        Field<Dictionary<string, StudyRevision>>(product, "_revisions").Add(revision.Id, revision);
        Field<LocalRunWorkspace>(product, "_workspace").RetainStudy(revision);
        Assert.True(product.Snapshot().Actions.Single(action =>
            action.Kind == ActorActionKind.StartPreparation).Enabled);
        return product;
    }

    private static T Field<T>(ProductRoot product, string name) =>
        (T)typeof(ProductRoot).GetField(name,
            System.Reflection.BindingFlags.Instance |
            System.Reflection.BindingFlags.NonPublic)!.GetValue(product)!;

    private static void Set(ProductRoot product, string name, object value) =>
        typeof(ProductRoot).GetField(name,
            System.Reflection.BindingFlags.Instance |
            System.Reflection.BindingFlags.NonPublic)!.SetValue(product, value);

    private static Task<BoundaryOutcome<WorkspaceState>> Command(ProductRoot product,
        ActorActionKind kind, object data) => product.ExecuteAsync(new ActorCommand(kind,
            JsonSerializer.SerializeToElement(data), product.Snapshot().Revision),
            TestContext.Current.CancellationToken);

    private static async Task Until(Func<bool> done)
    {
        var deadline = DateTime.UtcNow + TimeSpan.FromSeconds(30);
        while (!done() && DateTime.UtcNow < deadline)
            await Task.Delay(20, TestContext.Current.CancellationToken);
        Assert.True(done(), "The controlled route did not reach its expected account state.");
    }
}
