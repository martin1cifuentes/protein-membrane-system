using System.Collections.Immutable;
using System.Net;
using System.Security.Cryptography;
using System.Text.Json;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class MembraneModelRouteTests
{
    [Fact]
    public async Task A_proposal_becomes_a_chosen_and_assessed_model_only_after_explicit_adoption()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker();
        var product = Product(directory.Path, worker);
        var initial = product.Snapshot();

        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 1.0), ("POPC", 1.0)));

        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal(initial.Study!.Id, proposed.Value!.Study!.Id);
        Assert.Equal("proposed", proposed.Value.Membrane!.Status);
        Assert.Empty(worker.MembraneRequests);
        Assert.Equal(0.15, proposed.Value.Study.Conditions.TargetNaClMolar);
        Assert.Equal(303, proposed.Value.Study.Conditions.OptionalTemperatureKelvin);

        var adopted = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value.Membrane.ModelId });

        Assert.True(adopted.Established, adopted.Reason);
        Assert.Equal("assessed", adopted.Value!.Membrane!.Status);
        Assert.NotEqual(initial.Study.Id, adopted.Value.Study!.Id);
        Assert.Equal(initial.Study.Number + 1, adopted.Value.Study.Number);
        Assert.Single(worker.MembraneRequests);
        Assert.Equal(adopted.Value.Study.Id, worker.MembraneRequests[0].Payload.StudyRevisionId);
        Assert.Equal(adopted.Value.Membrane.ModelId, worker.MembraneRequests[0].Payload.MembraneModelId);
        Assert.Equal("POPC", adopted.Value.Membrane.Upper.Single().SpeciesId);
        Assert.Equal("POPC", adopted.Value.Membrane.Lower.Single().SpeciesId);
        Assert.Equal("policy-popc", adopted.Value.Membrane.PolicyId);
        Assert.Equal("1.0", adopted.Value.Membrane.PolicyVersion);
        var wire = JsonSerializer.SerializeToElement(adopted.Value,
            new JsonSerializerOptions(JsonSerializerDefaults.Web));
        Assert.Equal("policy-popc", wire.GetProperty("membrane").GetProperty("policyId").GetString());
        Assert.Equal("1.0", wire.GetProperty("membrane").GetProperty("policyVersion").GetString());
        Assert.Null(adopted.Value.Placement);
        Assert.Empty(adopted.Value.Stages);
    }

    [Fact]
    public async Task Zero_fraction_species_does_not_remove_qualified_membrane_support()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker();
        var product = Product(directory.Path, worker);
        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 1.0), ("UNSELECTED", 0.0), ("POPC", 1.0)));
        Assert.True(proposed.Established, proposed.Reason);

        var adopted = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value!.Membrane!.ModelId });

        Assert.True(adopted.Established, adopted.Reason);
        Assert.Equal("assessed", adopted.Value!.Membrane!.Status);
        Assert.Single(worker.MembraneRequests);
        Assert.Equal(["POPC"], worker.MembraneRequests[0].Payload.SpeciesRepresentations.Select(item => item.SpeciesId));
    }

    [Fact]
    public async Task Unsupported_chosen_model_remains_identified_and_inspectable_without_substitution()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker { CombinedParameterizationObserved = false };
        var product = Product(directory.Path, worker);
        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 1.0), ("POPC", 1.0)));
        Assert.True(proposed.Established, proposed.Reason);
        var modelId = proposed.Value!.Membrane!.ModelId;

        var adopted = await Command(product, ActorActionKind.AdoptMembrane, new { modelId });
        Assert.True(adopted.Established, adopted.Reason);
        Assert.Equal("notEstablished", adopted.Value!.Membrane!.Status);
        Assert.Equal(modelId, adopted.Value.Membrane.ModelId);
        Assert.Contains(adopted.Value.Notices, notice => notice.SubjectId == modelId &&
            notice.Message.Contains("combination", StringComparison.OrdinalIgnoreCase) &&
            notice.Message.Contains("parameterization", StringComparison.OrdinalIgnoreCase));
        Assert.Null(adopted.Value.Placement);

        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject, new { subjectId = modelId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal(modelId, inspected.Value!.Inspection!.SubjectId);
        Assert.True(inspected.Value.Inspection.Annotations.Length > 0 ||
            inspected.Value.Inspection.Metrics.Length > 0,
            "The selected unsupported membrane needs its own connected explanation.");
    }

    [Fact]
    public async Task Current_assessed_model_exposes_membrane_local_evidence_in_connected_inspection()
    {
        using var directory = new TemporaryDirectory();
        var product = Product(directory.Path, new MembraneRouteWorker());
        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 1.0), ("POPC", 1.0)));
        var modelId = proposed.Value!.Membrane!.ModelId;
        var adopted = await Command(product, ActorActionKind.AdoptMembrane, new { modelId });
        Assert.Equal("assessed", adopted.Value!.Membrane!.Status);

        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject, new { subjectId = modelId });

        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal("membraneModel", inspected.Value!.Inspection!.RepresentationKind);
        Assert.Contains(inspected.Value.Inspection.Metrics,
            metric => metric.Value.Contains("POPC", StringComparison.Ordinal) ||
                metric.Name.Contains("representation", StringComparison.OrdinalIgnoreCase));
        Assert.Equal("POPC", inspected.Value.Membrane!.Upper.Single().SpeciesId);
        Assert.Equal("POPC", inspected.Value.Membrane.Lower.Single().SpeciesId);
        Assert.Null(inspected.Value.Placement);
    }

    [Fact]
    public async Task Consequential_membrane_change_advances_revision_without_transferring_the_old_assessment()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker();
        var product = Product(directory.Path, worker);
        var firstProposal = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 1.0), ("POPC", 1.0)));
        var firstModelId = firstProposal.Value!.Membrane!.ModelId;
        var original = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = firstModelId });
        Assert.Equal("assessed", original.Value!.Membrane!.Status);
        Assert.Single(worker.MembraneRequests);

        var changedProposal = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("EXACT-OTHER", 1.0), ("EXACT-OTHER", 1.0)));
        Assert.True(changedProposal.Established, changedProposal.Reason);
        Assert.Equal(original.Value.Study!.Id, changedProposal.Value!.Study!.Id);
        Assert.Equal("proposed", changedProposal.Value.Membrane!.Status);
        Assert.NotEqual(firstModelId, changedProposal.Value.Membrane.ModelId);

        var changed = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = changedProposal.Value.Membrane.ModelId });
        Assert.True(changed.Established, changed.Reason);
        Assert.NotEqual(original.Value.Study.Id, changed.Value!.Study!.Id);
        Assert.Equal(original.Value.Study.Number + 1, changed.Value.Study.Number);
        Assert.Equal("notEstablished", changed.Value.Membrane!.Status);
        Assert.Equal("EXACT-OTHER", changed.Value.Membrane.Upper.Single().SpeciesId);
        Assert.Equal("assessed", original.Value.Membrane.Status);
        Assert.Equal(firstModelId, original.Value.Membrane.ModelId);
        Assert.Single(worker.MembraneRequests);
        Assert.Null(changed.Value.Placement);
    }

    [Fact]
    public async Task Incoherent_proposal_does_not_change_study_or_chosen_membrane()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker();
        var product = Product(directory.Path, worker);
        var initial = product.Snapshot();

        var refused = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("POPC", 0.4), ("POPC", 1.0)));

        Assert.False(refused.Established);
        Assert.Contains("coherent", refused.Reason);
        Assert.Equal(initial.Study!.Id, product.Snapshot().Study!.Id);
        Assert.Null(product.Snapshot().Membrane);
        Assert.Empty(worker.MembraneRequests);
    }

    [Fact]
    public async Task Researcher_can_choose_an_exact_unsupported_species_and_receive_non_establishment()
    {
        using var directory = new TemporaryDirectory();
        var worker = new MembraneRouteWorker();
        var product = Product(directory.Path, worker);
        Assert.True(product.Snapshot().Actions.Single(item => item.Kind == ActorActionKind.ProposeMembrane).Enabled);

        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            Composition(("EXACT-OTHER", 1.0), ("EXACT-OTHER", 1.0)));
        Assert.True(proposed.Established, proposed.Reason);
        var modelId = proposed.Value!.Membrane!.ModelId;
        Assert.Equal("proposed", proposed.Value.Membrane.Status);

        var adopted = await Command(product, ActorActionKind.AdoptMembrane, new { modelId });

        Assert.True(adopted.Established, adopted.Reason);
        Assert.Equal("notEstablished", adopted.Value!.Membrane!.Status);
        Assert.Equal("EXACT-OTHER", adopted.Value.Membrane.Upper.Single().SpeciesId);
        Assert.Contains(adopted.Value.Notices, notice => notice.SubjectId == modelId &&
            notice.Message.Contains("EXACT-OTHER", StringComparison.Ordinal) &&
            notice.Message.Contains("representation", StringComparison.OrdinalIgnoreCase));
        Assert.Empty(worker.MembraneRequests);
        Assert.Null(adopted.Value.Placement);
    }

    private static object Composition((string Species, double Fraction) firstUpper,
        (string Species, double Fraction) secondUpper,
        (string Species, double Fraction) lower)
        => new { upper = new[] { Fraction(firstUpper), Fraction(secondUpper) },
            lower = new[] { Fraction(lower) }, scientificPurpose = Purpose };

    private static object Composition((string Species, double Fraction) upper,
        (string Species, double Fraction) lower)
        => new { upper = new[] { Fraction(upper) }, lower = new[] { Fraction(lower) },
            scientificPurpose = Purpose };

    private const string Purpose = "A deliberately simplified POPC reference membrane.";
    private static object Fraction((string Species, double Fraction) item) =>
        new { speciesId = item.Species, fraction = item.Fraction };

    private static ProductRoot Product(string directory, MembraneRouteWorker worker)
    {
        var policy = WriteCatalogue(directory);
        var http = new HttpClient(new NoNetworkHandler());
        return new ProductRoot(worker, new ExternalSourceExchange(http),
            Path.Combine(directory, "workspace"), policy, string.Empty, () => null);
    }

    private static string WriteCatalogue(string directory)
    {
        var forceFieldPath = Path.Combine(directory, "lipid21.xml");
        var coordinatePath = Path.Combine(directory, "POPC.pdb");
        File.WriteAllText(forceFieldPath, "<ForceField/>\n");
        File.WriteAllText(coordinatePath, "HETATM test-only POPC coordinate identity\n");
        var lipid = new MolecularRepresentation("POPC", "POPC-chemistry", "lipid",
            forceFieldPath, Hash(forceFieldPath), coordinatePath, Hash(coordinatePath),
            134, 0, 68, 1200, ImmutableArray.Create(1, 20), "Lipid21", "8.6.0",
            ImmutableArray.Create("Membrane-local representation only."),
            ImmutableArray.Create(new MolecularStereoCheck("tetrahedral",
                ImmutableArray.Create("O21", "C1", "C3", "HS"), "negative")));
        var policy = new MembraneSupportPolicy("policy-popc", "1.0", ImmutableArray.Create("POPC"),
            false, false, true, ImmutableArray.Create("identified local policy evidence"),
            ImmutableArray.Create("Whole-system compatibility remains separate."));
        var path = Path.Combine(directory, "catalogue.json");
        File.WriteAllText(path, JsonSerializer.Serialize(new
        {
            version = "membrane-test", evidenceReferences = new[] { "test catalogue identity" },
            lipids = new[] { lipid }, proteinChemicalStates = Array.Empty<object>(),
            proteinStructuralPolicies = Array.Empty<object>(), membranePolicies = new[] { policy },
            placementPolicies = Array.Empty<object>(), placementWitnesses = Array.Empty<object>(),
            preparationPolicies = Array.Empty<object>(), equilibrationQualifications = Array.Empty<object>(),
            ppmVersion = "", ppmExecutableSha256 = "",
            maximumSourceAtoms = 1000
        }, new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        return path;
    }

    private static string Hash(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

    private static Task<BoundaryOutcome<WorkspaceState>> Command(ProductRoot product,
        ActorActionKind kind, object data) => product.ExecuteAsync(new ActorCommand(kind,
            JsonSerializer.SerializeToElement(data), product.Snapshot().Revision),
            TestContext.Current.CancellationToken);

    private sealed class NoNetworkHandler : HttpMessageHandler
    {
        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request,
            CancellationToken cancellationToken) => throw new InvalidOperationException("No remote source was requested.");
    }

    private sealed class MembraneRouteWorker : IScientificWorkerExchange
    {
        public List<ScientificWorkRequest<MembraneAssessmentPayload>> MembraneRequests { get; } = [];
        public bool CombinedParameterizationObserved { get; init; } = true;

        public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
            ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken)
        {
            MembraneRequests.Add(request);
            var species = request.Payload.SpeciesRepresentations.Select(item => new SpeciesTemplateObservation(
                item.SpeciesId, item.ChemistryId, item.AtomCount, item.AtomCount,
                true, ImmutableArray<string>.Empty)).ToImmutableArray();
            var observations = new MembraneAssessmentObservations(species, ImmutableArray<string>.Empty,
                CombinedParameterizationObserved);
            return Task.FromResult(new WorkerResult<MembraneAssessmentObservations>(request.RequestId,
                request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Observed,
                ImmutableArray<WorkerArtifact>.Empty, observations,
                new ProviderIdentity("controlled OpenMM observations", "8.6.0"), null, null));
        }

        public Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
            ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken) => NotUsed<SourceInspectionObservations>();
        public Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
            ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken) => NotUsed<PreparationChangeObservations>();
        public Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
            ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken) => NotUsed<ProteinPreparationObservations>();
        public Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
            ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken) => NotUsed<PredictionRegionSummaryObservations>();
        public Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
            ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken) => NotUsed<PlacementObservations>();
        public Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
            ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken) => NotUsed<PlacementAdjustmentObservations>();
        public Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
            ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken) => NotUsed<PlacementMeasurementObservations>();
        public Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
            ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken) => NotUsed<ConstructionObservations>();
        public Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
            ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken) => NotUsed<MinimizationObservations>();
        public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
            ScientificWorkRequest<EquilibrationPayload> request,
            IProgress<EquilibrationWorkProgress>? progress, CancellationToken cancellationToken) => NotUsed<EquilibrationObservations>();
        public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
            ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken) => NotUsed<StageObservationObservations>();
        public Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
            ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken) => NotUsed<ExportVerificationObservations>();
        private static Task<WorkerResult<T>> NotUsed<T>() where T : class =>
            throw new InvalidOperationException("An unrelated scientific operation was invoked.");
    }
}
