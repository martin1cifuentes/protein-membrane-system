using System.Collections.Immutable;
using System.Net;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProteinPreparationOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinPreparation.ProteinPreparation;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed partial class ProteinPreparationRouteTests
{
    private static readonly ResidueAddress SourceResidue = new(0, "A", 1, "", "");
    private static readonly ResidueAddress SelectedResidue = SourceResidue with { CopyId = "A" };
    private const string GeometryReason = "No declared observation radius exists for element 'X'.";

    [Fact]
    public async Task Actor_can_select_an_exact_database_reference_without_a_discovery_candidate()
    {
        using var directory = new TemporaryDirectory();
        var requested = new List<Uri>();
        using var http = new HttpClient(new RespondingHandler(request =>
        {
            requested.Add(request.RequestUri!);
            return new HttpResponseMessage(HttpStatusCode.OK)
            {
                Content = new ByteArrayContent(Encoding.ASCII.GetBytes("data_exact_source\n#\n"))
            };
        }));
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("mmcif", ImmutableArray.Create(Model()), null))
        };
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            string.Empty, string.Empty, () => null);

        var invalid = await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "1ABC/other" });
        Assert.False(invalid.Established);
        Assert.Empty(requested);
        var selected = await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "1ABC" });

        Assert.True(selected.Established);
        Assert.Empty(selected.Value!.SourceCandidates);
        Assert.Equal("rcsb:1ABC", selected.Value.Study?.SelectedSourceId);
        Assert.Equal(SourceRouteKind.Rcsb, selected.Value.Study?.SelectedSourceKind);
        Assert.Single(selected.Value.SourceModels);
        Assert.Single(requested);
        Assert.Equal("https://files.rcsb.org/download/1ABC.cif", requested[0].AbsoluteUri);
    }

    [Fact]
    public async Task Failed_exact_download_is_a_retrieval_refusal_without_a_protein_finding()
    {
        using var directory = new TemporaryDirectory();
        using var http = new HttpClient(new RespondingHandler(_ =>
            new HttpResponseMessage(HttpStatusCode.ServiceUnavailable)));
        var product = new ProductRoot(new PreparationWorkerStub(),
            new ExternalSourceExchange(http), directory.Path,
            string.Empty, string.Empty, () => null);

        var outcome = await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "1ABC" });

        Assert.False(outcome.Established);
        Assert.Contains("could not be retrieved", outcome.Reason);
        Assert.Empty(outcome.Findings);
        Assert.Null(product.Snapshot().Study?.SelectedSourceId);
        Assert.Null(product.Snapshot().Protein);
    }

    [Fact]
    public async Task Source_inspection_requires_correlated_worker_identity_and_observations()
    {
        using var directory = new TemporaryDirectory();
        var path = System.IO.Path.Combine(directory.Path, "source.pdb");
        File.WriteAllText(path, "source coordinates");
        var source = new StructuralSource("upload:fixture", SourceRouteKind.Upload, "researcher upload",
            path, Hash(path), null, "fixture", null, UploadOriginKind.Experimental);
        var observations = new SourceInspectionObservations("pdb", ImmutableArray.Create(Model()), null);
        var wrongIdWorker = new PreparationWorkerStub
        {
            InspectSourceResponse = _ => Observed("different-request", null, observations)
        };
        var wrongId = await new ProteinPreparationOwner(wrongIdWorker).InspectSourceAsync(source,
            directory.Path, 1000, TestContext.Current.CancellationToken);
        Assert.False(wrongId.Established);
        Assert.Contains("not observed", wrongId.Reason);

        var emptyWorker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => new WorkerResult<SourceInspectionObservations>(
                request.RequestId, null, null, null, WorkerResultStanding.Observed,
                ImmutableArray<WorkerArtifact>.Empty, null, null, null, null)
        };
        var empty = await new ProteinPreparationOwner(emptyWorker).InspectSourceAsync(source,
            directory.Path, 1000, TestContext.Current.CancellationToken);
        Assert.False(empty.Established);
        Assert.Contains("not observed", empty.Reason);
    }

    [Fact]
    public async Task Declined_required_heavy_atom_is_visible_on_current_proposal_without_preparing_a_candidate()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(Model()), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(
                        ImmutableArray.Create(new AtomAddress(SelectedResidue, "CB")),
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 4,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);

        var uploadToken = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "source.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        var selected = await Command(product, ActorActionKind.SelectSource, new { uploadToken });
        Assert.True(selected.Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established);
        var proposal = Assert.Single(modeled.Value!.Protein!.Changes);
        Assert.Equal(PreparationChangeKind.HeavyAtom, proposal.Kind);
        var revisionId = modeled.Value.Study!.Id;

        var prematureApproval = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true, rationale = "Approve exact CB" });
        Assert.False(prematureApproval.Established);
        Assert.Contains("not been selected for inspection", prematureApproval.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = proposal.Id });
        Assert.True(inspected.Established);
        Assert.Equal(proposal.Id, inspected.Value!.Inspection?.SubjectId);
        Assert.Contains(inspected.Value.Actions, action =>
            action.Kind == ActorActionKind.ApprovePreparationChange &&
            action.SubjectId == proposal.Id && action.Enabled);

        var declined = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = false, rationale = "Do not add this atom" });
        Assert.True(declined.Established);
        Assert.Equal("declined", declined.Value!.Protein!.Status);
        Assert.Contains(proposal.Id, declined.Value.Protein.Summary);
        Assert.Contains(revisionId, declined.Value.Protein.Summary);
        Assert.Contains(declined.Value.Notices, notice => notice.SubjectId == proposal.Id &&
            notice.Message.Contains("no assessed prepared protein", StringComparison.Ordinal));
        Assert.Equal(0, worker.PrepareProteinCalls);
        Assert.DoesNotContain(declined.Value.Actions, action =>
            action.Kind == ActorActionKind.ProposePlacement && action.Enabled);

        var stale = await product.ExecuteAsync(new ActorCommand(ActorActionKind.ApprovePreparationChange,
            JsonSerializer.SerializeToElement(new { proposalId = proposal.Id, approve = true }),
            modeled.Value.Revision), TestContext.Current.CancellationToken);
        Assert.False(stale.Established);
        Assert.Contains("workspace changed", stale.Reason);

        var replacementToken = await product.UploadAsync(
            new MemoryStream(Encoding.ASCII.GetBytes("DIFFERENT ATOM\n")), "replacement.pdb",
            UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        var replaced = await Command(product, ActorActionKind.SelectSource,
            new { uploadToken = replacementToken });
        Assert.True(replaced.Established);
        Assert.NotEqual(revisionId, replaced.Value!.Study!.Id);
        Assert.Null(replaced.Value.Inspection);
        Assert.Null(replaced.Value.Protein);
        var oldFocus = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = proposal.Id });
        Assert.False(oldFocus.Established);
        Assert.Contains("exact inspection subject is unavailable", oldFocus.Reason);
        var oldDecision = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true, rationale = "Stale evidence" });
        Assert.False(oldDecision.Established);
        Assert.Contains("absent or stale", oldDecision.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Mismatched_force_field_asset_keeps_selected_subject_but_prevents_proposal_and_preparation()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        File.AppendAllText(System.IO.Path.Combine(directory.Path, "forcefield.xml"), "changed bytes");
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(Model()), null))
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "source.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        var selected = await Command(product, ActorActionKind.SelectSource, new { uploadToken = token });
        Assert.True(selected.Established);

        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });

        Assert.True(modeled.Established);
        Assert.Equal(0, modeled.Value!.Study?.ModelIndex);
        Assert.Equal("notEstablished", modeled.Value.Protein?.Status);
        Assert.Empty(modeled.Value.Protein!.Changes);
        Assert.Contains(modeled.Value.Notices, notice =>
            notice.Message.Contains("No qualified protein chemical-state", StringComparison.Ordinal));
        Assert.Equal(0, worker.InspectChangesCalls);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Absent_policy_catalogue_is_visible_and_cannot_start_protein_preparation()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(Model()), null))
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            string.Empty, string.Empty, () => null);
        Assert.Contains(product.Snapshot().Notices, notice =>
            notice.Message.Contains("No qualified local policy catalogue", StringComparison.Ordinal));
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "source.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        var selected = await Command(product, ActorActionKind.SelectSource, new { uploadToken = token });
        Assert.True(selected.Established);

        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });

        Assert.True(modeled.Established);
        Assert.Equal("notEstablished", modeled.Value!.Protein?.Status);
        Assert.Contains(modeled.Value.Notices, notice => notice.SubjectId == modeled.Value.Protein!.SubjectId &&
            notice.Message.Contains("No qualified protein chemical-state", StringComparison.Ordinal));
        Assert.Equal(0, worker.InspectChangesCalls);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    private static WorkerResult<ProteinPreparationObservations> CandidateWithUnavailableGeometry(
        ScientificWorkRequest<ProteinPreparationPayload> request, string directory, string sourceHash)
        => CandidateResponse(request, directory, sourceHash, ["N", "CA", "C", "O"],
            new ProteinGeometryObservations(ObservationStanding.Unavailable,
                ImmutableArray.Create(new ProteinGeometryKindObservation("covalentBond",
                    GeometryKindStanding.Unavailable, 0, 0, null, null, GeometryReason)),
                ImmutableArray<GeometryDistanceObservation>.Empty, ImmutableArray.Create(GeometryReason)),
            ImmutableArray.Create(GeometryReason));

    private static WorkerResult<ProteinPreparationObservations> CandidateResponse(
        ScientificWorkRequest<ProteinPreparationPayload> request, string directory, string sourceHash,
        string[] atomNames, ProteinGeometryObservations geometry,
        ImmutableArray<string> geometryWarnings = default,
        ImmutableArray<ResidueVariantChoice> actualVariants = default,
        bool correspondenceComplete = true)
    {
        var prepared = System.IO.Path.Combine(directory, "prepared.pdb");
        var graph = System.IO.Path.Combine(directory, "bonds.json");
        var mapping = System.IO.Path.Combine(directory, "correspondence.json");
        File.WriteAllText(prepared, "prepared coordinates");
        File.WriteAllText(graph, "{}");
        var preparedHash = Hash(prepared);
        var atomMapping = atomNames.Select((name, index) => new AtomCorrespondence(
            index, $"result:{index}:{name}", $"source:{index}:{name}", AtomOriginKind.Source,
            MoleculeRoleKind.Protein, name is "N" or "CA" or "C" or "O"
                ? AtomRoleKind.Backbone : AtomRoleKind.Sidechain,
            name == "N" ? "N" : name is "O" ? "O" : name == "SG" ? "S" : "C",
            SelectedResidue, null)).ToImmutableArray();
        File.WriteAllText(mapping, JsonSerializer.Serialize(
            new SourceToResultCorrespondence(sourceHash, preparedHash, atomMapping, correspondenceComplete),
            new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        var observations = new ProteinPreparationObservations(atomNames.Length, atomNames.Length, 1,
            ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
            ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
            ImmutableArray<AtomAddress>.Empty, ImmutableArray<AtomAddress>.Empty,
            ImmutableArray<AtomAddress>.Empty,
            actualVariants.IsDefault ? ImmutableArray<ResidueVariantChoice>.Empty : actualVariants,
            ImmutableArray<DisulfideBond>.Empty,
            geometryWarnings.IsDefault ? ImmutableArray<string>.Empty : geometryWarnings,
            atomNames.Length, geometry);
        return Observed(request.RequestId, request.Payload.StudyRevisionId, observations,
            ImmutableArray.Create(Artifact(prepared, "preparedPdb"),
                Artifact(graph, "preparedBondGraph"), Artifact(mapping, "correspondenceJson")));
    }

    private static async Task<BoundaryOutcome<WorkspaceState>> Command(
        ProductRoot product, ActorActionKind kind, object data) =>
        await product.ExecuteAsync(new ActorCommand(kind, JsonSerializer.SerializeToElement(data),
            product.Snapshot().Revision), TestContext.Current.CancellationToken);

    private static SourceModelObservation Model() => new(0,
        ImmutableArray.Create(new SourceChainObservation("A", 1, 4)),
        ImmutableArray<SourceAssemblyObservation>.Empty,
        ImmutableArray<SourcePartnerObservation>.Empty,
        ImmutableArray.Create(new SourceResidueObservation(SourceResidue, "ALA", SourceResidueKind.Protein,
            true, ImmutableArray<string>.Empty, ImmutableArray<string>.Empty, null,
            ObservationStanding.Unavailable, ObservationStanding.Unavailable,
            ImmutableArray<string>.Empty)), 4);

    private static (StudyRevision Revision, IntendedProteinModel Intended,
        SourceInspectionReport Inspection, PreparationProposalReport Proposals,
        ProteinChemicalStatePolicy Chemical) OwnerContext(string directory, string residueName, string[] atomNames)
    {
        var path = System.IO.Path.Combine(directory, "source.pdb");
        File.WriteAllText(path, "source coordinates");
        var source = new StructuralSource("upload:fixture", SourceRouteKind.Upload, "researcher upload",
            path, Hash(path), null, "fixture", null, UploadOriginKind.Experimental);
        var intended = new IntendedProteinModel("intended", source, 0, null,
            ImmutableArray.Create(new ChainSelection("A", "A")),
            ImmutableArray<PartnerSelection>.Empty, ImmutableArray<AlternateLocationChoice>.Empty);
        var revision = new StudyRevision("revision", 1, intended, null, null, FixedStudyConditions.Initial);
        var model = new SourceModelObservation(0,
            ImmutableArray.Create(new SourceChainObservation("A", 1, atomNames.Length)),
            ImmutableArray<SourceAssemblyObservation>.Empty, ImmutableArray<SourcePartnerObservation>.Empty,
            ImmutableArray.Create(new SourceResidueObservation(SourceResidue, residueName, SourceResidueKind.Protein,
                true, ImmutableArray<string>.Empty, ImmutableArray<string>.Empty, null,
                ObservationStanding.Unavailable, ObservationStanding.Unavailable, ImmutableArray<string>.Empty)),
            atomNames.Length);
        var inspection = new SourceInspectionReport(source, "pdb", ImmutableArray.Create(model),
            ImmutableArray<string>.Empty);
        var proposals = new PreparationProposalReport(revision.Id, intended.Id,
            ImmutableArray<PreparationChangeProposal>.Empty, null,
            ImmutableArray<PreviewChainCorrespondence>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<string>.Empty);
        var chemical = new ProteinChemicalStatePolicy("chemical", "1", "canonical-amino-acid-assembly", 2.5,
            ImmutableArray.Create("policy evidence"),
            ImmutableArray.Create(new ForceFieldAsset("ff", "1", "Amber19", "/unused/ff.xml", "abc")),
            ImmutableDictionary<string, string>.Empty.Add("CYS", "CYS"),
            ImmutableArray.Create("CYS", "HID"), ImmutableArray<string>.Empty);
        return (revision, intended, inspection, proposals, chemical);
    }

    private static ProteinGeometryObservations ObservedGeometry() => new(ObservationStanding.Observed,
        ImmutableArray.Create(new ProteinGeometryKindObservation("covalentBond", GeometryKindStanding.Observed,
            1, 1, 1.4, 1.4, null)), ImmutableArray<GeometryDistanceObservation>.Empty,
        ImmutableArray<string>.Empty);

    private static ProteinStructuralAssessmentPolicy StructuralPolicy() => new("structural", "1",
        "canonical-amino-acid-assembly", ImmutableArray.Create("policy evidence"),
        new ProteinGeometryMeasurementSpec(ImmutableArray.Create("covalentBond"),
            ImmutableDictionary<string, double>.Empty.Add("C", 1.7).Add("N", 1.5), 4.0, 3, 10),
        ImmutableArray.Create(new ProteinGeometryCriterion("covalentBond", 1.0, 2.0, false)),
        ImmutableArray<string>.Empty);

    private static string WritePolicyCatalogue(string directory)
    {
        var asset = System.IO.Path.Combine(directory, "forcefield.xml");
        File.WriteAllText(asset, "<ForceField/>");
        var chemical = new ProteinChemicalStatePolicy("chemical", "1", "canonical-amino-acid-assembly",
            2.5, ImmutableArray.Create("policy evidence"),
            ImmutableArray.Create(new ForceFieldAsset("ff", "1", "Amber19", asset, Hash(asset))),
            ImmutableDictionary<string, string>.Empty, ImmutableArray.Create("HID"), ImmutableArray<string>.Empty);
        var catalogue = System.IO.Path.Combine(directory, "catalogue.json");
        File.WriteAllText(catalogue, JsonSerializer.Serialize(new
        {
            version = "1", evidenceReferences = new[] { "catalogue evidence" },
            lipids = Array.Empty<object>(), proteinChemicalStates = new[] { chemical },
            proteinStructuralPolicies = new[] { StructuralPolicy() },
            membranePolicies = Array.Empty<object>(), placementPolicies = Array.Empty<object>(),
            placementWitnesses = Array.Empty<object>(), preparationPolicies = Array.Empty<object>(),
            equilibrationQualifications = Array.Empty<object>(), ppmVersion = "", ppmExecutableSha256 = "",
            maximumSourceAtoms = 1000
        }, new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        return catalogue;
    }

    private static WorkerArtifact Artifact(string path, string role) => new(role, path, Hash(path));
    private static string Hash(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

    private static WorkerResult<T> Observed<T>(string requestId, string? revision, T observations,
        ImmutableArray<WorkerArtifact> artifacts = default) where T : class =>
        new(requestId, revision, null, null, WorkerResultStanding.Observed,
            artifacts.IsDefault ? ImmutableArray<WorkerArtifact>.Empty : artifacts,
            observations, new ProviderIdentity("controlled observed worker", "1"), null, null);

    private sealed class NoNetworkHandler : HttpMessageHandler
    {
        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request,
            CancellationToken cancellationToken) => throw new InvalidOperationException("No remote source was requested.");
    }

    private sealed class RespondingHandler(Func<HttpRequestMessage, HttpResponseMessage> respond) : HttpMessageHandler
    {
        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request,
            CancellationToken cancellationToken) => Task.FromResult(respond(request));
    }

    private sealed class PreparationWorkerStub : IScientificWorkerExchange
    {
        public Func<ScientificWorkRequest<SourceInspectionPayload>, WorkerResult<SourceInspectionObservations>>?
            InspectSourceResponse { get; init; }
        public Func<ScientificWorkRequest<PreparationChangeInspectionPayload>, WorkerResult<PreparationChangeObservations>>?
            InspectChangesResponse { get; init; }
        public Func<ScientificWorkRequest<ProteinPreparationPayload>, WorkerResult<ProteinPreparationObservations>>?
            PrepareResponse { get; init; }
        public int PrepareProteinCalls { get; private set; }
        public int InspectChangesCalls { get; private set; }

        public Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
            ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken) =>
            Task.FromResult(InspectSourceResponse!(request));
        public Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
            ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken)
        {
            InspectChangesCalls++;
            return Task.FromResult(InspectChangesResponse!(request));
        }
        public Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
            ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken)
        {
            PrepareProteinCalls++;
            return Task.FromResult(PrepareResponse!(request));
        }
        public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
            ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken) => NotUsed<MembraneAssessmentObservations>();
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
