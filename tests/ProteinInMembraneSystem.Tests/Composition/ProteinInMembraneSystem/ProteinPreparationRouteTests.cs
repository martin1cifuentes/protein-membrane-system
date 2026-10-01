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
    [Fact]
    public void Preparation_review_roles_keep_the_existing_browser_wire_values()
    {
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        Assert.Equal("\"beforePreparation\"", JsonSerializer.Serialize(
            DecisionInspectionRelation.BeforePreparation, options));
        Assert.Equal("\"preparedResult\"", JsonSerializer.Serialize(
            DecisionInspectionRelation.PreparedResult, options));
        Assert.Equal("\"modelAssumption\"", JsonSerializer.Serialize(
            PreparationInformationRole.ModelAssumption, options));
        Assert.Equal("\"observed\"", JsonSerializer.Serialize(
            PreparationInformationRole.Observed, options));
    }

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
        Assert.Equal("rcsb:1ABC", selected.Value.Inspection?.SubjectId);
        Assert.Equal(selected.Value.Study?.Id, selected.Value.Inspection?.StudyRevisionId);
        Assert.NotNull(selected.Value.Inspection?.StructureUrl);
        Assert.Single(selected.Value.SourceModels);
        Assert.Single(requested);
        Assert.Equal("https://files.rcsb.org/download/1ABC.cif", requested[0].AbsoluteUri);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public async Task Explicit_partner_disposition_needs_no_prose_but_each_partner_must_be_decided(bool retain)
    {
        using var directory = new TemporaryDirectory();
        using var http = new HttpClient(new RespondingHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new ByteArrayContent(Encoding.ASCII.GetBytes("data_exact_source\n#\n"))
        }));
        var model = Model() with
        {
            Partners = ImmutableArray.Create(new SourcePartnerObservation("partner-1", "Bound cofactor",
                "heterogen", 4, "A", 1))
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("mmcif", ImmutableArray.Create(model), null))
        };
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            string.Empty, string.Empty, () => null);
        Assert.True((await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "1ABC" })).Established);
        var undecided = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<PartnerSelection>(), alternateLocations = Array.Empty<AlternateLocationChoice>() });
        Assert.False(undecided.Established);
        var missingChains = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, chains = Array.Empty<ChainSelection>(),
                partners = new[] { new { sourceId = "partner-1", retain } },
                alternateLocations = Array.Empty<AlternateLocationChoice>() });
        Assert.False(missingChains.Established);
        var selected = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, chains = new[] { new ChainSelection("A", "A") },
                partners = new[] { new { sourceId = "partner-1", retain } },
                alternateLocations = Array.Empty<AlternateLocationChoice>() });
        Assert.True(selected.Established, selected.Reason);
        Assert.Equal(retain, Assert.Single(selected.Value!.Study!.Partners).Retain);
        Assert.Null(Assert.Single(selected.Value.Study.Partners).Reason);
    }

    [Fact]
    public async Task Unsupported_retained_molecules_are_named_without_changing_selection_or_starting_preparation()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        var model = Model() with
        {
            Partners = ImmutableArray.Create(
                new SourcePartnerObservation("heme-a", "HEM", "nonpolymer", 43, "A", 142,
                    DisplayName: "Heme"),
                new SourcePartnerObservation("heme-b", "HEM", "nonpolymer", 43, "A", 145,
                    DisplayName: "Heme"),
                new SourcePartnerObservation("water-a", "HOH", "water", 1, "A", 143,
                    DisplayName: "Water"),
                new SourcePartnerObservation("other-a", "X1Z", "nonpolymer", 10, "A", 144,
                    DisplayName: "123"))
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null))
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "source.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource,
            new { uploadToken = token })).Established);
        var selectedPartners = new[]
        {
            new PartnerSelection("heme-a", true), new PartnerSelection("heme-b", true),
            new PartnerSelection("water-a", true),
            new PartnerSelection("other-a", true)
        };

        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, chains = new[] { new ChainSelection("A", "A") },
                partners = selectedPartners, alternateLocations = Array.Empty<AlternateLocationChoice>() });

        Assert.True(modeled.Established, modeled.Reason);
        Assert.Equal(selectedPartners, modeled.Value!.Study!.Partners);
        Assert.Equal("notEstablished", modeled.Value.Protein?.Status);
        Assert.Equal("blocked", modeled.Value.ProteinTask?.Standing);
        Assert.Equal("This preparation method does not yet support the selected molecules: " +
            "Heme (HEM) · Chain A · Residue 142, Heme (HEM) · Chain A · Residue 145, " +
            "X1Z · Chain A · Residue 144.",
            modeled.Value.ProteinTask?.Message);
        Assert.True(modeled.Value.ProteinTask?.ReviewMoleculeSelection);
        Assert.Contains(modeled.Value.Actions, action => action.Kind == ActorActionKind.StartProteinPreparation &&
            !action.Enabled);
        var refused = await Command(product, ActorActionKind.StartProteinPreparation, new { });
        Assert.False(refused.Established);
        Assert.Equal(selectedPartners, product.Snapshot().Study!.Partners);
        Assert.Equal(0, worker.InspectChangesCalls);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Reviewed_histidine_variant_needs_no_written_reason_or_inspection_click()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        var histidineModel = Model() with
        {
            Residues = ImmutableArray.Create(Model().Residues[0] with { Name = "HIS" })
        };
        ProteinPreparationPayload? submitted = null;
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(histidineModel), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(
                        ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 4,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            },
            PrepareResponse = request =>
            {
                submitted = request.Payload;
                return new WorkerResult<ProteinPreparationObservations>(request.RequestId,
                    request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
                    ImmutableArray<WorkerArtifact>.Empty, null, null,
                    "controlled-provider-failure", "The controlled provider made no candidate.");
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var uploadToken = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "histidine.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        var proposal = Assert.Single(modeled.Value!.Protein!.Changes);
        Assert.Equal(PreparationChangeKind.ResidueState, proposal.Kind);

        Assert.Contains(modeled.Value.Actions, available => available.Kind == ActorActionKind.ApprovePreparationChange &&
            available.SubjectId == proposal.Id && available.Enabled);
        Assert.NotEqual(proposal.Id, modeled.Value.Inspection?.SubjectId);
        var approved = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true });
        Assert.True(approved.Established, approved.Reason);
        Assert.Equal(1, worker.PrepareProteinCalls);
        var variant = Assert.Single(submitted!.ResidueVariants);
        Assert.Equal("HID", variant.Variant);
        Assert.Equal(proposal.Residue, variant.Residue);
        Assert.NotEqual("assessed", approved.Value!.Protein?.Status);
        Assert.Equal("failed", approved.Value.PreparationReview?.PreparationStanding);
        Assert.Equal("failed", approved.Value.ProteinTask?.Standing);
        Assert.True(approved.Value.ProteinTask?.RetryAvailable);
        Assert.Equal(1, approved.Value.PreparationReview?.ConfirmedCount);
        Assert.Contains("controlled provider made no candidate", approved.Value.PreparationReview?.PreparationMessage,
            StringComparison.OrdinalIgnoreCase);
        Assert.Contains(approved.Value.Actions, available => available.Kind == ActorActionKind.RetryProteinPreparation &&
            available.Enabled);
        var duplicate = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true });
        Assert.False(duplicate.Established);
        Assert.Equal(1, worker.PrepareProteinCalls);
        var retried = await Command(product, ActorActionKind.RetryProteinPreparation, new { });
        Assert.True(retried.Established, retried.Reason);
        Assert.Equal(2, worker.PrepareProteinCalls);
        Assert.Equal(1, retried.Value!.PreparationReview?.ConfirmedCount);
    }

    [Fact]
    public async Task Missing_exact_chain_copy_mapping_refuses_the_proposal_report_before_review()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path, "HID", "HIE", "HIP");
        var first = Model().Residues[0] with { Name = "HIS" };
        var second = first with { Address = SourceResidue with { Chain = "B" } };
        var model = Model() with
        {
            Chains = ImmutableArray.Create(new SourceChainObservation("A", 1, 4),
                new SourceChainObservation("B", 1, 4)),
            Residues = ImmutableArray.Create(first, second), AtomCount = 8
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 8,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "two-site.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A"), new ChainSelection("B", "B") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        Assert.Null(modeled.Value!.PreparationReview);
        Assert.Contains(modeled.Value.Notices, notice =>
            notice.Message.Contains("exact selected chain copies", StringComparison.OrdinalIgnoreCase));
        Assert.DoesNotContain(modeled.Value.Actions, action =>
            action.Kind == ActorActionKind.ApprovePreparationChange && action.Enabled);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Review_groups_exact_histidine_sites_and_keeps_independent_repair_separate()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path, "HID", "HIE", "HIP");
        var first = Model().Residues[0] with { Name = "HIS" };
        var second = first with { Address = SourceResidue with { Chain = "B" } };
        var model = Model() with
        {
            Chains = ImmutableArray.Create(new SourceChainObservation("A", 1, 4),
                new SourceChainObservation("B", 1, 4)),
            Residues = ImmutableArray.Create(first, second), AtomCount = 8
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(
                        ImmutableArray.Create(new AtomAddress(SelectedResidue, "CB")),
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 8,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A"),
                            new PreviewChainCorrespondence("B", "B", "B")), ObservedGeometry()),
                    ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var uploadToken = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "two-histidines.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A"), new ChainSelection("B", "B") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        var review = modeled.Value!.PreparationReview!;
        Assert.Equal(7, modeled.Value.Protein!.Changes.Length);
        Assert.Equal(3, review.Decisions.Length);
        Assert.Equal(3, review.RemainingCount);
        var siteA = Assert.Single(review.Decisions.Where(item => item.Kind == PreparationChangeKind.ResidueState &&
            item.Residue.Chain == "A"));
        var siteB = Assert.Single(review.Decisions.Where(item => item.Kind == PreparationChangeKind.ResidueState &&
            item.Residue.Chain == "B"));
        Assert.Equal(3, siteA.Options.Length);
        Assert.Equal(3, siteB.Options.Length);
        Assert.NotEqual(siteA.Id, siteB.Id);
        Assert.Single(review.Decisions.Where(item => item.Kind == PreparationChangeKind.HeavyAtom));

        var chosen = siteA.Options.Single(item => item.ProposedChange == "HIE");
        var confirmed = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = chosen.ProposalId, approve = true });
        Assert.True(confirmed.Established, confirmed.Reason);
        Assert.Equal(modeled.Value.Inspection?.SubjectId, confirmed.Value?.Inspection?.SubjectId);
        var competing = siteA.Options.Single(item => item.ProposedChange == "HID");
        var competingConfirmation = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = competing.ProposalId, approve = true });
        Assert.False(competingConfirmation.Established);
        Assert.Contains("alternative", competingConfirmation.Reason, StringComparison.OrdinalIgnoreCase);
        review = confirmed.Value!.PreparationReview!;
        Assert.Equal(1, review.ConfirmedCount);
        Assert.Equal(2, review.RemainingCount);
        siteA = review.Decisions.Single(item => item.Id == siteA.Id);
        siteB = review.Decisions.Single(item => item.Id == siteB.Id);
        Assert.Equal("confirmed", siteA.Standing);
        Assert.Equal(chosen.ProposalId, siteA.ChosenProposalId);
        Assert.Equal(2, siteA.Options.Count(item => item.Disposition == "notChosen"));
        Assert.Equal("pending", siteB.Standing);
        Assert.Equal(0, worker.PrepareProteinCalls);

        var rejected = siteB.Options.Single(item => item.ProposedChange == "HIE");
        var declined = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = rejected.ProposalId, approve = false });
        Assert.False(declined.Established);
        Assert.Contains("cannot be declined", declined.Reason, StringComparison.OrdinalIgnoreCase);
        siteB = product.Snapshot().PreparationReview!.Decisions.Single(item => item.Id == siteB.Id);
        Assert.Equal("pending", siteB.Standing);
        Assert.Equal(3, siteB.Options.Count(item => item.Disposition == "available"));

        var repair = product.Snapshot().PreparationReview!.Decisions.Single(item => item.Kind == PreparationChangeKind.HeavyAtom);
        var blocked = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = repair.Options[0].ProposalId, approve = false });
        Assert.True(blocked.Established, blocked.Reason);
        Assert.Equal("blocked", blocked.Value!.PreparationReview?.PreparationStanding);
        Assert.Equal("blocked", blocked.Value.PreparationReview?.Decisions.Single(item => item.Id == repair.Id).Standing);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Last_confirmation_exposes_real_preparation_activity_then_retains_failure_and_choice()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        var histidine = Model() with
        {
            Residues = ImmutableArray.Create(Model().Residues[0] with { Name = "HIS" })
        };
        var entered = new TaskCompletionSource<ScientificWorkRequest<ProteinPreparationPayload>>(
            TaskCreationOptions.RunContinuationsAsynchronously);
        var released = new TaskCompletionSource<WorkerResult<ProteinPreparationObservations>>(
            TaskCreationOptions.RunContinuationsAsynchronously);
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(histidine), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 4,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            },
            PrepareResponseAsync = (request, _) =>
            {
                entered.TrySetResult(request);
                return released.Task;
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "histidine.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        var option = Assert.Single(Assert.Single(modeled.Value!.PreparationReview!.Decisions).Options);
        Assert.True(product.Snapshot().PreparationReview!.Decisions[0].Options[0].StartsPreparationOnConfirmation);
        var confirmation = Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = option.ProposalId, approve = true });
        var request = await entered.Task.WaitAsync(TimeSpan.FromSeconds(5),
            TestContext.Current.CancellationToken);
        Assert.Equal("HID", Assert.Single(request.Payload.ResidueVariants).Variant);
        var preparing = product.Snapshot().PreparationReview!;
        Assert.Equal("preparing", preparing.PreparationStanding);
        Assert.Equal("preparing", product.Snapshot().ProteinTask?.Standing);
        Assert.Equal(1, preparing.ConfirmedCount);
        Assert.Equal(0, preparing.RemainingCount);
        Assert.Null(product.Snapshot().Protein?.AtomCount);
        released.SetResult(new WorkerResult<ProteinPreparationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
            ImmutableArray<WorkerArtifact>.Empty, null, null,
            "controlled-provider-failure", "The controlled provider made no candidate."));
        var failed = await confirmation;
        Assert.True(failed.Established, failed.Reason);
        Assert.Equal("failed", failed.Value!.PreparationReview?.PreparationStanding);
        Assert.Equal("failed", failed.Value.ProteinTask?.Standing);
        Assert.True(failed.Value.ProteinTask?.RetryAvailable);
        Assert.Equal(option.ProposalId, failed.Value.PreparationReview?.Decisions[0].ChosenProposalId);
        Assert.Contains("controlled provider made no candidate", failed.Value.PreparationReview?.PreparationMessage,
            StringComparison.OrdinalIgnoreCase);
        Assert.Equal(1, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Review_scope_keeps_model_copy_and_insertion_code_distinct()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path, "HID", "HIE", "HIP");
        var site = Model().Residues[0] with { Name = "HIS" };
        var model0 = Model() with
        {
            Assemblies = ImmutableArray.Create(new SourceAssemblyObservation("dimer",
                ImmutableArray.Create(new ChainSelection("A", "A1"), new ChainSelection("A", "A2")))),
            Residues = ImmutableArray.Create(site,
                site with { Address = SourceResidue with { InsertionCode = "X" } }),
            Chains = ImmutableArray.Create(new SourceChainObservation("A", 2, 8)), AtomCount = 8
        };
        var model1 = Model() with
        {
            Index = 1,
            Residues = ImmutableArray.Create(site with { Address = SourceResidue with { Model = 1 } })
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model0, model1), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 8,
                        request.Payload.ChainSelections.Select((chain, index) =>
                            new PreviewChainCorrespondence(chain.SourceChain, chain.CopyId,
                                ((char)('A' + index)).ToString())).ToImmutableArray(),
                        ObservedGeometry()),
                    ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "two-models.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var first = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = "dimer",
                chains = new[] { new ChainSelection("A", "A1"), new ChainSelection("A", "A2") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(first.Established, first.Reason);
        var firstReview = first.Value!.PreparationReview!;
        Assert.Equal(4, firstReview.Decisions.Length);
        Assert.Equal(12, first.Value.Protein!.Changes.Length);
        Assert.Equal(4, firstReview.Decisions.Select(item =>
            (item.Residue.Model, item.Residue.Chain, item.Residue.CopyId,
                item.Residue.Residue, item.Residue.InsertionCode)).Distinct().Count());
        Assert.Equal(new[] { "A1", "A2" }, firstReview.Decisions.Select(item => item.Residue.CopyId)
            .Distinct().OrderBy(item => item));
        Assert.Equal(new[] { "", "X" }, firstReview.Decisions.Select(item => item.Residue.InsertionCode)
            .Distinct().OrderBy(item => item));

        var second = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 1, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(second.Established, second.Reason);
        Assert.NotEqual(first.Value.Study!.Id, second.Value!.Study!.Id);
        Assert.Equal(1, Assert.Single(second.Value.PreparationReview!.Decisions).Residue.Model);
        Assert.Equal(3, second.Value.Protein!.Changes.Length);
    }

    [Fact]
    public async Task Duplicate_exact_site_alternatives_block_decision_and_preparation()
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path, "HID", "HIE", "HIP");
        var site = Model().Residues[0] with { Name = "HIS" };
        var model = Model() with { Residues = ImmutableArray.Create(site, site), AtomCount = 8 };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null)),
            InspectChangesResponse = request =>
            {
                var preview = System.IO.Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
                File.WriteAllText(preview, "selected coordinates");
                return Observed(request.RequestId, request.Payload.StudyRevisionId,
                    new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 8,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            }
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "duplicate-histidine.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        var review = modeled.Value!.PreparationReview!;
        var decision = Assert.Single(review.Decisions);
        Assert.Equal(6, decision.Options.Length);
        Assert.Equal("blocked", decision.Standing);
        Assert.Equal("blocked", review.PreparationStanding);
        Assert.Contains("duplicate alternatives", decision.Blocker);
        var option = decision.Options[0];
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = option.ProposalId });
        Assert.True(inspected.Established);
        Assert.Contains(inspected.Value!.Actions, action => action.Kind == ActorActionKind.ApprovePreparationChange &&
            action.SubjectId == option.ProposalId && !action.Enabled);
        var refused = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = option.ProposalId, approve = true });
        Assert.False(refused.Established);
        Assert.Contains("duplicate alternatives", refused.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Failed_preparation_without_review_choices_retains_retry_route()
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
                    new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                        ImmutableArray<PossibleDisulfideObservation>.Empty,
                        ObservationStanding.Observed, ImmutableArray<string>.Empty, 4,
                        ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                        ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview")));
            },
            PrepareResponse = request => new WorkerResult<ProteinPreparationObservations>(request.RequestId,
                request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
                ImmutableArray<WorkerArtifact>.Empty, null, null,
                "controlled-provider-failure", "The controlled provider made no candidate.")
        };
        using var http = new HttpClient(new NoNetworkHandler());
        var product = new ProductRoot(worker, new ExternalSourceExchange(http), directory.Path,
            catalogue, string.Empty, () => null);
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "no-choice.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new ChainSelection("A", "A") },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        Assert.Empty(modeled.Value!.PreparationReview!.Decisions);
        Assert.Equal("ready", modeled.Value.PreparationReview.PreparationStanding);
        Assert.Equal("ready", modeled.Value.ProteinTask?.Standing);
        Assert.Equal(0, worker.PrepareProteinCalls);
        Assert.Contains(modeled.Value.Actions, action => action.Kind == ActorActionKind.StartProteinPreparation &&
            action.Enabled);
        var started = await Command(product, ActorActionKind.StartProteinPreparation, new { });
        Assert.True(started.Established, started.Reason);
        Assert.Equal("failed", started.Value!.ProteinTask?.Standing);
        Assert.True(started.Value.ProteinTask?.RetryAvailable);
        Assert.Contains(started.Value.Actions, action => action.Kind == ActorActionKind.RetryProteinPreparation &&
            action.Enabled);
        var retried = await Command(product, ActorActionKind.RetryProteinPreparation, new { });
        Assert.True(retried.Established, retried.Reason);
        Assert.Equal(2, worker.PrepareProteinCalls);
        Assert.Equal("failed", retried.Value!.PreparationReview?.PreparationStanding);
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

        Assert.Contains(modeled.Value.Actions, action =>
            action.Kind == ActorActionKind.ApprovePreparationChange &&
            action.SubjectId == proposal.Id && action.Enabled);

        var declined = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = false });
        Assert.True(declined.Established);
        Assert.Equal("declined", declined.Value!.Protein!.Status);
        Assert.Contains(proposal.Id, declined.Value.Protein.Summary);
        Assert.Contains(revisionId, declined.Value.Protein.Summary);
        Assert.Contains(declined.Value.Notices, notice =>
            notice.ConditionKey.StartsWith("preparation-blocker-", StringComparison.Ordinal) &&
            notice.Message.Contains("declined", StringComparison.Ordinal) &&
            notice.AffectedAreas.Contains("placement"));
        Assert.Equal(0, worker.PrepareProteinCalls);
        Assert.DoesNotContain(declined.Value.Actions, action =>
            action.Kind == ActorActionKind.ProposePlacement && action.Enabled);
        var refusedPlacement = await Command(product, ActorActionKind.ProposePlacement,
            new { orientationRoute = "manual", startingPosition = "center",
                offsetXAngstrom = 0, offsetYAngstrom = 0, offsetZAngstrom = 0,
                rotationXDegrees = 0, rotationYDegrees = 0, rotationZDegrees = 0 });
        Assert.False(refusedPlacement.Established);
        Assert.Contains("assessed protein", refusedPlacement.Reason, StringComparison.OrdinalIgnoreCase);

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
        Assert.Equal(replaced.Value.Study.SelectedSourceId, replaced.Value.Inspection?.SubjectId);
        Assert.Equal(replaced.Value.Study.Id, replaced.Value.Inspection?.StudyRevisionId);
        Assert.Null(replaced.Value.Protein);
        Assert.DoesNotContain(replaced.Value.Notices, notice =>
            notice.ConditionKey.StartsWith("preparation-blocker-", StringComparison.Ordinal));
        var oldFocus = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = proposal.Id });
        Assert.False(oldFocus.Established);
        Assert.Contains("exact inspection subject is unavailable", oldFocus.Reason);
        var oldDecision = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true });
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
        var model = Model() with
        {
            Partners = ImmutableArray.Create(new SourcePartnerObservation("heme-a", "HEM", "nonpolymer",
                43, "A", 142, DisplayName: "Heme"))
        };
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null))
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
                partners = new[] { new PartnerSelection("heme-a", true) },
                alternateLocations = Array.Empty<object>() });

        Assert.True(modeled.Established);
        Assert.Equal(0, modeled.Value!.Study?.ModelIndex);
        Assert.Equal("notEstablished", modeled.Value.Protein?.Status);
        Assert.Empty(modeled.Value.Protein!.Changes);
        Assert.Contains(modeled.Value.Notices, notice =>
            notice.Message.Contains("No qualified protein chemical-state", StringComparison.Ordinal));
        Assert.False(modeled.Value.ProteinTask?.ReviewMoleculeSelection);
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
        Assert.False(modeled.Value.ProteinTask?.ReviewMoleculeSelection);
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

    private static string WritePolicyCatalogue(string directory, params string[] variants)
    {
        var asset = System.IO.Path.Combine(directory, "forcefield.xml");
        File.WriteAllText(asset, "<ForceField/>");
        var chemical = new ProteinChemicalStatePolicy("chemical", "1", "canonical-amino-acid-assembly",
            2.5, ImmutableArray.Create("policy evidence"),
            ImmutableArray.Create(new ForceFieldAsset("ff", "1", "Amber19", asset, Hash(asset))),
            ImmutableDictionary<string, string>.Empty,
            variants.Length == 0 ? ImmutableArray.Create("HID") : variants.ToImmutableArray(),
            ImmutableArray<string>.Empty);
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
        public Task<WorkerResult<ManualPlacementObservations>> PlaceManualAsync(
            ScientificWorkRequest<ManualPlacementPayload> request, CancellationToken cancellationToken) =>
            Task.FromException<WorkerResult<ManualPlacementObservations>>(new NotSupportedException());
        public Task<WorkerResult<SourcePreviewObservations>> PreviewSourceModelAsync(
            ScientificWorkRequest<SourcePreviewPayload> request, CancellationToken cancellationToken) =>
            Task.FromException<WorkerResult<SourcePreviewObservations>>(new NotSupportedException());
        public Func<ScientificWorkRequest<SourceInspectionPayload>, WorkerResult<SourceInspectionObservations>>?
            InspectSourceResponse { get; init; }
        public Func<ScientificWorkRequest<PreparationChangeInspectionPayload>, WorkerResult<PreparationChangeObservations>>?
            InspectChangesResponse { get; init; }
        public Func<ScientificWorkRequest<ProteinPreparationPayload>, WorkerResult<ProteinPreparationObservations>>?
            PrepareResponse { get; init; }
        public Func<ScientificWorkRequest<ProteinPreparationPayload>, CancellationToken,
            Task<WorkerResult<ProteinPreparationObservations>>>? PrepareResponseAsync { get; init; }
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
            return PrepareResponseAsync is not null
                ? PrepareResponseAsync(request, cancellationToken)
                : Task.FromResult(PrepareResponse!(request));
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
