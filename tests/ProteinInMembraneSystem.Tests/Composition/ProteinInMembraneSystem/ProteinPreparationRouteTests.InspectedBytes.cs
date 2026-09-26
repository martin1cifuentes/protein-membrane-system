using System.Collections.Immutable;
using System.Net;
using System.Text;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed partial class ProteinPreparationRouteTests
{
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public async Task Changed_or_missing_selected_preparation_preview_withholds_approval(bool delete)
    {
        using var directory = new TemporaryDirectory();
        var catalogue = WritePolicyCatalogue(directory.Path);
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(Model()), null)),
            InspectChangesResponse = request =>
            {
                var preview = Path.Combine(request.WorkingDirectory, "selected-preview.pdb");
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

        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes("ATOM\n")),
            "source.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        Assert.True((await Command(product, ActorActionKind.SelectSource, new { uploadToken = token })).Established);
        var modeled = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners = Array.Empty<object>(), alternateLocations = Array.Empty<object>() });
        Assert.True(modeled.Established, modeled.Reason);
        var proposal = Assert.Single(modeled.Value!.Protein!.Changes);
        var selected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = proposal.Id });
        Assert.True(selected.Established, selected.Reason);
        Assert.Contains(selected.Value!.Actions, action =>
            action.Kind == ActorActionKind.ApprovePreparationChange &&
            action.SubjectId == proposal.Id && action.Enabled);

        var previewPath = Assert.Single(Directory.EnumerateFiles(directory.Path,
            "selected-preview.pdb", SearchOption.AllDirectories));
        if (delete) File.Delete(previewPath);
        else File.AppendAllText(previewPath, "changed after inspection");

        Assert.DoesNotContain(product.Snapshot().Actions, action =>
            action.Kind == ActorActionKind.ApprovePreparationChange &&
            action.SubjectId == proposal.Id && action.Enabled);
        var refusal = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = true, rationale = "Approve exact CB" });
        Assert.False(refusal.Established);
        Assert.Contains("structure", refusal.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Equal(0, worker.PrepareProteinCalls);
        Assert.Equal("review", product.Snapshot().Protein!.Status);
        Assert.Contains(product.Snapshot().Actions, action =>
            action.Kind == ActorActionKind.DeclinePreparationChange &&
            action.SubjectId == proposal.Id && action.Enabled);
        var declined = await Command(product, ActorActionKind.ApprovePreparationChange,
            new { proposalId = proposal.Id, approve = false, rationale = "Decline the proposed change" });
        Assert.True(declined.Established, declined.Reason);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public async Task Changed_or_missing_selected_oriented_protein_withholds_placement_adoption(bool delete)
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var proposed = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(proposed.Established, proposed.Reason);
        var placement = Assert.IsType<PlacementAccount>(proposed.Value!.Placement);
        Assert.Equal("supported", placement.Status);
        var selected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = placement.ProposalId });
        Assert.True(selected.Established, selected.Reason);
        Assert.Contains(selected.Value!.Actions, action =>
            action.Kind == ActorActionKind.AdoptPlacement && action.Enabled);

        var orientedPath = Assert.Single(Directory.EnumerateFiles(directory.Path,
            "oriented.pdb", SearchOption.AllDirectories));
        if (delete) File.Delete(orientedPath);
        else File.AppendAllText(orientedPath, "changed after inspection");

        Assert.DoesNotContain(product.Snapshot().Actions, action =>
            action.Kind == ActorActionKind.AdoptPlacement && action.Enabled);
        var refusal = await Command(product, ActorActionKind.AdoptPlacement,
            new { proposalId = placement.ProposalId });
        Assert.False(refusal.Established);
        Assert.Contains("structure", refusal.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Null(product.Snapshot().Study!.AdoptedPlacementProposalId);
    }
}
