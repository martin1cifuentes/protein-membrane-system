using System.Collections.Immutable;
using System.Net;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed partial class ProteinPreparationRouteTests
{
    private const string PlacementSourceText = "ATOM      1  CA  ALA A   1       0.000   0.000  15.000  1.00 50.00           C\n" +
        "ATOM      2  CA  ALA A   2       0.000   0.000   0.000  1.00 50.00           C\n" +
        "ATOM      3  CA  ALA A   3       0.000   0.000 -15.000  1.00 50.00           C\nEND\n";
    private static readonly ImmutableArray<ResidueAddress> PlacementResidues =
        ImmutableArray.Create(new ResidueAddress(0, "A", 1, "", "A"),
            new ResidueAddress(0, "A", 2, "", "A"),
            new ResidueAddress(0, "A", 3, "", "A"));

    [Fact]
    public async Task Exact_source_and_bond_graph_witness_supports_measured_placement_without_a_prepared_coordinate_digest()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var before = product.Snapshot();
        Assert.Contains(before.Actions, action => action.Kind == ActorActionKind.ProposePlacement && action.Enabled);

        var proposed = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(proposed.Established, proposed.Reason);
        var account = Assert.IsType<PlacementAccount>(proposed.Value!.Placement);
        Assert.True(account.Status == "supported", account.Reason + "; witness=" + account.WitnessId +
            "; policy=" + account.PolicyId + "; " + string.Join("; ", account.Limitations) +
            "; evidence=" + string.Join("; ", account.Evidence.Select(item => item.Method + ":" + item.Bearing)));
        var wire = JsonSerializer.SerializeToElement(proposed.Value,
            new JsonSerializerOptions(JsonSerializerDefaults.Web));
        Assert.Equal("supported", wire.GetProperty("placement").GetProperty("status").GetString());
        Assert.Equal("placement-popc", account.PolicyId);
        Assert.Contains("Controlled POPC policy has a contextual source limit.", account.Limitations);
        Assert.Equal("source-topology", account.WitnessId);
        Assert.Equal(proposed.Value.Protein!.SubjectId, account.PreparedProteinId);
        Assert.Equal(proposed.Value.Membrane!.ModelId, account.MembraneModelId);
        Assert.Equal(PlacementPhysicalSide.Both, account.PhysicalSide);
        Assert.Equal(0, account.MidplaneAngstrom);
        Assert.Equal(20, account.ThicknessAngstrom);
        Assert.Contains(account.Evidence, evidence => evidence.Bearing == EvidenceBearing.Supports);
        Assert.Contains(account.ContactingRegions, region => region.Contains("A[A]:2", StringComparison.Ordinal));
        Assert.Single(worker.PpmRequests);
        Assert.Single(worker.MeasurementRequests);
        Assert.Equal(worker.PreparedHash, worker.PpmRequests[0].Payload.PreparedSha256);
        Assert.Equal(worker.OrientedHash, worker.MeasurementRequests[0].Payload.OrientedPdbSha256);

        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = account.ProposalId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal(account.ProposalId, inspected.Value!.Inspection!.SubjectId);
        Assert.NotEmpty(inspected.Value.Inspection.Annotations);
        Assert.NotEmpty(inspected.Value.Inspection.Metrics);
        var adopted = await Command(product, ActorActionKind.AdoptPlacement,
            new { proposalId = account.ProposalId });
        Assert.True(adopted.Established, adopted.Reason);
        Assert.NotEqual(before.Study!.Id, adopted.Value!.Study!.Id);
        Assert.Equal(account.ProposalId, adopted.Value.Study.AdoptedPlacementProposalId);
        Assert.Equal(account.ProposalId, adopted.Value.Placement!.ProposalId);
        Assert.Equal("supported", adopted.Value.Placement!.Status);
        Assert.DoesNotContain(adopted.Value.Actions, action =>
            action.Kind == ActorActionKind.AdoptPlacement && action.Enabled);
        var duplicate = await Command(product, ActorActionKind.AdoptPlacement,
            new { proposalId = account.ProposalId });
        Assert.False(duplicate.Established);
        Assert.Equal(adopted.Value.Study.Id, product.Snapshot().Study!.Id);
    }

    [Fact]
    public async Task Missing_measurement_retains_inspectable_position_with_not_established_reason()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker { FailMeasurement = true };
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);

        var proposed = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(proposed.Established, proposed.Reason);
        var account = Assert.IsType<PlacementAccount>(proposed.Value!.Placement);
        Assert.Equal("notEstablished", account.Status);
        Assert.Contains("measurement", account.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Equal(0, account.MidplaneAngstrom);
        Assert.Contains(account.Evidence, evidence => evidence.Method == "PPM orientation candidate");
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = account.ProposalId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.Equal(account.ProposalId, inspected.Value!.Inspection!.SubjectId);
        Assert.NotEmpty(inspected.Value.Inspection.Annotations);
        Assert.NotEmpty(inspected.Value.Inspection.Metrics);
        Assert.DoesNotContain(inspected.Value.Actions, action =>
            action.Kind == ActorActionKind.AdoptPlacement && action.Enabled);
    }

    [Fact]
    public async Task Chosen_unassessed_membrane_still_allows_an_inspectable_position_without_support()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var replacement = await Command(product, ActorActionKind.ProposeMembrane,
            new { upper = new[] { new { speciesId = "UNKNOWN", fraction = 1.0 } },
                lower = new[] { new { speciesId = "UNKNOWN", fraction = 1.0 } },
                scientificPurpose = "Unqualified explicit membrane choice" });
        var chosen = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = replacement.Value!.Membrane!.ModelId });
        Assert.Equal("notEstablished", chosen.Value!.Membrane!.Status);
        Assert.Equal("assessed", chosen.Value.Protein!.Status);

        var proposed = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal("notEstablished", proposed.Value!.Placement!.Status);
        Assert.Equal(chosen.Value.Membrane.ModelId, proposed.Value.Placement.MembraneModelId);
        Assert.Empty(worker.MeasurementRequests);
        var inspected = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = proposed.Value.Placement.ProposalId });
        Assert.True(inspected.Established, inspected.Reason);
        Assert.NotEmpty(inspected.Value!.Inspection!.Metrics);
    }

    [Fact]
    public async Task Failed_local_orientation_without_an_alternate_has_no_proposal_and_no_scientific_unsupported_verdict()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker { FailPpm = true };
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);

        var attempted = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(attempted.Established, attempted.Reason);
        Assert.Null(attempted.Value!.Placement);
        Assert.Single(worker.PpmRequests);
        Assert.Contains(attempted.Value.Notices, notice =>
            notice.Message.Contains("orientation failed", StringComparison.OrdinalIgnoreCase));
    }

    [Theory]
    [InlineData(0, "No placement support policy")]
    [InlineData(2, "ambiguous placement support policies")]
    public async Task Absent_or_ambiguous_placement_policy_cannot_promote_a_measured_candidate(
        int policyCopies, string expectedReason)
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker, policyCopies);
        await EstablishProteinAndMembrane(product);
        var proposed = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal("notEstablished", proposed.Value!.Placement!.Status);
        Assert.Contains(expectedReason, proposed.Value.Placement.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Single(worker.MeasurementRequests);
        Assert.DoesNotContain(proposed.Value.Actions, action =>
            action.Kind == ActorActionKind.AdoptPlacement && action.Enabled);
    }

    [Fact]
    public async Task Correcting_a_supported_position_withdraws_its_conclusion_until_new_measurement()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker { FailSecondMeasurement = true };
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var original = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.Equal("supported", original.Value!.Placement!.Status);
        var priorId = original.Value.Placement.ProposalId;

        var revised = await Command(product, ActorActionKind.RevisePlacement, new
        {
            proposalId = priorId, depthShiftAngstrom = 1.0, tiltAboutXDegrees = 0.0,
            tiltAboutYDegrees = 0.0, rotationAboutNormalDegrees = 0.0,
            rationale = "Reobserve the upper interface after a bounded shift."
        });
        Assert.True(revised.Established, revised.Reason);
        Assert.NotEqual(priorId, revised.Value!.Placement!.ProposalId);
        Assert.Equal("notEstablished", revised.Value.Placement.Status);
        Assert.Null(revised.Value.Placement.TiltDegrees);
        Assert.Equal("supported", original.Value.Placement.Status);
        Assert.Equal(2, worker.MeasurementRequests.Count);
        Assert.Single(worker.AdjustmentRequests);
        var old = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = priorId });
        Assert.False(old.Established);
    }

    [Fact]
    public async Task Changed_ppm_residue_library_disables_local_route_before_worker_invocation()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        File.AppendAllText(Path.Combine(directory.Path, "res.lib"), "changed");
        Assert.DoesNotContain(product.Snapshot().Actions, action =>
            action.Kind == ActorActionKind.ProposePlacement && action.Enabled);

        var refused = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.True(refused.Established, refused.Reason);
        Assert.Null(refused.Value!.Placement);
        Assert.Empty(worker.PpmRequests);
        Assert.Contains(refused.Value.Notices, notice => notice.Message.Contains("residue library", StringComparison.Ordinal));
    }

    [Fact]
    public async Task Qualified_exact_coordinate_opm_reference_can_propose_without_invoking_ppm()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var catalogue = WritePlacementCatalogue(directory.Path);
        var ppm = Path.Combine(directory.Path, "immers");
        File.AppendAllText(ppm, "unavailable");
        using var http = new HttpClient(new RespondingHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        { Content = new ByteArrayContent(Encoding.ASCII.GetBytes(PlacementSourceText)) }));
        var sources = new ControlledOpmExchange(http, () =>
            new OpmReferenceRecord("6QWR", "controlled OPM reference", worker.PreparedPath,
                worker.PreparedHash, null, ImmutableArray.Create(new ChainSelection("A", "A")),
                ImmutableArray<string>.Empty, ImmutableArray.Create("source-1", "source-2", "source-3"),
                "implicit symmetric POPC", 20, 8, "POPC", true, 0));
        var product = new ProductRoot(worker, sources, Path.Combine(directory.Path, "workspace"),
            catalogue, ppm, () => null);
        var selected = await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "6QWR" });
        Assert.True(selected.Established, selected.Reason);
        await EstablishSelectedProteinAndMembrane(product);
        Assert.Contains(product.Snapshot().Actions, action =>
            action.Kind == ActorActionKind.ProposePlacement && action.Enabled);

        var proposed = await Command(product, ActorActionKind.ProposePlacement, new
        {
            orientationRoute = "opm", topologyKind = "membrane-spanning", physicalSide = "both",
            biologicalSidedness = "extracellular upper"
        });
        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal("supported", proposed.Value!.Placement!.Status);
        Assert.Contains(proposed.Value.Placement.Evidence,
            evidence => evidence.Method == "OPM exact-coordinate orientation candidate" &&
                evidence.Bearing == EvidenceBearing.Context);
        Assert.Empty(worker.PpmRequests);
        Assert.Single(worker.MeasurementRequests);
    }

    [Fact]
    public async Task Contextual_opm_reference_and_fresh_ppm_position_keep_distinct_provenance()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var catalogue = WritePlacementCatalogue(directory.Path);
        var opmPath = Path.Combine(directory.Path, "opm-reference.pdb");
        File.WriteAllText(opmPath, PlacementSourceText + "REMARK independently oriented context\n");
        using var http = new HttpClient(new RespondingHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        { Content = new ByteArrayContent(Encoding.ASCII.GetBytes(PlacementSourceText)) }));
        var sources = new ControlledOpmExchange(http, () =>
            new OpmReferenceRecord("6QWR", "controlled OPM reference", opmPath,
                Hash(opmPath), null, ImmutableArray.Create(new ChainSelection("A", "A")),
                ImmutableArray<string>.Empty, ImmutableArray.Create("source-1", "source-2", "source-3"),
                "implicit symmetric POPC", 20, 8, "POPC", true, 0));
        var product = new ProductRoot(worker, sources, Path.Combine(directory.Path, "workspace"),
            catalogue, Path.Combine(directory.Path, "immers"), () => null);
        var selected = await Command(product, ActorActionKind.SelectSource,
            new { sourceKind = "rcsb", exactIdentifier = "6QWR" });
        Assert.True(selected.Established, selected.Reason);
        await EstablishSelectedProteinAndMembrane(product);

        var proposed = await Command(product, ActorActionKind.ProposePlacement, new
        {
            orientationRoute = "auto", topologyKind = "membrane-spanning", physicalSide = "both",
            biologicalSidedness = "extracellular upper", ppmNterminalSide = "out"
        });
        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal("supported", proposed.Value!.Placement!.Status);
        Assert.Contains(proposed.Value.Placement.Evidence,
            evidence => evidence.Method == "OPM oriented-structure reference" &&
                evidence.Bearing == EvidenceBearing.Context);
        Assert.Contains(proposed.Value.Placement.Evidence,
            evidence => evidence.Method == "PPM orientation candidate" &&
                evidence.Bearing == EvidenceBearing.Context);
        Assert.Single(worker.PpmRequests);
        Assert.Single(worker.MeasurementRequests);
    }

    [Fact]
    public async Task Consequential_membrane_change_withdraws_current_placement_without_rewriting_earlier_account()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var positioned = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.Equal("supported", positioned.Value!.Placement!.Status);
        var priorProposalId = positioned.Value.Placement.ProposalId;
        var replacement = await Command(product, ActorActionKind.ProposeMembrane,
            new { upper = new[] { new { speciesId = "UNKNOWN", fraction = 1.0 } },
                lower = new[] { new { speciesId = "UNKNOWN", fraction = 1.0 } },
                scientificPurpose = "A new exact bilayer premise" });
        Assert.True(replacement.Established, replacement.Reason);
        var changed = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = replacement.Value!.Membrane!.ModelId });
        Assert.True(changed.Established, changed.Reason);
        Assert.Null(changed.Value!.Placement);
        Assert.Equal("notEstablished", changed.Value.Membrane!.Status);
        Assert.Equal("supported", positioned.Value.Placement.Status);
        var stale = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = priorProposalId });
        Assert.False(stale.Established);
        Assert.Empty(worker.MeasurementRequests.Skip(1));
    }

    [Fact]
    public async Task New_exact_protein_source_withdraws_the_prior_placement_reliance()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var positioned = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.Equal("supported", positioned.Value!.Placement!.Status);
        var formerProposalId = positioned.Value.Placement.ProposalId;
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes(
            PlacementSourceText + "REMARK changed exact source\n")), "changed.pdb",
            UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);

        var changed = await Command(product, ActorActionKind.SelectSource, new { uploadToken = token });
        Assert.True(changed.Established, changed.Reason);
        Assert.Null(changed.Value!.Placement);
        Assert.NotEqual(positioned.Value.Study!.Id, changed.Value.Study!.Id);
        Assert.Equal("supported", positioned.Value.Placement.Status);
        var oldInspection = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = formerProposalId });
        Assert.False(oldInspection.Established);
    }

    [Fact]
    public async Task Changed_retained_partner_decision_withdraws_the_prior_placement_reliance()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker { IncludePartner = true };
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product, includePartner: true);
        var positioned = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.Equal("supported", positioned.Value!.Placement!.Status);
        Assert.False(positioned.Value.Study!.Partners.Single().Retain);
        var formerProposalId = positioned.Value.Placement.ProposalId;

        var changed = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners = new[] { new { sourceId = "partner-X", retain = true,
                    reason = "Retain the exact observed cofactor for this revision." } },
                alternateLocations = Array.Empty<object>() });
        Assert.True(changed.Established, changed.Reason);
        Assert.Null(changed.Value!.Placement);
        Assert.True(changed.Value.Study!.Partners.Single().Retain);
        Assert.NotEqual(positioned.Value.Study.Id, changed.Value.Study.Id);
        var oldInspection = await Command(product, ActorActionKind.SelectInspectionSubject,
            new { subjectId = formerProposalId });
        Assert.False(oldInspection.Established);
    }

    [Fact]
    public async Task Changed_biological_sidedness_creates_a_new_unwitnessed_proposal()
    {
        using var directory = new TemporaryDirectory();
        var worker = new PlacementRouteWorker();
        var product = PlacementProduct(directory.Path, worker);
        await EstablishProteinAndMembrane(product);
        var positioned = await Command(product, ActorActionKind.ProposePlacement, PlacementChoice());
        Assert.Equal("supported", positioned.Value!.Placement!.Status);
        var priorId = positioned.Value.Placement.ProposalId;

        var changed = await Command(product, ActorActionKind.ProposePlacement,
            new { orientationRoute = "ppm", topologyKind = "membrane-spanning", physicalSide = "both",
                biologicalSidedness = "periplasmic upper", ppmNterminalSide = "out" });
        Assert.True(changed.Established, changed.Reason);
        Assert.NotEqual(priorId, changed.Value!.Placement!.ProposalId);
        Assert.Equal("notEstablished", changed.Value.Placement.Status);
        Assert.Null(changed.Value.Placement.WitnessId);
        Assert.Equal("supported", positioned.Value.Placement.Status);
        Assert.Equal(2, worker.PpmRequests.Count);
        Assert.Equal(2, worker.MeasurementRequests.Count);
    }

    private static object PlacementChoice() => new
    {
        orientationRoute = "ppm", topologyKind = "membrane-spanning", physicalSide = "both",
        biologicalSidedness = "extracellular upper", ppmNterminalSide = "out"
    };

    private static ProductRoot PlacementProduct(string directory, PlacementRouteWorker worker,
        int placementPolicyCopies = 1)
    {
        var catalogue = WritePlacementCatalogue(directory, placementPolicyCopies);
        using var http = new HttpClient(new NoNetworkHandler());
        return new ProductRoot(worker, new ExternalSourceExchange(http),
            Path.Combine(directory, "workspace"), catalogue,
            Path.Combine(directory, "immers"), () => null);
    }

    private static async Task EstablishProteinAndMembrane(ProductRoot product, bool includePartner = false)
    {
        var token = await product.UploadAsync(new MemoryStream(Encoding.ASCII.GetBytes(PlacementSourceText)),
            "alkl.pdb", UploadOriginKind.Experimental, null, TestContext.Current.CancellationToken);
        var source = await Command(product, ActorActionKind.SelectSource, new { uploadToken = token });
        Assert.True(source.Established, source.Reason);
        await EstablishSelectedProteinAndMembrane(product, includePartner);
    }

    private static async Task EstablishSelectedProteinAndMembrane(ProductRoot product,
        bool includePartner = false)
    {
        object[] partners = includePartner
            ? [new { sourceId = "partner-X", retain = false,
                reason = "Exclude the observed cofactor from the selected construct." }]
            : [];
        var protein = await Command(product, ActorActionKind.SelectProteinModel,
            new { modelIndex = 0, biologicalAssemblyId = (string?)null,
                chains = new[] { new { sourceChain = "A", copyId = "A" } },
                partners, alternateLocations = Array.Empty<object>() });
        Assert.True(protein.Established, protein.Reason);
        Assert.Equal("assessed", protein.Value!.Protein!.Status);
        var proposed = await Command(product, ActorActionKind.ProposeMembrane,
            new { upper = new[] { new { speciesId = "POPC", fraction = 1.0 } },
                lower = new[] { new { speciesId = "POPC", fraction = 1.0 } },
                scientificPurpose = "Controlled exact POPC bilayer" });
        Assert.True(proposed.Established, proposed.Reason);
        var adopted = await Command(product, ActorActionKind.AdoptMembrane,
            new { modelId = proposed.Value!.Membrane!.ModelId });
        Assert.True(adopted.Established, adopted.Reason);
        Assert.Equal("assessed", adopted.Value!.Membrane!.Status);
        Assert.Equal("assessed", adopted.Value.Protein!.Status);
    }

    private static string WritePlacementCatalogue(string directory, int placementPolicyCopies = 1)
    {
        var ff = Path.Combine(directory, "protein-forcefield.xml");
        var lipidFf = Path.Combine(directory, "lipid-forcefield.xml");
        var lipidCoordinates = Path.Combine(directory, "popc.pdb");
        var ppm = Path.Combine(directory, "immers");
        var residueLibrary = Path.Combine(directory, "res.lib");
        File.WriteAllText(ff, "<ForceField/>");
        File.WriteAllText(lipidFf, "<ForceField/>");
        File.WriteAllText(lipidCoordinates, "POPC");
        File.WriteAllText(ppm, "identified local executable fixture");
        File.WriteAllText(residueLibrary, "identified local residue library fixture");
        var chemical = new ProteinChemicalStatePolicy("chemical", "1", "canonical-amino-acid-assembly",
            2.5, ImmutableArray.Create("identified protein chemistry"),
            ImmutableArray.Create(new ForceFieldAsset("ff", "1", "Amber19", ff, Hash(ff))),
            ImmutableDictionary<string, string>.Empty, ImmutableArray.Create("HID"), ImmutableArray<string>.Empty);
        var lipid = new MolecularRepresentation("POPC", "POPC-chemistry", "lipid", lipidFf,
            Hash(lipidFf), lipidCoordinates, Hash(lipidCoordinates), 3, 0, 68, 1200,
            ImmutableArray.Create(1), "Lipid21", "8.6", ImmutableArray<string>.Empty,
            ImmutableArray.Create(new MolecularStereoCheck("tetrahedral",
                ImmutableArray.Create("O21", "C1", "C3", "HS"), "negative")));
        var membranePolicy = new MembraneSupportPolicy("membrane-popc", "1",
            ImmutableArray.Create("POPC"), false, false, true,
            ImmutableArray.Create("identified membrane policy"), ImmutableArray<string>.Empty);
        var placementPolicy = new PlacementSupportPolicy("placement-popc", "1",
            ImmutableArray.Create(ProteinTopologyKind.MembraneSpanning), ImmutableArray.Create("POPC"),
            false, false, false,
            ImmutableArray.Create("Independently witnessed placement relationship"),
            ImmutableArray.Create(new PlacementMeasurementCriterion("atomsWithinCore", "atoms", 1, null)),
            3, null, ImmutableArray.Create("independent experimental topology"),
            ImmutableArray.Create("Controlled POPC policy has a contextual source limit."));
        var upper = new LeafletComposition(LeafletSide.Upper,
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var lower = new LeafletComposition(LeafletSide.Lower,
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var witness = new PlacementStructuralWitness("source-topology", "1",
            Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(PlacementSourceText))).ToLowerInvariant(),
            0, null, ImmutableArray.Create(new ChainSelection("A", "A")),
            ImmutableArray<string>.Empty, null,
            Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes("{}"))).ToLowerInvariant(),
            upper, lower, FixedStudyConditions.Initial, ProteinTopologyKind.MembraneSpanning,
            "extracellular upper", "independent experimental topology",
            ImmutableArray.Create("independent experimental topology"),
            ImmutableArray.Create(
                new PlacementResidueWitness(PlacementResidues[0], PlacementWitnessRoleKind.Topology,
                    "upper-water", "independent experimental topology"),
                new PlacementResidueWitness(PlacementResidues[2], PlacementWitnessRoleKind.Topology,
                    "lower-water", "independent experimental topology"),
                new PlacementResidueWitness(PlacementResidues[1], PlacementWitnessRoleKind.Contact,
                    "core-contact", "independent experimental topology"),
                new PlacementResidueWitness(PlacementResidues[0], PlacementWitnessRoleKind.Sidedness,
                    "upper-water", "independent experimental topology")), ImmutableArray<string>.Empty,
            ImmutableArray<PlacementPredictionResolution>.Empty);
        var catalogue = Path.Combine(directory, "catalogue.json");
        File.WriteAllText(catalogue, JsonSerializer.Serialize(new
        {
            version = "1", evidenceReferences = new[] { "controlled catalogue evidence" },
            lipids = new[] { lipid }, proteinChemicalStates = new[] { chemical },
            proteinStructuralPolicies = new[] { StructuralPolicy() },
            membranePolicies = new[] { membranePolicy },
            placementPolicies = Enumerable.Range(0, placementPolicyCopies).Select(index => index == 0
                ? placementPolicy : placementPolicy with { Id = $"placement-popc-{index}" }).ToArray(),
            placementWitnesses = new[] { witness }, preparationPolicies = Array.Empty<object>(),
            equilibrationQualifications = Array.Empty<object>(), ppmVersion = "2.0",
            ppmExecutableSha256 = Hash(ppm), ppmResidueLibraryPath = "res.lib",
            ppmResidueLibrarySha256 = Hash(residueLibrary),
            maximumSourceAtoms = 1000
        }, new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        return catalogue;
    }

    private sealed class PlacementRouteWorker : IScientificWorkerExchange
    {
        public bool IncludePartner { get; init; }
        public bool FailMeasurement { get; init; }
        public bool FailSecondMeasurement { get; init; }
        public bool FailPpm { get; init; }
        public string PreparedHash { get; private set; } = "";
        public string PreparedPath { get; private set; } = "";
        public string OrientedHash { get; private set; } = "";
        public List<ScientificWorkRequest<PlacementPayload>> PpmRequests { get; } = [];
        public List<ScientificWorkRequest<PlacementMeasurementPayload>> MeasurementRequests { get; } = [];
        public List<ScientificWorkRequest<PlacementAdjustmentPayload>> AdjustmentRequests { get; } = [];

        public Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
            ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken)
        {
            var sourceResidues = PlacementResidues.Select(address => new SourceResidueObservation(
                address with { CopyId = "" }, "ALA", SourceResidueKind.Protein, true,
                ImmutableArray<string>.Empty, ImmutableArray<string>.Empty, null,
                ObservationStanding.Unavailable, ObservationStanding.Unavailable,
                ImmutableArray<string>.Empty)).ToImmutableArray();
            var model = new SourceModelObservation(0,
                ImmutableArray.Create(new SourceChainObservation("A", 3, 3)),
                ImmutableArray<SourceAssemblyObservation>.Empty,
                IncludePartner
                    ? ImmutableArray.Create(new SourcePartnerObservation("partner-X", "Observed cofactor",
                        "cofactor", 1, "A", 4))
                    : ImmutableArray<SourcePartnerObservation>.Empty,
                sourceResidues, 3);
            return Task.FromResult(Observed(request.RequestId, null,
                new SourceInspectionObservations("pdb", ImmutableArray.Create(model), null)));
        }

        public Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
            ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken)
        {
            var preview = Path.Combine(request.WorkingDirectory, "selected.pdb");
            Directory.CreateDirectory(request.WorkingDirectory);
            File.WriteAllText(preview, PlacementSourceText);
            return Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId,
                new PreparationChangeObservations(ImmutableArray<AtomAddress>.Empty,
                    ImmutableArray<PossibleDisulfideObservation>.Empty,
                    ObservationStanding.Observed, ImmutableArray<string>.Empty, 3,
                    ImmutableArray.Create(new PreviewChainCorrespondence("A", "A", "A")),
                    ObservedGeometry()), ImmutableArray.Create(Artifact(preview, "selectedProteinPreview"))));
        }

        public Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
            ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken)
        {
            Directory.CreateDirectory(request.WorkingDirectory);
            var prepared = Path.Combine(request.WorkingDirectory, "prepared.pdb");
            PreparedPath = prepared;
            var graph = Path.Combine(request.WorkingDirectory, "graph.json");
            var mapping = Path.Combine(request.WorkingDirectory, "correspondence.json");
            File.WriteAllText(prepared, PlacementSourceText);
            File.WriteAllText(graph, "{}");
            PreparedHash = Hash(prepared);
            var correspondence = new SourceToResultCorrespondence(request.Payload.SourceSha256,
                PreparedHash, PlacementResidues.Select((address, index) => new AtomCorrespondence(index,
                    $"{index}:A:{address.Residue}::CA", $"source-{index + 1}", AtomOriginKind.Source,
                    MoleculeRoleKind.Protein, AtomRoleKind.Backbone, "C", address, null)).ToImmutableArray(), true);
            File.WriteAllText(mapping, JsonSerializer.Serialize(correspondence,
                new JsonSerializerOptions(JsonSerializerDefaults.Web)));
            var observation = new ProteinPreparationObservations(3, 3, 3,
                ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
                ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
                ImmutableArray<AtomAddress>.Empty, ImmutableArray<AtomAddress>.Empty,
                ImmutableArray<AtomAddress>.Empty, ImmutableArray<ResidueVariantChoice>.Empty,
                ImmutableArray<DisulfideBond>.Empty, ImmutableArray<string>.Empty, 3,
                ObservedGeometry());
            return Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId, observation,
                ImmutableArray.Create(Artifact(prepared, "preparedPdb"),
                    Artifact(graph, "preparedBondGraph"), Artifact(mapping, "correspondenceJson"))));
        }

        public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
            ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken) =>
            Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId,
                new MembraneAssessmentObservations(ImmutableArray.Create(
                    new SpeciesTemplateObservation("POPC", "POPC-chemistry", 3, 3, true,
                        ImmutableArray<string>.Empty)), ImmutableArray<string>.Empty, true)));

        public Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
            ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken)
        {
            PpmRequests.Add(request);
            if (FailPpm)
                return Task.FromResult(new WorkerResult<PlacementObservations>(request.RequestId,
                    request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
                    ImmutableArray<WorkerArtifact>.Empty, null, null, "orientationFailed",
                    "The local orientation failed without a corresponding output."));
            Directory.CreateDirectory(request.WorkingDirectory);
            var oriented = Path.Combine(request.WorkingDirectory, "oriented.pdb");
            File.Copy(request.Payload.PreparedPdbPath, oriented);
            OrientedHash = Hash(oriented);
            var observations = new PlacementObservations(0, 20, 8, 3,
                ImmutableArray<MeasuredValue>.Empty, ImmutableArray.Create("N:1", "O:1"),
                ImmutableArray<string>.Empty, "PPM 2.0 implicit symmetric DOPC membrane");
            return Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId,
                observations, ImmutableArray.Create(Artifact(oriented, "orientedPdb"))));
        }

        public Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
            ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken)
        {
            MeasurementRequests.Add(request);
            if (FailMeasurement || FailSecondMeasurement && MeasurementRequests.Count > 1)
                return Task.FromResult(new WorkerResult<PlacementMeasurementObservations>(request.RequestId,
                    request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
                    ImmutableArray<WorkerArtifact>.Empty, null, null, "measurementFailed",
                    "The exact placement measurement failed."));
            var observations = PlacementResidues.Select((address, index) => new PlacementResidueObservation(
                address, "A", address.Residue.ToString(), "", "ALA", 1,
                index == 0 ? 15 : index == 1 ? 0 : -15,
                index == 0 ? 15 : index == 1 ? 0 : -15,
                index == 0 ? 15 : index == 1 ? 0 : -15,
                index == 1 ? 1 : 0, index == 1 ? 1 : 0,
                index == 0 ? 1 : 0, index == 2 ? 1 : 0)).ToImmutableArray();
            return Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId,
                new PlacementMeasurementObservations(3, observations, 1, 1, 1, -15, 15,
                    ImmutableArray<string>.Empty)));
        }

        public Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
            ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken) =>
            NotUsed<PredictionRegionSummaryObservations>();
        public Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
            ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken)
        {
            AdjustmentRequests.Add(request);
            Directory.CreateDirectory(request.WorkingDirectory);
            var adjusted = Path.Combine(request.WorkingDirectory, "adjusted.pdb");
            File.Copy(request.Payload.OrientedPdbPath, adjusted);
            return Task.FromResult(Observed(request.RequestId, request.Payload.StudyRevisionId,
                new PlacementAdjustmentObservations(3, 3, 1, 0, 0, 0,
                    ImmutableArray<string>.Empty),
                ImmutableArray.Create(Artifact(adjusted, "adjustedPdb"))));
        }
        public Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
            ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken) =>
            NotUsed<ConstructionObservations>();
        public Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
            ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken) =>
            NotUsed<MinimizationObservations>();
        public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
            ScientificWorkRequest<EquilibrationPayload> request,
            IProgress<EquilibrationWorkProgress>? progress, CancellationToken cancellationToken) =>
            NotUsed<EquilibrationObservations>();
        public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
            ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken) =>
            NotUsed<StageObservationObservations>();
        public Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
            ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken) =>
            NotUsed<ExportVerificationObservations>();
        private static Task<WorkerResult<T>> NotUsed<T>() where T : class =>
            throw new InvalidOperationException("An unrelated scientific operation was invoked.");
    }

    private sealed class ControlledOpmExchange(HttpClient http, Func<OpmReferenceRecord> reference)
        : ExternalSourceExchange(http)
    {
        public override Task<OpmReferenceRecord?> TryRetrieveOpmReferenceAsync(string pdbAccession,
            string targetDirectory, CancellationToken cancellationToken) =>
            Task.FromResult<OpmReferenceRecord?>(reference());
    }
}
