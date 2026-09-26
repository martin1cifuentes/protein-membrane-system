using System.Collections.Immutable;
using System.Security.Cryptography;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using PlacementOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.PlacementAssessment.PlacementAssessment;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class PlacementAssessmentOwnerTests
{
    [Fact]
    public void Opm_context_is_corresponding_only_for_the_exact_identified_construct_and_membrane()
    {
        using var artifact = new IdentifiedArtifact();
        using var fixture = new PlacementFixture();
        var owner = new PlacementOwner(new PlacementWorkerStub());
        var reference = new OpmReferenceRecord("6QWR", "https://opm.example/6qwr", artifact.Path,
            artifact.Sha256, "assembly-1", fixture.Protein.Intended.Chains,
            ImmutableArray<string>.Empty, ImmutableArray.Create("atom-1", "atom-2", "atom-3"),
            "symmetric POPC", 22, 8, "POPC", true, 0);

        var matching = owner.ReviewOpmReference(fixture.Revision, fixture.Protein, fixture.Membrane.Intended,
            reference with { MembraneContext = "symmetric POPC" });
        Assert.True(matching.Established, matching.Reason);
        Assert.True(matching.Value!.CorrespondsToSelectedConstruct);
        Assert.True(matching.Value.MembraneContextApplicable);
        Assert.All(matching.Value.Evidence, evidence => Assert.Equal(EvidenceBearing.Context, evidence.Bearing));

        var mismatches = new (OpmReferenceRecord Reference, string Expected)[]
        {
            (reference with { PdbAccession = "OTHER" }, "accession"),
            (reference with { BiologicalAssemblyId = "assembly-2" }, "assembly"),
            (reference with { ChainCopies = ImmutableArray.Create(new ChainSelection("B", "other")) }, "chain"),
            (reference with { RetainedPartnerSourceIds = ImmutableArray.Create("partner-x") }, "partners"),
            (reference with { SourceAtomIds = ImmutableArray.Create("atom-1", "atom-2", "other") }, "atoms"),
            (reference with { OrientedCoordinateSha256 = new string('0', 64) }, "hash"),
            (reference with { OrientedCoordinatePath = artifact.Path + ".absent" }, "artifact")
        };
        foreach (var (changed, expected) in mismatches)
        {
            var reviewed = owner.ReviewOpmReference(fixture.Revision, fixture.Protein,
                fixture.Membrane.Intended, changed);
            Assert.True(reviewed.Established, reviewed.Reason);
            Assert.False(reviewed.Value!.CorrespondsToSelectedConstruct);
            Assert.Contains(reviewed.Value.Limitations, limitation =>
                limitation.Contains(expected, StringComparison.OrdinalIgnoreCase));
        }
        foreach (var changed in new[]
                 {
                     reference with { AssumedMembraneSpeciesId = "DOPC" },
                     reference with { AssumedMembraneSpeciesId = null },
                     reference with { ImplicitSymmetric = false }
                 })
        {
            var reviewed = owner.ReviewOpmReference(fixture.Revision, fixture.Protein,
                fixture.Membrane.Intended, changed);
            Assert.True(reviewed.Established, reviewed.Reason);
            Assert.True(reviewed.Value!.CorrespondsToSelectedConstruct);
            Assert.False(reviewed.Value.MembraneContextApplicable);
            Assert.Contains(reviewed.Value.Limitations, limitation =>
                limitation.Contains("membrane", StringComparison.OrdinalIgnoreCase));
        }
    }

    [Fact]
    public async Task Ppm_candidate_requires_correlated_complete_hashed_oriented_observation()
    {
        using var artifact = new IdentifiedArtifact();
        using var fixture = new PlacementFixture();
        var ownerWorker = new PlacementWorkerStub
        {
            PpmReply = request => Observed(request, artifact)
        };
        var owner = new PlacementOwner(ownerWorker);
        var valid = await Propose(owner, fixture);
        Assert.True(valid.Established, valid.Reason);
        Assert.Equal(fixture.Protein.Id, valid.Value!.PreparedProteinId);
        Assert.Equal(fixture.Membrane.Intended.Id, valid.Value.MembraneModelId);
        Assert.Equal(3, valid.Value.OrientedProtein.AtomCount);
        Assert.Equal(PlacementPhysicalSide.Both, valid.Value.PhysicalSide);
        Assert.Contains("symmetric DOPC", valid.Value.Evidence.Single().Applicability,
            StringComparison.OrdinalIgnoreCase);
        Assert.Single(ownerWorker.PpmRequests);
        Assert.Equal(fixture.Protein.Molecule.CoordinateSha256,
            ownerWorker.PpmRequests[0].Payload.PreparedSha256);

        var replies = new (string Name, Func<ScientificWorkRequest<PlacementPayload>, WorkerResult<PlacementObservations>> Reply)[]
        {
            ("request", request => Observed(request, artifact) with { RequestId = "wrong" }),
            ("revision", request => Observed(request, artifact) with { StudyRevisionId = "wrong" }),
            ("failed", request => Observed(request, artifact) with { Standing = WorkerResultStanding.Failed }),
            ("stopped", request => Observed(request, artifact) with { Standing = WorkerResultStanding.Stopped }),
            ("observations", request => Observed(request, artifact) with { Observations = null }),
            ("artifact", request => Observed(request, artifact) with { Artifacts = ImmutableArray<WorkerArtifact>.Empty }),
            ("artifact bytes", request => Observed(request, artifact) with { Artifacts = ImmutableArray.Create(
                new WorkerArtifact("orientedPdb", artifact.Path, new string('0', 64))) }),
            ("atom count", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { AlignedSourceAtomCount = 2 } }),
            ("midplane", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { MidplaneAngstrom = double.NaN } }),
            ("tilt", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { TiltDegrees = double.NaN } }),
            ("thickness", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { ThicknessAngstrom = 0 } }),
            ("implicit context", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { AssumedMembrane = "unknown" } }),
            ("paired planes", request => Observed(request, artifact) with { Observations =
                Observed(request, artifact).Observations! with { PlaneMarkerIds = ImmutableArray.Create("N:1") } })
        };
        foreach (var (name, reply) in replies)
        {
            var worker = new PlacementWorkerStub { PpmReply = reply };
            var result = await Propose(new PlacementOwner(worker), fixture);
            Assert.False(result.Established, name);
            Assert.Single(worker.PpmRequests);
        }
        var invalidSide = await owner.ProposeWithPpmAsync(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, ProteinTopologyKind.MembraneSpanning, PlacementPhysicalSide.Upper,
            "extracellular upper", PpmNterminalSide.Out, "/identified/immers", "2.0",
            new string('a', 64), "/tmp/placement-unit", TestContext.Current.CancellationToken);
        Assert.False(invalidSide.Established);
        Assert.Single(ownerWorker.PpmRequests);
        var unidentified = await owner.ProposeWithPpmAsync(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, ProteinTopologyKind.MembraneSpanning, PlacementPhysicalSide.Both,
            "extracellular upper", PpmNterminalSide.Out, "", "2.0", new string('a', 64),
            "/tmp/placement-unit", TestContext.Current.CancellationToken,
            "/identified/res.lib", new string('f', 64));
        Assert.False(unidentified.Established);
        Assert.Single(ownerWorker.PpmRequests);
    }

    [Fact]
    public void Opm_only_proposal_needs_the_exact_prepared_bytes_and_observed_frame_bounds()
    {
        using var fixture = new PlacementFixture();
        using var artifact = new IdentifiedArtifact();
        fixture.Protein = fixture.Protein with { Molecule = fixture.Protein.Molecule with
        { CoordinatePath = artifact.Path, CoordinateSha256 = artifact.Sha256 } };
        var reference = new OpmReferenceRecord("6QWR", "controlled exact OPM reference",
            artifact.Path, artifact.Sha256, "assembly-1", fixture.Protein.Intended.Chains,
            ImmutableArray<string>.Empty, ImmutableArray.Create("atom-1", "atom-2", "atom-3"),
            "implicit symmetric POPC", 20, 8, "POPC", true, 0);
        var owner = new PlacementOwner(new PlacementWorkerStub());
        var review = owner.ReviewOpmReference(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, reference);
        Assert.True(review.Value!.CorrespondsToSelectedConstruct);
        Assert.True(review.Value.MembraneContextApplicable);
        var proposed = owner.ProposeFromOpmReference(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, reference, review.Value,
            ProteinTopologyKind.MembraneSpanning, PlacementPhysicalSide.Both,
            "extracellular upper");
        Assert.True(proposed.Established, proposed.Reason);
        Assert.Equal(artifact.Sha256, proposed.Value!.OrientedProtein.CoordinateSha256);
        Assert.Equal(0, proposed.Value.MidplaneAngstrom);
        Assert.Equal(20, proposed.Value.ThicknessAngstrom);
        Assert.All(proposed.Value.Evidence, item => Assert.Equal(EvidenceBearing.Context, item.Bearing));
        foreach (var invalid in new[]
                 {
                     reference with { MidplaneAngstrom = null },
                     reference with { HydrophobicThicknessAngstrom = double.NaN },
                     reference with { OrientedCoordinateSha256 = new string('0', 64) },
                     reference with { AssumedMembraneSpeciesId = "DOPC" }
                 })
            Assert.False(owner.ProposeFromOpmReference(fixture.Revision, fixture.Protein,
                fixture.Membrane.Intended, invalid, review.Value,
                ProteinTopologyKind.MembraneSpanning, PlacementPhysicalSide.Both,
                "extracellular upper").Established);
    }

    [Fact]
    public async Task Independent_witness_and_policy_can_support_exact_measured_spanning_relationship()
    {
        using var fixture = new PlacementFixture();
        var worker = new PlacementWorkerStub { MeasurementReply = request => Measured(request, fixture) };
        var owner = new PlacementOwner(worker);
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.Equal(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        var result = Assess(owner, fixture, measured.Value);
        Assert.Equal(AssessmentStanding.Supported, result.Standing);
        Assert.Single(worker.MeasurementRequests);
        Assert.Equal(fixture.Proposal.Id, worker.MeasurementRequests[0].Payload.ProposalId);
        Assert.Equal(3, worker.MeasurementRequests[0].Payload.ResidueAddressesInOrder.Length);
        Assert.Equal(["A", "A", "A"], worker.MeasurementRequests[0].Payload.OutputChainIdsInOrder);
        Assert.Equal(["0:A:1::CA", "1:A:2::CA", "2:A:3::CA"],
            worker.MeasurementRequests[0].Payload.ExpectedResultAtomIdsInOrder);
        Assert.Equal(fixture.Proposal.OrientedProtein.CoordinateSha256,
            worker.MeasurementRequests[0].Payload.OrientedPdbSha256);
    }

    [Fact]
    public void Spanning_assembly_witness_must_cover_each_selected_chain_and_reject_a_reversed_chain()
    {
        using var fixture = new PlacementFixture();
        var selectedB = new ChainSelection("B", "copy-B");
        var intended = fixture.Protein.Intended with
        { Chains = fixture.Protein.Intended.Chains.Add(selectedB) };
        var bResidues = fixture.Residues.Select(address => address with
        { Chain = selectedB.SourceChain, CopyId = selectedB.CopyId }).ToImmutableArray();
        var bAtoms = bResidues.Select((address, index) => new AtomCorrespondence(
            index + 3, $"{index + 3}:B:{address.Residue}::CA", $"atom-B-{index + 1}",
            AtomOriginKind.Source, MoleculeRoleKind.Protein, AtomRoleKind.Backbone,
            "C", address, null)).ToImmutableArray();
        var protein = fixture.Protein with
        {
            Intended = intended,
            Molecule = fixture.Protein.Molecule with { AtomCount = 6 },
            Correspondence = fixture.Protein.Correspondence with
            { Atoms = fixture.Protein.Correspondence.Atoms.AddRange(bAtoms) }
        };
        var revision = fixture.Revision with { IntendedProtein = intended };
        var proposal = fixture.Proposal with
        { OrientedProtein = fixture.Proposal.OrientedProtein with { AtomCount = 6 } };
        var bMarkers = fixture.Witness.Residues.Select(marker => marker with
        { Residue = marker.Residue with { Chain = selectedB.SourceChain, CopyId = selectedB.CopyId } })
            .ToImmutableArray();
        var completeWitness = fixture.Witness with
        {
            ChainCopies = intended.Chains,
            Residues = fixture.Witness.Residues.AddRange(bMarkers)
        };
        static PlacementResidueObservation ObservedResidue(ResidueAddress address, int index) =>
            new(address, address.Chain, address.Residue.ToString(), "", "ALA", 1,
                index switch { 0 => 15, 1 => 0, _ => -15 },
                index switch { 0 => 15, 1 => 0, _ => -15 },
                index switch { 0 => 15, 1 => 0, _ => -15 },
                index == 1 ? 1 : 0, index == 1 ? 1 : 0,
                index == 0 ? 1 : 0, index == 2 ? 1 : 0);
        var residues = fixture.Residues.Select(ObservedResidue)
            .Concat(bResidues.Select(ObservedResidue)).ToImmutableArray();
        var report = new PlacementMeasurementReport(proposal.Id,
            ImmutableArray.Create(new MeasuredValue("atomsWithinCore", 2, "atoms", "oriented protein")),
            residues, ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<string>.Empty);
        var evidence = ImmutableArray.Create(new ScientificEvidence("witnessed", proposal.Id,
            "independent published topology", "Independently witnessed placement relationship",
            "Both chain copies are positioned", $"Selected membrane {fixture.Membrane.Id}",
            "Evidence applies to this exact assembly and bilayer.", EvidenceBearing.Supports));
        var owner = new PlacementOwner(new PlacementWorkerStub());
        AssessmentStanding Standing(PlacementStructuralWitness witness, PlacementMeasurementReport measured) =>
            owner.Assess(revision, protein, fixture.Membrane, proposal, fixture.Policy,
                measured, witness, evidence, ImmutableArray<ScientificFinding>.Empty).Standing;

        Assert.Equal(AssessmentStanding.Supported, Standing(completeWitness, report));
        var missingBMarkers = completeWitness with { Residues = fixture.Witness.Residues };
        Assert.Equal(AssessmentStanding.NotEstablished, Standing(missingBMarkers, report));
        var reversedB = report with { Residues = residues
            .SetItem(3, residues[3] with
                { MinZAngstrom = -15, MaxZAngstrom = -15, MeanZAngstrom = -15,
                  AtomsAboveCore = 0, AtomsBelowCore = 1 })
            .SetItem(5, residues[5] with
                { MinZAngstrom = 15, MaxZAngstrom = 15, MeanZAngstrom = 15,
                  AtomsAboveCore = 1, AtomsBelowCore = 0 }) };
        Assert.Equal(AssessmentStanding.NotEstablished, Standing(completeWitness, reversedB));

        var emptySelectionProtein = protein with
        { Intended = intended with { Chains = ImmutableArray<ChainSelection>.Empty } };
        var emptySelectionWitness = completeWitness with
        { ChainCopies = ImmutableArray<ChainSelection>.Empty };
        Assert.Equal(AssessmentStanding.NotEstablished,
            owner.Assess(revision, emptySelectionProtein, fixture.Membrane, proposal,
                fixture.Policy, report, emptySelectionWitness, evidence,
                ImmutableArray<ScientificFinding>.Empty).Standing);
    }

    [Fact]
    public async Task Pure_lipid_core_reframe_keeps_the_exact_PPM_proposal_and_distinct_boundary_provenance()
    {
        using var fixture = new PlacementFixture();
        var citation = "https://doi.org/10.1016/j.bbamem.2011.07.022";
        var review = "review:popc-reference-to-target";
        var frame = new PlacementPureLipidCoreFrame("popc-core-reference", "1", "POPC",
            28.8, 0.6, 303.15, "protein-free fully hydrated H2O/D2O vesicles; no stated NaCl treatment",
            citation, fixture.Membrane.Intended.Conditions, 1.5, review,
            "The 0.15 K and salt/protein-context differences are disclosed for policy review.");
        var policy = fixture.Policy with
        {
            EvidenceReferences = fixture.Policy.EvidenceReferences.Add(citation).Add(review),
            PureLipidCoreFrame = frame
        };
        var owner = new PlacementOwner(new PlacementWorkerStub());
        var proposed = owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal, policy);

        Assert.True(proposed.Established, proposed.Reason);
        Assert.NotEqual(fixture.Proposal.Id, proposed.Value!.Id);
        Assert.Equal(fixture.Proposal.OrientedProtein.CoordinateSha256,
            proposed.Value.OrientedProtein.CoordinateSha256);
        Assert.Equal(20, fixture.Proposal.ThicknessAngstrom);
        Assert.Equal(28.8, proposed.Value.ThicknessAngstrom);
        Assert.Equal(1.5, proposed.Value.MidplaneAngstrom);
        Assert.Contains(proposed.Value.Evidence, item =>
            item.Method == "PPM orientation candidate" && item.Bearing == EvidenceBearing.Context &&
            item.Observation.Contains(fixture.Proposal.Id, StringComparison.Ordinal));
        Assert.Contains(proposed.Value.Evidence, item =>
            item.Method == "pure-lipid hydrocarbon core reference" &&
            item.Source == citation && item.Bearing == EvidenceBearing.Context &&
            item.Uncertainty.Contains(review, StringComparison.Ordinal));
        Assert.All(proposed.Value.Evidence, item => Assert.Equal(proposed.Value.Id, item.SubjectId));
        Assert.DoesNotContain(fixture.Proposal.Evidence, item =>
            item.Method == "pure-lipid hydrocarbon core reference");
        var measuringWorker = new PlacementWorkerStub
        { MeasurementReply = request => Measured(request, fixture) };
        var measurement = await new PlacementOwner(measuringWorker).MeasureAgainstMembraneAsync(
            fixture.Revision, fixture.Protein, fixture.Membrane, proposed.Value,
            policy, fixture.Witness, "/tmp/placement-unit", TestContext.Current.CancellationToken);
        Assert.True(measurement.Established, measurement.Reason);
        Assert.Single(measuringWorker.MeasurementRequests);
        Assert.Equal(-12.9, measuringWorker.MeasurementRequests[0].Payload.CoreLowerZAngstrom, 6);
        Assert.Equal(15.9, measuringWorker.MeasurementRequests[0].Payload.CoreUpperZAngstrom, 6);

        PlacementSupportPolicy WithFrame(PlacementPureLipidCoreFrame changed) =>
            policy with { PureLipidCoreFrame = changed };
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            policy with { PureLipidCoreFrame = null }).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal, WithFrame(frame with { SpeciesId = "DOPC" })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal, WithFrame(frame with { ReferenceCitation = "" })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { SourceToTargetReviewRationale = "" })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { SourceToTargetReviewReference = "unidentified-review" })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { IntendedConditions = frame.IntendedConditions with
                { OptionalTemperatureKelvin = 310 } })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { HydrocarbonThicknessAngstrom = double.NaN })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { ReferenceTemperatureKelvin = double.NaN })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended, fixture.Proposal,
            WithFrame(frame with { MidplaneOffsetFromPpmAngstrom = double.PositiveInfinity })).Established);
        Assert.False(owner.ReframePpmProposalWithPolicy(fixture.Revision, fixture.Protein,
            fixture.Membrane.Intended,
            fixture.Proposal with { OrientedProtein = fixture.Proposal.OrientedProtein with
                { CoordinateSha256 = new string('a', 64) } }, policy).Established);
    }

    [Fact]
    public async Task Measurement_refuses_stale_or_swapped_prepared_atom_and_chain_correspondence()
    {
        using var fixture = new PlacementFixture();
        var changedChainWorker = new PlacementWorkerStub { MeasurementReply = request =>
        {
            var ordinary = Measured(request, fixture);
            return ordinary with { Observations = ordinary.Observations! with
            { Residues = ordinary.Observations.Residues.SetItem(0,
                ordinary.Observations.Residues[0] with { OutputChainId = "B" }) } };
        }};
        var wrongChain = await Measure(new PlacementOwner(changedChainWorker), fixture);
        Assert.False(wrongChain.Established);
        Assert.Single(changedChainWorker.MeasurementRequests);
        File.AppendAllText(fixture.Proposal.OrientedProtein.CoordinatePath, "changed");
        var staleWorker = new PlacementWorkerStub();
        var stale = await Measure(new PlacementOwner(staleWorker), fixture);
        Assert.False(stale.Established);
        Assert.Empty(staleWorker.MeasurementRequests);
    }

    [Fact]
    public async Task Independent_exact_witness_can_establish_explicit_membrane_transfer_when_implicit_ppm_transfer_is_disallowed()
    {
        using var fixture = new PlacementFixture();
        fixture.Policy = fixture.Policy with { AllowsTransferFromPpmDopc = false };
        var worker = new PlacementWorkerStub { MeasurementReply = request => Measured(request, fixture) };
        var owner = new PlacementOwner(worker);
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.Equal(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        Assert.Equal(AssessmentStanding.Supported, Assess(owner, fixture, measured.Value).Standing);
        Assert.Equal(AssessmentStanding.NotEstablished,
            owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Proposal,
                fixture.Policy, measured.Value, null, measured.Value.Evidence,
                ImmutableArray<ScientificFinding>.Empty).Standing);
    }

    [Fact]
    public async Task One_surface_witness_needs_the_selected_interface_and_physical_water_side()
    {
        using var fixture = new PlacementFixture();
        fixture.Proposal = fixture.Proposal with
        {
            TopologyKind = ProteinTopologyKind.OneSurfaceAssociated,
            PhysicalSide = PlacementPhysicalSide.Upper
        };
        fixture.Policy = fixture.Policy with
        { CoveredTopologyKinds = ImmutableArray.Create(ProteinTopologyKind.OneSurfaceAssociated) };
        fixture.Witness = fixture.Witness with
        {
            TopologyKind = ProteinTopologyKind.OneSurfaceAssociated,
            Residues = ImmutableArray.Create(
                new PlacementResidueWitness(fixture.Residues[0], PlacementWitnessRoleKind.Topology,
                    "upper-water", "independent published topology"),
                new PlacementResidueWitness(fixture.Residues[1], PlacementWitnessRoleKind.Contact,
                    "upper-interface", "independent published topology"),
                new PlacementResidueWitness(fixture.Residues[2], PlacementWitnessRoleKind.Sidedness,
                    "upper-water", "independent published topology"))
        };
        var worker = new PlacementWorkerStub { MeasurementReply = request =>
        {
            var result = Measured(request, fixture);
            var residues = result.Observations!.Residues;
            var onUpperSide = residues.SetItem(1, residues[1] with
                { MinZAngstrom = 10, MeanZAngstrom = 10, MaxZAngstrom = 10 })
                .SetItem(2, residues[2] with
                { MinZAngstrom = 15, MeanZAngstrom = 15, MaxZAngstrom = 15,
                  AtomsBelowCore = 0, AtomsAboveCore = 1 });
            return result with { Observations = result.Observations with
            { Residues = onUpperSide, AtomsAboveCore = 2, AtomsBelowCore = 0,
              ProteinZMinAngstrom = 10, ProteinZMaxAngstrom = 15 } };
        }};
        var owner = new PlacementOwner(worker);
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.Equal(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        Assert.Equal(AssessmentStanding.Supported, Assess(owner, fixture, measured.Value).Standing);
        fixture.Proposal = fixture.Proposal with { PhysicalSide = PlacementPhysicalSide.Lower };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
    }

    [Fact]
    public async Task Predicted_placement_needs_local_confidence_and_both_directional_pae_regions()
    {
        using var fixture = new PlacementFixture();
        var source = fixture.Protein.Intended.Source;
        var asset = new PredictionEvidenceAsset("AF-exact", "1", 1, 3,
            "https://prediction.example/model", source.Sha256, "https://prediction.example/pae",
            "/identified/pae.json", new string('f', 64), PaeAcquisitionStanding.Available, null);
        fixture.Protein = fixture.Protein with
        {
            Intended = fixture.Protein.Intended with
            { Source = source with { Kind = SourceRouteKind.AlphaFold, Prediction = asset } },
            Prediction = new PredictionEvidenceObservations("AF-exact", source.Sha256,
                fixture.Residues.Select(address => new PredictedResidueConfidence(
                    address with { CopyId = "" }, 95, PredictionObservationStanding.Observed, null))
                    .ToImmutableArray(),
                PredictionObservationStanding.Observed, null, 3,
                "/identified/pae-map.json", new string('e', 64), ImmutableArray<string>.Empty)
        };
        fixture.Policy = fixture.Policy with { PredictionCriterion = new PlacementPredictionCriterion(
            80, 1, 5, 1, false, 20) };
        var goodDirection = new DirectionalPredictionSummary(3, 3, 1, 3, 2,
            ImmutableArray<PredictionPairObservation>.Empty);
        var worker = new PlacementWorkerStub
        {
            MeasurementReply = request => Measured(request, fixture),
            PredictionReply = request => new WorkerResult<PredictionRegionSummaryObservations>(
                request.RequestId, request.Payload.StudyRevisionId, null, null,
                WorkerResultStanding.Observed, ImmutableArray<WorkerArtifact>.Empty,
                new PredictionRegionSummaryObservations("AF-exact", PredictionObservationStanding.Observed,
                    null, goodDirection, goodDirection), new ProviderIdentity("controlled PAE", "1"), null, null)
        };
        var owner = new PlacementOwner(worker);
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.Equal(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        Assert.Equal(AssessmentStanding.Supported, Assess(owner, fixture, measured.Value).Standing);
        Assert.Single(worker.PredictionRequests);
        Assert.Equal("AF-exact", worker.PredictionRequests[0].Payload.RecordId);

        var highPae = measured.Value with { Prediction = measured.Value.Prediction! with
        { SecondAlignedOnFirst = goodDirection with { MaximumAngstrom = 12 } } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, highPae).Standing);
        var highFirstDirection = measured.Value with { Prediction = measured.Value.Prediction! with
        { FirstAlignedOnSecond = goodDirection with { MaximumAngstrom = 12 } } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, highFirstDirection).Standing);
        var incompletePairs = measured.Value with { Prediction = measured.Value.Prediction! with
        { FirstAlignedOnSecond = goodDirection with { ValidPairCount = 2 } } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, incompletePairs).Standing);
        var completeLocal = fixture.Protein.Prediction!.LocalConfidence;
        fixture.Protein = fixture.Protein with { Prediction = fixture.Protein.Prediction with
        { LocalConfidence = completeLocal.SetItem(0, completeLocal[0] with { PLddt = 60 }) } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
        fixture.Protein = fixture.Protein with { Prediction = fixture.Protein.Prediction! with
        { LocalConfidence = fixture.Protein.Prediction.LocalConfidence.RemoveAt(0) } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
        fixture.Protein = fixture.Protein with { Prediction = fixture.Protein.Prediction! with
        { LocalConfidence = completeLocal } };
        var wrongRecord = measured.Value with { Prediction = measured.Value.Prediction! with { RecordId = "AF-other" } };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, wrongRecord).Standing);
        fixture.Policy = fixture.Policy with { PredictionCriterion = fixture.Policy.PredictionCriterion! with
        { IndependentWitnessCanResolvePredictionLimitations = true } };
        fixture.Witness = fixture.Witness with { PredictionResolutions = ImmutableArray.Create(
            new PlacementPredictionResolution(fixture.Residues,
                true, false, "independent published topology")) };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, highPae).Standing);
        fixture.Protein = fixture.Protein with { Prediction = fixture.Protein.Prediction! with
        { LocalConfidence = completeLocal.SetItem(0, completeLocal[0] with { PLddt = 60 }) } };
        fixture.Witness = fixture.Witness with { PredictionResolutions = ImmutableArray.Create(
            new PlacementPredictionResolution(fixture.Residues,
                false, true, "independent published topology")) };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
        fixture.Protein = fixture.Protein with { Prediction = fixture.Protein.Prediction! with
        { LocalConfidence = completeLocal } };
        fixture.Witness = fixture.Witness with { PredictionResolutions = ImmutableArray.Create(
            new PlacementPredictionResolution(fixture.Residues,
                true, true, "independent published topology")) };
        Assert.Equal(AssessmentStanding.Supported, Assess(owner, fixture, highPae).Standing);
        fixture.Witness = fixture.Witness with { PredictionResolutions = ImmutableArray.Create(
            new PlacementPredictionResolution(fixture.Residues.RemoveAt(0),
                true, true, "independent published topology")) };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, highPae).Standing);
    }

    [Fact]
    public async Task Duplicate_leaflet_witness_cannot_match_the_distinct_selected_mixture()
    {
        using var fixture = new PlacementFixture(mixed: true);
        var duplicate = new LeafletComposition(LeafletSide.Upper,
            ImmutableArray.Create(new LipidFraction("POPC", .5), new LipidFraction("POPC", .5)));
        fixture.Witness = fixture.Witness with { Upper = duplicate };
        var worker = new PlacementWorkerStub { MeasurementReply = request => Measured(request, fixture) };
        var owner = new PlacementOwner(worker);
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.NotEqual(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
    }

    [Fact]
    public async Task Insufficient_or_contradictory_evidence_has_distinct_standing()
    {
        using var fixture = new PlacementFixture();
        var worker = new PlacementWorkerStub { MeasurementReply = request => Measured(request, fixture) };
        var owner = new PlacementOwner(worker);
        var measured = (await Measure(owner, fixture)).Value!;
        Assert.Equal(AssessmentStanding.NotEstablished,
            owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Proposal,
                null, measured, fixture.Witness, measured.Evidence,
                ImmutableArray<ScientificFinding>.Empty).Standing);
        Assert.Equal(AssessmentStanding.NotEstablished,
            Assess(owner, fixture, measured with { ProposalId = "other" }).Standing);
        Assert.Equal(AssessmentStanding.NotEstablished,
            Assess(owner, fixture, measured, witness: fixture.Witness with
            { PreparedBondGraphSha256 = "wrong" }).Standing);
        var contradiction = new ScientificEvidence("contradiction", fixture.Proposal.Id, "independent source",
            "opposed topology", "Observed contact contradicts the proposed physical side.",
            fixture.Membrane.Intended.Id, "Exact observed region.", EvidenceBearing.Contradicts);
        Assert.Equal(AssessmentStanding.Unsupported,
            Assess(owner, fixture, measured, additional: ImmutableArray.Create(contradiction)).Standing);
        var refusals = new (PlacementSupportPolicy Policy, PlacementStructuralWitness Witness)[]
        {
            (fixture.Policy with { GeometryCriteria = ImmutableArray<PlacementMeasurementCriterion>.Empty }, fixture.Witness),
            (fixture.Policy with { GeometryCriteria = ImmutableArray.Create(
                new PlacementMeasurementCriterion("atomsWithinCore", "atoms", 2, null)) }, fixture.Witness),
            (fixture.Policy with { CoveredSpeciesIds = ImmutableArray.Create("DOPC") }, fixture.Witness),
            (fixture.Policy, fixture.Witness with { Upper = fixture.Witness.Upper with
                { Fractions = ImmutableArray.Create(new LipidFraction("DOPC", 1)) } }),
            (fixture.Policy, fixture.Witness with { ChainCopies = ImmutableArray.Create(new ChainSelection("B", "other")) }),
            (fixture.Policy, fixture.Witness with { RetainedPartnerSourceIds = ImmutableArray.Create("unselected") }),
            (fixture.Policy, fixture.Witness with { Residues = fixture.Witness.Residues.SetItem(2,
                fixture.Witness.Residues[2] with { Residue = new ResidueAddress(1, "A", 99, "", "copy-A") }) })
        };
        foreach (var (policy, witness) in refusals)
            Assert.Equal(AssessmentStanding.NotEstablished,
                owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Proposal,
                    policy, measured, witness, measured.Evidence,
                    ImmutableArray<ScientificFinding>.Empty).Standing);
    }

    [Fact]
    public async Task Retained_partner_cannot_gain_support_without_an_exact_partner_to_region_mapping()
    {
        using var fixture = new PlacementFixture();
        fixture.Protein = fixture.Protein with { Intended = fixture.Protein.Intended with
        { Partners = ImmutableArray.Create(new PartnerSelection("partner-X", true, "Retain cofactor.")) } };
        fixture.Witness = fixture.Witness with
        { RetainedPartnerSourceIds = ImmutableArray.Create("partner-X") };
        var owner = new PlacementOwner(new PlacementWorkerStub
        { MeasurementReply = request => Measured(request, fixture) });
        var measured = await Measure(owner, fixture);
        Assert.True(measured.Established, measured.Reason);
        Assert.NotEqual(EvidenceBearing.Supports, measured.Value!.Evidence.Single().Bearing);
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, measured.Value).Standing);
    }

    [Fact]
    public async Task Later_material_challenge_withholds_current_support_and_disqualification_is_unsupported()
    {
        using var fixture = new PlacementFixture();
        var owner = new PlacementOwner(new PlacementWorkerStub
        { MeasurementReply = request => Measured(request, fixture) });
        var measured = (await Measure(owner, fixture)).Value!;
        var prior = Assess(owner, fixture, measured);
        Assert.Equal(AssessmentStanding.Supported, prior.Standing);
        var finding = new ScientificFinding("later", fixture.Proposal.Id, measured.Evidence[0].Id,
            "New spatial finding", "Reconsider the placement premise", FindingDisposition.Challenges,
            true, DateTimeOffset.UtcNow);
        var challenged = owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane,
            fixture.Proposal, fixture.Policy, measured, fixture.Witness,
            measured.Evidence, ImmutableArray.Create(finding));
        Assert.Equal(AssessmentStanding.NotEstablished, challenged.Standing);
        Assert.Equal(AssessmentStanding.Supported, prior.Standing);
        var disqualified = owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane,
            fixture.Proposal, fixture.Policy, measured, fixture.Witness,
            measured.Evidence, ImmutableArray.Create(finding with
            { Disposition = FindingDisposition.Disqualifies }));
        Assert.Equal(AssessmentStanding.Unsupported, disqualified.Standing);
    }

    [Fact]
    public async Task Exact_witness_with_observed_opposite_sides_is_unsupported_but_missing_mapping_is_not_established()
    {
        using var fixture = new PlacementFixture();
        var oppositeWorker = new PlacementWorkerStub { MeasurementReply = request =>
        {
            var ordinary = Measured(request, fixture);
            var residues = ordinary.Observations!.Residues;
            return ordinary with { Observations = ordinary.Observations with
            {
                Residues = residues.SetItem(0, residues[0] with
                    { MinZAngstrom = -15, MaxZAngstrom = -15, MeanZAngstrom = -15,
                      AtomsAboveCore = 0, AtomsBelowCore = 1 })
                    .SetItem(2, residues[2] with
                    { MinZAngstrom = 15, MaxZAngstrom = 15, MeanZAngstrom = 15,
                      AtomsAboveCore = 1, AtomsBelowCore = 0 })
            }};
        }};
        var owner = new PlacementOwner(oppositeWorker);
        var opposite = await Measure(owner, fixture);
        Assert.True(opposite.Established, opposite.Reason);
        Assert.Contains(opposite.Value!.Evidence, item => item.Bearing == EvidenceBearing.Contradicts);
        Assert.Equal(AssessmentStanding.Unsupported, Assess(owner, fixture, opposite.Value).Standing);

        var missing = opposite.Value with { Residues = opposite.Value.Residues.RemoveAt(0),
            Evidence = ImmutableArray<ScientificEvidence>.Empty };
        Assert.Equal(AssessmentStanding.NotEstablished, Assess(owner, fixture, missing).Standing);
        var wrongSourceWitness = fixture.Witness with { SourceCoordinateSha256 = new string('0', 64) };
        var wrongSource = await owner.MeasureAgainstMembraneAsync(fixture.Revision, fixture.Protein,
            fixture.Membrane, fixture.Proposal, fixture.Policy, wrongSourceWitness,
            "/tmp/placement-unit", TestContext.Current.CancellationToken);
        Assert.True(wrongSource.Established, wrongSource.Reason);
        Assert.DoesNotContain(wrongSource.Value!.Evidence, item => item.Bearing == EvidenceBearing.Contradicts);
    }

    [Fact]
    public async Task Corrected_position_is_a_new_proposal_and_requires_fresh_measurement()
    {
        using var artifact = new IdentifiedArtifact();
        using var fixture = new PlacementFixture();
        var worker = new PlacementWorkerStub { AdjustmentReply = request => new WorkerResult<PlacementAdjustmentObservations>(
            request.RequestId, request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Observed,
            ImmutableArray.Create(new WorkerArtifact("adjustedPdb", artifact.Path, artifact.Sha256)),
            new PlacementAdjustmentObservations(3, 3, 1, 2, 0, 5, ImmutableArray<string>.Empty),
            new ProviderIdentity("controlled rigid transformation", "1"), null, null) };
        var owner = new PlacementOwner(worker);
        var revised = await owner.ReviseProposalAsync(fixture.Revision, fixture.Proposal,
            1, 2, 0, 5, "Move to inspect the upper loop contact.", "/tmp/placement-unit",
            TestContext.Current.CancellationToken);
        Assert.True(revised.Established, revised.Reason);
        Assert.NotEqual(fixture.Proposal.Id, revised.Value!.Id);
        Assert.Equal(fixture.Proposal.OrientedProtein.AtomCount, revised.Value.OrientedProtein.AtomCount);
        Assert.Null(revised.Value.TiltDegrees);
        Assert.Equal(1, revised.Value.MidplaneAngstrom);
        Assert.Equal(AssessmentStanding.NotEstablished,
            owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane,
                revised.Value, fixture.Policy, null, fixture.Witness,
                ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty).Standing);
        var bad = await owner.ReviseProposalAsync(fixture.Revision, fixture.Proposal,
            double.NaN, 0, 0, 0, "", "/tmp/placement-unit", TestContext.Current.CancellationToken);
        Assert.False(bad.Established);
        Assert.Single(worker.AdjustmentRequests);
    }

    private static Task<BoundaryOutcome<PlacementProposal>> Propose(PlacementOwner owner,
        PlacementFixture fixture) => owner.ProposeWithPpmAsync(fixture.Revision, fixture.Protein,
        fixture.Membrane.Intended, ProteinTopologyKind.MembraneSpanning, PlacementPhysicalSide.Both,
        "extracellular upper", PpmNterminalSide.Out, "/identified/immers", "2.0",
        new string('a', 64), "/tmp/placement-unit", TestContext.Current.CancellationToken,
        "/identified/res.lib", new string('f', 64));

    private static Task<BoundaryOutcome<PlacementMeasurementReport>> Measure(PlacementOwner owner,
        PlacementFixture fixture) => owner.MeasureAgainstMembraneAsync(fixture.Revision, fixture.Protein,
        fixture.Membrane, fixture.Proposal, fixture.Policy, fixture.Witness,
        "/tmp/placement-unit", TestContext.Current.CancellationToken);

    private static AssessedProteinMembranePlacement Assess(PlacementOwner owner, PlacementFixture fixture,
        PlacementMeasurementReport measurement, PlacementSupportPolicy? policy = null,
        PlacementStructuralWitness? witness = null, ImmutableArray<ScientificEvidence> additional = default) =>
        owner.Assess(fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Proposal,
            policy ?? fixture.Policy, measurement, witness ?? fixture.Witness,
            additional.IsDefault ? measurement.Evidence : additional,
            ImmutableArray<ScientificFinding>.Empty);

    private static WorkerResult<PlacementObservations> Observed(
        ScientificWorkRequest<PlacementPayload> request, IdentifiedArtifact artifact) =>
        new(request.RequestId, request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Observed,
            ImmutableArray.Create(new WorkerArtifact("orientedPdb", artifact.Path, artifact.Sha256)),
            new PlacementObservations(0, 20, 8, 3, ImmutableArray<MeasuredValue>.Empty,
                ImmutableArray.Create("N:1", "O:1"), ImmutableArray<string>.Empty,
                "PPM 2.0 implicit symmetric DOPC membrane"),
            new ProviderIdentity("PPM", "2.0"), null, null);

    private static WorkerResult<PlacementMeasurementObservations> Measured(
        ScientificWorkRequest<PlacementMeasurementPayload> request, PlacementFixture fixture)
    {
        var residues = fixture.Residues.Select((address, index) => new PlacementResidueObservation(
            address, "A", address.Residue.ToString(), "", "ALA", 1,
            index switch { 0 => 15, 1 => 0, _ => -15 },
            index switch { 0 => 15, 1 => 0, _ => -15 },
            index switch { 0 => 15, 1 => 0, _ => -15 },
            index == 1 ? 1 : 0, index == 1 ? 1 : 0,
            index == 0 ? 1 : 0, index == 2 ? 1 : 0)).ToImmutableArray();
        return new WorkerResult<PlacementMeasurementObservations>(request.RequestId,
            request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Observed,
            ImmutableArray<WorkerArtifact>.Empty,
            new PlacementMeasurementObservations(3, residues, 1, 1, 1, -15, 15,
                ImmutableArray<string>.Empty), new ProviderIdentity("controlled measurement", "1"), null, null);
    }

    private sealed class PlacementFixture : IDisposable
    {
        private readonly IdentifiedArtifact _orientedArtifact = new();
        public ImmutableArray<ResidueAddress> Residues { get; } =
            ImmutableArray.Create(new ResidueAddress(1, "A", 1, "", "copy-A"),
                new ResidueAddress(1, "A", 2, "", "copy-A"),
                new ResidueAddress(1, "A", 3, "", "copy-A"));
        public StudyRevision Revision { get; }
        public AssessedPreparedProtein Protein { get; set; }
        public AssessedMembraneModel Membrane { get; }
        public PlacementProposal Proposal { get; set; }
        public PlacementSupportPolicy Policy { get; set; }
        public PlacementStructuralWitness Witness { get; set; }

        public PlacementFixture(bool mixed = false)
        {
            var source = new StructuralSource("6QWR", SourceRouteKind.Rcsb, "RCSB experimental",
                "/identified/6qwr.pdb", new string('b', 64), "6QWR", "model 1");
            var intended = new IntendedProteinModel("intended", source, 1, "assembly-1",
                ImmutableArray.Create(new ChainSelection("A", "copy-A")),
                ImmutableArray<PartnerSelection>.Empty, ImmutableArray<AlternateLocationChoice>.Empty);
            var atoms = Residues.Select((address, index) => new AtomCorrespondence(index,
                $"{index}:A:{address.Residue}::CA", "atom-" + (index + 1), AtomOriginKind.Source,
                MoleculeRoleKind.Protein, AtomRoleKind.Backbone, "C", address, null)).ToImmutableArray();
            var molecule = new MolecularArtifact("prepared", "/identified/prepared.pdb", new string('c', 64),
                "/identified/prepared-topology.json", null, null, 3, null, new string('d', 64));
            Protein = new AssessedPreparedProtein("protein", "revision", intended, molecule, "chemical-policy",
                ImmutableArray<ResidueVariantChoice>.Empty, ImmutableArray<PreparationChangeProposal>.Empty,
                new SourceToResultCorrespondence(source.Id, molecule.Id, atoms, true),
                ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty,
                ImmutableArray<string>.Empty);
            var upper = new LeafletComposition(LeafletSide.Upper, mixed
                ? ImmutableArray.Create(new LipidFraction("POPC", .5), new LipidFraction("DOPC", .5))
                : ImmutableArray.Create(new LipidFraction("POPC", 1)));
            var lower = new LeafletComposition(LeafletSide.Lower,
                ImmutableArray.Create(new LipidFraction("POPC", 1)));
            var membrane = new MembraneModel("membrane", upper, lower,
                FixedStudyConditions.Initial, "identified planar bilayer");
            Membrane = new AssessedMembraneModel("assessed-membrane", "revision", membrane,
                ImmutableArray<MolecularRepresentation>.Empty, ImmutableArray<ScientificEvidence>.Empty,
                ImmutableArray<string>.Empty);
            Revision = new StudyRevision("revision", 3, intended, membrane, null, FixedStudyConditions.Initial);
            Proposal = new PlacementProposal("proposal", Protein.Id, membrane.Id,
                ProteinTopologyKind.MembraneSpanning,
                new MolecularArtifact("oriented", _orientedArtifact.Path, _orientedArtifact.Sha256,
                    molecule.TopologyPath, null, null, 3, null, molecule.TopologySha256),
                0, 20, 8, PlacementPhysicalSide.Both, "extracellular upper",
                ImmutableArray<string>.Empty,
                ImmutableArray.Create(new ScientificEvidence("ppm", "proposal", "PPM 2.0",
                    "PPM orientation candidate", "midplane 0 Å", "symmetric DOPC",
                    "Implicit symmetric membrane only.", EvidenceBearing.Context)),
                ImmutableArray<string>.Empty);
            Policy = new PlacementSupportPolicy("placement-policy", "1",
                ImmutableArray.Create(ProteinTopologyKind.MembraneSpanning),
                mixed ? ImmutableArray.Create("POPC", "DOPC") : ImmutableArray.Create("POPC"),
                mixed, mixed, true,
                ImmutableArray.Create("Independently witnessed placement relationship"),
                ImmutableArray.Create(new PlacementMeasurementCriterion("atomsWithinCore", "atoms", 1, null)),
                3, null, ImmutableArray.Create("independent published topology"),
                ImmutableArray<string>.Empty);
            Witness = new PlacementStructuralWitness("witness", "1", source.Sha256, 1, "assembly-1",
                intended.Chains, ImmutableArray<string>.Empty, null, molecule.TopologySha256,
                upper, lower, membrane.Conditions, ProteinTopologyKind.MembraneSpanning,
                "extracellular upper", "independent published topology",
                ImmutableArray.Create("independent published topology"),
                ImmutableArray.Create(
                    new PlacementResidueWitness(Residues[0], PlacementWitnessRoleKind.Topology,
                        "upper-water", "independent published topology"),
                    new PlacementResidueWitness(Residues[2], PlacementWitnessRoleKind.Topology,
                        "lower-water", "independent published topology"),
                    new PlacementResidueWitness(Residues[1], PlacementWitnessRoleKind.Contact,
                        "core-contact", "independent published topology"),
                    new PlacementResidueWitness(Residues[0], PlacementWitnessRoleKind.Sidedness,
                        "upper-water", "independent published topology")),
                ImmutableArray<string>.Empty);
        }

        public void Dispose() => _orientedArtifact.Dispose();
    }

    private sealed class PlacementWorkerStub : IPlacementAssessmentWork
    {
        public Func<ScientificWorkRequest<PlacementPayload>, WorkerResult<PlacementObservations>>? PpmReply { get; init; }
        public Func<ScientificWorkRequest<PlacementAdjustmentPayload>, WorkerResult<PlacementAdjustmentObservations>>? AdjustmentReply { get; init; }
        public Func<ScientificWorkRequest<PlacementMeasurementPayload>, WorkerResult<PlacementMeasurementObservations>>? MeasurementReply { get; init; }
        public Func<ScientificWorkRequest<PredictionRegionSummaryPayload>, WorkerResult<PredictionRegionSummaryObservations>>? PredictionReply { get; init; }
        public List<ScientificWorkRequest<PlacementPayload>> PpmRequests { get; } = [];
        public List<ScientificWorkRequest<PlacementAdjustmentPayload>> AdjustmentRequests { get; } = [];
        public List<ScientificWorkRequest<PlacementMeasurementPayload>> MeasurementRequests { get; } = [];
        public List<ScientificWorkRequest<PredictionRegionSummaryPayload>> PredictionRequests { get; } = [];
        public Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
            ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken)
        {
            PpmRequests.Add(request);
            return Task.FromResult(PpmReply!(request));
        }
        public Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
            ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken)
        {
            AdjustmentRequests.Add(request);
            return Task.FromResult(AdjustmentReply!(request));
        }
        public Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
            ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken)
        {
            MeasurementRequests.Add(request);
            return Task.FromResult(MeasurementReply!(request));
        }
        public Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
            ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken)
        {
            PredictionRequests.Add(request);
            return Task.FromResult(PredictionReply!(request));
        }
    }

    private sealed class IdentifiedArtifact : IDisposable
    {
        public string Path { get; } = System.IO.Path.Combine(System.IO.Path.GetTempPath(),
            "placement-owner-" + Guid.NewGuid().ToString("N") + ".pdb");
        public string Sha256 { get; }
        public IdentifiedArtifact()
        {
            File.WriteAllText(Path, "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C\nEND\n");
            Sha256 = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(Path)));
        }
        public void Dispose() => File.Delete(Path);
    }
}
