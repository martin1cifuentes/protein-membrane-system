using System.Collections.Immutable;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ProteinPreparationOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinPreparation.ProteinPreparation;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed partial class ProteinPreparationRouteTests
{
    [Fact]
    public async Task Prediction_confidence_requires_the_exact_archive_record_and_never_inherits_upload_provenance()
    {
        using var directory = new TemporaryDirectory();
        var sourcePath = System.IO.Path.Combine(directory.Path, "prediction.cif");
        File.WriteAllText(sourcePath, "source coordinates");
        var digest = Hash(sourcePath);
        var asset = new PredictionEvidenceAsset("AF-P12345-F1", "4", 1, 1,
            "https://alphafold.ebi.ac.uk/files/prediction.cif", digest,
            null, null, null, PaeAcquisitionStanding.Unavailable, "No PAE was supplied.");
        var predicted = new StructuralSource("alphafold:AF-P12345-F1", SourceRouteKind.AlphaFold,
            "identified AlphaFold DB prediction", sourcePath, digest, "P12345", "AF-P12345-F1", asset);
        var mismatched = new PredictionEvidenceObservations("AF-P12345-F2", digest,
            ImmutableArray.Create(new PredictedResidueConfidence(SourceResidue, 93,
                PredictionObservationStanding.Observed, null)),
            PredictionObservationStanding.Unavailable, "No PAE was supplied.", null, null, null,
            ImmutableArray<string>.Empty);
        var worker = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("mmcif", ImmutableArray.Create(Model()), mismatched))
        };
        var outcome = await new ProteinPreparationOwner(worker).InspectSourceAsync(predicted,
            directory.Path, 1000, TestContext.Current.CancellationToken);
        Assert.False(outcome.Established);
        Assert.Contains("same prediction", outcome.Reason);

        var upload = predicted with
        {
            Id = "upload:fixture", Kind = SourceRouteKind.Upload, Prediction = null,
            UploadProvenance = UploadOriginKind.Predicted
        };
        var uploadWithArchiveConfidence = await new ProteinPreparationOwner(worker).InspectSourceAsync(upload,
            directory.Path, 1000, TestContext.Current.CancellationToken);
        Assert.False(uploadWithArchiveConfidence.Established);
        Assert.Contains("cannot be attributed", uploadWithArchiveConfidence.Reason);

        var noArchiveConfidence = new PreparationWorkerStub
        {
            InspectSourceResponse = request => Observed(request.RequestId, null,
                new SourceInspectionObservations("mmcif", ImmutableArray.Create(Model()), null))
        };
        var declaredUpload = await new ProteinPreparationOwner(noArchiveConfidence).InspectSourceAsync(upload,
            directory.Path, 1000, TestContext.Current.CancellationToken);
        Assert.True(declaredUpload.Established);
        Assert.Null(declaredUpload.Value!.Prediction);
        Assert.Equal(UploadOriginKind.Predicted, declaredUpload.Value.Source.UploadProvenance);
    }

    [Fact]
    public async Task Candidate_with_incomplete_atom_correspondence_cannot_be_promoted()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "ALA", ["N", "CA", "C", "O"]);
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateResponse(request, directory.Path,
                context.Intended.Source.Sha256, ["N", "CA", "C", "O"], ObservedGeometry(),
                correspondenceComplete: false)
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(),
            context.Proposals, ImmutableArray<ResearcherDecision>.Empty, directory.Path,
            TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains("atom-for-atom", result.Reason);
        Assert.Null(result.Diagnostic);
    }

    [Fact]
    public async Task Measured_out_of_policy_candidate_geometry_remains_unqualified_with_a_material_finding()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "ALA", ["N", "CA", "C", "O"]);
        var measured = new ProteinGeometryObservations(ObservationStanding.Observed,
            ImmutableArray.Create(new ProteinGeometryKindObservation("covalentBond",
                GeometryKindStanding.Observed, 1, 1, 2.4, 2.4, null)),
            ImmutableArray<GeometryDistanceObservation>.Empty, ImmutableArray<string>.Empty);
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateResponse(request, directory.Path,
                context.Intended.Source.Sha256, ["N", "CA", "C", "O"], measured)
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(),
            context.Proposals, ImmutableArray<ResearcherDecision>.Empty, directory.Path,
            TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        var diagnostic = Assert.IsType<ProteinPreparationDiagnostic>(result.Diagnostic);
        Assert.Equal(ObservationStanding.Observed, diagnostic.Geometry.Standing);
        Assert.Contains(diagnostic.Findings, finding =>
            finding.SubjectId == diagnostic.Candidate.Id && finding.Material &&
            finding.Disposition == FindingDisposition.Disqualifies);
    }

    [Fact]
    public async Task Histidine_needs_one_current_reviewed_variant_with_a_recorded_rationale()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "HIS", ["N", "CA", "C", "O", "CB", "CG"]);
        var policy = context.Chemical with { PermittedVariants = ImmutableArray.Create("HID", "HIE", "HIP") };
        var hid = new PreparationChangeProposal("his-hid", context.Revision.Id, context.Intended.Id,
            SelectedResidue, PreparationChangeKind.ResidueState, "HID", "Review this site",
            ImmutableArray<string>.Empty, true);
        var hie = hid with { Id = "his-hie", ProposedChange = "HIE" };
        var proposals = context.Proposals with { Changes = ImmutableArray.Create(hid, hie) };
        var worker = new PreparationWorkerStub();
        var owner = new ProteinPreparationOwner(worker);

        var absent = await owner.PrepareAsync(context.Revision, context.Intended, context.Inspection,
            policy, StructuralPolicy(), proposals, ImmutableArray<ResearcherDecision>.Empty,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(absent.Established);
        Assert.Contains("Exactly one reviewed", absent.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        var noRationale = ImmutableArray.Create(new ResearcherDecision("choose-hid", context.Revision.Id,
            hid.Id, ResearcherDecisionKind.ApprovePreparationChange, ResearcherDecisionValue.Approved,
            DateTimeOffset.UtcNow, string.Empty));
        var unreasoned = await owner.PrepareAsync(context.Revision, context.Intended, context.Inspection,
            policy, StructuralPolicy(), proposals, noRationale,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(unreasoned.Established);
        Assert.Contains("recorded rationale", unreasoned.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        var conflicting = noRationale.SetItem(0, noRationale[0] with { Rationale = "At this site" })
            .Add(new ResearcherDecision("choose-hie", context.Revision.Id, hie.Id,
                ResearcherDecisionKind.ApprovePreparationChange, ResearcherDecisionValue.Approved,
                DateTimeOffset.UtcNow, "Competing assignment"));
        var conflict = await owner.PrepareAsync(context.Revision, context.Intended, context.Inspection,
            policy, StructuralPolicy(), proposals, conflicting,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(conflict.Established);
        Assert.Contains("Exactly one reviewed", conflict.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        ProteinPreparationPayload? submitted = null;
        var observedWorker = new PreparationWorkerStub
        {
            PrepareResponse = request =>
            {
                submitted = request.Payload;
                return ControlledFailure(request);
            }
        };
        var reasoned = await new ProteinPreparationOwner(observedWorker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, policy, StructuralPolicy(), proposals,
            ImmutableArray.Create(conflicting[0]), directory.Path, TestContext.Current.CancellationToken);
        Assert.False(reasoned.Established);
        Assert.Equal(1, observedWorker.PrepareProteinCalls);
        var variant = Assert.Single(submitted!.ResidueVariants);
        Assert.Equal("HID", variant.Variant);
        Assert.Equal("choose-hid", variant.DecisionId);
    }

    [Fact]
    public async Task Possible_disulfide_requires_an_explicit_reviewed_bond_decision()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "CYS", ["N", "CA", "C", "O", "CB", "SG"]);
        var second = new ResidueAddress(0, "A", 2, "", "A");
        var inspectedModel = context.Inspection.Models[0] with
        {
            Chains = ImmutableArray.Create(new SourceChainObservation("A", 2, 12)),
            Residues = context.Inspection.Models[0].Residues.Add(
                context.Inspection.Models[0].Residues[0] with
                { Address = second with { CopyId = string.Empty } }),
            AtomCount = 12
        };
        var inspection = context.Inspection with { Models = ImmutableArray.Create(inspectedModel) };
        var policy = context.Chemical with
        { PermittedVariants = context.Chemical.PermittedVariants.Add("CYX") };
        var bond = new PreparationChangeProposal("possible-cys-pair", context.Revision.Id,
            context.Intended.Id, SelectedResidue, PreparationChangeKind.Disulfide,
            "A:1–A:2", "SG proximity suggests a possible bond", ImmutableArray<string>.Empty,
            true, second);
        var proposals = context.Proposals with { Changes = ImmutableArray.Create(bond) };
        var worker = new PreparationWorkerStub();
        var owner = new ProteinPreparationOwner(worker);

        var absent = await owner.PrepareAsync(context.Revision, context.Intended, inspection,
            policy, StructuralPolicy(), proposals, ImmutableArray<ResearcherDecision>.Empty,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(absent.Established);
        Assert.Contains("explicitly approved or declined", absent.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        var withoutRationale = ImmutableArray.Create(new ResearcherDecision("approve-bond",
            context.Revision.Id, bond.Id, ResearcherDecisionKind.ApprovePreparationChange,
            ResearcherDecisionValue.Approved, DateTimeOffset.UtcNow, string.Empty));
        var unreasoned = await owner.PrepareAsync(context.Revision, context.Intended, inspection,
            policy, StructuralPolicy(), proposals, withoutRationale,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(unreasoned.Established);
        Assert.Contains("recorded rationale", unreasoned.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        var conflicting = withoutRationale.SetItem(0,
            withoutRationale[0] with { Rationale = "Site evidence supports a bond" })
            .Add(new ResearcherDecision("decline-bond", context.Revision.Id, bond.Id,
                ResearcherDecisionKind.ApprovePreparationChange, ResearcherDecisionValue.Declined,
                DateTimeOffset.UtcNow, "Keep these cysteines separate"));
        var conflict = await owner.PrepareAsync(context.Revision, context.Intended, inspection,
            policy, StructuralPolicy(), proposals, conflicting,
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(conflict.Established);
        Assert.Contains("conflicting", conflict.Reason);
        Assert.Equal(0, worker.PrepareProteinCalls);

        ProteinPreparationPayload? submitted = null;
        var observedWorker = new PreparationWorkerStub
        {
            PrepareResponse = request =>
            {
                submitted = request.Payload;
                return ControlledFailure(request);
            }
        };
        var approved = await new ProteinPreparationOwner(observedWorker).PrepareAsync(context.Revision,
            context.Intended, inspection, policy, StructuralPolicy(), proposals,
            ImmutableArray.Create(withoutRationale[0] with { Rationale = "Site evidence supports a bond" }),
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(approved.Established);
        Assert.Equal(1, observedWorker.PrepareProteinCalls);
        var approvedPair = Assert.Single(submitted!.ApprovedDisulfides);
        Assert.Equal(SelectedResidue, approvedPair.First);
        Assert.Equal(second, approvedPair.Second);
        Assert.Equal(["CYX", "CYX"], submitted.ResidueVariants.Select(item => item.Variant));

        var mismatchedWorker = new PreparationWorkerStub
        {
            PrepareResponse = request => Observed(request.RequestId, request.Payload.StudyRevisionId,
                new ProteinPreparationObservations(12, 12, 2,
                    ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
                    ImmutableArray<ResidueAddress>.Empty, ImmutableArray<ResidueAddress>.Empty,
                    ImmutableArray<AtomAddress>.Empty, ImmutableArray<AtomAddress>.Empty,
                    ImmutableArray<AtomAddress>.Empty,
                    ImmutableArray.Create(new ResidueVariantChoice(SelectedResidue, "CYX", "approve-bond"),
                        new ResidueVariantChoice(second, "CYX", "approve-bond")),
                    ImmutableArray<DisulfideBond>.Empty, ImmutableArray<string>.Empty,
                    12, ObservedGeometry()))
        };
        var bondMismatch = await new ProteinPreparationOwner(mismatchedWorker).PrepareAsync(
            context.Revision, context.Intended, inspection, policy, StructuralPolicy(), proposals,
            ImmutableArray.Create(withoutRationale[0] with { Rationale = "Site evidence supports a bond" }),
            directory.Path, TestContext.Current.CancellationToken);
        Assert.False(bondMismatch.Established);
        Assert.Contains("disulfide bonds do not match", bondMismatch.Reason);
    }

    private static WorkerResult<ProteinPreparationObservations> ControlledFailure(
        ScientificWorkRequest<ProteinPreparationPayload> request) =>
        new(request.RequestId, request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Failed,
            ImmutableArray<WorkerArtifact>.Empty, null, null, "controlled-provider-failure",
            "The controlled provider made no candidate.");

    [Fact]
    public async Task Unavailable_actual_candidate_geometry_keeps_an_inspectable_unqualified_candidate_and_reason()
    {
        using var directory = new TemporaryDirectory();
        var sourcePath = System.IO.Path.Combine(directory.Path, "source.pdb");
        File.WriteAllText(sourcePath, "source coordinates");
        var source = new StructuralSource("upload:fixture", SourceRouteKind.Upload, "researcher upload",
            sourcePath, Hash(sourcePath), null, "fixture", null, UploadOriginKind.Experimental);
        var intended = new IntendedProteinModel("intended", source, 0, null,
            ImmutableArray.Create(new ChainSelection("A", "A")),
            ImmutableArray<PartnerSelection>.Empty, ImmutableArray<AlternateLocationChoice>.Empty);
        var revision = new StudyRevision("revision", 1, intended, null, null, FixedStudyConditions.Initial);
        var inspection = new SourceInspectionReport(source, "pdb", ImmutableArray.Create(Model()),
            ImmutableArray<string>.Empty);
        var proposals = new PreparationProposalReport(revision.Id, intended.Id,
            ImmutableArray<PreparationChangeProposal>.Empty, null,
            ImmutableArray<PreviewChainCorrespondence>.Empty,
            ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<string>.Empty);
        var chemical = new ProteinChemicalStatePolicy("chemical", "1", "canonical-amino-acid-assembly", 2.5,
            ImmutableArray.Create("policy evidence"),
            ImmutableArray.Create(new ForceFieldAsset("ff", "1", "Amber19", "/unused/ff.xml", "abc")),
            ImmutableDictionary<string, string>.Empty, ImmutableArray.Create("HID"), ImmutableArray<string>.Empty);
        var structural = StructuralPolicy();
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateWithUnavailableGeometry(request, directory.Path, source.Sha256)
        };

        var outcome = await new ProteinPreparationOwner(worker).PrepareAsync(revision, intended,
            inspection, chemical, structural, proposals, ImmutableArray<ResearcherDecision>.Empty,
            directory.Path, TestContext.Current.CancellationToken);

        Assert.False(outcome.Established);
        var diagnostic = Assert.IsType<ProteinPreparationDiagnostic>(outcome.Diagnostic);
        Assert.Contains(GeometryReason, outcome.Reason);
        Assert.Equal(ObservationStanding.Unavailable, diagnostic.Geometry.Standing);
        Assert.Contains(diagnostic.Findings, finding => finding.SubjectId == diagnostic.Candidate.Id &&
            finding.Meaning.Contains(GeometryReason, StringComparison.Ordinal));
        Assert.Equal(1, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Approved_heavy_atom_missing_from_actual_candidate_is_not_promoted()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "ALA", ["N", "CA", "C", "O"]);
        var proposal = new PreparationChangeProposal("add-cb", context.Revision.Id, context.Intended.Id,
            SelectedResidue, PreparationChangeKind.HeavyAtom, "CB", "Observed missing CB",
            ImmutableArray<string>.Empty, true);
        var report = context.Proposals with { Changes = ImmutableArray.Create(proposal) };
        var decisions = ImmutableArray.Create(new ResearcherDecision("approve-cb", context.Revision.Id,
            proposal.Id, ResearcherDecisionKind.ApprovePreparationChange,
            ResearcherDecisionValue.Approved, DateTimeOffset.UtcNow, "Approved exact CB"));
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateResponse(request, directory.Path,
                context.Intended.Source.Sha256, ["N", "CA", "C", "O"], ObservedGeometry())
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(), report,
            decisions, directory.Path, TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains("heavy-atom", result.Reason);
        Assert.Equal(1, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Actual_heavy_atom_without_an_exact_approval_is_not_promoted()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "ALA", ["N", "CA", "C", "O"]);
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request =>
            {
                var candidate = CandidateResponse(request, directory.Path,
                    context.Intended.Source.Sha256, ["N", "CA", "C", "O", "CB"], ObservedGeometry());
                return candidate with
                {
                    Observations = candidate.Observations! with
                    {
                        AddedHeavyAtoms = ImmutableArray.Create(new AtomAddress(SelectedResidue, "CB"))
                    }
                };
            }
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(),
            context.Proposals, ImmutableArray<ResearcherDecision>.Empty, directory.Path,
            TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains("heavy-atom additions do not match", result.Reason);
        Assert.Equal(1, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Required_explicit_variant_missing_from_actual_candidate_is_not_promoted()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "CYS", ["N", "CA", "C", "O", "CB", "SG"]);
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateResponse(request, directory.Path,
                context.Intended.Source.Sha256, ["N", "CA", "C", "O", "CB", "SG"], ObservedGeometry())
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(),
            context.Proposals, ImmutableArray<ResearcherDecision>.Empty, directory.Path,
            TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains("residue state", result.Reason);
        Assert.Equal(1, worker.PrepareProteinCalls);
    }

    [Fact]
    public async Task Corresponding_explicit_variant_and_measured_candidate_can_establish_only_a_prepared_protein()
    {
        using var directory = new TemporaryDirectory();
        var context = OwnerContext(directory.Path, "CYS", ["N", "CA", "C", "O", "CB", "SG"]);
        var worker = new PreparationWorkerStub
        {
            PrepareResponse = request => CandidateResponse(request, directory.Path,
                context.Intended.Source.Sha256, ["N", "CA", "C", "O", "CB", "SG"], ObservedGeometry(),
                actualVariants: ImmutableArray.Create(new ResidueVariantChoice(SelectedResidue, "CYS", null)))
        };

        var result = await new ProteinPreparationOwner(worker).PrepareAsync(context.Revision,
            context.Intended, context.Inspection, context.Chemical, StructuralPolicy(),
            context.Proposals, ImmutableArray<ResearcherDecision>.Empty, directory.Path,
            TestContext.Current.CancellationToken);

        var protein = Assert.IsType<AssessedPreparedProtein>(result.Value);
        Assert.Equal(context.Revision.Id, protein.StudyRevisionId);
        Assert.Equal(context.Intended.Id, protein.Intended.Id);
        Assert.Contains(protein.Limitations, limitation => limitation.Contains("placement", StringComparison.Ordinal));
    }
}
