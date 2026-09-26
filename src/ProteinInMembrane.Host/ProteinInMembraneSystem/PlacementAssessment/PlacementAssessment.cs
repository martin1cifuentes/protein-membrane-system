using System.Collections.Immutable;
using System.Security.Cryptography;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.PlacementAssessment;

/// <summary>Interprets an orientation for one exact protein and membrane pair.</summary>
public sealed class PlacementAssessment
{
    private readonly IPlacementAssessmentWork _worker;

    public PlacementAssessment(IPlacementAssessmentWork worker) => _worker = worker;

    public BoundaryOutcome<OpmReferenceReview> ReviewOpmReference(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        MembraneModel membrane,
        OpmReferenceRecord reference)
    {
        if (revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Id != protein.StudyRevisionId || revision.Membrane?.Id != membrane.Id)
            return BoundaryOutcome<OpmReferenceReview>.Unavailable("The OPM reference must be reviewed against one exact current protein and membrane choice.");
        var limitations = ImmutableArray.CreateBuilder<string>();
        if (!string.Equals(protein.Intended.Source.Accession, reference.PdbAccession,
                StringComparison.OrdinalIgnoreCase))
            limitations.Add("The OPM PDB accession differs from the selected source.");
        if (protein.Intended.BiologicalAssemblyId != reference.BiologicalAssemblyId ||
            !protein.Intended.Chains.ToHashSet().SetEquals(reference.ChainCopies))
            limitations.Add("The OPM assembly or exact retained chain copies differ from the selected protein.");
        var retainedPartners = protein.Intended.Partners.Where(item => item.Retain)
            .Select(item => item.SourceId).ToHashSet(StringComparer.Ordinal);
        if (!retainedPartners.SetEquals(reference.RetainedPartnerSourceIds))
            limitations.Add("The OPM retained partners differ from the selected construct.");
        var sourceAtoms = protein.Correspondence.Atoms.Where(atom => atom.SourceAtomId is not null)
            .Select(atom => atom.SourceAtomId!).ToArray();
        if (sourceAtoms.Length == 0 || reference.SourceAtomIds.IsDefaultOrEmpty ||
            sourceAtoms.Length != reference.SourceAtomIds.Length ||
            sourceAtoms.Distinct(StringComparer.Ordinal).Count() != sourceAtoms.Length ||
            reference.SourceAtomIds.Distinct(StringComparer.Ordinal).Count() != reference.SourceAtomIds.Length ||
            !sourceAtoms.ToHashSet(StringComparer.Ordinal).SetEquals(reference.SourceAtomIds))
            limitations.Add("The OPM oriented coordinate atoms cannot be mapped exactly to the retained source atoms.");
        if (string.IsNullOrWhiteSpace(reference.SourceUrl) ||
            string.IsNullOrWhiteSpace(reference.OrientedCoordinateSha256) ||
            !File.Exists(reference.OrientedCoordinatePath))
            limitations.Add("The identified OPM oriented-coordinate artifact is unavailable.");
        else
        {
            try
            {
                using var stream = File.OpenRead(reference.OrientedCoordinatePath);
                var actualHash = Convert.ToHexString(SHA256.HashData(stream));
                if (!actualHash.Equals(reference.OrientedCoordinateSha256,
                        StringComparison.OrdinalIgnoreCase))
                    limitations.Add("The OPM coordinate bytes no longer match their identified source hash.");
            }
            catch (IOException)
            {
                limitations.Add("The identified OPM oriented-coordinate artifact could not be read.");
            }
        }
        var corresponds = limitations.Count == 0;
        var upper = membrane.Upper.Fractions.Where(item => item.Fraction > 0).ToArray();
        var lower = membrane.Lower.Fractions.Where(item => item.Fraction > 0).ToArray();
        var membraneContextApplicable = reference.ImplicitSymmetric &&
            !string.IsNullOrWhiteSpace(reference.MembraneContext) &&
            !string.IsNullOrWhiteSpace(reference.AssumedMembraneSpeciesId) &&
            upper.Length == 1 && lower.Length == 1 &&
            upper[0].SpeciesId == reference.AssumedMembraneSpeciesId &&
            lower[0].SpeciesId == reference.AssumedMembraneSpeciesId &&
            Math.Abs(upper[0].Fraction - 1) <= 1e-9 &&
            Math.Abs(lower[0].Fraction - 1) <= 1e-9;
        if (!membraneContextApplicable)
            limitations.Add("The OPM implicit membrane context is unidentified or differs from the selected explicit membrane.");
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            protein.Id, reference.SourceUrl, "OPM oriented-structure reference",
            $"PDB {reference.PdbAccession}; membrane context {reference.MembraneContext}; " +
            $"hydrophobic thickness/depth {reference.HydrophobicThicknessAngstrom?.ToString("G6") ?? "unavailable"} Å; " +
            $"tilt {reference.TiltDegrees?.ToString("G6") ?? "unavailable"}°",
            $"Selected protein {protein.Id}; chosen membrane {membrane.Id}; exact construct correspondence {corresponds}; " +
            $"implicit membrane assumption consistent with selected species {membraneContextApplicable}",
            "This reference is contextual; it neither positions the prepared construct nor proves support for the chosen explicit lipid mixture.",
            EvidenceBearing.Context));
        return BoundaryOutcome<OpmReferenceReview>.Success(new OpmReferenceReview(
            protein.Id, membrane.Id, corresponds, evidence, limitations.ToImmutable(), membraneContextApplicable));
    }

    public BoundaryOutcome<PlacementProposal> ProposeFromOpmReference(
        StudyRevision revision, AssessedPreparedProtein protein, MembraneModel membrane,
        OpmReferenceRecord reference, OpmReferenceReview review,
        ProteinTopologyKind topologyKind, PlacementPhysicalSide physicalSide,
        string? biologicalSidedness)
    {
        if (revision.Id != protein.StudyRevisionId || revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Id || review.PreparedProteinId != protein.Id ||
            review.MembraneModelId != membrane.Id || !review.CorrespondsToSelectedConstruct ||
            !review.MembraneContextApplicable ||
            ReviewOpmReference(revision, protein, membrane, reference).Value is not
                { CorrespondsToSelectedConstruct: true, MembraneContextApplicable: true } ||
            !((topologyKind == ProteinTopologyKind.MembraneSpanning && physicalSide == PlacementPhysicalSide.Both) ||
              (topologyKind == ProteinTopologyKind.OneSurfaceAssociated &&
                  physicalSide is PlacementPhysicalSide.Upper or PlacementPhysicalSide.Lower)) ||
            reference.MidplaneAngstrom is not double midplane || !double.IsFinite(midplane) ||
            reference.HydrophobicThicknessAngstrom is not double thickness ||
            !double.IsFinite(thickness) || thickness <= 0 ||
            reference.TiltDegrees is double tilt && !double.IsFinite(tilt) ||
            !reference.OrientedCoordinateSha256.Equals(protein.Molecule.CoordinateSha256,
                StringComparison.OrdinalIgnoreCase) ||
            !ArtifactMatches(new WorkerArtifact("orientedPdb", reference.OrientedCoordinatePath,
                reference.OrientedCoordinateSha256)) ||
            !ArtifactMatches(new WorkerArtifact("preparedPdb", protein.Molecule.CoordinatePath,
                protein.Molecule.CoordinateSha256)))
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The OPM reference does not identify these exact prepared coordinates in its observed membrane frame.");
        var proposalId = Guid.NewGuid().ToString("N");
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            proposalId, reference.SourceUrl, "OPM exact-coordinate orientation candidate",
            $"The prepared coordinate digest equals the independently identified OPM oriented digest; " +
            $"midplane {midplane:G6} Å; hydrophobic thickness {thickness:G6} Å; tilt {reference.TiltDegrees?.ToString("G6") ?? "unavailable"}°",
            $"Prepared protein {protein.Id}; selected membrane {membrane.Id}; " +
            $"implicit {reference.AssumedMembraneSpeciesId} reference context",
            "This is a contextual position of already identical prepared coordinates; it does not establish support for the selected explicit bilayer.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(proposalId,
            protein.Id, membrane.Id, topologyKind,
            protein.Molecule with { Id = proposalId }, midplane, thickness, reference.TiltDegrees,
            physicalSide, biologicalSidedness, ImmutableArray<string>.Empty, evidence,
            ImmutableArray.Create("OPM implicit membrane reference; exact selected-bilayer support requires independent assessment.")));
    }

    public async Task<BoundaryOutcome<PlacementProposal>> ProposeWithPpmAsync(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        MembraneModel membrane,
        ProteinTopologyKind topologyKind,
        PlacementPhysicalSide physicalSide,
        string? biologicalSidedness,
        PpmNterminalSide ppmNterminalSide,
        string ppmExecutablePath,
        string ppmVersion,
        string ppmExecutableSha256,
        string workingDirectory,
        CancellationToken cancellationToken,
        string ppmResidueLibraryPath = "",
        string ppmResidueLibrarySha256 = "")
    {
        if (protein.StudyRevisionId != revision.Id || revision.Membrane?.Id != membrane.Id ||
            revision.IntendedProtein?.Id != protein.Intended.Id)
            return BoundaryOutcome<PlacementProposal>.Unavailable("Protein and membrane must belong to the exact current study revision.");
        if (!((topologyKind == ProteinTopologyKind.MembraneSpanning && physicalSide == PlacementPhysicalSide.Both) ||
              (topologyKind == ProteinTopologyKind.OneSurfaceAssociated &&
                  physicalSide is PlacementPhysicalSide.Upper or PlacementPhysicalSide.Lower)) ||
            !Enum.IsDefined(ppmNterminalSide) ||
            string.IsNullOrWhiteSpace(ppmExecutablePath) || string.IsNullOrWhiteSpace(ppmVersion) ||
            string.IsNullOrWhiteSpace(ppmExecutableSha256) || ppmExecutableSha256.Length != 64 ||
            !ppmExecutableSha256.All(Uri.IsHexDigit) ||
            string.IsNullOrWhiteSpace(ppmResidueLibraryPath) ||
            ppmResidueLibrarySha256.Length != 64 || !ppmResidueLibrarySha256.All(Uri.IsHexDigit))
            return BoundaryOutcome<PlacementProposal>.Unavailable("A declared topology, physical side, and identified PPM installation are required.");

        var request = new ScientificWorkRequest<PlacementPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new PlacementPayload(revision.Id, protein.Id, protein.Molecule.CoordinatePath,
                protein.Molecule.CoordinateSha256, ppmExecutablePath, ppmVersion,
                ppmExecutableSha256, topologyKind,
                ppmNterminalSide, ppmResidueLibraryPath, ppmResidueLibrarySha256));
        var result = await _worker.PlacePpmAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementProposal>.Unavailable(result.FailureMessage ?? "PPM orientation was not observed.");
        var observed = result.Observations;
        var oriented = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "orientedPdb");
        if (oriented is null || !ArtifactMatches(oriented) ||
            observed.AlignedSourceAtomCount != protein.Molecule.AtomCount ||
            observed.AssumedMembrane != "PPM 2.0 implicit symmetric DOPC membrane" ||
            !PairedPlaneMarkers(observed.PlaneMarkerIds) ||
            observed.TiltDegrees is double tilt && !double.IsFinite(tilt) ||
            observed.MidplaneAngstrom is null || observed.ThicknessAngstrom is null ||
            !double.IsFinite(observed.MidplaneAngstrom.Value) ||
            !double.IsFinite(observed.ThicknessAngstrom.Value) || observed.ThicknessAngstrom <= 0)
            return BoundaryOutcome<PlacementProposal>.Unavailable("PPM did not supply a corresponding positioned molecular model and membrane-boundary account.");

        var proposalId = Guid.NewGuid().ToString("N");
        var molecule = new MolecularArtifact(
            proposalId, oriented.Path, oriented.Sha256, protein.Molecule.TopologyPath, null, null,
            protein.Molecule.AtomCount, null, protein.Molecule.TopologySha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), proposalId, "PPM " + ppmVersion,
            "PPM orientation candidate",
            $"Midplane {observed.MidplaneAngstrom.Value:G6} Å; thickness {observed.ThicknessAngstrom.Value:G6} Å",
            $"Prepared protein {protein.Id}; assumed membrane {observed.AssumedMembrane}",
            "Implicit symmetric-membrane orientation is not itself support for the chosen explicit bilayer or biological sidedness.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(
            proposalId, protein.Id, membrane.Id, topologyKind, molecule,
            observed.MidplaneAngstrom, observed.ThicknessAngstrom, observed.TiltDegrees,
            physicalSide, biologicalSidedness, ImmutableArray<string>.Empty,
            evidence, observed.InterpretationWarnings));
    }

    public BoundaryOutcome<PlacementProposal> ReframePpmProposalWithPolicy(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        MembraneModel membrane,
        PlacementProposal ppmProposal,
        PlacementSupportPolicy policy)
    {
        var frame = policy.PureLipidCoreFrame;
        var ppmEvidence = ppmProposal.Evidence.IsDefaultOrEmpty ? null : ppmProposal.Evidence.FirstOrDefault(item =>
            item.SubjectId == ppmProposal.Id && item.Method == "PPM orientation candidate" &&
            item.Bearing == EvidenceBearing.Context &&
            item.Applicability?.Contains("DOPC", StringComparison.Ordinal) == true);
        if (frame is null ||
            revision.Id != protein.StudyRevisionId || revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Id || revision.Conditions != membrane.Conditions ||
            ppmProposal.PreparedProteinId != protein.Id || ppmProposal.MembraneModelId != membrane.Id ||
            ppmProposal.MidplaneAngstrom is not double ppmMidplane || !double.IsFinite(ppmMidplane) ||
            ppmProposal.ThicknessAngstrom is not double ppmThickness || !double.IsFinite(ppmThickness) ||
            ppmThickness <= 0 || ppmProposal.Evidence.Length != 1 || ppmEvidence is null ||
            ppmProposal.OrientedProtein.AtomCount != protein.Molecule.AtomCount ||
            ppmProposal.OrientedProtein.TopologySha256 != protein.Molecule.TopologySha256 ||
            !ArtifactMatches(new WorkerArtifact("orientedPdb", ppmProposal.OrientedProtein.CoordinatePath,
                ppmProposal.OrientedProtein.CoordinateSha256)) ||
            string.IsNullOrWhiteSpace(policy.Id) || string.IsNullOrWhiteSpace(policy.Version) ||
            policy.EvidenceReferences.IsDefaultOrEmpty ||
            !policy.CoveredTopologyKinds.Contains(ppmProposal.TopologyKind) ||
            !policy.CoveredSpeciesIds.Contains(frame.SpeciesId) ||
            string.IsNullOrWhiteSpace(frame.Id) || string.IsNullOrWhiteSpace(frame.Version) ||
            string.IsNullOrWhiteSpace(frame.SpeciesId) ||
            !double.IsFinite(frame.HydrocarbonThicknessAngstrom) ||
            frame.HydrocarbonThicknessAngstrom <= 0 ||
            !double.IsFinite(frame.ThicknessUncertaintyAngstrom) ||
            frame.ThicknessUncertaintyAngstrom <= 0 ||
            !double.IsFinite(frame.ReferenceTemperatureKelvin) || frame.ReferenceTemperatureKelvin <= 0 ||
            !double.IsFinite(frame.MidplaneOffsetFromPpmAngstrom) ||
            !double.IsFinite(ppmMidplane + frame.MidplaneOffsetFromPpmAngstrom) ||
            string.IsNullOrWhiteSpace(frame.ReferenceCondition) ||
            string.IsNullOrWhiteSpace(frame.ReferenceCitation) ||
            string.IsNullOrWhiteSpace(frame.SourceToTargetReviewReference) ||
            string.IsNullOrWhiteSpace(frame.SourceToTargetReviewRationale) ||
            frame.SourceToTargetReviewReference == frame.ReferenceCitation ||
            !Uri.TryCreate(frame.ReferenceCitation, UriKind.Absolute, out var citationUri) ||
            citationUri.Scheme != Uri.UriSchemeHttps ||
            !policy.EvidenceReferences.Contains(frame.ReferenceCitation) ||
            !policy.EvidenceReferences.Contains(frame.SourceToTargetReviewReference) ||
            frame.IntendedConditions != membrane.Conditions ||
            !PureSymmetricSpecies(membrane, frame.SpeciesId))
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "An exact PPM orientation and a separately identified, reviewed pure-lipid hydrocarbon frame for this chosen membrane are required.");

        var proposalId = Guid.NewGuid().ToString("N");
        var evidence = ImmutableArray.Create(ppmEvidence with
            {
                Id = Guid.NewGuid().ToString("N"), SubjectId = proposalId,
                Observation = $"From PPM proposal {ppmProposal.Id}; {ppmEvidence.Observation}"
            }).Add(new ScientificEvidence(Guid.NewGuid().ToString("N"), proposalId,
            frame.ReferenceCitation, "pure-lipid hydrocarbon core reference",
            $"2D_C {frame.HydrocarbonThicknessAngstrom:G6} ± " +
                $"{frame.ThicknessUncertaintyAngstrom:G6} Å; proposed midplane " +
                $"{ppmMidplane + frame.MidplaneOffsetFromPpmAngstrom:G6} Å, offset " +
                $"{frame.MidplaneOffsetFromPpmAngstrom:G6} Å from PPM midplane {ppmMidplane:G6} Å",
            $"Policy {policy.Id} version {policy.Version}; frame {frame.Id} version {frame.Version}; " +
                $"pure {frame.SpeciesId}; chosen membrane {membrane.Id}; intended conditions " +
                $"pH {membrane.Conditions.NominalPh:G6}, NaCl {membrane.Conditions.TargetNaClMolar:G6} M, " +
                $"temperature {membrane.Conditions.OptionalTemperatureKelvin:G6} K",
            $"Reference {frame.ReferenceTemperatureKelvin:G6} K, {frame.ReferenceCondition}. " +
                $"Context review {frame.SourceToTargetReviewReference}: {frame.SourceToTargetReviewRationale}. " +
                "This is an intended pure-lipid core frame, not a measured local protein-containing bilayer or placement support.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(ppmProposal with
        {
            Id = proposalId,
            MidplaneAngstrom = ppmMidplane + frame.MidplaneOffsetFromPpmAngstrom,
            ThicknessAngstrom = frame.HydrocarbonThicknessAngstrom,
            OrientedProtein = ppmProposal.OrientedProtein with { Id = proposalId },
            Evidence = evidence,
            Limitations = (ppmProposal.Limitations.IsDefault
                ? ImmutableArray<string>.Empty : ppmProposal.Limitations).Add(
                "The intended pure-lipid hydrocarbon frame is a sourced proposal assumption; " +
                "PPM's implicit DOPC boundaries and any protein-local POPC boundary remain distinct.")
        });
    }

    private static bool PureSymmetricSpecies(MembraneModel membrane, string speciesId) =>
        membrane.Upper.PhysicalSide == LeafletSide.Upper &&
        membrane.Lower.PhysicalSide == LeafletSide.Lower &&
        membrane.Upper.Fractions.Length == 1 && membrane.Lower.Fractions.Length == 1 &&
        membrane.Upper.Fractions[0].SpeciesId == speciesId &&
        membrane.Lower.Fractions[0].SpeciesId == speciesId &&
        double.IsFinite(membrane.Upper.Fractions[0].Fraction) &&
        double.IsFinite(membrane.Lower.Fractions[0].Fraction) &&
        Math.Abs(membrane.Upper.Fractions[0].Fraction - 1) <= 1e-9 &&
        Math.Abs(membrane.Lower.Fractions[0].Fraction - 1) <= 1e-9;

    public async Task<BoundaryOutcome<PlacementProposal>> ReviseProposalAsync(
        StudyRevision revision,
        PlacementProposal source,
        double depthShiftAngstrom,
        double tiltAboutXDegrees,
        double tiltAboutYDegrees,
        double rotationAboutNormalDegrees,
        string rationale,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        if (revision.Membrane?.Id != source.MembraneModelId ||
            !new[] { depthShiftAngstrom, tiltAboutXDegrees, tiltAboutYDegrees, rotationAboutNormalDegrees }.All(double.IsFinite) ||
            string.IsNullOrWhiteSpace(rationale) ||
            !ArtifactMatches(new WorkerArtifact("orientedPdb", source.OrientedProtein.CoordinatePath,
                source.OrientedProtein.CoordinateSha256)))
            return BoundaryOutcome<PlacementProposal>.Unavailable("A finite adjustment, its rationale, and the corresponding current membrane are required.");
        var request = new ScientificWorkRequest<PlacementAdjustmentPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new PlacementAdjustmentPayload(revision.Id, source.Id, source.OrientedProtein.CoordinatePath,
                depthShiftAngstrom, tiltAboutXDegrees, tiltAboutYDegrees,
                rotationAboutNormalDegrees, rationale, source.OrientedProtein.CoordinateSha256));
        var result = await _worker.AdjustPlacementAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementProposal>.Unavailable(result.FailureMessage ?? "The adjusted orientation was not observed.");
        var observed = result.Observations;
        var adjusted = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "adjustedPdb");
        if (adjusted is null || !ArtifactMatches(adjusted) ||
            observed.SourceAtomCount != source.OrientedProtein.AtomCount ||
            observed.AdjustedAtomCount != observed.SourceAtomCount ||
            !observed.GeometryWarnings.IsDefaultOrEmpty)
            return BoundaryOutcome<PlacementProposal>.Unavailable("The proposed adjustment did not preserve a corresponding molecular structure.");

        var proposalId = Guid.NewGuid().ToString("N");
        var artifact = new MolecularArtifact(proposalId, adjusted.Path, adjusted.Sha256,
            source.OrientedProtein.TopologyPath, null, null, observed.AdjustedAtomCount,
            null, source.OrientedProtein.TopologySha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), proposalId, "researcher adjustment",
            "Observed rigid-body placement adjustment", rationale,
            $"Derived from proposal {source.Id}",
            "The adjustment requires fresh assessment; earlier orientation energies and support do not transfer.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(
            proposalId, source.PreparedProteinId, source.MembraneModelId,
            source.TopologyKind, artifact,
            source.MidplaneAngstrom is null ? null : source.MidplaneAngstrom + observed.AppliedDepthShiftAngstrom,
            source.ThicknessAngstrom,
            // A scalar tilt cannot be updated by adding one Euler rotation. The
            // adjusted construct requires a fresh axis-based tilt observation.
            null,
            source.PhysicalSide, source.BiologicalSidedness, source.ContactingRegions,
            evidence, ImmutableArray.Create("Manual placement adjustment; new support assessment required.")));
    }

    public async Task<BoundaryOutcome<PlacementMeasurementReport>> MeasureAgainstMembraneAsync(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        AssessedMembraneModel membrane,
        PlacementProposal proposal,
        PlacementSupportPolicy? policy,
        PlacementStructuralWitness? witness,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        if (revision.Id != protein.StudyRevisionId || revision.Id != membrane.StudyRevisionId ||
            revision.Membrane?.Id != membrane.Intended.Id || proposal.PreparedProteinId != protein.Id ||
            proposal.MembraneModelId != membrane.Intended.Id || proposal.MidplaneAngstrom is not double midplane ||
            proposal.ThicknessAngstrom is not double thickness || !double.IsFinite(midplane) ||
            !double.IsFinite(thickness) || thickness <= 0 ||
            !ArtifactMatches(new WorkerArtifact("orientedPdb", proposal.OrientedProtein.CoordinatePath,
                proposal.OrientedProtein.CoordinateSha256)))
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable("Exact corresponding placement and membrane geometry are required for measurement.");

        var orderedAtoms = protein.Correspondence.Atoms.OrderBy(atom => atom.ResultAtomIndex).ToArray();
        if (!protein.Correspondence.Complete || orderedAtoms.Length != protein.Molecule.AtomCount ||
            orderedAtoms.Where((atom, index) => atom.ResultAtomIndex != index ||
                !ExactResultAtomId(atom.ResultAtomId, index)).Any())
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable(
                "The positioned artifact lacks exact ordered prepared-atom identities for measurement.");
        var sourceResidues = orderedAtoms
            .OrderBy(atom => atom.ResultAtomIndex)
            .Select(atom => atom.SourceResidue)
            .Where(address => address is not null)
            .Cast<ResidueAddress>()
            .Distinct()
            .ToImmutableArray();
        if (sourceResidues.IsDefaultOrEmpty)
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable("The oriented protein has no preserved source-residue identities for placement measurement.");
        var outputChainIds = ImmutableArray.CreateBuilder<string>();
        var outputResidueIds = ImmutableArray.CreateBuilder<string>();
        var outputInsertionCodes = ImmutableArray.CreateBuilder<string>();
        foreach (var residue in sourceResidues)
        {
            var outputResidues = orderedAtoms.Where(atom => atom.SourceResidue == residue)
                .Select(atom => atom.ResultAtomId.Split(':'))
                .Select(parts => (Chain: parts[1], Residue: parts[2], Insertion: parts[3]))
                .Distinct().ToArray();
            if (outputResidues.Length != 1)
                return BoundaryOutcome<PlacementMeasurementReport>.Unavailable(
                    "The prepared source residue has no unique output-chain and residue identity for placement measurement.");
            outputChainIds.Add(outputResidues[0].Chain);
            outputResidueIds.Add(outputResidues[0].Residue);
            outputInsertionCodes.Add(outputResidues[0].Insertion);
        }

        var request = new ScientificWorkRequest<PlacementMeasurementPayload>(Guid.NewGuid().ToString("N"),
            workingDirectory, new PlacementMeasurementPayload(revision.Id, proposal.Id,
                proposal.OrientedProtein.CoordinatePath, midplane,
                midplane - thickness / 2.0, midplane + thickness / 2.0, sourceResidues,
                outputChainIds.ToImmutable(), orderedAtoms.Select(atom => atom.ResultAtomId).ToImmutableArray(),
                proposal.OrientedProtein.CoordinateSha256));
        var result = await _worker.MeasurePlacementAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable(result.FailureMessage ?? "Placement geometry was not observed.");
        var observed = result.Observations;
        if (observed.AtomCount != protein.Molecule.AtomCount ||
            observed.AtomsWithinCore + observed.AtomsAboveCore + observed.AtomsBelowCore != observed.AtomCount ||
            observed.Residues.Length != sourceResidues.Length ||
            !observed.Residues.Select(item => item.Address).ToHashSet().SetEquals(sourceResidues) ||
            sourceResidues.Where((address, index) => !observed.Residues.Any(item =>
                item.Address == address && item.OutputChainId == outputChainIds[index] &&
                item.OutputResidueId == outputResidueIds[index] &&
                item.OutputInsertionCode == outputInsertionCodes[index])).Any() ||
            observed.Residues.Sum(residue => residue.AtomCount) != observed.AtomCount ||
            observed.Residues.Any(residue => residue.AtomsWithinCore + residue.AtomsAboveCore +
                residue.AtomsBelowCore != residue.AtomCount || residue.AtomCount <= 0 ||
                !double.IsFinite(residue.MinZAngstrom) || !double.IsFinite(residue.MeanZAngstrom) ||
                !double.IsFinite(residue.MaxZAngstrom) ||
                residue.MinZAngstrom > residue.MeanZAngstrom ||
                residue.MeanZAngstrom > residue.MaxZAngstrom) ||
            !double.IsFinite(observed.ProteinZMinAngstrom) || !double.IsFinite(observed.ProteinZMaxAngstrom))
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable("The oriented protein's atom and residue geometry is inconsistent.");

        var chargedCore = observed.Residues.Count(residue => residue.AtomsWithinCore > 0 &&
            residue.Name is "ASP" or "GLU" or "LYS" or "ARG");
        var measurements = ImmutableArray.Create(
            new MeasuredValue("atomsWithinCore", observed.AtomsWithinCore, "atoms", "oriented protein"),
            new MeasuredValue("atomsAboveCore", observed.AtomsAboveCore, "atoms", "oriented protein"),
            new MeasuredValue("atomsBelowCore", observed.AtomsBelowCore, "atoms", "oriented protein"),
            new MeasuredValue("proteinSpan", observed.ProteinZMaxAngstrom - observed.ProteinZMinAngstrom, "Å", "oriented protein"),
            new MeasuredValue("coreOccupancyFraction", (double)observed.AtomsWithinCore / observed.AtomCount, "fraction", "oriented protein"),
            new MeasuredValue("chargedResiduesWithCoreAtoms", chargedCore, "residues", "oriented protein"));
        var reportForWitness = new PlacementMeasurementReport(proposal.Id, measurements,
            observed.Residues, ImmutableArray<ScientificEvidence>.Empty, observed.Limitations);
        var witnessed = WitnessMatches(protein, membrane, proposal, policy, witness, reportForWitness);
        var contradicted = !witnessed && observed.Limitations.IsDefaultOrEmpty &&
            HasAttributableWitnessContradiction(protein, membrane, proposal, policy, witness, reportForWitness);
        var prediction = await SummarizePredictionForPlacementAsync(revision, protein, proposal,
            observed.Residues, witness, policy, workingDirectory, cancellationToken);
        var predictionAdequate = PredictionAdequate(protein, observed.Residues, witness,
            policy?.PredictionCriterion, prediction, witnessed);
        var applicable = policy is not null && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.GeometryCriteria.IsDefaultOrEmpty && witnessed && predictionAdequate &&
            !policy.EvidenceReferences.IsDefaultOrEmpty &&
            policy.CoveredTopologyKinds.Contains(proposal.TopologyKind) &&
            membrane.Intended.Upper.Fractions.Concat(membrane.Intended.Lower.Fractions)
                .Where(fraction => fraction.Fraction > 0)
                .All(fraction => policy.CoveredSpeciesIds.Contains(fraction.SpeciesId)) &&
            observed.Limitations.IsDefaultOrEmpty &&
            policy.GeometryCriteria.All(criterion =>
                measurements.Any(measurement => measurement.Name == criterion.MeasurementName &&
                    measurement.Unit == criterion.Unit &&
                    (criterion.Minimum is null || measurement.Value >= criterion.Minimum) &&
                    (criterion.Maximum is null || measurement.Value <= criterion.Maximum)));
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            proposal.Id, result.Provider?.Name ?? "local structural measurement",
            applicable ? "Independently witnessed placement relationship" :
                contradicted ? "Contradictory measured placement relationship" : "Oriented-protein placement geometry",
            string.Join("; ", measurements.Select(value => $"{value.Name}={value.Value:G6} {value.Unit}")),
            $"Prepared protein {protein.Id}; membrane {membrane.Id}; proposal {proposal.Id}; policy {policy?.Id ?? "unavailable"}; witness {witness?.Id ?? "unavailable"} version {witness?.Version ?? "unavailable"}",
            contradicted ? "An exact independent witness and measured region identify a material opposite-side or absent-contact relationship." :
                applicable ? "Support is bounded to the policy's evidenced topology, membrane class and predicted-model uncertainty gate; explicit packing may still fail." :
                "Raw geometry, prediction uncertainty or unmet policy criteria do not establish support for the chosen membrane.",
            contradicted ? EvidenceBearing.Contradicts :
                applicable ? EvidenceBearing.Supports : EvidenceBearing.Context));
        return BoundaryOutcome<PlacementMeasurementReport>.Success(new PlacementMeasurementReport(
            proposal.Id, measurements, observed.Residues, evidence, observed.Limitations, prediction));
    }

    private async Task<PredictionRegionSummaryObservations?> SummarizePredictionForPlacementAsync(
        StudyRevision revision, AssessedPreparedProtein protein, PlacementProposal proposal,
        ImmutableArray<PlacementResidueObservation> residues, PlacementStructuralWitness? witness,
        PlacementSupportPolicy? policy, string workingDirectory, CancellationToken cancellationToken)
    {
        if (protein.Intended.Source.Kind != SourceRouteKind.AlphaFold || protein.Intended.Source.Prediction is not { } asset ||
            protein.Prediction is not { PaeStanding: PredictionObservationStanding.Observed } prediction ||
            string.IsNullOrWhiteSpace(asset.PaePath) || string.IsNullOrWhiteSpace(asset.PaeSha256) ||
            string.IsNullOrWhiteSpace(prediction.PaeMappingPath) ||
            string.IsNullOrWhiteSpace(prediction.PaeMappingSha256) ||
            policy?.PredictionCriterion is not { MaximumReportedPaePairs: > 0 } criterion || witness is null)
            return null;
        var contacting = ContactingPredictionResidues(residues, witness);
        var sidedness = witness.Residues.Where(item => item.Role is PlacementWitnessRoleKind.Sidedness or PlacementWitnessRoleKind.Topology)
            .Select(item => item.Residue with { CopyId = string.Empty }).Distinct().ToImmutableArray();
        if (contacting.IsDefaultOrEmpty || sidedness.IsDefaultOrEmpty)
            return null;
        var request = new ScientificWorkRequest<PredictionRegionSummaryPayload>(Guid.NewGuid().ToString("N"),
            workingDirectory, new PredictionRegionSummaryPayload(revision.Id, asset.RecordId,
                asset.CoordinateSha256, asset.PaePath, asset.PaeSha256,
                prediction.PaeMappingPath, prediction.PaeMappingSha256,
                contacting, sidedness, criterion.MaximumReportedPaePairs));
        var result = await _worker.SummarizePredictionEvidenceAsync(request, cancellationToken);
        return result.RequestId == request.RequestId && result.StudyRevisionId == revision.Id &&
            result.Standing == WorkerResultStanding.Observed && result.Observations is { } summary &&
            summary.RecordId == asset.RecordId ? summary : null;
    }

    private static bool PredictionAdequate(AssessedPreparedProtein protein,
        ImmutableArray<PlacementResidueObservation> residues, PlacementStructuralWitness? witness,
        PlacementPredictionCriterion? criterion, PredictionRegionSummaryObservations? summary,
        bool witnessed)
    {
        if (protein.Intended.Source.Kind == SourceRouteKind.Upload)
            // Actor-declared provenance is not a mapped prediction-evidence asset.
            // A predicted or unknown upload cannot silently take the experimental route.
            return protein.Intended.Source.UploadProvenance == UploadOriginKind.Experimental;
        if (protein.Intended.Source.Kind != SourceRouteKind.AlphaFold) return true;
        if (protein.Intended.Source.Prediction is not { } predictionAsset ||
            protein.Prediction is not { } observedPrediction ||
            observedPrediction.RecordId != predictionAsset.RecordId ||
            observedPrediction.CoordinateSha256 != predictionAsset.CoordinateSha256 ||
            summary is not null && summary.RecordId != predictionAsset.RecordId)
            return false;
        if (criterion is null || !double.IsFinite(criterion.MinimumLocalConfidence) ||
            criterion.MinimumLocalConfidence < 0 || criterion.MinimumLocalConfidence > 100 ||
            !double.IsFinite(criterion.MinimumLocalCoverageFraction) ||
            criterion.MinimumLocalCoverageFraction <= 0 || criterion.MinimumLocalCoverageFraction > 1 ||
            !double.IsFinite(criterion.MaximumDirectionalPaeAngstrom) ||
            criterion.MaximumDirectionalPaeAngstrom < 0 ||
            !double.IsFinite(criterion.MinimumPaePairCoverageFraction) ||
            criterion.MinimumPaePairCoverageFraction <= 0 || criterion.MinimumPaePairCoverageFraction > 1)
            return false;
        var affected = ContactingPredictionResidues(residues, witness)
            .Concat(witness?.Residues.Where(item => item.Role is PlacementWitnessRoleKind.Sidedness or PlacementWitnessRoleKind.Topology)
                .Select(item => item.Residue with { CopyId = string.Empty }) ?? Enumerable.Empty<ResidueAddress>())
            .Distinct().ToHashSet();
        if (affected.Count == 0) return false;
        var local = protein.Prediction?.LocalConfidence
            .Where(item => affected.Contains(item.Residue))
            .GroupBy(item => item.Residue).ToDictionary(group => group.Key, group => group.ToArray());
        var observedLocal = local is null ? Array.Empty<PredictedResidueConfidence>() :
            local.Values.Where(group => group.Length == 1 && group[0].Standing == PredictionObservationStanding.Observed &&
                group[0].PLddt is double value && double.IsFinite(value) && value is >= 0 and <= 100)
                .Select(group => group[0]).ToArray();
        var localAdequate = (double)observedLocal.Length / affected.Count >=
            criterion.MinimumLocalCoverageFraction &&
            observedLocal.All(item => item.PLddt >= criterion.MinimumLocalConfidence);
        static bool DirectionAdequate(DirectionalPredictionSummary? direction,
            PlacementPredictionCriterion rule) =>
            direction is { PossiblePairCount: > 0, ValidPairCount: > 0, MaximumAngstrom: double maximum } &&
            direction.ValidPairCount <= direction.PossiblePairCount &&
            (double)direction.ValidPairCount / direction.PossiblePairCount >=
                rule.MinimumPaePairCoverageFraction &&
            double.IsFinite(maximum) && maximum <= rule.MaximumDirectionalPaeAngstrom;
        var relativeAdequate = summary is { Standing: PredictionObservationStanding.Observed } &&
            DirectionAdequate(summary.FirstAlignedOnSecond, criterion) &&
            DirectionAdequate(summary.SecondAlignedOnFirst, criterion);
        if (localAdequate && relativeAdequate) return true;
        if (!criterion.IndependentWitnessCanResolvePredictionLimitations || !witnessed || witness is null ||
            witness.PredictionResolutions.IsDefaultOrEmpty) return false;
        var exactResolutions = witness.PredictionResolutions.Where(resolution =>
            witness.EvidenceReferences.Contains(resolution.EvidenceReference) &&
            affected.IsSubsetOf(resolution.ResolvedResidues
                .Select(item => item with { CopyId = string.Empty }).ToHashSet())).ToArray();
        return (localAdequate || exactResolutions.Any(item => item.LocalStructureResolved)) &&
            (relativeAdequate || exactResolutions.Any(item => item.RelativePositionResolved));
    }

    private static ImmutableArray<ResidueAddress> ContactingPredictionResidues(
        ImmutableArray<PlacementResidueObservation> residues, PlacementStructuralWitness? witness)
    {
        var core = residues.Where(item => item.AtomsWithinCore > 0)
            .Select(item => item.Address with { CopyId = string.Empty });
        var witnessedContacts = witness?.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Contact &&
                (item.ExpectedRegion is "core-contact" or "upper-interface" or "lower-interface") &&
                residues.Any(observed => observed.Address == item.Residue && observed.AtomCount > 0))
            .Select(item => item.Residue with { CopyId = string.Empty }) ??
            Enumerable.Empty<ResidueAddress>();
        return core.Concat(witnessedContacts).Distinct().ToImmutableArray();
    }

    public BoundaryOutcome<InspectionSubject> ProposalInspectionSubject(
        StudyRevision revision,
        PlacementProposal proposal,
        PlacementMeasurementReport? report,
        string orientedProteinUrl,
        string? measurementIssue = null)
    {
        if (revision.Membrane?.Id != proposal.MembraneModelId ||
            string.IsNullOrWhiteSpace(orientedProteinUrl) ||
            proposal.Evidence.IsDefaultOrEmpty ||
            proposal.Evidence.Any(item => item.SubjectId != proposal.Id))
            return BoundaryOutcome<InspectionSubject>.Unavailable("This placement has no corresponding observed spatial and numerical account.");
        if (report is null)
        {
            var evidence = proposal.Evidence;
            if (!string.IsNullOrWhiteSpace(measurementIssue))
                evidence = evidence.Add(new ScientificEvidence(Guid.NewGuid().ToString("N"), proposal.Id,
                    "Placement Assessment", "Unavailable placement measurement", measurementIssue,
                    $"Prepared protein {proposal.PreparedProteinId}; membrane {proposal.MembraneModelId}; proposal {proposal.Id}",
                    "The positioned proposal is inspectable, but its contact and side relationship was not observed.",
                    EvidenceBearing.Unknown));
            var basis = proposal.Evidence[0];
            var proposedAnnotations = ImmutableArray.Create(new InspectionAnnotation(Guid.NewGuid().ToString("N"),
                proposal.Id, "proposed membrane plane",
                "The oriented protein is shown against the proposal's intended membrane bounds; contact remains unmeasured.",
                basis.Id, null));
            var proposedMetrics = ImmutableArray.CreateBuilder<InspectionMetric>();
            if (proposal.MidplaneAngstrom is double midplane && double.IsFinite(midplane))
                proposedMetrics.Add(new InspectionMetric("proposed midplane", midplane.ToString("G17"), "Å", proposal.Id, basis.Id));
            if (proposal.ThicknessAngstrom is double thickness && double.IsFinite(thickness))
                proposedMetrics.Add(new InspectionMetric("proposed thickness", thickness.ToString("G17"), "Å", proposal.Id, basis.Id));
            if (!string.IsNullOrWhiteSpace(measurementIssue))
                proposedMetrics.Add(new InspectionMetric("measurement", measurementIssue, null, proposal.Id, evidence[^1].Id));
            return BoundaryOutcome<InspectionSubject>.Success(new InspectionSubject(proposal.Id,
                revision.Id, orientedProteinUrl, "oriented-protein-with-proposed-membrane-bounds",
                ImmutableArray<string>.Empty, evidence, ImmutableArray<ScientificFinding>.Empty,
                proposedAnnotations, proposedMetrics.ToImmutable(), null));
        }
        if (report.ProposalId != proposal.Id || report.Evidence.IsDefaultOrEmpty ||
            report.Evidence.Any(item => item.SubjectId != proposal.Id))
            return BoundaryOutcome<InspectionSubject>.Unavailable("This placement has no corresponding observed spatial and numerical account.");
        var evidenceId = report.Evidence[0].Id;
        var anchors = new[]
        {
            (Name: "upper side", Residue: report.Residues.FirstOrDefault(item => item.AtomsAboveCore > 0)),
            (Name: "membrane core", Residue: report.Residues.FirstOrDefault(item => item.AtomsWithinCore > 0)),
            (Name: "lower side", Residue: report.Residues.FirstOrDefault(item => item.AtomsBelowCore > 0))
        };
        var annotations = ImmutableArray.CreateBuilder<InspectionAnnotation>();
        foreach (var anchor in anchors)
        {
            if (anchor.Residue is null ||
                !int.TryParse(anchor.Residue.OutputResidueId,
                    System.Globalization.NumberStyles.Integer,
                    System.Globalization.CultureInfo.InvariantCulture, out var sequenceId))
                continue;
            var residue = anchor.Residue;
            annotations.Add(new InspectionAnnotation(Guid.NewGuid().ToString("N"),
                $"{residue.Address.Chain}:{residue.Address.Residue}{residue.Address.InsertionCode}",
                anchor.Name, $"Observed {anchor.Name} relative to the proposal's membrane core bounds.",
                evidenceId, new StructureFocus(residue.OutputChainId, sequenceId,
                    string.IsNullOrWhiteSpace(residue.OutputInsertionCode) ? null : residue.OutputInsertionCode, null)));
        }
        if (annotations.Count == 0)
            return BoundaryOutcome<InspectionSubject>.Unavailable("No observed residue can spatially anchor this placement account.");
        var metrics = report.Measurements.Select(measurement => new InspectionMetric(measurement.Name,
            measurement.Value.ToString("G17", System.Globalization.CultureInfo.InvariantCulture),
            measurement.Unit, annotations[0].SubjectPartId, evidenceId)).ToImmutableArray();
        return BoundaryOutcome<InspectionSubject>.Success(new InspectionSubject(
            proposal.Id, revision.Id, orientedProteinUrl, "oriented-protein-with-proposed-membrane-bounds",
            ImmutableArray<string>.Empty, proposal.Evidence.AddRange(report.Evidence), ImmutableArray<ScientificFinding>.Empty,
            annotations.ToImmutable(), metrics, null));
    }

    public AssessedProteinMembranePlacement Assess(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        AssessedMembraneModel? membrane,
        PlacementProposal proposal,
        PlacementSupportPolicy? policy,
        PlacementMeasurementReport? measurement,
        PlacementStructuralWitness? witness,
        ImmutableArray<ScientificEvidence> additionalEvidence,
        ImmutableArray<ScientificFinding> findings)
    {
        var reviewedEvidence = additionalEvidence.IsDefault
            ? ImmutableArray<ScientificEvidence>.Empty : additionalEvidence;
        var allEvidence = proposal.Evidence.AddRange(reviewedEvidence);
        var reason = "Placement support is not established for these corresponding inputs.";
        var standing = AssessmentStanding.NotEstablished;

        if (proposal.PreparedProteinId != protein.Id || protein.StudyRevisionId != revision.Id ||
            revision.Membrane?.Id != proposal.MembraneModelId ||
            membrane?.StudyRevisionId != revision.Id || membrane.Intended.Id != proposal.MembraneModelId)
            reason = "Prepared protein, membrane, placement and study revision do not correspond.";
        else if (findings.Any(finding => finding.Material && finding.Disposition == FindingDisposition.Disqualifies &&
                     finding.SubjectId == proposal.Id) ||
                 allEvidence.Any(evidence => evidence.SubjectId == proposal.Id && evidence.Bearing == EvidenceBearing.Contradicts))
        {
            standing = AssessmentStanding.Unsupported;
            reason = "A material, attributable finding contradicts the proposed relationship.";
        }
        else if (policy is null || string.IsNullOrWhiteSpace(policy.Id) ||
                 string.IsNullOrWhiteSpace(policy.Version) || policy.EvidenceReferences.IsDefaultOrEmpty ||
                 policy.RequiredEvidenceMethods.IsDefaultOrEmpty)
            reason = "No versioned, evidence-backed placement support policy is available.";
        else if (protein.Intended.Partners.Any(item => item.Retain))
            reason = "Retained partners lack an exact source-partner-to-prepared-residue map and independently witnessed physical-side relationship.";
        else if (measurement is null || measurement.ProposalId != proposal.Id ||
                 !WitnessMatches(protein, membrane, proposal, policy, witness, measurement) ||
                 !PredictionAdequate(protein, measurement.Residues, witness,
                     policy.PredictionCriterion, measurement.Prediction, true) ||
                 !reviewedEvidence.Any(evidence => evidence.SubjectId == proposal.Id &&
                     evidence.Method == "Independently witnessed placement relationship" &&
                     evidence.Bearing == EvidenceBearing.Supports))
            reason = "Exact structural, topology, contact and sidedness witnesses were not observed and met for this placement.";
        else if (!policy.CoveredTopologyKinds.Contains(proposal.TopologyKind) ||
                 membrane.Intended.Upper.Fractions.Concat(membrane.Intended.Lower.Fractions)
                     .Where(fraction => fraction.Fraction > 0).Any(fraction => !policy.CoveredSpeciesIds.Contains(fraction.SpeciesId)) ||
                 !policy.AllowsMixtures && (membrane.Intended.Upper.Fractions.Count(f => f.Fraction > 0) > 1 ||
                                            membrane.Intended.Lower.Fractions.Count(f => f.Fraction > 0) > 1) ||
                 !policy.AllowsAsymmetry && !membrane.Intended.Upper.Fractions.SequenceEqual(membrane.Intended.Lower.Fractions))
            reason = "The selected topology or explicit membrane is outside the declared placement-policy scope.";
        else if (!policy.AllowsTransferFromPpmDopc && proposal.Evidence.Any(evidence =>
                     evidence.Method == "PPM orientation candidate") &&
                 !reviewedEvidence.Any(evidence => evidence.Bearing == EvidenceBearing.Supports &&
                     evidence.Applicability.Contains(membrane.Id, StringComparison.Ordinal)))
            reason = "The implicit-membrane orientation has not been shown applicable to this explicit bilayer.";
        else if (allEvidence.Any(evidence => evidence.SubjectId != proposal.Id) ||
                 policy.RequiredEvidenceMethods.Any(method => !allEvidence.Any(evidence =>
                     evidence.Method == method && evidence.Bearing == EvidenceBearing.Supports &&
                     evidence.Applicability.Contains(membrane.Id, StringComparison.Ordinal))))
            reason = "Required subject-specific placement, contact or sidedness evidence is absent or inapplicable.";
        else if (findings.Any(finding => finding.Material && finding.SubjectId == proposal.Id &&
                     finding.Disposition == FindingDisposition.Challenges))
            reason = "A later material finding requires reassessment of the placement premise.";
        else
        {
            standing = AssessmentStanding.Supported;
            reason = "Applicable evidence meets the identified placement policy without a material contradiction.";
        }

        return new AssessedProteinMembranePlacement(
            Guid.NewGuid().ToString("N"), revision.Id, proposal,
            standing, reason,
            findings.IsDefault ? ImmutableArray<ScientificFinding>.Empty : findings,
            DateTimeOffset.UtcNow);
    }

    private static bool WitnessMatches(
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        PlacementProposal proposal, PlacementSupportPolicy? policy,
        PlacementStructuralWitness? witness, PlacementMeasurementReport report)
    {
        if (!WitnessScopeMatches(protein, membrane, proposal, policy, witness, report))
            return false;
        var midplane = proposal.MidplaneAngstrom!.Value;
        var thickness = proposal.ThicknessAngstrom!.Value;
        var lower = midplane - thickness / 2.0;
        var upper = midplane + thickness / 2.0;
        foreach (var marker in witness!.Residues)
        {
            var observed = report.Residues.First(item => item.Address == marker.Residue);
            var match = marker.ExpectedRegion switch
            {
                "upper-water" => observed.MinZAngstrom > upper,
                "lower-water" => observed.MaxZAngstrom < lower,
                "core-contact" => observed.AtomsWithinCore > 0,
                "upper-interface" => Math.Abs(observed.MeanZAngstrom - upper) <= policy!.InterfaceBandAngstrom,
                "lower-interface" => Math.Abs(observed.MeanZAngstrom - lower) <= policy!.InterfaceBandAngstrom,
                _ => false
            };
            if (!match) return false;
        }

        var topology = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Topology).ToArray();
        var contact = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Contact).ToArray();
        var sidedness = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Sidedness).ToArray();
        if (proposal.TopologyKind == ProteinTopologyKind.MembraneSpanning)
            return proposal.PhysicalSide == PlacementPhysicalSide.Both &&
                witness.ChainCopies.All(chain =>
                    topology.Any(item => OnSelectedChain(item, chain) &&
                        item.ExpectedRegion == "upper-water") &&
                    topology.Any(item => OnSelectedChain(item, chain) &&
                        item.ExpectedRegion == "lower-water") &&
                    contact.Any(item => OnSelectedChain(item, chain) &&
                        item.ExpectedRegion == "core-contact") &&
                    sidedness.Any(item => OnSelectedChain(item, chain) &&
                        item.ExpectedRegion is "upper-water" or "lower-water"));
        if (proposal.TopologyKind == ProteinTopologyKind.OneSurfaceAssociated)
        {
            var side = proposal.PhysicalSide;
            var waterRegion = side switch
            {
                PlacementPhysicalSide.Upper => "upper-water",
                PlacementPhysicalSide.Lower => "lower-water",
                _ => null
            };
            var interfaceRegion = side == PlacementPhysicalSide.Upper ? "upper-interface" : "lower-interface";
            var oppositeWaterRegion = side == PlacementPhysicalSide.Upper ? "lower-water" : "upper-water";
            return waterRegion is not null &&
                topology.Any(item => item.ExpectedRegion == waterRegion) &&
                contact.Any(item => item.ExpectedRegion == interfaceRegion) &&
                sidedness.Any(item => item.ExpectedRegion == waterRegion) &&
                !witness.Residues.Any(item => item.ExpectedRegion == oppositeWaterRegion);
        }
        return false;
    }

    private static bool OnSelectedChain(PlacementResidueWitness marker, ChainSelection chain) =>
        marker.Residue.Chain == chain.SourceChain && marker.Residue.CopyId == chain.CopyId;

    private static bool WitnessScopeMatches(
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        PlacementProposal proposal, PlacementSupportPolicy? policy,
        PlacementStructuralWitness? witness, PlacementMeasurementReport report)
    {
        if (policy is null || witness is null || report.ProposalId != proposal.Id ||
            string.IsNullOrWhiteSpace(witness.Id) || string.IsNullOrWhiteSpace(witness.Version) ||
            string.IsNullOrWhiteSpace(witness.Source) || witness.EvidenceReferences.IsDefaultOrEmpty ||
            !witness.Limitations.IsDefaultOrEmpty || witness.Residues.IsDefaultOrEmpty ||
            policy.EvidenceReferences.IsDefaultOrEmpty || policy.GeometryCriteria.IsDefaultOrEmpty ||
            policy.GeometryCriteria.Any(criterion =>
                !report.Measurements.Any(measurement =>
                    measurement.Name == criterion.MeasurementName && measurement.Unit == criterion.Unit &&
                    double.IsFinite(measurement.Value) &&
                    (criterion.Minimum is null || measurement.Value >= criterion.Minimum) &&
                    (criterion.Maximum is null || measurement.Value <= criterion.Maximum))) ||
            !policy.CoveredTopologyKinds.Contains(proposal.TopologyKind) ||
            membrane.Intended.Upper.Fractions.Concat(membrane.Intended.Lower.Fractions)
                .Where(item => item.Fraction > 0)
                .Any(item => !policy.CoveredSpeciesIds.Contains(item.SpeciesId)) ||
            !policy.AllowsMixtures && (membrane.Intended.Upper.Fractions.Count(item => item.Fraction > 0) > 1 ||
                                       membrane.Intended.Lower.Fractions.Count(item => item.Fraction > 0) > 1) ||
            !policy.AllowsAsymmetry && !SameSymmetricLeaflets(membrane.Intended.Upper, membrane.Intended.Lower) ||
            !double.IsFinite(policy.InterfaceBandAngstrom) || policy.InterfaceBandAngstrom <= 0 ||
            witness.SourceCoordinateSha256 != protein.Intended.Source.Sha256 ||
            witness.SourceModelIndex != protein.Intended.ModelIndex ||
            witness.BiologicalAssemblyId != protein.Intended.BiologicalAssemblyId ||
            witness.ChainCopies.IsDefaultOrEmpty ||
            !witness.ChainCopies.ToHashSet().SetEquals(protein.Intended.Chains) ||
            !witness.RetainedPartnerSourceIds.ToHashSet(StringComparer.Ordinal).SetEquals(
                protein.Intended.Partners.Where(item => item.Retain).Select(item => item.SourceId)) ||
            protein.Intended.Partners.Any(item => item.Retain) ||
            (witness.PreparedCoordinateSha256 is not null &&
                witness.PreparedCoordinateSha256 != protein.Molecule.CoordinateSha256) ||
            string.IsNullOrWhiteSpace(witness.PreparedBondGraphSha256) ||
            witness.PreparedBondGraphSha256 != protein.Molecule.TopologySha256 ||
            !protein.Correspondence.Complete ||
            !SameLeaflet(witness.Upper, membrane.Intended.Upper) ||
            !SameLeaflet(witness.Lower, membrane.Intended.Lower) ||
            witness.Conditions != membrane.Intended.Conditions ||
            witness.TopologyKind != proposal.TopologyKind ||
            string.IsNullOrWhiteSpace(proposal.BiologicalSidedness) ||
            witness.BiologicalSidedness != proposal.BiologicalSidedness ||
            proposal.MidplaneAngstrom is not double midplane ||
            proposal.ThicknessAngstrom is not double thickness ||
            !double.IsFinite(midplane) || !double.IsFinite(thickness) || thickness <= 0 ||
            witness.Residues.Any(item => !Enum.IsDefined(item.Role)) ||
            !witness.Residues.Any(item => item.Role == PlacementWitnessRoleKind.Topology) ||
            !witness.Residues.Any(item => item.Role == PlacementWitnessRoleKind.Contact) ||
            !witness.Residues.Any(item => item.Role == PlacementWitnessRoleKind.Sidedness) ||
            report.Residues.IsDefaultOrEmpty ||
            report.Residues.Select(item => item.Address).Distinct().Count() != report.Residues.Length ||
            witness.Residues.Any(item => !witness.EvidenceReferences.Contains(item.EvidenceReference) ||
                !protein.Correspondence.Atoms.Any(atom => atom.SourceResidue == item.Residue) ||
                !report.Residues.Any(observed => observed.Address == item.Residue && observed.AtomCount > 0)))
            return false;
        return true;
    }

    private static bool HasAttributableWitnessContradiction(
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        PlacementProposal proposal, PlacementSupportPolicy? policy,
        PlacementStructuralWitness? witness, PlacementMeasurementReport report)
    {
        if (!WitnessScopeMatches(protein, membrane, proposal, policy, witness, report)) return false;
        var lower = proposal.MidplaneAngstrom!.Value - proposal.ThicknessAngstrom!.Value / 2;
        var upper = proposal.MidplaneAngstrom!.Value + proposal.ThicknessAngstrom!.Value / 2;
        return witness!.Residues.Any(marker =>
        {
            var observed = report.Residues.First(item => item.Address == marker.Residue);
            return marker.ExpectedRegion switch
            {
                "upper-water" => observed.MaxZAngstrom < lower,
                "lower-water" => observed.MinZAngstrom > upper,
                "core-contact" => observed.MaxZAngstrom < lower || observed.MinZAngstrom > upper,
                "upper-interface" => observed.MaxZAngstrom < lower,
                "lower-interface" => observed.MinZAngstrom > upper,
                _ => false
            };
        });
    }

    private static bool SameSymmetricLeaflets(LeafletComposition upper, LeafletComposition lower) =>
        SameDistinctFractions(upper.Fractions, lower.Fractions);

    private static bool SameLeaflet(LeafletComposition witness, LeafletComposition selected) =>
        witness.PhysicalSide == selected.PhysicalSide &&
        SameDistinctFractions(witness.Fractions, selected.Fractions);

    private static bool SameDistinctFractions(ImmutableArray<LipidFraction> first,
        ImmutableArray<LipidFraction> second)
    {
        if (first.IsDefaultOrEmpty || second.IsDefaultOrEmpty || first.Length != second.Length ||
            first.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() != first.Length ||
            second.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() != second.Length)
            return false;
        var orderedFirst = first.OrderBy(item => item.SpeciesId, StringComparer.Ordinal).ToArray();
        var orderedSecond = second.OrderBy(item => item.SpeciesId, StringComparer.Ordinal).ToArray();
        return orderedFirst.Zip(orderedSecond).All(pair =>
            pair.First.SpeciesId == pair.Second.SpeciesId &&
            double.IsFinite(pair.First.Fraction) && double.IsFinite(pair.Second.Fraction) &&
            Math.Abs(pair.First.Fraction - pair.Second.Fraction) <= 1e-9);
    }

    private static bool ArtifactMatches(WorkerArtifact artifact)
    {
        if (string.IsNullOrWhiteSpace(artifact.Path) || artifact.Sha256.Length != 64 ||
            !artifact.Sha256.All(Uri.IsHexDigit) || !File.Exists(artifact.Path)) return false;
        try
        {
            using var stream = File.OpenRead(artifact.Path);
            return Convert.ToHexString(SHA256.HashData(stream))
                .Equals(artifact.Sha256, StringComparison.OrdinalIgnoreCase);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException) { return false; }
    }

    private static bool PairedPlaneMarkers(ImmutableArray<string> markers)
    {
        if (markers.IsDefaultOrEmpty || markers.Length % 2 != 0) return false;
        for (var index = 0; index < markers.Length; index += 2)
        {
            if (!markers[index].StartsWith("N:", StringComparison.Ordinal) ||
                !markers[index + 1].StartsWith("O:", StringComparison.Ordinal) ||
                markers[index].Length <= 2 ||
                markers[index][2..] != markers[index + 1][2..]) return false;
        }
        return markers.Distinct(StringComparer.Ordinal).Count() == markers.Length;
    }

    private static bool ExactResultAtomId(string? resultAtomId, int index)
    {
        if (string.IsNullOrWhiteSpace(resultAtomId)) return false;
        var parts = resultAtomId.Split(':');
        return parts.Length == 5 &&
            int.TryParse(parts[0], out var actualIndex) && actualIndex == index &&
            !string.IsNullOrWhiteSpace(parts[1]) && int.TryParse(parts[2], out _) &&
            !string.IsNullOrWhiteSpace(parts[4]);
    }
}
