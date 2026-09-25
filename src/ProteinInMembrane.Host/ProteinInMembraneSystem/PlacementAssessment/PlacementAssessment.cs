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
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            protein.Id, reference.SourceUrl, "OPM oriented-structure reference",
            $"PDB {reference.PdbAccession}; membrane context {reference.MembraneContext}; " +
            $"hydrophobic thickness/depth {reference.HydrophobicThicknessAngstrom?.ToString("G6") ?? "unavailable"} Å; " +
            $"tilt {reference.TiltDegrees?.ToString("G6") ?? "unavailable"}°",
            $"Selected protein {protein.Id}; chosen membrane {membrane.Id}; exact construct correspondence {corresponds}",
            "This reference is contextual; it neither positions the prepared construct nor proves support for the chosen explicit lipid mixture.",
            EvidenceBearing.Context));
        return BoundaryOutcome<OpmReferenceReview>.Success(new OpmReferenceReview(
            protein.Id, membrane.Id, corresponds, evidence, limitations.ToImmutable()));
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
        CancellationToken cancellationToken)
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
            !ppmExecutableSha256.All(Uri.IsHexDigit))
            return BoundaryOutcome<PlacementProposal>.Unavailable("A declared topology, physical side, and identified PPM installation are required.");

        var request = new ScientificWorkRequest<PlacementPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new PlacementPayload(revision.Id, protein.Id, protein.Molecule.CoordinatePath,
                protein.Molecule.CoordinateSha256, ppmExecutablePath, ppmVersion,
                ppmExecutableSha256, topologyKind,
                ppmNterminalSide));
        var result = await _worker.PlacePpmAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementProposal>.Unavailable(result.FailureMessage ?? "PPM orientation was not observed.");
        var observed = result.Observations;
        var oriented = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "orientedPdb");
        if (oriented is null || observed.AlignedSourceAtomCount != protein.Molecule.AtomCount ||
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
            string.IsNullOrWhiteSpace(rationale))
            return BoundaryOutcome<PlacementProposal>.Unavailable("A finite adjustment, its rationale, and the corresponding current membrane are required.");
        var request = new ScientificWorkRequest<PlacementAdjustmentPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new PlacementAdjustmentPayload(revision.Id, source.Id, source.OrientedProtein.CoordinatePath,
                depthShiftAngstrom, tiltAboutXDegrees, tiltAboutYDegrees,
                rotationAboutNormalDegrees, rationale));
        var result = await _worker.AdjustPlacementAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementProposal>.Unavailable(result.FailureMessage ?? "The adjusted orientation was not observed.");
        var observed = result.Observations;
        var adjusted = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "adjustedPdb");
        if (adjusted is null || observed.SourceAtomCount != source.OrientedProtein.AtomCount ||
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
            !double.IsFinite(thickness) || thickness <= 0)
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable("Exact corresponding placement and membrane geometry are required for measurement.");

        var sourceResidues = protein.Correspondence.Atoms
            .OrderBy(atom => atom.ResultAtomIndex)
            .Select(atom => atom.SourceResidue)
            .Where(address => address is not null)
            .Cast<ResidueAddress>()
            .Distinct()
            .ToImmutableArray();
        if (sourceResidues.IsDefaultOrEmpty)
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable("The oriented protein has no preserved source-residue identities for placement measurement.");

        var request = new ScientificWorkRequest<PlacementMeasurementPayload>(Guid.NewGuid().ToString("N"),
            workingDirectory, new PlacementMeasurementPayload(revision.Id, proposal.Id,
                proposal.OrientedProtein.CoordinatePath, midplane,
                midplane - thickness / 2.0, midplane + thickness / 2.0, sourceResidues));
        var result = await _worker.MeasurePlacementAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementMeasurementReport>.Unavailable(result.FailureMessage ?? "Placement geometry was not observed.");
        var observed = result.Observations;
        if (observed.AtomCount != protein.Molecule.AtomCount ||
            observed.AtomsWithinCore + observed.AtomsAboveCore + observed.AtomsBelowCore != observed.AtomCount ||
            observed.Residues.Sum(residue => residue.AtomCount) != observed.AtomCount ||
            observed.Residues.Any(residue => residue.AtomsWithinCore + residue.AtomsAboveCore +
                residue.AtomsBelowCore != residue.AtomCount) ||
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
            applicable ? "Independently witnessed placement relationship" : "Oriented-protein placement geometry",
            string.Join("; ", measurements.Select(value => $"{value.Name}={value.Value:G6} {value.Unit}")),
            $"Prepared protein {protein.Id}; membrane {membrane.Id}; proposal {proposal.Id}; policy {policy?.Id ?? "unavailable"}; witness {witness?.Id ?? "unavailable"} version {witness?.Version ?? "unavailable"}",
            applicable ? "Support is bounded to the policy's evidenced topology, membrane class and predicted-model uncertainty gate; explicit packing may still fail." :
                "Raw geometry, prediction uncertainty or unmet policy criteria do not establish support for the chosen membrane.",
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
        PlacementMeasurementReport report,
        string orientedProteinUrl)
    {
        if (revision.Membrane?.Id != proposal.MembraneModelId || report.ProposalId != proposal.Id ||
            string.IsNullOrWhiteSpace(orientedProteinUrl) || report.Evidence.IsDefaultOrEmpty ||
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
            ImmutableArray<string>.Empty, report.Evidence, ImmutableArray<ScientificFinding>.Empty,
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
        if (policy is null || witness is null || report.ProposalId != proposal.Id ||
            string.IsNullOrWhiteSpace(witness.Id) || string.IsNullOrWhiteSpace(witness.Version) ||
            string.IsNullOrWhiteSpace(witness.Source) || witness.EvidenceReferences.IsDefaultOrEmpty ||
            !witness.Limitations.IsDefaultOrEmpty || witness.Residues.IsDefaultOrEmpty ||
            !double.IsFinite(policy.InterfaceBandAngstrom) || policy.InterfaceBandAngstrom <= 0 ||
            witness.SourceCoordinateSha256 != protein.Intended.Source.Sha256 ||
            witness.SourceModelIndex != protein.Intended.ModelIndex ||
            witness.BiologicalAssemblyId != protein.Intended.BiologicalAssemblyId ||
            !witness.ChainCopies.ToHashSet().SetEquals(protein.Intended.Chains) ||
            !witness.RetainedPartnerSourceIds.ToHashSet(StringComparer.Ordinal).SetEquals(
                protein.Intended.Partners.Where(item => item.Retain).Select(item => item.SourceId)) ||
            (witness.PreparedCoordinateSha256 is not null &&
                witness.PreparedCoordinateSha256 != protein.Molecule.CoordinateSha256) ||
            (witness.PreparedBondGraphSha256 is not null &&
                witness.PreparedBondGraphSha256 != protein.Molecule.TopologySha256) ||
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
            !witness.Residues.Any(item => item.Role == PlacementWitnessRoleKind.Sidedness))
            return false;

        var lower = midplane - thickness / 2.0;
        var upper = midplane + thickness / 2.0;
        foreach (var marker in witness.Residues)
        {
            if (!witness.EvidenceReferences.Contains(marker.EvidenceReference)) return false;
            if (!protein.Correspondence.Atoms.Any(atom => atom.SourceResidue == marker.Residue)) return false;
            var observed = report.Residues.FirstOrDefault(item => item.Address == marker.Residue);
            if (observed is null || observed.AtomCount <= 0) return false;
            var match = marker.ExpectedRegion switch
            {
                "upper-water" => observed.MinZAngstrom > upper,
                "lower-water" => observed.MaxZAngstrom < lower,
                "core-contact" => observed.AtomsWithinCore > 0,
                "upper-interface" => Math.Abs(observed.MeanZAngstrom - upper) <= policy.InterfaceBandAngstrom,
                "lower-interface" => Math.Abs(observed.MeanZAngstrom - lower) <= policy.InterfaceBandAngstrom,
                _ => false
            };
            if (!match) return false;
        }

        var topology = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Topology).ToArray();
        var contact = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Contact).ToArray();
        var sidedness = witness.Residues.Where(item => item.Role == PlacementWitnessRoleKind.Sidedness).ToArray();
        if (proposal.TopologyKind == ProteinTopologyKind.MembraneSpanning)
            return proposal.PhysicalSide == PlacementPhysicalSide.Both &&
                topology.Any(item => item.ExpectedRegion == "upper-water") &&
                topology.Any(item => item.ExpectedRegion == "lower-water") &&
                contact.Any(item => item.ExpectedRegion == "core-contact") &&
                sidedness.Any(item => item.ExpectedRegion is "upper-water" or "lower-water");
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

    private static bool SameLeaflet(LeafletComposition witness, LeafletComposition selected) =>
        witness.PhysicalSide == selected.PhysicalSide && witness.Fractions.Length == selected.Fractions.Length &&
        witness.Fractions.All(item => selected.Fractions.Any(other => other.SpeciesId == item.SpeciesId &&
            Math.Abs(other.Fraction - item.Fraction) <= 1e-9));
}
