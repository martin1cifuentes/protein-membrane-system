using System.Collections.Immutable;
using System.Globalization;
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
        if (reference.BiologicalAssemblyId is not null &&
            protein.Intended.BiologicalAssemblyId != reference.BiologicalAssemblyId ||
            !reference.ChainCopies.IsDefaultOrEmpty &&
            !protein.Intended.Chains.ToHashSet().SetEquals(reference.ChainCopies))
            limitations.Add("The OPM assembly or exact retained chain copies differ from the selected protein.");
        var retainedPartners = protein.Intended.Partners.Where(item => item.Retain)
            .Select(item => item.SourceId).ToHashSet(StringComparer.Ordinal);
        if (!reference.RetainedPartnerSourceIds.IsDefaultOrEmpty &&
            !retainedPartners.SetEquals(reference.RetainedPartnerSourceIds))
            limitations.Add("The OPM retained partners differ from the selected construct.");
        var sourceAtoms = protein.Correspondence.Atoms.Where(atom => atom.SourceAtomId is not null)
            .Select(atom => atom.SourceAtomId!).ToArray();
        if (sourceAtoms.Length == 0 ||
            !reference.SourceAtomIds.IsDefaultOrEmpty &&
            (sourceAtoms.Length != reference.SourceAtomIds.Length ||
             sourceAtoms.Distinct(StringComparer.Ordinal).Count() != sourceAtoms.Length ||
             reference.SourceAtomIds.Distinct(StringComparer.Ordinal).Count() != reference.SourceAtomIds.Length ||
             !sourceAtoms.ToHashSet(StringComparer.Ordinal).SetEquals(reference.SourceAtomIds)))
            limitations.Add("The OPM oriented coordinate atoms cannot be mapped exactly to the retained source atoms.");
        if (!TryOpmAlignment(protein, reference, out _, out var mappingIssue))
            limitations.Add(mappingIssue);
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
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            protein.Id, reference.SourceUrl, "OPM oriented-structure reference",
            $"PDB {reference.PdbAccession}; membrane context {reference.MembraneContext}; " +
            $"hydrophobic thickness/depth {reference.HydrophobicThicknessAngstrom?.ToString("G6") ?? "unavailable"} Å; " +
            $"tilt {reference.TiltDegrees?.ToString("G6") ?? "unavailable"}°",
            $"Selected protein {protein.Id}; chosen membrane {membrane.Id}; exact construct correspondence {corresponds}; " +
            $"provider membrane assumption matches selected pure species {membraneContextApplicable}",
            "The provider membrane assumption is reference provenance; the selected construct needs an independently verified rigid mapping.",
            EvidenceBearing.Context));
        return BoundaryOutcome<OpmReferenceReview>.Success(new OpmReferenceReview(
            protein.Id, membrane.Id, corresponds, evidence, limitations.ToImmutable(), membraneContextApplicable));
    }

    public BoundaryOutcome<PlacementProposal> ProposeFromOpmReference(
        StudyRevision revision, AssessedPreparedProtein protein, MembraneModel membrane,
        PlacementSupportPolicy framePolicy,
        OpmReferenceRecord reference, OpmReferenceReview review,
        ProteinTopologyKind topologyKind, PlacementPhysicalSide physicalSide,
        string? biologicalSidedness)
    {
        if (revision.Id != protein.StudyRevisionId || revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Id || review.PreparedProteinId != protein.Id ||
            review.MembraneModelId != membrane.Id || !review.CorrespondsToSelectedConstruct ||
            ReviewOpmReference(revision, protein, membrane, reference).Value is not
                { CorrespondsToSelectedConstruct: true } ||
            !ValidFramePolicy(framePolicy, topologyKind, membrane) ||
            !Enum.IsDefined(topologyKind) || !Enum.IsDefined(physicalSide) ||
            reference.MidplaneAngstrom is not double midplane || !double.IsFinite(midplane) ||
            reference.HydrophobicThicknessAngstrom is not double thickness ||
            !double.IsFinite(thickness) || thickness <= 0 ||
            reference.TiltDegrees is double tilt && !double.IsFinite(tilt) ||
            !TryOpmAlignment(protein, reference, out var transform, out _))
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The OPM reference has no verified rigid mapping from this exact prepared construct into its membrane frame.");
        var proposalId = Guid.NewGuid().ToString("N");
        var output = Path.Combine(Path.GetDirectoryName(reference.OrientedCoordinatePath)!,
            $"opm-proposal-{proposalId}.pdb");
        try
        {
            var transformed = File.ReadAllLines(protein.Molecule.CoordinatePath);
            for (var i = 0; i < transformed.Length; i++)
            {
                var line = transformed[i];
                if (!IsPdbAtom(line)) continue;
                if (!TryPdbPoint(line, out var point))
                    return BoundaryOutcome<PlacementProposal>.Unavailable(
                        "The prepared PDB has an uninterpretable atom coordinate.");
                var oriented = transform!(point);
                var moved = oriented with { Z = oriented.Z - midplane };
                if (!double.IsFinite(moved.X) || !double.IsFinite(moved.Y) ||
                    !double.IsFinite(moved.Z) ||
                    new[] { moved.X, moved.Y, moved.Z }.Any(value => Math.Abs(value) >= 10000))
                    return BoundaryOutcome<PlacementProposal>.Unavailable(
                        "The mapped OPM orientation exceeds the PDB coordinate representation.");
                transformed[i] = line[..30] +
                    string.Format(CultureInfo.InvariantCulture, "{0,8:F3}{1,8:F3}{2,8:F3}",
                        moved.X, moved.Y, moved.Z) + line[54..];
            }
            File.WriteAllLines(output, transformed);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        {
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The mapped OPM proposal could not be written as a separate prepared-construct artifact.");
        }
        using var outputStream = File.OpenRead(output);
        var outputHash = Convert.ToHexString(SHA256.HashData(outputStream)).ToLowerInvariant();
        var molecule = new MolecularArtifact(proposalId, output, outputHash,
            protein.Molecule.TopologyPath, null, null, protein.Molecule.AtomCount,
            null, protein.Molecule.TopologySha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            proposalId, reference.SourceUrl, "OPM mapped rigid orientation candidate",
            $"Every retained source heavy atom mapped within 0.05 Å under one proper rigid transform; " +
            $"OPM midplane {midplane:G6} Å and hydrophobic thickness {thickness:G6} Å; " +
            $"chosen outer leaflet envelope ±{framePolicy.OuterLeafletEnvelopeAngstrom:G6} Å; " +
            $"tilt {reference.TiltDegrees?.ToString("G6") ?? "unavailable"}°",
            $"Prepared protein {protein.Id}; selected membrane {membrane.Id}; " +
            $"OPM reference context {reference.MembraneContext}; source digest {reference.OrientedCoordinateSha256}",
            "The OPM membrane assumption describes the reference; this candidate keeps the selected prepared atoms and bonds.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(proposalId,
            protein.Id, membrane.Id, topologyKind,
            molecule, 0, 2 * framePolicy.OuterLeafletEnvelopeAngstrom, reference.TiltDegrees,
            physicalSide, biologicalSidedness, ImmutableArray<string>.Empty, evidence,
            ImmutableArray<string>.Empty, FramePolicyId: framePolicy.Id,
            FramePolicyVersion: framePolicy.Version));
    }

    private readonly record struct PdbPoint(double X, double Y, double Z);
    private readonly record struct PdbAtom(string Key, string Element, PdbPoint Position);

    private static bool ValidFramePolicy(PlacementSupportPolicy policy,
        ProteinTopologyKind topology, MembraneModel membrane) =>
        !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
        policy.OuterLeafletEnvelopeAngstrom is double envelope &&
        double.IsFinite(envelope) && envelope is > 0 and <= 80 &&
        policy.CoveredTopologyKinds.Contains(topology) &&
        membrane.Upper.Fractions.Concat(membrane.Lower.Fractions)
            .Where(fraction => fraction.Fraction > 0)
            .All(fraction => policy.CoveredSpeciesIds.Contains(fraction.SpeciesId));

    private static bool IsPdbAtom(string line) =>
        line.StartsWith("ATOM  ", StringComparison.Ordinal) ||
        line.StartsWith("HETATM", StringComparison.Ordinal);

    private static bool TryPdbPoint(string line, out PdbPoint point)
    {
        point = default;
        if (line.Length < 54 ||
            !double.TryParse(line.Substring(30, 8), NumberStyles.Float,
                CultureInfo.InvariantCulture, out var x) ||
            !double.TryParse(line.Substring(38, 8), NumberStyles.Float,
                CultureInfo.InvariantCulture, out var y) ||
            !double.TryParse(line.Substring(46, 8), NumberStyles.Float,
                CultureInfo.InvariantCulture, out var z) ||
            !double.IsFinite(x) || !double.IsFinite(y) || !double.IsFinite(z))
            return false;
        point = new PdbPoint(x, y, z);
        return true;
    }

    private static bool TryPdbAtom(string line, out PdbAtom atom)
    {
        atom = default;
        if (!IsPdbAtom(line) || line.Length < 54 ||
            !int.TryParse(line.Substring(22, 4), NumberStyles.Integer,
                CultureInfo.InvariantCulture, out var residue) ||
            !TryPdbPoint(line, out var point)) return false;
        var chain = line.Substring(21, 1).Trim();
        var insertion = line.Substring(26, 1).Trim();
        var name = line.Substring(12, 4).Trim();
        var element = line.Length >= 78 ? line.Substring(76, 2).Trim() :
            name.TrimStart('0', '1', '2', '3').FirstOrDefault().ToString();
        if (name.Length == 0 || element.Length == 0) return false;
        atom = new PdbAtom($"{chain}:{residue}:{insertion}:{name}", element, point);
        return true;
    }

    private static bool TryWriteShiftedPdb(string source, string output, double zShift,
        out string sha256)
    {
        sha256 = string.Empty;
        if (!double.IsFinite(zShift)) return false;
        try
        {
            var lines = File.ReadAllLines(source);
            var atomCount = 0;
            for (var index = 0; index < lines.Length; index++)
            {
                var line = lines[index];
                if (!IsPdbAtom(line)) continue;
                if (!TryPdbPoint(line, out var point) ||
                    !double.IsFinite(point.Z + zShift) ||
                    Math.Abs(point.Z + zShift) >= 10000) return false;
                lines[index] = line[..46] +
                    string.Format(CultureInfo.InvariantCulture, "{0,8:F3}", point.Z + zShift) +
                    line[54..];
                atomCount++;
            }
            if (atomCount == 0) return false;
            File.WriteAllLines(output, lines);
            using var stream = File.OpenRead(output);
            sha256 = Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
            return true;
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        { return false; }
    }

    private static PdbPoint Add(PdbPoint a, PdbPoint b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    private static PdbPoint Sub(PdbPoint a, PdbPoint b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    private static PdbPoint Scale(PdbPoint a, double factor) =>
        new(a.X * factor, a.Y * factor, a.Z * factor);
    private static double Dot(PdbPoint a, PdbPoint b) =>
        a.X * b.X + a.Y * b.Y + a.Z * b.Z;
    private static PdbPoint Cross(PdbPoint a, PdbPoint b) =>
        new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z,
            a.X * b.Y - a.Y * b.X);
    private static double Length(PdbPoint a) => Math.Sqrt(Dot(a, a));

    private static bool TryOpmAlignment(AssessedPreparedProtein protein,
        OpmReferenceRecord reference, out Func<PdbPoint, PdbPoint>? transform,
        out string issue)
    {
        transform = null;
        issue = "The OPM oriented coordinate atoms cannot be mapped as one rigid transform to the exact selected prepared construct.";
        if (protein.Intended.Source.Kind != SourceRouteKind.Rcsb ||
            string.IsNullOrWhiteSpace(reference.SourceUrl) ||
            !protein.Correspondence.Complete ||
            !string.Equals(protein.Correspondence.ResultId,
                protein.Molecule.CoordinateSha256, StringComparison.OrdinalIgnoreCase) ||
            protein.Correspondence.Atoms.IsDefault ||
            protein.Correspondence.Atoms.Length != protein.Molecule.AtomCount ||
            !ArtifactMatches(new WorkerArtifact("preparedPdb", protein.Molecule.CoordinatePath,
                protein.Molecule.CoordinateSha256)) ||
            !ArtifactMatches(new WorkerArtifact("opmPdb", reference.OrientedCoordinatePath,
                reference.OrientedCoordinateSha256)))
        {
            issue = "The identified prepared or OPM coordinate bytes and correspondence are unavailable or changed.";
            return false;
        }
        string[] preparedLines;
        string[] opmLines;
        try
        {
            preparedLines = File.ReadAllLines(protein.Molecule.CoordinatePath);
            opmLines = File.ReadAllLines(reference.OrientedCoordinatePath);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        { return false; }
        var prepared = new List<PdbAtom>();
        foreach (var line in preparedLines.Where(IsPdbAtom))
        {
            if (!TryPdbAtom(line, out var atom)) return false;
            prepared.Add(atom);
        }
        if (prepared.Count != protein.Molecule.AtomCount) return false;
        var hasModels = opmLines.Any(line => line.StartsWith("MODEL ", StringComparison.Ordinal));
        if (!hasModels && protein.Intended.ModelIndex != 0) return false;
        var selectedModel = !hasModels;
        var opm = new List<PdbAtom>();
        var markers = new List<double>();
        foreach (var line in opmLines)
        {
            if (line.StartsWith("MODEL ", StringComparison.Ordinal))
            {
                selectedModel = int.TryParse(line[6..].Trim(), out var model) &&
                    model == protein.Intended.ModelIndex + 1;
                continue;
            }
            if (line.StartsWith("ENDMDL", StringComparison.Ordinal))
            {
                selectedModel = false;
                continue;
            }
            if (!selectedModel || !IsPdbAtom(line)) continue;
            if (line.Length >= 20 && line.Substring(17, 3) == "DUM")
            {
                if (!TryPdbPoint(line, out var marker)) return false;
                markers.Add(marker.Z);
                continue;
            }
            if (!TryPdbAtom(line, out var atom)) return false;
            if (!string.Equals(atom.Element, "H", StringComparison.OrdinalIgnoreCase))
                opm.Add(atom);
        }
        if (markers.Count < 2 || reference.MidplaneAngstrom is not double midplane ||
            reference.HydrophobicThicknessAngstrom is not double thickness ||
            !double.IsFinite(midplane) || !double.IsFinite(thickness) || thickness <= 0 ||
            Math.Abs(midplane - (markers.Min() + markers.Max()) / 2) > 0.05 ||
            Math.Abs(thickness - (markers.Max() - markers.Min())) > 0.05)
            return false;
        var opmByKey = opm.GroupBy(atom => atom.Key).ToDictionary(group => group.Key,
            group => group.ToArray(), StringComparer.Ordinal);
        var matched = new List<(PdbPoint Prepared, PdbPoint Oriented)>();
        var used = new HashSet<string>(StringComparer.Ordinal);
        foreach (var mapping in protein.Correspondence.Atoms)
        {
            if (mapping.ResultAtomIndex < 0 || mapping.ResultAtomIndex >= prepared.Count)
                return false;
            var atom = prepared[mapping.ResultAtomIndex];
            if (mapping.SourceAtomId is null ||
                string.Equals(atom.Element, "H", StringComparison.OrdinalIgnoreCase))
                continue;
            var fields = mapping.SourceAtomId.Split(':');
            if (fields.Length != 6 || !int.TryParse(fields[0], out var model) ||
                model != protein.Intended.ModelIndex ||
                !int.TryParse(fields[3], out var residue)) return false;
            var key = $"{fields[1]}:{residue}:{fields[4]}:{fields[5]}";
            if (!used.Add(key) || !opmByKey.TryGetValue(key, out var candidates) ||
                candidates.Length != 1 ||
                !string.Equals(atom.Element, candidates[0].Element,
                    StringComparison.OrdinalIgnoreCase)) return false;
            matched.Add((atom.Position, candidates[0].Position));
        }
        if (matched.Count < 3 || opmByKey.Count != matched.Count) return false;
        var source0 = matched[0].Prepared;
        var target0 = matched[0].Oriented;
        var farthest = matched.Select((pair, index) =>
            (Index: index, Distance: Dot(Sub(pair.Prepared, source0),
                Sub(pair.Prepared, source0)))).MaxBy(item => item.Distance);
        if (farthest.Distance < 1) return false;
        var sourceAxis = Scale(Sub(matched[farthest.Index].Prepared, source0),
            1 / Math.Sqrt(farthest.Distance));
        var targetAxisVector = Sub(matched[farthest.Index].Oriented, target0);
        if (Math.Abs(Length(targetAxisVector) - Math.Sqrt(farthest.Distance)) > 0.05)
            return false;
        var targetAxis = Scale(targetAxisVector, 1 / Length(targetAxisVector));
        var third = matched.Select((pair, index) =>
        {
            var difference = Sub(pair.Prepared, source0);
            var perpendicular = Sub(difference, Scale(sourceAxis, Dot(difference, sourceAxis)));
            return (Index: index, Distance: Dot(perpendicular, perpendicular));
        }).MaxBy(item => item.Distance);
        if (third.Distance < 1) return false;
        var sourceSecondVector = Sub(matched[third.Index].Prepared, source0);
        var sourceSecond = Scale(Sub(sourceSecondVector,
            Scale(sourceAxis, Dot(sourceSecondVector, sourceAxis))),
            1 / Math.Sqrt(third.Distance));
        var targetSecondVector = Sub(matched[third.Index].Oriented, target0);
        var targetPerpendicular = Sub(targetSecondVector,
            Scale(targetAxis, Dot(targetSecondVector, targetAxis)));
        if (Math.Abs(Length(targetPerpendicular) - Math.Sqrt(third.Distance)) > 0.05)
            return false;
        var targetSecond = Scale(targetPerpendicular, 1 / Length(targetPerpendicular));
        var sourceThird = Cross(sourceAxis, sourceSecond);
        var targetThird = Cross(targetAxis, targetSecond);
        PdbPoint Apply(PdbPoint point)
        {
            var relative = Sub(point, source0);
            return Add(target0, Add(Scale(targetAxis, Dot(relative, sourceAxis)),
                Add(Scale(targetSecond, Dot(relative, sourceSecond)),
                    Scale(targetThird, Dot(relative, sourceThird)))));
        }
        if (matched.Any(pair => Length(Sub(Apply(pair.Prepared), pair.Oriented)) > 0.05))
            return false;
        transform = Apply;
        issue = string.Empty;
        return true;
    }

    public async Task<BoundaryOutcome<PlacementProposal>> ProposeWithPpmAsync(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        MembraneModel membrane,
        PlacementSupportPolicy framePolicy,
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
        if (!ValidFramePolicy(framePolicy, topologyKind, membrane) ||
            !((topologyKind == ProteinTopologyKind.MembraneSpanning && physicalSide == PlacementPhysicalSide.Both) ||
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
        var framedPath = Path.Combine(workingDirectory, $"ppm-framed-{proposalId}.pdb");
        if (!TryWriteShiftedPdb(oriented.Path, framedPath,
                -observed.MidplaneAngstrom.Value, out var framedSha256))
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The PPM orientation could not be placed in the chosen membrane frame without changing the construct.");
        var molecule = new MolecularArtifact(
            proposalId, framedPath, framedSha256, protein.Molecule.TopologyPath, null, null,
            protein.Molecule.AtomCount, null, protein.Molecule.TopologySha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), proposalId, "PPM " + ppmVersion,
            "PPM orientation candidate",
            $"PPM midplane {observed.MidplaneAngstrom.Value:G6} Å and thickness " +
            $"{observed.ThicknessAngstrom.Value:G6} Å; chosen outer leaflet envelope " +
            $"±{framePolicy.OuterLeafletEnvelopeAngstrom:G6} Å",
            $"Prepared protein {protein.Id}; assumed membrane {observed.AssumedMembrane}; " +
            $"frame policy {framePolicy.Id} version {framePolicy.Version}",
            "Implicit symmetric-membrane orientation is not itself support for the chosen explicit bilayer or biological sidedness.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(
            proposalId, protein.Id, membrane.Id, topologyKind, molecule,
            0, 2 * framePolicy.OuterLeafletEnvelopeAngstrom, observed.TiltDegrees,
            physicalSide, biologicalSidedness, ImmutableArray<string>.Empty,
            evidence, observed.InterpretationWarnings,
            FramePolicyId: framePolicy.Id, FramePolicyVersion: framePolicy.Version));
    }

    public async Task<BoundaryOutcome<PlacementProposal>> ProposeManualAsync(
        StudyRevision revision, AssessedPreparedProtein protein, MembraneModel membrane,
        IReadOnlyDictionary<string, MolecularRepresentation> frameSpecies,
        PlacementSupportPolicy framePolicy,
        PlacementStartingPosition startingPosition,
        double offsetXAngstrom, double offsetYAngstrom, double offsetZAngstrom,
        double rotationXDegrees, double rotationYDegrees, double rotationZDegrees,
        int maximumAtomCount, string workingDirectory, CancellationToken cancellationToken)
    {
        if (revision.Id != protein.StudyRevisionId ||
            revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Id || !Enum.IsDefined(startingPosition) ||
            !ValidFramePolicy(framePolicy, ProteinTopologyKind.Unclassified, membrane) ||
            !new[] { offsetXAngstrom, offsetYAngstrom, offsetZAngstrom,
                rotationXDegrees, rotationYDegrees, rotationZDegrees }.All(double.IsFinite) ||
            maximumAtomCount <= 0 || protein.Molecule.AtomCount > maximumAtomCount)
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "Choose a current prepared protein, membrane and finite starting transform within the atom limit.");
        var speciesIds = membrane.Upper.Fractions.Concat(membrane.Lower.Fractions)
            .Where(fraction => fraction.Fraction > 0).Select(fraction => fraction.SpeciesId)
            .Distinct(StringComparer.Ordinal).ToHashSet(StringComparer.Ordinal);
        var species = speciesIds.Where(frameSpecies.ContainsKey).Select(id => frameSpecies[id]).ToImmutableArray();
        if (species.Length != speciesIds.Count ||
            !species.Any(item => item.Category == "lipid") ||
            !ArtifactMatches(new WorkerArtifact("preparedPdb", protein.Molecule.CoordinatePath,
                protein.Molecule.CoordinateSha256)))
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The exact prepared coordinates or chosen membrane frame species are unavailable.");
        var request = new ScientificWorkRequest<ManualPlacementPayload>(Guid.NewGuid().ToString("N"),
            workingDirectory, new ManualPlacementPayload(revision.Id, protein.Id,
                protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256,
                protein.Molecule.AtomCount, species, startingPosition,
                offsetXAngstrom, offsetYAngstrom, offsetZAngstrom,
                rotationXDegrees, rotationYDegrees, rotationZDegrees, maximumAtomCount,
                framePolicy.OuterLeafletEnvelopeAngstrom));
        var result = await _worker.PlaceManualAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PlacementProposal>.Unavailable(result.FailureMessage ??
                "The user-defined position was not observed.");
        var observed = result.Observations;
        var oriented = result.Artifacts.FirstOrDefault(item => item.Role == "orientedPdb");
        if (oriented is null || !ArtifactMatches(oriented) ||
            observed.SourceAtomCount != protein.Molecule.AtomCount ||
            observed.OrientedAtomCount != protein.Molecule.AtomCount ||
            !new[] { observed.AppliedTranslationXAngstrom, observed.AppliedTranslationYAngstrom,
                observed.AppliedTranslationZAngstrom, observed.HeadgroupBoundaryAngstrom,
                observed.MaximumRigidDeviationAngstrom }.All(double.IsFinite) ||
            Math.Abs(observed.HeadgroupBoundaryAngstrom -
                framePolicy.OuterLeafletEnvelopeAngstrom!.Value) > 0.001 ||
            observed.MaximumRigidDeviationAngstrom is < 0 or > 0.001 ||
            !observed.GeometryWarnings.IsDefaultOrEmpty)
            return BoundaryOutcome<PlacementProposal>.Unavailable(
                "The positioned artifact did not preserve the complete exact construct and declared rigid transform.");
        var id = Guid.NewGuid().ToString("N");
        var transform = new PlacementTransform(startingPosition, offsetXAngstrom, offsetYAngstrom,
            offsetZAngstrom, rotationXDegrees, rotationYDegrees, rotationZDegrees,
            observed.AppliedTranslationXAngstrom, observed.AppliedTranslationYAngstrom,
            observed.AppliedTranslationZAngstrom, observed.HeadgroupBoundaryAngstrom);
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"), id,
            result.Provider?.Name ?? "local rigid transform", "User-defined starting position",
            $"{startingPosition}; translation ({transform.AppliedTranslationXAngstrom:G6}, " +
            $"{transform.AppliedTranslationYAngstrom:G6}, {transform.AppliedTranslationZAngstrom:G6}) Å; " +
            $"rotation X/Y/Z ({rotationXDegrees:G6}, {rotationYDegrees:G6}, {rotationZDegrees:G6})°; " +
            $"maximum written-coordinate deviation {observed.MaximumRigidDeviationAngstrom:G6} Å",
            $"Prepared protein {protein.Id}; intended membrane {membrane.Id}; " +
            $"declared outer leaflet envelope {framePolicy.OuterLeafletEnvelopeAngstrom:G6} Å; " +
            $"frame policy {framePolicy.Id} version {framePolicy.Version}; " +
            $"selected species {string.Join(", ", species.Select(item => item.SpeciesId))}",
            "This is an adjustable geometric starting position for construction.", EvidenceBearing.Context));
        var molecule = new MolecularArtifact(id, oriented.Path, oriented.Sha256,
            protein.Molecule.TopologyPath, null, null, protein.Molecule.AtomCount,
            null, protein.Molecule.TopologySha256);
        var side = startingPosition switch
        {
            PlacementStartingPosition.Upper => PlacementPhysicalSide.Upper,
            PlacementStartingPosition.Lower => PlacementPhysicalSide.Lower,
            _ => PlacementPhysicalSide.Both
        };
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(id, protein.Id,
            membrane.Id, ProteinTopologyKind.Unclassified, molecule, 0,
            2 * framePolicy.OuterLeafletEnvelopeAngstrom, null, side, null,
            ImmutableArray<string>.Empty, evidence, ImmutableArray<string>.Empty, transform,
            framePolicy.Id, framePolicy.Version));
    }

    public async Task<BoundaryOutcome<PlacementProposal>> ReviseProposalAsync(
        StudyRevision revision,
        PlacementProposal source,
        double depthShiftAngstrom,
        double tiltAboutXDegrees,
        double tiltAboutYDegrees,
        double rotationAboutNormalDegrees,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        var adjustment = new[] { depthShiftAngstrom, tiltAboutXDegrees, tiltAboutYDegrees,
            rotationAboutNormalDegrees };
        if (revision.Membrane?.Id != source.MembraneModelId)
            return BoundaryOutcome<PlacementProposal>.Unavailable("The placement correction needs the corresponding current membrane.");
        if (!adjustment.All(double.IsFinite) || !adjustment.Any(value => value != 0))
            return BoundaryOutcome<PlacementProposal>.Unavailable("Change at least one finite depth, tilt, or rotation value.");
        if (!ArtifactMatches(new WorkerArtifact("orientedPdb", source.OrientedProtein.CoordinatePath,
                source.OrientedProtein.CoordinateSha256)))
            return BoundaryOutcome<PlacementProposal>.Unavailable("The identified source placement structure is missing or changed.");
        var request = new ScientificWorkRequest<PlacementAdjustmentPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new PlacementAdjustmentPayload(revision.Id, source.Id, source.OrientedProtein.CoordinatePath,
                depthShiftAngstrom, tiltAboutXDegrees, tiltAboutYDegrees,
                rotationAboutNormalDegrees, source.OrientedProtein.CoordinateSha256));
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
        var adjustmentSummary = $"Depth {observed.AppliedDepthShiftAngstrom:G6} Å; " +
            $"tilt X {observed.AppliedTiltAboutXDegrees:G6}°, tilt Y {observed.AppliedTiltAboutYDegrees:G6}°; " +
            $"normal rotation {observed.AppliedRotationAboutNormalDegrees:G6}°.";
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), proposalId, "researcher adjustment",
            "Observed rigid-body placement adjustment", adjustmentSummary,
            $"Derived from proposal {source.Id}",
            "The adjustment requires fresh assessment; earlier orientation energies and support do not transfer.",
            EvidenceBearing.Context));
        return BoundaryOutcome<PlacementProposal>.Success(new PlacementProposal(
            proposalId, source.PreparedProteinId, source.MembraneModelId,
            source.TopologyKind, artifact,
            source.MidplaneAngstrom,
            source.ThicknessAngstrom,
            // A scalar tilt cannot be updated by adding one Euler rotation. The
            // adjusted construct requires a fresh axis-based tilt observation.
            null,
            source.PhysicalSide, source.BiologicalSidedness, source.ContactingRegions,
            evidence, ImmutableArray.Create("Manual placement adjustment; new support assessment required."),
            FramePolicyId: source.FramePolicyId,
            FramePolicyVersion: source.FramePolicyVersion));
    }

    public async Task<BoundaryOutcome<PlacementMeasurementReport>> MeasureAgainstMembraneAsync(
        StudyRevision revision,
        AssessedPreparedProtein protein,
        MembraneModel membrane,
        PlacementProposal proposal,
        PlacementSupportPolicy? policy,
        PlacementStructuralWitness? witness,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        if (revision.Id != protein.StudyRevisionId ||
            revision.Membrane?.Id != membrane.Id || proposal.PreparedProteinId != protein.Id ||
            proposal.MembraneModelId != membrane.Id || proposal.MidplaneAngstrom is not double midplane ||
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
        var prediction = await SummarizePredictionForPlacementAsync(revision, protein, proposal,
            observed.Residues, witness, policy, workingDirectory, cancellationToken);
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            proposal.Id, result.Provider?.Name ?? "local structural measurement",
            "Measured placement geometry",
            string.Join("; ", measurements.Select(value => $"{value.Name}={value.Value:G6} {value.Unit}")),
            $"Prepared protein {protein.Id}; membrane {membrane.Id}; proposal {proposal.Id}; " +
            $"observed residue count {observed.Residues.Length}",
            "Atom positions relative to the intended slab were measured; no biological topology judgment was made.",
            EvidenceBearing.Context));
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
        MembraneModel? membrane,
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
        var reason = "The exact position has not been checked.";
        var standing = AssessmentStanding.NotEstablished;

        if (proposal.PreparedProteinId != protein.Id || protein.StudyRevisionId != revision.Id ||
            revision.Membrane?.Id != proposal.MembraneModelId ||
            membrane?.Id != proposal.MembraneModelId)
            reason = "Prepared protein, membrane, placement and study revision do not correspond.";
        else if (policy is null || !ValidFramePolicy(policy, proposal.TopologyKind, membrane!) ||
                 proposal.FramePolicyId != policy.Id ||
                 proposal.FramePolicyVersion != policy.Version ||
                 proposal.MidplaneAngstrom is not double midplane ||
                 Math.Abs(midplane) > 0.001 ||
                 proposal.ThicknessAngstrom is not double thickness ||
                 Math.Abs(thickness - 2 * policy.OuterLeafletEnvelopeAngstrom!.Value) > 0.001)
            reason = "The proposed position no longer corresponds to the chosen membrane frame policy.";
        else if (findings.Any(finding => finding.Material && finding.Disposition == FindingDisposition.Disqualifies &&
                     finding.SubjectId == proposal.Id) ||
                 allEvidence.Any(evidence => evidence.SubjectId == proposal.Id && evidence.Bearing == EvidenceBearing.Contradicts))
        {
            standing = AssessmentStanding.Unsupported;
            reason = "A material, attributable finding contradicts the proposed relationship.";
        }
        else if (measurement is null || measurement.ProposalId != proposal.Id ||
                 proposal.OrientedProtein.AtomCount != protein.Molecule.AtomCount ||
                 proposal.OrientedProtein.TopologySha256 != protein.Molecule.TopologySha256 ||
                 measurement.Residues.IsDefaultOrEmpty || !measurement.Limitations.IsDefaultOrEmpty ||
                 !measurement.Evidence.Any(evidence => evidence.SubjectId == proposal.Id &&
                     evidence.Method == "Measured placement geometry" && evidence.Bearing == EvidenceBearing.Context) ||
                 !new[] { "atomsWithinCore", "atomsAboveCore", "atomsBelowCore", "proteinSpan" }
                     .All(name => measurement.Measurements.Any(item => item.Name == name &&
                         double.IsFinite(item.Value))) ||
                 allEvidence.Any(evidence => evidence.SubjectId != proposal.Id))
            reason = "The complete placed construct or its membrane-frame measurement could not be checked.";
        else if (findings.Any(finding => finding.Material && finding.SubjectId == proposal.Id &&
                     finding.Disposition == FindingDisposition.Challenges))
            reason = "A later material finding requires rechecking this exact position.";
        else
        {
            standing = AssessmentStanding.Supported;
            reason = "The exact construct and chosen membrane frame passed the technical position checks.";
        }

        return new AssessedProteinMembranePlacement(
            Guid.NewGuid().ToString("N"), revision.Id, proposal,
            standing, reason,
            findings.IsDefault ? ImmutableArray<ScientificFinding>.Empty : findings,
            DateTimeOffset.UtcNow);
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
