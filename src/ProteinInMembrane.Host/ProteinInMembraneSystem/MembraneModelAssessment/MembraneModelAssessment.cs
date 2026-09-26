using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.MembraneModelAssessment;

/// <summary>Assesses one chosen planar bilayer locally, not its protein-specific placement.</summary>
public sealed class MembraneModelAssessment
{
    private readonly IMembraneModelAssessmentWork _worker;

    public MembraneModelAssessment(IMembraneModelAssessmentWork worker) => _worker = worker;

    public async Task<BoundaryOutcome<AssessedMembraneModel>> AssessAsync(
        StudyRevision revision,
        MembraneModel chosen,
        IReadOnlyDictionary<string, MolecularRepresentation> catalogue,
        MembraneSupportPolicy? policy,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(revision);
        ArgumentNullException.ThrowIfNull(chosen);
        ArgumentNullException.ThrowIfNull(catalogue);
        if (revision.Membrane is not { } adopted || !SameChoice(adopted, chosen))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The selected bilayer does not belong to this study revision.");
        if (chosen.Upper.PhysicalSide != LeafletSide.Upper)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The upper leaflet has the wrong physical side identity.");
        if (chosen.Lower.PhysicalSide != LeafletSide.Lower)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The lower leaflet has the wrong physical side identity.");
        var upperProblem = LeafletProblem(chosen.Upper, "upper");
        if (upperProblem is not null)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(upperProblem);
        var lowerProblem = LeafletProblem(chosen.Lower, "lower");
        if (lowerProblem is not null)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(lowerProblem);
        var selectedSpecies = chosen.Upper.Fractions.Concat(chosen.Lower.Fractions)
            .Where(fraction => fraction.Fraction > 0)
            .Select(fraction => fraction.SpeciesId)
            .Distinct(StringComparer.Ordinal)
            .ToImmutableArray();
        var unrepresented = selectedSpecies.Where(species => !catalogue.ContainsKey(species)).ToArray();
        if (unrepresented.Length > 0)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                $"Selected species {string.Join(", ", unrepresented)} has no qualified molecular representation.");
        if (!HasPositiveLipid(chosen.Upper, catalogue))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The upper physical leaflet has no positive-fraction lipid species; a sterol-only leaflet cannot establish a planar bilayer.");
        if (!HasPositiveLipid(chosen.Lower, catalogue))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The lower physical leaflet has no positive-fraction lipid species; a sterol-only leaflet cannot establish a planar bilayer.");
        if (policy is null || string.IsNullOrWhiteSpace(policy.Id) || string.IsNullOrWhiteSpace(policy.Version) ||
            policy.EvidenceReferences.IsDefaultOrEmpty || !policy.CoversArbitraryCoherentFractions)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                $"No evidence-backed membrane-local policy covers selected species {string.Join(", ", selectedSpecies)} and this fraction class.");
        var outsidePolicy = selectedSpecies.Where(species => !policy.CoveredSpeciesIds.Contains(species)).ToArray();
        if (outsidePolicy.Length > 0)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                $"Selected species {string.Join(", ", outsidePolicy)} is outside the declared membrane-policy scope.");
        if ((!policy.AllowsMixedLeaflets &&
             (chosen.Upper.Fractions.Count(fraction => fraction.Fraction > 0) > 1 ||
              chosen.Lower.Fractions.Count(fraction => fraction.Fraction > 0) > 1)) ||
            (!policy.AllowsAsymmetricLeaflets && !SameFractions(chosen.Upper, chosen.Lower)))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The chosen mixture or asymmetry is outside the declared membrane-policy scope.");

        var representations = selectedSpecies.Select(species => catalogue[species]).ToImmutableArray();
        if (selectedSpecies.Where((species, index) => representations[index].SpeciesId != species).Any() ||
            representations.Any(rep =>
            rep.AtomCount <= 0 || string.IsNullOrWhiteSpace(rep.ChemistryId) ||
            string.IsNullOrWhiteSpace(rep.TemplatePath) || string.IsNullOrWhiteSpace(rep.TemplateSha256) ||
            string.IsNullOrWhiteSpace(rep.CoordinateTemplatePath) || string.IsNullOrWhiteSpace(rep.CoordinateTemplateSha256) ||
            rep.ForceFieldFamily != "Lipid21" || string.IsNullOrWhiteSpace(rep.ForceFieldVersion) ||
            rep.StereoChecks.IsDefaultOrEmpty))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The catalogue does not identify exact Lipid21 coordinates, stereochemistry and parameter templates for each selected species.");
        var assetGroups = representations.GroupBy(rep => rep.TemplatePath, StringComparer.Ordinal).ToArray();
        if (assetGroups.Any(group => group.Select(rep => (rep.TemplateSha256, rep.ForceFieldVersion,
                rep.ForceFieldFamily)).Distinct().Count() != 1))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("One molecular parameter path has conflicting identified bytes or version.");
        var distinctAssets = assetGroups.Select(group =>
        {
            var rep = group.First();
            return new ForceFieldAsset(rep.SpeciesId, rep.ForceFieldVersion,
                rep.ForceFieldFamily, rep.TemplatePath, rep.TemplateSha256);
        }).ToImmutableArray();

        var representedUpper = chosen.Upper with
        { Fractions = chosen.Upper.Fractions.Where(item => item.Fraction > 0).ToImmutableArray() };
        var representedLower = chosen.Lower with
        { Fractions = chosen.Lower.Fractions.Where(item => item.Fraction > 0).ToImmutableArray() };
        var request = new ScientificWorkRequest<MembraneAssessmentPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new MembraneAssessmentPayload(revision.Id, chosen.Id,
                representedUpper, representedLower, representations,
                distinctAssets));
        var result = await _worker.AssessMembraneAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(result.FailureMessage ?? "Exact membrane representation was not observed.");

        var observed = result.Observations;
        if (observed.Species.IsDefault)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The worker did not report membrane species observations.");
        if (observed.Species.Length != representations.Length ||
            observed.Species.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
                representations.Length ||
            observed.Species.Any(item => !selectedSpecies.Contains(item.SpeciesId)))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                $"The worker reported membrane species {string.Join(", ", observed.Species.Select(item => item.SpeciesId))} " +
                $"instead of one observation for each exact selected species {string.Join(", ", selectedSpecies)}.");
        foreach (var item in observed.Species)
        {
            var representation = catalogue[item.SpeciesId];
            if (item.ChemistryId != representation.ChemistryId)
                return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                    $"Selected species {item.SpeciesId} has a different observed chemical identity than its qualified representation.");
            if (item.CoordinateAtomCount != item.ParameterAtomCount ||
                item.CoordinateAtomCount != representation.AtomCount)
                return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                    $"Selected species {item.SpeciesId} has incompatible coordinate and parameter atom counts.");
            if (item.Warnings.IsDefault)
                return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                    $"Selected species {item.SpeciesId} has no complete coordinate, bond and stereochemical observation.");
            if (!item.AtomIdentityAndBondMatch || !item.Warnings.IsEmpty)
                return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                    $"Selected species {item.SpeciesId} does not establish exact coordinate, bond, stereochemical and parameter correspondence." +
                    (item.Warnings.IsEmpty ? string.Empty : $" {string.Join(" ", item.Warnings)}"));
        }
        if (observed.CombinationWarnings.IsDefault)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                "The membrane-local combination has no complete parameterization observation.");
        if (observed.CombinedParameterizationObserved != true)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                observed.CombinedParameterizationObserved == false
                    ? "The selected membrane-local combination failed combined parameterization." +
                      (observed.CombinationWarnings.IsEmpty ? string.Empty :
                          $" {string.Join(" ", observed.CombinationWarnings)}")
                    : "The selected membrane-local combination was not observed for parameterization.");
        if (!observed.CombinationWarnings.IsEmpty)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(
                $"The selected membrane-local combination has unresolved parameterization warnings. {string.Join(" ", observed.CombinationWarnings)}");

        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), chosen.Id, result.Provider?.Name ?? "local scientific worker",
            "Exact membrane representation check",
            $"{representations.Length} selected molecular species matched their exact atom and bond identities, " +
            "declared stereochemical geometry, parameter templates and combined membrane-local parameterization",
            $"Membrane model {chosen.Id}; policy {policy.Id} version {policy.Version}",
            "Membrane-local support only; protein placement and complete-system parameterization remain separate.",
            EvidenceBearing.Supports));
        return BoundaryOutcome<AssessedMembraneModel>.Success(new AssessedMembraneModel(
            Guid.NewGuid().ToString("N"), revision.Id, chosen, representations, evidence,
            policy.Limitations.IsDefault ? ImmutableArray<string>.Empty : policy.Limitations,
            policy.Id, policy.Version));
    }

    private static string? LeafletProblem(LeafletComposition leaflet, string side)
    {
        if (leaflet.Fractions.IsDefaultOrEmpty)
            return $"The {side} physical leaflet needs identified species fractions.";
        if (leaflet.Fractions.Any(fraction => string.IsNullOrWhiteSpace(fraction.SpeciesId)) ||
            leaflet.Fractions.GroupBy(fraction => fraction.SpeciesId, StringComparer.Ordinal)
                .Any(group => group.Count() > 1))
            return $"The {side} physical leaflet needs unambiguous species; an identity is absent or repeated.";
        if (leaflet.Fractions.Any(fraction => !double.IsFinite(fraction.Fraction)))
            return $"The {side} physical leaflet has a nonfinite fraction.";
        if (leaflet.Fractions.Any(fraction => fraction.Fraction < 0))
            return $"The {side} physical leaflet has a negative fraction.";
        if (Math.Abs(leaflet.Fractions.Sum(fraction => fraction.Fraction) - 1.0) > 1e-9)
            return $"The {side} physical leaflet fractions do not sum to one.";
        return null;
    }

    private static bool HasPositiveLipid(LeafletComposition leaflet,
        IReadOnlyDictionary<string, MolecularRepresentation> catalogue) =>
        leaflet.Fractions.Any(fraction => fraction.Fraction > 0 &&
            catalogue[fraction.SpeciesId].Category == "lipid");

    private static bool SameFractions(LeafletComposition left, LeafletComposition right)
    {
        var presentLeft = left.Fractions.Where(item => item.Fraction > 0).ToArray();
        var presentRight = right.Fractions.Where(item => item.Fraction > 0).ToArray();
        return presentLeft.Length == presentRight.Length && presentLeft.All(item =>
            presentRight.Any(other => other.SpeciesId == item.SpeciesId &&
                Math.Abs(other.Fraction - item.Fraction) <= 1e-9));
    }

    private static bool SameChoice(MembraneModel adopted, MembraneModel chosen)
    {
        return adopted.Id == chosen.Id && adopted.Conditions == chosen.Conditions &&
            adopted.ScientificPurpose == chosen.ScientificPurpose &&
            adopted.Upper.PhysicalSide == chosen.Upper.PhysicalSide &&
            adopted.Lower.PhysicalSide == chosen.Lower.PhysicalSide &&
            SameDeclaredFractions(adopted.Upper.Fractions, chosen.Upper.Fractions) &&
            SameDeclaredFractions(adopted.Lower.Fractions, chosen.Lower.Fractions);
    }

    private static bool SameDeclaredFractions(ImmutableArray<LipidFraction> first,
        ImmutableArray<LipidFraction> second) =>
        !first.IsDefault && !second.IsDefault && first.Length == second.Length &&
        first.OrderBy(item => item.SpeciesId, StringComparer.Ordinal).ThenBy(item => item.Fraction)
            .SequenceEqual(second.OrderBy(item => item.SpeciesId, StringComparer.Ordinal)
                .ThenBy(item => item.Fraction));
}
