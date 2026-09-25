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
        if (revision.Membrane?.Id != chosen.Id)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The selected bilayer does not belong to this study revision.");
        if (chosen.Upper.PhysicalSide != LeafletSide.Upper || chosen.Lower.PhysicalSide != LeafletSide.Lower ||
            !Coherent(chosen.Upper) || !Coherent(chosen.Lower))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("Both physical leaflets need finite nonnegative fractions summing to one and unambiguous species.");
        if (policy is null || string.IsNullOrWhiteSpace(policy.Id) || string.IsNullOrWhiteSpace(policy.Version) ||
            policy.EvidenceReferences.IsDefaultOrEmpty || !policy.CoversArbitraryCoherentFractions)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("No evidence-backed membrane-local policy covers this fraction class.");
        if ((!policy.AllowsMixedLeaflets &&
             (chosen.Upper.Fractions.Count(fraction => fraction.Fraction > 0) > 1 ||
              chosen.Lower.Fractions.Count(fraction => fraction.Fraction > 0) > 1)) ||
            (!policy.AllowsAsymmetricLeaflets && !SameFractions(chosen.Upper, chosen.Lower)))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The chosen mixture or asymmetry is outside the declared membrane-policy scope.");

        var selectedSpecies = chosen.Upper.Fractions.Concat(chosen.Lower.Fractions)
            .Where(fraction => fraction.Fraction > 0)
            .Select(fraction => fraction.SpeciesId)
            .Distinct(StringComparer.Ordinal)
            .ToImmutableArray();
        if (selectedSpecies.Any(species => !policy.CoveredSpeciesIds.Contains(species) || !catalogue.ContainsKey(species)))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("A selected lipid or sterol has no covered, qualified molecular representation.");
        var representations = selectedSpecies.Select(species => catalogue[species]).ToImmutableArray();
        if (selectedSpecies.Where((species, index) => representations[index].SpeciesId != species).Any() ||
            representations.Any(rep =>
            rep.AtomCount <= 0 || string.IsNullOrWhiteSpace(rep.ChemistryId) ||
            string.IsNullOrWhiteSpace(rep.TemplatePath) || string.IsNullOrWhiteSpace(rep.TemplateSha256) ||
            string.IsNullOrWhiteSpace(rep.CoordinateTemplatePath) || string.IsNullOrWhiteSpace(rep.CoordinateTemplateSha256) ||
            rep.ForceFieldFamily != "Lipid21" || string.IsNullOrWhiteSpace(rep.ForceFieldVersion)))
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The catalogue does not identify exact Lipid21 coordinates and parameter templates for each selected species.");
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

        var request = new ScientificWorkRequest<MembraneAssessmentPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new MembraneAssessmentPayload(revision.Id, chosen.Id,
                chosen.Upper, chosen.Lower, representations,
                distinctAssets));
        var result = await _worker.AssessMembraneAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable(result.FailureMessage ?? "Exact membrane representation was not observed.");

        var observed = result.Observations;
        if (observed.Species.Length != representations.Length ||
            observed.Species.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
                representations.Length ||
            observed.Species.Any(item => !selectedSpecies.Contains(item.SpeciesId) ||
                item.ChemistryId != catalogue[item.SpeciesId].ChemistryId ||
                !item.AtomIdentityAndBondMatch || item.CoordinateAtomCount != item.ParameterAtomCount ||
                item.CoordinateAtomCount != catalogue[item.SpeciesId].AtomCount || !item.Warnings.IsDefaultOrEmpty) ||
            observed.CombinedParameterizationObserved != true || !observed.CombinationWarnings.IsDefaultOrEmpty)
            return BoundaryOutcome<AssessedMembraneModel>.Unavailable("The observed coordinates, templates or membrane-local combination do not establish exact support.");

        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), chosen.Id, result.Provider?.Name ?? "local scientific worker",
            "Exact membrane representation check",
            $"{representations.Length} selected molecular species matched their coordinate and parameter templates",
            $"Membrane model {chosen.Id}; policy {policy.Id} version {policy.Version}",
            "Membrane-local support only; protein placement and complete-system parameterization remain separate.",
            EvidenceBearing.Supports));
        return BoundaryOutcome<AssessedMembraneModel>.Success(new AssessedMembraneModel(
            Guid.NewGuid().ToString("N"), revision.Id, chosen, representations, evidence,
            policy.Limitations.IsDefault ? ImmutableArray<string>.Empty : policy.Limitations));
    }

    private static bool Coherent(LeafletComposition leaflet)
    {
        if (leaflet.Fractions.IsDefaultOrEmpty ||
            leaflet.Fractions.GroupBy(fraction => fraction.SpeciesId, StringComparer.Ordinal).Any(group => group.Count() > 1) ||
            leaflet.Fractions.Any(fraction => string.IsNullOrWhiteSpace(fraction.SpeciesId) ||
                !double.IsFinite(fraction.Fraction) || fraction.Fraction < 0))
            return false;
        return Math.Abs(leaflet.Fractions.Sum(fraction => fraction.Fraction) - 1.0) <= 1e-9;
    }

    private static bool SameFractions(LeafletComposition left, LeafletComposition right)
        => left.Fractions.Length == right.Fractions.Length && left.Fractions.All(item =>
            right.Fractions.Any(other => other.SpeciesId == item.SpeciesId &&
                Math.Abs(other.Fraction - item.Fraction) <= 1e-9));
}
