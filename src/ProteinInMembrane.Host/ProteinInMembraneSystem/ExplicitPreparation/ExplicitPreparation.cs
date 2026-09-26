using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation;

/// <summary>Coordinates one identified attempt from assessed inputs to an eligible explicit system.</summary>
public sealed class ExplicitPreparation
{
    private const double MoleculesPerMolarAngstromCubed = 6.02214076e-4;
    private readonly IExplicitConstructionWork _worker;
    private readonly Minimization _minimization;
    private readonly OptionalEquilibrationProcedure _equilibration;

    public ExplicitPreparation(IExplicitConstructionWork worker,
        IMinimizationWork minimizationWork, IOptionalEquilibrationWork equilibrationWork)
    {
        _worker = worker;
        _minimization = new Minimization(minimizationWork);
        _equilibration = new OptionalEquilibrationProcedure(equilibrationWork);
    }

    public Task<StageOperationResult> MinimizeAsync(ConstructedExplicitSystem source,
        ApplicablePreparationPolicy policy, string workingDirectory, IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
        => _minimization.RunAsync(source, policy, workingDirectory, progress, cancellationToken);

    public Task<StageOperationResult> RequestOptionalEquilibrationAsync(CompletedStage minimized,
        ApplicablePreparationPolicy policy, string workingDirectory, IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
        => _equilibration.RunAsync(minimized, policy, workingDirectory, progress, cancellationToken);

    public static bool PolicyScopeMatches(PreparationPolicyScope? scope, StudyRevision revision,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane, PlacementProposal proposal)
    {
        if (scope is null || !HashLike(scope.SourceCoordinateSha256) ||
            !HashLike(scope.PreparedBondGraphSha256) || !HashLike(scope.ResidueVariantsSha256) ||
            scope.SourceModelIndex < 0 || string.IsNullOrWhiteSpace(scope.ChemicalStatePolicyId) ||
            string.IsNullOrWhiteSpace(scope.ChemicalStatePolicyVersion) ||
            string.IsNullOrWhiteSpace(scope.ProteinStructuralPolicyId) ||
            string.IsNullOrWhiteSpace(scope.ProteinStructuralPolicyVersion) ||
            string.IsNullOrWhiteSpace(scope.MembraneSupportPolicyId) ||
            string.IsNullOrWhiteSpace(scope.MembraneSupportPolicyVersion) ||
            scope.ChainCopies.IsDefaultOrEmpty || scope.RetainedPartnerSourceIds.IsDefault ||
            scope.ChainCopies.Distinct().Count() != scope.ChainCopies.Length ||
            scope.RetainedPartnerSourceIds.Any(string.IsNullOrWhiteSpace) ||
            scope.RetainedPartnerSourceIds.Distinct(StringComparer.Ordinal).Count() !=
                scope.RetainedPartnerSourceIds.Length ||
            !Enum.IsDefined(scope.TopologyKind) ||
            scope.Upper?.PhysicalSide != LeafletSide.Upper ||
            scope.Lower?.PhysicalSide != LeafletSide.Lower ||
            !SameLeaflet(scope.Upper, membrane.Intended.Upper) ||
            !SameLeaflet(scope.Lower, membrane.Intended.Lower))
            return false;
        var intended = protein.Intended;
        if (intended.Partners.IsDefault) return false;
        var retained = intended.Partners.Where(item => item.Retain).Select(item => item.SourceId).ToArray();
        return scope.SourceCoordinateSha256.Equals(intended.Source.Sha256, StringComparison.OrdinalIgnoreCase) &&
            scope.SourceModelIndex == intended.ModelIndex &&
            scope.BiologicalAssemblyId == intended.BiologicalAssemblyId &&
            !intended.Chains.IsDefaultOrEmpty && intended.Chains.Distinct().Count() == intended.Chains.Length &&
            scope.ChainCopies.ToHashSet().SetEquals(intended.Chains) &&
            retained.Distinct(StringComparer.Ordinal).Count() == retained.Length &&
            scope.RetainedPartnerSourceIds.ToHashSet(StringComparer.Ordinal).SetEquals(retained) &&
            scope.PreparedBondGraphSha256.Equals(protein.Molecule.TopologySha256,
                StringComparison.OrdinalIgnoreCase) &&
            scope.ChemicalStatePolicyId == protein.ChemicalStatePolicyId &&
            scope.ChemicalStatePolicyVersion == protein.ChemicalStatePolicyVersion &&
            scope.ProteinStructuralPolicyId == protein.StructuralAssessmentPolicyId &&
            scope.ProteinStructuralPolicyVersion == protein.StructuralAssessmentPolicyVersion &&
            scope.MembraneSupportPolicyId == membrane.PolicyId &&
            scope.MembraneSupportPolicyVersion == membrane.PolicyVersion &&
            scope.ResidueVariantsSha256.Equals(PreparationPolicyFingerprint.ComputeResidueVariants(protein.ResidueVariants),
                StringComparison.OrdinalIgnoreCase) &&
            scope.Conditions == revision.Conditions && scope.Conditions == membrane.Intended.Conditions &&
            scope.TopologyKind == proposal.TopologyKind &&
            proposal.PreparedProteinId == protein.Id && proposal.MembraneModelId == membrane.Intended.Id &&
            revision.IntendedProtein?.Id == intended.Id && revision.Membrane?.Id == membrane.Intended.Id;
    }

    private static bool SameLeaflet(LeafletComposition? declared, LeafletComposition actual)
    {
        if (declared is null || declared.PhysicalSide != actual.PhysicalSide ||
            declared.Fractions.IsDefaultOrEmpty || actual.Fractions.IsDefaultOrEmpty ||
            declared.Fractions.Length != actual.Fractions.Length ||
            declared.Fractions.Any(item => string.IsNullOrWhiteSpace(item.SpeciesId) ||
                !double.IsFinite(item.Fraction) || item.Fraction < 0) ||
            actual.Fractions.Any(item => string.IsNullOrWhiteSpace(item.SpeciesId) ||
                !double.IsFinite(item.Fraction) || item.Fraction < 0) ||
            declared.Fractions.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
                declared.Fractions.Length ||
            actual.Fractions.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
                actual.Fractions.Length ||
            Math.Abs(declared.Fractions.Sum(item => item.Fraction) - 1) > 1e-9 ||
            Math.Abs(actual.Fractions.Sum(item => item.Fraction) - 1) > 1e-9)
            return false;
        return declared.Fractions.OrderBy(item => item.SpeciesId, StringComparer.Ordinal)
            .Zip(actual.Fractions.OrderBy(item => item.SpeciesId, StringComparer.Ordinal))
            .All(pair => pair.First.SpeciesId == pair.Second.SpeciesId &&
                Math.Abs(pair.First.Fraction - pair.Second.Fraction) <= 1e-9);
    }

    private static bool HashLike(string? hash) => hash is { Length: 64 } && hash.All(Uri.IsHexDigit);


    public async Task<PreparationStartResult> StartAsync(
        PreparationAttempt attempt,
        StudyRevision revision,
        AssessedPreparedProtein protein,
        AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy,
        string workingDirectory,
        Action<PreparationAttempt>? onAccepted,
        IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
    {
        var construction = policy.Construction;
        if (attempt.StudyRevisionId != revision.Id || attempt.ProteinId != protein.Id ||
            attempt.MembraneId != membrane.Id || attempt.PlacementId != placement.Id ||
            !PreparationPolicyFingerprint.Matches(attempt, policy) ||
            revision.Id != protein.StudyRevisionId || revision.Id != membrane.StudyRevisionId ||
            revision.Id != placement.StudyRevisionId || revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Intended.Id ||
            revision.AdoptedPlacementProposalId != placement.Proposal.Id ||
            placement.Standing != AssessmentStanding.Supported ||
            placement.Proposal.PreparedProteinId != protein.Id ||
            placement.Proposal.MembraneModelId != membrane.Intended.Id ||
            !PolicyScopeMatches(policy.Scope, revision, protein, membrane, placement.Proposal) ||
            !ValidConstructionPolicy(construction) || !ValidLocalStatePolicy(policy) ||
            !ValidStageProteinGeometryPolicy(policy) ||
            policy.EvidenceReferences.IsDefaultOrEmpty || policy.ForceFieldFiles.IsDefaultOrEmpty ||
            policy.MaximumMinimizationIterations <= 0 ||
            policy.FinalUnrestrainedRmsForceTargetKjMolNm != 10.0 ||
            policy.SystemSettings is not { NonbondedMethod: "PME", Constraints: "HBonds" } ||
            !double.IsFinite(policy.SystemSettings.NonbondedCutoffNanometers) ||
            policy.SystemSettings.NonbondedCutoffNanometers <= 0 ||
            !double.IsFinite(policy.SystemSettings.EwaldErrorTolerance) ||
            policy.SystemSettings.EwaldErrorTolerance is <= 0 or >= 1 ||
            construction.CoveredTopologyKinds.IsDefaultOrEmpty ||
            !construction.CoveredTopologyKinds.Contains(placement.Proposal.TopologyKind) ||
            membrane.SpeciesRepresentations.Length != 1 ||
            membrane.SpeciesRepresentations[0].SpeciesId != construction.LipidTypeArgument ||
            !PureSelectedLeaflet(membrane.Intended.Upper, construction.LipidTypeArgument) ||
            !PureSelectedLeaflet(membrane.Intended.Lower, construction.LipidTypeArgument) ||
            revision.Conditions.TargetNaClMolar != 0.15 ||
            placement.Proposal.MidplaneAngstrom is not double midplane ||
            !double.IsFinite(midplane))
            return NoAttempt("Corresponding supported inputs and an identified native construction policy are required.");

        var lipid = membrane.SpeciesRepresentations[0];
        var needed = new[] { lipid, policy.Water, policy.Sodium, policy.Chloride };
        if (policy.Water.SpeciesId != "HOH" || policy.Water.Category != "water" ||
            policy.Sodium.SpeciesId != "NA" || policy.Sodium.Category != "ion" ||
            policy.Chloride.SpeciesId != "CL" || policy.Chloride.Category != "ion" ||
            Math.Abs(lipid.NetChargeElementary) > 1e-8 ||
            Math.Abs(policy.Water.NetChargeElementary) > 1e-8 ||
            Math.Abs(policy.Sodium.NetChargeElementary - 1.0) > 1e-8 ||
            Math.Abs(policy.Chloride.NetChargeElementary + 1.0) > 1e-8 ||
            needed.Any(species => !construction.CoveredSpeciesIds.Contains(species.SpeciesId) ||
                species.AtomCount <= 0 || string.IsNullOrWhiteSpace(species.ForceFieldFamily) ||
                string.IsNullOrWhiteSpace(species.ForceFieldVersion) ||
                !Available(workingDirectory, species.TemplatePath) ||
                !Available(workingDirectory, species.CoordinateTemplatePath)) ||
            !File.Exists(protein.Molecule.CoordinatePath) ||
            !protein.Correspondence.Complete ||
            protein.Correspondence.ResultId != protein.Molecule.CoordinateSha256 ||
            protein.Correspondence.Atoms.Length != protein.Molecule.AtomCount ||
            !File.Exists(placement.Proposal.OrientedProtein.CoordinatePath) ||
            string.IsNullOrWhiteSpace(placement.Proposal.OrientedProtein.CoordinateSha256) ||
            !Available(workingDirectory, placement.Proposal.OrientedProtein.TopologyPath ?? "") ||
            !Available(workingDirectory, construction.NativePatchPath) ||
            !HashMatches(construction.NativePatchPath, construction.NativePatchSha256) ||
            (construction.NativePatchMode == "popc-62-109-deletion" &&
                !HashMatches(construction.NativeSourcePatchPath!, construction.NativeSourcePatchSha256!)) ||
            policy.ForceFieldFiles.Any(asset =>
                string.IsNullOrWhiteSpace(asset.Id) || string.IsNullOrWhiteSpace(asset.Version) ||
                string.IsNullOrWhiteSpace(asset.Family) || !HashMatches(asset.Path, asset.Sha256)) ||
            policy.ForceFieldFiles.GroupBy(asset => asset.Path, StringComparer.Ordinal)
                .Any(group => group.Select(asset => asset.Sha256)
                    .Distinct(StringComparer.OrdinalIgnoreCase).Count() != 1))
            return Refusal("A qualified molecular representation, native package patch, or exact parameter asset is unavailable before start.");

        onAccepted?.Invoke(attempt);
        progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
            "OpenMM is constructing one identified membrane, solvent and ion candidate."));
        var forceFieldFiles = policy.ForceFieldFiles.DistinctBy(asset => asset.Sha256,
            StringComparer.OrdinalIgnoreCase).ToImmutableArray();
        var request = new ScientificWorkRequest<ConstructionPayload>(Guid.NewGuid().ToString("N"), workingDirectory,
            new ConstructionPayload(revision.Id, attempt.Id,
                placement.Proposal.OrientedProtein.CoordinatePath,
                placement.Proposal.OrientedProtein.CoordinateSha256,
                protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256,
                protein.Correspondence,
                placement.Proposal.OrientedProtein.TopologyPath!,
                placement.Proposal.OrientedProtein.TopologySha256!,
                lipid, policy.Water, policy.Sodium, policy.Chloride,
                construction.NativePatchPath, construction.NativePatchSha256,
                construction.ProviderName, construction.ProviderVersion,
                construction.LipidTypeArgument, construction.PositiveIonArgument,
                construction.NegativeIonArgument, midplane / 10.0,
                construction.MinimumPaddingNanometers, revision.Conditions.TargetNaClMolar,
                forceFieldFiles, policy.SystemSettings, policy.LocalStateObservation,
                construction.MaximumAtomCount, construction.MaximumCellDimensionAngstrom,
                construction.NativePatchMode, construction.NativeSourcePatchPath,
                construction.NativeSourcePatchSha256, construction.RemovedNativeLipidResidueIds));
        using var bounded = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
        bounded.CancelAfter(TimeSpan.FromSeconds(construction.MaximumConstructionSeconds));
        try
        {
            var built = await _worker.ConstructSystemAsync(request, bounded.Token);
            if (cancellationToken.IsCancellationRequested)
                return Failed(attempt, WorkerResultStanding.Stopped, "Construction was stopped before a candidate was established.");
            if (bounded.IsCancellationRequested)
                return Failed(attempt, WorkerResultStanding.Failed,
                    "The native construction exceeded its declared execution bound.");
            if (built.RequestId != request.RequestId || built.StudyRevisionId != revision.Id ||
                built.AttemptId != attempt.Id || built.Standing != WorkerResultStanding.Observed ||
                built.Observations is null)
                return Failed(attempt, built.Standing,
                    built.FailureMessage ?? "The native candidate was not observed.");
            if (!TryDeriveNative(attempt, revision, protein, membrane, policy,
                    built.Observations, out var derivation, out var reason))
                return new PreparationStartResult(attempt, null, null,
                    State(attempt.Id, StageExecutionStanding.Failed, reason));
            var candidate = await AssessConstructedAsync(attempt, revision, protein, membrane,
                placement, policy, derivation!, built, cancellationToken);
            if (candidate is null)
                return new PreparationStartResult(attempt, derivation, null,
                    State(attempt.Id, StageExecutionStanding.Failed,
                        "Actual membership, periodic geometry, contacts, parameters or correspondence were not established."));
            return new PreparationStartResult(attempt, derivation, candidate,
                State(attempt.Id, StageExecutionStanding.ReadyForMinimization,
                    "The provider candidate is checked and ready for researcher review before minimization."));
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return Failed(attempt, WorkerResultStanding.Stopped,
                "Construction was stopped before a candidate was established.");
        }
        catch (OperationCanceledException) when (bounded.IsCancellationRequested)
        {
            return Failed(attempt, WorkerResultStanding.Failed,
                "The native construction exceeded its declared execution bound.");
        }
    }

    private static bool PureSelectedLeaflet(LeafletComposition leaflet, string species) =>
        leaflet.Fractions.Length == 1 && leaflet.Fractions[0].SpeciesId == species &&
        leaflet.Fractions[0].Fraction == 1.0;

    private static PreparationStartResult Refusal(string reason) =>
        new(null, null, null, State(string.Empty, StageExecutionStanding.ResourceRefused, reason));

    private static PreparationStartResult NoAttempt(string reason) =>
        new(null, null, null, State(string.Empty, StageExecutionStanding.Pending, reason));

    private static PreparationStartResult Failed(PreparationAttempt attempt, WorkerResultStanding standing,
        string reason, ConstructionDerivation? derivation = null) =>
        new(attempt, derivation, null, State(attempt.Id, standing switch
        {
            WorkerResultStanding.Stopped => StageExecutionStanding.Stopped,
            WorkerResultStanding.Unobserved => StageExecutionStanding.Unobserved,
            _ => StageExecutionStanding.Failed
        }, reason));

    private static StageExecutionState State(string attemptId, StageExecutionStanding standing, string message) =>
        new(attemptId, null, null, standing, message, null, DateTimeOffset.UtcNow);

    private static bool Available(string workingDirectory, string path) =>
        !string.IsNullOrWhiteSpace(path) &&
        File.Exists(Path.IsPathRooted(path) ? path : Path.Combine(workingDirectory, path));

    private static bool HashMatches(string path, string expected)
    {
        if (string.IsNullOrWhiteSpace(path) || expected is not { Length: 64 } ||
            !expected.All(Uri.IsHexDigit) || !File.Exists(path)) return false;
        try
        {
            using var stream = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(stream)).Equals(expected,
                StringComparison.OrdinalIgnoreCase);
        }
        catch (IOException) { return false; }
        catch (UnauthorizedAccessException) { return false; }
    }

    internal static bool ValidConstructionPolicy(ConstructionPolicy p) =>
        p.ProviderName == "OpenMM Modeller.addMembrane" &&
        !string.IsNullOrWhiteSpace(p.ProviderVersion) &&
        p.NativePatchSha256 is { Length: 64 } && p.NativePatchSha256.All(Uri.IsHexDigit) &&
        p.LipidTypeArgument is "DMPC" or "POPC" && p.PositiveIonArgument == "Na+" &&
        p.NegativeIonArgument == "Cl-" &&
        ((p.NativePatchMode is null or "installed") && p.NativeSourcePatchPath is null &&
            p.NativeSourcePatchSha256 is null && p.RemovedNativeLipidResidueIds is null ||
         p.NativePatchMode == "popc-62-109-deletion" && p.LipidTypeArgument == "POPC" &&
            !string.IsNullOrWhiteSpace(p.NativeSourcePatchPath) &&
            p.NativeSourcePatchSha256 is { Length: 64 } &&
            p.NativeSourcePatchSha256.All(Uri.IsHexDigit) &&
            p.RemovedNativeLipidResidueIds is { } removed && !removed.IsDefault &&
            removed.SequenceEqual(["62", "109"])) &&
        double.IsFinite(p.MinimumPaddingNanometers) && p.MinimumPaddingNanometers > 0 &&
        p.WaterMolarityForIonRounding == 55.4 &&
        p.MaximumAtomCount > 0 &&
        double.IsFinite(p.MaximumCellDimensionAngstrom) && p.MaximumCellDimensionAngstrom > 0 &&
        p.MaximumConstructionSeconds is > 0 and <= 86400 &&
        !string.IsNullOrWhiteSpace(p.ApproximationStatement);

    private static bool TryDeriveNative(PreparationAttempt attempt, StudyRevision revision,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        ApplicablePreparationPolicy policy, ConstructionObservations observed,
        out ConstructionDerivation? derivation, out string reason)
    {
        derivation = null;
        reason = "The same-invocation native counts and cell were not coherently observed.";
        var construction = policy.Construction;
        var cell = observed.ActualCellAngstrom;
        var gaps = observed.ProteinPeriodicImageGapsAngstrom;
        if (cell.Length != 3 || gaps.Length != 3 ||
            cell.Any(value => !double.IsFinite(value) || value <= 0 ||
                value > construction.MaximumCellDimensionAngstrom) ||
            gaps.Any(value => !double.IsFinite(value) ||
                value + policy.ExportCellLengthReadBackToleranceAngstrom <
                    20.0 * construction.MinimumPaddingNanometers) ||
            cell.Any(value => value <= 20.0 * policy.SystemSettings.NonbondedCutoffNanometers) ||
            observed.AtomCount <= protein.Molecule.AtomCount ||
            observed.AtomCount > construction.MaximumAtomCount ||
            observed.WaterCount <= 0 || observed.PositiveIonCount < 0 ||
            observed.NegativeIonCount < 0 ||
            observed.SpeciesCounts.IsDefaultOrEmpty ||
            !double.IsFinite(observed.ProteinNetChargeElementary) ||
            !double.IsFinite(observed.NetChargeElementary) ||
            Math.Abs(observed.NetChargeElementary) > 1e-5)
            return false;
        var lipidId = membrane.SpeciesRepresentations[0].SpeciesId;
        var counts = ImmutableArray.CreateBuilder<SpeciesCount>();
        foreach (var side in new[] { LeafletSide.Upper, LeafletSide.Lower })
        {
            var matches = observed.SpeciesCounts.Where(item =>
                item.Role == GeneratedComponentRoleKind.Lipid &&
                item.PhysicalSide == side && item.SpeciesId == lipidId).ToArray();
            if (matches.Length != 1 || matches[0].Count <= 0) return false;
            counts.Add(new SpeciesCount(side, lipidId, matches[0].Count, 1.0));
        }
        var generatedAtoms = (long)counts.Sum(item => item.Count) * membrane.SpeciesRepresentations[0].AtomCount +
            (long)observed.WaterCount * policy.Water.AtomCount +
            (long)observed.PositiveIonCount * policy.Sodium.AtomCount +
            (long)observed.NegativeIonCount * policy.Chloride.AtomCount;
        if (generatedAtoms + protein.Molecule.AtomCount != observed.AtomCount) return false;
        var charge = Math.Round(observed.ProteinNetChargeElementary);
        if (Math.Abs(observed.ProteinNetChargeElementary - charge) > 1e-5 ||
            Math.Abs(charge) > int.MaxValue / 2) return false;
        var neutralizers = (int)Math.Abs(charge);
        var expectedSodium = charge < 0 ? neutralizers : 0;
        var expectedChloride = charge > 0 ? neutralizers : 0;
        var saltPairs = Math.Min(observed.PositiveIonCount, observed.NegativeIonCount);
        if (observed.PositiveIonCount != saltPairs + expectedSodium ||
            observed.NegativeIonCount != saltPairs + expectedChloride) return false;
        var equivalent = (long)observed.WaterCount + observed.PositiveIonCount + observed.NegativeIonCount;
        if (equivalent <= neutralizers) return false;
        var nativeSaltPairs = Math.Floor(0.5 + (equivalent - neutralizers) *
            revision.Conditions.TargetNaClMolar / construction.WaterMolarityForIonRounding);
        if (!double.IsFinite(nativeSaltPairs) || nativeSaltPairs != saltPairs) return false;
        var estimatedAqueousVolume = equivalent /
            (construction.WaterMolarityForIonRounding * MoleculesPerMolarAngstromCubed);
        var estimatedMolar = construction.WaterMolarityForIonRounding * saltPairs / equivalent;
        if (!double.IsFinite(estimatedAqueousVolume) || estimatedAqueousVolume <= 0 ||
            !double.IsFinite(estimatedMolar)) return false;
        derivation = new ConstructionDerivation(attempt.Id, counts.ToImmutable(), cell,
            observed.WaterCount, observed.PositiveIonCount, observed.NegativeIonCount,
            observed.ProteinNetChargeElementary, revision.Conditions.TargetNaClMolar,
            estimatedMolar, estimatedAqueousVolume,
            ImmutableArray.Create(construction.ApproximationStatement,
                "OpenMM selected the finite cell and populations in this same construction invocation; water-equivalent salt molarity is an estimate, not an equilibrated concentration."),
            policy.Limitations);
        reason = string.Empty;
        return true;
    }

    private static async Task<ConstructedExplicitSystem?> AssessConstructedAsync(
        PreparationAttempt attempt, StudyRevision revision, AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy,
        ConstructionDerivation derivation, WorkerResult<ConstructionObservations> result,
        CancellationToken cancellationToken)
    {
        var observed = result.Observations!;
        if (result.Provider?.Name != policy.Construction.ProviderName ||
            result.Provider.Version != policy.Construction.ProviderVersion ||
            observed.NativePatchMode != (policy.Construction.NativePatchMode ?? "installed") ||
            !string.Equals(observed.NativeSourcePatchSha256,
                policy.Construction.NativeSourcePatchSha256, StringComparison.OrdinalIgnoreCase) ||
            !string.Equals(observed.NativePatchSha256, policy.Construction.NativePatchSha256,
                StringComparison.OrdinalIgnoreCase) ||
            observed.CorrespondedResultAtomCount != observed.AtomCount ||
            observed.WaterCount != derivation.WaterCount ||
            observed.PositiveIonCount != derivation.SodiumCount ||
            observed.NegativeIonCount != derivation.ChlorideCount ||
            observed.ActualCellAngstrom.Length != 3 ||
            observed.ActualCellAngstrom.Where((value, index) =>
                !double.IsFinite(value) || Math.Abs(value - derivation.CellAngstrom[index]) > 1e-6).Any() ||
            !observed.ProteinIdentityAndBondsPreserved ||
            !double.IsFinite(observed.MaximumProteinCoordinateDeviationAngstrom) ||
            observed.MaximumProteinCoordinateDeviationAngstrom < 0 ||
            observed.MaximumProteinCoordinateDeviationAngstrom > 1e-6 ||
            !double.IsFinite(observed.InitialPotentialEnergyKjMol) ||
            !observed.ContactWarnings.IsDefaultOrEmpty || !observed.GeometryWarnings.IsDefaultOrEmpty ||
            !observed.ParameterWarnings.IsDefaultOrEmpty ||
            result.Artifacts.IsDefault)
            return null;
        if (!LocalStateSupportsConstruction(observed.LocalState, policy))
            return null;
        if (observed.SpeciesCounts.IsDefaultOrEmpty ||
            observed.SpeciesCounts.Any(item => item.Count <= 0 ||
                item.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                item.Role is not (GeneratedComponentRoleKind.Lipid or GeneratedComponentRoleKind.Water or
                    GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon)) ||
            observed.SpeciesCounts.Select(item => (item.Role, item.PhysicalSide, item.SpeciesId))
                .Distinct().Count() != observed.SpeciesCounts.Length)
            return null;
        var expectedNames = new Dictionary<GeneratedComponentRoleKind, MolecularRepresentation>
        {
            [GeneratedComponentRoleKind.Lipid] = membrane.SpeciesRepresentations[0],
            [GeneratedComponentRoleKind.Water] = policy.Water,
            [GeneratedComponentRoleKind.PositiveIon] = policy.Sodium,
            [GeneratedComponentRoleKind.NegativeIon] = policy.Chloride
        };
        if (observed.SpeciesCounts.Any(item =>
                item.SpeciesId != expectedNames[item.Role].SpeciesId) ||
            observed.SpeciesCounts.Where(item => item.Role == GeneratedComponentRoleKind.Water)
                .Sum(item => item.Count) != observed.WaterCount ||
            observed.SpeciesCounts.Where(item => item.Role == GeneratedComponentRoleKind.PositiveIon)
                .Sum(item => item.Count) != observed.PositiveIonCount ||
            observed.SpeciesCounts.Where(item => item.Role == GeneratedComponentRoleKind.NegativeIon)
                .Sum(item => item.Count) != observed.NegativeIonCount ||
            derivation.LipidCounts.Any(item =>
                observed.SpeciesCounts.Count(observedItem =>
                    observedItem.Role == GeneratedComponentRoleKind.Lipid &&
                    observedItem.PhysicalSide == item.PhysicalSide &&
                    observedItem.SpeciesId == item.SpeciesId &&
                    observedItem.Count == item.Count) != 1))
            return null;
        var expectedGeneratedAtoms = observed.SpeciesCounts.ToDictionary(
            item => (item.Role, item.PhysicalSide, item.SpeciesId),
            item => (long)item.Count * expectedNames[item.Role].AtomCount);
        if (protein.Molecule.AtomCount + expectedGeneratedAtoms.Values.Sum() != observed.AtomCount)
            return null;

        var topologyCif = result.Artifacts.FirstOrDefault(item => item.Role == "topologyCif");
        var topologyJson = result.Artifacts.FirstOrDefault(item => item.Role == "topologyJson");
        var system = result.Artifacts.FirstOrDefault(item => item.Role == "systemXml");
        var state = result.Artifacts.FirstOrDefault(item => item.Role == "stateXml");
        var mapping = result.Artifacts.FirstOrDefault(item => item.Role == "correspondenceJson");
        if (topologyCif is null || topologyJson is null || system is null || state is null || mapping is null)
            return null;
        foreach (var artifact in new[] { topologyCif, topologyJson, system, state, mapping })
            if (!await ArtifactMatchesAsync(artifact, cancellationToken)) return null;
        SourceToResultCorrespondence? correspondence;
        try
        {
            correspondence = JsonSerializer.Deserialize<SourceToResultCorrespondence>(
                await File.ReadAllTextAsync(mapping.Path, cancellationToken),
                new JsonSerializerOptions(JsonSerializerDefaults.Web));
        }
        catch (JsonException) { return null; }
        if (correspondence is null || !correspondence.Complete ||
            correspondence.SourceId != placement.Proposal.OrientedProtein.CoordinateSha256 ||
            correspondence.ResultId != topologyCif.Sha256 || correspondence.Atoms.Length != observed.AtomCount ||
            correspondence.Atoms.Select(item => item.ResultAtomIndex).Distinct().Count() != observed.AtomCount ||
            correspondence.Atoms.Select(item => item.ResultAtomId)
                .Distinct(StringComparer.Ordinal).Count() != observed.AtomCount ||
            correspondence.Atoms.Any(item => item.ResultAtomIndex < 0 ||
                item.ResultAtomIndex >= observed.AtomCount || !Enum.IsDefined(item.MoleculeRole) ||
                !Enum.IsDefined(item.AtomRole) || !Enum.IsDefined(item.Role)) ||
            correspondence.Atoms.Count(item => item.MoleculeRole is MoleculeRoleKind.Protein or MoleculeRoleKind.RetainedPartner) != protein.Molecule.AtomCount)
            return null;
        var actualAtoms = correspondence.Atoms.OrderBy(item => item.ResultAtomIndex).ToArray();
        var preparedAtoms = protein.Correspondence.Atoms.OrderBy(item => item.ResultAtomIndex).ToArray();
        for (var index = 0; index < preparedAtoms.Length; index++)
        {
            var expected = preparedAtoms[index];
            var actual = actualAtoms[index];
            if (expected.ResultAtomIndex != index || actual.ResultAtomIndex != index ||
                actual.Role != expected.Role || actual.SourceAtomId != expected.SourceAtomId ||
                actual.SourceResidue != expected.SourceResidue ||
                actual.ApprovedChangeId != expected.ApprovedChangeId ||
                actual.MoleculeRole != expected.MoleculeRole || actual.AtomRole != expected.AtomRole ||
                actual.Element != expected.Element ||
                actual.MoleculeRole is not (MoleculeRoleKind.Protein or MoleculeRoleKind.RetainedPartner) ||
                actual.GeneratedSpeciesId is not null || actual.GeneratedComponentRole is not null ||
                actual.PhysicalSide is not null ||
                (actual.Role == AtomOriginKind.Source && (string.IsNullOrWhiteSpace(actual.SourceAtomId) ||
                    actual.SourceResidue is null)) ||
                (actual.Role == AtomOriginKind.Generated && actual.SourceAtomId is not null) ||
                actual.Role is not (AtomOriginKind.Source or AtomOriginKind.Generated))
                return null;
        }
        var generatedCounts = new Dictionary<(GeneratedComponentRoleKind Role, LeafletSide Side, string SpeciesId), long>();
        for (var index = preparedAtoms.Length; index < actualAtoms.Length; index++)
        {
            var atom = actualAtoms[index];
            if (atom.ResultAtomIndex != index || atom.Role != AtomOriginKind.Generated ||
                atom.SourceAtomId is not null || atom.SourceResidue is not null ||
                atom.ApprovedChangeId is not null ||
                atom.GeneratedComponentRole is not (GeneratedComponentRoleKind.Lipid or GeneratedComponentRoleKind.Water or GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon) ||
                atom.AtomRole is not (AtomRoleKind.Head or AtomRoleKind.Body) ||
                atom.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                string.IsNullOrWhiteSpace(atom.GeneratedSpeciesId) ||
                string.IsNullOrWhiteSpace(atom.Element) ||
                string.IsNullOrWhiteSpace(atom.ResultAtomId) ||
                !atom.ResultAtomId.StartsWith("component:", StringComparison.Ordinal) ||
                atom.MoleculeRole != (atom.GeneratedComponentRole is GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon
                    ? MoleculeRoleKind.Ion : atom.GeneratedComponentRole == GeneratedComponentRoleKind.Lipid
                        ? MoleculeRoleKind.Lipid : MoleculeRoleKind.Water))
                return null;
            var key = (atom.GeneratedComponentRole!.Value, atom.PhysicalSide!.Value, atom.GeneratedSpeciesId);
            generatedCounts[key] = generatedCounts.GetValueOrDefault(key) + 1;
        }
        if (generatedCounts.Count != expectedGeneratedAtoms.Count ||
            generatedCounts.Any(item => !expectedGeneratedAtoms.TryGetValue(item.Key, out var count) ||
                item.Value != count))
            return null;

        var molecule = new MolecularArtifact(topologyCif.Sha256, topologyCif.Path, topologyCif.Sha256,
            topologyJson.Path, system.Path, state.Path, observed.AtomCount,
            string.Join(" × ", derivation.CellAngstrom.Select(value => value.ToString("G17", System.Globalization.CultureInfo.InvariantCulture))),
            topologyJson.Sha256, mapping.Path, mapping.Sha256, system.Sha256, state.Sha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"), molecule.Id,
            result.Provider?.Name ?? "local scientific worker", "Whole-system construction and parameter assessment",
            $"{observed.AtomCount} mapped atoms; {observed.WaterCount} water; {observed.PositiveIonCount} sodium and {observed.NegativeIonCount} chloride ions",
            $"Attempt {attempt.Id}; policy {policy.Id}; placement {placement.Id}; " +
            $"OpenMM build {policy.Construction.ProviderVersion}; native patch SHA-256 {observed.NativePatchSha256}",
            "Same-invocation provider cell, populations, protein preservation, combined parameters and initial contacts were checked. " +
            "This is a candidate for researcher review before required minimization, not a completed preparation.",
            EvidenceBearing.Supports));
        var conditions = $"Fixed nominal pH {revision.Conditions.NominalPh:G6}; intended NaCl {derivation.IntendedNaClMolar:G6} M; " +
                         $"finite-cell estimated NaCl {derivation.EstimatedNaClMolar:G6} M; " +
                         "thermal equilibration has not been established.";
        return new ConstructedExplicitSystem(molecule.Id, attempt, molecule, derivation,
            correspondence, derivation.LipidCounts, evidence, ImmutableArray<ScientificFinding>.Empty, conditions,
            observed.LocalState, observed.ActualCellAngstrom);
    }

    private static async Task<bool> ArtifactMatchesAsync(WorkerArtifact artifact,
        CancellationToken cancellationToken)
    {
        if (string.IsNullOrWhiteSpace(artifact.Path) || artifact.Sha256 is not { Length: 64 } ||
            !artifact.Sha256.All(Uri.IsHexDigit) || !File.Exists(artifact.Path))
            return false;
        try
        {
            await using var stream = File.OpenRead(artifact.Path);
            var hash = await SHA256.HashDataAsync(stream, cancellationToken);
            return Convert.ToHexString(hash).Equals(artifact.Sha256,
                StringComparison.OrdinalIgnoreCase);
        }
        catch (IOException) { return false; }
        catch (UnauthorizedAccessException) { return false; }
    }

    private static bool ValidLocalStatePolicy(ApplicablePreparationPolicy policy)
    {
        var spec = policy.LocalStateObservation;
        if (spec is null || spec.ContactRolePairs.IsDefaultOrEmpty ||
            spec.RequiredMetricNames.IsDefaultOrEmpty || spec.AtomRadiusByElementAngstrom.IsEmpty ||
            spec.ContactRolePairs.Any(pair => !Enum.IsDefined(pair.FirstMoleculeRole) ||
                !Enum.IsDefined(pair.SecondMoleculeRole)) ||
            spec.ContactRolePairs.Distinct().Count() != spec.ContactRolePairs.Length ||
            spec.RequiredMetricNames.Distinct(StringComparer.Ordinal).Count() != spec.RequiredMetricNames.Length ||
            spec.RequiredMetricNames.Any(name => name is not (
                "minimumIntermolecularDistanceAngstrom" or
                "minimumIntermolecularHeavyAtomDistanceAngstrom" or "upperLipidHeadMeanZAngstrom" or
                "lowerLipidHeadMeanZAngstrom" or "leafletHeadSeparationAngstrom" or
                "proteinBilayerMidplaneOffsetAngstrom")) ||
            !spec.RequiredMetricNames.Contains("leafletHeadSeparationAngstrom") ||
            !spec.RequiredMetricNames.Contains("proteinBilayerMidplaneOffsetAngstrom") ||
            spec.AtomRadiusByElementAngstrom.Any(item => string.IsNullOrWhiteSpace(item.Key) ||
                !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(spec.ContactSearchRadiusAngstrom) ||
            spec.ContactSearchRadiusAngstrom <= 0 ||
            spec.MaximumReportedPairs <= 0 || !spec.UsePeriodicBoundary ||
            policy.ConstructionCriteria.IsDefaultOrEmpty ||
            policy.ContactCriteria.IsDefaultOrEmpty ||
            policy.ConstructionCriteria.Select(item => item.MeasurementName)
                .Distinct(StringComparer.Ordinal).Count() != policy.ConstructionCriteria.Length ||
            policy.ConstructionCriteria.Any(item => item.MeasurementName is
                "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom") ||
            !policy.ConstructionCriteria.Any(item =>
                item.MeasurementName == "minimumIntermolecularHeavyAtomDistanceAngstrom" &&
                item.Minimum is double minimum && double.IsFinite(minimum) && minimum > 0) ||
            !HasOrganizationCriterion(policy.ConstructionCriteria, "leafletHeadSeparationAngstrom",
                "bilayer", positiveMinimum: true) ||
            !HasOrganizationCriterion(policy.ConstructionCriteria, "proteinBilayerMidplaneOffsetAngstrom",
                "proteinVsBilayer", positiveMinimum: false) ||
            policy.ContactCriteria.Select(item => (item.StageKind, item.FirstMoleculeRole,
                item.SecondMoleculeRole)).Distinct().Count() != policy.ContactCriteria.Length ||
            policy.ContactCriteria.Any(item => item.MinimumPairsWithinSearchRadius < 0 ||
                !spec.ContactRolePairs.Any(pair => pair.FirstMoleculeRole == item.FirstMoleculeRole &&
                    pair.SecondMoleculeRole == item.SecondMoleculeRole) ||
                item.MinimumNearestDistanceAngstrom is double lower &&
                    (!double.IsFinite(lower) || lower <= 0) ||
                item.MaximumNearestDistanceAngstrom is double upper &&
                    (!double.IsFinite(upper) || upper <= 0 || upper > spec.ContactSearchRadiusAngstrom) ||
                item.MinimumNearestDistanceAngstrom is double minimum &&
                    item.MaximumNearestDistanceAngstrom is double maximum && minimum > maximum) ||
            !new StageKind?[] { null, StageKind.Minimization, StageKind.Equilibration }
                .Where(kind => kind is null || kind != StageKind.Equilibration ||
                    policy.OptionalEquilibration is not null)
                .All(kind => policy.ContactCriteria.Any(item => item.StageKind == kind &&
                    item.FirstMoleculeRole == MoleculeRoleKind.Protein && item.SecondMoleculeRole == MoleculeRoleKind.Lipid &&
                    item.MinimumPairsWithinSearchRadius > 0 &&
                    item.MaximumNearestDistanceAngstrom is double maximum &&
                    double.IsFinite(maximum) && maximum > 0)))
            return false;
        return policy.ConstructionCriteria.All(criterion =>
            spec.RequiredMetricNames.Contains(criterion.MeasurementName) &&
            !string.IsNullOrWhiteSpace(criterion.Unit) && !string.IsNullOrWhiteSpace(criterion.Scope) &&
            (criterion.Minimum is null || double.IsFinite(criterion.Minimum.Value)) &&
            (criterion.Maximum is null || double.IsFinite(criterion.Maximum.Value)));
    }

    private static bool HasOrganizationCriterion(ImmutableArray<LocalStateCriterion> criteria,
        string name, string scope, bool positiveMinimum)
    {
        var matching = criteria.Where(item => item.MeasurementName == name).ToArray();
        return matching.Length == 1 && matching[0].Unit == "angstrom" && matching[0].Scope == scope &&
               matching[0].Minimum is double minimum && double.IsFinite(minimum) &&
               (!positiveMinimum || minimum > 0) &&
               matching[0].Maximum is double maximum && double.IsFinite(maximum) &&
               maximum >= minimum;
    }

    private static bool ValidStageProteinGeometryPolicy(ApplicablePreparationPolicy policy)
    {
        var spec = policy.StageProteinGeometryMeasurement;
        var required = new[] { "covalentBond", "chainContinuity", "nonbondedDistance" };
        if (spec is null || spec.RequiredKinds.IsDefaultOrEmpty ||
            spec.RequiredKinds.Distinct(StringComparer.Ordinal).Count() != spec.RequiredKinds.Length ||
            required.Any(kind => !spec.RequiredKinds.Contains(kind)) ||
            spec.RequiredKinds.Any(kind => !required.Contains(kind)) ||
            spec.AtomRadiusByElementAngstrom.IsEmpty ||
            spec.AtomRadiusByElementAngstrom.Any(item => string.IsNullOrWhiteSpace(item.Key) ||
                !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(spec.NeighborSearchRadiusAngstrom) ||
            spec.NeighborSearchRadiusAngstrom <= 0 || spec.ExcludedBondHops < 0 ||
            spec.MaximumReportedPairs <= 0 || policy.StageProteinGeometryCriteria.IsDefaultOrEmpty)
            return false;
        return new[] { StageKind.Minimization, StageKind.Equilibration }
            .Where(kind => kind == StageKind.Minimization || policy.OptionalEquilibration is not null)
            .All(kind => spec.RequiredKinds.All(requiredKind =>
            {
                var criteria = policy.StageProteinGeometryCriteria.Where(item =>
                    item.StageKind == kind && item.Criterion?.Kind == requiredKind).ToArray();
                if (criteria.Length != 1)
                    return false;
                var criterion = criteria[0].Criterion;
                return (criterion.MinimumObservedAngstrom is not null ||
                        criterion.MaximumObservedAngstrom is not null) &&
                       (requiredKind != "covalentBond" || !criterion.AllowNotApplicable) &&
                       (criterion.MinimumObservedAngstrom is not double minimum || double.IsFinite(minimum)) &&
                       (criterion.MaximumObservedAngstrom is not double maximum || double.IsFinite(maximum)) &&
                       (criterion.MinimumObservedAngstrom is not double lower ||
                        criterion.MaximumObservedAngstrom is not double upper || lower <= upper);
            }));
    }

    private static bool LocalStateSupportsConstruction(LocalStateObservations? observed,
        ApplicablePreparationPolicy policy)
    {
        if (observed is null || observed.Standing != ObservationStanding.Observed ||
            observed.Measurements.IsDefault || observed.LocatedContacts.IsDefault ||
            observed.CoveredRolePairs.IsDefault || observed.RolePairMeasurements.IsDefault ||
            observed.CoveredRolePairs.Distinct().Count() !=
                observed.CoveredRolePairs.Length ||
            observed.RolePairMeasurements.Any(item => !Enum.IsDefined(item.FirstMoleculeRole) ||
                !Enum.IsDefined(item.SecondMoleculeRole)) ||
            observed.RolePairMeasurements.Select(item => (item.FirstMoleculeRole,
                item.SecondMoleculeRole)).Distinct().Count() != observed.RolePairMeasurements.Length ||
            policy.LocalStateObservation.ContactRolePairs.Any(pair =>
                !observed.CoveredRolePairs.Contains(pair) ||
                !observed.RolePairMeasurements.Any(item =>
                    item.FirstMoleculeRole == pair.FirstMoleculeRole &&
                    item.SecondMoleculeRole == pair.SecondMoleculeRole &&
                    item.PairsWithinSearchRadius >= 0 &&
                    (item.PairsWithinSearchRadius == 0 && item.MinimumDistanceAngstrom is null ||
                     item.PairsWithinSearchRadius > 0 && item.MinimumDistanceAngstrom is double distance &&
                     double.IsFinite(distance) && distance > 0))) ||
            observed.LocatedContacts.Any(contact =>
                contact.FirstAtomIndex < 0 || contact.SecondAtomIndex < 0 ||
                !Enum.IsDefined(contact.FirstMoleculeRole) || !Enum.IsDefined(contact.SecondMoleculeRole) ||
                !double.IsFinite(contact.DistanceAngstrom) || contact.DistanceAngstrom < 0 ||
                !double.IsFinite(contact.RadiusSumAngstrom) || contact.RadiusSumAngstrom <= 0))
            return false;
        foreach (var name in policy.LocalStateObservation.RequiredMetricNames)
        {
            var measurements = observed.Measurements.Where(item => item.Name == name).ToArray();
            if (measurements.Length != 1 || !double.IsFinite(measurements[0].Value))
                return false;
            var criterion = policy.ConstructionCriteria.FirstOrDefault(item => item.MeasurementName == name);
            if (criterion is not null &&
                (measurements[0].Unit != criterion.Unit || measurements[0].Scope != criterion.Scope ||
                 criterion.Minimum is double minimum && measurements[0].Value < minimum ||
                 criterion.Maximum is double maximum && measurements[0].Value > maximum))
                return false;
        }
        foreach (var criterion in policy.ContactCriteria.Where(item => item.StageKind is null))
        {
            var pair = observed.RolePairMeasurements.Single(item =>
                item.FirstMoleculeRole == criterion.FirstMoleculeRole &&
                item.SecondMoleculeRole == criterion.SecondMoleculeRole);
            if (pair.PairsWithinSearchRadius < criterion.MinimumPairsWithinSearchRadius ||
                criterion.MinimumNearestDistanceAngstrom is double minimum &&
                    (pair.MinimumDistanceAngstrom is not double distance || distance < minimum) ||
                criterion.MaximumNearestDistanceAngstrom is double maximum &&
                    (pair.MinimumDistanceAngstrom is not double distance2 || distance2 > maximum))
                return false;
        }
        return true;
    }
}
