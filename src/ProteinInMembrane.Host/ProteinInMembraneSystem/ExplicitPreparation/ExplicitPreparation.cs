using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation;

/// <summary>Coordinates one identified attempt from assessed inputs to an eligible explicit system.</summary>
public sealed class ExplicitPreparation
{
    private const double MoleculesPerMolarAngstromCubed = 6.02214076e-4;
    // The pinned Memgen 2026.3.25 source uses this value for integer salt rounding.
    private const double MemgenMoleculesPerMolarAngstromCubed = 6.02214086e-4;
    private const long MaximumDiagnosticFileBytes = 512L * 1024 * 1024;
    private const long MaximumDiagnosticTotalBytes = 2L * 1024 * 1024 * 1024;
    // A trial account may expose method files, never arbitrary worker output.
    private static readonly IReadOnlyDictionary<string, PreparationPhase> MemgenDiagnosticRoles =
        new Dictionary<string, PreparationPhase>(StringComparer.Ordinal)
        {
            ["providerOptions"] = PreparationPhase.ProviderPopulation,
            ["providerLog"] = PreparationPhase.ProviderPopulation,
            ["providerCommand"] = PreparationPhase.ProviderPopulation,
            ["providerStdout"] = PreparationPhase.ProviderPopulation,
            ["preparedChargeLeapInput"] = PreparationPhase.ProviderPopulation,
            ["preparedChargeLeapLog"] = PreparationPhase.ProviderPopulation,
            ["preparedChargeLeapTopology"] = PreparationPhase.ProviderPopulation,
            ["preparedChargeLeapCoordinates"] = PreparationPhase.ProviderPopulation,
            ["preparedChargePdb"] = PreparationPhase.ProviderPopulation,
            ["providerInput"] = PreparationPhase.ProviderPopulation,
            ["providerProtein"] = PreparationPhase.ProviderPopulation,
            ["providerPackmolInput"] = PreparationPhase.ProviderPacking,
            ["providerPackmolLog"] = PreparationPhase.ProviderPacking,
            ["providerPacked"] = PreparationPhase.ProviderPacking,
            ["providerPostCleanup"] = PreparationPhase.ProviderCleanup,
            ["providerLeapInput"] = PreparationPhase.AmberParameterization,
            ["providerLeapLog"] = PreparationPhase.AmberParameterization,
            ["providerParameterized"] = PreparationPhase.AmberParameterization,
            ["amberCoordinates"] = PreparationPhase.AmberParameterization,
            ["amberTopology"] = PreparationPhase.AmberParameterization,
            ["providerRestrainedInput"] = PreparationPhase.ProviderConditioningRestrained,
            ["providerRestrainedLog"] = PreparationPhase.ProviderConditioningRestrained,
            ["providerRestrainedRestart"] = PreparationPhase.ProviderConditioningRestrained,
            ["providerUnrestrainedInput"] = PreparationPhase.ProviderConditioningUnrestrained,
            ["providerUnrestrainedLog"] = PreparationPhase.ProviderConditioningUnrestrained,
            ["amberFinalRestart"] = PreparationPhase.ProviderConditioningUnrestrained
        };
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
        // A class-level recipe still binds this exact checked subject and pair.
        // A non-null scope is retained for historical, explicitly exact policies.
        if (scope is null)
            return protein.StudyRevisionId == revision.Id && membrane.StudyRevisionId == revision.Id &&
                revision.IntendedProtein?.Id == protein.Intended.Id &&
                revision.Membrane?.Id == membrane.Intended.Id &&
                proposal.PreparedProteinId == protein.Id &&
                proposal.MembraneModelId == membrane.Intended.Id &&
                !protein.Intended.Chains.IsDefaultOrEmpty &&
                protein.Correspondence.Complete &&
                protein.Correspondence.ResultId == protein.Molecule.CoordinateSha256;
        if (!HashLike(scope.SourceCoordinateSha256) ||
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
            policy.MaximumMinimizationIterations != 20000 ||
            policy.FinalUnrestrainedRmsForceTargetKjMolNm != 10.0 ||
            policy.SystemSettings is not
                { NonbondedMethod: "PME", NonbondedCutoffNanometers: 1.0,
                  Constraints: "HBonds", RigidWater: true,
                  EwaldErrorTolerance: 0.0005, SwitchDistanceNanometers: null,
                  UseDispersionCorrection: true, RemoveCMMotion: true,
                  HydrogenMassDaltons: null } ||
            construction.CoveredTopologyKinds.IsDefaultOrEmpty ||
            !construction.CoveredTopologyKinds.Contains(placement.Proposal.TopologyKind) ||
            membrane.SpeciesRepresentations.IsDefaultOrEmpty ||
            revision.Conditions.TargetNaClMolar != 0.15 ||
            placement.Proposal.MidplaneAngstrom is not double midplane ||
            !double.IsFinite(midplane) ||
            !HashMatches(protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256) ||
            !HashMatches(protein.Molecule.TopologyPath, protein.Molecule.TopologySha256) ||
            !HashMatches(placement.Proposal.OrientedProtein.CoordinatePath,
                placement.Proposal.OrientedProtein.CoordinateSha256) ||
            !HashMatches(placement.Proposal.OrientedProtein.TopologyPath,
                placement.Proposal.OrientedProtein.TopologySha256))
            return NoAttempt("Corresponding supported inputs and an identified complete construction policy are required.");

        if (membrane.SpeciesRepresentations.Any(item =>
                !HashMatches(item.TemplatePath, item.TemplateSha256) ||
                !HashMatches(item.CoordinateTemplatePath, item.CoordinateTemplateSha256)) ||
            new[] { policy.Water, policy.Sodium, policy.Chloride }.Any(item =>
                !HashMatches(item.TemplatePath, item.TemplateSha256) ||
                !HashMatches(item.CoordinateTemplatePath, item.CoordinateTemplateSha256)))
            return Refusal("A selected molecular representation changed before preparation admission.");

        if (construction.Route == ConstructionRouteKind.PackmolMemgen)
            return await StartMemgenAsync(attempt, revision, protein, membrane, placement, policy,
                workingDirectory, midplane, onAccepted, progress, cancellationToken);

        if (membrane.SpeciesRepresentations.Length != 1 ||
            membrane.SpeciesRepresentations[0].SpeciesId != construction.LipidTypeArgument ||
            !PureSelectedLeaflet(membrane.Intended.Upper, construction.LipidTypeArgument!) ||
            !PureSelectedLeaflet(membrane.Intended.Lower, construction.LipidTypeArgument!))
            return NoAttempt("The selected native recipe requires its exact qualified pure lipid on both leaflets.");

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
            (construction.NativePatchMode == "mapped-lipid21-zenodo-popc" &&
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
                construction.NativeSourcePatchSha256, construction.RemovedNativeLipidResidueIds,
                ProviderAssets: ImmutableArray<ProviderAsset>.Empty,
                SelectedSpeciesRepresentations: ImmutableArray<MolecularRepresentation>.Empty));
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
                State(attempt.Id, StageExecutionStanding.Running,
                    "The checked native construction is continuing to required minimization.",
                    PreparationPhase.HandoffChecks));
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

    private async Task<PreparationStartResult> StartMemgenAsync(
        PreparationAttempt attempt, StudyRevision revision, AssessedPreparedProtein protein,
        AssessedMembraneModel membrane, AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy, string workingDirectory, double midplane,
        Action<PreparationAttempt>? onAccepted, IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
    {
        var construction = policy.Construction;
        var settings = construction.Memgen!;
        var constructionBoundMessage = $"The complete construction trials exceeded the declared {construction.MaximumConstructionSeconds}-second attempt bound.";
        var targets = new[] { membrane.Intended.Upper, membrane.Intended.Lower };
        var selected = targets.SelectMany(side => side.Fractions)
            .Where(fraction => fraction.Fraction > 0).Select(fraction => fraction.SpeciesId)
            .ToHashSet(StringComparer.Ordinal);
        var representations = membrane.SpeciesRepresentations;
        if (targets[0].PhysicalSide != LeafletSide.Upper ||
            targets[1].PhysicalSide != LeafletSide.Lower ||
            targets.Any(side => side.Fractions.IsDefaultOrEmpty ||
                side.Fractions.Any(fraction => fraction.Fraction < 0 ||
                    !double.IsFinite(fraction.Fraction) || string.IsNullOrWhiteSpace(fraction.SpeciesId)) ||
                Math.Abs(side.Fractions.Sum(fraction => fraction.Fraction) - 1) > 1e-9 ||
                side.Fractions.Select(fraction => fraction.SpeciesId)
                    .Distinct(StringComparer.Ordinal).Count() != side.Fractions.Length) ||
            selected.Count == 0 || !selected.SetEquals(representations.Select(item => item.SpeciesId)) ||
            representations.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
                representations.Length ||
            representations.Any(item => !construction.CoveredSpeciesIds.Contains(item.SpeciesId) ||
                item.AtomCount <= 0 || Math.Abs(item.NetChargeElementary) > 1e-8 ||
                string.IsNullOrWhiteSpace(item.ChemistryId) ||
                item.ForceFieldFamily != "Lipid21") ||
            policy.Water.SpeciesId != "HOH" || policy.Water.Category != "water" ||
            policy.Sodium.SpeciesId != "NA" || policy.Sodium.Category != "ion" ||
            policy.Chloride.SpeciesId != "CL" || policy.Chloride.Category != "ion" ||
            Math.Abs(policy.Water.NetChargeElementary) > 1e-8 ||
            Math.Abs(policy.Sodium.NetChargeElementary - 1) > 1e-8 ||
            Math.Abs(policy.Chloride.NetChargeElementary + 1) > 1e-8)
            return NoAttempt("The selected exact composition or retained chemistry has no corresponding Memgen route.");

        if (!HashMatches(protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256) ||
            !protein.Correspondence.Complete ||
            protein.Correspondence.ResultId != protein.Molecule.CoordinateSha256 ||
            protein.Correspondence.Atoms.IsDefault ||
            protein.Correspondence.Atoms.Length != protein.Molecule.AtomCount ||
            !HashMatches(placement.Proposal.OrientedProtein.CoordinatePath,
                placement.Proposal.OrientedProtein.CoordinateSha256) ||
            !HashMatches(placement.Proposal.OrientedProtein.TopologyPath,
                placement.Proposal.OrientedProtein.TopologySha256) ||
            policy.ForceFieldFiles.Any(asset => string.IsNullOrWhiteSpace(asset.Id) ||
                string.IsNullOrWhiteSpace(asset.Version) ||
                !HashMatches(asset.Path, asset.Sha256)) ||
            construction.ProviderAssets.Any(asset => !HashMatches(asset.Path, asset.Sha256)))
            return Refusal("An exact prepared construct, positioned construct, provider asset or parameter asset is unavailable before start.");

        onAccepted?.Invoke(attempt);
        using var bounded = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
        bounded.CancelAfter(TimeSpan.FromSeconds(construction.MaximumConstructionSeconds));
        var trials = ImmutableArray.CreateBuilder<ConstructionTrialSummary>();
        var lateral = settings.LateralPaddingAngstrom;
        var aqueous = settings.AqueousPaddingAngstrom;
        var forceFieldFiles = policy.ForceFieldFiles.DistinctBy(asset => asset.Sha256,
            StringComparer.OrdinalIgnoreCase).ToImmutableArray();
        PreparationStartResult EndWithoutObservedTrial(string trialId, int trialIndex,
            ConstructionTrialStanding trialStanding, WorkerResultStanding resultStanding,
            string code, string message, PreparationPhase phase = PreparationPhase.ProviderPopulation)
        {
            var incomplete = new ConstructionTrialSummary(trialId, trialIndex, trialStanding,
                lateral, aqueous, ImmutableArray<SpeciesCount>.Empty,
                ImmutableArray<SpeciesCount>.Empty, ImmutableArray<SpeciesCount>.Empty,
                ImmutableArray<double>.Empty, ImmutableArray<double>.Empty, null, code, message,
                ImmutableArray<TrialDiagnosticArtifact>.Empty);
            trials.Add(incomplete);
            progress?.Report(State(attempt.Id, StageExecutionStanding.Running, message,
                phase, trialId, trialIndex, code, incomplete));
            return Failed(attempt, resultStanding, message, code, trials.ToImmutable(),
                trialId, trialIndex, phase);
        }
        PreparationStartResult InterruptedHandoff(ConstructionTrialSummary trial)
        {
            var stopped = cancellationToken.IsCancellationRequested;
            var code = stopped ? "stopped" : "resourceLimit";
            var message = stopped
                ? "Construction handoff checks were stopped before completion."
                : $"Construction handoff checks exceeded the declared {construction.MaximumConstructionSeconds}-second attempt bound.";
            var incomplete = trial with
            {
                Standing = stopped ? ConstructionTrialStanding.Stopped : ConstructionTrialStanding.Failed,
                FailureCode = code,
                Message = message,
                DiagnosticArtifacts = ImmutableArray<TrialDiagnosticArtifact>.Empty
            };
            trials.Add(incomplete);
            progress?.Report(State(attempt.Id, StageExecutionStanding.Running, message,
                PreparationPhase.HandoffChecks, trial.TrialId, trial.TrialIndex, code, incomplete));
            return Failed(attempt, stopped ? WorkerResultStanding.Stopped : WorkerResultStanding.Failed,
                message, code, trials.ToImmutable(), trial.TrialId, trial.TrialIndex,
                PreparationPhase.HandoffChecks);
        }
        for (var index = 0; index <= settings.MaximumGeometryRetries; index++)
        {
            if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
            {
                var stopped = cancellationToken.IsCancellationRequested;
                return Failed(attempt, stopped ? WorkerResultStanding.Stopped : WorkerResultStanding.Failed,
                    stopped ? "Construction was stopped before another provider trial began." :
                        constructionBoundMessage,
                    stopped ? "stopped" : "resourceLimit", trials.ToImmutable(),
                    phase: PreparationPhase.ProviderPopulation);
            }
            // Pinned Memgen clamps both z faces beyond the 23 Å leaflet envelope
            // by --dist_wat. This lower bound alone can make a retry impossible.
            var minimumZAngstrom = 2 * (settings.LeafletEnvelopeAngstrom + aqueous);
            if (minimumZAngstrom > construction.MaximumCellDimensionAngstrom)
                return Failed(attempt, WorkerResultStanding.Failed,
                    $"The next provider trial requires at least {minimumZAngstrom:G} Å in z, " +
                    $"above the declared {construction.MaximumCellDimensionAngstrom:G} Å cell bound.",
                    "resourceRefused", trials.ToImmutable(), phase: PreparationPhase.ProviderPopulation);
            var trialId = Guid.NewGuid().ToString("N");
            var trialDirectory = Path.Combine(workingDirectory, "trials", trialId);
            Directory.CreateDirectory(trialDirectory);
            progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                $"Memgen complete construction trial {index + 1} is running.",
                PreparationPhase.ProviderPopulation, trialId, index));
            var payload = new ConstructionPayload(revision.Id, attempt.Id,
                placement.Proposal.OrientedProtein.CoordinatePath,
                placement.Proposal.OrientedProtein.CoordinateSha256,
                protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256,
                protein.Correspondence,
                placement.Proposal.OrientedProtein.TopologyPath!,
                placement.Proposal.OrientedProtein.TopologySha256!,
                null, policy.Water, policy.Sodium, policy.Chloride,
                null, null, construction.ProviderName, construction.ProviderVersion,
                null, construction.PositiveIonArgument, construction.NegativeIonArgument,
                midplane / 10.0, lateral / 10.0, revision.Conditions.TargetNaClMolar,
                forceFieldFiles, policy.SystemSettings, policy.LocalStateObservation,
                construction.MaximumAtomCount, construction.MaximumCellDimensionAngstrom,
                Route: ConstructionRouteKind.PackmolMemgen,
                SaltConvention: SaltConventionKind.MemgenChargeCompensated,
                ProviderAssets: construction.ProviderAssets, Memgen: settings,
                TrialId: trialId, TrialIndex: index,
                TargetUpper: membrane.Intended.Upper, TargetLower: membrane.Intended.Lower,
                SelectedSpeciesRepresentations: representations,
                LateralPaddingAngstrom: lateral, AqueousPaddingAngstrom: aqueous);
            var request = new ScientificWorkRequest<ConstructionPayload>(
                Guid.NewGuid().ToString("N"), trialDirectory, payload);
            WorkerResult<ConstructionObservations> built;
            try
            {
                built = await _worker.ConstructSystemAsync(request, bounded.Token);
            }
            catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
            {
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Stopped, WorkerResultStanding.Stopped, "stopped",
                    "Construction was stopped before a checked handoff was established.",
                    PreparationPhase.ProviderPopulation);
            }
            catch (OperationCanceledException) when (bounded.IsCancellationRequested)
            {
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Failed, WorkerResultStanding.Failed, "resourceLimit",
                    constructionBoundMessage,
                    PreparationPhase.ProviderPopulation);
            }
            catch (OperationCanceledException)
            {
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Unobserved, WorkerResultStanding.Unobserved,
                    "uncorrelatedWorkerCancellation",
                    "The provider crossing was cancelled without an attributable stop or resource limit.");
            }
            if (cancellationToken.IsCancellationRequested)
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Stopped, WorkerResultStanding.Stopped, "stopped",
                    "Construction was stopped before a checked handoff was established.");
            if (bounded.IsCancellationRequested)
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Failed, WorkerResultStanding.Failed, "resourceLimit",
                    constructionBoundMessage);
            if (built.RequestId != request.RequestId || built.StudyRevisionId != revision.Id ||
                built.AttemptId != attempt.Id)
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Unobserved, WorkerResultStanding.Unobserved,
                    "uncorrelatedProviderResult",
                    "The provider result did not identify this exact attempt and request.",
                    PreparationPhase.ProviderPopulation);

            if (built.Standing != WorkerResultStanding.Observed || built.Observations is null)
            {
                var failedTrial = FailedTrialDetail(built.FailureDetails, trialId, index,
                    lateral, aqueous, built.FailureCode);
                if (failedTrial is not null)
                {
                    failedTrial = failedTrial with
                    {
                        DiagnosticArtifacts = built.Standing == WorkerResultStanding.Failed
                            ? await VerifiedMemgenDiagnosticsAsync(built.Artifacts, trialDirectory, bounded.Token)
                            : ImmutableArray<TrialDiagnosticArtifact>.Empty
                    };
                    trials.Add(failedTrial);
                    progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                        failedTrial.Message ?? "A complete provider trial failed its reported check.",
                        PreparationPhase.ProviderPopulation, trialId, index,
                        failedTrial.FailureCode, failedTrial));
                }
                var code = built.FailureCode;
                var axis = RetryAxis(code, built.FailureDetails);
                if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                {
                    var stopped = cancellationToken.IsCancellationRequested;
                    var terminalCode = stopped ? "stopped" : "resourceLimit";
                    var terminalMessage = stopped
                        ? "Construction was stopped before a checked handoff was established."
                        : constructionBoundMessage;
                    if (failedTrial is not null)
                    {
                        var terminalTrial = failedTrial with
                        {
                            Standing = stopped ? ConstructionTrialStanding.Stopped : ConstructionTrialStanding.Failed,
                            FailureCode = terminalCode,
                            Message = terminalMessage,
                            DiagnosticArtifacts = ImmutableArray<TrialDiagnosticArtifact>.Empty
                        };
                        trials[^1] = terminalTrial;
                        progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                            terminalMessage, PreparationPhase.ProviderPopulation, trialId, index,
                            terminalCode, terminalTrial));
                    }
                    return failedTrial is null
                        ? EndWithoutObservedTrial(trialId, index,
                            stopped ? ConstructionTrialStanding.Stopped : ConstructionTrialStanding.Failed,
                            stopped ? WorkerResultStanding.Stopped : WorkerResultStanding.Failed,
                            terminalCode, terminalMessage)
                        : Failed(attempt, stopped ? WorkerResultStanding.Stopped : WorkerResultStanding.Failed,
                            terminalMessage, terminalCode, trials.ToImmutable(), trialId, index);
                }
                if (failedTrial is null)
                {
                    var standing = built.Standing switch
                    {
                        WorkerResultStanding.Failed => ConstructionTrialStanding.Failed,
                        WorkerResultStanding.Stopped => ConstructionTrialStanding.Stopped,
                        _ => ConstructionTrialStanding.Unobserved
                    };
                    return EndWithoutObservedTrial(trialId, index, standing,
                        built.Standing == WorkerResultStanding.Observed
                            ? WorkerResultStanding.Unobserved : built.Standing,
                        code ?? "missingProviderTrialAccount",
                        built.FailureMessage ?? "The provider returned no attributable trial account.");
                }
                if (built.Standing == WorkerResultStanding.Failed && failedTrial is not null &&
                    axis is not null && index < settings.MaximumGeometryRetries)
                {
                    if (axis == "lateral")
                        lateral = settings.LateralPaddingAngstrom *
                            (lateral == settings.LateralPaddingAngstrom ? 2 : 4);
                    else
                        aqueous = settings.AqueousPaddingAngstrom *
                            (aqueous == settings.AqueousPaddingAngstrom ? 2 : 4);
                    continue;
                }
                return Failed(attempt, built.Standing,
                    built.FailureMessage ?? "The full provider did not establish a checked construction.",
                    code, trials.ToImmutable(), trialId, index);
            }

            var observed = built.Observations;
            var trial = observed.Trial;
            if (trial is null || trial.TrialId != trialId || trial.TrialIndex != index ||
                trial.Standing != ConstructionTrialStanding.Checked ||
                trial.LateralPaddingAngstrom != lateral ||
                trial.AqueousPaddingAngstrom != aqueous ||
                observed.Conditions is null || trial.Conditions is null ||
                !SameConditionAccount(trial.Conditions, observed.Conditions))
                return EndWithoutObservedTrial(trialId, index,
                    ConstructionTrialStanding.Unobserved, WorkerResultStanding.Unobserved,
                    "incompleteTrialAccount",
                    "The successful provider trial lacks its exact population and condition account.",
                    PreparationPhase.HandoffChecks);
            progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                "Provider conditioning completed; checking the exact Amber handoff.",
                PreparationPhase.HandoffChecks, trialId, index));
            var derived = TryDeriveMemgen(attempt, revision, protein, membrane, policy,
                    observed, trials.ToImmutable().Add(trial), out var derivation, out var cause,
                    out var clearanceAxis);
            if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                return InterruptedHandoff(trial);
            if (!derived)
            {
                var rejected = trial with { Standing = ConstructionTrialStanding.Failed,
                    FailureCode = cause,
                    Message = "Product handoff checks rejected the provider observation.",
                    DiagnosticArtifacts = await VerifiedMemgenDiagnosticsAsync(
                        observed.ProviderIntermediates, trialDirectory, bounded.Token) };
                if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                    return InterruptedHandoff(trial);
                trials.Add(rejected);
                progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                    rejected.Message, PreparationPhase.HandoffChecks, trialId, index,
                    cause, rejected));
                if (cause == "insufficientCellClearance" && index < settings.MaximumGeometryRetries)
                {
                    if (clearanceAxis == "lateral")
                        lateral = settings.LateralPaddingAngstrom *
                            (lateral == settings.LateralPaddingAngstrom ? 2 : 4);
                    else if (clearanceAxis == "aqueous")
                        aqueous = settings.AqueousPaddingAngstrom *
                            (aqueous == settings.AqueousPaddingAngstrom ? 2 : 4);
                    else return Failed(attempt, WorkerResultStanding.Failed,
                        "The insufficient-clearance observation did not identify a retry axis.",
                        cause, trials.ToImmutable(), trialId, index);
                    continue;
                }
                return Failed(attempt, WorkerResultStanding.Failed,
                    "The actual full-provider counts, chemistry or cell did not establish a faithful handoff.",
                    cause, trials.ToImmutable(), trialId, index);
            }
            ConstructedExplicitSystem? constructed;
            try
            {
                constructed = await AssessConstructedMemgenAsync(attempt, revision, protein,
                    membrane, placement, policy, derivation!, built, bounded.Token);
            }
            catch (OperationCanceledException) when (bounded.IsCancellationRequested)
            {
                return InterruptedHandoff(trial);
            }
            if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                return InterruptedHandoff(trial);
            if (constructed is null)
            {
                var rejected = trial with { Standing = ConstructionTrialStanding.Failed,
                    FailureCode = "invalidHandoff",
                    Message = "Amber artifacts, retained correspondence, parameters or contacts did not pass product handoff checks.",
                    DiagnosticArtifacts = await VerifiedMemgenDiagnosticsAsync(
                        observed.ProviderIntermediates, trialDirectory, bounded.Token) };
                if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                    return InterruptedHandoff(trial);
                trials.Add(rejected);
                progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                    rejected.Message, PreparationPhase.HandoffChecks, trialId, index,
                    rejected.FailureCode, rejected));
                return Failed(attempt, WorkerResultStanding.Failed,
                    "The conditioned Amber handoff failed correspondence, parameter or contact checks.",
                    "invalidHandoff", trials.ToImmutable(), trialId, index);
            }
            trial = trial with { DiagnosticArtifacts = await VerifiedMemgenDiagnosticsAsync(
                observed.ProviderIntermediates, trialDirectory, bounded.Token) };
            if (cancellationToken.IsCancellationRequested || bounded.IsCancellationRequested)
                return InterruptedHandoff(trial);
            derivation = derivation! with { Trials = trials.ToImmutable().Add(trial) };
            constructed = constructed with { Derivation = derivation };
            trials.Add(trial);
            progress?.Report(State(attempt.Id, StageExecutionStanding.Running,
                "The complete provider trial passed product handoff checks.",
                PreparationPhase.HandoffChecks, trialId, index, trial: trial));
            return new PreparationStartResult(attempt, derivation, constructed,
                State(attempt.Id, StageExecutionStanding.Running,
                    "The checked full-provider handoff is continuing to required minimization.",
                    PreparationPhase.HandoffChecks, trialId, index), trials.ToImmutable());
        }
        return Failed(attempt, WorkerResultStanding.Failed,
            "The declared full-provider geometry trials were exhausted without a checked handoff.",
            "geometryRetriesExhausted", trials.ToImmutable());
    }

    private static ConstructionTrialSummary? FailedTrialDetail(JsonElement? details,
        string trialId, int trialIndex, double lateral, double aqueous, string? failureCode)
    {
        if (details is not { ValueKind: JsonValueKind.Object } value ||
            !value.TryGetProperty("trial", out var encoded) || encoded.ValueKind != JsonValueKind.Object)
            return null;
        try
        {
            var trial = encoded.Deserialize<ConstructionTrialSummary>(
                new JsonSerializerOptions(JsonSerializerDefaults.Web));
            return trial is { Standing: ConstructionTrialStanding.Failed } &&
                trial.TrialId == trialId && trial.TrialIndex == trialIndex &&
                trial.LateralPaddingAngstrom == lateral && trial.AqueousPaddingAngstrom == aqueous &&
                trial.FailureCode == failureCode
                    ? trial with { DiagnosticArtifacts = ImmutableArray<TrialDiagnosticArtifact>.Empty } : null;
        }
        catch (JsonException) { return null; }
    }

    private static string? RetryAxis(string? failureCode, JsonElement? details)
    {
        if (failureCode == "zeroRoundedSpecies") return "lateral";
        if (failureCode == "chargeBudgetRefusal") return "aqueous";
        if (failureCode != "insufficientCellClearance" ||
            details is not { ValueKind: JsonValueKind.Object } value ||
            !value.TryGetProperty("axis", out var axis) || axis.ValueKind != JsonValueKind.String)
            return null;
        return axis.GetString() is "lateral" or "aqueous" ? axis.GetString() : null;
    }

    private static bool SameConditionAccount(ConstructionConditionAccount first,
        ConstructionConditionAccount second) =>
        first with { AqueousRegions = default } == second with { AqueousRegions = default } &&
        !first.AqueousRegions.IsDefault && !second.AqueousRegions.IsDefault &&
        first.AqueousRegions.SequenceEqual(second.AqueousRegions);

    private static PreparationStartResult Refusal(string reason) =>
        new(null, null, null, State(string.Empty, StageExecutionStanding.ResourceRefused, reason));

    private static PreparationStartResult NoAttempt(string reason) =>
        new(null, null, null, State(string.Empty, StageExecutionStanding.Failed, reason,
            failureCode: "preAdmissionRefused"));

    private static PreparationStartResult Failed(PreparationAttempt attempt, WorkerResultStanding standing,
        string reason, string? failureCode = null,
        ImmutableArray<ConstructionTrialSummary> trials = default,
        string? trialId = null, int? trialIndex = null,
        PreparationPhase? phase = null) =>
        new(attempt, null, null, State(attempt.Id, standing switch
        {
            WorkerResultStanding.Stopped => StageExecutionStanding.Stopped,
            WorkerResultStanding.Unobserved => StageExecutionStanding.Unobserved,
            WorkerResultStanding.Failed when failureCode is "resourceLimit" or "resourceRefused" =>
                StageExecutionStanding.ResourceRefused,
            _ => StageExecutionStanding.Failed
        }, reason, phase, trialId, trialIndex, failureCode), trials);

    private static StageExecutionState State(string attemptId, StageExecutionStanding standing,
        string message, PreparationPhase? phase = null, string? trialId = null,
        int? trialIndex = null, string? failureCode = null,
        ConstructionTrialSummary? trial = null) =>
        new(attemptId, null, null, standing, message, null, DateTimeOffset.UtcNow,
            phase, trialId, trialIndex, failureCode, trial);

    private static bool Available(string workingDirectory, string? path) =>
        !string.IsNullOrWhiteSpace(path) &&
        File.Exists(Path.IsPathRooted(path) ? path : Path.Combine(workingDirectory, path));

    private static bool HashMatches(string? path, string? expected)
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
        !string.IsNullOrWhiteSpace(p.ProviderVersion) &&
        p.PositiveIonArgument == "Na+" &&
        p.NegativeIonArgument == "Cl-" &&
        double.IsFinite(p.MinimumPaddingNanometers) && p.MinimumPaddingNanometers > 0 &&
        p.WaterMolarityForIonRounding == 55.4 &&
        p.MaximumAtomCount > 0 &&
        double.IsFinite(p.MaximumCellDimensionAngstrom) && p.MaximumCellDimensionAngstrom > 0 &&
        p.MaximumConstructionSeconds is > 0 and <= 86400 &&
        !string.IsNullOrWhiteSpace(p.ApproximationStatement) &&
        (p.Route switch
        {
            ConstructionRouteKind.NativeOpenMm => ValidNativeConstructionPolicy(p),
            ConstructionRouteKind.PackmolMemgen => ValidMemgenConstructionPolicy(p),
            _ => false
        });

    private static bool ValidNativeConstructionPolicy(ConstructionPolicy p) =>
        p.ProviderName == "OpenMM Modeller.addMembrane" &&
        p.SaltConvention == SaltConventionKind.NativeAddedSaltPlusNeutralization &&
        p.Memgen is null && p.NativePatchSha256 is { Length: 64 } &&
        p.NativePatchSha256.All(Uri.IsHexDigit) &&
        p.LipidTypeArgument is "DLPC" or "DLPE" or "DMPC" or "DOPC" or "DPPC" or
            "POPC" or "POPE" &&
        ((p.NativePatchMode is null or "installed") && p.NativeSourcePatchPath is null &&
            p.NativeSourcePatchSha256 is null && p.RemovedNativeLipidResidueIds is null ||
         (p.NativePatchMode == "mapped-lipid21-zenodo-popc" && p.LipidTypeArgument == "POPC" &&
          p.RemovedNativeLipidResidueIds is null && !string.IsNullOrWhiteSpace(p.NativeSourcePatchPath) &&
          HashLike(p.NativeSourcePatchSha256)));

    private static bool ValidMemgenConstructionPolicy(ConstructionPolicy p)
    {
        var m = p.Memgen;
        return p.ProviderName == "PACKMOL-Memgen" && p.ProviderVersion == "2026.3.25" &&
            p.SaltConvention == SaltConventionKind.MemgenChargeCompensated &&
            p.NativePatchPath is null && p.NativePatchSha256 is null &&
            p.LipidTypeArgument is null && p.NativePatchMode is null &&
            p.NativeSourcePatchPath is null && p.NativeSourcePatchSha256 is null &&
            p.RemovedNativeLipidResidueIds is null &&
            !p.ProviderAssets.IsDefaultOrEmpty && p.ProviderAssets.All(asset =>
                !string.IsNullOrWhiteSpace(asset.Id) && !string.IsNullOrWhiteSpace(asset.Version) &&
                !string.IsNullOrWhiteSpace(asset.Path) && HashLike(asset.Sha256)) &&
            p.ProviderAssets.Select(asset => asset.Id).Distinct(StringComparer.Ordinal).Count() ==
                p.ProviderAssets.Length &&
            p.CoveredSpeciesIds is { IsDefaultOrEmpty: false } &&
            p.CoveredSpeciesIds.ToHashSet(StringComparer.Ordinal).SetEquals(
                ["POPC", "POPE", "DLPC", "DLPE", "DMPC", "DOPC", "DPPC", "CHL1",
                    "HOH", "NA", "CL"]) &&
            p.MaximumAtomCount == 120000 && p.MaximumCellDimensionAngstrom == 180 &&
            p.MaximumConstructionSeconds == 3000 &&
            m is { Engine: "sander", ProteinForceField: "ff19SB", LipidForceField: "lipid21",
                WaterForceField: "tip3p", LocalPackingLoops: 20, TotalPackingLoops: 100,
                PackingOptimizerIterations: 20, ConditioningSteepestDescentSteps: 250,
                ConditioningConjugateGradientSteps: 250, MaximumGeometryRetries: 2,
                Preoriented: true, DoNotProtonate: true, DoNotTrim: true, DoNotCenterXy: true,
                RetainIntermediates: true, UseRatio: true, UsePbc: true, UseTightBox: true,
                UseSalt: true, Parameterize: true, Condition: true } &&
            m.LateralPaddingAngstrom == 15 && m.AqueousPaddingAngstrom == 17.5 &&
            m.LeafletEnvelopeAngstrom == 23 && m.PackingToleranceAngstrom == 2.0 &&
            m.ConditioningRestraintKcalMolAngstromSquared == 10;
    }

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
            Math.Abs(observed.NetChargeElementary) > 1e-4)
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
            policy.Limitations, Trials: ImmutableArray<ConstructionTrialSummary>.Empty);
        reason = string.Empty;
        return true;
    }

    private static bool TryDeriveMemgen(PreparationAttempt attempt, StudyRevision revision,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        ApplicablePreparationPolicy policy, ConstructionObservations observed,
        ImmutableArray<ConstructionTrialSummary> trials,
        out ConstructionDerivation? derivation, out string cause, out string? clearanceAxis)
    {
        derivation = null;
        cause = "invalidProviderAccount";
        clearanceAxis = null;
        var selected = observed.Trial;
        var account = observed.Conditions;
        if (selected is null || account is null || selected.Standing != ConstructionTrialStanding.Checked ||
            selected.ProposedLipidCounts.IsDefaultOrEmpty ||
            selected.AchievedLipidCounts.IsDefaultOrEmpty ||
            selected.CleanupRemovedLipidCounts.IsDefault ||
            selected.ProposedCellAngstrom.Length != 3 ||
            selected.ActualCellAngstrom.Length != 3 ||
            observed.ActualCellAngstrom.Length != 3 ||
            observed.ProteinPeriodicImageGapsAngstrom.Length != 3 ||
            !selected.ActualCellAngstrom.SequenceEqual(observed.ActualCellAngstrom) ||
            selected.ActualCellAngstrom.Any(value => !double.IsFinite(value) || value <= 0 ||
                value > policy.Construction.MaximumCellDimensionAngstrom) ||
            selected.ProposedCellAngstrom.Any(value => !double.IsFinite(value) || value <= 0) ||
            observed.ProteinPeriodicImageGapsAngstrom.Any(value => !double.IsFinite(value) || value < 0) ||
            selected.ActualCellAngstrom.Any(value =>
                value <= 20 * policy.SystemSettings.NonbondedCutoffNanometers) ||
            observed.AtomCount <= protein.Molecule.AtomCount ||
            observed.AtomCount > policy.Construction.MaximumAtomCount ||
            !double.IsFinite(observed.InitialPotentialEnergyKjMol) ||
            !double.IsFinite(observed.NetChargeElementary) ||
            Math.Abs(observed.NetChargeElementary) > 1e-4)
            return false;
        for (var axis = 0; axis < 3; axis++)
        {
            var required = 2 * (axis == 2 ? selected.AqueousPaddingAngstrom :
                selected.LateralPaddingAngstrom);
            if (observed.ProteinPeriodicImageGapsAngstrom[axis] +
                policy.ExportCellLengthReadBackToleranceAngstrom < required)
            {
                cause = "insufficientCellClearance";
                clearanceAxis = axis == 2 ? "aqueous" : "lateral";
                return false;
            }
        }

        var targets = new[] { membrane.Intended.Upper, membrane.Intended.Lower }
            .SelectMany(side => side.Fractions.Where(fraction => fraction.Fraction > 0)
                .Select(fraction => (side.PhysicalSide, fraction.SpeciesId, fraction.Fraction)))
            .ToDictionary(item => (item.PhysicalSide, item.SpeciesId), item => item.Fraction);
        var proposed = selected.ProposedLipidCounts;
        var achieved = selected.AchievedLipidCounts;
        var removed = selected.CleanupRemovedLipidCounts;
        if (proposed.Length != targets.Count || achieved.Length != targets.Count ||
            proposed.Select(item => (item.PhysicalSide, item.SpeciesId)).Distinct().Count() != proposed.Length ||
            achieved.Select(item => (item.PhysicalSide, item.SpeciesId)).Distinct().Count() != achieved.Length ||
            removed.Select(item => (item.PhysicalSide, item.SpeciesId)).Distinct().Count() != removed.Length ||
            proposed.Any(item => !targets.TryGetValue((item.PhysicalSide, item.SpeciesId), out var fraction) ||
                !double.IsFinite(item.IntendedFraction) ||
                Math.Abs(item.IntendedFraction - fraction) > 1e-12 || item.Count <= 0) ||
            achieved.Any(item => !targets.TryGetValue((item.PhysicalSide, item.SpeciesId), out var fraction) ||
                !double.IsFinite(item.IntendedFraction) ||
                Math.Abs(item.IntendedFraction - fraction) > 1e-12 || item.Count <= 0))
            return false;
        foreach (var current in achieved)
        {
            var starting = proposed.Single(item => item.PhysicalSide == current.PhysicalSide &&
                item.SpeciesId == current.SpeciesId).Count;
            var cleanup = removed.Where(item => item.PhysicalSide == current.PhysicalSide &&
                item.SpeciesId == current.SpeciesId).Sum(item => item.Count);
            if (cleanup < 0 || starting - current.Count != cleanup) return false;
        }
        if (removed.Any(item => !targets.TryGetValue((item.PhysicalSide, item.SpeciesId),
                out var fraction) || item.Count < 0 || !double.IsFinite(item.IntendedFraction) ||
                Math.Abs(item.IntendedFraction - fraction) > 1e-12))
            return false;

        var generated = observed.SpeciesCounts;
        if (generated.IsDefaultOrEmpty || generated.Any(item => item.Count < 0 ||
                item.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                item.Role is not (GeneratedComponentRoleKind.Lipid or GeneratedComponentRoleKind.Water or
                    GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon)) ||
            generated.Select(item => (item.Role, item.PhysicalSide, item.SpeciesId))
                .Distinct().Count() != generated.Length ||
            generated.Where(item => item.Role == GeneratedComponentRoleKind.Lipid)
                .Any(item => !targets.ContainsKey((item.PhysicalSide, item.SpeciesId))) ||
            achieved.Any(item => generated.Count(actual => actual.Role == GeneratedComponentRoleKind.Lipid &&
                actual.PhysicalSide == item.PhysicalSide && actual.SpeciesId == item.SpeciesId &&
                actual.Count == item.Count) != 1) ||
            generated.Any(item => item.Role switch
            {
                GeneratedComponentRoleKind.Water => item.SpeciesId != policy.Water.SpeciesId,
                GeneratedComponentRoleKind.PositiveIon => item.SpeciesId != policy.Sodium.SpeciesId,
                GeneratedComponentRoleKind.NegativeIon => item.SpeciesId != policy.Chloride.SpeciesId,
                _ => false
            }))
            return false;

        if (!MemgenConditionsCoherent(account, observed, generated, revision.Conditions.TargetNaClMolar,
                policy.Construction.WaterMolarityForIonRounding))
            return false;
        var representations = membrane.SpeciesRepresentations.ToDictionary(item => item.SpeciesId,
            StringComparer.Ordinal);
        var generatedAtoms = generated.Sum(item => (long)item.Count *
            (item.Role switch
            {
                GeneratedComponentRoleKind.Lipid => representations[item.SpeciesId].AtomCount,
                GeneratedComponentRoleKind.Water => policy.Water.AtomCount,
                GeneratedComponentRoleKind.PositiveIon => policy.Sodium.AtomCount,
                GeneratedComponentRoleKind.NegativeIon => policy.Chloride.AtomCount,
                _ => 0
            }));
        if (generatedAtoms + protein.Molecule.AtomCount != observed.AtomCount)
            return false;

        var pairMolar = Math.Min(account.SodiumAqueousMolar!.Value,
            account.ChlorideAqueousMolar!.Value);
        derivation = new ConstructionDerivation(attempt.Id, proposed,
            selected.ProposedCellAngstrom, account.FinalWaterCount,
            account.FinalSodiumCount, account.FinalChlorideCount,
            account.PreparedFormalChargeElementary, revision.Conditions.TargetNaClMolar,
            pairMolar, account.EstimatedAqueousVolumeAngstromCubed,
            ImmutableArray.Create(policy.Construction.ApproximationStatement,
                "Memgen derived population and cell from its selected ratio method; actual post-cleanup populations and aqueous estimates are reported separately."),
            policy.Limitations, trials, selected.TrialId, account, observed.AmberImport);
        cause = string.Empty;
        return true;
    }

    private static bool MemgenConditionsCoherent(ConstructionConditionAccount account,
        ConstructionObservations observed, ImmutableArray<ObservedSpeciesCount> generated,
        double nominalMolar, double waterMolarity)
    {
        if (account.SaltBranch is not ("chargeCompensated" or "neutralizationOnly") ||
            account.AqueousRegions.IsDefault || account.AqueousRegions.Length != 2 ||
            account.AqueousRegions.Select(region => region.Side).ToHashSet()
                .SetEquals([LeafletSide.Upper, LeafletSide.Lower]) == false ||
            account.AqueousRegions.Any(region =>
                !double.IsFinite(region.EstimatedVolumeAngstromCubed) ||
                region.EstimatedVolumeAngstromCubed <= 0 || region.FlooredNominalSaltCount < 0 ||
                region.ProviderGeneratedSodiumCount < 0 || region.ProviderGeneratedChlorideCount < 0) ||
            !double.IsFinite(account.EstimatedAqueousVolumeAngstromCubed) ||
            account.EstimatedAqueousVolumeAngstromCubed <= 0 ||
            Math.Abs(account.EstimatedAqueousVolumeAngstromCubed -
                account.AqueousRegions.Sum(region => region.EstimatedVolumeAngstromCubed)) > 1e-6 ||
            !double.IsFinite(account.PreparedFormalChargeElementary) ||
            !double.IsFinite(observed.ProteinNetChargeElementary) ||
            Math.Abs(account.PreparedFormalChargeElementary -
                observed.ProteinNetChargeElementary) > 1e-6 ||
            !double.IsFinite(account.ProviderResidueNameChargeElementary) ||
            !double.IsFinite(account.ChargePdbDeltaElementary) ||
            Math.Abs(account.PreparedFormalChargeElementary -
                account.ProviderResidueNameChargeElementary - account.ChargePdbDeltaElementary) > 1e-6 ||
            !double.IsFinite(account.FinalNetChargeElementary) ||
            Math.Abs(account.FinalNetChargeElementary) > 1e-4 ||
            Math.Abs(account.FinalNetChargeElementary - observed.NetChargeElementary) > 1e-6 ||
            account.RetainedWaterCount < 0 || account.ProviderGeneratedWaterCount < 0 ||
            account.FinalWaterCount <= 0 || account.RetainedSodiumCount < 0 ||
            account.RetainedChlorideCount < 0 || account.ProviderGeneratedSodiumCount < 0 ||
            account.ProviderGeneratedChlorideCount < 0 || account.LeapAddedSodiumCount < 0 ||
            account.LeapAddedChlorideCount < 0 || account.LeapRemovedGeneratedWaterCount < 0 ||
            account.LeapRemovedGeneratedSodiumCount < 0 ||
            account.LeapRemovedGeneratedChlorideCount < 0 ||
            account.LeapRemovedRetainedWaterCount != 0 ||
            account.LeapRemovedRetainedSodiumCount != 0 ||
            account.LeapRemovedRetainedChlorideCount != 0 ||
            account.FinalSodiumCount < 0 || account.FinalChlorideCount < 0 ||
            account.FinalWaterCount != observed.WaterCount ||
            account.FinalSodiumCount != observed.PositiveIonCount ||
            account.FinalChlorideCount != observed.NegativeIonCount ||
            account.FinalWaterCount != account.RetainedWaterCount +
                account.ProviderGeneratedWaterCount - account.LeapRemovedGeneratedWaterCount ||
            account.FinalSodiumCount != account.RetainedSodiumCount +
                account.ProviderGeneratedSodiumCount + account.LeapAddedSodiumCount -
                account.LeapRemovedGeneratedSodiumCount ||
            account.FinalChlorideCount != account.RetainedChlorideCount +
                account.ProviderGeneratedChlorideCount + account.LeapAddedChlorideCount -
                account.LeapRemovedGeneratedChlorideCount ||
            generated.Where(item => item.Role == GeneratedComponentRoleKind.Water)
                .Sum(item => item.Count) != account.FinalWaterCount - account.RetainedWaterCount ||
            generated.Where(item => item.Role == GeneratedComponentRoleKind.PositiveIon)
                .Sum(item => item.Count) != account.FinalSodiumCount - account.RetainedSodiumCount ||
            generated.Where(item => item.Role == GeneratedComponentRoleKind.NegativeIon)
                .Sum(item => item.Count) != account.FinalChlorideCount - account.RetainedChlorideCount ||
            account.AqueousRegions.Sum(region => region.ProviderGeneratedSodiumCount) !=
                account.ProviderGeneratedSodiumCount ||
            account.AqueousRegions.Sum(region => region.ProviderGeneratedChlorideCount) !=
                account.ProviderGeneratedChlorideCount)
            return false;
        foreach (var region in account.AqueousRegions)
        {
            // Memgen computes salt from its unrounded water-box volume, then
            // prints that volume with two decimal places. The logged value
            // identifies a 0.01 Å³ interval, not the exact floor operand.
            var logged = region.EstimatedVolumeAngstromCubed;
            var lowerFloor = Math.Floor(nominalMolar * Math.Max(0, logged - 0.005) *
                MemgenMoleculesPerMolarAngstromCubed);
            var upperFloor = Math.Floor(nominalMolar * (logged + 0.005) *
                MemgenMoleculesPerMolarAngstromCubed);
            if (region.FlooredNominalSaltCount < lowerFloor ||
                region.FlooredNominalSaltCount > upperFloor) return false;
        }
        // Memgen receives an integral residue-name estimate and an int-parsed
        // charge delta. The Amber preflight charge may carry small floating noise.
        var providerEstimate = Math.Round(account.ProviderResidueNameChargeElementary);
        var providerDelta = Math.Round(account.ChargePdbDeltaElementary);
        if (Math.Abs(account.ProviderResidueNameChargeElementary - providerEstimate) > 1e-4 ||
            Math.Abs(account.ChargePdbDeltaElementary - providerDelta) > 1e-4)
            return false;
        var expectedBranch = Math.Abs(providerEstimate + providerDelta) / 2 <
            account.AqueousRegions.Min(region => region.FlooredNominalSaltCount)
                ? "chargeCompensated" : "neutralizationOnly";
        if (account.SaltBranch != expectedBranch) return false;
        var volume = account.EstimatedAqueousVolumeAngstromCubed;
        var sodium = account.FinalSodiumCount / (MoleculesPerMolarAngstromCubed * volume);
        var chloride = account.FinalChlorideCount / (MoleculesPerMolarAngstromCubed * volume);
        var finiteWaterSodium = waterMolarity * account.FinalSodiumCount / account.FinalWaterCount;
        var finiteWaterChloride = waterMolarity * account.FinalChlorideCount / account.FinalWaterCount;
        return SameEstimate(account.SodiumAqueousMolar, sodium) &&
            SameEstimate(account.ChlorideAqueousMolar, chloride) &&
            SameEstimate(account.SodiumFiniteWaterMolar, finiteWaterSodium) &&
            SameEstimate(account.ChlorideFiniteWaterMolar, finiteWaterChloride);
    }

    private static bool SameEstimate(double? reported, double calculated) =>
        reported is double value && double.IsFinite(value) &&
        Math.Abs(value - calculated) <= 1e-8 * Math.Max(1, Math.Abs(calculated));

    private static async Task<ConstructedExplicitSystem?> AssessConstructedMemgenAsync(
        PreparationAttempt attempt, StudyRevision revision, AssessedPreparedProtein protein,
        AssessedMembraneModel membrane, AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy, ConstructionDerivation derivation,
        WorkerResult<ConstructionObservations> result, CancellationToken cancellationToken)
    {
        var observed = result.Observations!;
        var trial = observed.Trial!;
        var account = observed.Conditions!;
        if (result.Provider?.Name != policy.Construction.ProviderName ||
            result.Provider.Version != policy.Construction.ProviderVersion ||
            observed.NativePatchSha256 is not null || observed.NativePatchMode is not null ||
            observed.NativeSourcePatchSha256 is not null ||
            observed.CorrespondedResultAtomCount != observed.AtomCount ||
            !observed.ProteinIdentityAndBondsPreserved ||
            !double.IsFinite(observed.MaximumProteinCoordinateDeviationAngstrom) ||
            observed.MaximumProteinCoordinateDeviationAngstrom < 0 ||
            !observed.ContactWarnings.IsDefaultOrEmpty ||
            !observed.GeometryWarnings.IsDefaultOrEmpty ||
            !observed.ParameterWarnings.IsDefaultOrEmpty ||
            observed.OutputFrameMidplaneZAngstrom is not double outputMidplane ||
            !double.IsFinite(outputMidplane) ||
            !policy.LocalStateObservation.RequiredMetricNames.Contains("upperLipidHeadMeanZAngstrom") ||
            !policy.LocalStateObservation.RequiredMetricNames.Contains("lowerLipidHeadMeanZAngstrom") ||
            !LocalStateSupportsConstruction(observed.LocalState, policy) ||
            result.Artifacts.IsDefaultOrEmpty ||
            result.Artifacts.Select(item => item.Role).Distinct(StringComparer.Ordinal).Count() !=
                result.Artifacts.Length ||
            observed.ProviderIntermediates.IsDefaultOrEmpty)
            return null;
        var upperHead = observed.LocalState.Measurements.Single(item =>
            item.Name == "upperLipidHeadMeanZAngstrom");
        var lowerHead = observed.LocalState.Measurements.Single(item =>
            item.Name == "lowerLipidHeadMeanZAngstrom");
        var separation = observed.LocalState.Measurements.Single(item =>
            item.Name == "leafletHeadSeparationAngstrom");
        if (upperHead.Value <= outputMidplane || lowerHead.Value >= outputMidplane ||
            separation.Value <= 0 ||
            Math.Abs(separation.Value - (upperHead.Value - lowerHead.Value)) > 1e-6)
            return null;
        string[] handoffRoles = ["topologyCif", "topologyJson", "systemXml", "stateXml",
            "correspondenceJson"];
        string[] providerRoles = ["providerOptions", "providerLog", "providerCommand",
            "providerStdout", "preparedChargeLeapInput", "preparedChargeLeapLog",
            "preparedChargeLeapTopology", "preparedChargeLeapCoordinates",
            "preparedChargePdb",
            "providerInput", "providerProtein",
            "providerPackmolInput", "providerPackmolLog", "providerPacked",
            "providerLeapInput", "providerLeapLog", "providerParameterized",
            "amberCoordinates", "amberTopology",
            "providerRestrainedInput", "providerRestrainedLog", "providerRestrainedRestart",
            "providerUnrestrainedInput", "providerUnrestrainedLog", "amberFinalRestart"];
        var artifacts = result.Artifacts.ToDictionary(item => item.Role, StringComparer.Ordinal);
        if (handoffRoles.Concat(providerRoles).Any(role => !artifacts.ContainsKey(role)) ||
            observed.ProviderIntermediates.Select(item => item.Role)
                .Distinct(StringComparer.Ordinal).Count() != observed.ProviderIntermediates.Length ||
            providerRoles.Any(role => observed.ProviderIntermediates.Count(item =>
                item.Role == role && artifacts.TryGetValue(role, out var topLevel) &&
                item.Path == topLevel.Path && item.Sha256 == topLevel.Sha256) != 1) ||
            observed.ProviderIntermediates.Any(item =>
                !artifacts.TryGetValue(item.Role, out var topLevel) ||
                item.Path != topLevel.Path || item.Sha256 != topLevel.Sha256))
            return null;
        var imported = observed.AmberImport;
        if (imported is null ||
            !string.Equals(imported.PrmtopSha256, artifacts["amberTopology"].Sha256,
                StringComparison.OrdinalIgnoreCase) ||
            !string.Equals(imported.FinalRestartSha256, artifacts["amberFinalRestart"].Sha256,
                StringComparison.OrdinalIgnoreCase) ||
            imported.AmberAtomCount != observed.AtomCount ||
            imported.ImportedAtomCount != observed.AtomCount ||
            !double.IsFinite(imported.AmberNetChargeElementary) ||
            !double.IsFinite(imported.ImportedNetChargeElementary) ||
            Math.Abs(imported.AmberNetChargeElementary) > 1e-4 ||
            Math.Abs(imported.ImportedNetChargeElementary) > 1e-4 ||
            Math.Abs(imported.AmberNetChargeElementary - imported.ImportedNetChargeElementary) > 1e-6 ||
            !double.IsFinite(imported.MaximumCellVectorDeviationAngstrom) ||
            imported.MaximumCellVectorDeviationAngstrom < 0 ||
            imported.MaximumCellVectorDeviationAngstrom >
                policy.ExportCellLengthReadBackToleranceAngstrom ||
            imported.AmberCmapTermCount < 0 ||
            imported.ImportedCmapTermCount != imported.AmberCmapTermCount ||
            imported.ImportedForceKinds.IsDefaultOrEmpty ||
            !imported.ImportedForceKinds.Contains("NonbondedForce") ||
            imported.AmberCmapTermCount > 0 &&
                !imported.ImportedForceKinds.Contains("CMAPTorsionForce") ||
            !imported.AtomOrderPreserved || !imported.BondedTermsPreserved ||
            !imported.NonbondedTermsPreserved || !imported.ExclusionsPreserved ||
            !imported.UnitsPreserved || !imported.ParameterComparisonPerformed ||
            !HashLike(imported.ParameterCorrespondenceSha256))
            return null;
        foreach (var artifact in handoffRoles.Concat(providerRoles)
            .Concat(artifacts.ContainsKey("providerPostCleanup") ? ["providerPostCleanup"] : [])
            .Select(role => artifacts[role]))
            if (!await ArtifactMatchesAsync(artifact, cancellationToken)) return null;

        var coordinate = artifacts["topologyCif"];
        var topology = artifacts["topologyJson"];
        var system = artifacts["systemXml"];
        var state = artifacts["stateXml"];
        var mapping = artifacts["correspondenceJson"];
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
            correspondence.ResultId != coordinate.Sha256 ||
            correspondence.Atoms.IsDefault || correspondence.Atoms.Length != observed.AtomCount ||
            correspondence.Atoms.Select(atom => atom.ResultAtomIndex).Distinct().Count() !=
                observed.AtomCount ||
            correspondence.Atoms.Any(atom => atom.ResultAtomIndex < 0 ||
                atom.ResultAtomIndex >= observed.AtomCount ||
                string.IsNullOrWhiteSpace(atom.ResultAtomId) ||
                !Enum.IsDefined(atom.Role) || !Enum.IsDefined(atom.MoleculeRole) ||
                !Enum.IsDefined(atom.AtomRole)) ||
            correspondence.Atoms.Select(atom => atom.ResultAtomId)
                .Distinct(StringComparer.Ordinal).Count() != observed.AtomCount)
            return null;

        var expectedRetained = protein.Correspondence.Atoms;
        var actualRetained = correspondence.Atoms.Where(atom =>
            atom.GeneratedComponentRole is null).ToArray();
        if (actualRetained.Length != expectedRetained.Length ||
            expectedRetained.Any(atom => atom.GeneratedComponentRole is not null) ||
            !SameRetainedAtoms(expectedRetained, actualRetained))
            return null;

        var expectedGenerated = observed.SpeciesCounts.ToDictionary(
            item => (item.Role, item.PhysicalSide, item.SpeciesId),
            item => (long)item.Count * (item.Role switch
            {
                GeneratedComponentRoleKind.Lipid => membrane.SpeciesRepresentations
                    .Single(species => species.SpeciesId == item.SpeciesId).AtomCount,
                GeneratedComponentRoleKind.Water => policy.Water.AtomCount,
                GeneratedComponentRoleKind.PositiveIon => policy.Sodium.AtomCount,
                GeneratedComponentRoleKind.NegativeIon => policy.Chloride.AtomCount,
                _ => 0
            }));
        var actualGenerated = new Dictionary<(GeneratedComponentRoleKind Role,
            LeafletSide Side, string SpeciesId), long>();
        foreach (var atom in correspondence.Atoms.Where(item => item.GeneratedComponentRole is not null))
        {
            if (atom.Role != AtomOriginKind.Generated || atom.SourceAtomId is not null ||
                atom.SourceResidue is not null || atom.ApprovedChangeId is not null ||
                atom.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                atom.GeneratedComponentRole is not (GeneratedComponentRoleKind.Lipid or
                    GeneratedComponentRoleKind.Water or GeneratedComponentRoleKind.PositiveIon or
                    GeneratedComponentRoleKind.NegativeIon) ||
                string.IsNullOrWhiteSpace(atom.GeneratedSpeciesId) ||
                atom.MoleculeRole != (atom.GeneratedComponentRole is
                    GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon
                    ? MoleculeRoleKind.Ion : atom.GeneratedComponentRole == GeneratedComponentRoleKind.Lipid
                        ? MoleculeRoleKind.Lipid : MoleculeRoleKind.Water))
                return null;
            var key = (atom.GeneratedComponentRole.Value, atom.PhysicalSide.Value,
                atom.GeneratedSpeciesId);
            actualGenerated[key] = actualGenerated.GetValueOrDefault(key) + 1;
        }
        if (actualGenerated.Count != expectedGenerated.Count ||
            actualGenerated.Any(item => !expectedGenerated.TryGetValue(item.Key, out var expected) ||
                item.Value != expected))
            return null;

        var molecule = new MolecularArtifact(coordinate.Sha256, coordinate.Path, coordinate.Sha256,
            topology.Path, system.Path, state.Path, observed.AtomCount,
            string.Join(" × ", observed.ActualCellAngstrom.Select(value =>
                value.ToString("G17", System.Globalization.CultureInfo.InvariantCulture))),
            topology.Sha256, mapping.Path, mapping.Sha256, system.Sha256, state.Sha256);
        var constructedId = "constructed-" + attempt.Id;
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"),
            constructedId, policy.Construction.ProviderName,
            "Complete Memgen construction, conditioning and Amber handoff",
            $"{observed.AtomCount} mapped atoms; actual {trial.AchievedLipidCounts.Length} species/leaflet counts; " +
            $"{account.FinalWaterCount} water, {account.FinalSodiumCount} sodium and " +
            $"{account.FinalChlorideCount} chloride ions; salt branch {account.SaltBranch}.",
            $"Attempt {attempt.Id}; trial {trial.TrialId}; policy {policy.Id}; " +
            $"Memgen {policy.Construction.ProviderVersion}; Amber topology " +
            artifacts["amberTopology"].Sha256,
            "The conditioned handoff is checked for final OpenMM minimization; provider-internal " +
            "conditioning is not a completed minimized stage or an equilibrium claim.",
            EvidenceBearing.Context));
        var findings = account.SaltBranch == "neutralizationOnly"
            ? ImmutableArray.Create(new ScientificFinding(Guid.NewGuid().ToString("N"),
                constructedId, evidence[0].Id,
                "Memgen's charge-compensated salt branch added neutralization only; the nominal 0.15 M input did not produce background salt pairs.",
                "The actual Na/Cl populations and aqueous-volume estimates are reported for this attempt.",
                FindingDisposition.Challenges, true, DateTimeOffset.UtcNow))
            : ImmutableArray<ScientificFinding>.Empty;
        var conditions = $"Fixed nominal pH {revision.Conditions.NominalPh:G6}; Memgen nominal " +
            $"charge-compensated NaCl input {revision.Conditions.TargetNaClMolar:G6} M; " +
            $"salt branch {account.SaltBranch}; actual Na+ {account.SodiumAqueousMolar:G6} M and " +
            $"Cl- {account.ChlorideAqueousMolar:G6} M by estimated aqueous volume; " +
            $"finite-water estimates Na+ {account.SodiumFiniteWaterMolar:G6} M and " +
            $"Cl- {account.ChlorideFiniteWaterMolar:G6} M; neither denominator is measured bulk concentration.";
        return new ConstructedExplicitSystem(constructedId, attempt, molecule, derivation,
            correspondence, trial.AchievedLipidCounts, evidence, findings, conditions,
            observed.LocalState, observed.ActualCellAngstrom, outputMidplane,
            observed.MaximumProteinCoordinateDeviationAngstrom);
    }

    private static bool SameRetainedAtoms(ImmutableArray<AtomCorrespondence> expected,
        IReadOnlyCollection<AtomCorrespondence> actual)
    {
        static (AtomOriginKind, string?, ResidueAddress?, string?, MoleculeRoleKind,
            AtomRoleKind, string) Key(AtomCorrespondence atom) =>
            (atom.Role, atom.SourceAtomId, atom.SourceResidue, atom.ApprovedChangeId,
                atom.MoleculeRole, atom.AtomRole, atom.Element);
        var expectedCounts = expected.GroupBy(Key).ToDictionary(group => group.Key,
            group => group.Count());
        var actualCounts = actual.GroupBy(Key).ToDictionary(group => group.Key,
            group => group.Count());
        return expectedCounts.Count == actualCounts.Count &&
            expectedCounts.All(item => actualCounts.GetValueOrDefault(item.Key) == item.Value);
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
        var constructedId = "constructed-" + attempt.Id;
        var evidence = ImmutableArray.Create(new ScientificEvidence(Guid.NewGuid().ToString("N"), constructedId,
            result.Provider?.Name ?? "local scientific worker", "Whole-system construction and parameter assessment",
            $"{observed.AtomCount} mapped atoms; {observed.WaterCount} water; {observed.PositiveIonCount} sodium and {observed.NegativeIonCount} chloride ions",
            $"Attempt {attempt.Id}; policy {policy.Id}; placement {placement.Id}; " +
            $"OpenMM build {policy.Construction.ProviderVersion}; native patch SHA-256 {observed.NativePatchSha256}",
            "Same-invocation provider cell, populations, protein preservation, combined parameters and initial contacts were checked. " +
            "This checked construction intermediate is continuing to required minimization under the accepted operation; completion has not yet been established.",
            EvidenceBearing.Supports));
        var conditions = $"Fixed nominal pH {revision.Conditions.NominalPh:G6}; intended NaCl {derivation.IntendedNaClMolar:G6} M; " +
                         $"finite-cell estimated NaCl {derivation.EstimatedNaClMolar:G6} M; " +
                         "thermal equilibration has not been established.";
        return new ConstructedExplicitSystem(constructedId, attempt, molecule, derivation,
            correspondence, derivation.LipidCounts, evidence, ImmutableArray<ScientificFinding>.Empty, conditions,
            observed.LocalState, observed.ActualCellAngstrom,
            MaximumProteinCoordinateDeviationAngstrom: observed.MaximumProteinCoordinateDeviationAngstrom);
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

    private static async Task<ImmutableArray<TrialDiagnosticArtifact>> VerifiedMemgenDiagnosticsAsync(
        ImmutableArray<WorkerArtifact> reported, string trialDirectory, CancellationToken cancellationToken)
    {
        if (reported.IsDefaultOrEmpty) return ImmutableArray<TrialDiagnosticArtifact>.Empty;
        var result = ImmutableArray.CreateBuilder<TrialDiagnosticArtifact>();
        string root;
        try { root = Path.GetFullPath(trialDirectory); }
        catch (Exception exception) when (exception is ArgumentException or NotSupportedException or PathTooLongException)
        { return result.ToImmutable(); }
        long totalBytes = 0;
        // The named provider roles bound both the count and the kinds of files exposed.
        foreach (var group in reported.Where(item => !string.IsNullOrWhiteSpace(item.Role) &&
                         MemgenDiagnosticRoles.ContainsKey(item.Role))
                     .GroupBy(item => item.Role, StringComparer.Ordinal))
        {
            if (group.Count() != 1) continue;
            var artifact = group.Single();
            if (!HashLike(artifact.Sha256) || string.IsNullOrWhiteSpace(artifact.Path)) continue;
            string path;
            long length;
            try
            {
                path = Path.GetFullPath(Path.IsPathRooted(artifact.Path)
                    ? artifact.Path : Path.Combine(root, artifact.Path));
                var relative = Path.GetRelativePath(root, path);
                if (Path.IsPathRooted(relative) || relative == ".." ||
                    relative.StartsWith($"..{Path.DirectorySeparatorChar}", StringComparison.Ordinal) ||
                    relative == ".") continue;
                // A lexical child through a symlink can resolve outside the attempt.
                var part = Path.GetPathRoot(root)!;
                var linked = (File.GetAttributes(part) & FileAttributes.ReparsePoint) != 0;
                foreach (var segment in root[part.Length..].Split(Path.DirectorySeparatorChar,
                    StringSplitOptions.RemoveEmptyEntries))
                {
                    part = Path.Combine(part, segment);
                    if ((File.GetAttributes(part) & FileAttributes.ReparsePoint) != 0)
                    { linked = true; break; }
                }
                if (linked) continue;
                foreach (var segment in relative.Split(Path.DirectorySeparatorChar))
                {
                    part = Path.Combine(part, segment);
                    if ((File.GetAttributes(part) & FileAttributes.ReparsePoint) != 0)
                    { linked = true; break; }
                }
                if (linked || !File.Exists(path)) continue;
                length = new FileInfo(path).Length;
                if (length > MaximumDiagnosticFileBytes ||
                    length > MaximumDiagnosticTotalBytes - totalBytes) continue;
            }
            catch (Exception exception) when (exception is ArgumentException or NotSupportedException or
                PathTooLongException or IOException or UnauthorizedAccessException)
            { continue; }
            try
            {
                if (!await ArtifactMatchesAsync(artifact with { Path = path }, cancellationToken)) continue;
            }
            catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
            { return ImmutableArray<TrialDiagnosticArtifact>.Empty; }
            totalBytes += length;
            result.Add(new TrialDiagnosticArtifact(artifact.Role, artifact.Sha256.ToLowerInvariant(),
                Path.GetFileName(path), MemgenDiagnosticRoles[artifact.Role], LocalPath: path));
        }
        return result.ToImmutable();
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
            policy.ContactCriteria.IsDefault ||
            policy.ConstructionCriteria.Select(item => item.MeasurementName)
                .Distinct(StringComparer.Ordinal).Count() != policy.ConstructionCriteria.Length ||
            policy.ConstructionCriteria.Any(item => item.MeasurementName is
                "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom") ||
            !policy.ConstructionCriteria.Any(item =>
                item.MeasurementName == "minimumIntermolecularHeavyAtomDistanceAngstrom" &&
                item.Minimum is double minimum && double.IsFinite(minimum) && minimum == 1.5) ||
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
            policy.ContactCriteria.Any(item => item.StageKind == StageKind.Equilibration &&
                policy.OptionalEquilibration is null))
            return false;
        return policy.ConstructionCriteria.All(criterion =>
            spec.RequiredMetricNames.Contains(criterion.MeasurementName) &&
            !string.IsNullOrWhiteSpace(criterion.Unit) && !string.IsNullOrWhiteSpace(criterion.Scope) &&
            (criterion.Minimum is null || double.IsFinite(criterion.Minimum.Value)) &&
            (criterion.Maximum is null || double.IsFinite(criterion.Maximum.Value)));
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
