using System.Collections.Immutable;
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

    public async Task<PreparationStartResult> StartAsync(
        PreparationAttempt attempt,
        StudyRevision revision,
        AssessedPreparedProtein protein,
        AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy,
        string packmolExecutablePath,
        string packmolVersion,
        string packmolExecutableSha256,
        string workingDirectory,
        Action<PreparationAttempt>? onAccepted,
        IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
    {
        if (attempt.StudyRevisionId != revision.Id || attempt.ProteinId != protein.Id ||
            attempt.MembraneId != membrane.Id || attempt.PlacementId != placement.Id ||
            attempt.PolicyId != policy.Id || revision.Id != protein.StudyRevisionId || revision.Id != membrane.StudyRevisionId ||
            revision.Id != placement.StudyRevisionId || revision.IntendedProtein?.Id != protein.Intended.Id ||
            revision.Membrane?.Id != membrane.Intended.Id ||
            revision.AdoptedPlacementProposalId != placement.Proposal.Id ||
            placement.Standing != AssessmentStanding.Supported ||
            placement.Proposal.PreparedProteinId != protein.Id ||
            placement.Proposal.MembraneModelId != membrane.Intended.Id ||
            string.IsNullOrWhiteSpace(policy.Id) || string.IsNullOrWhiteSpace(policy.Version) ||
            policy.EvidenceReferences.IsDefaultOrEmpty || policy.Construction.EvidenceReferences.IsDefaultOrEmpty ||
            policy.Construction.CoveredTopologyKinds.IsDefaultOrEmpty ||
            !policy.Construction.CoveredTopologyKinds.Contains(placement.Proposal.TopologyKind) ||
            policy.ForceFieldFiles.IsDefaultOrEmpty ||
            policy.SystemSettings is null ||
            policy.SystemSettings.NonbondedMethod != "PME" ||
            policy.SystemSettings.Constraints != "HBonds" ||
            !double.IsFinite(policy.SystemSettings.NonbondedCutoffNanometers) ||
            policy.SystemSettings.NonbondedCutoffNanometers <= 0 ||
            !double.IsFinite(policy.SystemSettings.EwaldErrorTolerance) ||
            policy.SystemSettings.EwaldErrorTolerance <= 0 ||
            policy.SystemSettings.EwaldErrorTolerance >= 1 ||
            (policy.SystemSettings.SwitchDistanceNanometers is double switchDistance &&
                (!double.IsFinite(switchDistance) || switchDistance <= 0 ||
                 switchDistance >= policy.SystemSettings.NonbondedCutoffNanometers)) ||
            (policy.SystemSettings.HydrogenMassDaltons is double hydrogenMass &&
                (!double.IsFinite(hydrogenMass) || hydrogenMass <= 0)))
            return NoAttempt("Corresponding supported inputs and an identified, evidence-backed preparation policy are required.");
        var needed = membrane.SpeciesRepresentations.Append(policy.Water).Append(policy.Sodium).Append(policy.Chloride).ToArray();
        if (policy.Water.Category != "water" || policy.Sodium.Category != "ion" ||
            policy.Chloride.Category != "ion" ||
            Math.Abs(policy.Water.NetChargeElementary) > 1e-8 ||
            Math.Abs(policy.Sodium.NetChargeElementary - 1.0) > 1e-8 ||
            Math.Abs(policy.Chloride.NetChargeElementary + 1.0) > 1e-8 ||
            needed.Any(species => !policy.Construction.CoveredSpeciesIds.Contains(species.SpeciesId) ||
                string.IsNullOrWhiteSpace(species.ForceFieldFamily) ||
                string.IsNullOrWhiteSpace(species.ForceFieldVersion) || species.AtomCount <= 0 ||
                string.IsNullOrWhiteSpace(species.TemplateSha256) ||
                !Available(workingDirectory, species.CoordinateTemplatePath)) ||
            !File.Exists(protein.Molecule.CoordinatePath) ||
            !protein.Correspondence.Complete ||
            protein.Correspondence.ResultId != protein.Molecule.CoordinateSha256 ||
            protein.Correspondence.Atoms.Length != protein.Molecule.AtomCount ||
            protein.Correspondence.Atoms.Select(atom => atom.ResultAtomIndex).Distinct().Count() !=
                protein.Molecule.AtomCount ||
            !File.Exists(placement.Proposal.OrientedProtein.CoordinatePath) ||
            string.IsNullOrWhiteSpace(placement.Proposal.OrientedProtein.TopologyPath) ||
            string.IsNullOrWhiteSpace(placement.Proposal.OrientedProtein.TopologySha256) ||
            !Available(workingDirectory, placement.Proposal.OrientedProtein.TopologyPath) ||
            !File.Exists(packmolExecutablePath) ||
            string.IsNullOrWhiteSpace(packmolExecutableSha256) || packmolExecutableSha256.Length != 64 ||
            !packmolExecutableSha256.All(Uri.IsHexDigit) ||
            policy.ForceFieldFiles.Any(asset =>
                string.IsNullOrWhiteSpace(asset.Id) || string.IsNullOrWhiteSpace(asset.Version) ||
                string.IsNullOrWhiteSpace(asset.Family) || string.IsNullOrWhiteSpace(asset.Sha256) ||
                !Available(workingDirectory, asset.Path)) ||
            policy.ForceFieldFiles.GroupBy(asset => asset.Path, StringComparer.Ordinal)
                .Any(group => group.Select(asset => asset.Sha256).Distinct(StringComparer.OrdinalIgnoreCase).Count() != 1) ||
            policy.ForceFieldFiles.GroupBy(asset => asset.Sha256, StringComparer.OrdinalIgnoreCase)
                .Any(group => group.Select(asset => (asset.Family, asset.Version)).Distinct().Count() != 1))
            return Refusal("A qualified molecular representation, installed external tool, or exact coordinate/parameter asset is unavailable before start.");
        var forceFieldFiles = policy.ForceFieldFiles.DistinctBy(asset => asset.Sha256,
            StringComparer.OrdinalIgnoreCase).ToImmutableArray();
        var construction = policy.Construction;
        if (!ValidConstructionPolicy(construction) || !ValidLocalStatePolicy(policy) ||
            !ValidStageProteinGeometryPolicy(policy) ||
            placement.Proposal.MidplaneAngstrom is not double midplane ||
            placement.Proposal.ThicknessAngstrom is not double thickness ||
            !double.IsFinite(midplane) || !double.IsFinite(thickness) || thickness <= 0)
            return NoAttempt("The applicable construction geometry and parameter policy are incomplete.");

        onAccepted?.Invoke(attempt);
        var running = State(attempt.Id, StageExecutionStanding.Running, "The identified preparation attempt is deriving and constructing its explicit system.");
        progress?.Report(running);
        var half = thickness / 2.0;
        var regions = ImmutableArray.Create(
            new LeafletProjectionRegion(LeafletSide.Upper, midplane, midplane + half + construction.HeadRegionThicknessAngstrom),
            new LeafletProjectionRegion(LeafletSide.Lower, midplane - half - construction.HeadRegionThicknessAngstrom, midplane));
        try
        {
            var measureRequest = new ScientificWorkRequest<ConstructionInputPayload>(Guid.NewGuid().ToString("N"), workingDirectory,
                new ConstructionInputPayload(revision.Id, attempt.Id, protein.Id,
                    placement.Proposal.OrientedProtein.CoordinatePath,
                    placement.Proposal.OrientedProtein.TopologyPath!,
                    placement.Proposal.OrientedProtein.TopologySha256!,
                    forceFieldFiles, midplane, regions, construction.AtomRadiusByElementAngstrom,
                    construction.LateralClearanceAngstrom, construction.WaterMarginAngstrom,
                    construction.HeadRegionThicknessAngstrom, construction.ProjectionClearanceAngstrom,
                    construction.GridResolutionAngstrom, construction.VolumeGridResolutionAngstrom));
            var measured = await _worker.MeasureConstructionInputsAsync(measureRequest, cancellationToken);
            if (measured.RequestId != measureRequest.RequestId || measured.StudyRevisionId != revision.Id ||
                measured.AttemptId != attempt.Id ||
                measured.Standing != WorkerResultStanding.Observed || measured.Observations is null)
                return Failed(attempt, measured.Standing, measured.FailureMessage ?? "Construction inputs were not observed.");
            if (!TryDerive(attempt, revision, membrane, policy, measured.Observations, midplane, thickness,
                    protein.Molecule.AtomCount, out var derivation, out var reason))
                return new PreparationStartResult(attempt, null, null, State(attempt.Id, StageExecutionStanding.Failed, reason));
            if (derivation!.CellAngstrom.Any(length =>
                    policy.SystemSettings.NonbondedCutoffNanometers * 20 >= length))
                return new PreparationStartResult(attempt, derivation, null,
                    State(attempt.Id, StageExecutionStanding.Failed,
                        "The declared nonbonded cutoff is not smaller than half of every derived periodic cell dimension."));

            progress?.Report(State(attempt.Id, StageExecutionStanding.Running, "The exact intended molecules are being packed and parameterized."));
            var components = MakeComponents(derivation!, membrane, policy, midplane, thickness);
            var constructRequest = new ScientificWorkRequest<ConstructionPayload>(Guid.NewGuid().ToString("N"), workingDirectory,
                new ConstructionPayload(revision.Id, attempt.Id, placement.Proposal.OrientedProtein.CoordinatePath,
                    protein.Molecule.CoordinatePath, protein.Molecule.CoordinateSha256, protein.Correspondence,
                    placement.Proposal.OrientedProtein.TopologyPath!,
                    placement.Proposal.OrientedProtein.TopologySha256!,
                    derivation!.CellOriginAngstrom, derivation.CellAngstrom, components,
                    packmolExecutablePath, packmolVersion, packmolExecutableSha256, forceFieldFiles,
                    policy.SystemSettings,
                    policy.PackingToleranceAngstrom, construction.MaximumPackingAttempts,
                    policy.LocalStateObservation));
            var built = await _worker.ConstructSystemAsync(constructRequest, cancellationToken);
            if (built.RequestId != constructRequest.RequestId || built.StudyRevisionId != revision.Id ||
                built.AttemptId != attempt.Id || built.Standing != WorkerResultStanding.Observed || built.Observations is null)
                return Failed(attempt, built.Standing, built.FailureMessage ?? "The packed candidate was not observed.", derivation);
            var candidate = await AssessConstructedAsync(attempt, revision, protein, membrane, placement, policy,
                derivation, components, built, cancellationToken);
            if (candidate is null)
                return new PreparationStartResult(attempt, derivation, null,
                    State(attempt.Id, StageExecutionStanding.Failed,
                        "Actual membership, whole-system parameters, geometry, and correspondence were not established."));
            return new PreparationStartResult(attempt, derivation, candidate,
                State(attempt.Id, StageExecutionStanding.Completed, "An eligible explicit system was constructed; minimization is still required."));
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new PreparationStartResult(attempt, null, null,
                State(attempt.Id, StageExecutionStanding.Stopped, "The preparation attempt was stopped without a constructed system."));
        }
    }

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

    private static bool ValidConstructionPolicy(ConstructionPolicy p) =>
        p.MaximumAtomCount > 0 && p.MaximumPackingAttempts > 0 &&
        p.MaximumCellDimensionAngstrom > 0 && p.LateralClearanceAngstrom > 0 &&
        p.WaterMarginAngstrom > 0 && p.HeadRegionThicknessAngstrom > 0 &&
        p.WaterNumberDensityPerAngstromCubed > 0 && p.GridResolutionAngstrom > 0 &&
        p.VolumeGridResolutionAngstrom > 0 && p.ProjectionClearanceAngstrom >= 0 &&
        !p.AtomRadiusByElementAngstrom.IsEmpty &&
        p.AtomRadiusByElementAngstrom.All(pair => pair.Value > 0 && double.IsFinite(pair.Value));

    private static bool TryDerive(PreparationAttempt attempt, StudyRevision revision,
        AssessedMembraneModel membrane, ApplicablePreparationPolicy policy,
        ConstructionInputObservations measured, double midplane, double membraneThickness,
        int proteinAtomCount, out ConstructionDerivation? derivation, out string reason)
    {
        derivation = null;
        reason = "The policy-governed finite construction could not be derived.";
        var origin = measured.CandidateCellOriginAngstrom;
        var cell = measured.CandidateCellAngstrom;
        if (origin.Length != 3 || cell.Length != 3 ||
            origin.Any(value => !double.IsFinite(value)) ||
            cell.Any(value => !double.IsFinite(value) || value <= 0 ||
                value > policy.Construction.MaximumCellDimensionAngstrom) ||
            !double.IsFinite(measured.UpperOccludedAreaAngstromSquared) ||
            !double.IsFinite(measured.LowerOccludedAreaAngstromSquared) ||
            !double.IsFinite(measured.ProteinAqueousOccludedVolumeAngstromCubed) ||
            !double.IsFinite(measured.NetChargeElementary))
        {
            reason = "Measured coordinate extents, accessible regions, charge, or policy cell bound are invalid.";
            return false;
        }
        var xyArea = cell[0] * cell[1];
        if (measured.UpperOccludedAreaAngstromSquared < 0 || measured.LowerOccludedAreaAngstromSquared < 0 ||
            measured.UpperOccludedAreaAngstromSquared >= xyArea || measured.LowerOccludedAreaAngstromSquared >= xyArea)
        {
            reason = "The projected protein occupancy leaves no feasible area for both leaflets.";
            return false;
        }
        var representations = membrane.SpeciesRepresentations.ToDictionary(item => item.SpeciesId, StringComparer.Ordinal);
        var counts = ImmutableArray.CreateBuilder<SpeciesCount>();
        if (!TryApportion(membrane.Intended.Upper, xyArea - measured.UpperOccludedAreaAngstromSquared,
                representations, counts, out reason) ||
            !TryApportion(membrane.Intended.Lower, xyArea - measured.LowerOccludedAreaAngstromSquared,
                representations, counts, out reason))
            return false;

        var membraneBottom = midplane - membraneThickness / 2.0 - policy.Construction.HeadRegionThicknessAngstrom;
        var membraneTop = midplane + membraneThickness / 2.0 + policy.Construction.HeadRegionThicknessAngstrom;
        var cellBottom = origin[2];
        var cellTop = origin[2] + cell[2];
        var aqueousHeight = Math.Max(0, membraneBottom - cellBottom) + Math.Max(0, cellTop - membraneTop);
        if (membraneBottom - cellBottom < policy.Construction.WaterMarginAngstrom ||
            cellTop - membraneTop < policy.Construction.WaterMarginAngstrom)
        {
            reason = "The derived cell lacks the declared hydration margin on both sides of the bilayer.";
            return false;
        }
        var availableAqueousVolume = xyArea * aqueousHeight - measured.ProteinAqueousOccludedVolumeAngstromCubed;
        if (!double.IsFinite(availableAqueousVolume) || availableAqueousVolume <= 0)
        {
            reason = "The measured cell and membrane region leave no positive aqueous volume.";
            return false;
        }
        var solventPopulation = availableAqueousVolume * policy.Construction.WaterNumberDensityPerAngstromCubed;
        var intendedSaltPairs = revision.Conditions.TargetNaClMolar * MoleculesPerMolarAngstromCubed * availableAqueousVolume;
        if (!double.IsFinite(solventPopulation) || !double.IsFinite(intendedSaltPairs) ||
            solventPopulation > int.MaxValue || intendedSaltPairs > int.MaxValue || intendedSaltPairs < 0)
        {
            reason = "The finite solvent or salt population cannot be represented.";
            return false;
        }
        var lipidCharge = counts.Sum(item => representations[item.SpeciesId].NetChargeElementary * item.Count);
        var totalUnneutralizedCharge = measured.NetChargeElementary + lipidCharge;
        var integralCharge = Math.Round(totalUnneutralizedCharge, MidpointRounding.AwayFromZero);
        if (!double.IsFinite(totalUnneutralizedCharge) || Math.Abs(totalUnneutralizedCharge - integralCharge) > 1e-5 ||
            Math.Abs(integralCharge) > int.MaxValue / 2)
        {
            reason = "The exact protein–lipid charge cannot be neutralized by the supported monovalent-ion rule.";
            return false;
        }
        var saltPairs = checked((int)Math.Round(intendedSaltPairs, MidpointRounding.AwayFromZero));
        var sodium = saltPairs + (integralCharge < 0 ? checked((int)-integralCharge) : 0);
        var chloride = saltPairs + (integralCharge > 0 ? checked((int)integralCharge) : 0);
        var water = (int)Math.Round(solventPopulation, MidpointRounding.AwayFromZero) - sodium - chloride;
        if (water < 2)
        {
            reason = "The finite water region cannot accommodate the required neutralization and salt pairs.";
            return false;
        }
        var estimatedAtoms = proteinAtomCount > 0 && measured.ProteinXMaxAngstrom > measured.ProteinXMinAngstrom ?
            proteinAtomCount + counts.Sum(item => (long)item.Count * representations[item.SpeciesId].AtomCount) +
            (long)water * policy.Water.AtomCount + (long)sodium * policy.Sodium.AtomCount +
            (long)chloride * policy.Chloride.AtomCount : long.MaxValue;
        if (estimatedAtoms > policy.Construction.MaximumAtomCount)
        {
            reason = "The policy's evidence-backed resource bound is exceeded before packing.";
            return false;
        }
        var approximations = ImmutableArray.CreateBuilder<string>();
        approximations.Add(policy.Construction.ApproximationStatement);
        approximations.AddRange(measured.ApproximationWarnings);
        approximations.Add("Finite lipid and ion populations approximate intended fractions and 0.15 M NaCl; they are not achieved-equilibration measurements.");
        derivation = new ConstructionDerivation(attempt.Id, counts.ToImmutable(), origin, cell,
            water, sodium, chloride, measured.NetChargeElementary, revision.Conditions.TargetNaClMolar,
            saltPairs / (MoleculesPerMolarAngstromCubed * availableAqueousVolume),
            availableAqueousVolume, approximations.ToImmutable(), policy.Limitations);
        return true;
    }

    private static bool TryApportion(LeafletComposition leaflet, double availableArea,
        IReadOnlyDictionary<string, MolecularRepresentation> representations,
        ImmutableArray<SpeciesCount>.Builder output, out string reason)
    {
        reason = "The leaflet cannot be apportioned under the identified species footprints.";
        var fractions = leaflet.Fractions.Where(item => item.Fraction > 0).ToArray();
        if (fractions.Length == 0 || fractions.Any(item => !double.IsFinite(item.Fraction) ||
                !representations.TryGetValue(item.SpeciesId, out var representation) ||
                !double.IsFinite(representation.AreaPerMoleculeAngstromSquared) ||
                representation.AreaPerMoleculeAngstromSquared <= 0) ||
            Math.Abs(fractions.Sum(item => item.Fraction) - 1.0) > 1e-8)
            return false;
        var weightedFootprint = fractions.Sum(item => item.Fraction * representations[item.SpeciesId].AreaPerMoleculeAngstromSquared);
        var rawTotal = availableArea / weightedFootprint;
        if (!double.IsFinite(rawTotal) || rawTotal < fractions.Length || rawTotal > int.MaxValue)
            return false;
        var total = (int)Math.Round(rawTotal, MidpointRounding.AwayFromZero);
        var targets = fractions.Select(item => item.Fraction * total).ToArray();
        var allotted = targets.Select(target => Math.Max(1, (int)Math.Floor(target))).ToArray();
        if (allotted.Sum() > total)
            return false;
        var remaining = total - allotted.Sum();
        foreach (var index in Enumerable.Range(0, targets.Length)
                     .OrderByDescending(index => targets[index] - Math.Floor(targets[index])))
        {
            if (remaining == 0) break;
            var upper = Math.Max(1, (int)Math.Ceiling(targets[index]));
            if (allotted[index] < upper)
            {
                allotted[index]++;
                remaining--;
            }
        }
        if (remaining != 0 || allotted.Where((count, index) =>
                Math.Abs(count - targets[index]) >= 1.0).Any())
            return false;
        for (var index = 0; index < fractions.Length; index++)
            output.Add(new SpeciesCount(leaflet.PhysicalSide, fractions[index].SpeciesId,
                allotted[index], fractions[index].Fraction));
        reason = string.Empty;
        return true;
    }

    private static ImmutableArray<ConstructionComponent> MakeComponents(
        ConstructionDerivation derivation, AssessedMembraneModel membrane,
        ApplicablePreparationPolicy policy, double midplane, double thickness)
    {
        var origin = derivation.CellOriginAngstrom;
        var cell = derivation.CellAngstrom;
        var x0 = origin[0]; var x1 = origin[0] + cell[0];
        var y0 = origin[1]; var y1 = origin[1] + cell[1];
        var z0 = origin[2]; var z1 = origin[2] + cell[2];
        var membraneTop = midplane + thickness / 2.0 + policy.Construction.HeadRegionThicknessAngstrom;
        var membraneBottom = midplane - thickness / 2.0 - policy.Construction.HeadRegionThicknessAngstrom;
        var upper = new SpatialRegionAngstrom(x0, y0, midplane, x1, y1, membraneTop);
        var lower = new SpatialRegionAngstrom(x0, y0, membraneBottom, x1, y1, midplane);
        var upperHead = new SpatialRegionAngstrom(x0, y0, midplane + thickness / 2.0, x1, y1, membraneTop);
        var lowerHead = new SpatialRegionAngstrom(x0, y0, membraneBottom, x1, y1, midplane - thickness / 2.0);
        var waterUpper = new SpatialRegionAngstrom(x0, y0, membraneTop, x1, y1, z1);
        var waterLower = new SpatialRegionAngstrom(x0, y0, z0, x1, y1, membraneBottom);
        var species = membrane.SpeciesRepresentations.ToDictionary(item => item.SpeciesId, StringComparer.Ordinal);
        var components = ImmutableArray.CreateBuilder<ConstructionComponent>();
        foreach (var lipid in derivation.LipidCounts)
        {
            var representation = species[lipid.SpeciesId];
            var isUpper = lipid.PhysicalSide == LeafletSide.Upper;
            components.Add(new ConstructionComponent(GeneratedComponentRoleKind.Lipid, lipid.PhysicalSide,
                lipid.SpeciesId, representation.CoordinateTemplatePath, representation.CoordinateTemplateSha256,
                lipid.Count,
                isUpper ? upper : lower, isUpper ? upperHead : lowerHead,
                representation.HeadAtomIndices));
        }
        AddSplit(components, GeneratedComponentRoleKind.Water, policy.Water, derivation.WaterCount, waterUpper, waterLower);
        AddSplit(components, GeneratedComponentRoleKind.PositiveIon, policy.Sodium, derivation.SodiumCount, waterUpper, waterLower);
        AddSplit(components, GeneratedComponentRoleKind.NegativeIon, policy.Chloride, derivation.ChlorideCount, waterUpper, waterLower);
        return components.ToImmutable();
    }

    private static void AddSplit(ImmutableArray<ConstructionComponent>.Builder output, GeneratedComponentRoleKind kind,
        MolecularRepresentation species, int count, SpatialRegionAngstrom upper, SpatialRegionAngstrom lower)
    {
        var upperCount = count / 2 + count % 2;
        var lowerCount = count / 2;
        if (upperCount > 0)
            output.Add(new ConstructionComponent(kind, LeafletSide.Upper, species.SpeciesId,
                species.CoordinateTemplatePath, species.CoordinateTemplateSha256,
                upperCount, upper, null, ImmutableArray<int>.Empty));
        if (lowerCount > 0)
            output.Add(new ConstructionComponent(kind, LeafletSide.Lower, species.SpeciesId,
                species.CoordinateTemplatePath, species.CoordinateTemplateSha256,
                lowerCount, lower, null, ImmutableArray<int>.Empty));
    }

    private static async Task<ConstructedExplicitSystem?> AssessConstructedAsync(
        PreparationAttempt attempt, StudyRevision revision, AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy,
        ConstructionDerivation derivation, ImmutableArray<ConstructionComponent> components,
        WorkerResult<ConstructionObservations> result,
        CancellationToken cancellationToken)
    {
        var observed = result.Observations!;
        if (observed.AtomCount <= protein.Molecule.AtomCount ||
            observed.CorrespondedResultAtomCount != observed.AtomCount ||
            observed.WaterCount != derivation.WaterCount ||
            observed.PositiveIonCount != derivation.SodiumCount ||
            observed.NegativeIonCount != derivation.ChlorideCount ||
            observed.ActualCellAngstrom.Length != 3 ||
            observed.ActualCellAngstrom.Where((value, index) =>
                !double.IsFinite(value) || Math.Abs(value - derivation.CellAngstrom[index]) > 0.01).Any() ||
            !double.IsFinite(observed.NetChargeElementary) || Math.Abs(observed.NetChargeElementary) > 1e-5 ||
            !double.IsFinite(observed.InitialPotentialEnergyKjMol) ||
            !observed.ContactWarnings.IsDefaultOrEmpty || !observed.GeometryWarnings.IsDefaultOrEmpty ||
            !observed.ParameterWarnings.IsDefaultOrEmpty)
            return null;
        if (!LocalStateSupportsConstruction(observed.LocalState, policy))
            return null;
        if (components.IsDefaultOrEmpty || observed.SpeciesCounts.IsDefault ||
            components.Any(item => item.Count <= 0 || item.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                item.Role is not (GeneratedComponentRoleKind.Lipid or GeneratedComponentRoleKind.Water or GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon)) ||
            observed.SpeciesCounts.Any(item => item.Count <= 0 || item.PhysicalSide is not (LeafletSide.Upper or LeafletSide.Lower) ||
                item.Role is not (GeneratedComponentRoleKind.Lipid or GeneratedComponentRoleKind.Water or GeneratedComponentRoleKind.PositiveIon or GeneratedComponentRoleKind.NegativeIon)))
            return null;
        var expectedMolecules = components.GroupBy(item => (item.Role, item.PhysicalSide, item.SpeciesId))
            .ToDictionary(group => group.Key, group => group.Sum(item => item.Count));
        if (observed.SpeciesCounts.Length != expectedMolecules.Count ||
            observed.SpeciesCounts.Select(item => (item.Role, item.PhysicalSide, item.SpeciesId))
                .Distinct().Count() != observed.SpeciesCounts.Length ||
            observed.SpeciesCounts.Any(item =>
                !expectedMolecules.TryGetValue((item.Role, item.PhysicalSide, item.SpeciesId), out var count) ||
                item.Count != count))
            return null;
        var representations = membrane.SpeciesRepresentations.Append(policy.Water).Append(policy.Sodium)
            .Append(policy.Chloride).ToArray();
        if (representations.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() !=
            representations.Length)
            return null;
        long expectedAtomCount = protein.Molecule.AtomCount;
        var expectedGeneratedAtoms = new Dictionary<(GeneratedComponentRoleKind Role, LeafletSide Side, string SpeciesId), long>();
        foreach (var component in components)
        {
            var representation = representations.FirstOrDefault(item => item.SpeciesId == component.SpeciesId);
            if (representation is null || representation.AtomCount <= 0 ||
                component.TemplateCoordinateSha256 != representation.CoordinateTemplateSha256)
                return null;
            var atoms = (long)component.Count * representation.AtomCount;
            expectedAtomCount += atoms;
            var key = (component.Role, component.PhysicalSide, component.SpeciesId);
            expectedGeneratedAtoms[key] = expectedGeneratedAtoms.GetValueOrDefault(key) + atoms;
        }
        if (expectedAtomCount != observed.AtomCount ||
            components.Where(item => item.Role == GeneratedComponentRoleKind.Water).Sum(item => item.Count) != observed.WaterCount ||
            components.Where(item => item.Role == GeneratedComponentRoleKind.PositiveIon).Sum(item => item.Count) != observed.PositiveIonCount ||
            components.Where(item => item.Role == GeneratedComponentRoleKind.NegativeIon).Sum(item => item.Count) != observed.NegativeIonCount ||
            !derivation.LipidCounts.All(item =>
                expectedMolecules.TryGetValue((GeneratedComponentRoleKind.Lipid, item.PhysicalSide, item.SpeciesId), out var count) &&
                count == item.Count))
            return null;

        var topologyCif = result.Artifacts.FirstOrDefault(item => item.Role == "topologyCif");
        var topologyJson = result.Artifacts.FirstOrDefault(item => item.Role == "topologyJson");
        var system = result.Artifacts.FirstOrDefault(item => item.Role == "systemXml");
        var state = result.Artifacts.FirstOrDefault(item => item.Role == "stateXml");
        var mapping = result.Artifacts.FirstOrDefault(item => item.Role == "correspondenceJson");
        if (topologyCif is null || topologyJson is null || system is null || state is null || mapping is null ||
            !File.Exists(topologyCif.Path) || !File.Exists(topologyJson.Path) ||
            !File.Exists(system.Path) || !File.Exists(state.Path) || !File.Exists(mapping.Path))
            return null;
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
            $"Attempt {attempt.Id}; policy {policy.Id}; placement {placement.Id}",
            "Eligible for required minimization only; not yet a qualified completed preparation.", EvidenceBearing.Supports));
        var conditions = $"Fixed nominal pH {revision.Conditions.NominalPh:G6}; intended NaCl {derivation.IntendedNaClMolar:G6} M; " +
                         $"finite-cell estimated NaCl {derivation.EstimatedNaClMolar:G6} M; " +
                         "thermal equilibration has not been established.";
        return new ConstructedExplicitSystem(molecule.Id, attempt, molecule, derivation,
            correspondence, derivation.LipidCounts, evidence, ImmutableArray<ScientificFinding>.Empty, conditions,
            observed.LocalState);
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
                "minimumIntermolecularDistanceAngstrom" or "upperLipidHeadMeanZAngstrom" or
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
                item.MeasurementName == "minimumIntermolecularDistanceAngstrom" &&
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
