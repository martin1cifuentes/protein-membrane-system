using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation;

/// <summary>Runs only a separately requested, already applicable optional procedure.</summary>
public sealed class OptionalEquilibrationProcedure
{
    private readonly IOptionalEquilibrationWork _worker;

    public OptionalEquilibrationProcedure(IOptionalEquilibrationWork worker) => _worker = worker;

    public async Task<StageOperationResult> RunAsync(
        CompletedStage minimized,
        ApplicablePreparationPolicy policy,
        string workingDirectory,
        IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
    {
        var stageId = Guid.NewGuid().ToString("N");
        var protocol = policy.OptionalEquilibration;
        var resolvedObservables = protocol is null
            ? ImmutableArray<ResolvedEquilibrationObservable>.Empty
            : ResolveObservables(protocol, minimized.Correspondence);
        if (minimized.Kind != StageKind.Minimization || minimized.PolicyId != policy.Id ||
            minimized.Molecule.TopologyPath is null || minimized.Molecule.SystemXmlPath is null ||
            minimized.Molecule.StateXmlPath is null || minimized.Molecule.TopologySha256 is null ||
            minimized.Molecule.SystemXmlSha256 is null || minimized.Molecule.StateXmlSha256 is null ||
            minimized.Molecule.CorrespondencePath is null ||
            minimized.Molecule.CorrespondenceSha256 is null ||
            policy.LocalStateObservation is null || protocol is null ||
            !ValidProtocol(protocol, resolvedObservables) || !minimized.Correspondence.Complete)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "The source minimized stage or exact applicable optional procedure is unavailable."),
                ImmutableArray<ScientificFinding>.Empty);

        var proteinBackbone = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Protein && atom.AtomRole == AtomRoleKind.Backbone)
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        var lipidHeavy = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Lipid && atom.Element != "H")
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        if ((protocol.Stages.Any(stage => stage.ProteinRestraintKjMolNm2 > 0) && proteinBackbone.IsDefaultOrEmpty) ||
            (protocol.Stages.Any(stage => stage.LipidRestraintKjMolNm2 > 0) && lipidHeavy.IsDefaultOrEmpty) ||
            proteinBackbone.Concat(lipidHeavy).Any(index => index < 0 || index >= minimized.Molecule.AtomCount))
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "Declared restraint targets cannot be mapped to the exact minimized molecular state."),
                ImmutableArray<ScientificFinding>.Empty);

        progress?.Report(State(minimized.Attempt.Id, stageId, StageExecutionStanding.Running,
            "Optional equilibration is running."));
        var request = new ScientificWorkRequest<EquilibrationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new EquilibrationPayload(minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId,
                minimized.Id, minimized.Molecule.CoordinatePath, minimized.Molecule.CoordinateSha256,
                minimized.Molecule.TopologyPath, minimized.Molecule.TopologySha256,
                minimized.Molecule.SystemXmlPath, minimized.Molecule.SystemXmlSha256,
                minimized.Molecule.StateXmlPath, minimized.Molecule.StateXmlSha256,
                policy.Id, proteinBackbone, lipidHeavy,
                resolvedObservables, protocol));
        WorkerResult<EquilibrationObservations> result;
        try
        {
            result = await _worker.EquilibrateAsync(request, cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Optional equilibration stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        }
        if (result.RequestId != request.RequestId || result.StudyRevisionId != minimized.Attempt.StudyRevisionId ||
            result.AttemptId != minimized.Attempt.Id || result.StageId != stageId ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, result.Standing switch
                {
                    WorkerResultStanding.Stopped => StageExecutionStanding.Stopped,
                    WorkerResultStanding.Unobserved => StageExecutionStanding.Unobserved,
                    _ => StageExecutionStanding.Failed
                }, result.FailureMessage ?? "Optional equilibration did not establish an observed completed stage."),
                ImmutableArray<ScientificFinding>.Empty);

        var observed = result.Observations;
        var equilibratedState = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedStateXml");
        var equilibratedCoordinates = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedCif");
        var equilibratedSystem = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedSystemXml");
        if (equilibratedState is null || equilibratedCoordinates is null || equilibratedSystem is null ||
            observed.FinalAtomCount != minimized.Molecule.AtomCount ||
            observed.ExtensionsPerformed < 0 || observed.ExtensionsPerformed > protocol.MaximumExtensions ||
            observed.Windows.Length != protocol.Stages.Length + observed.ExtensionsPerformed ||
            observed.Windows.Where((window, index) =>
                window.Name != (index < protocol.Stages.Length ? protocol.Stages[index].Name :
                    $"{protocol.ExtensionWindow.Name}-{index - protocol.Stages.Length + 1}") ||
                window.RequestedSteps != (index < protocol.Stages.Length ? protocol.Stages[index].Steps : protocol.ExtensionWindow.Steps) ||
                window.CompletedSteps != window.RequestedSteps ||
                window.Measurements.IsDefaultOrEmpty ||
                window.Measurements.Any(value => !double.IsFinite(value.Value))).Any() ||
            !observed.FinalObservationUnrestrained ||
            !observed.NumericalWarnings.IsDefaultOrEmpty ||
            !ValidAdequacyResult(protocol, observed) ||
            observed.Termination != StageTermination.Completed)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "The complete declared optional procedure and trustworthy final state were not observed."),
                ImmutableArray<ScientificFinding>.Empty);

        var observationRequest = new ScientificWorkRequest<StageObservationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new StageObservationPayload(minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId,
                equilibratedCoordinates.Path, equilibratedCoordinates.Sha256,
                minimized.Molecule.TopologyPath, minimized.Molecule.TopologySha256,
                equilibratedSystem.Path, equilibratedSystem.Sha256,
                equilibratedState.Path, equilibratedState.Sha256,
                StageKind.Equilibration, minimized.Molecule.CorrespondencePath!,
                minimized.Molecule.CorrespondenceSha256!, policy.LocalStateObservation,
                policy.StageProteinGeometryMeasurement));
        var stageResult = await _worker.ObserveStageAsync(observationRequest, cancellationToken);
        if (stageResult.RequestId != observationRequest.RequestId ||
            stageResult.StudyRevisionId != minimized.Attempt.StudyRevisionId ||
            stageResult.AttemptId != minimized.Attempt.Id || stageResult.StageId != stageId ||
            stageResult.Standing != WorkerResultStanding.Observed || stageResult.Observations is null ||
            stageResult.Observations.AtomCount != minimized.Molecule.AtomCount ||
            !stageResult.Observations.AtomOrderMatched || !stageResult.Observations.BondsMatched ||
            !stageResult.Observations.NumericalWarnings.IsDefaultOrEmpty ||
            protocol.RequiredObservations.Any(name =>
                !stageResult.Observations.Measurements.Concat(
                    observed.Windows.SelectMany(window => window.Measurements)).Concat(
                    observed.Samples.SelectMany(sample => sample.Measurements))
                    .Any(value => value.Name == name && double.IsFinite(value.Value))))
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Unobserved,
                    "Optional procedure output exists, but attributable final-stage observation is unavailable."),
                ImmutableArray<ScientificFinding>.Empty);

        var molecule = minimized.Molecule with
        {
            Id = stageId,
            CoordinatePath = equilibratedCoordinates.Path,
            CoordinateSha256 = equilibratedCoordinates.Sha256,
            SystemXmlPath = equilibratedSystem.Path,
            SystemXmlSha256 = equilibratedSystem.Sha256,
            StateXmlPath = equilibratedState.Path,
            StateXmlSha256 = equilibratedState.Sha256
        };
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), stageId, result.Provider?.Name ?? "OpenMM",
            "Observed declared optional equilibration",
            $"{observed.Windows.Length} declared windows completed with an unrestrained final observation",
            $"Attempt {minimized.Attempt.Id}; source minimized stage {minimized.Id}; policy {policy.Id}",
            "Procedure completion does not prove global thermodynamic equilibrium or scientific qualification.",
            EvidenceBearing.Context));
        // The final sampled, policy-named observables remain available to the
        // preparation assessment's stage-specific absolute-value criteria;
        // trace sufficiency alone is not a suitable-state conclusion.
        var measurements = stageResult.Observations.Measurements.AddRange(observed.Samples[^1].Measurements)
            .AddRange(observed.Windows
            .Where(window => window.TemperatureKelvin.HasValue)
            .Select(window => new MeasuredValue("temperature", window.TemperatureKelvin!.Value,
                "K", window.Name)));
        var local = stageResult.Observations.LocalState;
        if (local is not null && !local.Measurements.IsDefault)
            measurements = measurements.AddRange(local.Measurements);
        var observation = new StageObservation(stageId, minimized.Attempt.Id, StageKind.Equilibration,
            measurements, evidence, observed.Termination, result.Provider?.Version ?? "unknown", DateTimeOffset.UtcNow,
            observed.ObservationAdequacy, observed.ObservationAssessments, observed.Samples,
            stageResult.Observations.LocalState, stageResult.Observations.ProteinGeometry);
        var findings = stageResult.Observations.ContactWarnings.Concat(stageResult.Observations.StructuralWarnings)
            .Concat(observed.Windows.SelectMany(window => window.Warnings))
            .Select(warning => new ScientificFinding(
                Guid.NewGuid().ToString("N"), stageId, evidence[0].Id, warning,
                "Requires stage-specific preparation assessment.", FindingDisposition.Challenges,
                true, DateTimeOffset.UtcNow)).ToImmutableArray();
        var completed = new CompletedStage(stageId, minimized.Attempt, StageKind.Equilibration,
            molecule, observation, minimized.Correspondence with { ResultId = stageId }, policy.Id,
            minimized.Id, findings, DateTimeOffset.UtcNow);
        var state = State(minimized.Attempt.Id, stageId, StageExecutionStanding.Completed,
            "Optional equilibration completed; its scientific assessment remains separate.", 1.0);
        progress?.Report(state);
        return new StageOperationResult(completed, state, findings);
    }

    private static bool ValidProtocol(EquilibrationProtocol protocol,
        ImmutableArray<ResolvedEquilibrationObservable> resolved)
    {
        if (string.IsNullOrWhiteSpace(protocol.Id) || protocol.RandomSeed <= 0 ||
            protocol.Stages.IsDefaultOrEmpty || protocol.MaximumExtensions < 0 ||
            protocol.MaximumSampleCount <= 0 || protocol.ExtensionWindow is null ||
            protocol.Observables.IsDefaultOrEmpty || resolved.IsDefaultOrEmpty ||
            protocol.SufficiencyRules.IsDefaultOrEmpty ||
            !double.IsFinite(protocol.TargetTemperatureKelvin) || protocol.TargetTemperatureKelvin != 303.0 ||
            protocol.RequiredObservations.IsDefaultOrEmpty || string.IsNullOrWhiteSpace(protocol.ComparisonBasis))
            return false;
        if (protocol.Stages.Append(protocol.ExtensionWindow).Any(stage => stage.Steps <= 0 || stage.ReportIntervalSteps <= 0 ||
            stage.TimestepPicoseconds <= 0 || !double.IsFinite(stage.TimestepPicoseconds) ||
            !double.IsFinite(stage.TemperatureKelvin) || stage.TemperatureKelvin <= 0 ||
            !double.IsFinite(stage.FrictionPerPicosecond) || stage.FrictionPerPicosecond <= 0 ||
            !double.IsFinite(stage.ProteinRestraintKjMolNm2) || stage.ProteinRestraintKjMolNm2 < 0 ||
            !double.IsFinite(stage.LipidRestraintKjMolNm2) || stage.LipidRestraintKjMolNm2 < 0 ||
            stage.PressureBar.HasValue && (!double.IsFinite(stage.PressureBar.Value) ||
                                           stage.BarostatFrequencySteps is null or <= 0)))
            return false;
        var last = protocol.Stages[^1];
        var extension = protocol.ExtensionWindow;
        if (last.ProteinRestraintKjMolNm2 != 0 || last.LipidRestraintKjMolNm2 != 0 ||
            extension.ProteinRestraintKjMolNm2 != 0 || extension.LipidRestraintKjMolNm2 != 0 ||
            last.TemperatureKelvin != protocol.TargetTemperatureKelvin ||
            extension.TemperatureKelvin != last.TemperatureKelvin ||
            extension.PressureBar != last.PressureBar || extension.PressureMode != last.PressureMode ||
            extension.BarostatFrequencySteps != last.BarostatFrequencySteps ||
            extension.SurfaceTensionBarNm != last.SurfaceTensionBarNm ||
            extension.FrictionPerPicosecond != last.FrictionPerPicosecond ||
            extension.TimestepPicoseconds != last.TimestepPicoseconds ||
            protocol.Stages.Any(stage => string.IsNullOrWhiteSpace(stage.Name)) ||
            string.IsNullOrWhiteSpace(extension.Name) ||
            protocol.Stages.Select(stage => stage.Name).Distinct(StringComparer.Ordinal).Count() != protocol.Stages.Length ||
            protocol.Stages.Any(stage => stage.Name == extension.Name ||
                stage.Name.StartsWith(extension.Name + "-", StringComparison.Ordinal)) ||
            protocol.Observables.Select(item => item.Name).Distinct(StringComparer.Ordinal).Count() != protocol.Observables.Length ||
            protocol.SufficiencyRules.Select(item => item.ObservableName).Distinct(StringComparer.Ordinal).Count() != protocol.SufficiencyRules.Length ||
            protocol.SufficiencyRules.Length != protocol.Observables.Length ||
            protocol.RequiredObservations.Any(name => !protocol.Observables.Any(item => item.Name == name)))
            return false;
        var declaredSamples = protocol.Stages.Sum(stage =>
            (long)Math.Ceiling((double)stage.Steps / stage.ReportIntervalSteps)) +
            (long)protocol.MaximumExtensions *
            (long)Math.Ceiling((double)extension.Steps / extension.ReportIntervalSteps);
        if (declaredSamples > protocol.MaximumSampleCount) return false;
        var categories = protocol.Observables.Select(item => item.Category).ToHashSet(StringComparer.Ordinal);
        if (!new[] { "system", "membraneOrganization", "proteinPlacement", "hydrationIons" }
                .All(categories.Contains)) return false;
        if (resolved.Length != protocol.Observables.Length) return false;
        foreach (var pair in protocol.Observables.Zip(resolved))
        {
            var observable = pair.First;
            var actual = pair.Second;
            if (string.IsNullOrWhiteSpace(observable.Name) || string.IsNullOrWhiteSpace(observable.Unit) ||
                string.IsNullOrWhiteSpace(observable.Scope) || observable.Name != actual.Name ||
                observable.Unit != actual.Unit || observable.Scope != actual.Scope ||
                observable.Category != actual.Category || observable.Method != actual.Method ||
                observable.DistanceCutoffAngstrom != actual.DistanceCutoffAngstrom ||
                (observable.ThirdSelector == "none" ? !actual.ThirdAtomIndices.IsEmpty :
                    actual.ThirdAtomIndices.IsDefaultOrEmpty))
                return false;
            var valid = observable.Category switch
            {
                "system" => (observable.Method is "potentialEnergy" or "kineticEnergy" or "temperature" or
                    "cellAreaXY" or "cellHeight") && observable.AtomSelector == "none" &&
                    observable.ComparisonSelector == "none" && observable.ThirdSelector == "none",
                "membraneOrganization" => observable.Method == "centroidSeparationZ" &&
                    observable.AtomSelector == "upper-lipid-head" &&
                    observable.ComparisonSelector == "lower-lipid-head" &&
                    observable.ThirdSelector == "none",
                "proteinPlacement" => observable.Method == "proteinMidplaneOffset" &&
                    observable.AtomSelector == "protein-backbone" &&
                    observable.ComparisonSelector == "upper-lipid-head" &&
                    observable.ThirdSelector == "lower-lipid-head",
                "hydrationIons" => observable.Method == "countWithinDistance" &&
                    (observable.AtomSelector is "water-oxygen" or "all-ions") &&
                    observable.ComparisonSelector == "protein-or-lipid-heavy" &&
                    observable.ThirdSelector == "none" &&
                    observable.DistanceCutoffAngstrom is double cutoff && double.IsFinite(cutoff) && cutoff > 0,
                _ => false
            };
            if (!valid) return false;
        }
        if (!protocol.Observables.Where(item => item.Category == "hydrationIons")
                .Any(item => item.AtomSelector == "water-oxygen") ||
            !protocol.Observables.Where(item => item.Category == "hydrationIons")
                .Any(item => item.AtomSelector == "all-ions")) return false;
        return protocol.SufficiencyRules.All(rule =>
            protocol.Observables.Any(item => item.Name == rule.ObservableName) &&
            rule.BlockSizeSamples > 0 && rule.MinimumEffectiveBlocks >= 3 &&
            double.IsFinite(rule.MaximumAbsoluteFirstVsLastBlockMeanDifference) &&
            rule.MaximumAbsoluteFirstVsLastBlockMeanDifference >= 0 &&
            double.IsFinite(rule.MaximumAbsoluteLagOneBlockCorrelation) &&
            rule.MaximumAbsoluteLagOneBlockCorrelation >= 0 &&
            rule.MaximumAbsoluteLagOneBlockCorrelation < 1);
    }

    private static ImmutableArray<ResolvedEquilibrationObservable> ResolveObservables(
        EquilibrationProtocol protocol, SourceToResultCorrespondence correspondence)
    {
        if (protocol.Observables.IsDefaultOrEmpty || correspondence.Atoms.IsDefaultOrEmpty ||
            correspondence.Atoms.Select(atom => atom.ResultAtomIndex).Distinct().Count() !=
                correspondence.Atoms.Length)
            return ImmutableArray<ResolvedEquilibrationObservable>.Empty;

        ImmutableArray<int> Select(string selector)
        {
            if (selector == "none") return ImmutableArray<int>.Empty;
            Func<AtomCorrespondence, bool>? predicate = selector switch
            {
                "protein-backbone" => atom => atom.MoleculeRole == MoleculeRoleKind.Protein && atom.AtomRole == AtomRoleKind.Backbone,
                "upper-lipid-head" => atom => atom.MoleculeRole == MoleculeRoleKind.Lipid && atom.AtomRole == AtomRoleKind.Head && atom.PhysicalSide == LeafletSide.Upper,
                "lower-lipid-head" => atom => atom.MoleculeRole == MoleculeRoleKind.Lipid && atom.AtomRole == AtomRoleKind.Head && atom.PhysicalSide == LeafletSide.Lower,
                "all-lipid-head" => atom => atom.MoleculeRole == MoleculeRoleKind.Lipid && atom.AtomRole == AtomRoleKind.Head,
                "water-oxygen" => atom => atom.MoleculeRole == MoleculeRoleKind.Water && atom.Element == "O",
                "all-ions" => atom => atom.MoleculeRole == MoleculeRoleKind.Ion,
                "protein-or-lipid-heavy" => atom => (atom.MoleculeRole is MoleculeRoleKind.Protein or MoleculeRoleKind.Lipid) && atom.Element != "H",
                _ => null
            };
            return predicate is null ? default : correspondence.Atoms.Where(predicate)
                .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        }

        var resolved = ImmutableArray.CreateBuilder<ResolvedEquilibrationObservable>();
        foreach (var declared in protocol.Observables)
        {
            var first = Select(declared.AtomSelector);
            var second = Select(declared.ComparisonSelector);
            var third = Select(declared.ThirdSelector);
            if (first.IsDefault || second.IsDefault || third.IsDefault ||
                (declared.AtomSelector != "none" && first.IsEmpty) ||
                (declared.ComparisonSelector != "none" && second.IsEmpty) ||
                (declared.ThirdSelector != "none" && third.IsEmpty))
                return ImmutableArray<ResolvedEquilibrationObservable>.Empty;
            resolved.Add(new ResolvedEquilibrationObservable(declared.Name, declared.Unit,
                declared.Scope, declared.Category, declared.Method, first, second, third,
                declared.DistanceCutoffAngstrom));
        }
        return resolved.ToImmutable();
    }

    private static bool ValidAdequacyResult(EquilibrationProtocol protocol, EquilibrationObservations observed)
    {
        if (observed.Samples.IsDefaultOrEmpty || observed.Samples.Length > protocol.MaximumSampleCount ||
            observed.ObservationAssessments.Length != protocol.Observables.Length ||
            observed.ObservationAssessments.Select(item => item.ObservableName)
                .Distinct(StringComparer.Ordinal).Count() != protocol.Observables.Length)
            return false;
        var tracedWindows = observed.Windows.Select(window => window.Name).ToHashSet(StringComparer.Ordinal);
        if (observed.Samples.Any(sample => !tracedWindows.Contains(sample.WindowName) ||
            sample.Step < 0 || sample.Measurements.Length != protocol.Observables.Length ||
            protocol.Observables.Any(declared => !sample.Measurements.Any(value =>
                value.Name == declared.Name && value.Unit == declared.Unit &&
                value.Scope == declared.Scope && double.IsFinite(value.Value))))) return false;
        foreach (var assessment in observed.ObservationAssessments)
        {
            var rule = protocol.SufficiencyRules.SingleOrDefault(item => item.ObservableName == assessment.ObservableName);
            var adequateSamples = observed.Samples.Count(sample =>
                sample.WindowName == protocol.Stages[^1].Name ||
                sample.WindowName.StartsWith(protocol.ExtensionWindow.Name + "-", StringComparison.Ordinal));
            if (rule is null || assessment.SampleCount != adequateSamples ||
                !double.IsFinite(assessment.EffectiveBlockCount) || assessment.EffectiveBlockCount < 0 ||
                (assessment.FirstVsLastBlockMeanDifference is double drift && !double.IsFinite(drift)) ||
                (assessment.LagOneBlockCorrelation is double correlation && !double.IsFinite(correlation)))
                return false;
            if (assessment.Sufficient &&
                (assessment.EffectiveBlockCount < rule.MinimumEffectiveBlocks ||
                 assessment.FirstVsLastBlockMeanDifference is not double observedDrift ||
                 Math.Abs(observedDrift) > rule.MaximumAbsoluteFirstVsLastBlockMeanDifference ||
                 assessment.LagOneBlockCorrelation is not double observedCorrelation ||
                 Math.Abs(observedCorrelation) > rule.MaximumAbsoluteLagOneBlockCorrelation))
                return false;
        }
        var allSufficient = observed.ObservationAssessments.All(item => item.Sufficient);
        return observed.ObservationAdequacy switch
        {
            EquilibrationObservationAdequacy.Adequate => allSufficient,
            EquilibrationObservationAdequacy.InsufficientAtBound => !allSufficient && observed.ExtensionsPerformed == protocol.MaximumExtensions,
            _ => false
        };
    }

    private static StageExecutionState State(string attemptId, string stageId,
        StageExecutionStanding standing, string message, double? progress = null)
        => new(attemptId, stageId, StageKind.Equilibration, standing, message, progress, DateTimeOffset.UtcNow);
}
