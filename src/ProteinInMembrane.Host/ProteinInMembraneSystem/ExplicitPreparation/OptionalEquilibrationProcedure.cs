using System.Collections.Immutable;
using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text.Json;

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
        var protocolSha256 = protocol is not null &&
            EquilibrationProtocolFingerprint.TryCompute(protocol, out var digest) ? digest : null;
        if (minimized.Kind != StageKind.Minimization ||
            minimized.Observation.Kind != StageKind.Minimization ||
            minimized.Observation.StageId != minimized.Id ||
            minimized.Observation.AttemptId != minimized.Attempt.Id ||
            minimized.Molecule.Id != minimized.Id ||
            minimized.Correspondence.ResultId != minimized.Id ||
            minimized.PolicyId != policy.Id ||
            !PreparationPolicyFingerprint.Matches(minimized.Attempt, policy) ||
            minimized.Molecule.TopologyPath is null || minimized.Molecule.SystemXmlPath is null ||
            minimized.Molecule.StateXmlPath is null || minimized.Molecule.TopologySha256 is null ||
            minimized.Molecule.SystemXmlSha256 is null || minimized.Molecule.StateXmlSha256 is null ||
            minimized.Molecule.CorrespondencePath is null ||
            minimized.Molecule.CorrespondenceSha256 is null ||
            policy.LocalStateObservation is null || protocol is null ||
            !ValidProtocol(protocol, resolvedObservables) || protocolSha256 is null ||
            !minimized.Correspondence.Complete)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "The source minimized stage or exact applicable optional procedure is unavailable."),
                ImmutableArray<ScientificFinding>.Empty);

        var plannedFrames = protocol.Stages.Sum(stage =>
            (long)Math.Ceiling((double)stage.Steps / stage.ReportIntervalSteps)) +
            (long)protocol.MaximumExtensions *
            (long)Math.Ceiling((double)protocol.ExtensionWindow.Steps /
                protocol.ExtensionWindow.ReportIntervalSteps);
        if (plannedFrames * minimized.Molecule.AtomCount * 24L > protocol.MaximumFrameBytes)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.ResourceRefused,
                    "The declared optional frame series exceeds its explicit storage budget."),
                ImmutableArray<ScientificFinding>.Empty);

        var proteinBackbone = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Protein && atom.AtomRole == AtomRoleKind.Backbone)
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        var proteinHeavy = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Protein && atom.Element != "H")
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        var proteinBackboneHeavy = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Protein &&
                atom.AtomRole == AtomRoleKind.Backbone && atom.Element != "H")
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        var lipidHeavy = minimized.Correspondence.Atoms
            .Where(atom => atom.MoleculeRole == MoleculeRoleKind.Lipid && atom.Element != "H")
            .Select(atom => atom.ResultAtomIndex).ToImmutableArray();
        bool HasProteinTarget(EquilibrationStageControl stage) => stage.ProteinRestraintSelector switch
        {
            "backbone" => !proteinBackbone.IsEmpty,
            "protein-heavy" => !proteinHeavy.IsEmpty,
            "backbone-heavy" or "protein-ca" => !proteinBackboneHeavy.IsEmpty,
            _ => false
        };
        if (protocol.Stages.Any(stage => stage.ProteinRestraintKjMolNm2 > 0 && !HasProteinTarget(stage)) ||
            (protocol.Stages.Any(stage => stage.LipidRestraintKjMolNm2 > 0) && lipidHeavy.IsDefaultOrEmpty) ||
            proteinBackbone.Concat(proteinHeavy).Concat(lipidHeavy)
                .Any(index => index < 0 || index >= minimized.Molecule.AtomCount))
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
                protocol.Id, proteinBackbone, proteinHeavy, lipidHeavy,
                resolvedObservables, protocol, protocolSha256));
        WorkerResult<EquilibrationObservations> result;
        try
        {
            result = await _worker.EquilibrateAsync(request,
                progress is null ? null : new EquilibrationProgressForwarder(
                    minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId, protocol, progress),
                cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Optional equilibration stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        }
        if (IsCorrelatedStop(result, request.RequestId, minimized.Attempt.StudyRevisionId,
                minimized.Attempt.Id, stageId, cancellationToken))
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Optional equilibration stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        var equilibrationCorrelated = HasExactContext(result, request.RequestId,
            minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId);
        if (!equilibrationCorrelated || result.Standing != WorkerResultStanding.Observed ||
            result.Observations is null)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, !equilibrationCorrelated
                    ? StageExecutionStanding.Unobserved : result.Standing switch
                {
                    WorkerResultStanding.Stopped or WorkerResultStanding.Unobserved => StageExecutionStanding.Unobserved,
                    _ => StageExecutionStanding.Failed
                }, result.FailureMessage ?? "Optional equilibration did not establish an observed completed stage."),
                ImmutableArray<ScientificFinding>.Empty);

        var observed = result.Observations;
        var equilibratedState = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedStateXml");
        var equilibratedCoordinates = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedCif");
        var equilibratedSystem = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedSystemXml");
        var equilibratedTopology = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "equilibratedTopologyJson");
        var framesManifest = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "sampledFramesManifest");
        var framesPositions = result.Artifacts.IsDefault ? null :
            result.Artifacts.FirstOrDefault(artifact => artifact.Role == "sampledPositionsF64");
        if (equilibratedState is null || equilibratedCoordinates is null || equilibratedSystem is null ||
            equilibratedTopology is null ||
            !ArtifactMatches(equilibratedState) || !ArtifactMatches(equilibratedCoordinates) ||
            !ArtifactMatches(equilibratedSystem) || !ArtifactMatches(equilibratedTopology) ||
            observed.FinalAtomCount != minimized.Molecule.AtomCount ||
            observed.Windows.IsDefault || observed.Windows.Any(window => window is null) ||
            observed.ExtensionsPerformed < 0 || observed.ExtensionsPerformed > protocol.MaximumExtensions ||
            observed.Windows.Length != protocol.Stages.Length + observed.ExtensionsPerformed ||
            observed.Windows.Where((window, index) =>
                window.Name != (index < protocol.Stages.Length ? protocol.Stages[index].Name :
                    $"{protocol.ExtensionWindow.Name}-{index - protocol.Stages.Length + 1}") ||
                window.RequestedSteps != (index < protocol.Stages.Length ? protocol.Stages[index].Steps : protocol.ExtensionWindow.Steps) ||
                window.CompletedSteps != window.RequestedSteps ||
                window.TemperatureKelvin is not double temperature ||
                !double.IsFinite(temperature) || temperature <= 0 ||
                !PressureMatches(window.PressureBar,
                    index < protocol.Stages.Length ? protocol.Stages[index].PressureBar :
                    protocol.ExtensionWindow.PressureBar) ||
                window.Measurements.IsDefaultOrEmpty ||
                window.Warnings.IsDefault ||
                window.Measurements.Any(value => !double.IsFinite(value.Value))).Any() ||
            !observed.FinalObservationUnrestrained ||
            observed.NumericalWarnings.IsDefault || !observed.NumericalWarnings.IsEmpty ||
            !ValidAdequacyResult(protocol, observed) ||
            observed.Termination != StageTermination.Completed)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "The complete declared optional procedure and trustworthy final state were not observed."),
                ImmutableArray<ScientificFinding>.Empty);

        EquilibrationFrameSeries? frameSeries;
        try
        {
            frameSeries = await ReadFrameSeriesAsync(framesManifest, framesPositions,
                request.Payload, equilibratedCoordinates, equilibratedTopology,
                equilibratedSystem, equilibratedState,
                observed.FinalAtomCount, observed.Samples, cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Frame read-back stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        }
        if (frameSeries is null)
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Unobserved,
                    "The optional procedure's sampled periodic frame series is missing or does not match its declared identity and byte budget."),
                ImmutableArray<ScientificFinding>.Empty);

        var observationRequest = new ScientificWorkRequest<StageObservationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new StageObservationPayload(minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId,
                equilibratedCoordinates.Path, equilibratedCoordinates.Sha256,
                equilibratedTopology.Path, equilibratedTopology.Sha256,
                equilibratedSystem.Path, equilibratedSystem.Sha256,
                equilibratedState.Path, equilibratedState.Sha256,
                StageKind.Equilibration, minimized.Molecule.CorrespondencePath!,
                minimized.Molecule.CorrespondenceSha256!, policy.LocalStateObservation,
                policy.StageProteinGeometryMeasurement));
        WorkerResult<StageObservationObservations> stageResult;
        try
        {
            stageResult = await _worker.ObserveStageAsync(observationRequest, cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Stage remeasurement was stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        }
        if (IsCorrelatedStop(stageResult, observationRequest.RequestId,
                minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId, cancellationToken))
            return new StageOperationResult(null,
                State(minimized.Attempt.Id, stageId, StageExecutionStanding.Stopped,
                    "Stage remeasurement was stopped; the completed minimized stage remains available."),
                ImmutableArray<ScientificFinding>.Empty);
        if (!HasExactContext(stageResult, observationRequest.RequestId,
                minimized.Attempt.StudyRevisionId, minimized.Attempt.Id, stageId) ||
            stageResult.Standing != WorkerResultStanding.Observed || stageResult.Observations is null ||
            stageResult.Observations.AtomCount != minimized.Molecule.AtomCount ||
            !stageResult.Observations.AtomOrderMatched || !stageResult.Observations.BondsMatched ||
            stageResult.Observations.Measurements.IsDefault ||
            stageResult.Observations.ContactWarnings.IsDefault ||
            stageResult.Observations.StructuralWarnings.IsDefault ||
            stageResult.Observations.NumericalWarnings.IsDefault ||
            !stageResult.Observations.NumericalWarnings.IsEmpty ||
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
            TopologyPath = equilibratedTopology.Path,
            TopologySha256 = equilibratedTopology.Sha256,
            SystemXmlPath = equilibratedSystem.Path,
            SystemXmlSha256 = equilibratedSystem.Sha256,
            StateXmlPath = equilibratedState.Path,
            StateXmlSha256 = equilibratedState.Sha256,
            // Pressure control can change the cell. The exact output vectors are
            // in the stage topology and State; an inherited source-cell label
            // would describe the wrong completed molecule.
            CellDescription = null
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
            stageResult.Observations.LocalState, stageResult.Observations.ProteinGeometry,
            observed.Windows, frameSeries);
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
            protocol.MaximumSampleCount <= 0 || protocol.MaximumSampleCount > 100_000 ||
            protocol.MaximumFrameBytes <= 0 || protocol.MaximumFrameBytes > 16L * 1024 * 1024 * 1024 ||
            protocol.ExtensionWindow is null ||
            protocol.Observables.IsDefaultOrEmpty || resolved.IsDefaultOrEmpty ||
            protocol.SufficiencyRules.IsDefaultOrEmpty ||
            !double.IsFinite(protocol.TargetTemperatureKelvin) || protocol.TargetTemperatureKelvin != 303.0 ||
            protocol.RequiredObservations.IsDefaultOrEmpty || string.IsNullOrWhiteSpace(protocol.ComparisonBasis))
            return false;
        if (protocol.Stages.Append(protocol.ExtensionWindow).Any(stage => stage.Steps <= 0 || stage.ReportIntervalSteps <= 0 ||
            stage.TimestepPicoseconds <= 0 || !double.IsFinite(stage.TimestepPicoseconds) ||
            !double.IsFinite(stage.TemperatureKelvin) || stage.TemperatureKelvin <= 0 ||
            (stage.InitialTemperatureKelvin is double initial &&
                (!double.IsFinite(initial) || initial <= 0 || initial == stage.TemperatureKelvin ||
                 stage.Steps < 2 || stage.PressureBar is not null)) ||
            !double.IsFinite(stage.FrictionPerPicosecond) || stage.FrictionPerPicosecond <= 0 ||
            stage.ProteinRestraintSelector is not
                ("backbone" or "protein-heavy" or "backbone-heavy" or "protein-ca") ||
            !double.IsFinite(stage.ProteinRestraintKjMolNm2) || stage.ProteinRestraintKjMolNm2 < 0 ||
            !double.IsFinite(stage.LipidRestraintKjMolNm2) || stage.LipidRestraintKjMolNm2 < 0 ||
            (stage.PressureBar is null
                ? !string.Equals(stage.PressureMode, "none", StringComparison.OrdinalIgnoreCase) ||
                  stage.BarostatFrequencySteps is not null || stage.SurfaceTensionBarNm is not null
                : !double.IsFinite(stage.PressureBar.Value) || stage.PressureBar.Value < 0 ||
                  stage.BarostatFrequencySteps is null or <= 0 ||
                  stage.SurfaceTensionBarNm is not double surfaceTension ||
                  !double.IsFinite(surfaceTension) ||
                  stage.PressureMode is null ||
                  stage.PressureMode.ToLowerInvariant() is not
                      ("xyisotropiczfree" or "xyisotropiczfixed" or "xyanisotropiczfree"))))
            return false;
        var last = protocol.Stages[^1];
        var extension = protocol.ExtensionWindow;
        if (last.ProteinRestraintKjMolNm2 != 0 || last.LipidRestraintKjMolNm2 != 0 ||
            extension.ProteinRestraintKjMolNm2 != 0 || extension.LipidRestraintKjMolNm2 != 0 ||
            last.InitialTemperatureKelvin is not null || extension.InitialTemperatureKelvin is not null ||
            last.TemperatureKelvin != protocol.TargetTemperatureKelvin ||
            extension.TemperatureKelvin != last.TemperatureKelvin ||
            extension.PressureBar != last.PressureBar || extension.PressureMode != last.PressureMode ||
            extension.BarostatFrequencySteps != last.BarostatFrequencySteps ||
            extension.SurfaceTensionBarNm != last.SurfaceTensionBarNm ||
            extension.FrictionPerPicosecond != last.FrictionPerPicosecond ||
            extension.TimestepPicoseconds != last.TimestepPicoseconds ||
            protocol.Stages.Any(stage => string.IsNullOrWhiteSpace(stage.Name) || stage.Name.Length > 128) ||
            string.IsNullOrWhiteSpace(extension.Name) || extension.Name.Length > 128 ||
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

    private static bool PressureMatches(double? observed, double? declared) =>
        declared is null ? observed is null :
        observed is double value && double.IsFinite(value) &&
        Math.Abs(value - declared.Value) <= 0.000000001 * Math.Max(1, Math.Abs(declared.Value));

    private static bool HasExactContext<TObservation>(WorkerResult<TObservation> result,
        string requestId, string revisionId, string attemptId, string stageId)
        where TObservation : class =>
        result.RequestId == requestId && result.StudyRevisionId == revisionId &&
        result.AttemptId == attemptId && result.StageId == stageId;

    private static bool IsCorrelatedStop<TObservation>(WorkerResult<TObservation> result,
        string requestId, string revisionId, string attemptId, string stageId,
        CancellationToken cancellationToken)
        where TObservation : class =>
        result.Standing == WorkerResultStanding.Stopped && result.RequestId == requestId &&
        (HasExactContext(result, requestId, revisionId, attemptId, stageId) ||
         // The local exchange kills a cancelled one-request process and returns
         // this request-scoped stop before it has a worker terminal with IDs.
         (cancellationToken.IsCancellationRequested && result.FailureCode == "worker-stopped" &&
          result.StudyRevisionId is null && result.AttemptId is null && result.StageId is null));

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
            observed.Samples.Any(sample => sample is null || sample.Measurements.IsDefault ||
                sample.Measurements.Any(value => value is null)) ||
            observed.ObservationAssessments.IsDefault ||
            observed.ObservationAssessments.Any(item => item is null) ||
            observed.ObservationAssessments.Length != protocol.Observables.Length ||
            observed.ObservationAssessments.Select(item => item.ObservableName)
                .Distinct(StringComparer.Ordinal).Count() != protocol.Observables.Length)
            return false;
        var nextSample = 0;
        foreach (var window in observed.Windows.Select((value, index) =>
            (Value: value, Control: index < protocol.Stages.Length
                ? protocol.Stages[index] : protocol.ExtensionWindow)))
        {
            for (var step = Math.Min(window.Control.ReportIntervalSteps, window.Control.Steps);
                 step <= window.Control.Steps; step = (int)Math.Min(
                     (long)step + window.Control.ReportIntervalSteps, window.Control.Steps))
            {
                if (nextSample >= observed.Samples.Length ||
                    observed.Samples[nextSample].WindowName != window.Value.Name ||
                    observed.Samples[nextSample].Step != step)
                    return false;
                nextSample++;
                if (step == window.Control.Steps) break;
            }
        }
        if (nextSample != observed.Samples.Length) return false;
        var tracedWindows = observed.Windows.Select(window => window.Name).ToHashSet(StringComparer.Ordinal);
        if (observed.Samples.Any(sample => !tracedWindows.Contains(sample.WindowName) ||
            sample.Step < 0 || sample.Measurements.Length != protocol.Observables.Length ||
            protocol.Observables.Any(declared => !sample.Measurements.Any(value =>
                value.Name == declared.Name && value.Unit == declared.Unit &&
                value.Scope == declared.Scope && double.IsFinite(value.Value))))) return false;
        // The worker's sufficiency fields are claims about these retained
        // samples. Reconstruct the declared complete blocks before assessment.
        var adequacySamples = observed.Samples.Where(sample =>
            sample.WindowName == protocol.Stages[^1].Name ||
            sample.WindowName.StartsWith(protocol.ExtensionWindow.Name + "-", StringComparison.Ordinal))
            .ToImmutableArray();
        var allSufficient = true;
        foreach (var assessment in observed.ObservationAssessments)
        {
            var rule = protocol.SufficiencyRules.SingleOrDefault(item => item.ObservableName == assessment.ObservableName);
            var calculated = rule is null ? null : CalculateAdequacy(rule, adequacySamples);
            if (calculated is null || assessment.SampleCount != calculated.SampleCount ||
                !NearlyEqual(assessment.EffectiveBlockCount, calculated.EffectiveBlockCount) ||
                !NearlyEqual(assessment.FirstVsLastBlockMeanDifference,
                    calculated.FirstVsLastBlockMeanDifference) ||
                !NearlyEqual(assessment.LagOneBlockCorrelation, calculated.LagOneBlockCorrelation) ||
                assessment.Sufficient != calculated.Sufficient)
                return false;
            allSufficient &= calculated.Sufficient;
        }
        return observed.ObservationAdequacy switch
        {
            EquilibrationObservationAdequacy.Adequate => allSufficient,
            EquilibrationObservationAdequacy.InsufficientAtBound => !allSufficient && observed.ExtensionsPerformed == protocol.MaximumExtensions,
            _ => false
        };
    }

    private static EquilibrationObservationAssessment? CalculateAdequacy(
        EquilibrationSufficiencyRule rule, ImmutableArray<EquilibrationSample> samples)
    {
        var blockCount = samples.Length / rule.BlockSizeSamples;
        if (blockCount < 3)
            return new EquilibrationObservationAssessment(rule.ObservableName, samples.Length,
                0, null, null, false);
        var blocks = new double[blockCount];
        for (var block = 0; block < blockCount; block++)
        {
            var sum = 0.0;
            for (var offset = 0; offset < rule.BlockSizeSamples; offset++)
            {
                var sample = samples[block * rule.BlockSizeSamples + offset];
                var value = sample.Measurements.First(item => item.Name == rule.ObservableName).Value;
                sum += value;
            }
            blocks[block] = sum / rule.BlockSizeSamples;
            if (!double.IsFinite(blocks[block])) return null;
        }
        var blockSum = 0.0;
        foreach (var block in blocks) blockSum += block;
        var mean = blockSum / blocks.Length;
        if (!double.IsFinite(mean)) return null;
        var variance = 0.0;
        foreach (var block in blocks) variance += (block - mean) * (block - mean);
        if (!double.IsFinite(variance)) return null;
        var correlation = 0.0;
        if (variance != 0)
        {
            var lagSum = 0.0;
            for (var index = 0; index < blocks.Length - 1; index++)
                lagSum += (blocks[index] - mean) * (blocks[index + 1] - mean);
            correlation = Math.Clamp(lagSum / variance, -1.0, 1.0);
        }
        var positive = Math.Max(0, correlation);
        var effective = blocks.Length * (1 - positive) / (1 + positive);
        var drift = blocks[^1] - blocks[0];
        if (!double.IsFinite(correlation) || !double.IsFinite(effective) || !double.IsFinite(drift))
            return null;
        var sufficient = effective >= rule.MinimumEffectiveBlocks &&
            Math.Abs(drift) <= rule.MaximumAbsoluteFirstVsLastBlockMeanDifference &&
            Math.Abs(correlation) <= rule.MaximumAbsoluteLagOneBlockCorrelation;
        return new EquilibrationObservationAssessment(rule.ObservableName, samples.Length,
            effective, drift, correlation, sufficient);
    }

    // Allow only transport-scale differences between Python and C# sums.
    private static bool NearlyEqual(double? reported, double? calculated) =>
        reported is null || calculated is null
            ? reported is null && calculated is null
            : double.IsFinite(reported.Value) && double.IsFinite(calculated.Value) &&
              Math.Abs(reported.Value - calculated.Value) <=
              0.000000000001 * Math.Max(1, Math.Max(Math.Abs(reported.Value), Math.Abs(calculated.Value)));

    private static StageExecutionState State(string attemptId, string stageId,
        StageExecutionStanding standing, string message, double? progress = null)
        => new(attemptId, stageId, StageKind.Equilibration, standing, message, progress, DateTimeOffset.UtcNow);

    private static bool ArtifactMatches(WorkerArtifact artifact)
    {
        if (string.IsNullOrWhiteSpace(artifact.Path) || artifact.Sha256 is not { Length: 64 } ||
            !artifact.Sha256.All(Uri.IsHexDigit))
            return false;
        try
        {
            using var stream = File.OpenRead(artifact.Path);
            return Convert.ToHexString(SHA256.HashData(stream))
                .Equals(artifact.Sha256, StringComparison.OrdinalIgnoreCase);
        }
        catch (IOException) { return false; }
        catch (UnauthorizedAccessException) { return false; }
        catch (ArgumentException) { return false; }
    }

    private static async Task<EquilibrationFrameSeries?> ReadFrameSeriesAsync(
        WorkerArtifact? manifestArtifact, WorkerArtifact? positionsArtifact,
        EquilibrationPayload payload, WorkerArtifact finalCoordinates,
        WorkerArtifact finalTopology, WorkerArtifact finalSystem, WorkerArtifact finalState,
        int expectedAtomCount,
        ImmutableArray<EquilibrationSample> samples,
        CancellationToken cancellationToken)
    {
        if (manifestArtifact is null || positionsArtifact is null ||
            !File.Exists(positionsArtifact.Path) ||
            !string.Equals(Path.GetDirectoryName(manifestArtifact.Path),
                Path.GetDirectoryName(positionsArtifact.Path), StringComparison.Ordinal) ||
            positionsArtifact.Sha256 is not { Length: 64 } ||
            !positionsArtifact.Sha256.All(Uri.IsHexDigit))
            return null;
        try
        {
            if (new FileInfo(manifestArtifact.Path).Length > 64L * 1024 * 1024 ||
                !ArtifactMatches(manifestArtifact))
                return null;
            await using var manifestStream = File.OpenRead(manifestArtifact.Path);
            using var document = await JsonDocument.ParseAsync(manifestStream,
                cancellationToken: cancellationToken);
            var manifest = document.RootElement;
            var fullAtomCount = manifest.GetProperty("atomCount").GetInt32();
            var frameBytes = fullAtomCount * 24L;
            var length = new FileInfo(positionsArtifact.Path).Length;
            var frames = manifest.GetProperty("frames");
            if (manifest.GetProperty("schemaVersion").GetString() !=
                    "protein-in-membrane.equilibration-frames.v1" ||
                manifest.GetProperty("encoding").GetString() != "float64-le-xyz-angstrom" ||
                manifest.GetProperty("studyRevisionId").GetString() != payload.StudyRevisionId ||
                manifest.GetProperty("attemptId").GetString() != payload.AttemptId ||
                manifest.GetProperty("stageId").GetString() != payload.StageId ||
                manifest.GetProperty("sourceMinimizedStageId").GetString() != payload.SourceMinimizedStageId ||
                manifest.GetProperty("validatedPolicyId").GetString() != payload.ValidatedPolicyId ||
                !string.Equals(manifest.GetProperty("protocolSha256").GetString(), payload.ProtocolSha256,
                    StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("sourceCoordinateSha256").GetString(),
                    payload.TopologyCifSha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("sourceTopologySha256").GetString(),
                    payload.TopologyJsonSha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("sourceSystemSha256").GetString(),
                    payload.SystemXmlSha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("sourceStateSha256").GetString(),
                    payload.MinimizedStateXmlSha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("finalCoordinateSha256").GetString(),
                    finalCoordinates.Sha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("finalTopologySha256").GetString(),
                    finalTopology.Sha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("finalSystemSha256").GetString(),
                    finalSystem.Sha256, StringComparison.OrdinalIgnoreCase) ||
                !string.Equals(manifest.GetProperty("finalStateSha256").GetString(),
                    finalState.Sha256, StringComparison.OrdinalIgnoreCase) ||
                manifest.GetProperty("maximumFrameBytes").GetInt64() != payload.Protocol.MaximumFrameBytes ||
                manifest.GetProperty("positionsFile").GetString() != Path.GetFileName(positionsArtifact.Path) ||
                !string.Equals(manifest.GetProperty("positionsSha256").GetString(),
                    positionsArtifact.Sha256, StringComparison.OrdinalIgnoreCase) ||
                manifest.GetProperty("positionsByteLength").GetInt64() != length ||
                fullAtomCount != expectedAtomCount || fullAtomCount <= 0 || frameBytes <= 0 ||
                manifest.GetProperty("frameCount").GetInt32() != samples.Length ||
                frames.ValueKind != JsonValueKind.Array || frames.GetArrayLength() != samples.Length ||
                samples.Length > payload.Protocol.MaximumSampleCount ||
                length != frameBytes * samples.Length || length > payload.Protocol.MaximumFrameBytes)
                return null;

            await using var positions = File.OpenRead(positionsArtifact.Path);
            using var wholeDigest = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            var buffer = new byte[64 * 1024];
            for (var frameIndex = 0; frameIndex < samples.Length; frameIndex++)
            {
                var frame = frames[frameIndex];
                var expectedDigest = frame.GetProperty("positionsSha256").GetString();
                var box = frame.GetProperty("boxVectorsAngstrom");
                if (frame.GetProperty("windowName").GetString() != samples[frameIndex].WindowName ||
                    frame.GetProperty("step").GetInt32() != samples[frameIndex].Step ||
                    frame.GetProperty("offsetBytes").GetInt64() != positions.Position ||
                    frame.GetProperty("lengthBytes").GetInt64() != frameBytes ||
                    expectedDigest is not { Length: 64 } || !expectedDigest.All(Uri.IsHexDigit) ||
                    !ValidFrameBox(box))
                    return null;
                using var frameDigest = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
                var remaining = frameBytes;
                while (remaining > 0)
                {
                    var read = (int)Math.Min(buffer.Length, remaining);
                    await positions.ReadExactlyAsync(buffer.AsMemory(0, read), cancellationToken);
                    for (var offset = 0; offset < read; offset += sizeof(double))
                        if (!double.IsFinite(BinaryPrimitives.ReadDoubleLittleEndian(
                                buffer.AsSpan(offset, sizeof(double)))))
                            return null;
                    frameDigest.AppendData(buffer, 0, read);
                    wholeDigest.AppendData(buffer, 0, read);
                    remaining -= read;
                }
                if (!Convert.ToHexString(frameDigest.GetHashAndReset()).Equals(expectedDigest,
                        StringComparison.OrdinalIgnoreCase))
                    return null;
            }
            if (positions.Position != length ||
                !Convert.ToHexString(wholeDigest.GetHashAndReset()).Equals(positionsArtifact.Sha256,
                    StringComparison.OrdinalIgnoreCase))
                return null;
            using var finalTopologyDocument = JsonDocument.Parse(
                await File.ReadAllBytesAsync(finalTopology.Path, cancellationToken));
            var finalBox = finalTopologyDocument.RootElement.GetProperty("boxVectorsAngstrom");
            var lastFrameBox = frames[samples.Length - 1].GetProperty("boxVectorsAngstrom");
            if (!ValidFrameBox(finalBox)) return null;
            for (var axis = 0; axis < 3; axis++)
                for (var component = 0; component < 3; component++)
                    if (Math.Abs(finalBox[axis][component].GetDouble() -
                            lastFrameBox[axis][component].GetDouble()) > 1e-8)
                        return null;
            return new EquilibrationFrameSeries(manifestArtifact.Path,
                manifestArtifact.Sha256, positionsArtifact.Path, positionsArtifact.Sha256,
                length, samples.Length);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException or
                                         JsonException or KeyNotFoundException or InvalidOperationException or
                                         ArgumentException or OverflowException)
        {
            return null;
        }
    }

    private static bool ValidFrameBox(JsonElement box)
    {
        if (box.ValueKind != JsonValueKind.Array || box.GetArrayLength() != 3)
            return false;
        var values = new double[3, 3];
        for (var axis = 0; axis < 3; axis++)
        {
            if (box[axis].ValueKind != JsonValueKind.Array || box[axis].GetArrayLength() != 3)
                return false;
            for (var component = 0; component < 3; component++)
                if (box[axis][component].ValueKind != JsonValueKind.Number ||
                    !box[axis][component].TryGetDouble(out values[axis, component]) ||
                    !double.IsFinite(values[axis, component]))
                    return false;
        }
        var determinant = values[0, 0] * (values[1, 1] * values[2, 2] - values[1, 2] * values[2, 1]) -
            values[0, 1] * (values[1, 0] * values[2, 2] - values[1, 2] * values[2, 0]) +
            values[0, 2] * (values[1, 0] * values[2, 1] - values[1, 1] * values[2, 0]);
        return double.IsFinite(determinant) && determinant > 0;
    }

    private sealed class EquilibrationProgressForwarder : IProgress<EquilibrationWorkProgress>
    {
        private readonly string _revisionId;
        private readonly string _attemptId;
        private readonly string _stageId;
        private readonly IProgress<StageExecutionState> _output;
        private readonly (string Name, int Steps)[] _windows;
        private readonly int[] _completed;
        private readonly long _totalSteps;
        private readonly object _gate = new();

        public EquilibrationProgressForwarder(string revisionId, string attemptId, string stageId,
            EquilibrationProtocol protocol, IProgress<StageExecutionState> output)
        {
            _revisionId = revisionId;
            _attemptId = attemptId;
            _stageId = stageId;
            _output = output;
            _windows = protocol.Stages.Select(stage => (stage.Name, stage.Steps))
                .Concat(Enumerable.Range(1, protocol.MaximumExtensions).Select(index =>
                    ($"{protocol.ExtensionWindow.Name}-{index}", protocol.ExtensionWindow.Steps)))
                .ToArray();
            _completed = new int[_windows.Length];
            _totalSteps = _windows.Sum(window => (long)window.Steps);
        }

        public void Report(EquilibrationWorkProgress observed)
        {
            if (observed.StudyRevisionId != _revisionId || observed.AttemptId != _attemptId ||
                observed.StageId != _stageId) return;
            lock (_gate)
            {
                var index = Array.FindIndex(_windows, window => window.Name == observed.Window);
                if (index < 0 || observed.RequestedSteps != _windows[index].Steps ||
                    observed.CompletedSteps <= _completed[index] ||
                    observed.CompletedSteps > observed.RequestedSteps ||
                    _completed.Skip(index + 1).Any(value => value > 0) ||
                    _completed.Take(index).Where((value, prior) => value != _windows[prior].Steps).Any())
                    return;
                _completed[index] = observed.CompletedSteps;
                var fraction = _totalSteps > 0 ? (double)_completed.Sum(value => (long)value) / _totalSteps : 0;
                _output.Report(State(_attemptId, _stageId, StageExecutionStanding.Running,
                    $"Optional equilibration: {observed.Window} {observed.CompletedSteps}/{observed.RequestedSteps} steps observed.",
                    fraction));
            }
        }
    }
}
