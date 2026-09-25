using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation;

/// <summary>Establishes one factual completed minimized stage, not qualification.</summary>
public sealed class Minimization
{
    private readonly IMinimizationWork _worker;

    public Minimization(IMinimizationWork worker) => _worker = worker;

    public async Task<StageOperationResult> RunAsync(
        ConstructedExplicitSystem source,
        ApplicablePreparationPolicy policy,
        string workingDirectory,
        IProgress<StageExecutionState>? progress,
        CancellationToken cancellationToken)
    {
        var stageId = Guid.NewGuid().ToString("N");
        var running = State(source.Attempt.Id, stageId, StageExecutionStanding.Running, "Required minimization is running.");
        if (source.Attempt.PolicyId != policy.Id ||
            source.Molecule.TopologyPath is null || source.Molecule.SystemXmlPath is null ||
            source.Molecule.StateXmlPath is null || source.Molecule.TopologySha256 is null ||
            source.Molecule.SystemXmlSha256 is null || source.Molecule.StateXmlSha256 is null ||
            source.Molecule.CorrespondencePath is null || source.Molecule.CorrespondenceSha256 is null ||
            policy.MaximumMinimizationIterations <= 0 ||
            !double.IsFinite(policy.FinalUnrestrainedRmsForceTargetKjMolNm) ||
            policy.FinalUnrestrainedRmsForceTargetKjMolNm != 10.0)
            return new StageOperationResult(null,
                State(source.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "The exact constructed system or declared required minimization policy is unavailable."),
                ImmutableArray<ScientificFinding>.Empty);

        progress?.Report(running);
        var request = new ScientificWorkRequest<MinimizationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new MinimizationPayload(source.Attempt.StudyRevisionId, source.Attempt.Id, stageId,
                source.Molecule.CoordinatePath, source.Molecule.CoordinateSha256,
                source.Molecule.TopologyPath, source.Molecule.TopologySha256,
                source.Molecule.SystemXmlPath, source.Molecule.SystemXmlSha256,
                source.Molecule.StateXmlPath, source.Molecule.StateXmlSha256,
                policy.MaximumMinimizationIterations,
                policy.FinalUnrestrainedRmsForceTargetKjMolNm));
        WorkerResult<MinimizationObservations> result;
        try
        {
            result = await _worker.MinimizeAsync(request, cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            return new StageOperationResult(null,
                State(source.Attempt.Id, stageId, StageExecutionStanding.Stopped, "Minimization was stopped without a completed stage."),
                ImmutableArray<ScientificFinding>.Empty);
        }
        if (result.RequestId != request.RequestId || result.StudyRevisionId != source.Attempt.StudyRevisionId ||
            result.AttemptId != source.Attempt.Id || result.StageId != stageId ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return new StageOperationResult(null, FailureState(source.Attempt.Id, stageId, result.Standing,
                result.FailureMessage ?? "Minimization did not establish an observed result."), ImmutableArray<ScientificFinding>.Empty);

        var observed = result.Observations;
        var minimizedState = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "minimizedStateXml");
        var minimizedCoordinates = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "minimizedCif");
        if (minimizedState is null || minimizedCoordinates is null ||
            !observed.FinalTreatmentUnrestrained ||
            observed.FinalAtomCount != source.Molecule.AtomCount ||
            !double.IsFinite(observed.InitialPotentialEnergyKjMol) ||
            !double.IsFinite(observed.FinalRmsForceKjMolNm) ||
            observed.FinalRmsForceKjMolNm > policy.FinalUnrestrainedRmsForceTargetKjMolNm ||
            !double.IsFinite(observed.FinalPotentialEnergyKjMol) ||
            !observed.NumericalWarnings.IsDefaultOrEmpty ||
            observed.Termination != StageTermination.Converged)
            return new StageOperationResult(null,
                State(source.Attempt.Id, stageId, StageExecutionStanding.Failed,
                    "Final unrestrained convergence, numerical finiteness or molecular correspondence was not observed."),
                ImmutableArray<ScientificFinding>.Empty);

        var observationRequest = new ScientificWorkRequest<StageObservationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new StageObservationPayload(source.Attempt.StudyRevisionId, source.Attempt.Id, stageId,
                minimizedCoordinates.Path, minimizedCoordinates.Sha256,
                source.Molecule.TopologyPath, source.Molecule.TopologySha256,
                source.Molecule.SystemXmlPath, source.Molecule.SystemXmlSha256,
                minimizedState.Path, minimizedState.Sha256,
                StageKind.Minimization, source.Molecule.CorrespondencePath,
                source.Molecule.CorrespondenceSha256, policy.LocalStateObservation,
                policy.StageProteinGeometryMeasurement));
        var stageResult = await _worker.ObserveStageAsync(observationRequest, cancellationToken);
        if (stageResult.RequestId != observationRequest.RequestId ||
            stageResult.StudyRevisionId != source.Attempt.StudyRevisionId ||
            stageResult.AttemptId != source.Attempt.Id || stageResult.StageId != stageId ||
            stageResult.Standing != WorkerResultStanding.Observed || stageResult.Observations is null ||
            stageResult.Observations.AtomCount != source.Molecule.AtomCount ||
            !stageResult.Observations.AtomOrderMatched || !stageResult.Observations.BondsMatched ||
            !stageResult.Observations.NumericalWarnings.IsDefaultOrEmpty)
            return new StageOperationResult(null,
                State(source.Attempt.Id, stageId, StageExecutionStanding.Unobserved,
                    "The minimized coordinates exist, but attributable stage observation is unavailable."),
                ImmutableArray<ScientificFinding>.Empty);

        var molecule = source.Molecule with
        {
            Id = stageId,
            CoordinatePath = minimizedCoordinates.Path,
            CoordinateSha256 = minimizedCoordinates.Sha256,
            StateXmlPath = minimizedState.Path,
            StateXmlSha256 = minimizedState.Sha256
        };
        var measurements = stageResult.Observations.Measurements.AddRange(ImmutableArray.Create(
            new MeasuredValue("initialPotentialEnergy", observed.InitialPotentialEnergyKjMol, "kJ/mol", "minimization"),
            new MeasuredValue("finalPotentialEnergy", observed.FinalPotentialEnergyKjMol, "kJ/mol", "minimization"),
            new MeasuredValue("finalRmsForce", observed.FinalRmsForceKjMolNm, "kJ mol^-1 nm^-1", "final unrestrained")
        ));
        var local = stageResult.Observations.LocalState;
        if (local is not null && !local.Measurements.IsDefault)
            measurements = measurements.AddRange(local.Measurements);
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), stageId, result.Provider?.Name ?? "OpenMM",
            "Observed final unrestrained minimization",
            $"RMS force {observed.FinalRmsForceKjMolNm:G6} kJ mol^-1 nm^-1; termination {observed.Termination}",
            $"Attempt {source.Attempt.Id}; policy {policy.Id}",
            "Completion does not establish suitable membrane phase, thermal equilibration or scientific qualification.",
            EvidenceBearing.Context));
        var observation = new StageObservation(stageId, source.Attempt.Id, StageKind.Minimization,
            measurements, evidence, observed.Termination, result.Provider?.Version ?? "unknown", DateTimeOffset.UtcNow,
            LocalState: stageResult.Observations.LocalState,
            ProteinGeometry: stageResult.Observations.ProteinGeometry);
        var findings = stageResult.Observations.ContactWarnings.Concat(stageResult.Observations.StructuralWarnings)
            .Select(warning => new ScientificFinding(
                Guid.NewGuid().ToString("N"), stageId, evidence[0].Id, warning,
                "Requires stage-specific preparation assessment.", FindingDisposition.Challenges,
                true, DateTimeOffset.UtcNow))
            .ToImmutableArray();
        var completed = new CompletedStage(stageId, source.Attempt, StageKind.Minimization,
            molecule, observation, source.Correspondence with { ResultId = stageId },
            policy.Id, null, findings, DateTimeOffset.UtcNow);
        var state = State(source.Attempt.Id, stageId, StageExecutionStanding.Completed, "Required minimization completed.", 1.0);
        progress?.Report(state);
        return new StageOperationResult(completed, state, findings);
    }

    private static StageExecutionState State(string attemptId, string stageId, StageExecutionStanding standing,
        string message, double? progress = null)
        => new(attemptId, stageId, StageKind.Minimization, standing, message, progress, DateTimeOffset.UtcNow);

    private static StageExecutionState FailureState(string attemptId, string stageId,
        WorkerResultStanding standing, string message)
        => State(attemptId, stageId, standing switch
        {
            WorkerResultStanding.Stopped => StageExecutionStanding.Stopped,
            WorkerResultStanding.Unobserved => StageExecutionStanding.Unobserved,
            _ => StageExecutionStanding.Failed
        }, message);
}
