using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace;

/// <summary>
/// Retains the actual local standing while this host survives. It neither runs a
/// stage nor infers completion from an unfinished artifact or a reopened view.
/// </summary>
public sealed class LocalRunWorkspace
{
    private readonly object _gate = new();
    private StudyRevision? _currentStudy;
    private ImmutableDictionary<string, StudyRevision> _knownStudies =
        ImmutableDictionary<string, StudyRevision>.Empty;
    private PreparationAttempt? _currentAttempt;
    private StageExecutionState? _currentExecution;
    private ImmutableDictionary<string, CompletedStage> _stages = ImmutableDictionary<string, CompletedStage>.Empty;
    private ImmutableDictionary<string, PreparationAssessmentResult> _assessments = ImmutableDictionary<string, PreparationAssessmentResult>.Empty;
    private ImmutableDictionary<string, ScientificFinding> _findings = ImmutableDictionary<string, ScientificFinding>.Empty;

    public void RetainStudy(StudyRevision study)
    {
        ArgumentNullException.ThrowIfNull(study);
        lock (_gate)
        {
            if (_currentStudy is not null && study.Number < _currentStudy.Number)
                throw new InvalidOperationException("An older study cannot replace the current revision.");
            if (_currentStudy is not null && study.Number == _currentStudy.Number &&
                study.Id != _currentStudy.Id ||
                _knownStudies.TryGetValue(study.Id, out var known) && known.Number != study.Number)
                throw new InvalidOperationException("One study revision number or identity cannot be relabelled.");
            _currentStudy = study;
            _knownStudies = _knownStudies.SetItem(study.Id, study);
        }
    }

    public void RetainAttempt(PreparationAttempt attempt)
    {
        ArgumentNullException.ThrowIfNull(attempt);
        lock (_gate)
        {
            if (!_knownStudies.ContainsKey(attempt.StudyRevisionId))
                throw new InvalidOperationException("The attempt's originating study revision is not retained.");
            if (_currentAttempt?.Id == attempt.Id ||
                _currentAttempt is not null && _currentExecution?.Standing is
                    StageExecutionStanding.Pending or StageExecutionStanding.Running or
                    StageExecutionStanding.ReadyForMinimization)
                throw new InvalidOperationException("A surviving unfinished attempt cannot be replaced by reattachment.");
            _currentAttempt = attempt;
            _currentExecution = new StageExecutionState(
                attempt.Id, null, null, StageExecutionStanding.Pending,
                "Preparation attempt admitted.", null, DateTimeOffset.UtcNow);
        }
    }

    public void RetainExecution(StageExecutionState execution)
    {
        ArgumentNullException.ThrowIfNull(execution);
        lock (_gate)
        {
            if (_currentAttempt?.Id != execution.AttemptId)
                throw new InvalidOperationException("Execution standing must identify the surviving attempt.");
            if (execution.Standing == StageExecutionStanding.Completed &&
                (execution.StageId is null || !_stages.TryGetValue(execution.StageId, out var completed) ||
                 completed.Attempt.Id != execution.AttemptId || completed.Kind != execution.Kind))
                throw new InvalidOperationException("Completion requires the exact retained completed stage.");
            if (_currentExecution is { } previous)
            {
                var previousTerminal = previous.Standing is StageExecutionStanding.Completed or
                    StageExecutionStanding.Stopped or StageExecutionStanding.Failed or
                    StageExecutionStanding.ResourceRefused or StageExecutionStanding.Unobserved;
                if (previousTerminal && execution.Standing is StageExecutionStanding.Pending or
                    StageExecutionStanding.Running)
                {
                    var newOptionalStage = execution.Kind == StageKind.Equilibration &&
                        _stages.Values.Any(stage => stage.Attempt.Id == execution.AttemptId &&
                            stage.Kind == StageKind.Minimization) &&
                        (previous.StageId != execution.StageId || previous.Kind != execution.Kind);
                    if (!newOptionalStage)
                        throw new InvalidOperationException("A terminal stage cannot return to pending progress.");
                }
                if (previous.Kind == StageKind.Equilibration && execution.Kind == StageKind.Minimization ||
                    previous.Kind == StageKind.Minimization && execution.Kind is null ||
                    previous.Standing == StageExecutionStanding.ReadyForMinimization &&
                        execution.Kind is null && execution.Standing == StageExecutionStanding.Running ||
                    previous.StageId is not null && execution.StageId is not null &&
                        previous.StageId != execution.StageId && previous.Standing is
                            StageExecutionStanding.Pending or StageExecutionStanding.Running)
                    throw new InvalidOperationException("Execution progress does not follow the retained stage.");
            }
            _currentExecution = execution;
        }
    }

    public void RetainCompletedStage(CompletedStage stage)
    {
        ArgumentNullException.ThrowIfNull(stage);
        lock (_gate)
        {
            if (_currentAttempt?.Id != stage.Attempt.Id)
                throw new InvalidOperationException("The completed stage does not belong to the current attempt.");
            if (_currentAttempt.StudyRevisionId != stage.Attempt.StudyRevisionId ||
                _currentAttempt.PolicyId != stage.Attempt.PolicyId ||
                _currentExecution?.Standing != StageExecutionStanding.Running ||
                _currentExecution.Kind != stage.Kind ||
                _currentExecution.StageId != stage.Id ||
                stage.Observation.StageId != stage.Id || stage.Observation.AttemptId != stage.Attempt.Id ||
                stage.Observation.Kind != stage.Kind ||
                stage.Molecule.Id != stage.Id || stage.Molecule.AtomCount <= 0 ||
                stage.Correspondence.ResultId != stage.Id || !stage.Correspondence.Complete ||
                stage.Correspondence.Atoms.IsDefault ||
                stage.Correspondence.Atoms.Length != stage.Molecule.AtomCount ||
                stage.Correspondence.Atoms.Where((atom, index) =>
                    atom.ResultAtomIndex != index || string.IsNullOrWhiteSpace(atom.ResultAtomId)).Any())
                throw new InvalidOperationException("Incomplete or mismatched material cannot be retained as a completed stage.");
            if (_stages.TryGetValue(stage.Id, out var earlier) && earlier != stage)
                throw new InvalidOperationException("A completed stage's identity cannot be relabelled.");
            _stages = _stages.SetItem(stage.Id, stage);
            _currentExecution = new StageExecutionState(
                stage.Attempt.Id, stage.Id, stage.Kind, StageExecutionStanding.Completed,
                "Completed stage established.", 1.0, stage.CompletedAt);
        }
    }

    public void RetainAssessment(PreparationAssessmentResult assessment)
    {
        ArgumentNullException.ThrowIfNull(assessment);
        lock (_gate)
        {
            if (!_stages.ContainsKey(assessment.StageId))
                throw new InvalidOperationException("Assessment cannot be attached to an unfinished or unknown stage.");
            _assessments = _assessments.SetItem(assessment.StageId, assessment);
        }
    }

    public void RetainFinding(ScientificFinding finding)
    {
        ArgumentNullException.ThrowIfNull(finding);
        lock (_gate)
            _findings = _findings.SetItem(finding.Id, finding);
    }

    public LocalWorkspaceSnapshot Snapshot()
    {
        lock (_gate)
            return new LocalWorkspaceSnapshot(
                _currentStudy,
                _currentAttempt,
                _currentExecution,
                _stages.Values.OrderBy(stage => stage.CompletedAt).ToImmutableArray(),
                _assessments.Values.OrderBy(assessment => assessment.AssessedAt).ToImmutableArray(),
                _findings.Values.OrderBy(finding => finding.EstablishedAt).ToImmutableArray(),
                DateTimeOffset.UtcNow);
    }
}
