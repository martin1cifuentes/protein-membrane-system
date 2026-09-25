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
            _currentStudy = study;
        }
    }

    public void RetainAttempt(PreparationAttempt attempt)
    {
        ArgumentNullException.ThrowIfNull(attempt);
        lock (_gate)
        {
            if (_currentStudy?.Id != attempt.StudyRevisionId)
                throw new InvalidOperationException("The attempt does not belong to the current study revision.");
            if (_currentAttempt is not null && _currentExecution?.Standing is StageExecutionStanding.Pending or StageExecutionStanding.Running)
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
            if (_currentExecution?.StageId is not null && execution.StageId == _currentExecution.StageId &&
                _currentExecution.Standing is StageExecutionStanding.Completed or StageExecutionStanding.Stopped or
                    StageExecutionStanding.Failed or StageExecutionStanding.ResourceRefused or StageExecutionStanding.Unobserved &&
                execution.Standing is StageExecutionStanding.Pending or StageExecutionStanding.Running)
                throw new InvalidOperationException("A terminal stage cannot return to pending progress.");
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
            if (stage.Observation.StageId != stage.Id || stage.Observation.AttemptId != stage.Attempt.Id ||
                stage.Correspondence.ResultId != stage.Molecule.Id || !stage.Correspondence.Complete)
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
