using System.Collections.Immutable;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using WorkspaceOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class LocalRunWorkspaceOwnerTests
{
    [Fact]
    public void Reopened_snapshot_keeps_the_surviving_attempt_and_original_revision()
    {
        var workspace = new WorkspaceOwner();
        var origin = Revision("origin", 1);
        var current = Revision("current", 2);
        var attempt = Attempt("attempt", origin.Id);
        workspace.RetainStudy(origin);
        workspace.RetainAttempt(attempt);
        var pending = workspace.Snapshot();
        Assert.Equal(origin.Id, pending.CurrentStudy?.Id);
        Assert.Equal(attempt.Id, pending.CurrentAttempt?.Id);
        Assert.Equal(StageExecutionStanding.Pending, pending.CurrentExecution?.Standing);
        Assert.Empty(pending.CompletedStages);

        var running = State(attempt.Id, null, null, StageExecutionStanding.Running, 0.4);
        workspace.RetainExecution(running);
        workspace.RetainStudy(current);
        var reopened = workspace.Snapshot();
        Assert.Equal(current.Id, reopened.CurrentStudy?.Id);
        Assert.Equal(origin.Id, reopened.CurrentAttempt?.StudyRevisionId);
        Assert.Equal(running, reopened.CurrentExecution);
        Assert.Empty(reopened.CompletedStages);
        Assert.Throws<InvalidOperationException>(() => workspace.RetainAttempt(Attempt("replacement", current.Id)));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainExecution(
            State("foreign-attempt", null, null, StageExecutionStanding.Completed, 1)));
        Assert.Equal(running, workspace.Snapshot().CurrentExecution);

        var candidate = State(attempt.Id, null, null, StageExecutionStanding.ReadyForMinimization, 1);
        workspace.RetainExecution(candidate);
        Assert.Equal(StageExecutionStanding.ReadyForMinimization, workspace.Snapshot().CurrentExecution?.Standing);
        Assert.Throws<InvalidOperationException>(() => workspace.RetainAttempt(Attempt("replacement", current.Id)));
        Assert.Equal(attempt.Id, workspace.Snapshot().CurrentAttempt?.Id);
    }

    [Fact]
    public void Admission_after_study_advances_keeps_the_known_attempt_origin()
    {
        var workspace = new WorkspaceOwner();
        var origin = Revision("origin", 1);
        workspace.RetainStudy(origin);
        workspace.RetainStudy(Revision("current", 2));
        workspace.RetainAttempt(Attempt("admitted-after-revision", origin.Id));

        var reopened = workspace.Snapshot();
        Assert.Equal("current", reopened.CurrentStudy?.Id);
        Assert.Equal(origin.Id, reopened.CurrentAttempt?.StudyRevisionId);
        Assert.Equal(StageExecutionStanding.Pending, reopened.CurrentExecution?.Standing);
        Assert.Throws<InvalidOperationException>(() =>
            workspace.RetainAttempt(Attempt("unknown-origin", "never-retained")));
        Assert.Equal("admitted-after-revision", workspace.Snapshot().CurrentAttempt?.Id);
    }

    [Theory]
    [InlineData(StageExecutionStanding.Stopped)]
    [InlineData(StageExecutionStanding.Failed)]
    [InlineData(StageExecutionStanding.ResourceRefused)]
    [InlineData(StageExecutionStanding.Unobserved)]
    public void Terminal_work_without_a_stage_cannot_regress_or_appear_completed(StageExecutionStanding standing)
    {
        var workspace = Started("attempt");
        workspace.RetainExecution(State("attempt", null, null, StageExecutionStanding.Running, 0.2));
        workspace.RetainExecution(State("attempt", null, null, standing, null));

        Assert.Throws<InvalidOperationException>(() => workspace.RetainExecution(
            State("attempt", null, null, StageExecutionStanding.Running, 0.8)));
        var reopened = workspace.Snapshot();
        Assert.Equal(standing, reopened.CurrentExecution?.Standing);
        Assert.Empty(reopened.CompletedStages);
    }

    [Fact]
    public void Only_matching_complete_stage_can_be_retained_and_later_work_keeps_it_separate()
    {
        var workspace = Started("attempt");
        var attempt = workspace.Snapshot().CurrentAttempt!;
        var running = State(attempt.Id, "minimized", StageKind.Minimization, StageExecutionStanding.Running, 0.8);
        workspace.RetainExecution(running);
        var stage = Stage(attempt);

        Assert.Throws<InvalidOperationException>(() => workspace.RetainExecution(
            State(attempt.Id, stage.Id, stage.Kind, StageExecutionStanding.Completed, 1)));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainCompletedStage(
            stage with { Observation = stage.Observation with { AttemptId = "foreign-attempt" } }));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainCompletedStage(
            stage with { Observation = stage.Observation with { StageId = "foreign-stage" } }));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainCompletedStage(
            stage with { Correspondence = stage.Correspondence with { ResultId = "foreign-result" } }));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainCompletedStage(
            stage with { Correspondence = stage.Correspondence with { Complete = false } }));
        Assert.Throws<InvalidOperationException>(() => workspace.RetainCompletedStage(
            stage with { Attempt = attempt with { StudyRevisionId = "foreign-revision" } }));
        Assert.Empty(workspace.Snapshot().CompletedStages);
        Assert.Equal(running, workspace.Snapshot().CurrentExecution);

        workspace.RetainCompletedStage(stage);
        var assessment = Assessment(stage.Id);
        workspace.RetainAssessment(assessment);
        var completed = workspace.Snapshot();
        Assert.Equal(stage.Id, Assert.Single(completed.CompletedStages).Id);
        Assert.Equal(assessment.Id, Assert.Single(completed.Assessments).Id);
        Assert.Equal(StageExecutionStanding.Completed, completed.CurrentExecution?.Standing);
        Assert.Throws<InvalidOperationException>(() => workspace.RetainExecution(running));

        workspace.RetainExecution(State(attempt.Id, null, StageKind.Equilibration,
            StageExecutionStanding.Pending, null));
        workspace.RetainExecution(State(attempt.Id, "optional", StageKind.Equilibration,
            StageExecutionStanding.Stopped, null));
        var afterOptionalStop = workspace.Snapshot();
        Assert.Equal(StageExecutionStanding.Stopped, afterOptionalStop.CurrentExecution?.Standing);
        Assert.Equal(stage.Id, Assert.Single(afterOptionalStop.CompletedStages).Id);
        Assert.Equal(assessment.Id, Assert.Single(afterOptionalStop.Assessments).Id);
    }

    private static WorkspaceOwner Started(string attemptId)
    {
        var workspace = new WorkspaceOwner();
        workspace.RetainStudy(Revision("origin", 1));
        workspace.RetainAttempt(Attempt(attemptId, "origin"));
        return workspace;
    }

    private static StudyRevision Revision(string id, long number) =>
        new(id, number, null, null, null, FixedStudyConditions.Initial);

    private static PreparationAttempt Attempt(string id, string revisionId) =>
        new(id, revisionId, "protein", "membrane", "placement", "policy",
            DateTimeOffset.UtcNow, "v1", "fingerprint", ImmutableArray<ForceFieldAsset>.Empty,
            "OpenMM test", "patch-hash");

    private static StageExecutionState State(string attemptId, string? stageId, StageKind? kind,
        StageExecutionStanding standing, double? progress) =>
        new(attemptId, stageId, kind, standing, standing.ToString(), progress, DateTimeOffset.UtcNow);

    private static CompletedStage Stage(PreparationAttempt attempt)
    {
        const string id = "minimized";
        var molecule = new MolecularArtifact(id, "stage.cif", "digest", null, null, null,
            1, "periodic cell");
        var observation = new StageObservation(id, attempt.Id, StageKind.Minimization,
            ImmutableArray<MeasuredValue>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            StageTermination.Converged, "OpenMM test", DateTimeOffset.UtcNow);
        var atom = new AtomCorrespondence(0, "result-atom", "source-atom", AtomOriginKind.Source,
            MoleculeRoleKind.Protein, AtomRoleKind.Backbone, "C", null, null);
        var correspondence = new SourceToResultCorrespondence("source", id,
            ImmutableArray.Create(atom), true);
        return new CompletedStage(id, attempt, StageKind.Minimization, molecule, observation,
            correspondence, attempt.PolicyId, null, ImmutableArray<ScientificFinding>.Empty,
            DateTimeOffset.UtcNow);
    }

    private static PreparationAssessmentResult Assessment(string stageId) =>
        new("assessment", stageId, PreparationQualification.Indeterminate,
            "One exact completed stage assessment.", ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<ScientificFinding>.Empty, ImmutableArray<string>.Empty,
            DateTimeOffset.UtcNow, true);
}
