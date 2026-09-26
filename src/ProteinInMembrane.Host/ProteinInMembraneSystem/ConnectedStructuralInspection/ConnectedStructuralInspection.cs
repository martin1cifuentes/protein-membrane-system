using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ConnectedStructuralInspection;

/// <summary>
/// Keeps the spatial and numerical account on one exact subject. Presentation
/// and viewpoint choices cannot change a scientific result or study revision.
/// </summary>
public sealed class ConnectedStructuralInspection
{
    private InspectionSubject? _subject;
    private InspectionAccount? _account;

    public InspectionAccount? Current => _account;

    /// <summary>Invalidates a selection when its scientific study context changes.</summary>
    public void Clear()
    {
        _subject = null;
        _account = null;
    }

    public BoundaryOutcome<InspectionAccount> Select(StudyRevision revision, InspectionSubject subject)
    {
        ArgumentNullException.ThrowIfNull(revision);
        ArgumentNullException.ThrowIfNull(subject);
        if (revision.Id != subject.StudyRevisionId)
            return BoundaryOutcome<InspectionAccount>.Unavailable("This inspection subject belongs to a different study revision.");
        var evidence = subject.Evidence.IsDefault ? ImmutableArray<ScientificEvidence>.Empty : subject.Evidence;
        var findings = subject.Findings.IsDefault ? ImmutableArray<ScientificFinding>.Empty : subject.Findings;
        var annotations = subject.Annotations.IsDefault
            ? ImmutableArray<InspectionAnnotation>.Empty : subject.Annotations;
        var metrics = subject.Metrics.IsDefault ? ImmutableArray<InspectionMetric>.Empty : subject.Metrics;
        if (evidence.Any(item => item.SubjectId != subject.Id || string.IsNullOrWhiteSpace(item.Id)) ||
            evidence.Select(item => item.Id).Distinct(StringComparer.Ordinal).Count() != evidence.Length ||
            findings.Any(item => item.SubjectId != subject.Id || string.IsNullOrWhiteSpace(item.Id)) ||
            subject.Assessment is not null && subject.Assessment.StageId != subject.Id)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The supplied evidence, finding, or assessment does not identify this exact subject.");

        var evidenceIds = evidence.Select(item => item.Id).ToHashSet(StringComparer.Ordinal);
        if (findings.Any(item => !evidenceIds.Contains(item.EvidenceId)) ||
            annotations.Any(item => item.EvidenceId is not null && !evidenceIds.Contains(item.EvidenceId)) ||
            metrics.Any(item => item.EvidenceId is not null && !evidenceIds.Contains(item.EvidenceId)) ||
            subject.Assessment is { } assessment && !assessment.Evidence.IsDefault &&
                assessment.Evidence.Any(item => item.SubjectId != subject.Id))
            return BoundaryOutcome<InspectionAccount>.Unavailable("A spatial mark or numerical value lacks corresponding evidence for this subject.");

        _subject = subject;
        _account = new InspectionAccount(
            subject.Id, subject.StructureUrl, subject.RepresentationKind,
            subject.OmittedMolecules.IsDefault ? ImmutableArray<string>.Empty : subject.OmittedMolecules,
            null, null, annotations, metrics, revision.Id, revision.Number,
            evidence, findings, subject.Assessment);
        return BoundaryOutcome<InspectionAccount>.Success(_account);
    }

    public BoundaryOutcome<InspectionAccount> Focus(string subjectId, string annotationId)
    {
        if (_subject?.Id != subjectId || _account is null)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The requested inspection subject is no longer selected.");
        var annotation = _account.Annotations.FirstOrDefault(item => item.Id == annotationId);
        if (annotation is null)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The requested part does not belong to this subject's account.");
        if (annotation.GeometryFocus is null || string.IsNullOrWhiteSpace(_account.StructureUrl))
            return BoundaryOutcome<InspectionAccount>.Unavailable("This part has no verified spatial focus in the selected structure.");
        _account = _account with { FocusId = annotation.SubjectPartId, Focus = annotation.GeometryFocus };
        return BoundaryOutcome<InspectionAccount>.Success(_account);
    }

    public BoundaryOutcome<InspectionAccount> ClearFocus(string subjectId)
    {
        if (_subject?.Id != subjectId || _account is null)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The requested inspection subject is no longer selected.");
        _account = _account with { FocusId = null, Focus = null };
        return BoundaryOutcome<InspectionAccount>.Success(_account);
    }

    public bool HasRequiredEvidenceForApproval(
        string subjectId,
        string expectedStudyRevisionId,
        ImmutableArray<string> requiredEvidenceIds,
        out string reason)
    {
        if (_subject?.Id != subjectId || _account is null ||
            _account.StudyRevisionId != expectedStudyRevisionId)
        {
            reason = "The exact current proposal has not been selected for inspection.";
            return false;
        }
        if (requiredEvidenceIds.IsDefaultOrEmpty || string.IsNullOrWhiteSpace(_account.StructureUrl))
        {
            reason = "Required spatial and numerical evidence is not defined or available for this proposal.";
            return false;
        }
        foreach (var id in requiredEvidenceIds)
        {
            if (!_account.Evidence.Any(item => item.Id == id) ||
                !_account.Annotations.Any(annotation => annotation.EvidenceId == id &&
                    annotation.GeometryFocus is not null &&
                    _account.Metrics.Any(metric => metric.EvidenceId == id &&
                        metric.SubjectPartId == annotation.SubjectPartId &&
                        !string.IsNullOrWhiteSpace(metric.Value))))
            {
                reason = $"Evidence {id} is not connected in both the spatial and numerical account.";
                return false;
            }
        }
        reason = string.Empty;
        return true;
    }
}
