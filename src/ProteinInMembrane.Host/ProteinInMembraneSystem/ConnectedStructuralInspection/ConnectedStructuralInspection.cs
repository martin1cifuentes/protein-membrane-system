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
        if (subject.Evidence.Any(evidence => evidence.SubjectId != subject.Id) ||
            subject.Findings.Any(finding => finding.SubjectId != subject.Id) ||
            subject.Assessment is not null && subject.Assessment.StageId != subject.Id)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The supplied evidence, finding, or assessment does not identify this exact subject.");

        var evidenceIds = subject.Evidence.Select(evidence => evidence.Id).ToHashSet(StringComparer.Ordinal);
        if (subject.Annotations.Any(annotation => annotation.EvidenceId is not null && !evidenceIds.Contains(annotation.EvidenceId)) ||
            subject.Metrics.Any(metric => metric.EvidenceId is not null && !evidenceIds.Contains(metric.EvidenceId)))
            return BoundaryOutcome<InspectionAccount>.Unavailable("A spatial mark or numerical value lacks corresponding evidence for this subject.");

        _subject = subject;
        _account = new InspectionAccount(
            subject.Id, subject.StructureUrl, subject.RepresentationKind,
            subject.OmittedMolecules.IsDefault ? ImmutableArray<string>.Empty : subject.OmittedMolecules,
            null, null,
            subject.Annotations.IsDefault ? ImmutableArray<InspectionAnnotation>.Empty : subject.Annotations,
            subject.Metrics.IsDefault ? ImmutableArray<InspectionMetric>.Empty : subject.Metrics);
        return BoundaryOutcome<InspectionAccount>.Success(_account);
    }

    public BoundaryOutcome<InspectionAccount> Focus(string subjectId, string annotationId)
    {
        if (_subject?.Id != subjectId || _account is null)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The requested inspection subject is no longer selected.");
        var annotation = _subject.Annotations.FirstOrDefault(item => item.Id == annotationId);
        if (annotation is null)
            return BoundaryOutcome<InspectionAccount>.Unavailable("The requested part does not belong to this subject's account.");
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
        ImmutableArray<string> requiredEvidenceIds,
        out string reason)
    {
        if (_subject?.Id != subjectId || _account is null)
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
            if (!_subject.Evidence.Any(evidence => evidence.Id == id) ||
                !_account.Annotations.Any(annotation => annotation.EvidenceId == id) ||
                !_account.Metrics.Any(metric => metric.EvidenceId == id))
            {
                reason = $"Evidence {id} is not connected in both the spatial and numerical account.";
                return false;
            }
        }
        reason = string.Empty;
        return true;
    }
}
