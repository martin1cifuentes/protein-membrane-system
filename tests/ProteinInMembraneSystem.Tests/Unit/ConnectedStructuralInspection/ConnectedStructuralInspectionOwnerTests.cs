using System.Collections.Immutable;
using System.Text.Json;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using InspectionOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ConnectedStructuralInspection.ConnectedStructuralInspection;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class ConnectedStructuralInspectionOwnerTests
{
    [Fact]
    public void Geometry_stays_with_its_exact_subject_and_supporting_evidence()
    {
        var inspection = new InspectionOwner();
        var revision = Revision("current", 1);
        var residue = new ResidueAddress(0, "A", 42, "", "A");
        var geometry = new ProteinGeometryObservations(ObservationStanding.Observed,
            ImmutableArray.Create(new ProteinGeometryKindObservation("covalentBond",
                GeometryKindStanding.Observed, 1, 1, 1.4999999999, 1.4999999999, null)),
            ImmutableArray.Create(new GeometryDistanceObservation("covalentBond",
                new AtomAddress(residue, "C"), new AtomAddress(residue with { Residue = 43 }, "N"),
                1.4999999999, 1.5)),
            ImmutableArray.Create("One identified pair was measured."));
        var subject = Subject("protein", revision.Id, "geometry-evidence", "A:42", 42) with
        { Geometry = new InspectionGeometryAccount(geometry, ImmutableArray.Create("geometry-evidence")) };

        var selected = inspection.Select(revision, subject);

        Assert.True(selected.Established, selected.Reason);
        Assert.Same(geometry, selected.Value!.Geometry?.Observations);
        Assert.Equal("geometry-evidence", Assert.Single(selected.Value.Geometry!.EvidenceIds));
        using var json = JsonDocument.Parse(JsonSerializer.Serialize(selected.Value,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        Assert.Equal(1.4999999999, json.RootElement.GetProperty("geometry")
            .GetProperty("observations").GetProperty("locatedDistances")[0]
            .GetProperty("distanceAngstrom").GetDouble());
        var invalid = inspection.Select(revision, subject with
        { Geometry = subject.Geometry with { EvidenceIds = ImmutableArray.Create("other-evidence") } });
        Assert.False(invalid.Established);
        Assert.Equal("protein", inspection.Current?.SubjectId);
        Assert.Same(geometry, inspection.Current?.Geometry?.Observations);
    }

    [Fact]
    public void Selection_keeps_the_exact_origin_evidence_finding_assessment_and_focused_part()
    {
        var inspection = new InspectionOwner();
        var firstRevision = Revision("revision-one", 1);
        var secondRevision = Revision("revision-two", 2);
        var first = Subject("stage-one", firstRevision.Id, "evidence-one", "A:42", 42,
            assessment: Assessment("stage-one")) with
        { OmittedMolecules = ImmutableArray.Create("water molecules absent from this preview") };
        var second = Subject("stage-two", secondRevision.Id, "evidence-two", "B:7", 7);

        var selected = inspection.Select(firstRevision, first);
        Assert.True(selected.Established, selected.Reason);
        Assert.Equal(firstRevision.Id, selected.Value!.StudyRevisionId);
        Assert.Equal(firstRevision.Number, selected.Value.StudyRevisionNumber);
        Assert.Equal("stage-one", selected.Value.SubjectId);
        Assert.Equal("water molecules absent from this preview", Assert.Single(selected.Value.OmittedMolecules));
        Assert.Equal("evidence-one", Assert.Single(selected.Value.Evidence).Id);
        Assert.Equal("finding-evidence-one", Assert.Single(selected.Value.Findings).Id);
        Assert.Equal("stage-one", selected.Value.Assessment?.StageId);
        Assert.Equal("Å", Assert.Single(selected.Value.Metrics).Unit);

        var focused = inspection.Focus("stage-one", first.Annotations[0].Id);
        Assert.True(focused.Established, focused.Reason);
        Assert.Equal("A:42", focused.Value!.FocusId);
        Assert.Equal(42, focused.Value.Focus?.AuthSeqId);
        Assert.Same(selected.Value.Assessment, focused.Value.Assessment);
        Assert.Equal(firstRevision.Id, focused.Value.StudyRevisionId);
        Assert.Equal(firstRevision.Id, first.StudyRevisionId);

        var selectedSecond = inspection.Select(secondRevision, second);
        Assert.True(selectedSecond.Established, selectedSecond.Reason);
        Assert.Equal("stage-two", selectedSecond.Value!.SubjectId);
        Assert.Null(selectedSecond.Value.FocusId);
        Assert.Null(selectedSecond.Value.Assessment);
        Assert.DoesNotContain(selectedSecond.Value.Evidence, item => item.Id == "evidence-one");
        Assert.False(inspection.Focus("stage-one", first.Annotations[0].Id).Established);
        Assert.False(inspection.Focus("stage-two", first.Annotations[0].Id).Established);
        Assert.Equal("stage-two", inspection.Current?.SubjectId);
    }

    [Fact]
    public void Foreign_revision_evidence_finding_assessment_or_link_cannot_replace_selection()
    {
        var inspection = new InspectionOwner();
        var revision = Revision("current", 3);
        var valid = Subject("proposal", revision.Id, "evidence", "A:10", 10);
        Assert.True(inspection.Select(revision, valid).Established);
        var foreignEvidence = Evidence("foreign", "other-subject");
        var cases = new[]
        {
            valid with { StudyRevisionId = "older" },
            valid with { Evidence = ImmutableArray.Create(foreignEvidence) },
            valid with { Findings = ImmutableArray.Create(Finding("finding", "other-subject", "evidence")) },
            valid with { Findings = ImmutableArray.Create(Finding("finding", valid.Id, "missing-evidence")) },
            valid with { Assessment = Assessment("other-stage") },
            valid with { Annotations = ImmutableArray.Create(valid.Annotations[0] with { EvidenceId = "missing-evidence" }) },
            valid with { Metrics = ImmutableArray.Create(valid.Metrics[0] with { EvidenceId = "missing-evidence" }) },
        };

        foreach (var invalid in cases)
        {
            var refused = inspection.Select(revision, invalid);
            Assert.False(refused.Established);
            Assert.Equal("proposal", inspection.Current?.SubjectId);
            Assert.Equal(revision.Id, inspection.Current?.StudyRevisionId);
        }
    }

    [Fact]
    public void An_annotation_without_verified_structure_and_location_cannot_claim_spatial_focus()
    {
        var inspection = new InspectionOwner();
        var revision = Revision("current", 5);
        var located = Subject("subject", revision.Id, "evidence", "A:10", 10);
        var unlocated = located with
        { Annotations = ImmutableArray.Create(located.Annotations[0] with { GeometryFocus = null }) };
        Assert.True(inspection.Select(revision, unlocated).Established);
        Assert.False(inspection.Focus(located.Id, unlocated.Annotations[0].Id).Established);
        Assert.Null(inspection.Current?.FocusId);

        var unrenderable = located with { StructureUrl = null };
        Assert.True(inspection.Select(revision, unrenderable).Established);
        Assert.False(inspection.Focus(located.Id, located.Annotations[0].Id).Established);
        Assert.Null(inspection.Current?.FocusId);
        Assert.Equal("evidence", Assert.Single(inspection.Current!.Evidence).Id);
    }

    private static StudyRevision Revision(string id, long number) =>
        new(id, number, null, null, null, FixedStudyConditions.Initial);

    private static ScientificEvidence Evidence(string id, string subjectId) =>
        new(id, subjectId, "observed source", "measured distance", "2.4 Å", "one selected residue",
            "bounded coordinate precision", EvidenceBearing.Context);

    private static ScientificFinding Finding(string id, string subjectId, string evidenceId) =>
        new(id, subjectId, evidenceId, "located observation", "review this region",
            FindingDisposition.Context, false, DateTimeOffset.UtcNow);

    private static PreparationAssessmentResult Assessment(string stageId) =>
        new("assessment-" + stageId, stageId, PreparationCheckStanding.ChecksIncomplete,
            "An exact stage assessment.", ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<ScientificFinding>.Empty, ImmutableArray.Create("Scope is limited."),
            DateTimeOffset.UtcNow, true);

    private static InspectionSubject Subject(string id, string revisionId, string evidenceId,
        string partId, int sequence, PreparationAssessmentResult? assessment = null)
    {
        var evidence = Evidence(evidenceId, id);
        return new InspectionSubject(id, revisionId, "/api/structures/exact?format=mmcif", "completedStage",
            ImmutableArray<string>.Empty, ImmutableArray.Create(evidence),
            ImmutableArray.Create(Finding("finding-" + evidenceId, id, evidenceId)),
            ImmutableArray.Create(new InspectionAnnotation("annotation-" + evidenceId, partId,
                "located residue", "exact selected part", evidenceId,
                new StructureFocus(partId[..1], sequence, null, null))),
            ImmutableArray.Create(new InspectionMetric("separation", "2.4", "Å", partId, evidenceId)),
            assessment);
    }
}
