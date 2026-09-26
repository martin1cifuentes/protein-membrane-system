using System.Collections.Immutable;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using AssessmentOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment.PreparationAssessment;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class PreparationAssessmentOwnerTests
{
    [Fact]
    public async Task Exact_stage_observations_meet_only_their_declared_stage_criteria()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var result = fixture.Assess();

        Assert.Equal(PreparationQualification.QualifiedPrepared, result.Qualification);
        Assert.True(result.CurrentlyApplicable);
        Assert.Equal(fixture.Stage.Id, result.StageId);
        Assert.Equal(StageKind.Minimization, fixture.Stage.Kind);
        var wrongStagePolicy = fixture.ForPolicy(fixture.Policy with { AssessmentCriteria =
            fixture.Policy.AssessmentCriteria.Select(item => item with
                { StageKind = StageKind.Equilibration }).ToImmutableArray() });
        Assert.Equal(PreparationQualification.Indeterminate, wrongStagePolicy.Assess().Qualification);
    }

    [Fact]
    public async Task Observed_stage_violation_is_not_qualified_while_missing_basis_is_indeterminate()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var violated = fixture.Stage with { Observation = fixture.Stage.Observation with
        { Measurements = EditMeasurement(fixture.Stage.Observation.Measurements,
            "minimumIntermolecularDistanceAngstrom", item => item with { Value = 0.5 }) } };
        Assert.Equal(PreparationQualification.NotQualified,
            fixture.Assess(stage: violated).Qualification);
        var missing = fixture.Stage with { Observation = fixture.Stage.Observation with
        { Measurements = fixture.Stage.Observation.Measurements.Remove(
            fixture.Stage.Observation.Measurements.Single(item =>
                item.Name == "minimumIntermolecularDistanceAngstrom")) } };
        Assert.Equal(PreparationQualification.Indeterminate,
            fixture.Assess(stage: missing).Qualification);
        var noPositiveCriteria = fixture.ForPolicy(fixture.Policy with
        { AssessmentCriteria = ImmutableArray<PreparationAssessmentCriterion>.Empty });
        Assert.Equal(PreparationQualification.Indeterminate,
            noPositiveCriteria.Assess().Qualification);
        Assert.True(noPositiveCriteria.Assess().CurrentlyApplicable);
        Assert.Equal(fixture.Stage.Id, fixture.Stage.Observation.StageId);
    }

    [Fact]
    public async Task Contact_relative_organization_and_protein_geometry_each_drive_their_own_negative_assessment()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var cases = new (string Name, Func<CompletedStage, CompletedStage> Change)[]
        {
            ("protein lipid contact", stage => stage with { Observation = stage.Observation with
                { LocalState = stage.Observation.LocalState! with { RolePairMeasurements =
                    ImmutableArray.Create(new LocalRolePairMeasurement(MoleculeRoleKind.Protein,
                        MoleculeRoleKind.Lipid, 1, 0.5)) } } }),
            ("relative leaflet order", stage => stage with { Observation = stage.Observation with
                { Measurements = EditMeasurement(stage.Observation.Measurements,
                    "leafletHeadSeparationAngstrom", item => item with { Value = 5 }) } }),
            ("relative protein offset", stage => stage with { Observation = stage.Observation with
                { Measurements = EditMeasurement(stage.Observation.Measurements,
                    "proteinBilayerMidplaneOffsetAngstrom", item => item with { Value = 30 }) } }),
            ("protein covalent geometry", stage => stage with { Observation = stage.Observation with
                { ProteinGeometry = stage.Observation.ProteinGeometry! with { Kinds =
                    stage.Observation.ProteinGeometry.Kinds.SetItem(0,
                        stage.Observation.ProteinGeometry.Kinds[0] with
                        { MinimumDistanceAngstrom = 0.5 }) } } }),
            ("protein chain continuity", stage => stage with { Observation = stage.Observation with
                { ProteinGeometry = stage.Observation.ProteinGeometry! with { Kinds =
                    stage.Observation.ProteinGeometry.Kinds.SetItem(1,
                        stage.Observation.ProteinGeometry.Kinds[1] with
                        { MaximumDistanceAngstrom = 3 }) } } }),
            ("protein nonbonded separation", stage => stage with { Observation = stage.Observation with
                { ProteinGeometry = stage.Observation.ProteinGeometry! with { Kinds =
                    stage.Observation.ProteinGeometry.Kinds.SetItem(2,
                        stage.Observation.ProteinGeometry.Kinds[2] with
                        { MinimumDistanceAngstrom = 0.5 }) } } })
        };
        foreach (var (name, change) in cases)
        {
            var result = fixture.Assess(stage: change(fixture.Stage));
            Assert.True(result.Qualification == PreparationQualification.NotQualified,
                name + ": observed defect was not interpreted negatively");
        }
    }

    [Fact]
    public async Task Missing_wrong_or_duplicate_stage_measurements_cannot_become_positive_evidence()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var cases = new (string Name, Func<CompletedStage, CompletedStage> Change)[]
        {
            ("wrong unit", stage => stage with { Observation = stage.Observation with
                { Measurements = EditMeasurement(stage.Observation.Measurements,
                    "leafletHeadSeparationAngstrom", item => item with { Unit = "nm" }) } }),
            ("wrong scope", stage => stage with { Observation = stage.Observation with
                { Measurements = EditMeasurement(stage.Observation.Measurements,
                    "leafletHeadSeparationAngstrom", item => item with { Scope = "absolute-z" }) } }),
            ("duplicate", stage => stage with { Observation = stage.Observation with
                { Measurements = stage.Observation.Measurements.Add(
                    stage.Observation.Measurements.Single(item =>
                        item.Name == "leafletHeadSeparationAngstrom")) } }),
            ("nonfinite", stage => stage with { Observation = stage.Observation with
                { Measurements = EditMeasurement(stage.Observation.Measurements,
                    "leafletHeadSeparationAngstrom", item => item with { Value = double.NaN }) } }),
            ("unavailable local state", stage => stage with { Observation = stage.Observation with
                { LocalState = stage.Observation.LocalState! with
                    { Standing = ObservationStanding.Unavailable } } }),
            ("global water ion contact without protein lipid pair", stage => stage with
            { Observation = stage.Observation with { LocalState = stage.Observation.LocalState! with
                { CoveredRolePairs = ImmutableArray.Create(new LocalContactRolePair(
                    MoleculeRoleKind.Water, MoleculeRoleKind.Ion)),
                  RolePairMeasurements = ImmutableArray.Create(new LocalRolePairMeasurement(
                    MoleculeRoleKind.Water, MoleculeRoleKind.Ion, 1, 2)) } } }),
            ("unavailable protein geometry", stage => stage with { Observation = stage.Observation with
                { ProteinGeometry = stage.Observation.ProteinGeometry! with
                    { Standing = ObservationStanding.Unavailable } } })
        };
        foreach (var (name, change) in cases)
        {
            var result = fixture.Assess(stage: change(fixture.Stage));
            Assert.True(result.Qualification == PreparationQualification.Indeterminate,
                name + ": incomplete evidence qualified the stage");
        }
    }

    [Fact]
    public async Task Absolute_leaflet_position_or_wrong_policy_and_placement_cannot_qualify_the_stage()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var absolute = fixture.ForPolicy(fixture.Policy with { AssessmentCriteria =
            ImmutableArray.Create(
                new PreparationAssessmentCriterion(StageKind.Minimization,
                    "upperLipidHeadMeanZAngstrom", "angstrom", "absolute-z", 0, 30),
                new PreparationAssessmentCriterion(StageKind.Minimization,
                    "lowerLipidHeadMeanZAngstrom", "angstrom", "absolute-z", -30, 0)) });
        var absoluteAssessment = absolute.Assess();
        Assert.Equal(PreparationQualification.Indeterminate, absoluteAssessment.Qualification);
        Assert.Contains("absolute", absoluteAssessment.Reason, StringComparison.OrdinalIgnoreCase);

        var owner = new AssessmentOwner();
        var wrongPlacement = owner.Assess(fixture.Stage, fixture.Constructed, fixture.Protein,
            fixture.Membrane, fixture.Placement with { Id = "other-placement" }, fixture.Policy,
            ImmutableArray<ScientificFinding>.Empty);
        Assert.False(wrongPlacement.CurrentlyApplicable);
        Assert.Equal(PreparationQualification.Indeterminate, wrongPlacement.Qualification);
        var wrongPolicy = fixture.Assess(policy: fixture.Policy with { Id = "other-policy" });
        Assert.False(wrongPolicy.CurrentlyApplicable);
        Assert.Equal(PreparationQualification.Indeterminate, wrongPolicy.Qualification);
    }

    [Fact]
    public async Task Wrong_stage_identity_or_observation_cannot_be_qualified()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        foreach (var changed in new[]
                 {
                     fixture.Stage with { Observation = fixture.Stage.Observation with
                         { StageId = "other" } },
                     fixture.Stage with { Observation = fixture.Stage.Observation with
                         { Kind = StageKind.Equilibration } },
                     fixture.Stage with { Observation = fixture.Stage.Observation with
                         { Termination = StageTermination.MaxIterations } },
                     fixture.Stage with { Correspondence = fixture.Stage.Correspondence with
                         { ResultId = "other" } },
                     fixture.Stage with { Correspondence = fixture.Stage.Correspondence with
                         { SourceId = "other" } },
                     fixture.Stage with { Correspondence = fixture.Stage.Correspondence with
                         { Atoms = fixture.Stage.Correspondence.Atoms.SetItem(0,
                             fixture.Stage.Correspondence.Atoms[0] with
                             { SourceAtomId = "other" }) } },
                     fixture.Stage with { Attempt = fixture.Stage.Attempt with
                         { StudyRevisionId = "other" } }
                 })
        {
            var assessment = fixture.Assess(stage: changed);
            Assert.Equal(PreparationQualification.Indeterminate, assessment.Qualification);
            Assert.False(assessment.CurrentlyApplicable);
        }
    }

    [Fact]
    public async Task Later_material_findings_change_current_assessment_without_rewriting_stage_history()
    {
        using var basis = new ConstructionFixture();
        var fixture = await AssessmentFixture.Create(basis);
        var original = fixture.Assess();
        var unrelated = Finding("other-attempt", FindingDisposition.Disqualifies, true);
        Assert.Equal(PreparationQualification.QualifiedPrepared,
            fixture.Assess(findings: ImmutableArray.Create(unrelated)).Qualification);
        var challenge = Finding(fixture.Placement.Proposal.Id, FindingDisposition.Challenges, true);
        var challenged = fixture.Assess(findings: ImmutableArray.Create(challenge));
        Assert.Equal(PreparationQualification.Indeterminate, challenged.Qualification);
        Assert.Contains(challenge, challenged.Findings);
        var contradiction = Finding(fixture.Constructed.Id, FindingDisposition.Disqualifies, true);
        Assert.Equal(PreparationQualification.NotQualified,
            fixture.Assess(findings: ImmutableArray.Create(contradiction)).Qualification);
        Assert.Equal(PreparationQualification.QualifiedPrepared, original.Qualification);
        Assert.Empty(fixture.Stage.Findings);
    }

    private static ScientificFinding Finding(string subject, FindingDisposition disposition, bool material) =>
        new(Guid.NewGuid().ToString("N"), subject, "evidence", "later observation",
            "review current interpretation", disposition, material, DateTimeOffset.UtcNow);

    private static ImmutableArray<MeasuredValue> EditMeasurement(
        ImmutableArray<MeasuredValue> measurements, string name,
        Func<MeasuredValue, MeasuredValue> edit)
    {
        var selected = measurements.Single(item => item.Name == name);
        return measurements.Replace(selected, edit(selected));
    }
}

internal sealed record AssessmentFixture(CompletedStage Stage, ConstructedExplicitSystem Constructed,
    AssessedPreparedProtein Protein, AssessedMembraneModel Membrane,
    AssessedProteinMembranePlacement Placement, ApplicablePreparationPolicy Policy)
{
    public static async Task<AssessmentFixture> Create(ConstructionFixture basis)
    {
        var worker = new ConstructionWorker(basis);
        var started = await new ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation(
            worker, worker, worker).StartAsync(basis.Attempt, basis.Revision, basis.Protein,
            basis.Membrane, basis.Placement, basis.Policy, basis.Directory, null, null,
            TestContext.Current.CancellationToken);
        var constructed = Assert.IsType<ConstructedExplicitSystem>(started.Constructed);
        var local = Assert.IsType<LocalStateObservations>(constructed.LocalState);
        var criteria = ImmutableArray.Create(
            new PreparationAssessmentCriterion(StageKind.Minimization, "leafletHeadSeparationAngstrom",
                "angstrom", "bilayer", 10, 40),
            new PreparationAssessmentCriterion(StageKind.Minimization, "proteinBilayerMidplaneOffsetAngstrom",
                "angstrom", "proteinVsBilayer", -20, 20),
            new PreparationAssessmentCriterion(StageKind.Minimization, "minimumIntermolecularDistanceAngstrom",
                "angstrom", "allMolecules", 1, 6));
        var policy = basis.Policy with { AssessmentCriteria = criteria };
        var attempt = basis.Attempt with
        { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy) };
        constructed = constructed with { Attempt = attempt };
        var geometryKinds = ImmutableArray.Create(
            new ProteinGeometryKindObservation("covalentBond", GeometryKindStanding.Observed,
                1, 1, 1.5, 1.5, null),
            new ProteinGeometryKindObservation("chainContinuity", GeometryKindStanding.Observed,
                1, 1, 1.5, 1.5, null),
            new ProteinGeometryKindObservation("nonbondedDistance", GeometryKindStanding.Observed,
                1, 1, 2, 2, null));
        var geometry = new ProteinGeometryObservations(ObservationStanding.Observed,
            geometryKinds, ImmutableArray<GeometryDistanceObservation>.Empty,
            ImmutableArray<string>.Empty);
        var stageId = "minimized-result";
        var observation = new StageObservation(stageId, attempt.Id, StageKind.Minimization,
            local.Measurements, ImmutableArray<ScientificEvidence>.Empty,
            StageTermination.Converged, "OpenMM 8.6.0", DateTimeOffset.UtcNow,
            LocalState: local, ProteinGeometry: geometry);
        var stage = new CompletedStage(stageId, attempt, StageKind.Minimization,
            constructed.Molecule with { Id = stageId, CoordinateSha256 = stageId }, observation,
            constructed.Correspondence with { ResultId = stageId }, policy.Id, null,
            ImmutableArray<ScientificFinding>.Empty, DateTimeOffset.UtcNow);
        return new AssessmentFixture(stage, constructed, basis.Protein, basis.Membrane,
            basis.Placement, policy);
    }

    public PreparationAssessmentResult Assess(CompletedStage? stage = null,
        ApplicablePreparationPolicy? policy = null,
        ImmutableArray<ScientificFinding> findings = default) =>
        new AssessmentOwner().Assess(stage ?? Stage, Constructed, Protein, Membrane,
            Placement, policy ?? Policy, findings);

    public AssessmentFixture ForPolicy(ApplicablePreparationPolicy policy)
    {
        var attempt = Stage.Attempt with
        { PolicyVersion = policy.Version,
            PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy),
            ForceFieldFiles = policy.ForceFieldFiles };
        return this with { Policy = policy, Stage = Stage with { Attempt = attempt },
            Constructed = Constructed with { Attempt = attempt } };
    }
}
