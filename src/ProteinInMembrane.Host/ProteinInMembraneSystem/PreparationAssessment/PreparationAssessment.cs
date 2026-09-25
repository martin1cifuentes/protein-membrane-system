using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment;

/// <summary>Judges an exact completed stage under its declared scientific policy.</summary>
public sealed class PreparationAssessment
{
    public PreparationAssessmentResult Assess(
        CompletedStage stage,
        ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein,
        AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy,
        ImmutableArray<ScientificFinding> laterFindings)
    {
        ArgumentNullException.ThrowIfNull(stage);
        ArgumentNullException.ThrowIfNull(constructed);
        ArgumentNullException.ThrowIfNull(policy);

        var findings = stage.Findings.AddRange(laterFindings.IsDefault
            ? ImmutableArray<ScientificFinding>.Empty : laterFindings);
        var evidence = stage.Observation.Evidence;
        // A later finding can invalidate this stage through the exact input it used,
        // even when the finding is not addressed to the stage itself. Do not apply
        // findings about other attempts or historical study subjects.
        var applicableSubjects = new HashSet<string>(StringComparer.Ordinal)
        {
            stage.Id, stage.Attempt.Id, constructed.Id, protein.Id, membrane.Id,
            membrane.Intended.Id, placement.Id, placement.Proposal.Id
        };
        PreparationQualification qualification;
        string reason;
        var currentlyApplicable = true;
        var proteinGeometry = AssessStageProteinGeometry(stage, policy);

        if (stage.Attempt.Id != constructed.Attempt.Id || stage.PolicyId != policy.Id ||
            stage.Attempt.ProteinId != protein.Id || stage.Attempt.MembraneId != membrane.Id ||
            stage.Attempt.PlacementId != placement.Id ||
            placement.Standing != AssessmentStanding.Supported ||
            !stage.Correspondence.Complete || stage.Correspondence.ResultId != stage.Molecule.Id ||
            stage.Observation.StageId != stage.Id || stage.Observation.AttemptId != stage.Attempt.Id)
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "The completed stage no longer has an applicable, corresponding scientific basis.";
            currentlyApplicable = false;
        }
        else if (findings.Any(finding => finding.Material && applicableSubjects.Contains(finding.SubjectId) &&
                     finding.Disposition == FindingDisposition.Disqualifies))
        {
            qualification = PreparationQualification.NotQualified;
            reason = "A material finding about this stage or one of its exact inputs disqualifies this prepared result.";
        }
        else if (findings.Any(finding => finding.Material && applicableSubjects.Contains(finding.SubjectId) &&
                     finding.Disposition == FindingDisposition.Challenges))
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "A material finding challenges this stage or one of its exact inputs; required interpretation is unresolved.";
        }
        else if (proteinGeometry == PreparationQualification.NotQualified)
        {
            qualification = PreparationQualification.NotQualified;
            reason = "The completed stage's measured protein intramolecular geometry falls outside its applicable stage-specific condition.";
        }
        else if (stage.Kind == StageKind.Equilibration &&
                 (stage.Observation.ObservationAdequacy != EquilibrationObservationAdequacy.Adequate ||
                  stage.Observation.EquilibrationAssessments.IsDefaultOrEmpty ||
                  stage.Observation.EquilibrationAssessments.Any(item => !item.Sufficient)))
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "The completed optional procedure did not establish sufficient, policy-declared time-series observations; its actual state remains reviewable and exportable.";
        }
        else if (!LocalStateObserved(stage.Observation.LocalState, policy.LocalStateObservation))
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "The completed stage has no complete, attributable local contact and protein–bilayer observation under the applicable policy.";
        }
        else if (!ContactRulesAvailable(stage.Kind, policy))
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "No applicable exact protein–lipid contact interpretation was declared before this stage.";
        }
        else if (!ContactRulesMet(stage.Kind, stage.Observation.LocalState!, policy))
        {
            qualification = PreparationQualification.NotQualified;
            reason = "An observed protein–lipid or other declared molecular contact falls outside the stage's qualified contact condition.";
        }
        else if (proteinGeometry != PreparationQualification.QualifiedPrepared)
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "The completed stage lacks corresponding measured protein geometry or applicable stage-specific integrity criteria.";
        }
        else if (string.IsNullOrWhiteSpace(policy.Version) || policy.EvidenceReferences.IsDefaultOrEmpty ||
                 policy.AssessmentCriteria.IsDefaultOrEmpty)
        {
            qualification = PreparationQualification.Indeterminate;
            reason = "No versioned, evidence-backed qualification criteria are declared for this stage.";
        }
        else
        {
            var criteria = policy.AssessmentCriteria.Where(criterion => criterion.StageKind == stage.Kind).ToImmutableArray();
            if (criteria.IsDefaultOrEmpty)
            {
                qualification = PreparationQualification.Indeterminate;
                reason = "The policy has no criteria for this kind of completed stage.";
            }
            else
            {
                if (!HasOrganizationCriterion(criteria, "leafletHeadSeparationAngstrom",
                        "bilayer", positiveMinimum: true) ||
                    !HasOrganizationCriterion(criteria, "proteinBilayerMidplaneOffsetAngstrom",
                        "proteinVsBilayer", positiveMinimum: false) ||
                    criteria.Any(criterion => criterion.MeasurementName is
                        "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom"))
                {
                    qualification = PreparationQualification.Indeterminate;
                    reason = "The stage policy lacks qualified relative leaflet-order and protein–bilayer organization bounds, or relies on absolute leaflet z positions.";
                    return new PreparationAssessmentResult(Guid.NewGuid().ToString("N"), stage.Id,
                        qualification, reason, evidence, findings,
                        policy.Limitations.IsDefault ? ImmutableArray<string>.Empty : policy.Limitations,
                        DateTimeOffset.UtcNow, currentlyApplicable);
                }
                var requiredValues = criteria.Select(criterion => new
                {
                    Criterion = criterion,
                    Value = stage.Observation.Measurements.FirstOrDefault(value =>
                        value.Name == criterion.MeasurementName && value.Unit == criterion.Unit &&
                        value.Scope == criterion.Scope)
                }).ToImmutableArray();
                if (requiredValues.Any(entry => entry.Value is null || !double.IsFinite(entry.Value.Value)))
                {
                    qualification = PreparationQualification.Indeterminate;
                    reason = "A required stage observation is missing, mismatched in unit/scope, or non-finite.";
                }
                else if (requiredValues.Any(entry =>
                             entry.Criterion.Minimum.HasValue && entry.Value!.Value < entry.Criterion.Minimum.Value ||
                             entry.Criterion.Maximum.HasValue && entry.Value!.Value > entry.Criterion.Maximum.Value))
                {
                    qualification = PreparationQualification.NotQualified;
                    reason = "An observed stage value falls outside the declared qualification condition.";
                }
                else
                {
                    qualification = PreparationQualification.QualifiedPrepared;
                    reason = "All applicable declared stage criteria are supported by corresponding observations.";
                }
            }
        }

        return new PreparationAssessmentResult(
            Guid.NewGuid().ToString("N"), stage.Id, qualification, reason,
            evidence, findings, policy.Limitations.IsDefault ? ImmutableArray<string>.Empty : policy.Limitations,
            DateTimeOffset.UtcNow, currentlyApplicable);
    }

    private static bool LocalStateObserved(LocalStateObservations? local,
        LocalStateObservationSpec? spec)
    {
        if (local is null || spec is null || local.Standing != ObservationStanding.Observed ||
            local.Measurements.IsDefault || local.LocatedContacts.IsDefault ||
            local.CoveredRolePairs.IsDefault || local.RolePairMeasurements.IsDefault ||
            local.RolePairMeasurements.Any(item => !Enum.IsDefined(item.FirstMoleculeRole) ||
                !Enum.IsDefined(item.SecondMoleculeRole)) ||
            spec.ContactRolePairs.IsDefaultOrEmpty ||
            spec.RequiredMetricNames.IsDefaultOrEmpty ||
            !spec.RequiredMetricNames.Contains("leafletHeadSeparationAngstrom") ||
            !spec.RequiredMetricNames.Contains("proteinBilayerMidplaneOffsetAngstrom") ||
            spec.ContactRolePairs.Any(pair => !local.CoveredRolePairs.Contains(pair) ||
                local.RolePairMeasurements.Count(item =>
                    item.FirstMoleculeRole == pair.FirstMoleculeRole &&
                    item.SecondMoleculeRole == pair.SecondMoleculeRole &&
                    item.PairsWithinSearchRadius >= 0 &&
                    (item.PairsWithinSearchRadius == 0 && item.MinimumDistanceAngstrom is null ||
                     item.PairsWithinSearchRadius > 0 && item.MinimumDistanceAngstrom is double distance &&
                     double.IsFinite(distance) && distance > 0)) != 1) ||
            local.LocatedContacts.Any(contact => contact.FirstAtomIndex < 0 ||
                contact.SecondAtomIndex < 0 || !Enum.IsDefined(contact.FirstMoleculeRole) ||
                !Enum.IsDefined(contact.SecondMoleculeRole) || !double.IsFinite(contact.DistanceAngstrom) ||
                !double.IsFinite(contact.RadiusSumAngstrom)))
            return false;
        return spec.RequiredMetricNames.All(name => local.Measurements.Count(value =>
            value.Name == name && double.IsFinite(value.Value)) == 1);
    }

    private static bool HasOrganizationCriterion(ImmutableArray<PreparationAssessmentCriterion> criteria,
        string name, string scope, bool positiveMinimum)
    {
        var matching = criteria.Where(item => item.MeasurementName == name).ToArray();
        return matching.Length == 1 && matching[0].Unit == "angstrom" && matching[0].Scope == scope &&
               matching[0].Minimum is double minimum && double.IsFinite(minimum) &&
               (!positiveMinimum || minimum > 0) &&
               matching[0].Maximum is double maximum && double.IsFinite(maximum) &&
               maximum >= minimum;
    }

    private static PreparationQualification AssessStageProteinGeometry(CompletedStage stage,
        ApplicablePreparationPolicy policy)
    {
        var spec = policy.StageProteinGeometryMeasurement;
        var observed = stage.Observation.ProteinGeometry;
        var required = new[] { "covalentBond", "chainContinuity", "nonbondedDistance" };
        if (spec is null || spec.RequiredKinds.IsDefaultOrEmpty ||
            required.Any(kind => !spec.RequiredKinds.Contains(kind)) ||
            spec.RequiredKinds.Distinct(StringComparer.Ordinal).Count() != spec.RequiredKinds.Length ||
            spec.RequiredKinds.Any(kind => !required.Contains(kind)) ||
            spec.AtomRadiusByElementAngstrom.IsEmpty ||
            spec.AtomRadiusByElementAngstrom.Any(item => string.IsNullOrWhiteSpace(item.Key) ||
                !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(spec.NeighborSearchRadiusAngstrom) ||
            spec.NeighborSearchRadiusAngstrom <= 0 || spec.ExcludedBondHops < 0 ||
            spec.MaximumReportedPairs <= 0 || policy.StageProteinGeometryCriteria.IsDefaultOrEmpty ||
            observed is null || observed.Standing != ObservationStanding.Observed || observed.Kinds.IsDefault ||
            observed.LocatedDistances.IsDefault ||
            observed.Kinds.Select(item => item.Kind).Distinct(StringComparer.Ordinal).Count() != observed.Kinds.Length ||
            observed.LocatedDistances.Any(item => !double.IsFinite(item.DistanceAngstrom) ||
                item.DistanceAngstrom < 0))
            return PreparationQualification.Indeterminate;
        foreach (var kindName in spec.RequiredKinds)
        {
            var rules = policy.StageProteinGeometryCriteria.Where(item =>
                item.StageKind == stage.Kind && item.Criterion?.Kind == kindName).ToArray();
            if (rules.Length != 1)
                return PreparationQualification.Indeterminate;
            var criterion = rules[0].Criterion;
            if ((criterion.MinimumObservedAngstrom is null && criterion.MaximumObservedAngstrom is null) ||
                kindName == "covalentBond" && criterion.AllowNotApplicable ||
                criterion.MinimumObservedAngstrom is double lower && !double.IsFinite(lower) ||
                criterion.MaximumObservedAngstrom is double upper && !double.IsFinite(upper) ||
                criterion.MinimumObservedAngstrom is double minimum &&
                    criterion.MaximumObservedAngstrom is double maximum && minimum > maximum)
                return PreparationQualification.Indeterminate;
            var measurements = observed.Kinds.Where(item => item.Kind == kindName).ToArray();
            if (measurements.Length != 1)
                return PreparationQualification.Indeterminate;
            var measure = measurements[0];
            if (measure.Standing == GeometryKindStanding.NotApplicable && criterion.AllowNotApplicable &&
                measure.EligibleCount == 0 && measure.MeasuredCount == 0 &&
                measure.MinimumDistanceAngstrom is null && measure.MaximumDistanceAngstrom is null)
                continue;
            if (measure.Standing != GeometryKindStanding.Observed || measure.EligibleCount <= 0 ||
                measure.MeasuredCount != measure.EligibleCount ||
                measure.MinimumDistanceAngstrom is not double observedMinimum ||
                measure.MaximumDistanceAngstrom is not double observedMaximum ||
                !double.IsFinite(observedMinimum) || !double.IsFinite(observedMaximum) ||
                observedMinimum > observedMaximum)
                return PreparationQualification.Indeterminate;
            if (criterion.MinimumObservedAngstrom is double supportedMinimum &&
                    observedMinimum < supportedMinimum ||
                criterion.MaximumObservedAngstrom is double supportedMaximum &&
                    observedMaximum > supportedMaximum)
                return PreparationQualification.NotQualified;
        }
        return PreparationQualification.QualifiedPrepared;
    }

    private static bool ContactRulesAvailable(StageKind kind, ApplicablePreparationPolicy policy)
    {
        var rules = policy.ContactCriteria;
        if (rules.IsDefaultOrEmpty ||
            rules.Where(item => item.StageKind == kind)
                .Select(item => (item.FirstMoleculeRole, item.SecondMoleculeRole))
                .Distinct().Count() != rules.Count(item => item.StageKind == kind) ||
            rules.Where(item => item.StageKind == kind).Any(item =>
                item.MinimumPairsWithinSearchRadius < 0 ||
                !policy.LocalStateObservation.ContactRolePairs.Any(pair =>
                    pair.FirstMoleculeRole == item.FirstMoleculeRole &&
                    pair.SecondMoleculeRole == item.SecondMoleculeRole) ||
                item.MinimumNearestDistanceAngstrom is double minimum &&
                    (!double.IsFinite(minimum) || minimum <= 0) ||
                item.MaximumNearestDistanceAngstrom is double maximum &&
                    (!double.IsFinite(maximum) || maximum <= 0 ||
                     maximum > policy.LocalStateObservation.ContactSearchRadiusAngstrom) ||
                item.MinimumNearestDistanceAngstrom is double minimum2 &&
                    item.MaximumNearestDistanceAngstrom is double maximum2 && minimum2 > maximum2))
            return false;
        return rules.Any(item => item.StageKind == kind &&
            item.FirstMoleculeRole == MoleculeRoleKind.Protein && item.SecondMoleculeRole == MoleculeRoleKind.Lipid &&
            item.MinimumPairsWithinSearchRadius > 0 &&
            item.MaximumNearestDistanceAngstrom is double maximum &&
            double.IsFinite(maximum) && maximum > 0 &&
            maximum <= policy.LocalStateObservation.ContactSearchRadiusAngstrom);
    }

    private static bool ContactRulesMet(StageKind kind, LocalStateObservations observed,
        ApplicablePreparationPolicy policy)
    {
        foreach (var rule in policy.ContactCriteria.Where(item => item.StageKind == kind))
        {
            var pair = observed.RolePairMeasurements.FirstOrDefault(item =>
                item.FirstMoleculeRole == rule.FirstMoleculeRole &&
                item.SecondMoleculeRole == rule.SecondMoleculeRole);
            if (pair is null || pair.PairsWithinSearchRadius < rule.MinimumPairsWithinSearchRadius ||
                rule.MinimumNearestDistanceAngstrom is double minimum &&
                    (pair.MinimumDistanceAngstrom is not double distance || distance < minimum) ||
                rule.MaximumNearestDistanceAngstrom is double maximum &&
                    (pair.MinimumDistanceAngstrom is not double distance2 || distance2 > maximum))
                return false;
        }
        return true;
    }
}
