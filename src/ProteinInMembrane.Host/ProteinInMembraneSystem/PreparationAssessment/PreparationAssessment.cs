using System.Collections.Immutable;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment;

/// <summary>Interprets declared technical checks of one exact completed stage.</summary>
public sealed class PreparationAssessment
{
    public PreparationAssessmentResult Assess(CompletedStage stage, ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy,
        ImmutableArray<ScientificFinding> laterFindings)
    {
        ArgumentNullException.ThrowIfNull(stage);
        ArgumentNullException.ThrowIfNull(constructed);
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(membrane);
        ArgumentNullException.ThrowIfNull(placement);
        ArgumentNullException.ThrowIfNull(policy);

        var evidence = ImmutableArray.CreateBuilder<ScientificEvidence>();
        if (!stage.Observation.Evidence.IsDefault) evidence.AddRange(stage.Observation.Evidence);
        var findings = ImmutableArray.CreateBuilder<ScientificFinding>();
        if (!stage.Findings.IsDefault) findings.AddRange(stage.Findings);
        if (!laterFindings.IsDefault) findings.AddRange(laterFindings);
        var issues = new List<string>();
        var missing = new List<string>();
        var applicable = Corresponds(stage, constructed, protein, membrane, placement, policy);
        var relevantSubjects = new HashSet<string>(StringComparer.Ordinal)
        {
            stage.Id, stage.Attempt.Id, constructed.Id, protein.Id,
            membrane.Id, membrane.Intended.Id, placement.Id, placement.Proposal.Id
        };
        foreach (var finding in findings.Where(item => item.Material &&
                     relevantSubjects.Contains(item.SubjectId)))
        {
            if (finding.Disposition == FindingDisposition.Disqualifies)
                issues.Add($"Material finding: {finding.Meaning}");
            else if (finding.Disposition == FindingDisposition.Challenges)
                missing.Add($"Material finding requiring a check: {finding.Meaning}");
        }

        void Report(string check, string detail, bool issue)
        {
            var description = $"{check}: {detail}";
            (issue ? issues : missing).Add(description);
            var observation = new ScientificEvidence(Guid.NewGuid().ToString("N"), stage.Id,
                "Preparation Assessment", check, detail,
                $"Stage {stage.Id}; attempt {stage.Attempt.Id}; policy {policy.Id} {policy.Version}",
                issue ? "Observed technical condition." : "Required observation or rule unavailable.",
                issue ? EvidenceBearing.Contradicts : EvidenceBearing.Unknown);
            evidence.Add(observation);
            findings.Add(new ScientificFinding(Guid.NewGuid().ToString("N"), stage.Id,
                observation.Id, description,
                issue ? "An observed technical issue remains for this stage." :
                    "This required technical check remains incomplete.",
                issue ? FindingDisposition.Disqualifies : FindingDisposition.Challenges,
                true, DateTimeOffset.UtcNow));
        }

        if (!applicable)
            Report("Exact correspondence", "Stage, attempt, inputs, policy or molecular identity do not form one applicable account.", false);
        else
        {
            CheckLocal(stage, policy, Report);
            CheckGeometry(stage, policy, Report);
            CheckStageValues(stage, policy, Report);
            if (policy.Construction.Route == ConstructionRouteKind.PackmolMemgen)
            {
                var saltBranch = constructed.Derivation.Conditions?.SaltBranch;
                if (saltBranch == "neutralizationOnly")
                    Report("Nominal NaCl treatment",
                        "The provider added neutralization only; no background salt pairs were established.", true);
                else if (saltBranch != "chargeCompensated")
                    Report("Nominal NaCl treatment", "The provider's actual salt branch is unavailable or unrecognized.", false);
            }
            if (stage.Kind == StageKind.Equilibration &&
                (stage.Observation.ObservationAdequacy != EquilibrationObservationAdequacy.Adequate ||
                 stage.Observation.EquilibrationAssessments.IsDefaultOrEmpty ||
                 stage.Observation.EquilibrationAssessments.Any(item => !item.Sufficient)))
                Report("Optional procedure", "The declared time-series observations are absent or insufficient.", false);
        }

        // A known failure takes precedence, while each unavailable check remains visible.
        var standing = issues.Count > 0 ? PreparationCheckStanding.IssuesFound :
            missing.Count > 0 ? PreparationCheckStanding.ChecksIncomplete : PreparationCheckStanding.ChecksPassed;
        var reason = standing switch
        {
            PreparationCheckStanding.IssuesFound => "Issues found: " + string.Join(" ", issues) +
                (missing.Count == 0 ? string.Empty : " Checks incomplete: " + string.Join(" ", missing)),
            PreparationCheckStanding.ChecksIncomplete => "Checks incomplete: " + string.Join(" ", missing),
            _ => "All declared technical checks have corresponding observations and passed."
        };
        return new PreparationAssessmentResult(Guid.NewGuid().ToString("N"), stage.Id,
            standing, reason, evidence.ToImmutable(), findings.ToImmutable(),
            policy.Limitations.IsDefault ? ImmutableArray<string>.Empty : policy.Limitations,
            DateTimeOffset.UtcNow, applicable);
    }

    private static bool Corresponds(CompletedStage stage, ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy)
    {
        var attempt = stage.Attempt;
        return attempt.Id == constructed.Attempt.Id &&
            attempt.StudyRevisionId == constructed.Attempt.StudyRevisionId &&
            attempt.StudyRevisionId == protein.StudyRevisionId &&
            attempt.StudyRevisionId == membrane.StudyRevisionId &&
            attempt.StudyRevisionId == placement.StudyRevisionId &&
            stage.PolicyId == policy.Id && PreparationPolicyFingerprint.Matches(attempt, policy) &&
            attempt.ProteinId == protein.Id && attempt.MembraneId == membrane.Id &&
            attempt.PlacementId == placement.Id && placement.Standing == AssessmentStanding.Supported &&
            stage.Correspondence.Complete && constructed.Correspondence.Complete &&
            stage.Correspondence.ResultId == stage.Molecule.Id &&
            stage.Correspondence.SourceId == constructed.Correspondence.SourceId &&
            !stage.Correspondence.Atoms.IsDefault && !constructed.Correspondence.Atoms.IsDefault &&
            stage.Correspondence.Atoms.SequenceEqual(constructed.Correspondence.Atoms) &&
            stage.Molecule.AtomCount == constructed.Molecule.AtomCount &&
            stage.Molecule.Id == stage.Id && stage.Observation.StageId == stage.Id &&
            stage.Observation.AttemptId == attempt.Id && stage.Observation.Kind == stage.Kind &&
            (stage.Kind == StageKind.Minimization && stage.Observation.Termination == StageTermination.Converged ||
             stage.Kind == StageKind.Equilibration && stage.Observation.Termination == StageTermination.Completed);
    }

    private static void CheckLocal(CompletedStage stage, ApplicablePreparationPolicy policy,
        Action<string, string, bool> report)
    {
        var spec = policy.LocalStateObservation;
        var local = stage.Observation.LocalState;
        if (spec is null || spec.ContactRolePairs.IsDefaultOrEmpty || spec.RequiredMetricNames.IsDefaultOrEmpty ||
            spec.RequiredMetricNames.Distinct(StringComparer.Ordinal).Count() != spec.RequiredMetricNames.Length)
        {
            report("Local-state policy", "Required molecular coverage is not declared coherently.", false);
            return;
        }
        if (local is null || local.Standing != ObservationStanding.Observed)
        {
            report("Stage local state", local?.UnavailableReason ?? "No observed completed-stage local state.", false);
            return;
        }
        if (local.Measurements.IsDefault || local.LocatedContacts.IsDefault ||
            local.CoveredRolePairs.IsDefault || local.RolePairMeasurements.IsDefault ||
            local.CoveredRolePairs.Distinct().Count() != local.CoveredRolePairs.Length ||
            local.RolePairMeasurements.Select(item => (item.FirstMoleculeRole,
                item.SecondMoleculeRole)).Distinct().Count() != local.RolePairMeasurements.Length ||
            local.LocatedContacts.Any(item => item.FirstAtomIndex < 0 || item.SecondAtomIndex < 0 ||
                !Enum.IsDefined(item.FirstMoleculeRole) || !Enum.IsDefined(item.SecondMoleculeRole) ||
                !double.IsFinite(item.DistanceAngstrom) || item.DistanceAngstrom < 0 ||
                !double.IsFinite(item.RadiusSumAngstrom) || item.RadiusSumAngstrom <= 0))
            report("Stage local state", "Measurement or located-contact coverage is incomplete or malformed.", false);

        foreach (var name in spec.RequiredMetricNames)
        {
            var values = local.Measurements.IsDefault ? Array.Empty<MeasuredValue>() :
                local.Measurements.Where(item => item.Name == name).ToArray();
            var expectedScope = name switch
            {
                "minimumIntermolecularDistanceAngstrom" or
                    "minimumIntermolecularHeavyAtomDistanceAngstrom" => "wholeSystem",
                "upperLipidHeadMeanZAngstrom" => "upperLeaflet",
                "lowerLipidHeadMeanZAngstrom" => "lowerLeaflet",
                "leafletHeadSeparationAngstrom" => "bilayer",
                "proteinBilayerMidplaneOffsetAngstrom" => "proteinVsBilayer",
                _ => null
            };
            if (expectedScope is null || values.Length != 1 || !double.IsFinite(values[0].Value) ||
                values[0].Unit != "angstrom" || values[0].Scope != expectedScope)
                report($"Local measurement {name}",
                    "One finite observation with the declared angstrom unit and molecular scope is required.", false);
        }
        foreach (var pair in spec.ContactRolePairs)
        {
            var values = local.RolePairMeasurements.IsDefault ? Array.Empty<LocalRolePairMeasurement>() :
                local.RolePairMeasurements.Where(item => item.FirstMoleculeRole == pair.FirstMoleculeRole &&
                    item.SecondMoleculeRole == pair.SecondMoleculeRole).ToArray();
            if (local.CoveredRolePairs.IsDefault || !local.CoveredRolePairs.Contains(pair) ||
                values.Length != 1 || !ValidPair(values[0]))
                report($"Contact coverage {pair.FirstMoleculeRole}–{pair.SecondMoleculeRole}",
                    "The declared role pair lacks one valid count and nearest distance.", false);
        }
        if (policy.ContactCriteria.IsDefault)
        {
            report("Stage contact criteria", "The technical contact-rule account is unavailable.", false);
            return;
        }
        foreach (var rule in policy.ContactCriteria.Where(item => item.StageKind == stage.Kind))
        {
            var check = $"Contact {rule.FirstMoleculeRole}–{rule.SecondMoleculeRole}";
            if (rule.MinimumPairsWithinSearchRadius < 0 ||
                rule.MinimumNearestDistanceAngstrom is double minimum && (!double.IsFinite(minimum) || minimum <= 0) ||
                rule.MaximumNearestDistanceAngstrom is double maximum && (!double.IsFinite(maximum) || maximum <= 0) ||
                rule.MinimumNearestDistanceAngstrom is double lower &&
                    rule.MaximumNearestDistanceAngstrom is double upper && lower > upper)
            {
                report(check, "The declared contact bound is invalid.", false);
                continue;
            }
            var pair = new LocalContactRolePair(rule.FirstMoleculeRole, rule.SecondMoleculeRole);
            var values = local.RolePairMeasurements.IsDefault ? Array.Empty<LocalRolePairMeasurement>() :
                local.RolePairMeasurements.Where(item => item.FirstMoleculeRole == pair.FirstMoleculeRole &&
                    item.SecondMoleculeRole == pair.SecondMoleculeRole).ToArray();
            if (values.Length != 1 || !ValidPair(values[0]) || local.CoveredRolePairs.IsDefault ||
                !local.CoveredRolePairs.Contains(pair))
            {
                report(check, "The corresponding role-pair observation is unavailable or malformed.", false);
                continue;
            }
            var observed = values[0];
            if (observed.PairsWithinSearchRadius < rule.MinimumPairsWithinSearchRadius ||
                rule.MinimumNearestDistanceAngstrom is double min &&
                    (observed.MinimumDistanceAngstrom is not double distance || distance < min) ||
                rule.MaximumNearestDistanceAngstrom is double max &&
                    (observed.MinimumDistanceAngstrom is not double distance2 || distance2 > max))
                report(check, "Measured count or nearest distance violates its declared bound.", true);
        }
    }

    private static bool ValidPair(LocalRolePairMeasurement item) =>
        item.PairsWithinSearchRadius >= 0 &&
        (item.PairsWithinSearchRadius == 0 && item.MinimumDistanceAngstrom is null ||
         item.PairsWithinSearchRadius > 0 && item.MinimumDistanceAngstrom is double distance &&
             double.IsFinite(distance) && distance > 0);

    private static void CheckStageValues(CompletedStage stage, ApplicablePreparationPolicy policy,
        Action<string, string, bool> report)
    {
        if (policy.AssessmentCriteria.IsDefault)
        {
            report("Stage measurement criteria", "The declared criterion account is unavailable.", false);
            return;
        }
        foreach (var criterion in policy.AssessmentCriteria.Where(item => item.StageKind == stage.Kind))
        {
            var check = $"Stage measurement {criterion.MeasurementName}";
            if (string.IsNullOrWhiteSpace(criterion.MeasurementName) ||
                string.IsNullOrWhiteSpace(criterion.Unit) || string.IsNullOrWhiteSpace(criterion.Scope) ||
                criterion.Minimum is null && criterion.Maximum is null ||
                criterion.Minimum is double minimum && !double.IsFinite(minimum) ||
                criterion.Maximum is double maximum && !double.IsFinite(maximum) ||
                criterion.Minimum is double lower && criterion.Maximum is double upper && lower > upper)
            {
                report(check, "The declared measurement bound is invalid.", false);
                continue;
            }
            var values = stage.Observation.Measurements.IsDefault ? Array.Empty<MeasuredValue>() :
                stage.Observation.Measurements.Where(item => item.Name == criterion.MeasurementName).ToArray();
            if (values.Length != 1 || values[0].Unit != criterion.Unit ||
                values[0].Scope != criterion.Scope || !double.IsFinite(values[0].Value))
            {
                report(check, "One finite observation with the declared unit and scope is required.", false);
                continue;
            }
            if (criterion.Minimum is double min && values[0].Value < min ||
                criterion.Maximum is double max && values[0].Value > max)
                report(check, $"Observed {values[0].Value:G17} {criterion.Unit} violates the declared bound.", true);
        }
    }

    private static void CheckGeometry(CompletedStage stage, ApplicablePreparationPolicy policy,
        Action<string, string, bool> report)
    {
        var spec = policy.StageProteinGeometryMeasurement;
        var observed = stage.Observation.ProteinGeometry;
        if (spec is null || spec.RequiredKinds.IsDefaultOrEmpty ||
            spec.RequiredKinds.Distinct(StringComparer.Ordinal).Count() != spec.RequiredKinds.Length ||
            !new[] { "covalentBond", "chainContinuity", "nonbondedDistance" }.All(spec.RequiredKinds.Contains) ||
            spec.RequiredKinds.Any(kind => kind is not ("covalentBond" or "chainContinuity" or "nonbondedDistance")) ||
            spec.AtomRadiusByElementAngstrom.IsEmpty ||
            spec.AtomRadiusByElementAngstrom.Any(item => string.IsNullOrWhiteSpace(item.Key) ||
                !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(spec.NeighborSearchRadiusAngstrom) || spec.NeighborSearchRadiusAngstrom <= 0 ||
            spec.ExcludedBondHops < 0 || spec.MaximumReportedPairs <= 0 ||
            policy.StageProteinGeometryCriteria.IsDefault)
        {
            report("Protein geometry policy", "Required bond, chain and nonbonded checks are not declared coherently.", false);
            return;
        }
        if (observed is null || observed.Standing != ObservationStanding.Observed ||
            observed.Kinds.IsDefault || observed.LocatedDistances.IsDefault ||
            observed.LocatedDistances.Any(item => !double.IsFinite(item.DistanceAngstrom) ||
                item.DistanceAngstrom < 0))
        {
            report("Stage protein geometry", "A complete stage-specific geometry observation is unavailable.", false);
            return;
        }
        foreach (var kind in spec.RequiredKinds)
        {
            var check = $"Protein geometry {kind}";
            var criteria = policy.StageProteinGeometryCriteria.Where(item =>
                item.StageKind == stage.Kind && item.Criterion?.Kind == kind).ToArray();
            if (criteria.Length != 1 || criteria[0].Criterion is not { } criterion ||
                criterion.MinimumObservedAngstrom is null && criterion.MaximumObservedAngstrom is null ||
                kind == "covalentBond" && criterion.AllowNotApplicable ||
                criterion.MinimumObservedAngstrom is double minimum && !double.IsFinite(minimum) ||
                criterion.MaximumObservedAngstrom is double maximum && !double.IsFinite(maximum) ||
                criterion.MinimumObservedAngstrom is double lower &&
                    criterion.MaximumObservedAngstrom is double upper && lower > upper)
            {
                report(check, "The stage-specific geometry condition is missing or invalid.", false);
                continue;
            }
            var values = observed.Kinds.Where(item => item.Kind == kind).ToArray();
            if (values.Length != 1)
            {
                report(check, "Exactly one geometry observation is required.", false);
                continue;
            }
            var value = values[0];
            if (value.Standing == GeometryKindStanding.NotApplicable && criterion.AllowNotApplicable &&
                value.EligibleCount == 0 && value.MeasuredCount == 0 &&
                value.MinimumDistanceAngstrom is null && value.MaximumDistanceAngstrom is null)
                continue;
            if (value.Standing != GeometryKindStanding.Observed || value.EligibleCount <= 0 ||
                value.MeasuredCount != value.EligibleCount ||
                value.MinimumDistanceAngstrom is not double observedMin ||
                value.MaximumDistanceAngstrom is not double observedMax ||
                !double.IsFinite(observedMin) || !double.IsFinite(observedMax) || observedMin > observedMax)
            {
                report(check, value.UnavailableReason ?? "The exact geometry measurement is incomplete.", false);
                continue;
            }
            if (criterion.MinimumObservedAngstrom is double min && observedMin < min ||
                criterion.MaximumObservedAngstrom is double max && observedMax > max)
                report(check, $"Observed distance range {observedMin:G6}–{observedMax:G6} Å violates its declared bound.", true);
        }
    }
}
