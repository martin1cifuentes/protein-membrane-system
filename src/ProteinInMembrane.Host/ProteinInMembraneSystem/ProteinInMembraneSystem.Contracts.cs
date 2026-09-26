using System.Buffers;
using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem;

// These values belong to the product root. A provider artifact is an observation,
// not an assessed protein, supported placement, completed stage, or qualification.
public sealed record StudyRevision(
    string Id,
    long Number,
    IntendedProteinModel? IntendedProtein,
    MembraneModel? Membrane,
    string? AdoptedPlacementProposalId,
    FixedStudyConditions Conditions);

public sealed record FixedStudyConditions(double NominalPh, double TargetNaClMolar, double OptionalTemperatureKelvin)
{
    public static FixedStudyConditions Initial { get; } = new(7.0, 0.15, 303.0);
}

[JsonConverter(typeof(JsonStringEnumConverter<ResearcherDecisionKind>))]
public enum ResearcherDecisionKind
{
    [JsonStringEnumMemberName("approvePreparationChange")] ApprovePreparationChange = 1,
    [JsonStringEnumMemberName("adoptMembrane")] AdoptMembrane,
    [JsonStringEnumMemberName("adoptPlacement")] AdoptPlacement
}

[JsonConverter(typeof(JsonStringEnumConverter<ResearcherDecisionValue>))]
public enum ResearcherDecisionValue
{
    [JsonStringEnumMemberName("approved")] Approved = 1,
    [JsonStringEnumMemberName("declined")] Declined,
    [JsonStringEnumMemberName("adopted")] Adopted
}

public sealed record ResearcherDecision(
    string Id,
    string StudyRevisionId,
    string SubjectId,
    ResearcherDecisionKind Kind,
    ResearcherDecisionValue ChosenValue,
    DateTimeOffset At,
    string Rationale);

// These are product-owned choices, distinct from open provenance descriptions and provider data.
[JsonConverter(typeof(JsonStringEnumConverter<SourceRouteKind>))]
public enum SourceRouteKind
{
    [JsonStringEnumMemberName("rcsb")] Rcsb = 1,
    [JsonStringEnumMemberName("alphafold")] AlphaFold,
    [JsonStringEnumMemberName("upload")] Upload
}

[JsonConverter(typeof(JsonStringEnumConverter<UploadOriginKind>))]
public enum UploadOriginKind
{
    [JsonStringEnumMemberName("predicted")] Predicted = 1,
    [JsonStringEnumMemberName("experimental")] Experimental,
    [JsonStringEnumMemberName("unknown")] Unknown
}

public sealed record StructuralSource(
    string Id,
    SourceRouteKind Kind,
    string Provenance,
    string CoordinatePath,
    string Sha256,
    string? Accession,
    string? SourceModelDescription,
    PredictionEvidenceAsset? Prediction = null,
    UploadOriginKind? UploadProvenance = null,
    string? UploadProvenanceNote = null);

/// <summary>An identified prediction record; coordinate and PAE assets must belong to the same record.</summary>
[JsonConverter(typeof(JsonStringEnumConverter<PaeAcquisitionStanding>))]
public enum PaeAcquisitionStanding
{
    [JsonStringEnumMemberName("Available")] Available = 1,
    [JsonStringEnumMemberName("Unavailable")] Unavailable
}

[JsonConverter(typeof(JsonStringEnumConverter<PredictionObservationStanding>))]
public enum PredictionObservationStanding
{
    [JsonStringEnumMemberName("Observed")] Observed = 1,
    [JsonStringEnumMemberName("Unavailable")] Unavailable,
    [JsonStringEnumMemberName("Unmapped")] Unmapped
}

public sealed record PredictionEvidenceAsset(
    string RecordId,
    string? ModelVersion,
    int? SequenceStart,
    int? SequenceEnd,
    string CoordinateUrl,
    string CoordinateSha256,
    string? PaeUrl,
    string? PaePath,
    string? PaeSha256,
    PaeAcquisitionStanding PaeStanding,
    string? PaeReason);

public sealed record CandidateSourceRecord(
    string Id,
    string Label,
    SourceRouteKind Kind,
    string Provenance,
    ImmutableArray<string> Limitations,
    StructuralSource? ResolvedSource);

public sealed record SourceInspectionReport(
    StructuralSource Source,
    string Format,
    ImmutableArray<SourceModelObservation> Models,
    ImmutableArray<string> Limitations,
    PredictionEvidenceObservations? Prediction = null);

public sealed record ResidueAddress(
    int Model,
    string Chain,
    int Residue,
    string InsertionCode,
    string CopyId);
public sealed record AtomAddress(ResidueAddress Residue, string AtomName);

public sealed record ChainSelection(string SourceChain, string CopyId);
public sealed record PartnerSelection(string SourceId, bool Retain, string Reason);
public sealed record AlternateLocationChoice(ResidueAddress Residue, string Altloc, string DecisionId);
public sealed record ResidueVariantChoice(ResidueAddress Residue, string Variant, string? DecisionId);
public sealed record HeavyAtomApproval(ResidueAddress Residue, string AtomName, string DecisionId);
public sealed record DisulfideChoice(ResidueAddress First, ResidueAddress Second, string DecisionId);
public sealed record DisulfideBond(ResidueAddress First, ResidueAddress Second);

public sealed record IntendedProteinModel(
    string Id,
    StructuralSource Source,
    int ModelIndex,
    string? BiologicalAssemblyId,
    ImmutableArray<ChainSelection> Chains,
    ImmutableArray<PartnerSelection> Partners,
    ImmutableArray<AlternateLocationChoice> AlternateLocations);

[JsonConverter(typeof(JsonStringEnumConverter<PreparationChangeKind>))]
public enum PreparationChangeKind
{
    [JsonStringEnumMemberName("alternateLocation")] AlternateLocation = 1,
    [JsonStringEnumMemberName("residueState")] ResidueState,
    [JsonStringEnumMemberName("heavyAtom")] HeavyAtom,
    [JsonStringEnumMemberName("disulfide")] Disulfide
}

public sealed record PreparationChangeProposal(
    string Id,
    string StudyRevisionId,
    string IntendedProteinId,
    ResidueAddress Residue,
    PreparationChangeKind Kind,
    string ProposedChange,
    string Rationale,
    ImmutableArray<string> Limitations,
    bool ApprovalRequired,
    ResidueAddress? PartnerResidue = null);
public sealed record PreparationProposalReport(
    string StudyRevisionId,
    string IntendedProteinId,
    ImmutableArray<PreparationChangeProposal> Changes,
    WorkerArtifact? Preview,
    ImmutableArray<PreviewChainCorrespondence> PreviewChains,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<string> UnresolvedQuestions,
    ProteinGeometryObservations? SourceGeometry = null);

public interface IObservedDiagnostic { }
/// <summary>Observed but unpromoted prepared candidate; not an Assessed Prepared Protein.</summary>
public sealed record ProteinPreparationDiagnostic(
    MolecularArtifact Candidate,
    ProteinGeometryObservations Geometry,
    SourceToResultCorrespondence Correspondence,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings) : IObservedDiagnostic;

public sealed record PreviewChainCorrespondence(string SourceChain, string CopyId, string PreviewChain);

[JsonConverter(typeof(JsonStringEnumConverter<AtomOriginKind>))]
public enum AtomOriginKind
{
    [JsonStringEnumMemberName("source")] Source = 1,
    [JsonStringEnumMemberName("generated")] Generated
}

[JsonConverter(typeof(JsonStringEnumConverter<MoleculeRoleKind>))]
public enum MoleculeRoleKind
{
    [JsonStringEnumMemberName("protein")] Protein = 1,
    [JsonStringEnumMemberName("retainedPartner")] RetainedPartner,
    [JsonStringEnumMemberName("lipid")] Lipid,
    [JsonStringEnumMemberName("water")] Water,
    [JsonStringEnumMemberName("ion")] Ion
}

[JsonConverter(typeof(JsonStringEnumConverter<AtomRoleKind>))]
public enum AtomRoleKind
{
    [JsonStringEnumMemberName("backbone")] Backbone = 1,
    [JsonStringEnumMemberName("sidechain")] Sidechain,
    [JsonStringEnumMemberName("partnerAtom")] PartnerAtom,
    [JsonStringEnumMemberName("head")] Head,
    [JsonStringEnumMemberName("body")] Body
}

[JsonConverter(typeof(JsonStringEnumConverter<GeneratedComponentRoleKind>))]
public enum GeneratedComponentRoleKind
{
    [JsonStringEnumMemberName("lipid")] Lipid = 1,
    [JsonStringEnumMemberName("water")] Water,
    [JsonStringEnumMemberName("positiveIon")] PositiveIon,
    [JsonStringEnumMemberName("negativeIon")] NegativeIon
}

public static class MoleculeRoleTokens
{
    public static string WireToken(this MoleculeRoleKind role) => role switch
    {
        MoleculeRoleKind.Protein => "protein",
        MoleculeRoleKind.RetainedPartner => "retainedPartner",
        MoleculeRoleKind.Lipid => "lipid",
        MoleculeRoleKind.Water => "water",
        MoleculeRoleKind.Ion => "ion",
        _ => throw new ArgumentOutOfRangeException(nameof(role))
    };

    public static string PairKey(MoleculeRoleKind first, MoleculeRoleKind second) =>
        $"{first.WireToken()}|{second.WireToken()}";

    public static bool TryParsePairKey(string? key, out LocalContactRolePair? pair)
    {
        pair = null;
        if (key is null) return false;
        var separator = key.IndexOf('|');
        if (separator <= 0 || separator == key.Length - 1 || key.IndexOf('|', separator + 1) >= 0 ||
            !TryParseRole(key[..separator], out var first) ||
            !TryParseRole(key[(separator + 1)..], out var second)) return false;
        pair = new LocalContactRolePair(first, second);
        return true;
    }

    private static bool TryParseRole(string token, out MoleculeRoleKind role)
    {
        role = token switch
        {
            "protein" => MoleculeRoleKind.Protein,
            "retainedPartner" => MoleculeRoleKind.RetainedPartner,
            "lipid" => MoleculeRoleKind.Lipid,
            "water" => MoleculeRoleKind.Water,
            "ion" => MoleculeRoleKind.Ion,
            _ => default
        };
        return Enum.IsDefined(role);
    }
}

public sealed class CoveredRolePairsJsonConverter : JsonConverter<ImmutableArray<LocalContactRolePair>>
{
    public override ImmutableArray<LocalContactRolePair> Read(ref Utf8JsonReader reader, Type typeToConvert,
        JsonSerializerOptions options)
    {
        if (reader.TokenType != JsonTokenType.StartArray) throw new JsonException("Covered role pairs must be an array.");
        var pairs = ImmutableArray.CreateBuilder<LocalContactRolePair>();
        while (reader.Read())
        {
            if (reader.TokenType == JsonTokenType.EndArray) return pairs.ToImmutable();
            if (reader.TokenType != JsonTokenType.String ||
                !MoleculeRoleTokens.TryParsePairKey(reader.GetString(), out var pair))
                throw new JsonException("An unrecognized covered role pair was observed.");
            pairs.Add(pair!);
        }
        throw new JsonException("The covered role-pair array is incomplete.");
    }

    public override void Write(Utf8JsonWriter writer, ImmutableArray<LocalContactRolePair> value,
        JsonSerializerOptions options)
    {
        writer.WriteStartArray();
        foreach (var pair in value)
            writer.WriteStringValue(MoleculeRoleTokens.PairKey(pair.FirstMoleculeRole, pair.SecondMoleculeRole));
        writer.WriteEndArray();
    }
}

public sealed record AtomCorrespondence(
    int ResultAtomIndex,
    string ResultAtomId,
    string? SourceAtomId,
    AtomOriginKind Role,
    MoleculeRoleKind MoleculeRole,
    AtomRoleKind AtomRole,
    string Element,
    ResidueAddress? SourceResidue,
    string? ApprovedChangeId,
    LeafletSide? PhysicalSide = null,
    string? GeneratedSpeciesId = null,
    GeneratedComponentRoleKind? GeneratedComponentRole = null);

public sealed record SourceToResultCorrespondence(
    string SourceId,
    string ResultId,
    ImmutableArray<AtomCorrespondence> Atoms,
    bool Complete);

public sealed record MolecularArtifact(
    string Id,
    string CoordinatePath,
    string CoordinateSha256,
    string? TopologyPath,
    string? SystemXmlPath,
    string? StateXmlPath,
    int AtomCount,
    string? CellDescription,
    string? TopologySha256 = null,
    string? CorrespondencePath = null,
    string? CorrespondenceSha256 = null,
    string? SystemXmlSha256 = null,
    string? StateXmlSha256 = null);

public sealed record MolecularRepresentation(
    string SpeciesId,
    string ChemistryId,
    string Category,
    string TemplatePath,
    string TemplateSha256,
    string CoordinateTemplatePath,
    string CoordinateTemplateSha256,
    int AtomCount,
    double NetChargeElementary,
    double AreaPerMoleculeAngstromSquared,
    double VolumeAngstromCubed,
    ImmutableArray<int> HeadAtomIndices,
    string ForceFieldFamily,
    string ForceFieldVersion,
    ImmutableArray<string> Limitations,
    ImmutableArray<MolecularStereoCheck> StereoChecks = default);

/// <summary>Expected geometry for the exact identified molecular configuration.</summary>
public sealed record MolecularStereoCheck(
    string Kind,
    ImmutableArray<string> AtomNames,
    string Expected);

public sealed record ScientificEvidence(
    string Id,
    string SubjectId,
    string Source,
    string Method,
    string Observation,
    string Applicability,
    string Uncertainty,
    EvidenceBearing Bearing);

public enum EvidenceBearing { Context, Supports, Contradicts, Unknown }
public enum FindingDisposition { Context, Supports, Challenges, Disqualifies }

public sealed record ScientificFinding(
    string Id,
    string SubjectId,
    string EvidenceId,
    string Meaning,
    string Consequence,
    FindingDisposition Disposition,
    bool Material,
    DateTimeOffset EstablishedAt);

public sealed record ProteinChemicalStatePolicy(
    string Id,
    string Version,
    string ApplicableChemistry,
    double DisulfideCandidateMaxSgDistanceAngstrom,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles,
    ImmutableDictionary<string, string> DefaultVariants,
    ImmutableArray<string> PermittedVariants,
    ImmutableArray<string> Limitations);

/// <summary>Bounded, evidence-backed protein geometry interpretation, independent of the worker.</summary>
public sealed record ProteinStructuralAssessmentPolicy(
    string Id,
    string Version,
    string ApplicableProteinClass,
    ImmutableArray<string> EvidenceReferences,
    ProteinGeometryMeasurementSpec Measurement,
    ImmutableArray<ProteinGeometryCriterion> Criteria,
    ImmutableArray<string> Limitations);
public sealed record ProteinGeometryMeasurementSpec(
    ImmutableArray<string> RequiredKinds,
    ImmutableDictionary<string, double> AtomRadiusByElementAngstrom,
    double NeighborSearchRadiusAngstrom,
    int ExcludedBondHops,
    int MaximumReportedPairs);
public sealed record ProteinGeometryCriterion(
    string Kind,
    double? MinimumObservedAngstrom,
    double? MaximumObservedAngstrom,
    bool AllowNotApplicable);
public sealed record GeometryDistanceObservation(
    string Kind,
    AtomAddress First,
    AtomAddress Second,
    double DistanceAngstrom,
    double? RadiusSumAngstrom);

[JsonConverter(typeof(JsonStringEnumConverter<ObservationStanding>))]
public enum ObservationStanding
{
    [JsonStringEnumMemberName("Observed")] Observed = 1,
    [JsonStringEnumMemberName("Unavailable")] Unavailable
}

[JsonConverter(typeof(JsonStringEnumConverter<GeometryKindStanding>))]
public enum GeometryKindStanding
{
    [JsonStringEnumMemberName("Observed")] Observed = 1,
    [JsonStringEnumMemberName("Unavailable")] Unavailable,
    [JsonStringEnumMemberName("NotApplicable")] NotApplicable
}

public sealed record ProteinGeometryKindObservation(
    string Kind,
    GeometryKindStanding Standing,
    int EligibleCount,
    int MeasuredCount,
    double? MinimumDistanceAngstrom,
    double? MaximumDistanceAngstrom,
    string? UnavailableReason);
public sealed record ProteinGeometryObservations(
    ObservationStanding Standing,
    ImmutableArray<ProteinGeometryKindObservation> Kinds,
    ImmutableArray<GeometryDistanceObservation> LocatedDistances,
    ImmutableArray<string> Limitations);

public sealed record ApplicablePreparationPolicy(
    string Id,
    string Version,
    string ApplicableMolecularClass,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles,
    int MaximumMinimizationIterations,
    double FinalUnrestrainedRmsForceTargetKjMolNm,
    double ExportCoordinateReadBackToleranceAngstrom,
    double ExportCellLengthReadBackToleranceAngstrom,
    double ExportCellAngleReadBackToleranceDegrees,
    MolecularDynamicsSystemSettings SystemSettings,
    ConstructionPolicy Construction,
    MolecularRepresentation Water,
    MolecularRepresentation Sodium,
    MolecularRepresentation Chloride,
    LocalStateObservationSpec LocalStateObservation,
    ImmutableArray<LocalStateCriterion> ConstructionCriteria,
    ImmutableArray<LocalRolePairCriterion> ContactCriteria,
    ImmutableArray<PreparationAssessmentCriterion> AssessmentCriteria,
    ProteinGeometryMeasurementSpec StageProteinGeometryMeasurement,
    ImmutableArray<StageProteinGeometryCriterion> StageProteinGeometryCriteria,
    EquilibrationProtocol? OptionalEquilibration,
    ImmutableArray<string> Limitations,
    PreparationPolicyScope? Scope = null);

/// <summary>The independently qualified molecular and condition scope of a preparation basis.</summary>
public sealed record PreparationPolicyScope(
    string SourceCoordinateSha256,
    int SourceModelIndex,
    string? BiologicalAssemblyId,
    ImmutableArray<ChainSelection> ChainCopies,
    ImmutableArray<string> RetainedPartnerSourceIds,
    string PreparedBondGraphSha256,
    string ChemicalStatePolicyId,
    string ChemicalStatePolicyVersion,
    string ResidueVariantsSha256,
    LeafletComposition Upper,
    LeafletComposition Lower,
    FixedStudyConditions Conditions,
    ProteinTopologyKind TopologyKind,
    string ProteinStructuralPolicyId = "",
    string ProteinStructuralPolicyVersion = "",
    string MembraneSupportPolicyId = "",
    string MembraneSupportPolicyVersion = "");
public sealed record ForceFieldAsset(
    string Id,
    string Version,
    string Family,
    string Path,
    string Sha256);

/// <summary>Declared force-field system construction choices, not worker defaults.</summary>
public sealed record MolecularDynamicsSystemSettings(
    string NonbondedMethod,
    double NonbondedCutoffNanometers,
    string Constraints,
    bool RigidWater,
    double EwaldErrorTolerance,
    double? SwitchDistanceNanometers,
    bool UseDispersionCorrection,
    bool RemoveCMMotion,
    double? HydrogenMassDaltons);

public sealed record ConstructionPolicy(
    string Id,
    string Version,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<ProteinTopologyKind> CoveredTopologyKinds,
    ImmutableArray<string> CoveredSpeciesIds,
    string ProviderName,
    string ProviderVersion,
    string NativePatchPath,
    string NativePatchSha256,
    string LipidTypeArgument,
    string PositiveIonArgument,
    string NegativeIonArgument,
    double MinimumPaddingNanometers,
    double WaterMolarityForIonRounding,
    int MaximumAtomCount,
    double MaximumCellDimensionAngstrom,
    int MaximumConstructionSeconds,
    string ApproximationStatement,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativePatchMode = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativeSourcePatchPath = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativeSourcePatchSha256 = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] ImmutableArray<string>? RemovedNativeLipidResidueIds = null);

public sealed record PreparationAssessmentCriterion(
    StageKind StageKind,
    string MeasurementName,
    string Unit,
    string Scope,
    double? Minimum,
    double? Maximum);

/// <summary>Stage-specific interpretation of measured intraprotein geometry.</summary>
public sealed record StageProteinGeometryCriterion(
    StageKind StageKind,
    ProteinGeometryCriterion Criterion);

/// <summary>Declared finite local measurements of actual constructed and completed states.</summary>
public sealed record LocalStateObservationSpec(
    ImmutableArray<LocalContactRolePair> ContactRolePairs,
    ImmutableDictionary<string, double> AtomRadiusByElementAngstrom,
    double ContactSearchRadiusAngstrom,
    int MaximumReportedPairs,
    bool UsePeriodicBoundary,
    ImmutableArray<string> RequiredMetricNames);
public sealed record LocalContactRolePair(MoleculeRoleKind FirstMoleculeRole, MoleculeRoleKind SecondMoleculeRole);
public sealed record LocalStateCriterion(
    string MeasurementName,
    string Unit,
    string Scope,
    double? Minimum,
    double? Maximum);
/// <summary>Qualified observed contact rule; null stage kind means construction.</summary>
public sealed record LocalRolePairCriterion(
    StageKind? StageKind,
    MoleculeRoleKind FirstMoleculeRole,
    MoleculeRoleKind SecondMoleculeRole,
    int MinimumPairsWithinSearchRadius,
    double? MinimumNearestDistanceAngstrom,
    double? MaximumNearestDistanceAngstrom);
public sealed record LocalRolePairMeasurement(
    MoleculeRoleKind FirstMoleculeRole,
    MoleculeRoleKind SecondMoleculeRole,
    int PairsWithinSearchRadius,
    double? MinimumDistanceAngstrom);
public sealed record LocalContactObservation(
    int FirstAtomIndex,
    int SecondAtomIndex,
    MoleculeRoleKind FirstMoleculeRole,
    MoleculeRoleKind SecondMoleculeRole,
    double DistanceAngstrom,
    double RadiusSumAngstrom);
public sealed record LocalStateObservations(
    ObservationStanding Standing,
    string? UnavailableReason,
    ImmutableArray<MeasuredValue> Measurements,
    ImmutableArray<LocalContactObservation> LocatedContacts,
    [property: JsonConverter(typeof(CoveredRolePairsJsonConverter))]
    ImmutableArray<LocalContactRolePair> CoveredRolePairs,
    ImmutableArray<string> Limitations,
    ImmutableArray<LocalRolePairMeasurement> RolePairMeasurements);

public sealed record EquilibrationProtocol(
    string Id,
    double TargetTemperatureKelvin,
    int RandomSeed,
    ImmutableArray<EquilibrationStageControl> Stages,
    EquilibrationStageControl ExtensionWindow,
    int MaximumExtensions,
    int MaximumSampleCount,
    ImmutableArray<string> RequiredObservations,
    ImmutableArray<EquilibrationObservable> Observables,
    ImmutableArray<EquilibrationSufficiencyRule> SufficiencyRules,
    string ComparisonBasis,
    long MaximumFrameBytes = 0);

public sealed record EquilibrationFrameSeries(
    string ManifestPath,
    string ManifestSha256,
    string PositionsPath,
    string PositionsSha256,
    long PositionsByteLength,
    int FrameCount);

public sealed record EquilibrationObservable(
    string Name,
    string Unit,
    string Scope,
    string Category,
    string Method,
    string AtomSelector,
    string ComparisonSelector,
    string ThirdSelector,
    double? DistanceCutoffAngstrom);
public sealed record ResolvedEquilibrationObservable(
    string Name,
    string Unit,
    string Scope,
    string Category,
    string Method,
    ImmutableArray<int> AtomIndices,
    ImmutableArray<int> ComparisonAtomIndices,
    ImmutableArray<int> ThirdAtomIndices,
    double? DistanceCutoffAngstrom);
public sealed record EquilibrationSufficiencyRule(
    string ObservableName,
    int BlockSizeSamples,
    int MinimumEffectiveBlocks,
    double MaximumAbsoluteFirstVsLastBlockMeanDifference,
    double MaximumAbsoluteLagOneBlockCorrelation);
public sealed record EquilibrationSample(
    string WindowName,
    int Step,
    ImmutableArray<MeasuredValue> Measurements);
public sealed record EquilibrationObservationAssessment(
    string ObservableName,
    int SampleCount,
    double EffectiveBlockCount,
    double? FirstVsLastBlockMeanDifference,
    double? LagOneBlockCorrelation,
    bool Sufficient);

public sealed record EquilibrationStageControl(
    string Name,
    int Steps,
    double TimestepPicoseconds,
    double TemperatureKelvin,
    double? PressureBar,
    string PressureMode,
    int? BarostatFrequencySteps,
    double? SurfaceTensionBarNm,
    double FrictionPerPicosecond,
    double ProteinRestraintKjMolNm2,
    double LipidRestraintKjMolNm2,
    int ReportIntervalSteps,
    double? InitialTemperatureKelvin = null,
    string ProteinRestraintSelector = "backbone");

public enum AssessmentStanding { Supported, Unsupported, NotEstablished }
public enum StageKind { Minimization, Equilibration }
public enum StageExecutionStanding { Pending, Running, ReadyForMinimization, Completed, Stopped, Failed, ResourceRefused, Unobserved }
public enum PreparationQualification { QualifiedPrepared, NotQualified, Indeterminate }

public sealed record AssessedPreparedProtein(
    string Id,
    string StudyRevisionId,
    IntendedProteinModel Intended,
    MolecularArtifact Molecule,
    string ChemicalStatePolicyId,
    ImmutableArray<ResidueVariantChoice> ResidueVariants,
    ImmutableArray<PreparationChangeProposal> Changes,
    SourceToResultCorrespondence Correspondence,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings,
    ImmutableArray<string> Limitations,
    PredictionEvidenceObservations? Prediction = null,
    ProteinGeometryObservations? Geometry = null,
    string? ChemicalStatePolicyVersion = null,
    string? StructuralAssessmentPolicyId = null,
    string? StructuralAssessmentPolicyVersion = null);

public sealed record LipidFraction(string SpeciesId, double Fraction);

[JsonConverter(typeof(JsonStringEnumConverter<ProteinTopologyKind>))]
public enum ProteinTopologyKind
{
    [JsonStringEnumMemberName("membrane-spanning")] MembraneSpanning = 1,
    [JsonStringEnumMemberName("one-surface-associated")] OneSurfaceAssociated
}

[JsonConverter(typeof(JsonStringEnumConverter<LeafletSide>))]
public enum LeafletSide
{
    [JsonStringEnumMemberName("upper")] Upper = 1,
    [JsonStringEnumMemberName("lower")] Lower
}

[JsonConverter(typeof(JsonStringEnumConverter<PlacementPhysicalSide>))]
public enum PlacementPhysicalSide
{
    [JsonStringEnumMemberName("upper")] Upper = 1,
    [JsonStringEnumMemberName("lower")] Lower,
    [JsonStringEnumMemberName("both")] Both
}

[JsonConverter(typeof(JsonStringEnumConverter<PpmNterminalSide>))]
public enum PpmNterminalSide
{
    [JsonStringEnumMemberName("in")] In = 1,
    [JsonStringEnumMemberName("out")] Out
}

public sealed record LeafletComposition(LeafletSide PhysicalSide, ImmutableArray<LipidFraction> Fractions);
public sealed record MembraneModel(
    string Id,
    LeafletComposition Upper,
    LeafletComposition Lower,
    FixedStudyConditions Conditions,
    string ScientificPurpose);

public sealed record AssessedMembraneModel(
    string Id,
    string StudyRevisionId,
    MembraneModel Intended,
    ImmutableArray<MolecularRepresentation> SpeciesRepresentations,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<string> Limitations,
    string PolicyId = "",
    string PolicyVersion = "");

public sealed record MembraneSupportPolicy(
    string Id,
    string Version,
    ImmutableArray<string> CoveredSpeciesIds,
    bool AllowsMixedLeaflets,
    bool AllowsAsymmetricLeaflets,
    bool CoversArbitraryCoherentFractions,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<string> Limitations);

public sealed record PlacementProposal(
    string Id,
    string PreparedProteinId,
    string MembraneModelId,
    ProteinTopologyKind TopologyKind,
    MolecularArtifact OrientedProtein,
    double? MidplaneAngstrom,
    double? ThicknessAngstrom,
    double? TiltDegrees,
    PlacementPhysicalSide PhysicalSide,
    string? BiologicalSidedness,
    ImmutableArray<string> ContactingRegions,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<string> Limitations);

public sealed record OpmReferenceRecord(
    string PdbAccession,
    string SourceUrl,
    string OrientedCoordinatePath,
    string OrientedCoordinateSha256,
    string? BiologicalAssemblyId,
    ImmutableArray<ChainSelection> ChainCopies,
    ImmutableArray<string> RetainedPartnerSourceIds,
    ImmutableArray<string> SourceAtomIds,
    string MembraneContext,
    double? HydrophobicThicknessAngstrom,
    double? TiltDegrees,
    string? AssumedMembraneSpeciesId = null,
    bool ImplicitSymmetric = false,
    double? MidplaneAngstrom = null);
public sealed record OpmReferenceReview(
    string PreparedProteinId,
    string MembraneModelId,
    bool CorrespondsToSelectedConstruct,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<string> Limitations,
    bool MembraneContextApplicable = false);

public sealed record PlacementSupportPolicy(
    string Id,
    string Version,
    ImmutableArray<ProteinTopologyKind> CoveredTopologyKinds,
    ImmutableArray<string> CoveredSpeciesIds,
    bool AllowsMixtures,
    bool AllowsAsymmetry,
    bool AllowsTransferFromPpmDopc,
    ImmutableArray<string> RequiredEvidenceMethods,
    ImmutableArray<PlacementMeasurementCriterion> GeometryCriteria,
    double InterfaceBandAngstrom,
    PlacementPredictionCriterion? PredictionCriterion,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<string> Limitations,
    PlacementPureLipidCoreFrame? PureLipidCoreFrame = null);
/// <summary>An identified pure-lipid hydrocarbon reference for an intended, still unproved membrane frame.</summary>
public sealed record PlacementPureLipidCoreFrame(
    string Id,
    string Version,
    string SpeciesId,
    double HydrocarbonThicknessAngstrom,
    double ThicknessUncertaintyAngstrom,
    double ReferenceTemperatureKelvin,
    string ReferenceCondition,
    string ReferenceCitation,
    FixedStudyConditions IntendedConditions,
    double MidplaneOffsetFromPpmAngstrom,
    string SourceToTargetReviewReference,
    string SourceToTargetReviewRationale);
/// <summary>An uncertainty gate for identified predicted models, never independent placement proof.</summary>
public sealed record PlacementPredictionCriterion(
    double MinimumLocalConfidence,
    double MinimumLocalCoverageFraction,
    double MaximumDirectionalPaeAngstrom,
    double MinimumPaePairCoverageFraction,
    bool IndependentWitnessCanResolvePredictionLimitations,
    int MaximumReportedPaePairs);
/// <summary>Independently sourced, exact-construct structural interpretation to test against positioned atoms.</summary>
public sealed record PlacementStructuralWitness(
    string Id,
    string Version,
    string SourceCoordinateSha256,
    int SourceModelIndex,
    string? BiologicalAssemblyId,
    ImmutableArray<ChainSelection> ChainCopies,
    ImmutableArray<string> RetainedPartnerSourceIds,
    string? PreparedCoordinateSha256,
    string? PreparedBondGraphSha256,
    LeafletComposition Upper,
    LeafletComposition Lower,
    FixedStudyConditions Conditions,
    ProteinTopologyKind TopologyKind,
    string BiologicalSidedness,
    string Source,
    ImmutableArray<string> EvidenceReferences,
    ImmutableArray<PlacementResidueWitness> Residues,
    ImmutableArray<string> Limitations,
    ImmutableArray<PlacementPredictionResolution> PredictionResolutions = default);
[JsonConverter(typeof(JsonStringEnumConverter<PlacementWitnessRoleKind>))]
public enum PlacementWitnessRoleKind
{
    [JsonStringEnumMemberName("topology")] Topology = 1,
    [JsonStringEnumMemberName("contact")] Contact,
    [JsonStringEnumMemberName("sidedness")] Sidedness
}

public sealed record PlacementResidueWitness(
    ResidueAddress Residue,
    PlacementWitnessRoleKind Role,
    string ExpectedRegion,
    string EvidenceReference);
/// <summary>Independent, exact-region evidence resolving a prediction limitation.</summary>
public sealed record PlacementPredictionResolution(
    ImmutableArray<ResidueAddress> ResolvedResidues,
    bool LocalStructureResolved,
    bool RelativePositionResolved,
    string EvidenceReference);
public sealed record PlacementMeasurementCriterion(
    string MeasurementName,
    string Unit,
    double? Minimum,
    double? Maximum);
public sealed record PlacementMeasurementReport(
    string ProposalId,
    ImmutableArray<MeasuredValue> Measurements,
    ImmutableArray<PlacementResidueObservation> Residues,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<string> Limitations,
    PredictionRegionSummaryObservations? Prediction = null);

public sealed record AssessedProteinMembranePlacement(
    string Id,
    string StudyRevisionId,
    PlacementProposal Proposal,
    AssessmentStanding Standing,
    string Reason,
    ImmutableArray<ScientificFinding> Findings,
    DateTimeOffset AssessedAt);

public sealed record PreparationAttempt(
    string Id,
    string StudyRevisionId,
    string ProteinId,
    string MembraneId,
    string PlacementId,
    string PolicyId,
    DateTimeOffset StartedAt,
    string PolicyVersion,
    string PolicyFingerprintSha256,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles,
    string ConstructionProviderVersion,
    string NativePatchSha256);

public sealed record SpeciesCount(LeafletSide PhysicalSide, string SpeciesId, int Count, double IntendedFraction);
public sealed record ConstructionDerivation(
    string AttemptId,
    ImmutableArray<SpeciesCount> LipidCounts,
    ImmutableArray<double> CellAngstrom,
    int WaterCount,
    int SodiumCount,
    int ChlorideCount,
    double ProteinNetChargeElementary,
    double IntendedNaClMolar,
    double EstimatedNaClMolar,
    double EstimatedAqueousVolumeAngstromCubed,
    ImmutableArray<string> Approximations,
    ImmutableArray<string> Limitations);

public sealed record ConstructedExplicitSystem(
    string Id,
    PreparationAttempt Attempt,
    MolecularArtifact Molecule,
    ConstructionDerivation Derivation,
    SourceToResultCorrespondence Correspondence,
    ImmutableArray<SpeciesCount> AchievedComposition,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings,
    string ConditionsTreatment,
    LocalStateObservations? LocalState = null,
    ImmutableArray<double> ActualCellAngstrom = default);

public sealed record PreparationStartResult(
    PreparationAttempt? Attempt,
    ConstructionDerivation? Derivation,
    ConstructedExplicitSystem? Constructed,
    StageExecutionState State);

[JsonConverter(typeof(JsonStringEnumConverter<EquilibrationObservationAdequacy>))]
public enum EquilibrationObservationAdequacy
{
    [JsonStringEnumMemberName("adequate")] Adequate = 1,
    [JsonStringEnumMemberName("insufficientAtBound")] InsufficientAtBound
}

[JsonConverter(typeof(JsonStringEnumConverter<StageTermination>))]
public enum StageTermination
{
    [JsonStringEnumMemberName("converged")] Converged = 1,
    [JsonStringEnumMemberName("max_iterations")] MaxIterations,
    [JsonStringEnumMemberName("completed")] Completed,
    [JsonStringEnumMemberName("unknown")] Unknown
}

public sealed record StageObservation(
    string StageId,
    string AttemptId,
    StageKind Kind,
    ImmutableArray<MeasuredValue> Measurements,
    ImmutableArray<ScientificEvidence> Evidence,
    StageTermination Termination,
    string ProviderVersion,
    DateTimeOffset ObservedAt,
    EquilibrationObservationAdequacy? ObservationAdequacy = null,
    ImmutableArray<EquilibrationObservationAssessment> EquilibrationAssessments = default,
    ImmutableArray<EquilibrationSample> EquilibrationSamples = default,
    LocalStateObservations? LocalState = null,
    ProteinGeometryObservations? ProteinGeometry = null,
    ImmutableArray<EquilibrationWindowObservation> EquilibrationWindows = default,
    EquilibrationFrameSeries? FrameSeries = null);

public sealed record MeasuredValue(string Name, double Value, string Unit, string Scope);

public sealed record CompletedStage(
    string Id,
    PreparationAttempt Attempt,
    StageKind Kind,
    MolecularArtifact Molecule,
    StageObservation Observation,
    SourceToResultCorrespondence Correspondence,
    string PolicyId,
    string? SourceStageId,
    ImmutableArray<ScientificFinding> Findings,
    DateTimeOffset CompletedAt);

public sealed record StageOperationResult(
    CompletedStage? CompletedStage,
    StageExecutionState State,
    ImmutableArray<ScientificFinding> Findings);

public sealed record StageExecutionState(
    string AttemptId,
    string? StageId,
    StageKind? Kind,
    StageExecutionStanding Standing,
    string Message,
    double? Progress,
    DateTimeOffset UpdatedAt);

public sealed record LocalWorkspaceSnapshot(
    StudyRevision? CurrentStudy,
    PreparationAttempt? CurrentAttempt,
    StageExecutionState? CurrentExecution,
    ImmutableArray<CompletedStage> CompletedStages,
    ImmutableArray<PreparationAssessmentResult> Assessments,
    ImmutableArray<ScientificFinding> Findings,
    DateTimeOffset ObservedAt);

public sealed record PreparationAssessmentResult(
    string Id,
    string StageId,
    PreparationQualification Qualification,
    string Reason,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings,
    ImmutableArray<string> Limitations,
    DateTimeOffset AssessedAt,
    bool CurrentlyApplicable);

public sealed record StructureFocus(string AuthAsymId, int AuthSeqId, string? InsertionCode, int? AuthAtomId);
public sealed record InspectionAnnotation(string Id, string SubjectPartId, string Label, string Meaning, string? EvidenceId, StructureFocus? GeometryFocus);
public sealed record InspectionMetric(string Name, string Value, string? Unit, string SubjectPartId, string? EvidenceId);
public sealed record InspectionAccount(
    string SubjectId,
    string? StructureUrl,
    string RepresentationKind,
    ImmutableArray<string> OmittedMolecules,
    string? FocusId,
    StructureFocus? Focus,
    ImmutableArray<InspectionAnnotation> Annotations,
    ImmutableArray<InspectionMetric> Metrics,
    string StudyRevisionId,
    long StudyRevisionNumber,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings,
    PreparationAssessmentResult? Assessment);

/// <summary>A verified coordinate-row identity for the currently selected inspection subject.</summary>
public sealed record InspectionAtomAccount(
    string SubjectId,
    string StudyRevisionId,
    string StructureToken,
    int AtomSiteIndex,
    AtomCorrespondence Atom);

public sealed record InspectionSubject(
    string Id,
    string StudyRevisionId,
    string? StructureUrl,
    string RepresentationKind,
    ImmutableArray<string> OmittedMolecules,
    ImmutableArray<ScientificEvidence> Evidence,
    ImmutableArray<ScientificFinding> Findings,
    ImmutableArray<InspectionAnnotation> Annotations,
    ImmutableArray<InspectionMetric> Metrics,
    PreparationAssessmentResult? Assessment);

public sealed record CompletedStageBundle(
    string StageId,
    string StudyRevisionId,
    string AttemptId,
    string AssessmentId,
    string BundlePath,
    string Sha256,
    long ByteLength,
    DateTimeOffset DeliveredAt);

/// <summary>Delivery standing for one stage and one currently applicable assessment.</summary>
public sealed record StageExportAccount(
    string StageId,
    string AssessmentId,
    string Status,
    string? Reason,
    string? Sha256,
    long? ByteLength);

public enum ExportDeliveryStanding { Available, Absent, Corrupt, IdentityChanged }
public sealed record ExportDeliveryObservation(
    ExportDeliveryStanding Standing,
    string? Reason,
    byte[]? Bytes,
    string? Sha256,
    long? ByteLength);

// Actor-facing account. JSON uses lower-camel property names through the host's
// serializer options; the browser does not infer scientific standing from shape.
public sealed record WorkspaceState(
    long Revision,
    StudyAccount? Study,
    ImmutableArray<SourceCandidateAccount> SourceCandidates,
    ImmutableArray<SourceModelObservation> SourceModels,
    ImmutableArray<LipidCatalogueAccount> AvailableLipids,
    ProteinAccount? Protein,
    MembraneAccount? Membrane,
    PlacementAccount? Placement,
    AttemptAccount? Attempt,
    ImmutableArray<StageAccount> Stages,
    InspectionAccount? Inspection,
    ImmutableArray<AvailableAction> Actions,
    ImmutableArray<WorkspaceNotice> Notices,
    PredictionEvidenceAccount? SourcePrediction = null);

public sealed record StudyAccount(
    string Id,
    long Number,
    string Summary,
    string? SelectedSourceId,
    int? ModelIndex,
    string? BiologicalAssemblyId,
    ImmutableArray<string> ChainIds,
    ImmutableArray<PartnerSelection> Partners,
    ImmutableArray<AlternateLocationChoice> AlternateLocations,
    FixedStudyConditions Conditions,
    SourceRouteKind? SelectedSourceKind = null,
    UploadOriginKind? UploadProvenance = null,
    string? UploadProvenanceNote = null,
    string? AdoptedPlacementProposalId = null);

public sealed record SourceCandidateAccount(
    string Id,
    string Label,
    SourceRouteKind Kind,
    string Provenance,
    ImmutableArray<string> Limitations);

public sealed record LipidCatalogueAccount(
    string SpeciesId,
    string DisplayName,
    string ChemistryId,
    ImmutableArray<string> Limitations);

/// <summary>Actor-visible prediction scope; local PAE files and axis mapping remain private.</summary>
public sealed record PredictionEvidenceAccount(
    string RecordId,
    ImmutableArray<PredictedResidueConfidence> LocalConfidence,
    PredictionObservationStanding PaeStanding,
    string? PaeReason,
    int? PaeAxisResidueCount,
    ImmutableArray<string> Limitations);

public sealed record ProteinAccount(
    string SubjectId,
    string Status,
    string Summary,
    int? AtomCount,
    ImmutableArray<PreparationChangeProposal> Changes,
    ImmutableArray<ScientificFinding> Findings,
    PredictionEvidenceAccount? Prediction = null,
    ProteinGeometryObservations? Geometry = null,
    ProteinGeometryObservations? SourceGeometry = null,
    string? CandidateId = null);

public sealed record MembraneAccount(
    string ModelId,
    string Status,
    string ScientificPurpose,
    ImmutableArray<LipidFraction> Upper,
    ImmutableArray<LipidFraction> Lower,
    ImmutableArray<string> Limitations,
    string? Reason = null,
    string? PolicyId = null,
    string? PolicyVersion = null,
    ImmutableArray<ScientificEvidence> Evidence = default,
    ImmutableArray<MembraneSpeciesSupportAccount> SpeciesSupport = default);

public sealed record MembraneSpeciesSupportAccount(
    string SpeciesId,
    string ChemistryId,
    string Category,
    string ForceFieldFamily,
    string ForceFieldVersion,
    string CoordinateSha256,
    string ParameterSha256,
    ImmutableArray<string> Limitations);

public sealed record PlacementAccount(
    string ProposalId,
    string Status,
    ProteinTopologyKind TopologyKind,
    double? DepthAngstrom,
    double? TiltDegrees,
    string? Sidedness,
    string Reason,
    ImmutableArray<ScientificEvidence> Evidence,
    PredictionRegionSummaryObservations? Prediction = null,
    string PreparedProteinId = "",
    string MembraneModelId = "",
    PlacementPhysicalSide PhysicalSide = default,
    double? MidplaneAngstrom = null,
    double? ThicknessAngstrom = null,
    ImmutableArray<string> ContactingRegions = default,
    ImmutableArray<string> Limitations = default,
    string? PolicyId = null,
    string? PolicyVersion = null,
    string? WitnessId = null);

public sealed record AttemptAccount(
    string AttemptId,
    string Status,
    StageKind? StageKind,
    double? Progress,
    string Message,
    string? StudyRevisionId = null,
    string? PolicyId = null,
    string? PolicyVersion = null,
    string? CurrentStageId = null,
    ConstructionDerivation? Derivation = null,
    ConstructedSystemAccount? Constructed = null);

public sealed record ConstructedSystemAccount(
    string SubjectId,
    string AttemptId,
    int AtomCount,
    ImmutableArray<SpeciesCount> AchievedComposition,
    ImmutableArray<double> CellAngstrom,
    int WaterCount,
    int SodiumCount,
    int ChlorideCount,
    string ConditionsTreatment,
    LocalStateObservations? LocalState);

public sealed record StageAccount(
    string StageId,
    string AttemptId,
    string StudyRevisionId,
    StageKind Kind,
    string Status,
    PreparationAssessmentResult? Assessment,
    string Summary,
    StageObservation? Observation = null,
    ConstructedSystemAccount? Constructed = null,
    StageExportAccount? Export = null,
    string? SourceStageId = null);

[JsonConverter(typeof(JsonStringEnumConverter<ActorActionKind>))]
public enum ActorActionKind
{
    [JsonStringEnumMemberName("searchSource")] SearchSource = 1,
    [JsonStringEnumMemberName("selectSource")] SelectSource,
    [JsonStringEnumMemberName("selectProteinModel")] SelectProteinModel,
    [JsonStringEnumMemberName("approvePreparationChange")] ApprovePreparationChange,
    // Availability for the negative decision; the actual command is ApprovePreparationChange with approve=false.
    [JsonStringEnumMemberName("declinePreparationChange")] DeclinePreparationChange,
    [JsonStringEnumMemberName("proposeMembrane")] ProposeMembrane,
    [JsonStringEnumMemberName("adoptMembrane")] AdoptMembrane,
    [JsonStringEnumMemberName("proposePlacement")] ProposePlacement,
    [JsonStringEnumMemberName("revisePlacement")] RevisePlacement,
    [JsonStringEnumMemberName("adoptPlacement")] AdoptPlacement,
    [JsonStringEnumMemberName("startPreparation")] StartPreparation,
    [JsonStringEnumMemberName("continueMinimization")] ContinueMinimization,
    [JsonStringEnumMemberName("stopAttempt")] StopAttempt,
    [JsonStringEnumMemberName("requestEquilibration")] RequestEquilibration,
    [JsonStringEnumMemberName("selectInspectionSubject")] SelectInspectionSubject,
    [JsonStringEnumMemberName("setInspectionFocus")] SetInspectionFocus,
    [JsonStringEnumMemberName("exportStage")] ExportStage
}

public sealed record AvailableAction(ActorActionKind Kind, string? SubjectId, bool Enabled, string? Reason);
public sealed record WorkspaceNotice(string Id, string Severity, string Message, string? SubjectId);
public sealed record ActorCommand(ActorActionKind Kind, JsonElement Data, long ExpectedRevision);

public sealed record BoundaryOutcome<T>(T? Value, string Reason, ImmutableArray<ScientificFinding> Findings,
    IObservedDiagnostic? Diagnostic = null)
    where T : class
{
    public bool Established => Value is not null;
    public static BoundaryOutcome<T> Success(T value) => new(value, string.Empty, ImmutableArray<ScientificFinding>.Empty);
    public static BoundaryOutcome<T> Unavailable(string reason, ImmutableArray<ScientificFinding> findings = default,
        IObservedDiagnostic? diagnostic = null)
        => new(null, reason, findings.IsDefault ? ImmutableArray<ScientificFinding>.Empty : findings, diagnostic);
}

// Owner-facing capabilities share one local scientific-worker process. Each
// semantic owner receives only its operations; the physical adapter maps them
// to named four-field line-JSON envelopes:
// requestId, operation, workingDirectory and payload. Attempt and stage identity
// live in the typed payload and are echoed in the terminal result. The adapter
// must confine reported artifacts to the requested working directory, verify
// their hashes, and report lost terminal communication as Unobserved.
public interface IProteinPreparationWork
{
    Task<WorkerResult<SourceInspectionObservations>> InspectSourceAsync(
        ScientificWorkRequest<SourceInspectionPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<PreparationChangeObservations>> InspectPreparationChangesAsync(
        ScientificWorkRequest<PreparationChangeInspectionPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<ProteinPreparationObservations>> PrepareProteinAsync(
        ScientificWorkRequest<ProteinPreparationPayload> request, CancellationToken cancellationToken);
}

public interface IMembraneModelAssessmentWork
{
    Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
        ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken);
}

public interface IPlacementAssessmentWork
{
    Task<WorkerResult<PredictionRegionSummaryObservations>> SummarizePredictionEvidenceAsync(
        ScientificWorkRequest<PredictionRegionSummaryPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<PlacementObservations>> PlacePpmAsync(
        ScientificWorkRequest<PlacementPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<PlacementAdjustmentObservations>> AdjustPlacementAsync(
        ScientificWorkRequest<PlacementAdjustmentPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<PlacementMeasurementObservations>> MeasurePlacementAsync(
        ScientificWorkRequest<PlacementMeasurementPayload> request, CancellationToken cancellationToken);
}

public interface IExplicitConstructionWork
{
    Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
        ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken);
}

public interface IMinimizationWork
{
    Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
        ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken);
    Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken);
}

public interface IOptionalEquilibrationWork
{
    Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
        ScientificWorkRequest<EquilibrationPayload> request,
        IProgress<EquilibrationWorkProgress>? progress,
        CancellationToken cancellationToken);
    Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken);
}

public interface ICompletedStageExportWork
{
    Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
        ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken);
}

// Composition alone holds the complete local process capability. Semantic
// children receive only the owner-facing interface they are entitled to use.
public interface IScientificWorkerExchange : IProteinPreparationWork,
    IMembraneModelAssessmentWork, IPlacementAssessmentWork, IExplicitConstructionWork,
    IMinimizationWork, IOptionalEquilibrationWork, ICompletedStageExportWork
{
}

public sealed record ScientificWorkRequest<TPayload>(string RequestId, string WorkingDirectory, TPayload Payload)
    where TPayload : class;

public enum WorkerResultStanding { Observed, Failed, Stopped, Unobserved }
public sealed record WorkerArtifact(string Role, string Path, string Sha256);
public sealed record ProviderIdentity(string Name, string Version);
public sealed record WorkerResult<TObservation>(
    string RequestId,
    string? StudyRevisionId,
    string? AttemptId,
    string? StageId,
    WorkerResultStanding Standing,
    ImmutableArray<WorkerArtifact> Artifacts,
    TObservation? Observations,
    ProviderIdentity? Provider,
    string? FailureCode,
    string? FailureMessage)
    where TObservation : class;

public sealed record SourceInspectionPayload(
    string SourcePath,
    string SourceSha256,
    int MaxAtoms,
    SourceRouteKind SourceKind,
    PredictionEvidenceAsset? Prediction);
public sealed record PredictedResidueConfidence(
    ResidueAddress Residue,
    double? PLddt,
    PredictionObservationStanding Standing,
    string? Reason);
public sealed record PredictionEvidenceObservations(
    string RecordId,
    string CoordinateSha256,
    ImmutableArray<PredictedResidueConfidence> LocalConfidence,
    PredictionObservationStanding PaeStanding,
    string? PaeReason,
    int? PaeAxisResidueCount,
    string? PaeMappingPath,
    string? PaeMappingSha256,
    ImmutableArray<string> Limitations);
public sealed record PredictionRegionSummaryPayload(
    string StudyRevisionId,
    string RecordId,
    string CoordinateSha256,
    string PaePath,
    string PaeSha256,
    string PaeMappingPath,
    string PaeMappingSha256,
    ImmutableArray<ResidueAddress> FirstRegion,
    ImmutableArray<ResidueAddress> SecondRegion,
    int MaximumReportedPairs);
public sealed record PredictionPairObservation(
    ResidueAddress First,
    ResidueAddress Second,
    double PredictedAlignedErrorAngstrom);
public sealed record DirectionalPredictionSummary(
    int PossiblePairCount,
    int ValidPairCount,
    double? MinimumAngstrom,
    double? MaximumAngstrom,
    double? MeanAngstrom,
    ImmutableArray<PredictionPairObservation> HighestErrorPairs);
public sealed record PredictionRegionSummaryObservations(
    string RecordId,
    PredictionObservationStanding Standing,
    string? Reason,
    DirectionalPredictionSummary? FirstAlignedOnSecond,
    DirectionalPredictionSummary? SecondAlignedOnFirst);
public sealed record SourceChainObservation(string Name, int ResidueCount, int AtomCount);
public sealed record SourceAssemblyObservation(string Name, ImmutableArray<ChainSelection> ChainCopies);
public sealed record SourcePartnerObservation(
    string SourceId,
    string Label,
    string Kind,
    int AtomCount,
    string? Chain,
    int? Residue);

[JsonConverter(typeof(JsonStringEnumConverter<SourceResidueKind>))]
public enum SourceResidueKind
{
    [JsonStringEnumMemberName("protein")] Protein = 1,
    [JsonStringEnumMemberName("solvent")] Solvent,
    [JsonStringEnumMemberName("heterogen")] Heterogen
}

public sealed record SourceResidueObservation(
    ResidueAddress Address,
    string Name,
    SourceResidueKind ResidueKind,
    bool BackboneHeavyAtomsComplete,
    ImmutableArray<string> AlternateLocations,
    ImmutableArray<string> MissingNonbackboneHeavyAtomNames,
    ResidueAddress? PossibleDisulfidePartner,
    ObservationStanding MissingAtomAssessmentStanding,
    ObservationStanding DisulfideAssessmentStanding,
    ImmutableArray<string> Limitations);
public sealed record SourceModelObservation(
    int Index,
    ImmutableArray<SourceChainObservation> Chains,
    ImmutableArray<SourceAssemblyObservation> Assemblies,
    ImmutableArray<SourcePartnerObservation> Partners,
    ImmutableArray<SourceResidueObservation> Residues,
    int AtomCount);
public sealed record SourceInspectionObservations(
    string SourceFormat,
    ImmutableArray<SourceModelObservation> Models,
    PredictionEvidenceObservations? Prediction);

public sealed record PreparationChangeInspectionPayload(
    string StudyRevisionId,
    double DisulfideCandidateMaxSgDistanceAngstrom,
    string SourcePath,
    string SourceSha256,
    int ModelIndex,
    string? AssemblyId,
    ImmutableArray<ChainSelection> ChainSelections,
    ImmutableArray<AlternateLocationChoice> AltlocChoices,
    ImmutableArray<SourcePartnerObservation> RetainedPartners,
    ProteinGeometryMeasurementSpec GeometrySpec);
public sealed record PossibleDisulfideObservation(
    ResidueAddress First,
    ResidueAddress Second,
    double DistanceAngstrom);
public sealed record PreparationChangeObservations(
    ImmutableArray<AtomAddress> MissingNonbackboneHeavyAtoms,
    ImmutableArray<PossibleDisulfideObservation> PossibleDisulfides,
    ObservationStanding AssessmentStanding,
    ImmutableArray<string> Limitations,
    int SelectedAtomCount,
    ImmutableArray<PreviewChainCorrespondence> PreviewChains,
    ProteinGeometryObservations Geometry);

public sealed record MembraneAssessmentPayload(
    string StudyRevisionId,
    string MembraneModelId,
    LeafletComposition Upper,
    LeafletComposition Lower,
    ImmutableArray<MolecularRepresentation> SpeciesRepresentations,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles);
public sealed record SpeciesTemplateObservation(
    string SpeciesId,
    string ChemistryId,
    int CoordinateAtomCount,
    int ParameterAtomCount,
    bool AtomIdentityAndBondMatch,
    ImmutableArray<string> Warnings);
public sealed record MembraneAssessmentObservations(
    ImmutableArray<SpeciesTemplateObservation> Species,
    ImmutableArray<string> CombinationWarnings,
    bool? CombinedParameterizationObserved);

public sealed record ProteinPreparationPayload(
    string StudyRevisionId,
    double NominalPh,
    double DisulfideCandidateMaxSgDistanceAngstrom,
    string SourcePath,
    string SourceSha256,
    int ModelIndex,
    string? AssemblyId,
    ImmutableArray<ChainSelection> ChainSelections,
    ImmutableArray<AlternateLocationChoice> AltlocChoices,
    ImmutableArray<SourcePartnerObservation> RetainedPartners,
    ImmutableArray<HeavyAtomApproval> ApprovedHeavyAtoms,
    ImmutableArray<DisulfideChoice> ApprovedDisulfides,
    ImmutableArray<ResidueVariantChoice> ResidueVariants,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles,
    ProteinGeometryMeasurementSpec GeometrySpec);

public sealed record ProteinPreparationObservations(
    int SourceAtomCount,
    int PreparedAtomCount,
    int RetainedResidueCount,
    ImmutableArray<ResidueAddress> MissingBackboneResidues,
    ImmutableArray<ResidueAddress> NoncanonicalResidues,
    ImmutableArray<ResidueAddress> UnresolvedAlternateLocations,
    ImmutableArray<ResidueAddress> UnparameterizedResidues,
    ImmutableArray<AtomAddress> AddedHeavyAtoms,
    ImmutableArray<AtomAddress> AddedHydrogens,
    ImmutableArray<AtomAddress> RemovedSourceHydrogens,
    ImmutableArray<ResidueVariantChoice> ActualResidueVariants,
    ImmutableArray<DisulfideBond> ActualDisulfides,
    ImmutableArray<string> GeometryWarnings,
    int CorrespondedResultAtomCount,
    ProteinGeometryObservations Geometry);

public sealed record PlacementPayload(
    string StudyRevisionId,
    string PreparedProteinId,
    string PreparedPdbPath,
    string PreparedSha256,
    string PpmExecutablePath,
    string PpmVersion,
    string PpmExecutableSha256,
    ProteinTopologyKind TopologyKind,
    PpmNterminalSide PpmNterminalSide,
    string PpmResidueLibraryPath = "",
    string PpmResidueLibrarySha256 = "");

public sealed record PlacementObservations(
    double? MidplaneAngstrom,
    double? ThicknessAngstrom,
    double? TiltDegrees,
    int AlignedSourceAtomCount,
    ImmutableArray<MeasuredValue> NumericalOutput,
    ImmutableArray<string> PlaneMarkerIds,
    ImmutableArray<string> InterpretationWarnings,
    string AssumedMembrane);

public sealed record PlacementAdjustmentPayload(
    string StudyRevisionId,
    string SourceProposalId,
    string OrientedPdbPath,
    double DepthShiftAngstrom,
    double TiltAboutXDegrees,
    double TiltAboutYDegrees,
    double RotationAboutNormalDegrees,
    string Rationale,
    string OrientedPdbSha256 = "");
public sealed record PlacementAdjustmentObservations(
    int SourceAtomCount,
    int AdjustedAtomCount,
    double AppliedDepthShiftAngstrom,
    double AppliedTiltAboutXDegrees,
    double AppliedTiltAboutYDegrees,
    double AppliedRotationAboutNormalDegrees,
    ImmutableArray<string> GeometryWarnings);
public sealed record PlacementMeasurementPayload(
    string StudyRevisionId,
    string ProposalId,
    string OrientedPdbPath,
    double MembraneMidplaneAngstrom,
    double CoreLowerZAngstrom,
    double CoreUpperZAngstrom,
    ImmutableArray<ResidueAddress> ResidueAddressesInOrder,
    ImmutableArray<string> OutputChainIdsInOrder = default,
    ImmutableArray<string> ExpectedResultAtomIdsInOrder = default,
    string OrientedPdbSha256 = "");
public sealed record PlacementResidueObservation(
    ResidueAddress Address,
    string OutputChainId,
    string OutputResidueId,
    string OutputInsertionCode,
    string Name,
    int AtomCount,
    double MinZAngstrom,
    double MaxZAngstrom,
    double MeanZAngstrom,
    int AtomsWithinCore,
    int BackboneAtomsWithinCore,
    int AtomsAboveCore,
    int AtomsBelowCore);
public sealed record PlacementMeasurementObservations(
    int AtomCount,
    ImmutableArray<PlacementResidueObservation> Residues,
    int AtomsWithinCore,
    int AtomsAboveCore,
    int AtomsBelowCore,
    double ProteinZMinAngstrom,
    double ProteinZMaxAngstrom,
    ImmutableArray<string> Limitations);

public sealed record ConstructionPayload(
    string StudyRevisionId,
    string AttemptId,
    string OrientedPdbPath,
    string OrientedPdbSha256,
    string PreparedPdbPath,
    string PreparedPdbSha256,
    SourceToResultCorrespondence PreparedCorrespondence,
    string PreparedBondGraphPath,
    string PreparedBondGraphSha256,
    MolecularRepresentation Lipid,
    MolecularRepresentation Water,
    MolecularRepresentation Sodium,
    MolecularRepresentation Chloride,
    string NativePatchPath,
    string NativePatchSha256,
    string ProviderName,
    string ProviderVersion,
    string LipidTypeArgument,
    string PositiveIonArgument,
    string NegativeIonArgument,
    double MembraneCenterZNanometers,
    double MinimumPaddingNanometers,
    double IonicStrengthMolar,
    ImmutableArray<ForceFieldAsset> ForceFieldFiles,
    MolecularDynamicsSystemSettings SystemSettings,
    LocalStateObservationSpec LocalObservationSpec,
    int MaximumAtomCount,
    double MaximumCellDimensionAngstrom,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativePatchMode = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativeSourcePatchPath = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] string? NativeSourcePatchSha256 = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] ImmutableArray<string>? RemovedNativeLipidResidueIds = null);

public sealed record ObservedSpeciesCount(GeneratedComponentRoleKind Role, LeafletSide PhysicalSide, string SpeciesId, int Count);
public sealed record ConstructionObservations(
    int AtomCount,
    ImmutableArray<ObservedSpeciesCount> SpeciesCounts,
    ImmutableArray<double> ActualCellAngstrom,
    ImmutableArray<double> ProteinPeriodicImageGapsAngstrom,
    int WaterCount,
    int PositiveIonCount,
    int NegativeIonCount,
    double NetChargeElementary,
    double ProteinNetChargeElementary,
    double InitialPotentialEnergyKjMol,
    int CorrespondedResultAtomCount,
    double MaximumProteinCoordinateDeviationAngstrom,
    bool ProteinIdentityAndBondsPreserved,
    string NativePatchSha256,
    ImmutableArray<string> ContactWarnings,
    ImmutableArray<string> GeometryWarnings,
    ImmutableArray<string> ParameterWarnings,
    LocalStateObservations LocalState,
    string NativePatchMode = "installed",
    string? NativeSourcePatchSha256 = null);

public sealed record MinimizationPayload(
    string StudyRevisionId,
    string AttemptId,
    string StageId,
    string TopologyCifPath,
    string TopologyCifSha256,
    string TopologyJsonPath,
    string TopologyJsonSha256,
    string SystemXmlPath,
    string SystemXmlSha256,
    string StateXmlPath,
    string StateXmlSha256,
    int MaxIterations,
    double RmsForceTargetKjMolNm);

public sealed record MinimizationObservations(
    double InitialPotentialEnergyKjMol,
    double FinalPotentialEnergyKjMol,
    double FinalRmsForceKjMolNm,
    int? Iterations,
    StageTermination Termination,
    bool FinalTreatmentUnrestrained,
    int FinalAtomCount,
    ImmutableArray<string> NumericalWarnings,
    double? FinalRawRmsForceKjMolNm = null,
    double? MaximumRelativeConstraintError = null,
    double? AppliedConstraintTolerance = null,
    string? FinalRmsForceMethod = null);

public sealed record EquilibrationPayload(
    string StudyRevisionId,
    string AttemptId,
    string StageId,
    string SourceMinimizedStageId,
    string TopologyCifPath,
    string TopologyCifSha256,
    string TopologyJsonPath,
    string TopologyJsonSha256,
    string SystemXmlPath,
    string SystemXmlSha256,
    string MinimizedStateXmlPath,
    string MinimizedStateXmlSha256,
    string ValidatedPolicyId,
    ImmutableArray<int> ProteinBackboneAtomIndices,
    ImmutableArray<int> ProteinHeavyAtomIndices,
    ImmutableArray<int> LipidHeavyAtomIndices,
    ImmutableArray<ResolvedEquilibrationObservable> ResolvedObservables,
    EquilibrationProtocol Protocol,
    string ProtocolSha256 = "");

public sealed record EquilibrationWorkProgress(
    string StudyRevisionId,
    string AttemptId,
    string StageId,
    string Window,
    int CompletedSteps,
    int RequestedSteps);

public sealed record EquilibrationWindowObservation(
    string Name,
    int RequestedSteps,
    int CompletedSteps,
    double? TemperatureKelvin,
    double? PressureBar,
    ImmutableArray<MeasuredValue> Measurements,
    ImmutableArray<string> Warnings);

public sealed record EquilibrationObservations(
    ImmutableArray<EquilibrationWindowObservation> Windows,
    int FinalAtomCount,
    StageTermination Termination,
    bool FinalObservationUnrestrained,
    ImmutableArray<string> NumericalWarnings,
    int ExtensionsPerformed,
    EquilibrationObservationAdequacy ObservationAdequacy,
    ImmutableArray<EquilibrationSample> Samples,
    ImmutableArray<EquilibrationObservationAssessment> ObservationAssessments);

public sealed record StageObservationPayload(
    string StudyRevisionId,
    string AttemptId,
    string StageId,
    string TopologyCifPath,
    string TopologyCifSha256,
    string TopologyJsonPath,
    string TopologyJsonSha256,
    string SystemXmlPath,
    string SystemXmlSha256,
    string StateXmlPath,
    string StateXmlSha256,
    StageKind StageKind,
    string CorrespondencePath,
    string CorrespondenceSha256,
    LocalStateObservationSpec LocalObservationSpec,
    ProteinGeometryMeasurementSpec StageProteinGeometrySpec);

public sealed record StageObservationObservations(
    int AtomCount,
    bool AtomOrderMatched,
    bool BondsMatched,
    ImmutableArray<MeasuredValue> Measurements,
    ImmutableArray<string> ContactWarnings,
    ImmutableArray<string> StructuralWarnings,
    ImmutableArray<string> NumericalWarnings,
    LocalStateObservations LocalState,
    ProteinGeometryObservations ProteinGeometry);

public sealed record ExportVerificationPayload(
    string StudyRevisionId,
    string AttemptId,
    string StageId,
    string TopologyCifPath,
    string TopologyCifSha256,
    string TopologyJsonPath,
    string TopologyJsonSha256,
    string SystemXmlPath,
    string SystemXmlSha256,
    string StateXmlPath,
    string StateXmlSha256,
    double CoordinateReadBackToleranceAngstrom,
    double CellLengthReadBackToleranceAngstrom,
    double CellAngleReadBackToleranceDegrees);
public sealed record ExportVerificationObservations(
    int SourceAtomCount,
    int ExportedAtomCount,
    int CorrespondingElementAndResidueCount,
    bool AtomOrderMatched,
    bool BondsMatched,
    bool CellMatched,
    bool ReadBackMatched,
    double CoordinateMaxDeviationAngstrom,
    ImmutableArray<string> Warnings);

/// <summary>Stable identity of the complete attempt-bound preparation basis.</summary>
public static class PreparationPolicyFingerprint
{
    public static string Compute(ApplicablePreparationPolicy policy)
    {
        ArgumentNullException.ThrowIfNull(policy);
        // Optional stereo descriptors may be absent in older water/ion records. Treat that
        // representation and an explicit empty array identically in the policy identity.
        var normalized = policy with
        {
            Water = Normalize(policy.Water),
            Sodium = Normalize(policy.Sodium),
            Chloride = Normalize(policy.Chloride)
        };
        var encoded = JsonSerializer.SerializeToElement(normalized, new JsonSerializerOptions(JsonSerializerDefaults.Web));
        return ComputeCanonical(encoded);
    }

    private static MolecularRepresentation Normalize(MolecularRepresentation representation) =>
        representation is null ? null! :
        representation with { StereoChecks = representation.StereoChecks.IsDefault
            ? ImmutableArray<MolecularStereoCheck>.Empty : representation.StereoChecks };

    public static string ComputeResidueVariants(ImmutableArray<ResidueVariantChoice> variants)
    {
        if (variants.IsDefault || variants.Any(item => item.Residue is null ||
                string.IsNullOrWhiteSpace(item.Variant)) ||
            variants.Select(item => item.Residue).Distinct().Count() != variants.Length)
            return string.Empty;
        var ordered = variants.OrderBy(item => item.Residue.Model)
            .ThenBy(item => item.Residue.Chain, StringComparer.Ordinal)
            .ThenBy(item => item.Residue.Residue)
            .ThenBy(item => item.Residue.InsertionCode, StringComparer.Ordinal)
            .ThenBy(item => item.Residue.CopyId, StringComparer.Ordinal)
            .Select(item => new
            {
                item.Residue.Model, item.Residue.Chain, item.Residue.Residue,
                item.Residue.InsertionCode, item.Residue.CopyId, item.Variant
            });
        return ComputeCanonical(JsonSerializer.SerializeToElement(ordered,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)));
    }

    private static string ComputeCanonical(JsonElement encoded)
    {
        var bytes = new ArrayBufferWriter<byte>();
        using (var writer = new Utf8JsonWriter(bytes))
            WriteCanonical(writer, encoded);
        return Convert.ToHexString(SHA256.HashData(bytes.WrittenSpan)).ToLowerInvariant();
    }

    public static bool Matches(PreparationAttempt attempt, ApplicablePreparationPolicy policy) =>
        attempt.PolicyId == policy.Id &&
        !string.IsNullOrWhiteSpace(attempt.PolicyVersion) &&
        attempt.PolicyVersion == policy.Version &&
        attempt.PolicyFingerprintSha256.Length == 64 &&
        attempt.PolicyFingerprintSha256.All(Uri.IsHexDigit) &&
        attempt.PolicyFingerprintSha256.Equals(Compute(policy), StringComparison.OrdinalIgnoreCase) &&
        attempt.ConstructionProviderVersion == policy.Construction.ProviderVersion &&
        attempt.NativePatchSha256.Equals(policy.Construction.NativePatchSha256,
            StringComparison.OrdinalIgnoreCase) &&
        !attempt.ForceFieldFiles.IsDefault && !policy.ForceFieldFiles.IsDefault &&
        attempt.ForceFieldFiles.SequenceEqual(policy.ForceFieldFiles);

    private static void WriteCanonical(Utf8JsonWriter writer, JsonElement element)
    {
        if (element.ValueKind == JsonValueKind.Object)
        {
            writer.WriteStartObject();
            foreach (var property in element.EnumerateObject().OrderBy(item => item.Name, StringComparer.Ordinal))
            {
                writer.WritePropertyName(property.Name);
                WriteCanonical(writer, property.Value);
            }
            writer.WriteEndObject();
            return;
        }
        if (element.ValueKind == JsonValueKind.Array)
        {
            writer.WriteStartArray();
            foreach (var child in element.EnumerateArray()) WriteCanonical(writer, child);
            writer.WriteEndArray();
            return;
        }
        element.WriteTo(writer);
    }
}
