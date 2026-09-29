using System.Collections.Immutable;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ConstructionOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation;
using ExportOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport.CompletedStageExport;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class CompletedStageExportOwnerTests
{
    [Fact]
    public async Task Checked_plan_heavy_atom_marker_needs_exact_plan_and_approved_change()
    {
        using var basis = new ConstructionFixture();
        var assessed = await AssessmentFixture.Create(basis);
        var residue = Assert.IsType<ResidueAddress>(assessed.Protein.Correspondence.Atoms[0].SourceResidue);
        const string preparationRevisionId = "earlier-preparation-revision";
        var change = new PreparationChangeProposal("oxt-change", preparationRevisionId,
            assessed.Protein.Intended.Id, residue, PreparationChangeKind.HeavyAtom, "OXT",
            "Controlled planned terminal atom", ImmutableArray<string>.Empty, true);
        var plan = new PreparationPlanProvenance(new string('a', 64), "controlled ff19SB",
            "1", 7, 1, assessed.Protein.Molecule.CoordinateSha256,
            ImmutableArray<RecommendedStateChoice>.Empty,
            ImmutableArray.Create(new AtomAddress(residue, "OXT")),
            ImmutableArray<AtomAddress>.Empty);
        var preparedAtom = assessed.Protein.Correspondence.Atoms[0] with
        {
            ResultAtomId = "0:A:1::OXT", SourceAtomId = null,
            Role = AtomOriginKind.Generated, MoleculeRole = MoleculeRoleKind.Protein,
            Element = "O", ApprovedChangeId = "recommendation-plan-only"
        };
        var protein = assessed.Protein with
        {
            PreparationStudyRevisionId = preparationRevisionId,
            Changes = ImmutableArray.Create(change), RecommendationPlan = plan,
            Correspondence = assessed.Protein.Correspondence with
            { Atoms = assessed.Protein.Correspondence.Atoms.SetItem(0, preparedAtom) }
        };
        Assert.NotEqual(protein.StudyRevisionId, protein.PreparationStudyRevisionId);
        var constructedAtom = assessed.Constructed.Correspondence.Atoms[0] with
        {
            SourceAtomId = null, SourceResidue = residue, Role = AtomOriginKind.Generated,
            MoleculeRole = MoleculeRoleKind.Protein, Element = "O",
            ApprovedChangeId = "recommendation-plan-only"
        };
        var atoms = assessed.Constructed.Correspondence.Atoms.SetItem(0, constructedAtom);
        var constructed = assessed.Constructed with
        { Correspondence = assessed.Constructed.Correspondence with { Atoms = atoms } };
        var topologyPath = Path.Combine(basis.Directory, "plan-marker-stage-topology.json");
        File.WriteAllText(topologyPath, JsonSerializer.Serialize(new
        {
            formatVersion = 1,
            chains = new[] { new { id = "A" } },
            residues = new[] { new { chainIndex = 0, id = residue.Residue.ToString(),
                name = "GLU", insertionCode = "" } },
            atoms = Enumerable.Range(0, assessed.Stage.Molecule.AtomCount).Select(index => new
            {
                residueIndex = 0, id = (index + 1).ToString(),
                name = index == 0 ? "OXT" : "C", element = index == 0 ? "O" : "C"
            }).ToArray(),
            bonds = Array.Empty<object>(),
            boxVectorsAngstrom = new[] { new[] { 80.0, 0.0, 0.0 },
                new[] { 0.0, 81.0, 0.0 }, new[] { 0.0, 0.0, 100.0 } }
        }));
        constructed = constructed with { Molecule = constructed.Molecule with
        { TopologyPath = topologyPath, TopologySha256 = ConstructionFixture.Hash(topologyPath) } };
        var stage = assessed.Stage with
        {
            Molecule = assessed.Stage.Molecule with
            {
                CoordinateSha256 = ConstructionFixture.Hash(assessed.Stage.Molecule.CoordinatePath!),
                TopologyPath = topologyPath, TopologySha256 = ConstructionFixture.Hash(topologyPath)
            },
            Observation = assessed.Stage.Observation with
            {
                EquilibrationAssessments = ImmutableArray<EquilibrationObservationAssessment>.Empty,
                EquilibrationSamples = ImmutableArray<EquilibrationSample>.Empty,
                EquilibrationWindows = ImmutableArray<EquilibrationWindowObservation>.Empty
            },
            Correspondence = assessed.Stage.Correspondence with { Atoms = atoms }
        };
        var assessment = (assessed with { Stage = stage, Constructed = constructed,
            Protein = protein }).Assess();
        Assert.Equal(stage.Correspondence.Atoms.Length, stage.Molecule.AtomCount);
        Assert.True(stage.Correspondence.Atoms.SequenceEqual(constructed.Correspondence.Atoms));
        Assert.Equal(stage.Correspondence.SourceId, constructed.Correspondence.SourceId);
        Assert.Equal(stage.Correspondence.ResultId, stage.Molecule.Id);
        Assert.True(assessment.CurrentlyApplicable);
        Assert.Equal(PreparationPolicyFingerprint.Compute(assessed.Policy),
            stage.Attempt.PolicyFingerprintSha256);
        var approved = ImmutableArray.Create(new ResearcherDecision("approved-oxt",
            preparationRevisionId, change.Id, ResearcherDecisionKind.ApprovePreparationChange,
            ResearcherDecisionValue.Approved, DateTimeOffset.UtcNow, null));
        var worker = new ControlledExportReadBack(constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);

        async Task<BoundaryOutcome<CompletedStageBundle>> Export(
            CompletedStage selected, AssessedPreparedProtein selectedProtein,
            ImmutableArray<ResearcherDecision> selectedDecisions, string folder,
            ConstructedExplicitSystem? selectedConstructed = null) =>
            await owner.ExportAsync(selected, assessment, basis.Revision,
                selectedProtein, assessed.Membrane, assessed.Placement,
                selectedConstructed ?? constructed, constructed.Derivation, assessed.Policy,
                selectedDecisions, null, null, null,
                Path.Combine(basis.Directory, folder), TestContext.Current.CancellationToken);

        var accepted = await Export(stage, protein, approved, "plan-marker-accepted");
        Assert.True(accepted.Established, accepted.Reason);
        Assert.Single(worker.Requests);

        // Memgen can return a retained atom at another final index. Its checked
        // provenance multiset and hash-verified final atom name still identify OXT.
        var movedSource = assessed.Protein.Correspondence.Atoms[1] with
        { ResultAtomIndex = 0, ResultAtomId = "0:A:2::CA" };
        var movedOxt = preparedAtom with { ResultAtomIndex = 1,
            ResultAtomId = "1:A:1::OXT" };
        var reorderedProtein = protein with { Correspondence = protein.Correspondence with
        { Atoms = protein.Correspondence.Atoms.SetItem(0, movedSource)
            .SetItem(1, movedOxt) } };
        Assert.True((await Export(stage, reorderedProtein, approved,
            "plan-marker-reordered-retained")).Established);
        var wrongSourceAtom = constructedAtom with
        { SourceResidue = residue with { Residue = residue.Residue + 1 } };
        var wrongSourceAtoms = atoms.SetItem(0, wrongSourceAtom);
        var wrongSourceConstructed = constructed with { Correspondence =
            constructed.Correspondence with { Atoms = wrongSourceAtoms } };
        var wrongSourceStage = stage with { Correspondence =
            stage.Correspondence with { Atoms = wrongSourceAtoms } };
        Assert.False((await Export(wrongSourceStage, protein, approved,
            "plan-marker-wrong-retained-identity", wrongSourceConstructed)).Established);

        var wrongProtein = protein with { Correspondence = protein.Correspondence with
        { Atoms = protein.Correspondence.Atoms.SetItem(0, preparedAtom with
            { ApprovedChangeId = "recommendation-plan-unchecked" }) } };
        var wrong = await Export(stage, wrongProtein, approved, "plan-marker-wrong");
        Assert.False(wrong.Established);
        Assert.Contains("approved preparation change", wrong.Reason ?? "", StringComparison.Ordinal);
        var unapproved = await Export(stage, protein, ImmutableArray<ResearcherDecision>.Empty,
            "plan-marker-unapproved");
        Assert.False(unapproved.Established);
        Assert.Contains("approved preparation change", unapproved.Reason ?? "", StringComparison.Ordinal);
        var foreignChange = change with { StudyRevisionId = "foreign-preparation-revision" };
        var foreignProtein = protein with { Changes = ImmutableArray.Create(foreignChange) };
        var foreignApproval = ImmutableArray.Create(approved[0] with
        { StudyRevisionId = foreignChange.StudyRevisionId });
        Assert.False((await Export(stage, foreignProtein, foreignApproval,
            "plan-marker-foreign-revision-pair")).Established);

        var staleCandidate = protein with { RecommendationPlan = plan with
        { CheckedCandidateSha256 = new string('b', 64) } };
        Assert.False((await Export(stage, staleCandidate, approved,
            "plan-marker-stale-candidate")).Established);
        var wrongProposal = protein with { Changes = ImmutableArray.Create(change with
            { ProposedChange = "OD1" }) };
        Assert.False((await Export(stage, wrongProposal, approved,
            "plan-marker-wrong-proposal")).Established);
        var wrongIntended = protein with { Changes = ImmutableArray.Create(change with
            { IntendedProteinId = "other-intended-protein" }) };
        Assert.False((await Export(stage, wrongIntended, approved,
            "plan-marker-wrong-intended")).Established);
        var malformedDigest = protein with { RecommendationPlan = plan with { PlanSha256 = "stale" } };
        Assert.False((await Export(stage, malformedDigest, approved,
            "plan-marker-malformed-digest")).Established);

        var wrongNameTopology = Path.Combine(basis.Directory, "plan-marker-wrong-name.json");
        File.WriteAllText(wrongNameTopology, File.ReadAllText(topologyPath).Replace(
            "\"name\":\"OXT\"", "\"name\":\"OD1\"", StringComparison.Ordinal));
        Assert.Contains("\"name\":\"OD1\"", File.ReadAllText(wrongNameTopology),
            StringComparison.Ordinal);
        var wrongNameStage = stage with { Molecule = stage.Molecule with
        { TopologyPath = wrongNameTopology,
            TopologySha256 = ConstructionFixture.Hash(wrongNameTopology) } };
        var wrongNameConstructed = constructed with { Molecule = constructed.Molecule with
        { TopologyPath = wrongNameTopology,
            TopologySha256 = ConstructionFixture.Hash(wrongNameTopology) } };
        Assert.False((await Export(wrongNameStage, protein, approved,
            "plan-marker-wrong-final-name", wrongNameConstructed)).Established);

        // Provider chain and residue labels can differ from the source after
        // mapping, while the checked retained atom index remains identical.
        var remappedTopology = Path.Combine(basis.Directory, "plan-marker-remapped.json");
        File.WriteAllText(remappedTopology, File.ReadAllText(topologyPath).Replace(
            "\"chains\":[{\"id\":\"A\"}]", "\"chains\":[{\"id\":\"B\"}]",
            StringComparison.Ordinal).Replace("\"id\":\"1\",\"name\":\"GLU\"",
                "\"id\":\"999\",\"name\":\"GLU\"", StringComparison.Ordinal));
        Assert.Contains("\"chains\":[{\"id\":\"B\"}]",
            File.ReadAllText(remappedTopology), StringComparison.Ordinal);
        var remappedStage = stage with { Molecule = stage.Molecule with
        { TopologyPath = remappedTopology,
            TopologySha256 = ConstructionFixture.Hash(remappedTopology) } };
        var remappedConstructed = constructed with { Molecule = constructed.Molecule with
        { TopologyPath = remappedTopology,
            TopologySha256 = ConstructionFixture.Hash(remappedTopology) } };
        Assert.True((await Export(remappedStage, protein, approved,
            "plan-marker-remapped-labels", remappedConstructed)).Established);

        var nullMarkerAtoms = atoms.SetItem(0, constructedAtom with { ApprovedChangeId = null });
        var nullMarkerConstructed = constructed with { Correspondence =
            constructed.Correspondence with { Atoms = nullMarkerAtoms } };
        var nullMarkerStage = stage with { Correspondence = stage.Correspondence with
        { Atoms = nullMarkerAtoms } };
        Assert.False((await Export(nullMarkerStage, protein, approved,
            "plan-marker-null-heavy", nullMarkerConstructed)).Established);

        // Manual changes use an individual decision, but must bind to the
        // preparation revision even when this prepared protein was carried forward.
        const string manualDecisionId = "manual-oxt-approval";
        var manualAtom = constructedAtom with { ApprovedChangeId = manualDecisionId };
        var manualAtoms = atoms.SetItem(0, manualAtom);
        var manualConstructed = constructed with { Correspondence =
            constructed.Correspondence with { Atoms = manualAtoms } };
        var manualStage = stage with { Correspondence = stage.Correspondence with
        { Atoms = manualAtoms } };
        var manualProtein = protein with
        {
            RecommendationPlan = null,
            Correspondence = protein.Correspondence with { Atoms =
                protein.Correspondence.Atoms.SetItem(0, preparedAtom with
                { ApprovedChangeId = manualDecisionId }) }
        };
        var manualApproval = ImmutableArray.Create(new ResearcherDecision(manualDecisionId,
            preparationRevisionId, change.Id, ResearcherDecisionKind.ApprovePreparationChange,
            ResearcherDecisionValue.Approved, DateTimeOffset.UtcNow, null));
        var manualAccepted = await Export(manualStage, manualProtein, manualApproval,
            "manual-marker-accepted", manualConstructed);
        Assert.True(manualAccepted.Established, manualAccepted.Reason);
        var wrongManualNameStage = manualStage with { Molecule = wrongNameStage.Molecule };
        Assert.False((await Export(wrongManualNameStage, manualProtein, manualApproval,
            "manual-marker-wrong-final-name", manualConstructed with
            { Molecule = wrongNameConstructed.Molecule })).Established);
        var wrongManualResidue = manualProtein with { Changes = ImmutableArray.Create(change with
        { Residue = residue with { Residue = residue.Residue + 1 } }) };
        Assert.False((await Export(manualStage, wrongManualResidue, manualApproval,
            "manual-marker-wrong-approved-residue", manualConstructed)).Established);
        var foreignManual = manualProtein with { Changes = ImmutableArray.Create(foreignChange) };
        var foreignManualApproval = ImmutableArray.Create(manualApproval[0] with
        { StudyRevisionId = foreignChange.StudyRevisionId });
        Assert.False((await Export(manualStage, foreignManual, foreignManualApproval,
            "manual-marker-foreign-revision-pair", manualConstructed)).Established);
        Assert.Equal(4, worker.Requests.Count);
    }

    [Fact]
    public async Task Established_technical_assessment_is_exported_without_a_new_scientific_judgment()
    {
        using var basis = new ConstructionFixture();
        var assessed = await AssessmentFixture.Create(basis);
        var stage = assessed.Stage with
        {
            Molecule = assessed.Stage.Molecule with
            { CoordinateSha256 = ConstructionFixture.Hash(assessed.Stage.Molecule.CoordinatePath!) },
            Observation = assessed.Stage.Observation with
            {
                EquilibrationAssessments = ImmutableArray<EquilibrationObservationAssessment>.Empty,
                EquilibrationSamples = ImmutableArray<EquilibrationSample>.Empty,
                EquilibrationWindows = ImmutableArray<EquilibrationWindowObservation>.Empty
            }
        };
        var assessment = assessed.Assess(stage: stage);
        Assert.Equal(PreparationCheckStanding.ChecksPassed, assessment.CheckStanding);
        Assert.True(assessment.CurrentlyApplicable);

        var worker = new ControlledExportReadBack(assessed.Constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);
        var result = await owner.ExportAsync(stage, assessment, basis.Revision,
            assessed.Protein, assessed.Membrane, assessed.Placement,
            assessed.Constructed, assessed.Constructed.Derivation, assessed.Policy,
            ImmutableArray<ResearcherDecision>.Empty, null, null, null,
            Path.Combine(basis.Directory, "checked-export"), TestContext.Current.CancellationToken);

        Assert.True(result.Established, result.Reason);
        Assert.Single(worker.Requests);
        Assert.Equal(stage.Id, worker.Requests[0].Payload.StageId);
        using var archive = ZipFile.OpenRead(result.Value!.BundlePath);
        using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
        var declared = manifest.RootElement.GetProperty("assessment");
        Assert.Equal(assessment.Id, declared.GetProperty("id").GetString());
        Assert.Equal(stage.Id, declared.GetProperty("stageId").GetString());
        Assert.Equal("checksPassed", declared.GetProperty("checkStanding").GetString());
        Assert.Equal(assessment.Reason, declared.GetProperty("reason").GetString());
        var membrane = manifest.RootElement.GetProperty("membrane").GetProperty("intended");
        Assert.Equal(assessed.Membrane.Intended.ScientificPurpose,
            membrane.GetProperty("scientificPurpose").GetString());
        Assert.Equal("DMPC", membrane.GetProperty("upper").GetProperty("fractions")[0]
            .GetProperty("speciesId").GetString());
        Assert.Equal(assessed.Protein.Changes.Length,
            manifest.RootElement.GetProperty("preparedProtein").GetProperty("changes").GetArrayLength());
        Assert.Equal(stage.Id, manifest.RootElement.GetProperty("lineage").GetProperty("sourceToResult")
            .GetProperty("resultId").GetString());

        var withoutPurpose = assessed.Membrane with
        {
            Intended = assessed.Membrane.Intended with { ScientificPurpose = null }
        };
        var newResult = await owner.ExportAsync(stage, assessment, basis.Revision,
            assessed.Protein, withoutPurpose, assessed.Placement,
            assessed.Constructed, assessed.Constructed.Derivation, assessed.Policy,
            ImmutableArray<ResearcherDecision>.Empty, null, null, null,
            Path.Combine(basis.Directory, "no-purpose-export"), TestContext.Current.CancellationToken);
        Assert.True(newResult.Established, newResult.Reason);
        using var newArchive = ZipFile.OpenRead(newResult.Value!.BundlePath);
        using var newManifest = JsonDocument.Parse(newArchive.GetEntry("manifest.json")!.Open());
        Assert.Equal(JsonValueKind.Null, newManifest.RootElement.GetProperty("membrane")
            .GetProperty("intended").GetProperty("scientificPurpose").ValueKind);
        Assert.Equal("checksPassed", newManifest.RootElement.GetProperty("assessment")
            .GetProperty("checkStanding").GetString());

        var numberedIntended = assessed.Protein.Intended with { SourceModelId = "5" };
        var numberedProtein = assessed.Protein with { Intended = numberedIntended };
        var numberedRevision = basis.Revision with { IntendedProtein = numberedIntended };
        var numberedResult = await owner.ExportAsync(stage, assessment, numberedRevision,
            numberedProtein, assessed.Membrane, assessed.Placement,
            assessed.Constructed, assessed.Constructed.Derivation, assessed.Policy,
            ImmutableArray<ResearcherDecision>.Empty, null, null, null,
            Path.Combine(basis.Directory, "numbered-model-export"), TestContext.Current.CancellationToken);
        Assert.True(numberedResult.Established, numberedResult.Reason);
        using var numberedArchive = ZipFile.OpenRead(numberedResult.Value!.BundlePath);
        using var numberedManifest = JsonDocument.Parse(numberedArchive.GetEntry("manifest.json")!.Open());
        Assert.Equal("5", numberedManifest.RootElement.GetProperty("preparedProtein")
            .GetProperty("sourceModelId").GetString());
        Assert.Equal(5, numberedManifest.RootElement.GetProperty("lineage")
            .GetProperty("sourceModelNumber").GetInt32());
    }

    [Theory]
    [InlineData("changedAmberDigest")]
    [InlineData("changedSourceAtomMap")]
    [InlineData("changedAuthorizedMethod")]
    [InlineData("changedImportPrmtopDigest")]
    [InlineData("changedImportRestartDigest")]
    [InlineData("changedParameterStanding")]
    public async Task Memgen_export_refuses_changed_construction_provenance_before_delivery(
        string mutation)
    {
        using var basis = new ConstructionFixture();
        var assessed = await AssessmentFixture.Create(basis);
        var amberTopology = Path.Combine(basis.Directory, "selected-trial.top");
        var amberRestart = Path.Combine(basis.Directory, "selected-trial.restrt");
        File.WriteAllText(amberTopology, "controlled selected Amber topology bytes");
        File.WriteAllText(amberRestart, "controlled selected conditioned restart bytes");
        var asset = new ProviderAsset("controlled-memgen", "2026.3.25",
            basis.NativePatchPath, basis.NativePatchSha);
        var settings = new MemgenConstructionSettings("sander", "ff19SB", "lipid21", "tip3p",
            15, 17.5, 23, 20, 100, 20, 2, 250, 250, 10, 2,
            true, true, true, true, true, true, true, true, true, true, true);
        var policy = assessed.Policy with { Construction = assessed.Policy.Construction with
        {
            Route = ConstructionRouteKind.PackmolMemgen,
            SaltConvention = SaltConventionKind.MemgenChargeCompensated,
            ProviderName = "PACKMOL-Memgen", ProviderVersion = "2026.3.25",
            NativePatchPath = null, NativePatchSha256 = null, LipidTypeArgument = null,
            ProviderAssets = [asset], Memgen = settings
        } };
        var attempt = assessed.Stage.Attempt with
        {
            Route = ConstructionRouteKind.PackmolMemgen,
            SaltConvention = SaltConventionKind.MemgenChargeCompensated,
            NativePatchSha256 = null,
            ConstructionProviderVersion = policy.Construction.ProviderVersion,
            ProviderAssets = [asset],
            PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy)
        };
        var trial = new ConstructionTrialSummary("selected-trial", 0,
            ConstructionTrialStanding.Checked, 15, 17.5,
            assessed.Constructed.AchievedComposition, assessed.Constructed.AchievedComposition,
            ImmutableArray<SpeciesCount>.Empty, [80, 81, 100], [80, 81, 100],
            null, null, null, [
                new TrialDiagnosticArtifact("amberTopology", ConstructionFixture.Hash(amberTopology),
                    Path.GetFileName(amberTopology), PreparationPhase.AmberParameterization,
                    LocalPath: amberTopology),
                new TrialDiagnosticArtifact("amberFinalRestart", ConstructionFixture.Hash(amberRestart),
                    Path.GetFileName(amberRestart), PreparationPhase.ProviderConditioningUnrestrained,
                    LocalPath: amberRestart)
            ]);
        // This controlled premise supplies an already checked import account;
        // the owner test verifies its exact binding, not parameter science.
        var parameterProof = Path.Combine(basis.Directory, "selected-parameter-comparison.json");
        File.WriteAllText(parameterProof, JsonSerializer.Serialize(new
        {
            kind = "controlled Amber-to-OpenMM parameter comparison",
            atomCount = assessed.Constructed.Molecule.AtomCount,
            selectedPrmtopSha256 = ConstructionFixture.Hash(amberTopology),
            selectedRestartSha256 = ConstructionFixture.Hash(amberRestart),
            atomOrder = true, bonded = true, nonbonded = true, exclusions = true, units = true
        }));
        var import = new AmberImportAccount(ConstructionFixture.Hash(amberTopology),
            ConstructionFixture.Hash(amberRestart), assessed.Constructed.Molecule.AtomCount,
            assessed.Constructed.Molecule.AtomCount, 0, 0, 0, 0, 0,
            ["HarmonicBondForce", "NonbondedForce"], true, true, true, true, true, true,
            ConstructionFixture.Hash(parameterProof));
        var derivation = assessed.Constructed.Derivation with
        { Trials = [trial], SelectedTrialId = trial.TrialId, AmberImport = import };
        var constructed = assessed.Constructed with { Attempt = attempt, Derivation = derivation,
            MaximumProteinCoordinateDeviationAngstrom = 2.125 };
        var correspondencePath = Path.Combine(basis.Directory, "selected-trial-correspondence.json");
        File.WriteAllText(correspondencePath, JsonSerializer.Serialize(
            constructed.Correspondence, new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        var correspondenceSha = ConstructionFixture.Hash(correspondencePath);
        constructed = constructed with { Molecule = constructed.Molecule with
        { CorrespondencePath = correspondencePath, CorrespondenceSha256 = correspondenceSha } };
        var stage = assessed.Stage with
        {
            Attempt = attempt,
            Molecule = assessed.Stage.Molecule with
            {
                CoordinateSha256 = ConstructionFixture.Hash(assessed.Stage.Molecule.CoordinatePath!),
                CorrespondencePath = correspondencePath, CorrespondenceSha256 = correspondenceSha
            },
            Observation = assessed.Stage.Observation with
            {
                EquilibrationAssessments = ImmutableArray<EquilibrationObservationAssessment>.Empty,
                EquilibrationSamples = ImmutableArray<EquilibrationSample>.Empty,
                EquilibrationWindows = ImmutableArray<EquilibrationWindowObservation>.Empty
            }
        };
        var assessment = new ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment
            .PreparationAssessment().Assess(stage, constructed, assessed.Protein,
                assessed.Membrane, assessed.Placement, policy,
                ImmutableArray<ScientificFinding>.Empty);
        Assert.True(assessment.CurrentlyApplicable, assessment.Reason);
        var worker = new ControlledExportReadBack(constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);

        Task<BoundaryOutcome<CompletedStageBundle>> Export(CompletedStage selectedStage,
            ConstructedExplicitSystem selectedConstructed, string suffix) =>
            owner.ExportAsync(selectedStage, assessment, basis.Revision, assessed.Protein,
                assessed.Membrane, assessed.Placement, selectedConstructed,
                selectedConstructed.Derivation, policy,
                ImmutableArray<ResearcherDecision>.Empty, null, null, null,
                Path.Combine(basis.Directory, "memgen-export-" + suffix),
                TestContext.Current.CancellationToken);

        var baseline = await Export(stage, constructed, "unchanged-" + mutation);
        Assert.True(baseline.Established, baseline.Reason);
        Assert.Single(worker.Requests);
        Assert.True(File.Exists(baseline.Value!.BundlePath));
        using (var archive = ZipFile.OpenRead(baseline.Value.BundlePath))
        using (var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open()))
        {
            Assert.Equal(2.125, manifest.RootElement.GetProperty("construction")
                .GetProperty("maximumProteinCoordinateDeviationAngstrom").GetDouble());
            var recordedImport = manifest.RootElement.GetProperty("construction")
                .GetProperty("derivation").GetProperty("amberImport");
            Assert.Equal(import.PrmtopSha256,
                recordedImport.GetProperty("prmtopSha256").GetString());
            Assert.Equal(import.FinalRestartSha256,
                recordedImport.GetProperty("finalRestartSha256").GetString());
            Assert.Equal(import.ParameterCorrespondenceSha256,
                recordedImport.GetProperty("parameterCorrespondenceSha256").GetString());
            Assert.True(recordedImport.GetProperty("parameterComparisonPerformed").GetBoolean());
            Assert.Equal(import.AmberAtomCount,
                recordedImport.GetProperty("amberAtomCount").GetInt32());
        }

        if (mutation == "changedAmberDigest")
        {
            var altered = trial with { DiagnosticArtifacts = trial.DiagnosticArtifacts.SetItem(0,
                trial.DiagnosticArtifacts[0] with { Sha256 = new string('f', 64) }) };
            constructed = constructed with { Derivation = derivation with { Trials = [altered] } };
        }
        else if (mutation == "changedSourceAtomMap")
        {
            var atoms = constructed.Correspondence.Atoms;
            var altered = atoms.SetItem(0, atoms[0] with { SourceAtomId = atoms[1].SourceAtomId })
                .SetItem(1, atoms[1] with { SourceAtomId = atoms[0].SourceAtomId });
            constructed = constructed with { Correspondence = constructed.Correspondence with
            { Atoms = altered } };
            stage = stage with { Correspondence = stage.Correspondence with { Atoms = altered } };
        }
        else if (mutation == "changedAuthorizedMethod")
        {
            var altered = attempt with { Route = ConstructionRouteKind.NativeOpenMm };
            constructed = constructed with { Attempt = altered };
            stage = stage with { Attempt = altered };
        }
        else
        {
            var alteredImport = mutation switch
            {
                "changedImportPrmtopDigest" => import with { PrmtopSha256 = new string('e', 64) },
                "changedImportRestartDigest" => import with { FinalRestartSha256 = new string('d', 64) },
                "changedParameterStanding" => import with { ParameterComparisonPerformed = false },
                _ => throw new InvalidOperationException($"Unknown controlled mutation: {mutation}")
            };
            constructed = constructed with { Derivation = derivation with
            { AmberImport = alteredImport } };
        }

        var refused = await Export(stage, constructed, mutation);
        Assert.False(refused.Established,
            $"{mutation} delivered a Memgen bundle with changed construction provenance.");
        Assert.Single(worker.Requests);
        var refusedDirectory = Path.Combine(basis.Directory, "memgen-export-" + mutation);
        Assert.False(Directory.Exists(refusedDirectory) &&
            Directory.EnumerateFiles(refusedDirectory, "*.zip").Any());
    }

    [Fact]
    public async Task Write_destination_failure_publishes_nothing_and_retry_exports_same_stage()
    {
        if (!OperatingSystem.IsLinux()) return; // POSIX directory permissions are the fault injector.
        using var fixture = await ExportFixture.CreateAsync();
        var worker = new ControlledExportReadBack(fixture.Constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);
        var output = Path.Combine(fixture.Directory, "exports",
            fixture.Minimized.Id + "-" + fixture.MinimizedAssessment.Id);
        // The controlled worker first writes a valid verification artifact,
        // then makes only its output directory unwritable. The owner reaches
        // the ZIP write and must refuse publication without changing science.
        worker.AfterReadBack = directory =>
        {
            if (OperatingSystem.IsLinux()) File.SetUnixFileMode(directory,
                UnixFileMode.UserRead | UnixFileMode.UserExecute);
        };
        BoundaryOutcome<CompletedStageBundle> failed;
        try
        {
            failed = await fixture.ExportAsync(owner, fixture.Minimized,
                fixture.MinimizedAssessment, fixture.Policy);
        }
        finally
        {
            if (Directory.Exists(output)) File.SetUnixFileMode(output,
                UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);
        }
        Assert.False(failed.Established);
        Assert.Single(worker.Requests); // Verification crossed before the write fault.
        Assert.Contains("could not be delivered", failed.Reason ?? "", StringComparison.Ordinal);
        Assert.Empty(Directory.GetFiles(output, "*.zip"));
        Assert.Empty(Directory.GetFiles(output, "*.tmp"));
        Assert.Equal(fixture.Minimized.Id, fixture.MinimizedAssessment.StageId);

        var retried = await fixture.ExportAsync(owner, fixture.Minimized,
            fixture.MinimizedAssessment, fixture.Policy);
        Assert.True(retried.Established, retried.Reason);
        Assert.Equal(fixture.Minimized.Id, retried.Value!.StageId);
        Assert.Equal(fixture.MinimizedAssessment.Id, retried.Value.AssessmentId);
        Assert.Equal(2, worker.Requests.Count);
        Assert.True(File.Exists(retried.Value.BundlePath));
    }

    [Fact]
    public async Task Later_completed_stage_exports_its_own_cell_protocol_and_assessment_without_changing_minimized_bundle()
    {
        using var fixture = await ExportFixture.CreateAsync();
        var worker = new ControlledExportReadBack(fixture.Constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);

        var minimized = await fixture.ExportAsync(owner, fixture.Minimized,
            fixture.MinimizedAssessment, fixture.Policy);
        Assert.True(minimized.Established, minimized.Reason);
        var minimizedBytes = File.ReadAllBytes(minimized.Value!.BundlePath);

        var equilibrated = await fixture.ExportAsync(owner, fixture.Equilibrated,
            fixture.EquilibratedAssessment, fixture.Policy);
        Assert.True(equilibrated.Established, equilibrated.Reason);
        Assert.NotEqual(minimized.Value.BundlePath, equilibrated.Value!.BundlePath);
        Assert.Equal(minimizedBytes, File.ReadAllBytes(minimized.Value.BundlePath));
        Assert.Equal(2, worker.Requests.Count);
        Assert.Equal(fixture.Minimized.Id, worker.Requests[0].Payload.StageId);
        Assert.Equal(fixture.Equilibrated.Id, worker.Requests[1].Payload.StageId);
        Assert.Equal(fixture.Equilibrated.Molecule.TopologySha256,
            worker.Requests[1].Payload.TopologyJsonSha256);
        Assert.NotEqual(fixture.Minimized.Molecule.TopologySha256,
            fixture.Equilibrated.Molecule.TopologySha256);

        using var archive = ZipFile.OpenRead(equilibrated.Value.BundlePath);
        foreach (var (entryName, sourcePath) in new[]
                 {
                     ("structure.cif", fixture.Equilibrated.Molecule.CoordinatePath!),
                     ("topology.json", fixture.Equilibrated.Molecule.TopologyPath!),
                     ("system.xml", fixture.Equilibrated.Molecule.SystemXmlPath!),
                     ("state.xml", fixture.Equilibrated.Molecule.StateXmlPath!)
                 })
        {
            using var entryStream = archive.GetEntry(entryName)!.Open();
            using var bytes = new MemoryStream();
            await entryStream.CopyToAsync(bytes, TestContext.Current.CancellationToken);
            Assert.Equal(File.ReadAllBytes(sourcePath), bytes.ToArray());
        }
        using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
        var root = manifest.RootElement;
        Assert.Equal(fixture.Equilibrated.Id,
            root.GetProperty("stage").GetProperty("id").GetString());
        Assert.Equal(fixture.Minimized.Id,
            root.GetProperty("stage").GetProperty("sourceStageId").GetString());
        Assert.Equal(fixture.EquilibratedAssessment.Id,
            root.GetProperty("assessment").GetProperty("id").GetString());
        Assert.Equal("checksIncomplete",
            root.GetProperty("assessment").GetProperty("checkStanding").GetString());
        var optional = root.GetProperty("optionalEquilibration");
        Assert.Equal(fixture.Minimized.Id,
            optional.GetProperty("sourceMinimizedStageId").GetString());
        Assert.Equal(fixture.Protocol.Id,
            optional.GetProperty("protocolId").GetString());
        Assert.Equal(EquilibrationProtocolFingerprint.Compute(fixture.Protocol),
            optional.GetProperty("protocolFingerprintSha256").GetString());
        Assert.Equal(fixture.Protocol.Stages[0].Steps,
            optional.GetProperty("declaredProtocol").GetProperty("stages")[0]
                .GetProperty("steps").GetInt32());
        Assert.Equal(fixture.Protocol.Stages[0].PressureBar,
            optional.GetProperty("declaredProtocol").GetProperty("stages")[0]
                .GetProperty("pressureBar").GetDouble());
        Assert.Equal(fixture.Protocol.Stages[0].Steps,
            optional.GetProperty("observedProcedure").GetProperty("windows")[0]
                .GetProperty("completedSteps").GetInt32());
        Assert.Equal("insufficientAtBound", root.GetProperty("observation")
            .GetProperty("observationAdequacy").GetString());
        var cell = root.GetProperty("finalCell");
        Assert.Equal(81.5, cell.GetProperty("boxVectorsAngstrom")[0][0].GetDouble());
        Assert.Equal(80.5, cell.GetProperty("boxVectorsAngstrom")[1][1].GetDouble());
        Assert.Equal(104, cell.GetProperty("boxVectorsAngstrom")[2][2].GetDouble());
        Assert.Equal(new[] { 81.5, 80.5, 104.0 }, cell.GetProperty("lengthsAngstrom")
            .EnumerateArray().Select(item => item.GetDouble()).ToArray());
        Assert.All(cell.GetProperty("anglesDegrees").EnumerateArray(),
            angle => Assert.Equal(90, angle.GetDouble(), 10));
        Assert.Equal(fixture.Equilibrated.Molecule.TopologySha256,
            root.GetProperty("molecularIdentity").GetProperty("topologySha256").GetString());

        using var earlier = ZipFile.OpenRead(minimized.Value.BundlePath);
        using var earlierManifest = JsonDocument.Parse(earlier.GetEntry("manifest.json")!.Open());
        Assert.Equal(fixture.Minimized.Id,
            earlierManifest.RootElement.GetProperty("stage").GetProperty("id").GetString());
        Assert.Equal(JsonValueKind.Null,
            earlierManifest.RootElement.GetProperty("optionalEquilibration").ValueKind);
        Assert.Equal(JsonValueKind.Null,
            earlierManifest.RootElement.GetProperty("finalCell").ValueKind);
    }

    [Fact]
    public async Task Equilibrated_export_refuses_missing_lineage_protocol_or_stage_specific_cell_before_worker()
    {
        using var fixture = await ExportFixture.CreateAsync();
        var worker = new ControlledExportReadBack(fixture.Constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);
        async Task Refused(CompletedStage stage, ApplicablePreparationPolicy? policy = null)
        {
            var before = worker.Requests.Count;
            var result = await fixture.ExportAsync(owner, stage,
                fixture.EquilibratedAssessment, policy ?? fixture.Policy);
            Assert.False(result.Established);
            Assert.Equal(before, worker.Requests.Count);
        }

        await Refused(fixture.Equilibrated with { SourceStageId = null });
        await Refused(fixture.Equilibrated with { SourceStageId = fixture.Equilibrated.Id });
        await Refused(fixture.Equilibrated with { Observation = fixture.Equilibrated.Observation with
            { Termination = StageTermination.Unknown } });
        await Refused(fixture.Equilibrated with { Observation = fixture.Equilibrated.Observation with
            { EquilibrationWindows = ImmutableArray<EquilibrationWindowObservation>.Empty } });
        await Refused(fixture.Equilibrated, fixture.Policy with { OptionalEquilibration = null });
        await Refused(fixture.Equilibrated, fixture.Policy with
            { OptionalEquilibration = fixture.Protocol with { RandomSeed = 99 } });
        await Refused(fixture.Equilibrated with { Molecule = fixture.Equilibrated.Molecule with
            { TopologyPath = Path.Combine(fixture.Directory, "missing.json") } });

        File.WriteAllText(fixture.Equilibrated.Molecule.TopologyPath!,
            "{\"boxVectorsAngstrom\":[[82,0,0],[0,80.5,0],[0,0,104]]}");
        await Refused(fixture.Equilibrated);
        var changed = fixture.Equilibrated with { Molecule = fixture.Equilibrated.Molecule with
            { TopologySha256 = ConstructionFixture.Hash(fixture.Equilibrated.Molecule.TopologyPath!) } };
        var accepted = await fixture.ExportAsync(owner, changed,
            fixture.EquilibratedAssessment, fixture.Policy);
        Assert.True(accepted.Established, accepted.Reason);
        Assert.Single(worker.Requests);
        using (var archive = ZipFile.OpenRead(accepted.Value!.BundlePath))
        using (var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open()))
            Assert.Equal(82, manifest.RootElement.GetProperty("finalCell")
                .GetProperty("lengthsAngstrom")[0].GetDouble());

        foreach (var invalidCell in new[]
                 {
                     "{}",
                     "{\"boxVectorsAngstrom\":[[0,0,0],[0,80.5,0],[0,0,104]]}",
                     "{\"boxVectorsAngstrom\":[[-82,0,0],[0,80.5,0],[0,0,104]]}",
                     "{\"boxVectorsAngstrom\":[[82,0,0],[0,80.5,0],[0,0,1e400]]}",
                     "{\"boxVectorsAngstrom\":[[82,0,0],[0,80.5,0],[0,104]]}"
                 })
        {
            File.WriteAllText(fixture.Equilibrated.Molecule.TopologyPath!, invalidCell);
            await Refused(fixture.Equilibrated with { Molecule = fixture.Equilibrated.Molecule with
                { TopologySha256 = ConstructionFixture.Hash(fixture.Equilibrated.Molecule.TopologyPath!) } });
        }
    }

    private sealed class ExportFixture : IDisposable
    {
        private readonly ConstructionFixture _basis = new();
        public string Directory => _basis.Directory;
        public ApplicablePreparationPolicy Policy { get; private set; } = null!;
        public EquilibrationProtocol Protocol { get; private set; } = null!;
        public ConstructedExplicitSystem Constructed { get; private set; } = null!;
        public CompletedStage Minimized { get; private set; } = null!;
        public CompletedStage Equilibrated { get; private set; } = null!;
        public PreparationAssessmentResult MinimizedAssessment { get; private set; } = null!;
        public PreparationAssessmentResult EquilibratedAssessment { get; private set; } = null!;

        public static async Task<ExportFixture> CreateAsync()
        {
            var fixture = new ExportFixture();
            try { await fixture.InitializeAsync(); return fixture; }
            catch { fixture.Dispose(); throw; }
        }

        private async Task InitializeAsync()
        {
            var control = new EquilibrationStageControl("unrestrained", 5, 0.002, 303,
                1, "semiisotropic", 25, 0, 1, 0, 0, 1);
            Protocol = new EquilibrationProtocol("optional-protocol", 303, 47,
                ImmutableArray.Create(control), control with { Name = "extension" }, 0, 5,
                ImmutableArray.Create("potential"),
                ImmutableArray.Create(new EquilibrationObservable("potential", "kJ/mol", "system",
                    "system", "potentialEnergy", "none", "none", "none", null)),
                ImmutableArray.Create(new EquilibrationSufficiencyRule("potential", 1, 3, 1, 0.9)),
                "controlled declared comparison");
            // A policy that admits an optional stage must declare the same
            // stage-specific contact and intraprotein geometry checks as its
            // minimized predecessor before construction can bind the attempt.
            Policy = _basis.Policy with
            {
                OptionalEquilibration = Protocol,
                ContactCriteria = _basis.Policy.ContactCriteria.Add(
                    _basis.Policy.ContactCriteria.Single(item =>
                        item.StageKind == StageKind.Minimization) with
                    { StageKind = StageKind.Equilibration }),
                StageProteinGeometryCriteria = _basis.Policy.StageProteinGeometryCriteria.AddRange(
                    _basis.Policy.StageProteinGeometryCriteria.Select(item => item with
                    { StageKind = StageKind.Equilibration }))
            };
            var attempt = _basis.Attempt with
            { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(Policy) };
            var worker = new ConstructionWorker(_basis);
            var owner = new ConstructionOwner(worker, worker, worker);
            var started = await owner.StartAsync(attempt, _basis.Revision, _basis.Protein,
                _basis.Membrane, _basis.Placement, Policy, Directory, null, null,
                CancellationToken.None);
            Assert.True(started.State.Standing == StageExecutionStanding.Running,
                started.State.Message);
            Constructed = Assert.IsType<ConstructedExplicitSystem>(started.Constructed);
            Minimized = Stage("minimized-stage", StageKind.Minimization, null,
                Constructed.Molecule.TopologyPath!, Constructed.Molecule.TopologySha256!,
                StageTermination.Converged, null);
            var finalTopology = Write("equilibrated-topology.json",
                "{\"formatVersion\":1,\"boxVectorsAngstrom\":[[81.5,0,0],[0,80.5,0],[0,0,104]]}");
            Equilibrated = Stage("equilibrated-stage", StageKind.Equilibration, Minimized.Id,
                finalTopology, ConstructionFixture.Hash(finalTopology),
                StageTermination.Completed, EquilibrationObservationAdequacy.InsufficientAtBound);
            MinimizedAssessment = Assessment(Minimized, "minimized-assessment",
                PreparationCheckStanding.ChecksIncomplete);
            EquilibratedAssessment = Assessment(Equilibrated, "equilibrated-assessment",
                PreparationCheckStanding.ChecksIncomplete);
        }

        private CompletedStage Stage(string id, StageKind kind, string? sourceId,
            string topologyPath, string topologySha, StageTermination termination,
            EquilibrationObservationAdequacy? adequacy)
        {
            var coordinate = Write(id + ".cif", "controlled " + id + " coordinates");
            var system = Write(id + "-system.xml", "controlled " + id + " parameters");
            var state = Write(id + "-state.xml", "controlled " + id + " state");
            var molecule = Constructed.Molecule with
            {
                Id = id, CoordinatePath = coordinate,
                CoordinateSha256 = ConstructionFixture.Hash(coordinate),
                TopologyPath = topologyPath, TopologySha256 = topologySha,
                SystemXmlPath = system, SystemXmlSha256 = ConstructionFixture.Hash(system),
                StateXmlPath = state, StateXmlSha256 = ConstructionFixture.Hash(state),
                CellDescription = kind == StageKind.Equilibration ? null : "80 × 81 × 100 Å"
            };
            var windows = kind == StageKind.Equilibration
                ? ImmutableArray.Create(new EquilibrationWindowObservation("unrestrained", 5, 5,
                    303, 1, ImmutableArray.Create(new MeasuredValue("potential", -100,
                        "kJ/mol", "system")), ImmutableArray<string>.Empty))
                : ImmutableArray<EquilibrationWindowObservation>.Empty;
            var sample = new EquilibrationSample("unrestrained", 5,
                ImmutableArray.Create(new MeasuredValue("potential", -100, "kJ/mol", "system")));
            var observation = new StageObservation(id, Constructed.Attempt.Id, kind,
                ImmutableArray<MeasuredValue>.Empty, ImmutableArray<ScientificEvidence>.Empty,
                termination, "controlled OpenMM", DateTimeOffset.UtcNow, adequacy,
                kind == StageKind.Equilibration
                    ? ImmutableArray.Create(new EquilibrationObservationAssessment("potential", 5,
                        3, 0, 0, false)) : ImmutableArray<EquilibrationObservationAssessment>.Empty,
                kind == StageKind.Equilibration ? ImmutableArray.Create(sample) :
                    ImmutableArray<EquilibrationSample>.Empty,
                EquilibrationWindows: windows);
            return new CompletedStage(id, Constructed.Attempt, kind, molecule, observation,
                Constructed.Correspondence with { ResultId = id }, Policy.Id, sourceId,
                ImmutableArray<ScientificFinding>.Empty, DateTimeOffset.UtcNow);
        }

        private static PreparationAssessmentResult Assessment(CompletedStage stage, string id,
            PreparationCheckStanding checkStanding) => new(id, stage.Id, checkStanding,
            "Controlled stage-specific assessment", ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<ScientificFinding>.Empty, ImmutableArray<string>.Empty,
            DateTimeOffset.UtcNow, true);

        private string Write(string name, string value)
        {
            var path = Path.Combine(Directory, name);
            File.WriteAllText(path, value);
            return path;
        }

        public Task<BoundaryOutcome<CompletedStageBundle>> ExportAsync(ExportOwner owner,
            CompletedStage stage, PreparationAssessmentResult assessment,
            ApplicablePreparationPolicy policy)
        {
            string? optionalProtocolSha256 = null;
            if (stage.Kind == StageKind.Equilibration && policy.OptionalEquilibration is { } protocol &&
                EquilibrationProtocolFingerprint.TryCompute(protocol, out var digest))
                optionalProtocolSha256 = digest;
            return owner.ExportAsync(stage, assessment,
                _basis.Revision, _basis.Protein, _basis.Membrane, _basis.Placement,
                Constructed, Constructed.Derivation, policy,
                ImmutableArray<ResearcherDecision>.Empty, null, null, optionalProtocolSha256,
                Path.Combine(Directory, "exports", stage.Id + "-" + assessment.Id),
                CancellationToken.None);
        }

        public void Dispose() => _basis.Dispose();
    }

    private sealed class ControlledExportReadBack(int atomCount) : ICompletedStageExportWork
    {
        public List<ScientificWorkRequest<ExportVerificationPayload>> Requests { get; } = [];
        public Action<string>? AfterReadBack { get; set; }

        public Task<WorkerResult<ExportVerificationObservations>> VerifyExportAsync(
            ScientificWorkRequest<ExportVerificationPayload> request, CancellationToken cancellationToken)
        {
            Requests.Add(request);
            System.IO.Directory.CreateDirectory(request.WorkingDirectory);
            var output = Path.Combine(request.WorkingDirectory, "read-back.cif");
            File.Copy(request.Payload.TopologyCifPath, output, overwrite: true);
            var digest = ConstructionFixture.Hash(output);
            var after = AfterReadBack;
            AfterReadBack = null;
            after?.Invoke(request.WorkingDirectory);
            return Task.FromResult(new WorkerResult<ExportVerificationObservations>(
                request.RequestId, request.Payload.StudyRevisionId, request.Payload.AttemptId,
                request.Payload.StageId, WorkerResultStanding.Observed,
                ImmutableArray.Create(new WorkerArtifact("stageMmcif", output, digest)),
                new ExportVerificationObservations(atomCount, atomCount, atomCount,
                    true, true, true, true, 0, ImmutableArray<string>.Empty),
                new ProviderIdentity("controlled read-back", "1"), null, null));
        }
    }
}
