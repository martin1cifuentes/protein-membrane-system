using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ConstructionOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class ExplicitPreparationOwnerTests
{
    [Fact]
    public async Task One_native_invocation_supplies_the_checked_candidate_and_reviewable_proposal()
    {
        using var fixture = new ConstructionFixture();
        var worker = new ConstructionWorker(fixture);
        PreparationAttempt? accepted = null;
        var result = await Start(fixture, worker, onAccepted: value => accepted = value);

        Assert.Equal(fixture.Attempt.Id, accepted?.Id);
        Assert.True(result.State.Standing == StageExecutionStanding.Running, result.State.Message);
        Assert.Single(worker.ConstructionRequests);
        var request = worker.ConstructionRequests.Single();
        Assert.Equal(fixture.Attempt.Id, request.Payload.AttemptId);
        Assert.Equal(fixture.NativePatchPath, request.Payload.NativePatchPath);
        Assert.Equal(fixture.NativePatchSha, request.Payload.NativePatchSha256);
        Assert.Equal(fixture.Policy.Construction.ProviderVersion, request.Payload.ProviderVersion);
        Assert.Equal("DMPC", request.Payload.LipidTypeArgument);
        Assert.Equal("Na+", request.Payload.PositiveIonArgument);
        Assert.Equal("Cl-", request.Payload.NegativeIonArgument);
        Assert.Equal(fixture.Revision.Conditions.TargetNaClMolar, request.Payload.IonicStrengthMolar);
        Assert.Equal(fixture.Protein.Molecule.CoordinateSha256, request.Payload.PreparedPdbSha256);
        Assert.Equal(fixture.Placement.Proposal.OrientedProtein.CoordinateSha256,
            request.Payload.OrientedPdbSha256);
        Assert.Equal(fixture.Policy.ForceFieldFiles, request.Payload.ForceFieldFiles);
        Assert.Empty(request.Payload.ProviderAssets);
        Assert.Empty(request.Payload.SelectedSpeciesRepresentations);
        Assert.Equal("nativeOpenMm", JsonSerializer.SerializeToElement(request.Payload,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)).GetProperty("route").GetString());

        var proposed = Assert.IsType<ConstructionDerivation>(result.Derivation);
        var constructed = Assert.IsType<ConstructedExplicitSystem>(result.Constructed);
        Assert.Equal(fixture.Attempt.Id, proposed.AttemptId);
        Assert.Equal(fixture.Attempt.Id, constructed.Attempt.Id);
        Assert.Equal(new[] { (LeafletSide.Upper, "DMPC", 2), (LeafletSide.Lower, "DMPC", 3) },
            proposed.LipidCounts.Select(item => (item.PhysicalSide, item.SpeciesId, item.Count)));
        Assert.Equal(367, proposed.WaterCount);
        Assert.Equal(2, proposed.SodiumCount);
        Assert.Equal(1, proposed.ChlorideCount);
        Assert.Equal(-1, proposed.ProteinNetChargeElementary);
        Assert.Equal([80.0, 81.0, 100.0], proposed.CellAngstrom);
        Assert.Equal(proposed.CellAngstrom, constructed.ActualCellAngstrom);
        Assert.Equal(proposed.LipidCounts, constructed.AchievedComposition);
        Assert.Equal(0.0, constructed.MaximumProteinCoordinateDeviationAngstrom);
        Assert.Equal(55.4 / 370, proposed.EstimatedNaClMolar, 10);
        Assert.Contains(proposed.Approximations, text => text.Contains("same construction invocation"));
        Assert.Equal(constructed.Molecule.AtomCount, constructed.Correspondence.Atoms.Length);
        Assert.Equal(constructed.Molecule.Id, constructed.Correspondence.ResultId);
        Assert.Equal(3, constructed.Correspondence.Atoms.Count(item =>
            item.MoleculeRole == MoleculeRoleKind.Protein));
        Assert.Equal(ObservationStanding.Observed, constructed.LocalState?.Standing);
        Assert.Equal(1.1, constructed.LocalState!.Measurements.Single(item =>
            item.Name == "minimumIntermolecularDistanceAngstrom").Value);
        Assert.Equal(1.8, constructed.LocalState.Measurements.Single(item =>
            item.Name == "minimumIntermolecularHeavyAtomDistanceAngstrom").Value);
        Assert.Contains(constructed.Evidence, item => item.Bearing == EvidenceBearing.Supports &&
            item.Observation.Contains("mapped atoms", StringComparison.Ordinal));
        Assert.Contains("thermal equilibration has not been established", constructed.ConditionsTreatment);
    }

    [Fact]
    public async Task Native_output_changes_the_reviewable_counts_without_a_second_derivation_run()
    {
        using var fixture = new ConstructionFixture();
        var worker = new ConstructionWorker(fixture) { UpperLipidCount = 4, LowerLipidCount = 2,
            WaterCount = 368 };
        var result = await Start(fixture, worker);

        Assert.True(result.State.Standing == StageExecutionStanding.Running, result.State.Message);
        Assert.Single(worker.ConstructionRequests);
        Assert.Equal([4, 2], result.Derivation!.LipidCounts.Select(item => item.Count));
        Assert.Equal(368, result.Derivation.WaterCount);
        Assert.Equal(result.Derivation.LipidCounts, result.Constructed!.AchievedComposition);
    }

    [Fact]
    public async Task Identical_coordinate_bytes_in_separate_attempts_have_distinct_constructed_subjects()
    {
        using var fixture = new ConstructionFixture();
        var worker = new ConstructionWorker(fixture);
        var first = Assert.IsType<ConstructedExplicitSystem>((await Start(fixture, worker)).Constructed);
        var laterAttempt = fixture.Attempt with { Id = "later-attempt" };
        var second = Assert.IsType<ConstructedExplicitSystem>((await new ConstructionOwner(worker, worker, worker)
            .StartAsync(laterAttempt, fixture.Revision, fixture.Protein, fixture.Membrane,
                fixture.Placement, fixture.Policy, fixture.Directory, null, null,
                TestContext.Current.CancellationToken)).Constructed);

        Assert.Equal(first.Molecule.Id, second.Molecule.Id);
        Assert.Equal("constructed-" + fixture.Attempt.Id, first.Id);
        Assert.Equal("constructed-" + laterAttempt.Id, second.Id);
        Assert.All(first.Evidence, evidence => Assert.Equal(first.Id, evidence.SubjectId));
        Assert.All(second.Evidence, evidence => Assert.Equal(second.Id, evidence.SubjectId));
    }

    [Theory]
    [InlineData("resourceRefused")]
    [InlineData("resourceLimit")]
    public async Task Typed_memgen_capacity_refusal_keeps_the_attempt_and_unobserved_trial(
        string failureCode)
    {
        using var fixture = new ConstructionFixture();
        var (policy, attempt) = MemgenRoute(fixture);
        var worker = new RefusingMemgenWorker { FailureCode = failureCode };
        var progress = new CapturingStageProgress();
        var unused = new ConstructionWorker(fixture);
        var result = await new ConstructionOwner(worker, unused, unused).StartAsync(attempt,
            fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
            fixture.Directory, null, progress, TestContext.Current.CancellationToken);

        Assert.Equal(attempt.Id, result.Attempt?.Id);
        Assert.Equal(StageExecutionStanding.ResourceRefused, result.State.Standing);
        Assert.Equal(failureCode, result.State.FailureCode);
        Assert.Null(result.Constructed);
        var request = Assert.Single(worker.Requests);
        AssertUnobservedTrial(result, progress, request, ConstructionTrialStanding.Failed, failureCode);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public async Task Memgen_stop_records_the_active_unobserved_trial_whether_worker_throws_or_returns(
        bool workerThrows)
    {
        using var fixture = new ConstructionFixture();
        var (policy, attempt) = MemgenRoute(fixture);
        using var stop = new CancellationTokenSource();
        var worker = new RefusingMemgenWorker
        {
            BeforeResult = stop.Cancel,
            ThrowWhenCancelled = workerThrows
        };
        var progress = new CapturingStageProgress();
        var unused = new ConstructionWorker(fixture);
        var result = await new ConstructionOwner(worker, unused, unused).StartAsync(attempt,
            fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
            fixture.Directory, null, progress, stop.Token);

        Assert.Equal(attempt.Id, result.Attempt?.Id);
        Assert.Equal(StageExecutionStanding.Stopped, result.State.Standing);
        Assert.Equal("stopped", result.State.FailureCode);
        Assert.Null(result.Constructed);
        var request = Assert.Single(worker.Requests);
        AssertUnobservedTrial(result, progress, request, ConstructionTrialStanding.Stopped, "stopped");
    }

    [Fact]
    public async Task Memgen_stop_on_retry_keeps_the_prior_attributable_trial_and_current_request()
    {
        using var fixture = new ConstructionFixture();
        var (policy, attempt) = MemgenRoute(fixture);
        using var stop = new CancellationTokenSource();
        var worker = new RefusingMemgenWorker
        {
            RetryOnceWithAttributableTrial = true,
            BeforeResult = stop.Cancel,
            ThrowWhenCancelled = true
        };
        var progress = new CapturingStageProgress();
        var unused = new ConstructionWorker(fixture);
        var result = await new ConstructionOwner(worker, unused, unused).StartAsync(attempt,
            fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
            fixture.Directory, null, progress, stop.Token);

        Assert.Equal(StageExecutionStanding.Stopped, result.State.Standing);
        Assert.Equal(2, worker.Requests.Count);
        Assert.Equal(2, result.Trials.Length);
        Assert.Equal(worker.Requests[0].Payload.TrialId, result.Trials[0].TrialId);
        Assert.Equal("zeroRoundedSpecies", result.Trials[0].FailureCode);
        Assert.Equal(ConstructionTrialStanding.Failed, result.Trials[0].Standing);
        Assert.Equal(worker.Requests[1].Payload.TrialId, result.Trials[1].TrialId);
        Assert.Equal(1, result.Trials[1].TrialIndex);
        Assert.Equal(30, result.Trials[1].LateralPaddingAngstrom);
        Assert.Equal(17.5, result.Trials[1].AqueousPaddingAngstrom);
        Assert.Equal(ConstructionTrialStanding.Stopped, result.Trials[1].Standing);
        Assert.Equal("stopped", result.Trials[1].FailureCode);
        Assert.Empty(result.Trials[1].ProposedLipidCounts);
        Assert.Empty(result.Trials[1].ActualCellAngstrom);
        Assert.Equal(result.Trials, progress.States.Where(state => state.Trial is not null)
            .Select(state => state.Trial!).ToImmutableArray());
    }

    [Theory]
    [InlineData("zeroRoundedSpecies", 30.0, 60.0, 17.5, 17.5)]
    [InlineData("chargeBudgetRefusal", 15.0, 15.0, 35.0, 70.0)]
    public async Task Memgen_attributable_failures_retry_only_within_the_declared_cell_bound(
        string failureCode, double secondLateral, double thirdLateral,
        double secondAqueous, double thirdAqueous)
    {
        using var fixture = new ConstructionFixture();
        var (policy, attempt) = MemgenRoute(fixture);
        var worker = new RefusingMemgenWorker
        {
            AttributableFailureCodes = [failureCode, failureCode, failureCode]
        };
        var progress = new CapturingStageProgress();
        var unused = new ConstructionWorker(fixture);
        var result = await new ConstructionOwner(worker, unused, unused).StartAsync(attempt,
            fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
            fixture.Directory, null, progress, TestContext.Current.CancellationToken);

        Assert.Equal(attempt.Id, result.Attempt?.Id);
        Assert.Null(result.Constructed);
        Assert.Null(result.Derivation);
        var aqueousLimitRefusal = failureCode == "chargeBudgetRefusal";
        var dispatched = aqueousLimitRefusal ? 2 : 3;
        Assert.Equal(aqueousLimitRefusal ? StageExecutionStanding.ResourceRefused :
            StageExecutionStanding.Failed, result.State.Standing);
        Assert.Equal(aqueousLimitRefusal ? "resourceRefused" : failureCode,
            result.State.FailureCode);
        if (aqueousLimitRefusal)
            Assert.Contains("186 Å", result.State.Message);
        Assert.Equal(3000, policy.Construction.MaximumConstructionSeconds);
        Assert.Equal(dispatched, worker.Requests.Count);
        Assert.Equal(dispatched, worker.Tokens.Count);
        Assert.All(worker.Tokens, token => Assert.Equal(worker.Tokens[0], token));
        Assert.Equal(new[] { 15.0, secondLateral, thirdLateral }.Take(dispatched),
            worker.Requests.Select(request => request.Payload.LateralPaddingAngstrom!.Value));
        Assert.Equal(new[] { 17.5, secondAqueous, thirdAqueous }.Take(dispatched),
            worker.Requests.Select(request => request.Payload.AqueousPaddingAngstrom!.Value));
        Assert.Equal(dispatched, worker.Requests.Select(request => request.Payload.TrialId).Distinct().Count());
        Assert.Equal(new[] { 0, 1, 2 }.Take(dispatched),
            worker.Requests.Select(request => request.Payload.TrialIndex!.Value));
        Assert.All(worker.Requests, request =>
        {
            Assert.Equal(attempt.Id, request.Payload.AttemptId);
            Assert.Equal(fixture.Revision.Id, request.Payload.StudyRevisionId);
            Assert.Equal(policy.Construction.ProviderVersion, request.Payload.ProviderVersion);
        });
        Assert.Equal(dispatched, result.Trials.Length);
        for (var index = 0; index < result.Trials.Length; index++)
        {
            var trial = result.Trials[index];
            Assert.Equal(worker.Requests[index].Payload.TrialId, trial.TrialId);
            Assert.Equal(index, trial.TrialIndex);
            Assert.Equal(ConstructionTrialStanding.Failed, trial.Standing);
            Assert.Equal(failureCode, trial.FailureCode);
            Assert.Empty(trial.ProposedLipidCounts);
            Assert.Empty(trial.AchievedLipidCounts);
            Assert.Empty(trial.ActualCellAngstrom);
        }
        Assert.Equal(result.Trials, progress.States.Where(state => state.Trial is not null)
            .Select(state => state.Trial!).ToImmutableArray());
    }

    [Fact]
    public async Task Attributable_Memgen_failure_retains_only_verified_trial_diagnostics()
    {
        using var fixture = new ConstructionFixture();
        var (policy, attempt) = MemgenRoute(fixture);
        var outside = Path.Combine(fixture.Directory, "outside-provider-options.json");
        File.WriteAllText(outside, "outside the trial");
        var worker = new RefusingMemgenWorker
        {
            AttributableFailureCodes = ["zeroRoundedSpecies", "zeroRoundedSpecies", "zeroRoundedSpecies"],
            FailureArtifacts = request =>
            {
                var inside = Path.Combine(request.WorkingDirectory, "memgen-trial-0");
                System.IO.Directory.CreateDirectory(inside);
                var log = Path.Combine(inside, "packmol.log");
                var changed = Path.Combine(inside, "provider-stdout.log");
                var unknown = Path.Combine(inside, "unlisted.txt");
                File.WriteAllText(log, "the observed Packmol refusal");
                File.WriteAllText(changed, "different bytes");
                File.WriteAllText(unknown, "not an allowed diagnostic role");
                return ImmutableArray.Create(
                    new WorkerArtifact("providerPackmolLog", log, ConstructionFixture.Hash(log)),
                    new WorkerArtifact("providerOptions", outside, ConstructionFixture.Hash(outside)),
                    new WorkerArtifact("providerStdout", changed, new string('0', 64)),
                    new WorkerArtifact("arbitraryOutput", unknown, ConstructionFixture.Hash(unknown)));
            }
        };
        var progress = new CapturingStageProgress();
        var unused = new ConstructionWorker(fixture);
        var result = await new ConstructionOwner(worker, unused, unused).StartAsync(attempt,
            fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
            fixture.Directory, null, progress, TestContext.Current.CancellationToken);

        Assert.Equal(3, result.Trials.Length);
        foreach (var trial in result.Trials)
        {
            var diagnostic = Assert.Single(trial.DiagnosticArtifacts);
            Assert.Equal("providerPackmolLog", diagnostic.Role);
            Assert.Equal(PreparationPhase.ProviderPacking, diagnostic.Phase);
            Assert.Equal("packmol.log", diagnostic.FileName);
            Assert.Equal(ConstructionFixture.Hash(diagnostic.LocalPath!), diagnostic.Sha256);
            Assert.Null(diagnostic.DownloadUrl);
            Assert.Null(diagnostic.StructureUrl);
            Assert.DoesNotContain("localPath", JsonSerializer.Serialize(diagnostic,
                new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        }
        Assert.Equal(result.Trials, progress.States.Where(state => state.Trial is not null)
            .Select(state => state.Trial!).ToImmutableArray());
    }

    [Fact]
    public async Task Checked_Memgen_intermediates_require_current_digest_and_nonlinked_trial_path()
    {
        using var fixture = new ConstructionFixture();
        var trial = Path.Combine(fixture.Directory, "trial");
        System.IO.Directory.CreateDirectory(trial);
        var packed = Path.Combine(trial, "system.pdb");
        var restart = Path.Combine(trial, "system_min.restrt");
        File.WriteAllText(packed, "packed coordinates");
        File.WriteAllText(restart, "final unrestrained restart");
        var linked = Path.Combine(trial, "linked.pdb");
        File.CreateSymbolicLink(linked, packed);
        var reported = ImmutableArray.Create(
            new WorkerArtifact("providerPacked", packed, ConstructionFixture.Hash(packed)),
            new WorkerArtifact("amberFinalRestart", restart, ConstructionFixture.Hash(restart)),
            new WorkerArtifact("providerProtein", linked, ConstructionFixture.Hash(packed)));
        var method = typeof(ConstructionOwner).GetMethod("VerifiedMemgenDiagnosticsAsync",
            System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.NonPublic)!;
        async Task<ImmutableArray<TrialDiagnosticArtifact>> Verify(
            ImmutableArray<WorkerArtifact> artifacts) =>
            await Assert.IsAssignableFrom<Task<ImmutableArray<TrialDiagnosticArtifact>>>(
                method.Invoke(null, [artifacts, trial, TestContext.Current.CancellationToken]));

        var accepted = await Verify(reported);
        Assert.Equal(["providerPacked", "amberFinalRestart"], accepted.Select(item => item.Role));
        Assert.Equal(PreparationPhase.ProviderPacking, accepted[0].Phase);
        Assert.Equal(PreparationPhase.ProviderConditioningUnrestrained, accepted[1].Phase);
        Assert.All(accepted, item => Assert.StartsWith(trial, item.LocalPath,
            StringComparison.Ordinal));
        File.AppendAllText(restart, "tampered after report");
        var changed = await Verify(reported);
        Assert.Equal("providerPacked", Assert.Single(changed).Role);
    }

    [Fact]
    public void Unequal_Memgen_aqueous_regions_at_the_strict_charge_boundary_require_neutralization_only()
    {
        const double lowerVolume = 15000;
        const double upperVolume = 26000;
        const double volume = lowerVolume + upperVolume;
        const double molarFactor = 6.02214076e-4;
        var regions = ImmutableArray.Create(
            new AqueousRegionAccount(LeafletSide.Lower, lowerVolume, 1, 0, 1),
            new AqueousRegionAccount(LeafletSide.Upper, upperVolume, 2, 0, 1));
        // |Q|/2 equals the smaller nominal salt count, so Memgen adds only
        // neutralizers even though the two aqueous volumes differ.
        var account = new ConstructionConditionAccount(
            0, 1933, 1933, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 2, 2, 2, 0, 0, volume,
            0, 2 / (molarFactor * volume), 0, 55.4 * 2 / 1933,
            "neutralizationOnly", regions);
        var generated = ImmutableArray.Create(
            new ObservedSpeciesCount(GeneratedComponentRoleKind.Water, LeafletSide.Upper,
                "HOH", 1933),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.NegativeIon, LeafletSide.Upper,
                "CL", 1),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.NegativeIon, LeafletSide.Lower,
                "CL", 1));
        var observed = new ConstructionObservations(0, generated, [80, 80, 100], [30, 30, 60],
            1933, 0, 2, 0, 2, -100, 0, 0, true, null,
            ImmutableArray<string>.Empty, ImmutableArray<string>.Empty,
            ImmutableArray<string>.Empty, null!);
        var validator = typeof(ConstructionOwner).GetMethod("MemgenConditionsCoherent",
            System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.NonPublic)!;
        bool Coherent(ConstructionConditionAccount value, ConstructionObservations actual,
            ImmutableArray<ObservedSpeciesCount> species) =>
            Assert.IsType<bool>(validator.Invoke(null, [value, actual, species, 0.15, 55.4]));

        Assert.True(Coherent(account, observed, generated));
        Assert.False(Coherent(account with { SaltBranch = "chargeCompensated" }, observed, generated));
        Assert.False(Coherent(account with { AqueousRegions = regions.SetItem(0,
            regions[0] with { FlooredNominalSaltCount = 2 }) }, observed, generated));
        // One salt threshold lies inside the printed 11070.26 Å³ interval:
        // Memgen's hidden unrounded volume can floor to one even though
        // recomputing from the two-decimal log would floor to zero.
        const double loggedBoundaryVolume = 11070.26;
        var boundaryRegions = regions.SetItem(0, regions[0] with
        {
            EstimatedVolumeAngstromCubed = loggedBoundaryVolume,
            FlooredNominalSaltCount = 1
        });
        var boundaryVolume = loggedBoundaryVolume + upperVolume;
        var boundary = account with
        {
            AqueousRegions = boundaryRegions,
            EstimatedAqueousVolumeAngstromCubed = boundaryVolume,
            ChlorideAqueousMolar = 2 / (molarFactor * boundaryVolume)
        };
        Assert.True(Coherent(boundary, observed, generated));
        Assert.False(Coherent(boundary with { AqueousRegions = boundaryRegions.SetItem(0,
            boundaryRegions[0] with { FlooredNominalSaltCount = 2 }) }, observed, generated));
        var belowBoundaryRegions = boundaryRegions.SetItem(0, boundaryRegions[0] with
        { EstimatedVolumeAngstromCubed = 11070.25 });
        var belowBoundaryVolume = 11070.25 + upperVolume;
        Assert.False(Coherent(boundary with
        {
            AqueousRegions = belowBoundaryRegions,
            EstimatedAqueousVolumeAngstromCubed = belowBoundaryVolume,
            ChlorideAqueousMolar = 2 / (molarFactor * belowBoundaryVolume)
        }, observed, generated));
        // LEaP's measured charge can differ slightly from the integral charge
        // Memgen receives; that noise cannot turn strict equality into a pass.
        var noisy = account with { PreparedFormalChargeElementary = 1.9999999,
            ChargePdbDeltaElementary = -0.0000001 };
        Assert.True(Coherent(noisy, observed with { ProteinNetChargeElementary = 1.9999999 },
            generated));
        Assert.False(Coherent(noisy with { SaltBranch = "chargeCompensated" },
            observed with { ProteinNetChargeElementary = 1.9999999 }, generated));
        var negativeRegions = regions.Select(region => region with
        { ProviderGeneratedSodiumCount = 1, ProviderGeneratedChlorideCount = 0 }).ToImmutableArray();
        var negative = account with { PreparedFormalChargeElementary = -2,
            ProviderResidueNameChargeElementary = -2,
            ProviderGeneratedSodiumCount = 2, ProviderGeneratedChlorideCount = 0,
            FinalSodiumCount = 2, FinalChlorideCount = 0,
            SodiumAqueousMolar = 2 / (molarFactor * volume), ChlorideAqueousMolar = 0,
            SodiumFiniteWaterMolar = 55.4 * 2 / 1933, ChlorideFiniteWaterMolar = 0,
            AqueousRegions = negativeRegions };
        var negativeSpecies = ImmutableArray.Create(
            generated[0],
            new ObservedSpeciesCount(GeneratedComponentRoleKind.PositiveIon, LeafletSide.Upper,
                "NA", 1),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.PositiveIon, LeafletSide.Lower,
                "NA", 1));
        var negativeObserved = observed with { ProteinNetChargeElementary = -2,
            PositiveIonCount = 2, NegativeIonCount = 0 };
        Assert.True(Coherent(negative, negativeObserved, negativeSpecies));
        Assert.False(Coherent(negative with { SaltBranch = "chargeCompensated" },
            negativeObserved, negativeSpecies));
    }

    private static (ApplicablePreparationPolicy Policy, PreparationAttempt Attempt) MemgenRoute(
        ConstructionFixture fixture)
    {
        var assets = ImmutableArray.Create(new ProviderAsset("controlled-provider", "2026.3.25",
            fixture.NativePatchPath, fixture.NativePatchSha));
        var settings = new MemgenConstructionSettings("sander", "ff19SB", "lipid21", "tip3p",
            15, 17.5, 23, 20, 100, 20, 2, 250, 250, 10, 2,
            true, true, true, true, true, true, true, true, true, true, true);
        var construction = fixture.Policy.Construction with
        {
            Route = ConstructionRouteKind.PackmolMemgen,
            SaltConvention = SaltConventionKind.MemgenChargeCompensated,
            ProviderName = "PACKMOL-Memgen",
            ProviderVersion = "2026.3.25",
            NativePatchPath = null,
            NativePatchSha256 = null,
            LipidTypeArgument = null,
            CoveredSpeciesIds = ImmutableArray.Create("POPC", "POPE", "DLPC", "DLPE", "DMPC",
                "DOPC", "DPPC", "CHL1", "HOH", "NA", "CL"),
            MaximumAtomCount = 120000,
            MaximumCellDimensionAngstrom = 180,
            MaximumConstructionSeconds = 3000,
            ProviderAssets = assets,
            Memgen = settings
        };
        var policy = fixture.Policy with { Construction = construction };
        var attempt = fixture.Attempt with
        {
            Route = ConstructionRouteKind.PackmolMemgen,
            SaltConvention = SaltConventionKind.MemgenChargeCompensated,
            NativePatchSha256 = null,
            ConstructionProviderVersion = construction.ProviderVersion,
            ProviderAssets = assets,
            PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy)
        };
        return (policy, attempt);
    }

    [Fact]
    public async Task Changed_selected_lipid_bytes_refuse_native_and_Memgen_before_an_attempt()
    {
        using var fixture = new ConstructionFixture();
        var selected = fixture.Membrane.SpeciesRepresentations[0];
        var membrane = fixture.Membrane with
        {
            SpeciesRepresentations = fixture.Membrane.SpeciesRepresentations.SetItem(0,
                selected with { CoordinateTemplateSha256 = new string('0', 64) })
        };
        var (memgenPolicy, memgenAttempt) = MemgenRoute(fixture);
        foreach (var (policy, attempt) in new[]
                 { (fixture.Policy, fixture.Attempt), (memgenPolicy, memgenAttempt) })
        {
            var worker = new ConstructionWorker(fixture);
            var accepted = false;
            var result = await new ConstructionOwner(worker, worker, worker).StartAsync(
                attempt, fixture.Revision, fixture.Protein, membrane, fixture.Placement,
                policy, fixture.Directory, _ => accepted = true, null,
                TestContext.Current.CancellationToken);
            Assert.Null(result.Attempt);
            Assert.False(accepted);
            Assert.Empty(worker.ConstructionRequests);
            Assert.Equal(StageExecutionStanding.ResourceRefused, result.State.Standing);
            Assert.Equal("A selected molecular representation changed before preparation admission.",
                result.State.Message);
        }
    }

    private static void AssertUnobservedTrial(PreparationStartResult result,
        CapturingStageProgress progress, ScientificWorkRequest<ConstructionPayload> request,
        ConstructionTrialStanding standing, string failureCode)
    {
        var trial = Assert.Single(result.Trials);
        Assert.Equal(request.Payload.TrialId, trial.TrialId);
        Assert.Equal(request.Payload.TrialIndex, trial.TrialIndex);
        Assert.Equal(15, trial.LateralPaddingAngstrom);
        Assert.Equal(17.5, trial.AqueousPaddingAngstrom);
        Assert.Equal(standing, trial.Standing);
        Assert.Equal(failureCode, trial.FailureCode);
        Assert.Empty(trial.ProposedLipidCounts);
        Assert.Empty(trial.AchievedLipidCounts);
        Assert.Empty(trial.CleanupRemovedLipidCounts);
        Assert.Empty(trial.ProposedCellAngstrom);
        Assert.Empty(trial.ActualCellAngstrom);
        Assert.Null(trial.Conditions);
        Assert.True(trial.DiagnosticArtifacts.IsDefaultOrEmpty);
        Assert.Equal(trial, Assert.Single(progress.States, state => state.Trial is not null).Trial);
    }

    [Fact]
    public async Task Missing_or_incompatible_prestart_basis_refuses_without_invoking_the_worker()
    {
        using var fixture = new ConstructionFixture();
        var cases = new (string Name, bool MissingRepresentation,
            Func<ConstructionFixture, (PreparationAttempt Attempt,
            StudyRevision Revision, AssessedPreparedProtein Protein, AssessedMembraneModel Membrane,
            AssessedProteinMembranePlacement Placement, ApplicablePreparationPolicy Policy)> Change)[]
        {
            ("revision", false, f => (f.Attempt, f.Revision with { Id = "other" }, f.Protein,
                f.Membrane, f.Placement, f.Policy)),
            ("adopted placement", false, f => (f.Attempt,
                f.Revision with { AdoptedPlacementProposalId = null }, f.Protein,
                f.Membrane, f.Placement, f.Policy)),
            ("placement standing", false, f => (f.Attempt, f.Revision, f.Protein, f.Membrane,
                f.Placement with { Standing = AssessmentStanding.NotEstablished }, f.Policy)),
            ("protein correspondence", true, f => (f.Attempt, f.Revision,
                f.Protein with { Correspondence = f.Protein.Correspondence with { Complete = false } },
                f.Membrane, f.Placement, f.Policy)),
            ("policy fingerprint", false, f => (f.Attempt with
                { PolicyFingerprintSha256 = new string('0', 64) }, f.Revision, f.Protein,
                f.Membrane, f.Placement, f.Policy)),
            ("provider version binding", false, f => (f.Attempt with
                { ConstructionProviderVersion = "other" }, f.Revision, f.Protein,
                f.Membrane, f.Placement, f.Policy)),
            ("patch digest binding", false, f => (f.Attempt with
                { NativePatchSha256 = new string('0', 64) }, f.Revision, f.Protein,
                f.Membrane, f.Placement, f.Policy)),
            ("chemical state", false, f => (f.Attempt, f.Revision,
                f.Protein with { ChemicalStatePolicyVersion = "other" }, f.Membrane,
                f.Placement, f.Policy)),
            ("physical leaflet", false, f => (f.Attempt, f.Revision, f.Protein,
                f.Membrane with { Intended = f.Membrane.Intended with { Upper =
                    f.Membrane.Intended.Upper with { Fractions = ImmutableArray.Create(
                        new LipidFraction("POPC", 1)) } } }, f.Placement, f.Policy))
        };
        foreach (var (name, missingRepresentation, change) in cases)
        {
            var (attempt, revision, protein, membrane, placement, policy) = change(fixture);
            var worker = new ConstructionWorker(fixture);
            var result = await new ConstructionOwner(worker, worker, worker).StartAsync(attempt,
                revision, protein, membrane, placement, policy, fixture.Directory, null, null,
                TestContext.Current.CancellationToken);
            Assert.Null(result.Attempt);
            Assert.Null(result.Constructed);
            Assert.Empty(worker.ConstructionRequests);
            Assert.Equal(string.Empty, result.State.AttemptId);
            Assert.True(result.State.Standing == (missingRepresentation
                ? StageExecutionStanding.ResourceRefused : StageExecutionStanding.Failed),
                name + ": " + result.State.Message);
            if (!missingRepresentation)
                Assert.Equal("preAdmissionRefused", result.State.FailureCode);
            Assert.Equal(missingRepresentation
                    ? "A qualified molecular representation, native package patch, or exact parameter asset is unavailable before start."
                    : "Corresponding supported inputs and an identified complete construction policy are required.",
                result.State.Message);
        }

        using var mixed = new ConstructionFixture(mixed: true);
        var mixedWorker = new ConstructionWorker(mixed);
        var mixedResult = await Start(mixed, mixedWorker);
        Assert.Equal(StageExecutionStanding.Failed, mixedResult.State.Standing);
        Assert.Equal("preAdmissionRefused", mixedResult.State.FailureCode);
        Assert.Equal("The selected native recipe requires its exact qualified pure lipid on both leaflets.",
            mixedResult.State.Message);
        Assert.Empty(mixedWorker.ConstructionRequests);
    }

    [Fact]
    public async Task Incomplete_native_policy_or_missing_exact_assets_refuse_before_the_worker()
    {
        using var fixture = new ConstructionFixture();
        var cases = new (string Name, bool MissingAsset,
            Func<ApplicablePreparationPolicy, ApplicablePreparationPolicy> Change)[]
        {
            ("missing evidence", false, p => p with { EvidenceReferences = ImmutableArray<string>.Empty }),
            ("wrong provider", false, p => p with { Construction = p.Construction with
                { ProviderName = "unqualified provider" } }),
            ("wrong lipid argument", false, p => p with { Construction = p.Construction with
                { LipidTypeArgument = "POPC" } }),
            ("unbounded time", false, p => p with { Construction = p.Construction with
                { MaximumConstructionSeconds = 0 } }),
            ("missing padding", false, p => p with { Construction = p.Construction with
                { MinimumPaddingNanometers = 0 } }),
            ("wrong water convention", false, p => p with { Construction = p.Construction with
                { WaterMolarityForIonRounding = 55.5 } }),
            ("missing package patch", true, p => p with { Construction = p.Construction with
                { NativePatchPath = p.Construction.NativePatchPath + ".absent" } }),
            ("changed package patch", true, p => p with { Construction = p.Construction with
                { NativePatchSha256 = new string('0', 64) } }),
            ("changed force field", true, p => p with { ForceFieldFiles = p.ForceFieldFiles.SetItem(0,
                p.ForceFieldFiles[0] with { Sha256 = new string('0', 64) }) }),
            ("missing catalogue reference", true, p => p with { Water = p.Water with
                { CoordinateTemplatePath = p.Water.CoordinateTemplatePath + ".absent" } }),
            ("missing heavy contact criterion", false, p => p with { ConstructionCriteria =
                p.ConstructionCriteria.RemoveAt(0) }),
            ("wrong final force target", false, p => p with
                { FinalUnrestrainedRmsForceTargetKjMolNm = 11 })
        };
        foreach (var (name, missingAsset, change) in cases)
        {
            var policy = change(fixture.Policy);
            var attempt = fixture.Attempt with
            { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy),
                ForceFieldFiles = policy.ForceFieldFiles,
                ConstructionProviderVersion = policy.Construction.ProviderVersion,
                NativePatchSha256 = policy.Construction.NativePatchSha256 };
            var worker = new ConstructionWorker(fixture);
            var result = await new ConstructionOwner(worker, worker, worker).StartAsync(attempt,
                fixture.Revision, fixture.Protein, fixture.Membrane, fixture.Placement, policy,
                fixture.Directory, null, null, TestContext.Current.CancellationToken);
            Assert.Null(result.Attempt);
            Assert.Null(result.Constructed);
            Assert.Empty(worker.ConstructionRequests);
            Assert.Equal(missingAsset ? StageExecutionStanding.ResourceRefused :
                StageExecutionStanding.Failed, result.State.Standing);
            if (!missingAsset)
                Assert.Equal("preAdmissionRefused", result.State.FailureCode);
            Assert.Equal(name == "missing catalogue reference"
                    ? "A selected molecular representation changed before preparation admission."
                    : missingAsset
                    ? "A qualified molecular representation, native package patch, or exact parameter asset is unavailable before start."
                    : name == "wrong lipid argument"
                        ? "The selected native recipe requires its exact qualified pure lipid on both leaflets."
                        : "Corresponding supported inputs and an identified complete construction policy are required.",
                result.State.Message);
        }
    }

    [Fact]
    public async Task Native_counts_ions_cell_and_bounded_resource_must_be_coherent()
    {
        using var fixture = new ConstructionFixture();
        var cases = new (string Name, Func<ConstructionObservations, ConstructionObservations> Change)[]
        {
            ("missing upper lipid", o => o with { SpeciesCounts = o.SpeciesCounts.RemoveAt(0) }),
            ("wrong leaflet species", o => o with { SpeciesCounts = o.SpeciesCounts.SetItem(0,
                o.SpeciesCounts[0] with { SpeciesId = "POPC" }) }),
            ("zero lower lipid", o => o with { SpeciesCounts = o.SpeciesCounts.SetItem(1,
                o.SpeciesCounts[1] with { Count = 0 }) }),
            ("missing water", o => o with { WaterCount = 0 }),
            ("wrong salt pair count", o => o with { PositiveIonCount = 3 }),
            ("wrong neutralizer sign", o => o with { PositiveIonCount = 1, NegativeIonCount = 2 }),
            ("non-neutral system", o => o with { NetChargeElementary = 1 }),
            ("nonintegral protein charge", o => o with { ProteinNetChargeElementary = -1.2 }),
            ("nonfinite cell", o => o with { ActualCellAngstrom = [80.0, double.NaN, 100.0] }),
            ("oversized cell", o => o with { ActualCellAngstrom = [101.0, 81.0, 100.0] }),
            ("cutoff exceeds half cell", o => o with { ActualCellAngstrom = [20.0, 81.0, 100.0] }),
            ("protein periodic gap", o => o with { ProteinPeriodicImageGapsAngstrom = [1.0, 40.0, 40.0] }),
            ("atom bound", o => o with { AtomCount = 100001 })
        };
        foreach (var (name, change) in cases)
        {
            var worker = new ConstructionWorker(fixture) { ChangeConstruction = result => result with
                { Observations = change(result.Observations!) } };
            var result = await Start(fixture, worker);
            Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
            Assert.Null(result.Constructed);
            Assert.True(result.Derivation is null, name + ": unexpected derivation");
            Assert.Single(worker.ConstructionRequests);
            Assert.Equal("The same-invocation native counts and cell were not coherently observed.",
                result.State.Message);
        }
    }

    [Fact]
    public async Task Native_artifacts_identity_parameters_and_local_contacts_are_required_for_promotion()
    {
        using var fixture = new ConstructionFixture();
        var cases = new (string Name, Func<WorkerResult<ConstructionObservations>,
            WorkerResult<ConstructionObservations>> Change)[]
        {
            ("wrong provider", r => r with { Provider = new ProviderIdentity("other", "1") }),
            ("wrong patch bytes", r => r with { Observations = r.Observations! with
                { NativePatchSha256 = new string('0', 64) } }),
            ("unreported solvent count", r => r with { Observations = r.Observations! with
                { SpeciesCounts = r.Observations.SpeciesCounts.RemoveAt(2) } }),
            ("protein identity", r => r with { Observations = r.Observations! with
                { ProteinIdentityAndBondsPreserved = false } }),
            ("protein motion", r => r with { Observations = r.Observations! with
                { MaximumProteinCoordinateDeviationAngstrom = 0.1 } }),
            ("unparameterized state", r => r with { Observations = r.Observations! with
                { ParameterWarnings = ["missing lipid parameters"] } }),
            ("nonfinite energy", r => r with { Observations = r.Observations! with
                { InitialPotentialEnergyKjMol = double.NaN } }),
            ("heavy overlap", r => r with { Observations = r.Observations! with
                { LocalState = ChangeMetric(r.Observations!.LocalState,
                    "minimumIntermolecularHeavyAtomDistanceAngstrom", 1.4) } }),
            ("inverted leaflets", r => r with { Observations = r.Observations! with
                { LocalState = ChangeMetric(r.Observations!.LocalState,
                    "leafletHeadSeparationAngstrom", -20) } }),
            ("large protein offset", r => r with { Observations = r.Observations! with
                { LocalState = ChangeMetric(r.Observations!.LocalState,
                    "proteinBilayerMidplaneOffsetAngstrom", 30) } }),
            ("no protein lipid contact", r => r with { Observations = r.Observations! with
                { LocalState = r.Observations.LocalState with { RolePairMeasurements =
                    [new LocalRolePairMeasurement(MoleculeRoleKind.Protein, MoleculeRoleKind.Lipid,
                        0, null)] } } }),
            ("missing role coverage", r => r with { Observations = r.Observations! with
                { LocalState = r.Observations.LocalState with
                    { CoveredRolePairs = ImmutableArray<LocalContactRolePair>.Empty } } }),
            ("unavailable observation", r => r with { Observations = r.Observations! with
                { LocalState = r.Observations.LocalState with
                    { Standing = ObservationStanding.Unavailable } } }),
            ("missing artifact", r => r with { Artifacts = r.Artifacts.RemoveAt(0) }),
            ("incorrect artifact hash", r => r with { Artifacts = r.Artifacts.SetItem(0,
                r.Artifacts[0] with { Sha256 = new string('0', 64) }) }),
            ("incorrect mapped atom count", r => r with { Observations = r.Observations! with
                { CorrespondedResultAtomCount = r.Observations.AtomCount - 1 } })
        };
        foreach (var (name, change) in cases)
        {
            var worker = new ConstructionWorker(fixture) { ChangeConstruction = change };
            var result = await Start(fixture, worker);
            Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
            Assert.NotNull(result.Derivation);
            Assert.Null(result.Constructed);
            Assert.Single(worker.ConstructionRequests);
            Assert.Equal("Actual membership, periodic geometry, contacts, parameters or correspondence were not established.",
                result.State.Message);
        }
    }

    [Fact]
    public async Task Incorrect_correspondence_cannot_promote_an_otherwise_complete_candidate()
    {
        using var fixture = new ConstructionFixture();
        var worker = new ConstructionWorker(fixture) { ChangeConstruction = result =>
        {
            var artifact = result.Artifacts.Single(item => item.Role == "correspondenceJson");
            var mapping = JsonSerializer.Deserialize<SourceToResultCorrespondence>(
                File.ReadAllText(artifact.Path), new JsonSerializerOptions(JsonSerializerDefaults.Web))!;
            var wrong = mapping with { Atoms = mapping.Atoms.SetItem(0,
                mapping.Atoms[0] with { SourceAtomId = "different-source-atom" }) };
            File.WriteAllText(artifact.Path, JsonSerializer.Serialize(wrong,
                new JsonSerializerOptions(JsonSerializerDefaults.Web)));
            return result with { Artifacts = result.Artifacts.Select(item =>
                item.Role == "correspondenceJson" ? item with
                { Sha256 = ConstructionFixture.Hash(item.Path) } : item).ToImmutableArray() };
        } };
        var result = await Start(fixture, worker);
        Assert.NotNull(result.Derivation);
        Assert.Null(result.Constructed);
        Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
    }

    [Fact]
    public async Task Failure_unobserved_and_wrong_correlation_never_establish_a_candidate()
    {
        using var fixture = new ConstructionFixture();
        var cases = new (string Name, StageExecutionStanding Standing, string Reason,
            Func<WorkerResult<ConstructionObservations>,
            WorkerResult<ConstructionObservations>> Change)[]
        {
            ("worker failed", StageExecutionStanding.Failed, "Native construction failed.",
                r => r with { Standing = WorkerResultStanding.Failed,
                FailureMessage = "Native construction failed." }),
            ("worker unobserved", StageExecutionStanding.Unobserved,
                "The native candidate was not observed.",
                r => r with { Standing = WorkerResultStanding.Unobserved,
                Observations = null }),
            ("request ID", StageExecutionStanding.Failed,
                "The native candidate was not observed.", r => r with { RequestId = "other" }),
            ("attempt ID", StageExecutionStanding.Failed,
                "The native candidate was not observed.", r => r with { AttemptId = "other" }),
            ("revision ID", StageExecutionStanding.Failed,
                "The native candidate was not observed.", r => r with { StudyRevisionId = "other" })
        };
        foreach (var (name, standing, reason, change) in cases)
        {
            var worker = new ConstructionWorker(fixture) { ChangeConstruction = change };
            var result = await Start(fixture, worker);
            Assert.Null(result.Constructed);
            Assert.Null(result.Derivation);
            Assert.Single(worker.ConstructionRequests);
            Assert.Equal(standing, result.State.Standing);
            Assert.Equal(reason, result.State.Message);
        }
    }

    private static Task<PreparationStartResult> Start(ConstructionFixture fixture,
        ConstructionWorker worker, Action<PreparationAttempt>? onAccepted = null) =>
        new ConstructionOwner(worker, worker, worker).StartAsync(fixture.Attempt, fixture.Revision,
            fixture.Protein, fixture.Membrane, fixture.Placement, fixture.Policy, fixture.Directory,
            onAccepted, null, TestContext.Current.CancellationToken);

    private static LocalStateObservations ChangeMetric(LocalStateObservations local,
        string name, double value) => local with { Measurements = local.Measurements.Select(item =>
            item.Name == name ? item with { Value = value } : item).ToImmutableArray() };
}

internal sealed class ConstructionFixture : IDisposable
{
    public string Directory { get; } = Path.Combine(Path.GetTempPath(),
        "pims-construction-" + Guid.NewGuid().ToString("N"));
    public StudyRevision Revision { get; }
    public AssessedPreparedProtein Protein { get; }
    public AssessedMembraneModel Membrane { get; }
    public AssessedProteinMembranePlacement Placement { get; }
    public ApplicablePreparationPolicy Policy { get; }
    public PreparationAttempt Attempt { get; }
    public string NativePatchPath { get; }
    public string NativePatchSha { get; }

    public ConstructionFixture(bool mixed = false, bool retainedPartner = false)
    {
        System.IO.Directory.CreateDirectory(Directory);
        string Write(string name, string content)
        { var path = Path.Combine(Directory, name); File.WriteAllText(path, content); return path; }
        var preparedPath = Write("prepared.pdb", "three identified source atoms");
        var orientedPath = Write("oriented.pdb", "three positioned source atoms");
        var graphPath = Write("graph.json", "{}");
        var sourcePath = Write("source.pdb", "three source atoms");
        NativePatchPath = Write("DMPC.pdb", "controlled installed OpenMM DMPC patch bytes");
        NativePatchSha = Hash(NativePatchPath);
        var proteinFf = Write("protein.xml", "protein parameters");
        var lipidFf = Write("lipid.xml", "lipid parameters");
        var waterFf = Write("water.xml", "water and ion parameters");
        var dm = Write("dmpc.cif", "DMPC catalogue reference");
        var popc = Write("popc.cif", "POPC catalogue reference");
        var cholesterol = Write("chl1.cif", "CHL1 catalogue reference");
        var water = Write("water.cif", "water catalogue reference");
        var sodium = Write("sodium.cif", "sodium catalogue reference");
        var chloride = Write("chloride.cif", "chloride catalogue reference");
        var forceFields = ImmutableArray.Create(
            new ForceFieldAsset("protein", "1", "ff19SB", proteinFf, Hash(proteinFf)),
            new ForceFieldAsset("lipid", "1", "Lipid21", lipidFf, Hash(lipidFf)),
            new ForceFieldAsset("water-ion", "1", "TIP3P", waterFf, Hash(waterFf)));
        var dmRepresentation = Representation("DMPC", "lipid", lipidFf, dm, 3, 0, 60);
        var popcRepresentation = Representation("POPC", "lipid", lipidFf, popc, 3, 0, 80);
        var cholesterolRepresentation = Representation("CHL1", "sterol", lipidFf,
            cholesterol, 3, 0, 0);
        var waterRepresentation = Representation("HOH", "water", waterFf, water, 3, 0, 0);
        var sodiumRepresentation = Representation("NA", "ion", waterFf, sodium, 1, 1, 0);
        var chlorideRepresentation = Representation("CL", "ion", waterFf, chloride, 1, -1, 0);
        var conditions = FixedStudyConditions.Initial;
        var source = new StructuralSource("source", SourceRouteKind.Rcsb, "controlled source", sourcePath,
            Hash(sourcePath), "6QWR", "first model");
        var intended = new IntendedProteinModel("intended", source, 0, null,
            ImmutableArray.Create(new ChainSelection("A", "A")),
            retainedPartner ? ImmutableArray.Create(new PartnerSelection("partner", true,
                "controlled retained partner")) : ImmutableArray<PartnerSelection>.Empty,
            ImmutableArray<AlternateLocationChoice>.Empty);
        var addresses = Enumerable.Range(1, 3)
            .Select(number => new ResidueAddress(0, "A", number, "", "A")).ToArray();
        var preparedSha = Hash(preparedPath);
        var preparedAtoms = addresses.Select((address, index) => new AtomCorrespondence(index,
                $"{index}:A:{address.Residue}::CA", $"source-{index}", AtomOriginKind.Source,
                MoleculeRoleKind.Protein, AtomRoleKind.Backbone, "C", address, null)).ToImmutableArray();
        if (retainedPartner)
            preparedAtoms = preparedAtoms.SetItem(2, preparedAtoms[2] with
            {
                SourceAtomId = "partner-atom-2",
                SourceResidue = new ResidueAddress(0, "P", 1, "", "P"),
                MoleculeRole = MoleculeRoleKind.RetainedPartner,
                AtomRole = AtomRoleKind.Body
            });
        var correspondence = new SourceToResultCorrespondence(source.Sha256, preparedSha,
            preparedAtoms, true);
        Protein = new AssessedPreparedProtein("protein", "revision", intended,
            new MolecularArtifact("prepared", preparedPath, preparedSha, graphPath, null, null, 3, null,
                Hash(graphPath)), "chemical", ImmutableArray<ResidueVariantChoice>.Empty,
            ImmutableArray<PreparationChangeProposal>.Empty, correspondence,
            ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty,
            ImmutableArray<string>.Empty, ChemicalStatePolicyVersion: "1",
            StructuralAssessmentPolicyId: "structural-policy",
            StructuralAssessmentPolicyVersion: "1");
        var upper = new LeafletComposition(LeafletSide.Upper, mixed
            ? ImmutableArray.Create(new LipidFraction("DMPC", 0.6),
                new LipidFraction("POPC", 0.4), new LipidFraction("CHL1", 0))
            : ImmutableArray.Create(new LipidFraction("DMPC", 1)));
        var lower = new LeafletComposition(LeafletSide.Lower, mixed
            ? ImmutableArray.Create(new LipidFraction("DMPC", 0.25),
                new LipidFraction("POPC", 0.75), new LipidFraction("CHL1", 0))
            : ImmutableArray.Create(new LipidFraction("DMPC", 1)));
        var model = new MembraneModel("membrane", upper, lower, conditions, "controlled bilayer");
        Membrane = new AssessedMembraneModel("assessed-membrane", "revision", model,
            mixed ? ImmutableArray.Create(dmRepresentation, popcRepresentation, cholesterolRepresentation)
                : ImmutableArray.Create(dmRepresentation), ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<string>.Empty, PolicyId: "membrane-policy", PolicyVersion: "1");
        var proposal = new PlacementProposal("proposal", Protein.Id, model.Id,
            ProteinTopologyKind.MembraneSpanning,
            new MolecularArtifact("oriented", orientedPath, Hash(orientedPath), graphPath,
                null, null, 3, null, Hash(graphPath)), 0, 20, 8,
            PlacementPhysicalSide.Both, "upper extracellular",
            ImmutableArray<string>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<string>.Empty);
        Placement = new AssessedProteinMembranePlacement("placement", "revision", proposal,
            AssessmentStanding.Supported, "controlled independent witness", ImmutableArray<ScientificFinding>.Empty,
            DateTimeOffset.UtcNow);
        Revision = new StudyRevision("revision", 1, intended, model, proposal.Id, conditions);
        var scope = new PreparationPolicyScope(source.Sha256, 0, null, intended.Chains,
            retainedPartner ? ImmutableArray.Create("partner") : ImmutableArray<string>.Empty,
            Hash(graphPath), Protein.ChemicalStatePolicyId,
            Protein.ChemicalStatePolicyVersion!,
            PreparationPolicyFingerprint.ComputeResidueVariants(Protein.ResidueVariants), upper, lower,
            conditions, ProteinTopologyKind.MembraneSpanning,
            Protein.StructuralAssessmentPolicyId!, Protein.StructuralAssessmentPolicyVersion!,
            Membrane.PolicyId, Membrane.PolicyVersion);
        var radii = ImmutableDictionary<string, double>.Empty.Add("C", 1.7).Add("O", 1.52)
            .Add("H", 1.2).Add("Na", 2.27).Add("Cl", 1.75);
        var rolePair = new LocalContactRolePair(MoleculeRoleKind.Protein, MoleculeRoleKind.Lipid);
        var metrics = ImmutableArray.Create("minimumIntermolecularDistanceAngstrom",
            "minimumIntermolecularHeavyAtomDistanceAngstrom", "leafletHeadSeparationAngstrom",
            "proteinBilayerMidplaneOffsetAngstrom");
        var geometryKinds = ImmutableArray.Create("covalentBond", "chainContinuity", "nonbondedDistance");
        Policy = new ApplicablePreparationPolicy("policy", "1", "canonical-amino-acid-assembly",
            ImmutableArray.Create("controlled identified construction evidence"), forceFields,
            20000, 10, 0.01, 0.01, 0.01,
            new MolecularDynamicsSystemSettings("PME", 1, "HBonds", true, 0.0005, null, true, true, null),
            new ConstructionPolicy("construction", "1", ImmutableArray.Create("controlled native basis"),
                ImmutableArray.Create(ProteinTopologyKind.MembraneSpanning),
                mixed ? ImmutableArray.Create("DMPC", "POPC", "CHL1", "HOH", "NA", "CL")
                    : ImmutableArray.Create("DMPC", "HOH", "NA", "CL"),
                "OpenMM Modeller.addMembrane", "8.6.0.dev-c6173db", NativePatchPath, NativePatchSha,
                "DMPC", "Na+", "Cl-", 1, 55.4, 100000, 100, 3600,
                "Water-equivalent salt estimate for a native finite cell",
                ProviderAssets: ImmutableArray<ProviderAsset>.Empty),
            waterRepresentation, sodiumRepresentation, chlorideRepresentation,
            new LocalStateObservationSpec(ImmutableArray.Create(rolePair), radii, 6, 100, true, metrics),
            ImmutableArray.Create(
                new LocalStateCriterion("minimumIntermolecularHeavyAtomDistanceAngstrom",
                    "angstrom", "wholeSystem", 1.5, 6),
                new LocalStateCriterion("leafletHeadSeparationAngstrom", "angstrom", "bilayer", 10, 40),
                new LocalStateCriterion("proteinBilayerMidplaneOffsetAngstrom", "angstrom",
                    "proteinVsBilayer", -20, 20)),
            ImmutableArray.Create(
                new LocalRolePairCriterion(null, MoleculeRoleKind.Protein, MoleculeRoleKind.Lipid, 1, 1, 6),
                new LocalRolePairCriterion(StageKind.Minimization, MoleculeRoleKind.Protein,
                    MoleculeRoleKind.Lipid, 1, 1, 6)),
            ImmutableArray<PreparationAssessmentCriterion>.Empty,
            new ProteinGeometryMeasurementSpec(geometryKinds, radii, 5, 3, 100),
            geometryKinds.Select(kind => new StageProteinGeometryCriterion(StageKind.Minimization,
                new ProteinGeometryCriterion(kind, kind == "nonbondedDistance" ? 0.75 : 1.0,
                    kind == "nonbondedDistance" ? null : 2.5, false))).ToImmutableArray(),
            null, ImmutableArray.Create("No positive minimized-stage assessment basis"), scope);
        Attempt = new PreparationAttempt("attempt", Revision.Id, Protein.Id, Membrane.Id,
            Placement.Id, Policy.Id, DateTimeOffset.UtcNow, Policy.Version,
            PreparationPolicyFingerprint.Compute(Policy), Policy.ForceFieldFiles,
            Policy.Construction.ProviderVersion, NativePatchSha,
            ProviderAssets: ImmutableArray<ProviderAsset>.Empty);

        MolecularRepresentation Representation(string species, string category, string ff, string template,
            int atomCount, double charge, double area) =>
            new(species, species + "-chemistry", category, ff, Hash(ff), template, Hash(template),
                atomCount, charge, area, 0, category == "lipid" ? ImmutableArray.Create(0) :
                ImmutableArray<int>.Empty, category == "lipid" ? "Lipid21" : "TIP3P", "8.6",
                ImmutableArray<string>.Empty, ImmutableArray<MolecularStereoCheck>.Empty);
    }

    public static string Hash(string path) =>
        Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

    public void Dispose() => System.IO.Directory.Delete(Directory, true);
}

internal sealed class ConstructionWorker(ConstructionFixture fixture) :
    IExplicitConstructionWork, IMinimizationWork, IOptionalEquilibrationWork
{
    public List<ScientificWorkRequest<ConstructionPayload>> ConstructionRequests { get; } = [];
    public Func<WorkerResult<ConstructionObservations>, WorkerResult<ConstructionObservations>>?
        ChangeConstruction { get; init; }
    public int UpperLipidCount { get; init; } = 2;
    public int LowerLipidCount { get; init; } = 3;
    public int WaterCount { get; init; } = 367;
    public int SodiumCount { get; init; } = 2;
    public int ChlorideCount { get; init; } = 1;
    public double ProteinCharge { get; init; } = -1;
    public ImmutableArray<double> CellAngstrom { get; init; } = [80.0, 81.0, 100.0];

    public Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
        ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken)
    {
        ConstructionRequests.Add(request);
        cancellationToken.ThrowIfCancellationRequested();
        var speciesCounts = ImmutableArray.Create(
            new ObservedSpeciesCount(GeneratedComponentRoleKind.Lipid, LeafletSide.Upper,
                "DMPC", UpperLipidCount),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.Lipid, LeafletSide.Lower,
                "DMPC", LowerLipidCount),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.Water, LeafletSide.Upper,
                "HOH", WaterCount),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.PositiveIon, LeafletSide.Upper,
                "NA", SodiumCount),
            new ObservedSpeciesCount(GeneratedComponentRoleKind.NegativeIon, LeafletSide.Upper,
                "CL", ChlorideCount));
        var representations = new Dictionary<GeneratedComponentRoleKind, MolecularRepresentation>
        {
            [GeneratedComponentRoleKind.Lipid] = request.Payload.Lipid!,
            [GeneratedComponentRoleKind.Water] = request.Payload.Water,
            [GeneratedComponentRoleKind.PositiveIon] = request.Payload.Sodium,
            [GeneratedComponentRoleKind.NegativeIon] = request.Payload.Chloride
        };
        var atomCount = fixture.Protein.Molecule.AtomCount + speciesCounts.Sum(item =>
            item.Count * representations[item.Role].AtomCount);
        var topologyPath = Path.Combine(request.WorkingDirectory, "constructed.cif");
        var topologyJsonPath = Path.Combine(request.WorkingDirectory, "topology.json");
        var systemPath = Path.Combine(request.WorkingDirectory, "system.xml");
        var statePath = Path.Combine(request.WorkingDirectory, "state.xml");
        var mappingPath = Path.Combine(request.WorkingDirectory, "constructed-correspondence.json");
        File.WriteAllText(topologyPath, "controlled full system coordinates");
        File.WriteAllText(topologyJsonPath, "controlled topology and bonds");
        File.WriteAllText(systemPath, "controlled complete-system parameters");
        File.WriteAllText(statePath, "controlled native candidate state");
        var atoms = fixture.Protein.Correspondence.Atoms.Select(item => item with
            { ResultAtomId = "protein:" + item.ResultAtomId }).ToList();
        var index = atoms.Count;
        foreach (var species in speciesCounts)
            for (var molecule = 0; molecule < species.Count; molecule++)
                for (var atom = 0; atom < representations[species.Role].AtomCount; atom++)
                {
                    atoms.Add(new AtomCorrespondence(index, $"component:{index}", null,
                        AtomOriginKind.Generated,
                        species.Role is GeneratedComponentRoleKind.PositiveIon or
                            GeneratedComponentRoleKind.NegativeIon ? MoleculeRoleKind.Ion :
                            species.Role == GeneratedComponentRoleKind.Lipid ? MoleculeRoleKind.Lipid :
                            MoleculeRoleKind.Water,
                        species.Role == GeneratedComponentRoleKind.Lipid && atom == 0 ? AtomRoleKind.Head :
                            AtomRoleKind.Body,
                        species.Role == GeneratedComponentRoleKind.PositiveIon ? "Na" :
                            species.Role == GeneratedComponentRoleKind.NegativeIon ? "Cl" :
                            species.Role == GeneratedComponentRoleKind.Water ? "O" : "C",
                        null, null, species.PhysicalSide, species.SpeciesId, species.Role));
                    index++;
                }
        var mapping = new SourceToResultCorrespondence(
            fixture.Placement.Proposal.OrientedProtein.CoordinateSha256,
            ConstructionFixture.Hash(topologyPath), atoms.ToImmutableArray(), true);
        File.WriteAllText(mappingPath, JsonSerializer.Serialize(mapping,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)));
        var pair = new LocalContactRolePair(MoleculeRoleKind.Protein, MoleculeRoleKind.Lipid);
        var local = new LocalStateObservations(ObservationStanding.Observed, null,
            ImmutableArray.Create(
                new MeasuredValue("minimumIntermolecularDistanceAngstrom", 1.1,
                    "angstrom", "wholeSystem"),
                new MeasuredValue("minimumIntermolecularHeavyAtomDistanceAngstrom", 1.8,
                    "angstrom", "wholeSystem"),
                new MeasuredValue("leafletHeadSeparationAngstrom", 20, "angstrom", "bilayer"),
                new MeasuredValue("proteinBilayerMidplaneOffsetAngstrom", 0,
                    "angstrom", "proteinVsBilayer")),
            ImmutableArray<LocalContactObservation>.Empty, ImmutableArray.Create(pair),
            ImmutableArray<string>.Empty,
            ImmutableArray.Create(new LocalRolePairMeasurement(MoleculeRoleKind.Protein,
                MoleculeRoleKind.Lipid, 1, 2)));
        var observed = new ConstructionObservations(atomCount, speciesCounts, CellAngstrom,
            [40.0, 41.0, 60.0], WaterCount, SodiumCount, ChlorideCount, 0, ProteinCharge,
            -100, atomCount, 0, true, request.Payload.NativePatchSha256,
            ImmutableArray<string>.Empty, ImmutableArray<string>.Empty,
            ImmutableArray<string>.Empty, local, NativePatchMode: "installed");
        var result = new WorkerResult<ConstructionObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, null,
            WorkerResultStanding.Observed,
            ImmutableArray.Create(Artifact(topologyPath, "topologyCif"),
                Artifact(topologyJsonPath, "topologyJson"), Artifact(systemPath, "systemXml"),
                Artifact(statePath, "stateXml"), Artifact(mappingPath, "correspondenceJson")),
            observed, new ProviderIdentity(request.Payload.ProviderName,
                request.Payload.ProviderVersion), null, null);
        return Task.FromResult(ChangeConstruction?.Invoke(result) ?? result);
    }

    private static WorkerArtifact Artifact(string path, string role) =>
        new(role, path, ConstructionFixture.Hash(path));

    public Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
        ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Minimization is exercised by its owning test.");

    public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Stage observation is exercised by its owning test.");

    public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
        ScientificWorkRequest<EquilibrationPayload> request,
        IProgress<EquilibrationWorkProgress>? progress, CancellationToken cancellationToken) =>
        throw new InvalidOperationException("Optional equilibration is outside Slice 4.");
}

internal sealed class RefusingMemgenWorker : IExplicitConstructionWork
{
    public List<ScientificWorkRequest<ConstructionPayload>> Requests { get; } = [];
    public List<CancellationToken> Tokens { get; } = [];
    public string FailureCode { get; init; } = "resourceRefused";
    public Action? BeforeResult { get; init; }
    public bool ThrowWhenCancelled { get; init; }
    public bool RetryOnceWithAttributableTrial { get; init; }
    public ImmutableArray<string> AttributableFailureCodes { get; init; } = [];
    public Func<ScientificWorkRequest<ConstructionPayload>, ImmutableArray<WorkerArtifact>>?
        FailureArtifacts { get; init; }

    public Task<WorkerResult<ConstructionObservations>> ConstructSystemAsync(
        ScientificWorkRequest<ConstructionPayload> request, CancellationToken cancellationToken)
    {
        Requests.Add(request);
        Tokens.Add(cancellationToken);
        if (!AttributableFailureCodes.IsDefaultOrEmpty)
        {
            var code = AttributableFailureCodes[Requests.Count - 1];
            var trial = new ConstructionTrialSummary(request.Payload.TrialId!, request.Payload.TrialIndex!.Value,
                ConstructionTrialStanding.Failed, request.Payload.LateralPaddingAngstrom!.Value,
                request.Payload.AqueousPaddingAngstrom!.Value,
                ImmutableArray<SpeciesCount>.Empty, ImmutableArray<SpeciesCount>.Empty,
                ImmutableArray<SpeciesCount>.Empty, ImmutableArray<double>.Empty,
                ImmutableArray<double>.Empty, null, code,
                "The complete provider reported its typed population or charge refusal.",
                ImmutableArray<TrialDiagnosticArtifact>.Empty);
            return Task.FromResult(new WorkerResult<ConstructionObservations>(request.RequestId,
                request.Payload.StudyRevisionId, request.Payload.AttemptId, null,
                WorkerResultStanding.Failed,
                FailureArtifacts?.Invoke(request) ?? ImmutableArray<WorkerArtifact>.Empty, null,
                new ProviderIdentity("PACKMOL-Memgen", "2026.3.25"), code,
                trial.Message, JsonSerializer.SerializeToElement(new { trial },
                    new JsonSerializerOptions(JsonSerializerDefaults.Web))));
        }
        if (RetryOnceWithAttributableTrial && Requests.Count == 1)
        {
            var trial = new ConstructionTrialSummary(request.Payload.TrialId!, request.Payload.TrialIndex!.Value,
                ConstructionTrialStanding.Failed, request.Payload.LateralPaddingAngstrom!.Value,
                request.Payload.AqueousPaddingAngstrom!.Value,
                ImmutableArray<SpeciesCount>.Empty, ImmutableArray<SpeciesCount>.Empty,
                ImmutableArray<SpeciesCount>.Empty, ImmutableArray<double>.Empty,
                ImmutableArray<double>.Empty, null, "zeroRoundedSpecies",
                "The complete provider trial rounded a species to zero.",
                ImmutableArray<TrialDiagnosticArtifact>.Empty);
            return Task.FromResult(new WorkerResult<ConstructionObservations>(request.RequestId,
                request.Payload.StudyRevisionId, request.Payload.AttemptId, null,
                WorkerResultStanding.Failed, ImmutableArray<WorkerArtifact>.Empty, null,
                new ProviderIdentity("PACKMOL-Memgen", "2026.3.25"), "zeroRoundedSpecies",
                trial.Message, JsonSerializer.SerializeToElement(new { trial },
                    new JsonSerializerOptions(JsonSerializerDefaults.Web))));
        }
        BeforeResult?.Invoke();
        if (ThrowWhenCancelled)
            cancellationToken.ThrowIfCancellationRequested();
        return Task.FromResult(new WorkerResult<ConstructionObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, null,
            WorkerResultStanding.Failed, ImmutableArray<WorkerArtifact>.Empty, null,
            new ProviderIdentity("PACKMOL-Memgen", "2026.3.25"), FailureCode,
            "Actual provider result exceeds the authorized atom budget."));
    }
}

internal sealed class CapturingStageProgress : IProgress<StageExecutionState>
{
    public List<StageExecutionState> States { get; } = [];

    public void Report(StageExecutionState value) => States.Add(value);
}
