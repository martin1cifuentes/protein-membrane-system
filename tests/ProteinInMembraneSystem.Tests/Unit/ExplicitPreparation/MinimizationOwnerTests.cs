using System.Collections.Immutable;
using System.Security.Cryptography;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using MinimizationOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.Minimization;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class MinimizationOwnerTests
{
    [Fact]
    public async Task Observed_final_unrestrained_convergence_establishes_an_exact_completed_stage_without_qualifying_it()
    {
        using var fixture = new MinimizationFixture();
        var worker = new MinimizationWorker(fixture);
        var owner = new MinimizationOwner(worker);

        var result = await owner.RunAsync(fixture.Constructed, fixture.Policy, fixture.Directory,
            null, TestContext.Current.CancellationToken);

        Assert.Equal(StageExecutionStanding.Completed, result.State.Standing);
        var stage = Assert.IsType<CompletedStage>(result.CompletedStage);
        Assert.Equal(StageKind.Minimization, stage.Kind);
        Assert.Equal(fixture.Attempt.Id, stage.Attempt.Id);
        Assert.Equal(fixture.Constructed.Molecule.AtomCount, stage.Molecule.AtomCount);
        Assert.Equal(stage.Id, stage.Observation.StageId);
        Assert.Equal(stage.Id, stage.Correspondence.ResultId);
        Assert.Equal(StageTermination.Converged, stage.Observation.Termination);
        Assert.Contains(stage.Observation.Measurements, item =>
            item.Name == "finalRmsForce" && item.Value == 8.0);
        Assert.Contains(stage.Observation.Measurements, item =>
            item.Name == "finalRawRmsForce" && item.Value == 500.0);
        Assert.Contains(stage.Observation.Evidence, item =>
            item.Observation.Contains("constraint-tangent-per-particle-v1", StringComparison.Ordinal));
        Assert.All(stage.Observation.Evidence, item => Assert.Equal(EvidenceBearing.Context, item.Bearing));
        Assert.Single(worker.MinimizationRequests);
        Assert.Single(worker.ObservationRequests);
        Assert.Equal(fixture.Constructed.Molecule.CoordinateSha256,
            worker.MinimizationRequests[0].Payload.TopologyCifSha256);
        Assert.Equal(worker.MinimizationRequests[0].Payload.StageId,
            worker.ObservationRequests[0].Payload.StageId);
        Assert.Equal(fixture.MinimizedCoordinateSha,
            worker.ObservationRequests[0].Payload.TopologyCifSha256);
    }

    [Fact]
    public async Task A_provider_reply_does_not_complete_without_correlated_unrestrained_convergence()
    {
        using var fixture = new MinimizationFixture();
        var variants = new (string Name, Func<WorkerResult<MinimizationObservations>, WorkerResult<MinimizationObservations>> Change)[]
        {
            ("request", reply => reply with { RequestId = "other" }),
            ("revision", reply => reply with { StudyRevisionId = "other" }),
            ("attempt", reply => reply with { AttemptId = "other" }),
            ("stage", reply => reply with { StageId = "other" }),
            ("unobserved", reply => reply with { Standing = WorkerResultStanding.Unobserved }),
            ("failed", reply => reply with { Standing = WorkerResultStanding.Failed }),
            ("missing observation", reply => reply with { Observations = null }),
            ("restrained", reply => reply with { Observations = reply.Observations! with
                { FinalTreatmentUnrestrained = false } }),
            ("iteration limit", reply => reply with { Observations = reply.Observations! with
                { Termination = StageTermination.MaxIterations } }),
            ("high force", reply => reply with { Observations = reply.Observations! with
                { FinalRmsForceKjMolNm = 11.0 } }),
            ("missing projection method", reply => reply with { Observations = reply.Observations! with
                { FinalRmsForceMethod = null } }),
            ("wrong projection method", reply => reply with { Observations = reply.Observations! with
                { FinalRmsForceMethod = "raw-state-components" } }),
            ("missing raw diagnostic", reply => reply with { Observations = reply.Observations! with
                { FinalRawRmsForceKjMolNm = null } }),
            ("nonfinite raw diagnostic", reply => reply with { Observations = reply.Observations! with
                { FinalRawRmsForceKjMolNm = double.NaN } }),
            ("unobserved constraints", reply => reply with { Observations = reply.Observations! with
                { MaximumRelativeConstraintError = null } }),
            ("constraint tolerance exceeded", reply => reply with { Observations = reply.Observations! with
                { MaximumRelativeConstraintError = 0.00002 } }),
            ("unbounded constraint tolerance", reply => reply with { Observations = reply.Observations! with
                { AppliedConstraintTolerance = 0.1 } }),
            ("unknown termination", reply => reply with { Observations = reply.Observations! with
                { Termination = StageTermination.Unknown } }),
            ("lower energy without force convergence", reply => reply with
                { Observations = reply.Observations! with { FinalPotentialEnergyKjMol = -10000,
                    FinalRmsForceKjMolNm = 11.0 } }),
            ("wrong atoms", reply => reply with { Observations = reply.Observations! with
                { FinalAtomCount = 2 } }),
            ("nonfinite energy", reply => reply with { Observations = reply.Observations! with
                { FinalPotentialEnergyKjMol = double.NaN } }),
            ("nonfinite force", reply => reply with { Observations = reply.Observations! with
                { FinalRmsForceKjMolNm = double.NaN } }),
            ("numerical warning", reply => reply with { Observations = reply.Observations! with
                { NumericalWarnings = ImmutableArray.Create("diverged") } }),
            ("missing minimized coordinates", reply => reply with { Artifacts = reply.Artifacts
                .Where(item => item.Role != "minimizedCif").ToImmutableArray() }),
            ("missing minimized coordinate bytes", reply => reply with { Artifacts = reply.Artifacts
                .Select(item => item.Role == "minimizedCif" ? item with
                    { Path = item.Path + ".absent" } : item).ToImmutableArray() }),
            ("changed minimized coordinate digest", reply => reply with { Artifacts = reply.Artifacts
                .Select(item => item.Role == "minimizedCif" ? item with
                    { Sha256 = new string('0', 64) } : item).ToImmutableArray() }),
            ("changed minimized state digest", reply => reply with { Artifacts = reply.Artifacts
                .Select(item => item.Role == "minimizedStateXml" ? item with
                    { Sha256 = new string('0', 64) } : item).ToImmutableArray() })
        };
        foreach (var (name, change) in variants)
        {
            var worker = new MinimizationWorker(fixture) { ChangeMinimization = change };
            var result = await new MinimizationOwner(worker).RunAsync(fixture.Constructed,
                fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.True(result.CompletedStage is null, name + ": completed stage was promoted");
            Assert.True(result.State.Standing != StageExecutionStanding.Completed,
                name + ": operation was marked completed");
            Assert.True(worker.ObservationRequests.Count == 0,
                name + ": stage remeasurement was invoked");
        }
    }

    [Fact]
    public async Task Exact_force_target_is_inclusive_and_energy_decrease_alone_is_insufficient()
    {
        using var fixture = new MinimizationFixture();
        var boundary = new MinimizationWorker(fixture) { ChangeMinimization = reply => reply with
        { Observations = reply.Observations! with
            { FinalRmsForceKjMolNm = 10, FinalPotentialEnergyKjMol = -10000 } } };
        var result = await new MinimizationOwner(boundary).RunAsync(fixture.Constructed,
            fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.NotNull(result.CompletedStage);
        Assert.Equal(StageExecutionStanding.Completed, result.State.Standing);
    }

    [Fact]
    public async Task Completed_numerics_without_corresponding_stage_observation_remain_unobserved()
    {
        using var fixture = new MinimizationFixture();
        var variants = new (string Name, Func<WorkerResult<StageObservationObservations>, WorkerResult<StageObservationObservations>> Change)[]
        {
            ("request", reply => reply with { RequestId = "other" }),
            ("revision", reply => reply with { StudyRevisionId = "other" }),
            ("attempt", reply => reply with { AttemptId = "other" }),
            ("stage", reply => reply with { StageId = "other" }),
            ("unobserved", reply => reply with { Standing = WorkerResultStanding.Unobserved }),
            ("failed", reply => reply with { Standing = WorkerResultStanding.Failed }),
            ("atom order", reply => reply with { Observations = reply.Observations! with
                { AtomOrderMatched = false } }),
            ("bonds", reply => reply with { Observations = reply.Observations! with
                { BondsMatched = false } }),
            ("count", reply => reply with { Observations = reply.Observations! with
                { AtomCount = 2 } }),
            ("numerical warning", reply => reply with { Observations = reply.Observations! with
                { NumericalWarnings = ImmutableArray.Create("nonfinite") } })
        };
        foreach (var (name, change) in variants)
        {
            var worker = new MinimizationWorker(fixture) { ChangeObservation = change };
            var result = await new MinimizationOwner(worker).RunAsync(fixture.Constructed,
                fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.True(result.CompletedStage is null, name + ": completed stage was promoted");
            Assert.Equal(StageExecutionStanding.Unobserved, result.State.Standing);
            Assert.Single(worker.ObservationRequests);
        }
    }

    [Fact]
    public async Task Stopped_stage_remeasurement_does_not_report_an_unobserved_completion()
    {
        using var fixture = new MinimizationFixture();
        var worker = new MinimizationWorker(fixture) { ChangeObservation = reply => reply with
        { Standing = WorkerResultStanding.Stopped, Observations = null } };
        var result = await new MinimizationOwner(worker).RunAsync(fixture.Constructed,
            fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Null(result.CompletedStage);
        Assert.Equal(StageExecutionStanding.Stopped, result.State.Standing);
    }

    [Fact]
    public async Task A_stop_or_wrong_policy_withholds_the_completed_stage()
    {
        using var fixture = new MinimizationFixture();
        var stoppedWorker = new MinimizationWorker(fixture)
        { ChangeMinimization = reply => reply with { Standing = WorkerResultStanding.Stopped } };
        var stopped = await new MinimizationOwner(stoppedWorker).RunAsync(fixture.Constructed,
            fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Null(stopped.CompletedStage);
        Assert.Equal(StageExecutionStanding.Stopped, stopped.State.Standing);

        var wrongPolicyWorker = new MinimizationWorker(fixture);
        var wrongPolicy = await new MinimizationOwner(wrongPolicyWorker).RunAsync(fixture.Constructed,
            fixture.Policy with { Id = "other" }, fixture.Directory, null,
            TestContext.Current.CancellationToken);
        Assert.Null(wrongPolicy.CompletedStage);
        Assert.Empty(wrongPolicyWorker.MinimizationRequests);

        foreach (var changed in new[]
                 {
                     fixture.Policy with { Version = "2" },
                     fixture.Policy with { MaximumMinimizationIterations = 101 },
                     fixture.Policy with { Construction = fixture.Policy.Construction with
                         { NativePatchSha256 = new string('b', 64) } },
                     fixture.Policy with { ForceFieldFiles = ImmutableArray.Create(
                         new ForceFieldAsset("ff19SB", "19", "Amber19", "changed.xml", new string('a', 64))) }
                 })
        {
            var worker = new MinimizationWorker(fixture);
            var result = await new MinimizationOwner(worker).RunAsync(fixture.Constructed,
                changed, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Empty(worker.MinimizationRequests);
        }
        var alteredBindingWorker = new MinimizationWorker(fixture);
        var alteredBinding = await new MinimizationOwner(alteredBindingWorker).RunAsync(
            fixture.Constructed with { Attempt = fixture.Attempt with
                { PolicyFingerprintSha256 = new string('0', 64) } },
            fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Null(alteredBinding.CompletedStage);
        Assert.Empty(alteredBindingWorker.MinimizationRequests);
    }

    [Fact]
    public void Complete_policy_fingerprint_is_independent_of_dictionary_insertion_order()
    {
        using var fixture = new MinimizationFixture();
        var first = ImmutableDictionary<string, double>.Empty.Add("N", 1.5).Add("C", 1.7);
        var reverse = ImmutableDictionary<string, double>.Empty.Add("C", 1.7).Add("N", 1.5);
        var left = fixture.Policy with { LocalStateObservation = fixture.Policy.LocalStateObservation with
            { AtomRadiusByElementAngstrom = first } };
        var right = fixture.Policy with { LocalStateObservation = fixture.Policy.LocalStateObservation with
            { AtomRadiusByElementAngstrom = reverse } };
        Assert.Equal(PreparationPolicyFingerprint.Compute(left), PreparationPolicyFingerprint.Compute(right));
        Assert.NotEqual(PreparationPolicyFingerprint.Compute(left),
            PreparationPolicyFingerprint.Compute(right with { MaximumMinimizationIterations = 101 }));
        Assert.NotEqual(PreparationPolicyFingerprint.Compute(left),
            PreparationPolicyFingerprint.Compute(left with { Construction = left.Construction with
                { MinimumPaddingNanometers = 2 } }));
        Assert.NotEqual(PreparationPolicyFingerprint.Compute(left),
            PreparationPolicyFingerprint.Compute(left with { Construction = left.Construction with
                { NativePatchSha256 = new string('b', 64) } }));
    }
}

internal sealed class MinimizationFixture : IDisposable
{
    public string Directory { get; } = Path.Combine(Path.GetTempPath(), "pims-minimization-" + Guid.NewGuid().ToString("N"));
    public PreparationAttempt Attempt { get; }
    public ConstructedExplicitSystem Constructed { get; }
    public ApplicablePreparationPolicy Policy { get; }
    public string MinimizedCoordinatePath { get; }
    public string MinimizedCoordinateSha { get; }
    public string MinimizedStatePath { get; }
    public string MinimizedStateSha { get; }

    public MinimizationFixture()
    {
        System.IO.Directory.CreateDirectory(Directory);
        string Write(string name) { var path = Path.Combine(Directory, name); File.WriteAllText(path, name); return path; }
        static string Sha(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();
        var coordinates = Write("constructed.cif");
        var topology = Write("topology.json");
        var system = Write("system.xml");
        var state = Write("state.xml");
        var correspondence = Write("correspondence.json");
        MinimizedCoordinatePath = Write("minimized.cif");
        MinimizedCoordinateSha = Sha(MinimizedCoordinatePath);
        MinimizedStatePath = Write("minimized-state.xml");
        MinimizedStateSha = Sha(MinimizedStatePath);
        var artifact = new MolecularArtifact("constructed-1", coordinates, Sha(coordinates),
            topology, system, state, 3, "20 x 20 x 40", Sha(topology), correspondence,
            Sha(correspondence), Sha(system), Sha(state));
        var mapping = new SourceToResultCorrespondence("oriented-source", artifact.Id,
            ImmutableArray<AtomCorrespondence>.Empty, true);
        Policy = new ApplicablePreparationPolicy("preparation-policy-1", "1",
            "canonical-amino-acid-assembly", ImmutableArray.Create("controlled policy evidence"),
            ImmutableArray<ForceFieldAsset>.Empty, 100, 10, 0.01, 0.01, 0.01,
            new MolecularDynamicsSystemSettings("PME", 1, "HBonds", true, 0.0005, null, true, true, null),
            new ConstructionPolicy("construction-1", "1", ImmutableArray.Create("controlled evidence"),
                ImmutableArray.Create(ProteinTopologyKind.MembraneSpanning), ImmutableArray.Create("DMPC"),
                "OpenMM Modeller.addMembrane", "8.6.0.dev-c6173db", "patch.pdb",
                new string('a', 64), "DMPC", "Na+", "Cl-", 1, 55.4,
                100000, 100, 900, "controlled native derivation"),
            null!, null!, null!,
            new LocalStateObservationSpec(ImmutableArray<LocalContactRolePair>.Empty,
                ImmutableDictionary<string, double>.Empty, 5, 10, true, ImmutableArray<string>.Empty),
            ImmutableArray<LocalStateCriterion>.Empty, ImmutableArray<LocalRolePairCriterion>.Empty,
            ImmutableArray<PreparationAssessmentCriterion>.Empty,
            new ProteinGeometryMeasurementSpec(ImmutableArray<string>.Empty,
                ImmutableDictionary<string, double>.Empty, 5, 1, 10),
            ImmutableArray<StageProteinGeometryCriterion>.Empty, null, ImmutableArray<string>.Empty);
        Attempt = new PreparationAttempt("attempt-1", "revision-1", "protein-1", "membrane-1",
            "placement-1", Policy.Id, DateTimeOffset.UtcNow, Policy.Version,
            PreparationPolicyFingerprint.Compute(Policy), Policy.ForceFieldFiles,
            Policy.Construction.ProviderVersion, Policy.Construction.NativePatchSha256);
        Constructed = new ConstructedExplicitSystem(artifact.Id, Attempt, artifact,
            new ConstructionDerivation(Attempt.Id, ImmutableArray<SpeciesCount>.Empty,
                ImmutableArray.Create(20.0, 20.0, 40.0),
                0, 0, 0, 0, 0.15, 0, 0, ImmutableArray<string>.Empty, ImmutableArray<string>.Empty),
            mapping, ImmutableArray<SpeciesCount>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<ScientificFinding>.Empty, "untreated");
    }

    public void Dispose() => System.IO.Directory.Delete(Directory, true);
}

internal sealed class MinimizationWorker(MinimizationFixture fixture) : IMinimizationWork
{
    public List<ScientificWorkRequest<MinimizationPayload>> MinimizationRequests { get; } = [];
    public List<ScientificWorkRequest<StageObservationPayload>> ObservationRequests { get; } = [];
    public Func<WorkerResult<MinimizationObservations>, WorkerResult<MinimizationObservations>>? ChangeMinimization { get; init; }
    public Func<WorkerResult<StageObservationObservations>, WorkerResult<StageObservationObservations>>? ChangeObservation { get; init; }

    public Task<WorkerResult<MinimizationObservations>> MinimizeAsync(
        ScientificWorkRequest<MinimizationPayload> request, CancellationToken cancellationToken)
    {
        MinimizationRequests.Add(request);
        var result = new WorkerResult<MinimizationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed,
            ImmutableArray.Create(new WorkerArtifact("minimizedCif", fixture.MinimizedCoordinatePath,
                    fixture.MinimizedCoordinateSha),
                new WorkerArtifact("minimizedStateXml", fixture.MinimizedStatePath, fixture.MinimizedStateSha)),
            new MinimizationObservations(-100, -110, 8, 19, StageTermination.Converged,
                true, fixture.Constructed.Molecule.AtomCount, ImmutableArray<string>.Empty,
                FinalRawRmsForceKjMolNm: 500,
                MaximumRelativeConstraintError: 0.0000001,
                AppliedConstraintTolerance: 0.00001,
                FinalRmsForceMethod: "constraint-tangent-per-particle-v1"),
            new ProviderIdentity("OpenMM", "8.6"), null, null);
        return Task.FromResult(ChangeMinimization?.Invoke(result) ?? result);
    }

    public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken)
    {
        ObservationRequests.Add(request);
        var local = new LocalStateObservations(ObservationStanding.Unavailable,
            "local state not part of this factual completion fixture", ImmutableArray<MeasuredValue>.Empty,
            ImmutableArray<LocalContactObservation>.Empty, ImmutableArray<LocalContactRolePair>.Empty,
            ImmutableArray<string>.Empty, ImmutableArray<LocalRolePairMeasurement>.Empty);
        var geometry = new ProteinGeometryObservations(ObservationStanding.Unavailable,
            ImmutableArray<ProteinGeometryKindObservation>.Empty,
            ImmutableArray<GeometryDistanceObservation>.Empty, ImmutableArray<string>.Empty);
        var result = new WorkerResult<StageObservationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed, ImmutableArray<WorkerArtifact>.Empty,
            new StageObservationObservations(fixture.Constructed.Molecule.AtomCount, true, true,
                ImmutableArray<MeasuredValue>.Empty, ImmutableArray<string>.Empty,
                ImmutableArray<string>.Empty, ImmutableArray<string>.Empty, local, geometry),
            new ProviderIdentity("OpenMM", "8.6"), null, null);
        return Task.FromResult(ChangeObservation?.Invoke(result) ?? result);
    }
}
