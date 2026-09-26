using System.Collections.Immutable;
using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text.Json;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using EquilibrationOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.OptionalEquilibrationProcedure;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class OptionalEquilibrationOwnerTests
{
    [Fact]
    public async Task Completed_optional_stage_uses_exact_protocol_and_its_own_periodic_topology()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var worker = new OptionalEquilibrationWorker(fixture);
        var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized, fixture.Policy,
            fixture.Directory, null, TestContext.Current.CancellationToken);

        Assert.True(result.State.Standing == StageExecutionStanding.Completed, result.State.Message);
        var stage = Assert.IsType<CompletedStage>(result.CompletedStage);
        var sent = Assert.Single(worker.EquilibrationRequests).Payload;
        var observed = Assert.Single(worker.ObservationRequests).Payload;
        Assert.Equal(fixture.Protocol.Id, sent.ValidatedPolicyId);
        Assert.Equal(fixture.Minimized.Id, sent.SourceMinimizedStageId);
        Assert.Equal(fixture.Minimized.Molecule.TopologyPath, sent.TopologyJsonPath);
        Assert.Equal(fixture.EquilibratedTopologyPath, observed.TopologyJsonPath);
        Assert.Equal(fixture.EquilibratedTopologySha, observed.TopologyJsonSha256);
        Assert.Equal(fixture.EquilibratedTopologyPath, stage.Molecule.TopologyPath);
        Assert.Equal(fixture.EquilibratedTopologySha, stage.Molecule.TopologySha256);
        Assert.Null(stage.Molecule.CellDescription);
        Assert.Equal(fixture.Minimized.Id, stage.SourceStageId);
        Assert.Equal(fixture.Minimized.Attempt.Id, stage.Attempt.Id);
        Assert.Equal(stage.Id, stage.Correspondence.ResultId);
        Assert.Equal(StageKind.Equilibration, stage.Kind);
        Assert.Equal(EquilibrationObservationAdequacy.Adequate,
            stage.Observation.ObservationAdequacy);
        Assert.Equal("unrestrained", Assert.Single(stage.Observation.EquilibrationWindows).Name);
        var frames = Assert.IsType<EquilibrationFrameSeries>(stage.Observation.FrameSeries);
        Assert.Equal(3, frames.FrameCount);
        Assert.Equal(3 * fixture.Minimized.Molecule.AtomCount * 24, frames.PositionsByteLength);
        Assert.Equal(ConstructionFixture.Hash(frames.PositionsPath), frames.PositionsSha256);
        Assert.Equal(ConstructionFixture.Hash(frames.ManifestPath), frames.ManifestSha256);
        using var manifest = JsonDocument.Parse(File.ReadAllText(frames.ManifestPath));
        Assert.Equal(stage.Id, manifest.RootElement.GetProperty("stageId").GetString());
        Assert.Equal(fixture.Minimized.Id,
            manifest.RootElement.GetProperty("sourceMinimizedStageId").GetString());
        Assert.Equal(EquilibrationProtocolFingerprint.Compute(fixture.Protocol),
            manifest.RootElement.GetProperty("protocolSha256").GetString());
    }

    [Fact]
    public async Task Ramp_and_restraint_selectors_bind_exact_protocol_and_protein_atom_groups()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var observation = fixture.Protocol.Stages[0];
        var heat = observation with
        {
            Name = "heat", Steps = 2, InitialTemperatureKelvin = 100,
            ProteinRestraintSelector = "protein-heavy",
            ProteinRestraintKjMolNm2 = 2092, LipidRestraintKjMolNm2 = 2092
        };
        var release = observation with
        {
            Name = "ca-release", ProteinRestraintSelector = "protein-ca",
            ProteinRestraintKjMolNm2 = 418.4
        };
        var protocol = fixture.Protocol with
        {
            Stages = ImmutableArray.Create(heat, release, observation),
            MaximumSampleCount = 8,
            MaximumFrameBytes = 8 * 6 * 24
        };
        Assert.NotEqual(EquilibrationProtocolFingerprint.Compute(fixture.Protocol),
            EquilibrationProtocolFingerprint.Compute(protocol));
        Assert.NotEqual(EquilibrationProtocolFingerprint.Compute(protocol),
            EquilibrationProtocolFingerprint.Compute(protocol with
            { Stages = protocol.Stages.SetItem(0, heat with { InitialTemperatureKelvin = 110 }) }));
        Assert.NotEqual(EquilibrationProtocolFingerprint.Compute(protocol),
            EquilibrationProtocolFingerprint.Compute(protocol with
            { Stages = protocol.Stages.SetItem(1, release with
                { ProteinRestraintSelector = "backbone-heavy" }) }));

        var policy = fixture.Policy with { OptionalEquilibration = protocol };
        var sidechain = fixture.Minimized.Correspondence.Atoms[0] with
        { ResultAtomIndex = 5, ResultAtomId = "result:5", AtomRole = AtomRoleKind.Sidechain };
        var minimized = fixture.Minimized with
        {
            Attempt = fixture.Minimized.Attempt with
            { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy) },
            Molecule = fixture.Minimized.Molecule with { AtomCount = 6 },
            Correspondence = fixture.Minimized.Correspondence with
            { Atoms = fixture.Minimized.Correspondence.Atoms.Add(sidechain) }
        };
        var worker = new OptionalEquilibrationWorker(fixture);
        var result = await new EquilibrationOwner(worker).RunAsync(minimized, policy,
            fixture.Directory, null, TestContext.Current.CancellationToken);

        // This fixture's one-window worker response is deliberately too short for
        // the new three-window request; it must not establish a completed stage.
        Assert.Null(result.CompletedStage);
        Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
        var sent = Assert.Single(worker.EquilibrationRequests).Payload;
        Assert.Equal(new[] { 0 }, sent.ProteinBackboneAtomIndices);
        Assert.Equal(new[] { 0, 5 }, sent.ProteinHeavyAtomIndices);
        Assert.Equal(new[] { 1, 2 }, sent.LipidHeavyAtomIndices);
        Assert.Equal(100, sent.Protocol.Stages[0].InitialTemperatureKelvin);
        Assert.Equal("protein-heavy", sent.Protocol.Stages[0].ProteinRestraintSelector);
        Assert.Equal("protein-ca", sent.Protocol.Stages[1].ProteinRestraintSelector);
        Assert.Empty(worker.ObservationRequests);
    }

    [Fact]
    public async Task Ramp_with_pressure_or_unknown_restraint_selector_is_refused_before_worker()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var observation = fixture.Protocol.Stages[0];
        var invalidControls = new[]
        {
            observation with { Name = "invalid-ramp", InitialTemperatureKelvin = 100,
                PressureBar = 1, PressureMode = "xyisotropiczfree",
                BarostatFrequencySteps = 25, SurfaceTensionBarNm = 0 },
            observation with { Name = "invalid-selector", ProteinRestraintSelector = "unspecified" }
        };
        foreach (var invalid in invalidControls)
        {
            var protocol = fixture.Protocol with
            { Stages = ImmutableArray.Create(invalid, observation), MaximumSampleCount = 6 };
            var policy = fixture.Policy with { OptionalEquilibration = protocol };
            var minimized = fixture.Minimized with
            { Attempt = fixture.Minimized.Attempt with
                { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy) } };
            var worker = new OptionalEquilibrationWorker(fixture);
            var result = await new EquilibrationOwner(worker).RunAsync(minimized, policy,
                fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
            Assert.Empty(worker.EquilibrationRequests);
        }
    }

    [Fact]
    public async Task Frame_budget_is_bound_to_protocol_and_refused_before_worker()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var protocol = fixture.Protocol with
        { MaximumFrameBytes = 3 * fixture.Minimized.Molecule.AtomCount * 24 - 1 };
        Assert.NotEqual(EquilibrationProtocolFingerprint.Compute(fixture.Protocol),
            EquilibrationProtocolFingerprint.Compute(protocol));
        var policy = fixture.Policy with { OptionalEquilibration = protocol };
        var minimized = fixture.Minimized with { Attempt = fixture.Minimized.Attempt with
            { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy) } };
        var worker = new OptionalEquilibrationWorker(fixture);
        var result = await new EquilibrationOwner(worker).RunAsync(minimized, policy,
            fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Null(result.CompletedStage);
        Assert.Equal(StageExecutionStanding.ResourceRefused, result.State.Standing);
        Assert.Empty(worker.EquilibrationRequests);
    }

    [Fact]
    public async Task Rehashed_frame_data_or_manifest_identity_cannot_establish_stage()
    {
        foreach (var corruption in new[] { "rehashed frame data", "rehashed manifest step",
                     "rehashed final frame cell" })
        {
            using var fixture = new OptionalEquilibrationFixture();
            var worker = new OptionalEquilibrationWorker(fixture)
            {
                ChangeEquilibration = reply =>
                {
                    var artifacts = reply.Artifacts.ToArray();
                    var role = corruption == "rehashed frame data"
                        ? "sampledPositionsF64" : "sampledFramesManifest";
                    var index = Array.FindIndex(artifacts, item => item.Role == role);
                    Assert.True(index >= 0);
                    var path = artifacts[index].Path;
                    if (role == "sampledPositionsF64")
                    {
                        var bytes = File.ReadAllBytes(path);
                        bytes[0] ^= 1;
                        File.WriteAllBytes(path, bytes);
                    }
                    else
                    {
                        using var document = JsonDocument.Parse(File.ReadAllText(path));
                        var root = document.RootElement.Deserialize<Dictionary<string, JsonElement>>()!;
                        var frames = root["frames"].EnumerateArray()
                            .Select(item => item.Deserialize<Dictionary<string, JsonElement>>()!).ToArray();
                        if (corruption == "rehashed manifest step")
                            frames[1]["step"] = JsonSerializer.SerializeToElement(99);
                        else
                            frames[2]["boxVectorsAngstrom"] = JsonSerializer.SerializeToElement(
                                new[] { new[] { 31.0, 0.0, 0.0 }, new[] { 0.0, 30.0, 0.0 },
                                    new[] { 0.0, 0.0, 30.0 } });
                        root["frames"] = JsonSerializer.SerializeToElement(frames);
                        File.WriteAllText(path, JsonSerializer.Serialize(root));
                    }
                    artifacts[index] = artifacts[index] with { Sha256 = ConstructionFixture.Hash(path) };
                    return reply with { Artifacts = artifacts.ToImmutableArray() };
                }
            };
            var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
                fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Equal(StageExecutionStanding.Unobserved, result.State.Standing);
            Assert.Empty(worker.ObservationRequests);
        }
    }

    [Fact]
    public async Task Altered_attempt_policy_or_minimized_stage_cannot_start_optional_work()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var cases = new (string Name, CompletedStage Stage, ApplicablePreparationPolicy Policy)[]
        {
            ("attempt fingerprint", fixture.Minimized with { Attempt = fixture.Minimized.Attempt with
                { PolicyFingerprintSha256 = new string('0', 64) } }, fixture.Policy),
            ("different policy protocol", fixture.Minimized, fixture.Policy with
                { OptionalEquilibration = fixture.Protocol with { RandomSeed = 99 } }),
            ("invalid NVT pressure mode", fixture.Minimized with { Attempt = fixture.Minimized.Attempt with
                { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(fixture.Policy with
                    { OptionalEquilibration = fixture.Protocol with { Stages =
                        fixture.Protocol.Stages.SetItem(0, fixture.Protocol.Stages[0] with
                        { PressureMode = "constantVolume" }) } }) } }, fixture.Policy with
                { OptionalEquilibration = fixture.Protocol with { Stages =
                    fixture.Protocol.Stages.SetItem(0, fixture.Protocol.Stages[0] with
                    { PressureMode = "constantVolume" }) } }),
            ("wrong observed source", fixture.Minimized with { Observation =
                fixture.Minimized.Observation with { StageId = "other-stage" } }, fixture.Policy),
            ("wrong correspondence", fixture.Minimized with { Correspondence =
                fixture.Minimized.Correspondence with { ResultId = "other-stage" } }, fixture.Policy)
        };
        foreach (var (name, stage, policy) in cases)
        {
            var worker = new OptionalEquilibrationWorker(fixture);
            var result = await new EquilibrationOwner(worker).RunAsync(stage, policy,
                fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
            Assert.Empty(worker.EquilibrationRequests);
            Assert.Empty(worker.ObservationRequests);
        }
    }

    [Fact]
    public async Task Missing_or_altered_final_state_parts_cannot_be_observed_or_promoted()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var roles = new[] { "equilibratedCif", "equilibratedTopologyJson",
            "equilibratedSystemXml", "equilibratedStateXml" };
        foreach (var role in roles)
        {
            foreach (var change in new[] { "missing", "changed digest", "missing bytes" })
            {
                var worker = new OptionalEquilibrationWorker(fixture)
                {
                    ChangeEquilibration = reply => reply with
                    {
                        Artifacts = change switch
                        {
                            "missing" => reply.Artifacts.Where(item => item.Role != role).ToImmutableArray(),
                            "changed digest" => reply.Artifacts.Select(item => item.Role == role
                                ? item with { Sha256 = new string('0', 64) } : item).ToImmutableArray(),
                            _ => reply.Artifacts.Select(item => item.Role == role
                                ? item with { Path = item.Path + ".absent" } : item).ToImmutableArray()
                        }
                    }
                };
                var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
                    fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
                Assert.Null(result.CompletedStage);
                Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
                Assert.Empty(worker.ObservationRequests);
            }
        }
    }

    [Fact]
    public async Task Completed_trace_must_match_declared_window_samples_and_finite_controls()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var changes = new (string Name, Func<EquilibrationObservations,
            EquilibrationObservations> Change)[]
        {
            ("missing sample", value => value with { Samples = value.Samples.RemoveAt(1) }),
            ("wrong step", value => value with { Samples = value.Samples.SetItem(1,
                value.Samples[1] with { Step = 3 }) }),
            ("nonfinite sampled temperature", value => value with { Windows = value.Windows.SetItem(0,
                value.Windows[0] with { TemperatureKelvin = double.NaN }) }),
            ("wrong pressure treatment", value => value with { Windows = value.Windows.SetItem(0,
                value.Windows[0] with { PressureBar = 1.0 }) }),
            ("invented effective blocks", value => value with
            { ObservationAssessments = value.ObservationAssessments.SetItem(0,
                value.ObservationAssessments[0] with { EffectiveBlockCount = 4 }) }),
            ("fabricated drift within bound", value => value with
            { ObservationAssessments = value.ObservationAssessments.SetItem(0,
                value.ObservationAssessments[0] with { FirstVsLastBlockMeanDifference = 0.5 }) }),
            ("fabricated correlation within bound", value => value with
            { ObservationAssessments = value.ObservationAssessments.SetItem(0,
                value.ObservationAssessments[0] with { LagOneBlockCorrelation = 0.5 }) }),
            ("changed sample with stale statistics", value => value with
            { Samples = value.Samples.SetItem(2, value.Samples[2] with
                { Measurements = value.Samples[2].Measurements.SetItem(0,
                    value.Samples[2].Measurements[0] with { Value = 2 }) }) }),
            ("fabricated insufficiency at bound", value => value with
            {
                ObservationAdequacy = EquilibrationObservationAdequacy.InsufficientAtBound,
                ObservationAssessments = value.ObservationAssessments.SetItem(0,
                    value.ObservationAssessments[0] with { Sufficient = false })
            })
        };
        foreach (var (name, change) in changes)
        {
            var worker = new OptionalEquilibrationWorker(fixture)
            {
                ChangeEquilibration = reply => reply with
                { Observations = change(reply.Observations!) }
            };
            var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
                fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Equal(StageExecutionStanding.Failed, result.State.Standing);
            Assert.Empty(worker.ObservationRequests);
        }
    }

    [Fact]
    public async Task Recomputed_nonconstant_block_trace_can_establish_adequate_observation()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var worker = new OptionalEquilibrationWorker(fixture)
        {
            ChangeEquilibration = reply =>
            {
                var observed = reply.Observations!;
                // Blocks [1, 2, 1] have zero first-last drift, lag-one
                // correlation -2/3, and three effective blocks.
                var middle = observed.Samples[1];
                var samples = observed.Samples.SetItem(1, middle with
                { Measurements = middle.Measurements.SetItem(0,
                    middle.Measurements[0] with { Value = 2 }) });
                var assessments = observed.ObservationAssessments.SetItem(0,
                    observed.ObservationAssessments[0] with { LagOneBlockCorrelation = -2.0 / 3.0 });
                return reply with { Observations = observed with
                { Samples = samples, ObservationAssessments = assessments } };
            }
        };
        var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
            fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Equal(StageExecutionStanding.Completed, result.State.Standing);
        Assert.NotNull(result.CompletedStage);
        Assert.Equal(-2.0 / 3.0,
            result.CompletedStage.Observation.EquilibrationAssessments[0].LagOneBlockCorrelation!.Value,
            precision: 12);
    }

    [Fact]
    public async Task Truthful_insufficiency_before_three_complete_blocks_remains_a_completed_observation()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var protocol = fixture.Protocol with
        {
            SufficiencyRules = fixture.Protocol.SufficiencyRules.Select(rule =>
                rule with { BlockSizeSamples = 2 }).ToImmutableArray()
        };
        var policy = fixture.Policy with { OptionalEquilibration = protocol };
        var minimized = fixture.Minimized with
        {
            Attempt = fixture.Minimized.Attempt with
            { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(policy) }
        };
        var worker = new OptionalEquilibrationWorker(fixture)
        {
            ChangeEquilibration = reply =>
            {
                var observed = reply.Observations!;
                return reply with { Observations = observed with
                {
                    ObservationAdequacy = EquilibrationObservationAdequacy.InsufficientAtBound,
                    ObservationAssessments = observed.ObservationAssessments.Select(assessment =>
                        assessment with { EffectiveBlockCount = 0,
                            FirstVsLastBlockMeanDifference = null,
                            LagOneBlockCorrelation = null, Sufficient = false }).ToImmutableArray()
                } };
            }
        };
        var result = await new EquilibrationOwner(worker).RunAsync(minimized,
            policy, fixture.Directory, null, TestContext.Current.CancellationToken);
        Assert.Equal(StageExecutionStanding.Completed, result.State.Standing);
        var stage = Assert.IsType<CompletedStage>(result.CompletedStage);
        Assert.Equal(EquilibrationObservationAdequacy.InsufficientAtBound,
            stage.Observation.ObservationAdequacy);
        Assert.All(stage.Observation.EquilibrationAssessments, assessment =>
        {
            Assert.Equal(0.0, assessment.EffectiveBlockCount);
            Assert.Null(assessment.FirstVsLastBlockMeanDifference);
            Assert.Null(assessment.LagOneBlockCorrelation);
            Assert.False(assessment.Sufficient);
        });
    }

    [Fact]
    public async Task Only_correlated_monotonic_declared_step_progress_reaches_the_active_stage()
    {
        using var fixture = new OptionalEquilibrationFixture();
        var progress = new StageProgressCollector();
        var worker = new OptionalEquilibrationWorker(fixture)
        {
            EmitProgress = (request, receiver) =>
            {
                var observed = new EquilibrationWorkProgress(request.Payload.StudyRevisionId,
                    request.Payload.AttemptId, request.Payload.StageId, "unrestrained", 1, 3);
                receiver!.Report(observed);
                receiver.Report(observed with { AttemptId = "another-attempt", CompletedSteps = 2 });
                receiver.Report(observed with { StageId = "another-stage", CompletedSteps = 2 });
                receiver.Report(observed with { Window = "unplanned", CompletedSteps = 2 });
                receiver.Report(observed with { CompletedSteps = 1 });
                receiver.Report(observed with { CompletedSteps = 2 });
            }
        };
        var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
            fixture.Policy, fixture.Directory, progress, TestContext.Current.CancellationToken);

        Assert.NotNull(result.CompletedStage);
        var stepUpdates = progress.States.Where(state => state.Message.Contains("steps observed",
            StringComparison.Ordinal)).ToArray();
        Assert.Equal(2, stepUpdates.Length);
        Assert.All(stepUpdates, state =>
        {
            Assert.Equal(result.CompletedStage!.Id, state.StageId);
            Assert.Equal(fixture.Minimized.Attempt.Id, state.AttemptId);
            Assert.Equal(StageExecutionStanding.Running, state.Standing);
        });
        Assert.Equal(1.0 / 3.0, stepUpdates[0].Progress);
        Assert.Equal(2.0 / 3.0, stepUpdates[1].Progress);
        Assert.Equal(StageExecutionStanding.Completed, progress.States[^1].Standing);
    }

    [Fact]
    public async Task Unobserved_or_stopped_remeasurement_retains_only_the_source_stage()
    {
        using var fixture = new OptionalEquilibrationFixture();
        foreach (var standing in new[] { WorkerResultStanding.Stopped,
            WorkerResultStanding.Unobserved, WorkerResultStanding.Failed })
        {
            var worker = new OptionalEquilibrationWorker(fixture)
            {
                ChangeObservation = reply => reply with { Standing = standing }
            };
            var result = await new EquilibrationOwner(worker).RunAsync(fixture.Minimized,
                fixture.Policy, fixture.Directory, null, TestContext.Current.CancellationToken);
            Assert.Null(result.CompletedStage);
            Assert.Equal(standing == WorkerResultStanding.Stopped ? StageExecutionStanding.Stopped :
                StageExecutionStanding.Unobserved, result.State.Standing);
            Assert.Single(worker.ObservationRequests);
            Assert.Equal(fixture.Minimized.Id,
                Assert.Single(worker.EquilibrationRequests).Payload.SourceMinimizedStageId);
        }
    }

    [Fact]
    public async Task Local_exchange_cancellation_is_stopped_at_either_optional_boundary()
    {
        using var fixture = new OptionalEquilibrationFixture();
        using var equilibrationCancellation = CancellationTokenSource.CreateLinkedTokenSource(
            TestContext.Current.CancellationToken);
        var equilibrationWorker = new OptionalEquilibrationWorker(fixture)
        {
            ChangeEquilibration = reply =>
            {
                equilibrationCancellation.Cancel();
                return reply with
                {
                    StudyRevisionId = null, AttemptId = null, StageId = null,
                    Standing = WorkerResultStanding.Stopped, Observations = null,
                    Artifacts = ImmutableArray<WorkerArtifact>.Empty,
                    FailureCode = "worker-stopped"
                };
            }
        };
        var equilibration = await new EquilibrationOwner(equilibrationWorker).RunAsync(
            fixture.Minimized, fixture.Policy, fixture.Directory, null,
            equilibrationCancellation.Token);
        Assert.Null(equilibration.CompletedStage);
        Assert.Equal(StageExecutionStanding.Stopped, equilibration.State.Standing);
        Assert.Empty(equilibrationWorker.ObservationRequests);

        using var observationCancellation = CancellationTokenSource.CreateLinkedTokenSource(
            TestContext.Current.CancellationToken);
        var observationWorker = new OptionalEquilibrationWorker(fixture)
        {
            ChangeObservation = reply =>
            {
                observationCancellation.Cancel();
                return reply with
                {
                    StudyRevisionId = null, AttemptId = null, StageId = null,
                    Standing = WorkerResultStanding.Stopped, Observations = null,
                    Artifacts = ImmutableArray<WorkerArtifact>.Empty,
                    FailureCode = "worker-stopped"
                };
            }
        };
        var observation = await new EquilibrationOwner(observationWorker).RunAsync(
            fixture.Minimized, fixture.Policy, fixture.Directory, null,
            observationCancellation.Token);
        Assert.Null(observation.CompletedStage);
        Assert.Equal(StageExecutionStanding.Stopped, observation.State.Standing);
        Assert.Single(observationWorker.ObservationRequests);
    }

    [Fact]
    public async Task Foreign_stopped_reply_cannot_claim_an_optional_stop_even_when_token_is_canceled()
    {
        using var fixture = new OptionalEquilibrationFixture();
        foreach (var mismatch in new[] { "stage", "request" })
        {
            using var equilibrationCancellation = CancellationTokenSource.CreateLinkedTokenSource(
                TestContext.Current.CancellationToken);
            var equilibrationWorker = new OptionalEquilibrationWorker(fixture)
            {
                ChangeEquilibration = reply =>
                {
                    equilibrationCancellation.Cancel();
                    return mismatch == "stage"
                        ? reply with { Standing = WorkerResultStanding.Stopped,
                            StageId = "foreign-stage", FailureCode = "worker-stopped" }
                        : reply with { RequestId = "foreign-request", StudyRevisionId = null,
                            AttemptId = null, StageId = null, Standing = WorkerResultStanding.Stopped,
                            FailureCode = "worker-stopped" };
                }
            };
            var equilibration = await new EquilibrationOwner(equilibrationWorker).RunAsync(
                fixture.Minimized, fixture.Policy, fixture.Directory, null,
                equilibrationCancellation.Token);
            Assert.Null(equilibration.CompletedStage);
            Assert.Equal(StageExecutionStanding.Unobserved, equilibration.State.Standing);

            using var observationCancellation = CancellationTokenSource.CreateLinkedTokenSource(
                TestContext.Current.CancellationToken);
            var observationWorker = new OptionalEquilibrationWorker(fixture)
            {
                ChangeObservation = reply =>
                {
                    observationCancellation.Cancel();
                    return mismatch == "stage"
                        ? reply with { Standing = WorkerResultStanding.Stopped,
                            StageId = "foreign-stage", FailureCode = "worker-stopped" }
                        : reply with { RequestId = "foreign-request", StudyRevisionId = null,
                            AttemptId = null, StageId = null, Standing = WorkerResultStanding.Stopped,
                            FailureCode = "worker-stopped" };
                }
            };
            var observation = await new EquilibrationOwner(observationWorker).RunAsync(
                fixture.Minimized, fixture.Policy, fixture.Directory, null,
                observationCancellation.Token);
            Assert.Null(observation.CompletedStage);
            Assert.Equal(StageExecutionStanding.Unobserved, observation.State.Standing);
        }
    }
}

internal sealed class OptionalEquilibrationFixture : IDisposable
{
    private readonly MinimizationFixture _basis = new();
    public string Directory => _basis.Directory;
    public EquilibrationProtocol Protocol { get; }
    public ApplicablePreparationPolicy Policy { get; }
    public CompletedStage Minimized { get; }
    public string EquilibratedCoordinatePath { get; }
    public string EquilibratedCoordinateSha { get; }
    public string EquilibratedTopologyPath { get; }
    public string EquilibratedTopologySha { get; }
    public string EquilibratedSystemPath { get; }
    public string EquilibratedSystemSha { get; }
    public string EquilibratedStatePath { get; }
    public string EquilibratedStateSha { get; }

    public OptionalEquilibrationFixture()
    {
        var observables = ImmutableArray.Create(
            Observable("potential", "kJ/mol", "system", "system", "potentialEnergy", "none", "none"),
            Observable("separation", "angstrom", "bilayer", "membraneOrganization",
                "centroidSeparationZ", "upper-lipid-head", "lower-lipid-head"),
            Observable("offset", "angstrom", "proteinVsBilayer", "proteinPlacement",
                "proteinMidplaneOffset", "protein-backbone", "upper-lipid-head", "lower-lipid-head"),
            Observable("nearWater", "count", "proteinOrLipid", "hydrationIons",
                "countWithinDistance", "water-oxygen", "protein-or-lipid-heavy", cutoff: 5),
            Observable("nearIons", "count", "proteinOrLipid", "hydrationIons",
                "countWithinDistance", "all-ions", "protein-or-lipid-heavy", cutoff: 5));
        var rules = observables.Select(item => new EquilibrationSufficiencyRule(item.Name,
            1, 3, 1, 0.9)).ToImmutableArray();
        var stage = new EquilibrationStageControl("unrestrained", 3, 0.002, 303,
            null, "none", null, null, 1, 0, 0, 1);
        Protocol = new EquilibrationProtocol("validated-optional-protocol", 303, 47,
            ImmutableArray.Create(stage), stage with { Name = "extension" }, 0, 3,
            observables.Select(item => item.Name).ToImmutableArray(), observables, rules,
            "controlled exact-class comparison", 3 * 5 * 24);
        Policy = _basis.Policy with { OptionalEquilibration = Protocol };
        var attempt = _basis.Attempt with
        { PolicyFingerprintSha256 = PreparationPolicyFingerprint.Compute(Policy) };
        var atoms = ImmutableArray.Create(
            Atom(0, MoleculeRoleKind.Protein, AtomRoleKind.Backbone, "C"),
            Atom(1, MoleculeRoleKind.Lipid, AtomRoleKind.Head, "C", LeafletSide.Upper),
            Atom(2, MoleculeRoleKind.Lipid, AtomRoleKind.Head, "C", LeafletSide.Lower),
            Atom(3, MoleculeRoleKind.Water, AtomRoleKind.Body, "O"),
            Atom(4, MoleculeRoleKind.Ion, AtomRoleKind.Body, "Na"));
        var id = "minimized-stage";
        var molecule = _basis.Constructed.Molecule with
        {
            Id = id, AtomCount = atoms.Length,
            CoordinatePath = _basis.MinimizedCoordinatePath,
            CoordinateSha256 = _basis.MinimizedCoordinateSha,
            StateXmlPath = _basis.MinimizedStatePath,
            StateXmlSha256 = _basis.MinimizedStateSha
        };
        var observation = new StageObservation(id, attempt.Id, StageKind.Minimization,
            ImmutableArray<MeasuredValue>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            StageTermination.Converged, "controlled provider", DateTimeOffset.UtcNow);
        Minimized = new CompletedStage(id, attempt, StageKind.Minimization, molecule, observation,
            new SourceToResultCorrespondence("source", id, atoms, true), Policy.Id,
            null, ImmutableArray<ScientificFinding>.Empty, DateTimeOffset.UtcNow);
        EquilibratedCoordinatePath = Write("equilibrated.cif", out var coordinatesSha);
        EquilibratedCoordinateSha = coordinatesSha;
        EquilibratedTopologyPath = Write("equilibrated-topology.json", out var topologySha,
            "{\"boxVectorsAngstrom\":[[30,0,0],[0,30,0],[0,0,30]]}");
        EquilibratedTopologySha = topologySha;
        EquilibratedSystemPath = Write("equilibrated-system.xml", out var systemSha);
        EquilibratedSystemSha = systemSha;
        EquilibratedStatePath = Write("equilibrated-state.xml", out var stateSha);
        EquilibratedStateSha = stateSha;
    }

    private string Write(string filename, out string sha, string? contents = null)
    {
        var path = Path.Combine(Directory, filename);
        File.WriteAllText(path, contents ?? "controlled " + filename);
        sha = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();
        return path;
    }

    private static EquilibrationObservable Observable(string name, string unit, string scope,
        string category, string method, string first, string second, string third = "none",
        double? cutoff = null) =>
        new(name, unit, scope, category, method, first, second, third, cutoff);

    private static AtomCorrespondence Atom(int index, MoleculeRoleKind role, AtomRoleKind atomRole,
        string element, LeafletSide? side = null) =>
        new(index, "result:" + index, null, index == 0 ? AtomOriginKind.Source : AtomOriginKind.Generated,
            role, atomRole, element, null, null, side);

    public void Dispose() => _basis.Dispose();
}

internal sealed class OptionalEquilibrationWorker(OptionalEquilibrationFixture fixture) : IOptionalEquilibrationWork
{
    public List<ScientificWorkRequest<EquilibrationPayload>> EquilibrationRequests { get; } = [];
    public List<ScientificWorkRequest<StageObservationPayload>> ObservationRequests { get; } = [];
    public Func<WorkerResult<EquilibrationObservations>, WorkerResult<EquilibrationObservations>>?
        ChangeEquilibration { get; init; }
    public Func<WorkerResult<StageObservationObservations>, WorkerResult<StageObservationObservations>>?
        ChangeObservation { get; init; }
    public Action<ScientificWorkRequest<EquilibrationPayload>,
        IProgress<EquilibrationWorkProgress>?>? EmitProgress { get; init; }

    public Task<WorkerResult<EquilibrationObservations>> EquilibrateAsync(
        ScientificWorkRequest<EquilibrationPayload> request,
        IProgress<EquilibrationWorkProgress>? progress, CancellationToken cancellationToken)
    {
        EquilibrationRequests.Add(request);
        EmitProgress?.Invoke(request, progress);
        var values = fixture.Protocol.Observables.Select(item =>
            new MeasuredValue(item.Name, 1, item.Unit, item.Scope)).ToImmutableArray();
        var samples = ImmutableArray.Create(1, 2, 3).Select(step =>
            new EquilibrationSample("unrestrained", step, values)).ToImmutableArray();
        var assessments = fixture.Protocol.Observables.Select(item =>
            new EquilibrationObservationAssessment(item.Name, 3, 3, 0, 0, true)).ToImmutableArray();
        var observed = new EquilibrationObservations(
            ImmutableArray.Create(new EquilibrationWindowObservation("unrestrained", 3, 3, 303,
                null, values, ImmutableArray<string>.Empty)),
            fixture.Minimized.Molecule.AtomCount, StageTermination.Completed, true,
            ImmutableArray<string>.Empty, 0, EquilibrationObservationAdequacy.Adequate,
            samples, assessments);
        Directory.CreateDirectory(request.WorkingDirectory);
        var positionsPath = Path.Combine(request.WorkingDirectory,
            $"sampled-positions-{request.Payload.StageId}.f64");
        var manifestPath = Path.Combine(request.WorkingDirectory,
            $"sampled-frames-{request.Payload.StageId}.json");
        var atomCount = fixture.Minimized.Molecule.AtomCount;
        var frameBytes = atomCount * 24;
        var data = new byte[frameBytes * samples.Length];
        var frames = samples.Select((sample, index) =>
        {
            var offset = index * frameBytes;
            for (var atom = 0; atom < atomCount; atom++)
                for (var axis = 0; axis < 3; axis++)
                    BinaryPrimitives.WriteDoubleLittleEndian(
                        data.AsSpan(offset + atom * 24 + axis * 8, 8), sample.Step + atom + axis / 10.0);
            var digest = Convert.ToHexString(SHA256.HashData(data.AsSpan(offset, frameBytes)))
                .ToLowerInvariant();
            return new
            {
                windowName = sample.WindowName, step = sample.Step,
                offsetBytes = (long)offset, lengthBytes = (long)frameBytes,
                positionsSha256 = digest,
                boxVectorsAngstrom = new[]
                {
                    new[] { 30.0, 0.0, 0.0 }, new[] { 0.0, 30.0, 0.0 },
                    new[] { 0.0, 0.0, 30.0 }
                }
            };
        }).ToArray();
        File.WriteAllBytes(positionsPath, data);
        File.WriteAllText(manifestPath, JsonSerializer.Serialize(new
        {
            schemaVersion = "protein-in-membrane.equilibration-frames.v1",
            encoding = "float64-le-xyz-angstrom",
            studyRevisionId = request.Payload.StudyRevisionId,
            attemptId = request.Payload.AttemptId,
            stageId = request.Payload.StageId,
            sourceMinimizedStageId = request.Payload.SourceMinimizedStageId,
            validatedPolicyId = request.Payload.ValidatedPolicyId,
            protocolSha256 = request.Payload.ProtocolSha256,
            sourceCoordinateSha256 = request.Payload.TopologyCifSha256,
            sourceTopologySha256 = request.Payload.TopologyJsonSha256,
            sourceSystemSha256 = request.Payload.SystemXmlSha256,
            sourceStateSha256 = request.Payload.MinimizedStateXmlSha256,
            finalCoordinateSha256 = fixture.EquilibratedCoordinateSha,
            finalTopologySha256 = fixture.EquilibratedTopologySha,
            finalSystemSha256 = fixture.EquilibratedSystemSha,
            finalStateSha256 = fixture.EquilibratedStateSha,
            atomCount, frameCount = samples.Length,
            maximumFrameBytes = request.Payload.Protocol.MaximumFrameBytes,
            positionsFile = Path.GetFileName(positionsPath),
            positionsSha256 = ConstructionFixture.Hash(positionsPath),
            positionsByteLength = (long)data.Length,
            frames
        }));
        var reply = new WorkerResult<EquilibrationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed,
            ImmutableArray.Create(
                new WorkerArtifact("equilibratedCif", fixture.EquilibratedCoordinatePath,
                    fixture.EquilibratedCoordinateSha),
                new WorkerArtifact("equilibratedTopologyJson", fixture.EquilibratedTopologyPath,
                    fixture.EquilibratedTopologySha),
                new WorkerArtifact("equilibratedSystemXml", fixture.EquilibratedSystemPath,
                    fixture.EquilibratedSystemSha),
                new WorkerArtifact("equilibratedStateXml", fixture.EquilibratedStatePath,
                    fixture.EquilibratedStateSha),
                new WorkerArtifact("sampledPositionsF64", positionsPath,
                    ConstructionFixture.Hash(positionsPath)),
                new WorkerArtifact("sampledFramesManifest", manifestPath,
                    ConstructionFixture.Hash(manifestPath))),
            observed, new ProviderIdentity("OpenMM", "controlled provider"), null, null);
        return Task.FromResult(ChangeEquilibration?.Invoke(reply) ?? reply);
    }

    public Task<WorkerResult<StageObservationObservations>> ObserveStageAsync(
        ScientificWorkRequest<StageObservationPayload> request, CancellationToken cancellationToken)
    {
        ObservationRequests.Add(request);
        var local = new LocalStateObservations(ObservationStanding.Unavailable,
            "controlled fixture", ImmutableArray<MeasuredValue>.Empty,
            ImmutableArray<LocalContactObservation>.Empty,
            ImmutableArray<LocalContactRolePair>.Empty, ImmutableArray<string>.Empty,
            ImmutableArray<LocalRolePairMeasurement>.Empty);
        var geometry = new ProteinGeometryObservations(ObservationStanding.Unavailable,
            ImmutableArray<ProteinGeometryKindObservation>.Empty,
            ImmutableArray<GeometryDistanceObservation>.Empty, ImmutableArray<string>.Empty);
        var reply = new WorkerResult<StageObservationObservations>(request.RequestId,
            request.Payload.StudyRevisionId, request.Payload.AttemptId, request.Payload.StageId,
            WorkerResultStanding.Observed, ImmutableArray<WorkerArtifact>.Empty,
            new StageObservationObservations(fixture.Minimized.Molecule.AtomCount, true, true,
                ImmutableArray<MeasuredValue>.Empty, ImmutableArray<string>.Empty,
                ImmutableArray<string>.Empty, ImmutableArray<string>.Empty, local, geometry),
            new ProviderIdentity("OpenMM", "controlled provider"), null, null);
        return Task.FromResult(ChangeObservation?.Invoke(reply) ?? reply);
    }
}

internal sealed class StageProgressCollector : IProgress<StageExecutionState>
{
    public List<StageExecutionState> States { get; } = [];
    public void Report(StageExecutionState value) => States.Add(value);
}
