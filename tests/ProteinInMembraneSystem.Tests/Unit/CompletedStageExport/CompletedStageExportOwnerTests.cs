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
    public async Task Established_qualified_assessment_is_exported_without_a_new_scientific_judgment()
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
        Assert.Equal(PreparationQualification.QualifiedPrepared, assessment.Qualification);
        Assert.True(assessment.CurrentlyApplicable);

        var worker = new ControlledExportReadBack(assessed.Constructed.Molecule.AtomCount);
        var owner = new ExportOwner(worker);
        var result = await owner.ExportAsync(stage, assessment, basis.Revision,
            assessed.Protein, assessed.Membrane, assessed.Placement,
            assessed.Constructed, assessed.Constructed.Derivation, assessed.Policy,
            ImmutableArray<ResearcherDecision>.Empty, null, null, null,
            Path.Combine(basis.Directory, "qualified-export"), TestContext.Current.CancellationToken);

        Assert.True(result.Established, result.Reason);
        Assert.Single(worker.Requests);
        Assert.Equal(stage.Id, worker.Requests[0].Payload.StageId);
        using var archive = ZipFile.OpenRead(result.Value!.BundlePath);
        using var manifest = JsonDocument.Parse(archive.GetEntry("manifest.json")!.Open());
        var declared = manifest.RootElement.GetProperty("assessment");
        Assert.Equal(assessment.Id, declared.GetProperty("id").GetString());
        Assert.Equal(stage.Id, declared.GetProperty("stageId").GetString());
        Assert.Equal("QualifiedPrepared", declared.GetProperty("qualification").GetString());
        Assert.Equal(assessment.Reason, declared.GetProperty("reason").GetString());
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
        Assert.Equal("NotQualified",
            root.GetProperty("assessment").GetProperty("qualification").GetString());
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
            Assert.True(started.State.Standing == StageExecutionStanding.ReadyForMinimization,
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
                PreparationQualification.Indeterminate);
            EquilibratedAssessment = Assessment(Equilibrated, "equilibrated-assessment",
                PreparationQualification.NotQualified);
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
            PreparationQualification qualification) => new(id, stage.Id, qualification,
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
