using System.Collections.Immutable;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport;

/// <summary>Delivers one already completed stage with its own current assessment.</summary>
public sealed class CompletedStageExport
{
    private readonly ICompletedStageExportWork _worker;

    public CompletedStageExport(ICompletedStageExportWork worker) => _worker = worker;

    public async Task<BoundaryOutcome<CompletedStageBundle>> ExportAsync(
        CompletedStage stage,
        PreparationAssessmentResult assessment,
        StudyRevision originatingRevision,
        AssessedPreparedProtein protein,
        ApplicablePreparationPolicy policy,
        ImmutableArray<ResearcherDecision> decisions,
        ImmutableArray<ScientificFinding> currentFindings,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        if (assessment.StageId != stage.Id || !assessment.CurrentlyApplicable ||
            stage.Attempt.StudyRevisionId != originatingRevision.Id || stage.Observation.StageId != stage.Id ||
            stage.Observation.AttemptId != stage.Attempt.Id ||
            stage.PolicyId != policy.Id || stage.Attempt.PolicyId != policy.Id ||
            stage.Attempt.ProteinId != protein.Id ||
            !double.IsFinite(policy.ExportCoordinateReadBackToleranceAngstrom) ||
            policy.ExportCoordinateReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellLengthReadBackToleranceAngstrom) ||
            policy.ExportCellLengthReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellAngleReadBackToleranceDegrees) ||
            policy.ExportCellAngleReadBackToleranceDegrees < 0 ||
            stage.Correspondence.Atoms.IsDefaultOrEmpty || !stage.Correspondence.Complete ||
            stage.Correspondence.ResultId != stage.Molecule.Id ||
            stage.Correspondence.Atoms.Length != stage.Molecule.AtomCount ||
            stage.Correspondence.Atoms.Any(atom => atom.ApprovedChangeId is not null &&
                !decisions.Any(decision => decision.Id == atom.ApprovedChangeId &&
                    decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                    decision.ChosenValue == ResearcherDecisionValue.Approved &&
                    protein.Changes.Any(change => change.Id == decision.SubjectId &&
                        change.StudyRevisionId == decision.StudyRevisionId))) ||
            string.IsNullOrWhiteSpace(stage.Molecule.CoordinateSha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.TopologySha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.SystemXmlSha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.StateXmlSha256))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The selected completed stage has no current, corresponding assessment and complete molecular identity.");

        var molecule = stage.Molecule;
        if (string.IsNullOrWhiteSpace(molecule.TopologyPath) || string.IsNullOrWhiteSpace(molecule.SystemXmlPath) ||
            string.IsNullOrWhiteSpace(molecule.StateXmlPath) || !File.Exists(molecule.TopologyPath) ||
            !File.Exists(molecule.SystemXmlPath) || !File.Exists(molecule.StateXmlPath))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The completed stage's corresponding topology, System and State are unavailable.");

        var request = new ScientificWorkRequest<ExportVerificationPayload>(Guid.NewGuid().ToString("N"), workingDirectory,
            new ExportVerificationPayload(originatingRevision.Id, stage.Attempt.Id, stage.Id,
                molecule.CoordinatePath, molecule.CoordinateSha256,
                molecule.TopologyPath, molecule.TopologySha256!,
                molecule.SystemXmlPath, molecule.SystemXmlSha256!,
                molecule.StateXmlPath, molecule.StateXmlSha256!,
                policy.ExportCoordinateReadBackToleranceAngstrom,
                policy.ExportCellLengthReadBackToleranceAngstrom,
                policy.ExportCellAngleReadBackToleranceDegrees));
        var result = await _worker.VerifyExportAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != originatingRevision.Id ||
            result.AttemptId != stage.Attempt.Id || result.StageId != stage.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<CompletedStageBundle>.Unavailable(result.FailureMessage ?? "Export read-back was not observed.");

        var observed = result.Observations;
        var mmcif = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "stageMmcif");
        if (mmcif is null || !File.Exists(mmcif.Path) || observed.SourceAtomCount != molecule.AtomCount ||
            observed.ExportedAtomCount != molecule.AtomCount ||
            observed.CorrespondingElementAndResidueCount != molecule.AtomCount ||
            !observed.AtomOrderMatched || !observed.BondsMatched || !observed.CellMatched ||
            !observed.ReadBackMatched || !double.IsFinite(observed.CoordinateMaxDeviationAngstrom) ||
            observed.CoordinateMaxDeviationAngstrom > policy.ExportCoordinateReadBackToleranceAngstrom ||
            !observed.Warnings.IsDefaultOrEmpty)
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The all-atom coordinates, topology, state and cell did not pass mutual read-back verification.");

        Directory.CreateDirectory(workingDirectory);
        var bundlePath = Path.Combine(workingDirectory, $"completed-stage-{stage.Id}.zip");
        var temporaryPath = Path.Combine(workingDirectory, $"completed-stage-{stage.Id}-{Guid.NewGuid():N}.tmp");
        try
        {
            await using (var output = new FileStream(temporaryPath, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            using (var archive = new ZipArchive(output, ZipArchiveMode.Create, leaveOpen: false))
            {
                await AddVerifiedEntryAsync(archive, mmcif.Path, "structure.cif", mmcif.Sha256, cancellationToken);
                await AddVerifiedEntryAsync(archive, molecule.TopologyPath, "topology.json", molecule.TopologySha256!, cancellationToken);
                await AddVerifiedEntryAsync(archive, molecule.SystemXmlPath, "system.xml", molecule.SystemXmlSha256!, cancellationToken);
                await AddVerifiedEntryAsync(archive, molecule.StateXmlPath, "state.xml", molecule.StateXmlSha256!, cancellationToken);
                var manifest = archive.CreateEntry("manifest.json");
                await using var manifestStream = manifest.Open();
                await JsonSerializer.SerializeAsync(manifestStream, new
                {
                    schema = "protein-in-membrane.completed-stage.v1",
                    studyRevision = originatingRevision,
                    preparedProtein = new { protein.Id, protein.StudyRevisionId, protein.ChemicalStatePolicyId,
                        protein.Changes, protein.Molecule.CoordinateSha256, protein.Molecule.TopologySha256 },
                    attempt = stage.Attempt,
                    stage = new { stage.Id, stage.Kind, stage.PolicyId, stage.SourceStageId, stage.CompletedAt },
                    preparationPolicy = policy,
                    molecularIdentity = new { molecule.Id, molecule.AtomCount, molecule.CellDescription,
                        coordinatesSha256 = mmcif.Sha256, molecule.TopologySha256,
                        molecule.SystemXmlSha256, molecule.StateXmlSha256 },
                    correspondence = stage.Correspondence,
                    approvedPreparationDecisions = decisions.Where(decision =>
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved &&
                        protein.Changes.Any(change => change.Id == decision.SubjectId &&
                            change.StudyRevisionId == decision.StudyRevisionId)).ToArray(),
                    observation = stage.Observation,
                    assessment,
                    findings = currentFindings,
                    readBack = observed,
                    provider = result.Provider
                }, new JsonSerializerOptions(JsonSerializerDefaults.Web), cancellationToken);
            }
            if (File.Exists(bundlePath))
                return BoundaryOutcome<CompletedStageBundle>.Unavailable("This identified completed stage already has an export in the selected destination.");
            File.Move(temporaryPath, bundlePath);
            await using var bundle = File.OpenRead(bundlePath);
            var hash = Convert.ToHexString(await SHA256.HashDataAsync(bundle, cancellationToken)).ToLowerInvariant();
            return BoundaryOutcome<CompletedStageBundle>.Success(new CompletedStageBundle(
                stage.Id, originatingRevision.Id, stage.Attempt.Id, bundlePath, hash, DateTimeOffset.UtcNow));
        }
        catch (Exception exception) when (exception is IOException or InvalidDataException)
        {
            return BoundaryOutcome<CompletedStageBundle>.Unavailable($"The verified stage could not be delivered: {exception.Message}");
        }
        finally
        {
            if (File.Exists(temporaryPath)) File.Delete(temporaryPath);
        }
    }

    private static async Task AddVerifiedEntryAsync(ZipArchive archive, string path, string entryName,
        string expectedSha256, CancellationToken cancellationToken)
    {
        var entry = archive.CreateEntry(entryName);
        await using var source = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        await using var destination = entry.Open();
        using var digest = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        var buffer = new byte[64 * 1024];
        int read;
        while ((read = await source.ReadAsync(buffer, cancellationToken)) != 0)
        {
            digest.AppendData(buffer, 0, read);
            await destination.WriteAsync(buffer.AsMemory(0, read), cancellationToken);
        }
        var actualSha256 = Convert.ToHexString(digest.GetHashAndReset());
        if (!actualSha256.Equals(expectedSha256, StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException($"The {entryName} bytes being archived do not match the identified completed stage.");
    }
}
