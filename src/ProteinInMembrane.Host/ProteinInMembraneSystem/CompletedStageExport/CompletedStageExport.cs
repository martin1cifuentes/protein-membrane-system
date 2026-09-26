using System.Collections.Immutable;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport;

/// <summary>Delivers one already completed stage with its own current assessment.</summary>
public sealed class CompletedStageExport
{
    private const string Qualified6QwrSourceSha256 =
        "f4c1503a60321c0cfe513e8211e43e20e2199aa5f71f22429e0e0ae8c97b0779";
    private const string OpenMmTermsUri =
        "https://docs.openmm.org/latest/userguide/library/01_introduction.html#license";
    private readonly ICompletedStageExportWork _worker;

    public CompletedStageExport(ICompletedStageExportWork worker) => _worker = worker;

    public async Task<BoundaryOutcome<CompletedStageBundle>> ExportAsync(
        CompletedStage stage,
        PreparationAssessmentResult assessment,
        StudyRevision originatingRevision,
        AssessedPreparedProtein protein,
        AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement,
        ConstructedExplicitSystem constructed,
        ConstructionDerivation derivation,
        ApplicablePreparationPolicy policy,
        ImmutableArray<ResearcherDecision> decisions,
        string? ppmVersion,
        string? ppmExecutableSha256,
        string? optionalProtocolSha256,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        var equilibrated = stage.Kind == StageKind.Equilibration;
        if (stage.Kind is not (StageKind.Minimization or StageKind.Equilibration) ||
            stage.Observation.Kind != stage.Kind ||
            assessment.StageId != stage.Id || !assessment.CurrentlyApplicable ||
            stage.Attempt.StudyRevisionId != originatingRevision.Id || stage.Observation.StageId != stage.Id ||
            stage.Observation.AttemptId != stage.Attempt.Id ||
            stage.PolicyId != policy.Id || stage.Attempt.PolicyId != policy.Id ||
            stage.Attempt.PolicyVersion != policy.Version ||
            stage.Attempt.PolicyFingerprintSha256 != PreparationPolicyFingerprint.Compute(policy) ||
            stage.Attempt.ProteinId != protein.Id || protein.StudyRevisionId != originatingRevision.Id ||
            stage.Attempt.MembraneId != membrane.Id || membrane.StudyRevisionId != originatingRevision.Id ||
            stage.Attempt.PlacementId != placement.Id || placement.StudyRevisionId != originatingRevision.Id ||
            placement.Proposal.PreparedProteinId != protein.Id ||
            placement.Proposal.MembraneModelId != membrane.Intended.Id ||
            constructed.Attempt.Id != stage.Attempt.Id ||
            derivation.AttemptId != stage.Attempt.Id || constructed.Derivation != derivation ||
            stage.Attempt.NativePatchSha256 != policy.Construction.NativePatchSha256 ||
            stage.Attempt.ConstructionProviderVersion != policy.Construction.ProviderVersion ||
            !stage.Attempt.ForceFieldFiles.SequenceEqual(policy.ForceFieldFiles) ||
            !double.IsFinite(policy.ExportCoordinateReadBackToleranceAngstrom) ||
            policy.ExportCoordinateReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellLengthReadBackToleranceAngstrom) ||
            policy.ExportCellLengthReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellAngleReadBackToleranceDegrees) ||
            policy.ExportCellAngleReadBackToleranceDegrees < 0 ||
            stage.Correspondence.Atoms.IsDefaultOrEmpty || !stage.Correspondence.Complete ||
            stage.Correspondence.ResultId != stage.Molecule.Id ||
            stage.Correspondence.Atoms.Length != stage.Molecule.AtomCount ||
            stage.Correspondence.Atoms.Where((atom, index) => atom.ResultAtomIndex != index ||
                string.IsNullOrWhiteSpace(atom.ResultAtomId)).Any() ||
            stage.Correspondence.Atoms.Select(atom => atom.ResultAtomId)
                .Distinct(StringComparer.Ordinal).Count() != stage.Molecule.AtomCount ||
            stage.Correspondence.SourceId != constructed.Correspondence.SourceId ||
            stage.Correspondence.Atoms.Any(atom => atom.ApprovedChangeId is not null &&
                !decisions.Any(decision => decision.Id == atom.ApprovedChangeId &&
                    decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                    decision.ChosenValue == ResearcherDecisionValue.Approved &&
                    protein.Changes.Any(change => change.Id == decision.SubjectId &&
                        change.StudyRevisionId == decision.StudyRevisionId))) ||
            string.IsNullOrWhiteSpace(stage.Molecule.CoordinateSha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.CoordinatePath) ||
            string.IsNullOrWhiteSpace(stage.Molecule.TopologySha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.SystemXmlSha256) ||
            string.IsNullOrWhiteSpace(stage.Molecule.StateXmlSha256))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The selected completed stage has no current, corresponding assessment and complete molecular identity.");

        if (equilibrated &&
            (string.IsNullOrWhiteSpace(stage.SourceStageId) || stage.SourceStageId == stage.Id ||
             stage.Observation.Termination != StageTermination.Completed ||
             stage.Observation.ObservationAdequacy is null ||
             stage.Observation.EquilibrationSamples.IsDefaultOrEmpty ||
             stage.Observation.EquilibrationAssessments.IsDefaultOrEmpty ||
             stage.Observation.EquilibrationWindows.IsDefaultOrEmpty ||
             policy.OptionalEquilibration is null ||
             optionalProtocolSha256 is null || optionalProtocolSha256.Length != 64 ||
             !optionalProtocolSha256.All(Uri.IsHexDigit)))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable(
                "The selected equilibrated stage has no completed source-bound procedure and exact optional protocol.");

        var molecule = stage.Molecule;
        if (string.IsNullOrWhiteSpace(molecule.TopologyPath) || string.IsNullOrWhiteSpace(molecule.SystemXmlPath) ||
            string.IsNullOrWhiteSpace(molecule.StateXmlPath) || !File.Exists(molecule.TopologyPath) ||
            !File.Exists(molecule.SystemXmlPath) || !File.Exists(molecule.StateXmlPath))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The completed stage's corresponding topology, System and State are unavailable.");

        FinalCell? finalCell = null;
        if (equilibrated)
        {
            try
            {
                finalCell = await ReadFinalCellAsync(molecule.TopologyPath,
                    molecule.TopologySha256!, cancellationToken);
            }
            catch (Exception exception) when (exception is IOException or JsonException or
                                             InvalidDataException or UnauthorizedAccessException)
            {
                return BoundaryOutcome<CompletedStageBundle>.Unavailable(
                    $"The equilibrated stage's exact final cell is unavailable: {exception.Message}");
            }
        }

        cancellationToken.ThrowIfCancellationRequested();
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

        var bundlePath = Path.Combine(workingDirectory, $"completed-stage-{stage.Id}-{assessment.Id}.zip");
        var temporaryPath = Path.Combine(workingDirectory, $"completed-stage-{stage.Id}-{Guid.NewGuid():N}.tmp");
        try
        {
            Directory.CreateDirectory(workingDirectory);
            var contents = new[]
            {
                new BundleEntry("structure.cif", mmcif.Path, mmcif.Sha256),
                new BundleEntry("topology.json", molecule.TopologyPath, molecule.TopologySha256!),
                new BundleEntry("system.xml", molecule.SystemXmlPath, molecule.SystemXmlSha256!),
                new BundleEntry("state.xml", molecule.StateXmlPath, molecule.StateXmlSha256!)
            };
            await using (var output = new FileStream(temporaryPath, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            using (var archive = new ZipArchive(output, ZipArchiveMode.Create, leaveOpen: false))
            {
                foreach (var item in contents)
                    await AddVerifiedEntryAsync(archive, item.Path, item.Name, item.Sha256, cancellationToken);
                var manifest = archive.CreateEntry("manifest.json");
                await using var manifestStream = manifest.Open();
                await JsonSerializer.SerializeAsync(manifestStream, new
                {
                    schemaVersion = "protein-in-membrane.completed-stage.v1",
                    study = new { originatingRevision.Id, originatingRevision.Number,
                        originatingRevision.Conditions,
                        intendedProteinId = originatingRevision.IntendedProtein?.Id,
                        membraneModelId = originatingRevision.Membrane?.Id,
                        originatingRevision.AdoptedPlacementProposalId },
                    stage = new { stage.Id, stage.Kind, stage.PolicyId, stage.SourceStageId,
                        stage.CompletedAt, attemptId = stage.Attempt.Id,
                        moleculeId = molecule.Id, molecule.AtomCount,
                        molecule.CellDescription },
                    assessment,
                    attempt = new { stage.Attempt.Id, stage.Attempt.StudyRevisionId,
                        stage.Attempt.ProteinId, stage.Attempt.MembraneId,
                        stage.Attempt.PlacementId, stage.Attempt.PolicyId,
                        stage.Attempt.StartedAt, stage.Attempt.PolicyVersion,
                        stage.Attempt.PolicyFingerprintSha256,
                        stage.Attempt.ConstructionProviderVersion,
                        stage.Attempt.NativePatchSha256 },
                    preparedProtein = new
                    {
                        protein.Id, protein.StudyRevisionId,
                        source = new { protein.Intended.Source.Id, protein.Intended.Source.Kind,
                            protein.Intended.Source.Provenance,
                            accession = SourceAccession(protein, policy),
                            protein.Intended.Source.Sha256,
                            protein.Intended.Source.SourceModelDescription,
                            protein.Intended.Source.UploadProvenance,
                            protein.Intended.Source.UploadProvenanceNote },
                        protein.Intended.ModelIndex,
                        sourceModelNumber = protein.Intended.ModelIndex + 1,
                        protein.Intended.BiologicalAssemblyId,
                        protein.Intended.Chains, protein.Intended.Partners,
                        protein.Intended.AlternateLocations,
                        preparedCoordinateSha256 = protein.Molecule.CoordinateSha256,
                        preparedTopologySha256 = protein.Molecule.TopologySha256,
                        protein.ChemicalStatePolicyId, protein.ChemicalStatePolicyVersion,
                        protein.StructuralAssessmentPolicyId, protein.StructuralAssessmentPolicyVersion,
                        protein.ResidueVariants, protein.Changes, protein.Limitations,
                        protein.Evidence, protein.Findings
                    },
                    membrane = new { membrane.Id, membrane.StudyRevisionId,
                        membrane.Intended, membrane.PolicyId, membrane.PolicyVersion,
                        species = membrane.SpeciesRepresentations.Select(species => new
                        {
                            species.SpeciesId, species.ChemistryId, species.Category,
                            species.TemplateSha256, species.CoordinateTemplateSha256,
                            species.AtomCount, species.ForceFieldFamily, species.ForceFieldVersion,
                            species.Limitations
                        }).ToArray(), membrane.Evidence, membrane.Limitations },
                    placement = new { placement.Id, placement.StudyRevisionId,
                        placement.Standing, placement.Reason, placement.Findings,
                        proposal = new { placement.Proposal.Id, placement.Proposal.PreparedProteinId,
                            placement.Proposal.MembraneModelId, placement.Proposal.TopologyKind,
                            placement.Proposal.MidplaneAngstrom, placement.Proposal.ThicknessAngstrom,
                            placement.Proposal.TiltDegrees, placement.Proposal.PhysicalSide,
                            placement.Proposal.BiologicalSidedness, placement.Proposal.ContactingRegions,
                            placement.Proposal.Evidence, placement.Proposal.Limitations,
                            orientedCoordinateSha256 = placement.Proposal.OrientedProtein.CoordinateSha256,
                            orientedTopologySha256 = placement.Proposal.OrientedProtein.TopologySha256 } },
                    construction = new { constructed.Id, constructed.ConditionsTreatment,
                        constructed.AchievedComposition, constructed.ActualCellAngstrom,
                        constructed.Evidence, constructed.Findings,
                        constructed.LocalState, derivation,
                        constructedCoordinateSha256 = constructed.Molecule.CoordinateSha256,
                        constructedTopologySha256 = constructed.Molecule.TopologySha256,
                        provider = new { policy.Construction.ProviderName,
                            policy.Construction.ProviderVersion,
                            policy.Construction.NativePatchSha256,
                            policy.Construction.LipidTypeArgument,
                            policy.Construction.PositiveIonArgument,
                            policy.Construction.NegativeIonArgument,
                            policy.Construction.ApproximationStatement } },
                    forceFieldAssets = stage.Attempt.ForceFieldFiles.Select(asset => new
                    {
                        asset.Id, asset.Version, asset.Family, asset.Sha256
                    }).ToArray(),
                    preparationPolicy = new { policy.Id, policy.Version,
                        stage.Attempt.PolicyFingerprintSha256,
                        policy.EvidenceReferences, policy.SystemSettings,
                        policy.MaximumMinimizationIterations,
                        policy.FinalUnrestrainedRmsForceTargetKjMolNm,
                        policy.ExportCoordinateReadBackToleranceAngstrom,
                        policy.ExportCellLengthReadBackToleranceAngstrom,
                        policy.ExportCellAngleReadBackToleranceDegrees,
                        policy.Limitations },
                    optionalEquilibration = equilibrated ? new
                    {
                        sourceMinimizedStageId = stage.SourceStageId,
                        protocolId = policy.OptionalEquilibration!.Id,
                        protocolFingerprintSha256 = optionalProtocolSha256,
                        declaredProtocol = policy.OptionalEquilibration,
                        observedProcedure = new
                        {
                            stage.Observation.Termination,
                            stage.Observation.ObservationAdequacy,
                            windows = stage.Observation.EquilibrationWindows,
                            samples = stage.Observation.EquilibrationSamples,
                            assessments = stage.Observation.EquilibrationAssessments
                        }
                    } : null,
                    finalCell,
                    molecularIdentity = new { molecule.Id, molecule.AtomCount, molecule.CellDescription,
                        coordinatesSha256 = mmcif.Sha256, molecule.TopologySha256,
                        molecule.SystemXmlSha256, molecule.StateXmlSha256 },
                    correspondence = stage.Correspondence,
                    lineage = new
                    {
                        sourceId = protein.Intended.Source.Id,
                        sourceAccession = SourceAccession(protein, policy),
                        sourceCoordinateSha256 = protein.Intended.Source.Sha256,
                        sourceModelIndex = protein.Intended.ModelIndex,
                        sourceModelNumber = protein.Intended.ModelIndex + 1,
                        sourceBiologicalAssemblyId = protein.Intended.BiologicalAssemblyId,
                        intendedProteinId = protein.Intended.Id,
                        preparedProteinId = protein.Id,
                        preparedCoordinateSha256 = protein.Molecule.CoordinateSha256,
                        orientedCoordinateSha256 = placement.Proposal.OrientedProtein.CoordinateSha256,
                        constructedCoordinateSha256 = constructed.Molecule.CoordinateSha256,
                        completedCoordinateSha256 = molecule.CoordinateSha256,
                        approvedChangeIds = decisions.Where(decision =>
                            decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                            decision.ChosenValue == ResearcherDecisionValue.Approved &&
                            protein.Changes.Any(change => change.Id == decision.SubjectId))
                            .Select(decision => decision.SubjectId).ToArray(),
                        sourceToResult = stage.Correspondence
                    },
                    changes = protein.Changes,
                    approvedPreparationDecisions = decisions.Where(decision =>
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved &&
                        protein.Changes.Any(change => change.Id == decision.SubjectId &&
                            change.StudyRevisionId == decision.StudyRevisionId)).ToArray(),
                    observation = stage.Observation,
                    findings = assessment.Findings,
                    limitations = assessment.Limitations,
                    readBack = observed,
                    provider = result.Provider,
                    attribution = Attribution(protein, membrane, stage.Attempt, policy,
                        ppmVersion, ppmExecutableSha256),
                    artifacts = contents.Select(item => new
                    {
                        name = item.Name, sha256 = item.Sha256.ToLowerInvariant(),
                        byteLength = new FileInfo(item.Path).Length
                    }).ToArray()
                }, ManifestJson(), cancellationToken);
            }
            await ReadBackBundleAsync(temporaryPath, contents, stage.Id, assessment.Id, cancellationToken);
            string hash;
            await using (var bundle = File.OpenRead(temporaryPath))
                hash = Convert.ToHexString(await SHA256.HashDataAsync(bundle, cancellationToken)).ToLowerInvariant();
            var byteLength = new FileInfo(temporaryPath).Length;
            cancellationToken.ThrowIfCancellationRequested();
            // Same-directory replacement is atomic; a corrupted orphan from an interrupted
            // delivery can be repaired on retry without exposing a partial ZIP.
            File.Move(temporaryPath, bundlePath, overwrite: true);
            return BoundaryOutcome<CompletedStageBundle>.Success(new CompletedStageBundle(
                stage.Id, originatingRevision.Id, stage.Attempt.Id, assessment.Id,
                bundlePath, hash, byteLength, DateTimeOffset.UtcNow));
        }
        catch (Exception exception) when (exception is IOException or InvalidDataException or UnauthorizedAccessException)
        {
            return BoundaryOutcome<CompletedStageBundle>.Unavailable($"The verified stage could not be delivered: {exception.Message}");
        }
        finally
        {
            try
            {
                if (File.Exists(temporaryPath)) File.Delete(temporaryPath);
            }
            catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
            {
                // A leftover temporary file is never a published bundle.
            }
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

    private static async Task<FinalCell> ReadFinalCellAsync(string topologyPath,
        string expectedSha256, CancellationToken cancellationToken)
    {
        var bytes = await File.ReadAllBytesAsync(topologyPath, cancellationToken);
        var actualSha256 = Convert.ToHexString(SHA256.HashData(bytes));
        if (!actualSha256.Equals(expectedSha256, StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("The stage-specific topology bytes changed from their recorded digest.");
        using var document = JsonDocument.Parse(bytes);
        if (document.RootElement.ValueKind != JsonValueKind.Object ||
            !document.RootElement.TryGetProperty("boxVectorsAngstrom", out var box) ||
            box.ValueKind != JsonValueKind.Array || box.GetArrayLength() != 3)
            throw new InvalidDataException("The stage-specific topology has no three-vector periodic cell.");
        var vectors = new double[3][];
        for (var index = 0; index < 3; index++)
        {
            var vector = box[index];
            if (vector.ValueKind != JsonValueKind.Array || vector.GetArrayLength() != 3)
                throw new InvalidDataException("The stage-specific periodic cell has an invalid vector.");
            vectors[index] = new double[3];
            for (var axis = 0; axis < 3; axis++)
            {
                if (vector[axis].ValueKind != JsonValueKind.Number ||
                    !vector[axis].TryGetDouble(out var value) || !double.IsFinite(value))
                    throw new InvalidDataException("The stage-specific periodic cell has a nonfinite component.");
                vectors[index][axis] = value;
            }
        }
        static double Dot(double[] first, double[] second) =>
            first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
        var lengths = vectors.Select(vector => Math.Sqrt(Dot(vector, vector))).ToArray();
        if (lengths.Any(length => !double.IsFinite(length) || length <= 0))
            throw new InvalidDataException("The stage-specific periodic cell has a zero or nonfinite length.");
        static double Angle(double[] first, double[] second, double firstLength,
            double secondLength) => 180.0 / Math.PI * Math.Acos(Math.Clamp(
                Dot(first, second) / firstLength / secondLength, -1, 1));
        var angles = new[]
        {
            Angle(vectors[1], vectors[2], lengths[1], lengths[2]),
            Angle(vectors[0], vectors[2], lengths[0], lengths[2]),
            Angle(vectors[0], vectors[1], lengths[0], lengths[1])
        };
        var volume = Dot(vectors[0], [
            vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1],
            vectors[1][2] * vectors[2][0] - vectors[1][0] * vectors[2][2],
            vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0]
        ]);
        if (!double.IsFinite(volume) || volume <= 0 ||
            angles.Any(angle => !double.IsFinite(angle) || angle <= 0 || angle >= 180))
            throw new InvalidDataException("The stage-specific periodic cell is degenerate or inverted.");
        return new FinalCell(vectors, lengths, angles);
    }

    private static async Task ReadBackBundleAsync(string path, BundleEntry[] expected,
        string stageId, string assessmentId, CancellationToken cancellationToken)
    {
        await using var file = File.OpenRead(path);
        using var archive = new ZipArchive(file, ZipArchiveMode.Read, leaveOpen: false);
        if (archive.Entries.Count != expected.Length + 1 ||
            archive.Entries.Select(item => item.FullName).Distinct(StringComparer.Ordinal).Count() != archive.Entries.Count ||
            archive.GetEntry("manifest.json") is not { } manifest)
            throw new InvalidDataException("The completed-stage bundle did not read back as one complete archive.");
        foreach (var item in expected)
        {
            var entry = archive.GetEntry(item.Name);
            if (entry is null) throw new InvalidDataException($"The completed-stage bundle lacks {item.Name}.");
            await using var entryStream = entry.Open();
            var digest = Convert.ToHexString(await SHA256.HashDataAsync(entryStream, cancellationToken));
            if (!digest.Equals(item.Sha256, StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException($"The archived {item.Name} bytes changed during bundle read-back.");
        }
        await using var manifestStream = manifest.Open();
        using var document = await JsonDocument.ParseAsync(manifestStream, cancellationToken: cancellationToken);
        var root = document.RootElement;
        if (root.GetProperty("stage").GetProperty("id").GetString() != stageId ||
            root.GetProperty("assessment").GetProperty("id").GetString() != assessmentId ||
            root.GetProperty("assessment").GetProperty("stageId").GetString() != stageId ||
            root.GetProperty("artifacts").GetArrayLength() != expected.Length)
            throw new InvalidDataException("The bundle manifest did not preserve its stage and assessment identity.");
        foreach (var item in expected)
        {
            var declarations = root.GetProperty("artifacts").EnumerateArray()
                .Where(value => value.GetProperty("name").GetString() == item.Name)
                .Take(2).ToArray();
            if (declarations.Length != 1 ||
                !string.Equals(declarations[0].GetProperty("sha256").GetString(), item.Sha256,
                    StringComparison.OrdinalIgnoreCase) ||
                declarations[0].GetProperty("byteLength").GetInt64() != archive.GetEntry(item.Name)!.Length)
                throw new InvalidDataException($"The manifest does not bind the archived {item.Name} bytes.");
        }
    }

    private static object Attribution(AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        PreparationAttempt attempt, ApplicablePreparationPolicy policy,
        string? ppmVersion, string? ppmExecutableSha256)
    {
        var source = protein.Intended.Source;
        var exactQualified6Qwr = IsQualified6Qwr(protein, policy);
        var rcsb = source.Kind == SourceRouteKind.Rcsb && !string.IsNullOrWhiteSpace(source.Accession) ||
            exactQualified6Qwr;
        var sourceAccession = exactQualified6Qwr ? "6QWR" : source.Accession;
        var incorporated = new List<AttributionEntry>
        {
            new("protein source coordinates", sourceAccession ?? source.Id,
                source.SourceModelDescription ?? "identified source model",
                rcsb ? $"https://www.rcsb.org/structure/{sourceAccession}" : source.Provenance,
                source.Sha256,
                rcsb ? "RCSB PDB archive data: CC0" : "Source-specific rights were not established by this record",
                rcsb ? "https://www.rcsb.org/pages/usage-policy" : source.Provenance,
                exactQualified6Qwr
                    ? "10.2210/pdb6QWR/pdb; Schubeis et al., PNAS 2020, 10.1073/pnas.2002598117"
                    : sourceAccession ?? source.Provenance,
                "Derived coordinates are included; the original source file is identified by digest but is not bundled.",
                "incorporated and transformed"),
            new("native lipid coordinate patch", policy.Construction.LipidTypeArgument + " patch",
                policy.Construction.ProviderVersion,
                "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                attempt.NativePatchSha256,
                "OpenMM application layer MIT terms; no separate DMPC.pdb file terms are asserted",
                OpenMmTermsUri,
                "OpenMM application package, exact build identified by attempt construction provider version",
                "Lipid coordinates are derived into the bundle; the native source patch is identified, not bundled.",
                "incorporated and transformed")
        };
        var references = new List<AttributionEntry>();
        foreach (var representation in membrane.SpeciesRepresentations)
            references.Add(new AttributionEntry("assessed molecular reference", representation.SpeciesId,
                representation.ForceFieldVersion,
                "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                representation.CoordinateTemplateSha256,
                "OpenMM application package terms; no separate coordinate template terms are asserted",
                OpenMmTermsUri,
                representation.ForceFieldFamily,
                "One-molecule assessment template identified by digest. Construction used the native full membrane patch instead; this reference template was not inserted into the resulting system or bundled.",
                "used as assessment reference, not incorporated"));
        foreach (var asset in attempt.ForceFieldFiles)
            incorporated.Add(new AttributionEntry("force-field parameter contribution", asset.Id,
                asset.Version, "https://ambermd.org/AmberModels.php", asset.Sha256,
                "Amber developers state their force fields are public domain",
                "https://ambermd.org/index.php", ForceFieldCitation(asset.Family),
                "System XML includes derived parameter values; the source parameter file is not bundled.",
                "parameters incorporated into System XML"));
        var tools = new List<AttributionEntry>
        {
            new("construction tool", policy.Construction.ProviderName,
                attempt.ConstructionProviderVersion,
                "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                null,
                "OpenMM application layer MIT terms", OpenMmTermsUri,
                "OpenMM: Eastman et al., PLoS Computational Biology 2017, 10.1371/journal.pcbi.1005659",
                "The tool executable is not bundled and its code digest was not recorded. Its exact full build is stated; the native patch has its own digest under incorporated data.",
                "used to produce, not incorporated as code")
        };
        if (!string.IsNullOrWhiteSpace(ppmVersion) && !string.IsNullOrWhiteSpace(ppmExecutableSha256))
            tools.Add(new AttributionEntry("orientation tool", "PPM 2.0", ppmVersion,
                "https://cggit.cc.lehigh.edu/biomembhub/ppm2_server_code/-/tree/0de704acd10fbf7f4fbafb9bc925e01b06f38ceb",
                ppmExecutableSha256,
                "Pinned-source license not independently confirmed; code terms do not automatically apply to the derived orientation output",
                "https://cggit.cc.lehigh.edu/biomembhub/ppm2_server_code/-/blob/0de704acd10fbf7f4fbafb9bc925e01b06f38ceb/license.txt",
                "Lomize et al., Nucleic Acids Research 2012, 10.1093/nar/gkr703; PPM method 10.1021/ci200020k",
                "The pinned license URL was not independently retrievable. The PPM executable is not bundled; the oriented result is identified separately by digest.",
                "used to produce, not incorporated as code"));
        return new { incorporatedData = incorporated, referenceInputs = references, toolsUsed = tools };
    }

    private static string ForceFieldCitation(string family) => family switch
    {
        "Amber19-ff19SB" => "Tian et al., Journal of Chemical Theory and Computation 2020, 10.1021/acs.jctc.9b00591",
        "Lipid21" => "Dickson, Walker and Gould, Journal of Chemical Theory and Computation 2022, 10.1021/acs.jctc.1c01217",
        "Amber19-TIP3P-JC" => "Jorgensen et al., Journal of Chemical Physics 1983, 10.1063/1.445869; Joung and Cheatham, Journal of Physical Chemistry B 2008, 10.1021/jp8001614",
        _ => family
    };

    private static JsonSerializerOptions ManifestJson()
    {
        var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        options.Converters.Add(new JsonStringEnumConverter());
        return options;
    }

    private static string? SourceAccession(AssessedPreparedProtein protein,
        ApplicablePreparationPolicy policy)
    {
        return IsQualified6Qwr(protein, policy) ? "6QWR" : protein.Intended.Source.Accession;
    }

    private static bool IsQualified6Qwr(AssessedPreparedProtein protein,
        ApplicablePreparationPolicy policy) =>
        protein.Intended.Source.Sha256.Equals(Qualified6QwrSourceSha256,
            StringComparison.OrdinalIgnoreCase) &&
        policy.Scope?.SourceCoordinateSha256.Equals(Qualified6QwrSourceSha256,
            StringComparison.OrdinalIgnoreCase) == true &&
        protein.Intended.ModelIndex == 0;

    private sealed record BundleEntry(string Name, string Path, string Sha256);
    private sealed record FinalCell(double[][] BoxVectorsAngstrom,
        double[] LengthsAngstrom, double[] AnglesDegrees);
    private sealed record AttributionEntry(string Role, string Name, string Version,
        string SourceUri, string? Sha256, string Rights, string RightsUrl, string Citation,
        string UseLimitations, string Relationship);
}
