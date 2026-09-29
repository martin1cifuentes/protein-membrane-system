using System.Collections.Immutable;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport;

/// <summary>Delivers one already completed stage with its own current assessment.</summary>
public sealed class CompletedStageExport
{
    private const string OpenMmTermsUri =
        "https://docs.openmm.org/latest/userguide/library/01_introduction.html#license";
    private readonly ICompletedStageExportWork _worker;

    public CompletedStageExport(ICompletedStageExportWork worker) => _worker = worker;

    private sealed record FinalAtomIdentity(string Name, string Element);

    private static string? ExactApprovedHeavyAtomName(AtomCorrespondence atom,
        AssessedPreparedProtein protein, FinalAtomIdentity? finalAtom)
    {
        if (atom.ApprovedChangeId is null ||
            atom.Role != AtomOriginKind.Generated || atom.MoleculeRole != MoleculeRoleKind.Protein ||
            atom.SourceAtomId is not null || atom.SourceResidue is null || finalAtom is null ||
            atom.Element.Equals("H", StringComparison.OrdinalIgnoreCase) ||
            atom.Element.Equals("D", StringComparison.OrdinalIgnoreCase) ||
            finalAtom.Element != atom.Element ||
            string.IsNullOrWhiteSpace(finalAtom.Name) ||
            !protein.Correspondence.Complete ||
            protein.Correspondence.ResultId != protein.Molecule.CoordinateSha256 ||
            protein.Correspondence.Atoms.IsDefaultOrEmpty)
            return null;
        // Construction checks the retained atom multiset; Memgen may reorder
        // it. Exact source provenance plus the hash-verified final atom name
        // identifies one prepared atom without relying on display labels.
        var matching = protein.Correspondence.Atoms.Where(prepared =>
            prepared.ApprovedChangeId == atom.ApprovedChangeId &&
            prepared.Role == atom.Role && prepared.MoleculeRole == atom.MoleculeRole &&
            prepared.AtomRole == atom.AtomRole &&
            prepared.SourceAtomId == atom.SourceAtomId &&
            prepared.SourceResidue == atom.SourceResidue &&
            prepared.Element == atom.Element && prepared.ResultAtomId.EndsWith(
                ":" + finalAtom.Name, StringComparison.Ordinal)).ToArray();
        return matching.Length == 1 ? finalAtom.Name : null;
    }

    private static bool AuthorizedProteinChange(AtomCorrespondence atom,
        AssessedPreparedProtein protein, ImmutableArray<ResearcherDecision> decisions,
        FinalAtomIdentity? finalAtom)
    {
        if (atom.ApprovedChangeId is null)
            return atom.Role != AtomOriginKind.Generated ||
                atom.MoleculeRole != MoleculeRoleKind.Protein ||
                atom.Element.Equals("H", StringComparison.OrdinalIgnoreCase) ||
                atom.Element.Equals("D", StringComparison.OrdinalIgnoreCase);
        var atomName = ExactApprovedHeavyAtomName(atom, protein, finalAtom);
        if (atomName is null) return false;
        if (atom.ApprovedChangeId != "recommendation-plan-only")
            return decisions.Any(decision => decision.Id == atom.ApprovedChangeId &&
                decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                decision.ChosenValue == ResearcherDecisionValue.Approved &&
                protein.Changes.Any(change => change.Id == decision.SubjectId &&
                    change.Kind == PreparationChangeKind.HeavyAtom &&
                    change.Residue == atom.SourceResidue &&
                    change.ProposedChange == atomName &&
                    change.StudyRevisionId == decision.StudyRevisionId &&
                    change.StudyRevisionId == protein.PreparationStudyRevisionId &&
                    change.IntendedProteinId == protein.Intended.Id));

        var plan = protein.RecommendationPlan;
        if (plan is null || plan.PlanSha256.Length != 64 || !plan.PlanSha256.All(Uri.IsHexDigit) ||
            !plan.CheckedCandidateSha256.Equals(protein.Molecule.CoordinateSha256,
                StringComparison.OrdinalIgnoreCase) ||
            string.IsNullOrWhiteSpace(protein.PreparationStudyRevisionId) ||
            atom.SourceResidue is null)
            return false;

        // A checked plan's candidate uses this marker for a proposed heavy atom.
        // Bind it back to the exact prepared atom and the approved proposal;
        // the marker alone is never an authorization.
        if (plan.ProposedHeavyAtoms.IsDefaultOrEmpty ||
            plan.ProposedHeavyAtoms.Count(item => item.Residue == atom.SourceResidue &&
                item.AtomName == atomName) != 1)
            return false;
        var proposed = protein.Changes.Where(change =>
            change.Kind == PreparationChangeKind.HeavyAtom &&
            change.Residue == atom.SourceResidue && change.ProposedChange == atomName &&
            change.StudyRevisionId == protein.PreparationStudyRevisionId &&
            change.IntendedProteinId == protein.Intended.Id).ToArray();
        return proposed.Length == 1 && decisions.Count(decision =>
            decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
            decision.ChosenValue == ResearcherDecisionValue.Approved &&
            decision.SubjectId == proposed[0].Id &&
            decision.StudyRevisionId == proposed[0].StudyRevisionId) == 1;
    }

    private static async Task<Dictionary<int, FinalAtomIdentity>> ReadApprovedHeavyAtomsAsync(
        MolecularArtifact molecule, SourceToResultCorrespondence correspondence,
        CancellationToken cancellationToken)
    {
        var indices = correspondence.Atoms.Where(item =>
            item.ApprovedChangeId is not null && item.Role == AtomOriginKind.Generated &&
            item.MoleculeRole == MoleculeRoleKind.Protein &&
            !item.Element.Equals("H", StringComparison.OrdinalIgnoreCase) &&
            !item.Element.Equals("D", StringComparison.OrdinalIgnoreCase))
            .Select(item => item.ResultAtomIndex).ToArray();
        var result = new Dictionary<int, FinalAtomIdentity>();
        if (indices.Length == 0) return result;
        var bytes = await File.ReadAllBytesAsync(molecule.TopologyPath!, cancellationToken);
        if (!Convert.ToHexString(SHA256.HashData(bytes)).Equals(molecule.TopologySha256,
            StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("The stage topology changed from its recorded digest.");
        using var document = JsonDocument.Parse(bytes);
        var root = document.RootElement;
        if (!root.TryGetProperty("atoms", out var atoms) || atoms.ValueKind != JsonValueKind.Array ||
            atoms.GetArrayLength() != correspondence.Atoms.Length)
            throw new InvalidDataException("The stage topology lacks its exact atom index.");
        foreach (var index in indices)
        {
            var item = atoms[index];
            if (item.ValueKind != JsonValueKind.Object ||
                !item.TryGetProperty("name", out var name) || name.ValueKind != JsonValueKind.String ||
                !item.TryGetProperty("element", out var element) ||
                element.ValueKind != JsonValueKind.String)
                throw new InvalidDataException("An approved heavy atom has incomplete final topology identity.");
            result[index] = new FinalAtomIdentity(name.GetString() ?? "", element.GetString() ?? "");
        }
        return result;
    }

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
            stage.Kind == StageKind.Minimization && stage.Observation.Termination != StageTermination.Converged ||
            assessment.StageId != stage.Id || !assessment.CurrentlyApplicable ||
            stage.Attempt.StudyRevisionId != originatingRevision.Id || stage.Observation.StageId != stage.Id ||
            stage.Observation.AttemptId != stage.Attempt.Id ||
            stage.PolicyId != policy.Id ||
            !PreparationPolicyFingerprint.Matches(stage.Attempt, policy) ||
            !PreparationPolicyFingerprint.Matches(constructed.Attempt, policy) ||
            stage.Attempt.ProteinId != protein.Id || protein.StudyRevisionId != originatingRevision.Id ||
            string.IsNullOrWhiteSpace(protein.PreparationStudyRevisionId) ||
            protein.Intended.Id != originatingRevision.IntendedProtein?.Id ||
            stage.Attempt.MembraneId != membrane.Id || membrane.StudyRevisionId != originatingRevision.Id ||
            stage.Attempt.PlacementId != placement.Id || placement.StudyRevisionId != originatingRevision.Id ||
            placement.Proposal.PreparedProteinId != protein.Id ||
            placement.Proposal.MembraneModelId != membrane.Intended.Id ||
            constructed.Attempt.Id != stage.Attempt.Id ||
            stage.Kind == StageKind.Minimization &&
                stage.Molecule.TopologySha256 != constructed.Molecule.TopologySha256 ||
            derivation.AttemptId != stage.Attempt.Id || constructed.Derivation != derivation ||
            !double.IsFinite(policy.ExportCoordinateReadBackToleranceAngstrom) ||
            policy.ExportCoordinateReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellLengthReadBackToleranceAngstrom) ||
            policy.ExportCellLengthReadBackToleranceAngstrom < 0 ||
            !double.IsFinite(policy.ExportCellAngleReadBackToleranceDegrees) ||
            policy.ExportCellAngleReadBackToleranceDegrees < 0 ||
            stage.Correspondence.Atoms.IsDefaultOrEmpty || !stage.Correspondence.Complete ||
            stage.Molecule.Id != stage.Id ||
            stage.Correspondence.ResultId != stage.Molecule.Id ||
            stage.Correspondence.Atoms.Length != stage.Molecule.AtomCount ||
            stage.Correspondence.Atoms.Where((atom, index) => atom.ResultAtomIndex != index ||
                string.IsNullOrWhiteSpace(atom.ResultAtomId)).Any() ||
            stage.Correspondence.Atoms.Select(atom => atom.ResultAtomId)
                .Distinct(StringComparer.Ordinal).Count() != stage.Molecule.AtomCount ||
            stage.Correspondence.SourceId != constructed.Correspondence.SourceId ||
            constructed.Correspondence.Atoms.IsDefault ||
            !stage.Correspondence.Atoms.SequenceEqual(constructed.Correspondence.Atoms) ||
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

        if (policy.Construction.Route == ConstructionRouteKind.PackmolMemgen &&
            !await VerifiedMemgenConstructionProvenanceAsync(constructed, derivation,
                cancellationToken))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable(
                "The selected Memgen construction's checked Amber artifacts or retained atom correspondence are unavailable or changed.");

        var molecule = stage.Molecule;
        if (string.IsNullOrWhiteSpace(molecule.TopologyPath) || string.IsNullOrWhiteSpace(molecule.SystemXmlPath) ||
            string.IsNullOrWhiteSpace(molecule.StateXmlPath) || !File.Exists(molecule.TopologyPath) ||
            !File.Exists(molecule.SystemXmlPath) || !File.Exists(molecule.StateXmlPath))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable("The completed stage's corresponding topology, System and State are unavailable.");

        Dictionary<int, FinalAtomIdentity> approvedHeavyAtoms;
        try
        {
            approvedHeavyAtoms = await ReadApprovedHeavyAtomsAsync(molecule, stage.Correspondence,
                cancellationToken);
        }
        catch (Exception exception) when (exception is IOException or JsonException or
                                         InvalidDataException or UnauthorizedAccessException)
        {
            return BoundaryOutcome<CompletedStageBundle>.Unavailable(
                $"The completed stage's exact approved heavy atom identity is unavailable: {exception.Message}");
        }
        if (stage.Correspondence.Atoms.Any(atom => !AuthorizedProteinChange(atom, protein,
                decisions, approvedHeavyAtoms.GetValueOrDefault(atom.ResultAtomIndex))))
            return BoundaryOutcome<CompletedStageBundle>.Unavailable(
                "A generated protein atom does not correspond to an exact approved preparation change or checked recommendation plan.");

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
                            accession = protein.Intended.Source.Accession,
                            protein.Intended.Source.Sha256,
                            protein.Intended.Source.SourceModelDescription,
                            protein.Intended.Source.UploadProvenance,
                            protein.Intended.Source.UploadProvenanceNote },
                        protein.Intended.ModelIndex,
                        protein.Intended.SourceModelId,
                        sourceModelNumber = SourceModelNumber(protein.Intended),
                        protein.Intended.BiologicalAssemblyId,
                        protein.Intended.Chains, protein.Intended.Partners,
                        protein.Intended.AlternateLocations,
                        preparedCoordinateSha256 = protein.Molecule.CoordinateSha256,
                        preparedTopologySha256 = protein.Molecule.TopologySha256,
                        protein.ChemicalStatePolicyId, protein.ChemicalStatePolicyVersion,
                        protein.StructuralAssessmentPolicyId, protein.StructuralAssessmentPolicyVersion,
                        protein.ResidueVariants, protein.Changes, protein.Limitations,
                        protein.RecommendationPlan,
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
                        constructed.MaximumProteinCoordinateDeviationAngstrom,
                        constructed.Evidence, constructed.Findings,
                        constructed.LocalState, derivation,
                        constructedCoordinateSha256 = constructed.Molecule.CoordinateSha256,
                        constructedTopologySha256 = constructed.Molecule.TopologySha256,
                        provider = new { policy.Construction.ProviderName,
                            policy.Construction.ProviderVersion,
                            policy.Construction.Route, policy.Construction.SaltConvention,
                            assets = (policy.Construction.ProviderAssets.IsDefault
                                ? ImmutableArray<ProviderAsset>.Empty
                                : policy.Construction.ProviderAssets).Select(asset => new
                            {
                                asset.Id, asset.Version, asset.Sha256
                            }).ToArray(),
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
                        sourceAccession = protein.Intended.Source.Accession,
                        sourceCoordinateSha256 = protein.Intended.Source.Sha256,
                        sourceModelIndex = protein.Intended.ModelIndex,
                        sourceModelId = protein.Intended.SourceModelId,
                        sourceModelNumber = SourceModelNumber(protein.Intended),
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

    private static int SourceModelNumber(IntendedProteinModel intended) =>
        int.TryParse(intended.SourceModelId, System.Globalization.NumberStyles.None,
            System.Globalization.CultureInfo.InvariantCulture, out var observed) && observed > 0
            ? observed : intended.ModelIndex + 1;

    private static async Task<bool> VerifiedMemgenConstructionProvenanceAsync(
        ConstructedExplicitSystem constructed, ConstructionDerivation derivation,
        CancellationToken cancellationToken)
    {
        var molecule = constructed.Molecule;
        if (derivation.Trials.IsDefaultOrEmpty ||
            string.IsNullOrWhiteSpace(derivation.SelectedTrialId) ||
            derivation.AmberImport is not { } import ||
            import.AmberAtomCount != molecule.AtomCount ||
            import.ImportedAtomCount != molecule.AtomCount ||
            !import.AtomOrderPreserved || !import.BondedTermsPreserved ||
            !import.NonbondedTermsPreserved || !import.ExclusionsPreserved ||
            !import.UnitsPreserved || !import.ParameterComparisonPerformed ||
            import.ParameterCorrespondenceSha256 is not { Length: 64 } parameterSha ||
            !parameterSha.All(Uri.IsHexDigit) ||
            string.IsNullOrWhiteSpace(molecule.CorrespondencePath) ||
            molecule.CorrespondenceSha256 is not { Length: 64 } mappingSha ||
            !mappingSha.All(Uri.IsHexDigit) || !constructed.Correspondence.Complete)
            return false;
        var selected = derivation.Trials.Where(trial =>
            trial.TrialId == derivation.SelectedTrialId).Take(2).ToArray();
        if (selected.Length != 1 || selected[0].Standing != ConstructionTrialStanding.Checked ||
            selected[0].DiagnosticArtifacts.IsDefaultOrEmpty)
            return false;

        foreach (var role in new[] { "amberTopology", "amberFinalRestart" })
        {
            var artifacts = selected[0].DiagnosticArtifacts.Where(artifact =>
                artifact.Role == role).Take(2).ToArray();
            if (artifacts.Length != 1 || artifacts[0].LocalPath is not { Length: > 0 } path ||
                artifacts[0].Sha256.Length != 64 || !artifacts[0].Sha256.All(Uri.IsHexDigit) ||
                !artifacts[0].Sha256.Equals(role == "amberTopology" ? import.PrmtopSha256 :
                    import.FinalRestartSha256, StringComparison.OrdinalIgnoreCase) ||
                !await FileMatchesAsync(path, artifacts[0].Sha256,
                    cancellationToken))
                return false;
        }

        try
        {
            var bytes = await File.ReadAllBytesAsync(molecule.CorrespondencePath,
                cancellationToken);
            if (!Convert.ToHexString(SHA256.HashData(bytes)).Equals(mappingSha,
                    StringComparison.OrdinalIgnoreCase))
                return false;
            var mapped = JsonSerializer.Deserialize<SourceToResultCorrespondence>(bytes,
                new JsonSerializerOptions(JsonSerializerDefaults.Web));
            return mapped is { Complete: true } &&
                mapped.SourceId == constructed.Correspondence.SourceId &&
                mapped.ResultId == constructed.Correspondence.ResultId &&
                mapped.ResultId == molecule.CoordinateSha256 &&
                !mapped.Atoms.IsDefault &&
                mapped.Atoms.SequenceEqual(constructed.Correspondence.Atoms);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException or
                                         JsonException)
        {
            return false;
        }
    }

    private static async Task<bool> FileMatchesAsync(string path, string expectedSha256,
        CancellationToken cancellationToken)
    {
        try
        {
            await using var stream = File.OpenRead(path);
            var actual = Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken));
            return actual.Equals(expectedSha256, StringComparison.OrdinalIgnoreCase);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        {
            return false;
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
        var rcsb = source.Kind == SourceRouteKind.Rcsb && !string.IsNullOrWhiteSpace(source.Accession);
        var sourceAccession = source.Accession;
        var incorporated = new List<AttributionEntry>
        {
            new("protein source coordinates", sourceAccession ?? source.Id,
                source.SourceModelDescription ?? "identified source model",
                rcsb ? $"https://www.rcsb.org/structure/{sourceAccession}" : source.Provenance,
                source.Sha256,
                rcsb ? "RCSB PDB archive data: CC0" : "Source-specific rights were not established by this record",
                rcsb ? "https://www.rcsb.org/pages/usage-policy" : source.Provenance,
                sourceAccession ?? source.Provenance,
                "Derived coordinates are included; the original source file is identified by digest but is not bundled.",
                "incorporated and transformed"),
        };
        if (policy.Construction.Route == ConstructionRouteKind.PackmolMemgen)
        {
            foreach (var asset in (policy.Construction.ProviderAssets.IsDefault
                         ? ImmutableArray<ProviderAsset>.Empty
                         : policy.Construction.ProviderAssets).Where(asset =>
                asset.Id.EndsWith("pdbs.tar.gz", StringComparison.Ordinal) ||
                asset.Id.EndsWith("memgen.parm", StringComparison.Ordinal) ||
                asset.Id.Contains("leaprc.", StringComparison.Ordinal)))
                incorporated.Add(new AttributionEntry("provider molecular input asset",
                    asset.Id, asset.Version, "https://ambermd.org/",
                    asset.Sha256,
                    "The installed AmberTools asset's redistribution terms are not asserted by this record",
                    "https://ambermd.org/", "PACKMOL-Memgen 2026.3.25; AmberTools26",
                    "The identified provider asset supplied molecular templates, population data or parameter selection; its original file is not bundled.",
                    "incorporated through the constructed system"));
        }
        else if (policy.Construction.NativePatchMode == "mapped-lipid21-zenodo-popc")
        {
            const string record = "https://doi.org/10.5281/zenodo.14776136";
            const string rights = "https://creativecommons.org/licenses/by/4.0/";
            const string citation = "Frankel, Amber Lipid21 bilayer simulations, Zenodo v1 (2025), 10.5281/zenodo.14776136";
            incorporated.Add(new AttributionEntry("native lipid coordinate patch",
                "POPC one-residue Lipid21 conversion", "Zenodo v1 mapped for OpenMM 8.6",
                record, attempt.NativePatchSha256, "Creative Commons Attribution 4.0; adapted coordinates",
                rights, citation,
                "Exact atom permutation, bonds, coordinates and cell are identified by the mapped patch and source digests; neither source file is bundled in this export.",
                "incorporated and transformed"));
            incorporated.Add(new AttributionEntry("source bilayer coordinates", "POPC.gro",
                "Zenodo v1", record, policy.Construction.NativeSourcePatchSha256,
                "Creative Commons Attribution 4.0", rights, citation,
                "The original solvated periodic POPC coordinates supplied the mapped starting patch; the source GRO is identified by digest and is not bundled.",
                "incorporated and transformed"));
        }
        else
            incorporated.Add(new AttributionEntry("native lipid coordinate patch",
                policy.Construction.LipidTypeArgument + " patch",
                policy.Construction.ProviderVersion,
                "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                attempt.NativePatchSha256,
                "OpenMM application layer MIT terms; no separate patch-file terms are asserted",
                OpenMmTermsUri,
                "OpenMM application package, exact build identified by attempt construction provider version",
                "Lipid coordinates are derived into the bundle; the native source patch is identified, not bundled.",
                "incorporated and transformed"));
        var references = new List<AttributionEntry>();
        foreach (var representation in membrane.SpeciesRepresentations)
            references.Add(new AttributionEntry("assessed molecular reference", representation.SpeciesId,
                representation.ForceFieldVersion,
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "https://ambermd.org/" :
                    "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                representation.CoordinateTemplateSha256,
                "The identified asset's redistribution terms are not asserted by this record",
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "https://ambermd.org/" : OpenMmTermsUri,
                representation.ForceFieldFamily,
                "One-molecule assessment template identified by digest; its role in construction is reported by the selected recipe and provider-asset account.",
                "used as assessment reference"));
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
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "https://ambermd.org/" :
                    "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
                null,
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "Installed AmberTools terms; executable is not bundled" :
                    "OpenMM application layer MIT terms",
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "https://ambermd.org/" : OpenMmTermsUri,
                policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ?
                    "PACKMOL-Memgen, Journal of Chemical Information and Modeling 2019, 10.1021/acs.jcim.9b00269" :
                    "OpenMM: Eastman et al., PLoS Computational Biology 2017, 10.1371/journal.pcbi.1005659",
                "The tool executable is not bundled. The exact provider version and available asset digests are recorded with the construction account.",
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

    private sealed record BundleEntry(string Name, string Path, string Sha256);
    private sealed record FinalCell(double[][] BoxVectorsAngstrom,
        double[] LengthsAngstrom, double[] AnglesDegrees);
    private sealed record AttributionEntry(string Role, string Name, string Version,
        string SourceUri, string? Sha256, string Rights, string RightsUrl, string Citation,
        string UseLimitations, string Relationship);
}
