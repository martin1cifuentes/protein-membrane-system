using System.Collections.Immutable;
using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Serialization;
using ProteinPreparationBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinPreparation.ProteinPreparation;
using MembraneAssessmentBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.MembraneModelAssessment.MembraneModelAssessment;
using PlacementAssessmentBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.PlacementAssessment.PlacementAssessment;
using ExplicitPreparationBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation;
using PreparationAssessmentBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment.PreparationAssessment;
using InspectionBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.ConnectedStructuralInspection.ConnectedStructuralInspection;
using ExportBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport.CompletedStageExport;
using WorkspaceBoundary = ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem;

/// <summary>The installed package resource observed through the selected worker interpreter.</summary>
public sealed record ConstructionProviderInstallation(string FullVersion, string NativePatchPath,
    string NativePatchSha256, ImmutableArray<NativePatchInstallation> AdditionalPatches = default,
    string? MemgenVersion = null, string? MemgenHome = null);

/// <summary>An additional exact built-in patch from the same selected OpenMM installation.</summary>
public sealed record NativePatchInstallation(string SpeciesId, string Path, string Sha256);

/// <summary>
/// Owns the study, decisions, and handoffs among direct children. Provider
/// observations are never promoted here merely because a process succeeded.
/// </summary>
public sealed class ProteinInMembraneSystem
{
    private sealed record ProteinSelectionIssue(string Message, bool ReviewMoleculeSelection = false);

    private readonly object _gate = new();
    private readonly SemaphoreSlim _commands = new(1, 1);
    private readonly ExternalSourceExchange _sources;
    private readonly ProteinPreparationBoundary _proteinPreparation;
    private readonly MembraneAssessmentBoundary _membraneAssessment;
    private readonly PlacementAssessmentBoundary _placementAssessment;
    private readonly ExplicitPreparationBoundary _explicitPreparation;
    private readonly PreparationAssessmentBoundary _preparationAssessment = new();
    private readonly InspectionBoundary _inspection = new();
    private readonly ExportBoundary _export;
    private readonly WorkspaceBoundary _workspace = new();
    private readonly string _workspaceRoot;
    private readonly string _ppmExecutablePath;
    private readonly Func<ConstructionProviderInstallation?> _constructionProviderProbe;
    private readonly ConstructionProviderInstallation? _startupConstructionProvider;
    private ConstructionProviderInstallation? _lastObservedConstructionProvider;
    private readonly CatalogueDocument? _catalogue;
    private readonly string? _catalogueIssue;
    private readonly Dictionary<string, StructureBinding> _structures = new(StringComparer.Ordinal);
    private readonly Dictionary<string, DiagnosticBinding> _diagnostics = new(StringComparer.Ordinal);
    private readonly Dictionary<string, StructuralSource> _uploads = new(StringComparer.Ordinal);
    private readonly Dictionary<(string SourceId, int ModelIndex, string? AssemblyId), SourcePreviewAccount> _sourcePreviews = new();
    private readonly Dictionary<string, StudyRevision> _revisions = new(StringComparer.Ordinal);
    private readonly Dictionary<string, CompletedStage> _stages = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ConstructedExplicitSystem> _constructedByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ConstructionDerivation> _derivationByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ImmutableArray<ConstructionTrialSummary>> _trialsByAttempt = new(StringComparer.Ordinal);
    private readonly List<AttemptAccount> _priorAttempts = new();
    private readonly Dictionary<string, ApplicablePreparationPolicy> _policyByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedPreparedProtein> _proteinByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedMembraneModel> _membraneByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedProteinMembranePlacement> _placementByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ImmutableArray<ResearcherDecision>> _decisionsByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, PreparationAssessmentResult> _stageAssessments = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ImmutableArray<ScientificFinding>> _laterFindingsByStage = new(StringComparer.Ordinal);
    private readonly Dictionary<string, CompletedStageBundle> _bundles = new(StringComparer.Ordinal);
    private readonly Dictionary<string, StageExportAccount> _exportAccounts = new(StringComparer.Ordinal);
    private ImmutableArray<CandidateSourceRecord> _candidates = ImmutableArray<CandidateSourceRecord>.Empty;
    private ImmutableArray<string> _sourceSearchIssues = ImmutableArray<string>.Empty;
    private StudyRevision _study;
    private StructuralSource? _selectedSource;
    private SourceInspectionReport? _sourceInspection;
    private PreparationProposalReport? _preparationProposals;
    private RecommendedPreparationPlan? _preparationPlan;
    private RecommendedPreparationPlan? _authorizedPreparationPlan;
    private ImmutableArray<ResidueVariantChoice> _planOverrides = ImmutableArray<ResidueVariantChoice>.Empty;
    private bool _recommendationRunning;
    private long _recommendationGeneration;
    private string? _recommendationIssue;
    private ImmutableArray<ResearcherDecision> _decisions = ImmutableArray<ResearcherDecision>.Empty;
    private ImmutableArray<ResearcherDecision> _allDecisions = ImmutableArray<ResearcherDecision>.Empty;
    private AssessedPreparedProtein? _protein;
    private ProteinPreparationDiagnostic? _proteinDiagnostic;
    private bool _proteinSelectionRunning;
    private ProteinSelectionIssue? _proteinSelectionIssue;
    private bool _proteinPreparationRunning;
    private string? _proteinPreparationFailure;
    private bool _proteinPreparationRetryable;
    private string? _proteinPreparationObservationIssue;
    private AssessedMembraneModel? _membrane;
    private string? _membraneAssessmentReason;
    private bool _membraneAssessmentRunning;
    private bool _membraneAssessmentUnavailable;
    private PlacementProposal? _placementProposal;
    private OpmReferenceReview? _opmReview;
    private PlacementMeasurementReport? _placementMeasurement;
    private string? _placementMeasurementIssue;
    private PlacementStructuralWitness? _placementWitness;
    private AssessedProteinMembranePlacement? _placement;
    private bool _placementRunning;
    private bool _placementAssessing;
    private string? _placementOperationIssue;
    private PlacementRouteOutcomeAccount? _placementRouteOutcome;
    private PreparationAttempt? _currentAttempt;
    private StageExecutionState? _execution;
    private CancellationTokenSource? _attemptStop;
    private string? _stopRequestedAttemptId;
    private Task? _attemptTask;
    private long _revision;

    public event Action? Changed;

    public ProteinInMembraneSystem(
        IScientificWorkerExchange worker,
        ExternalSourceExchange sources,
        string workspaceRoot,
        string policyCataloguePath,
        string ppmExecutablePath,
        Func<ConstructionProviderInstallation?> constructionProviderProbe)
    {
        _sources = sources;
        _proteinPreparation = new ProteinPreparationBoundary(worker);
        _membraneAssessment = new MembraneAssessmentBoundary(worker);
        _placementAssessment = new PlacementAssessmentBoundary(worker);
        _explicitPreparation = new ExplicitPreparationBoundary(worker, worker, worker);
        if (worker is ScientificWorkerExchange localExchange)
            localExchange.ProgressObserved += ObserveConstructionProgress;
        _export = new ExportBoundary(worker);
        _workspaceRoot = Path.GetFullPath(workspaceRoot);
        Directory.CreateDirectory(_workspaceRoot);
        _ppmExecutablePath = ppmExecutablePath;
        _constructionProviderProbe = constructionProviderProbe;
        (_catalogue, _catalogueIssue) = ReadCatalogue(policyCataloguePath);
        _startupConstructionProvider = _catalogue is { PreparationPolicies: { IsEmpty: false } }
            ? _constructionProviderProbe() : null;
        _lastObservedConstructionProvider = _startupConstructionProvider;
        _study = new StudyRevision(NewId(), 1, null, null, null, FixedStudyConditions.Initial);
        _revisions.Add(_study.Id, _study);
        _workspace.RetainStudy(_study);
    }

    public WorkspaceState Snapshot()
    {
        lock (_gate) return SnapshotLocked();
    }

    /// <summary>Materializes an exact read-only source model or assembly for visual inspection.</summary>
    public async Task<BoundaryOutcome<SourcePreviewAccount>> PreviewSourceModelAsync(
        string sourceId, int modelIndex, string? assemblyId, CancellationToken cancellationToken)
    {
        SourceInspectionReport inspection;
        var key = (sourceId, modelIndex, assemblyId);
        lock (_gate)
        {
            if (_sourceInspection is null || _selectedSource?.Id != sourceId)
                return BoundaryOutcome<SourcePreviewAccount>.Unavailable("The requested source is no longer selected.");
            if (_sourcePreviews.TryGetValue(key, out var cached)) return BoundaryOutcome<SourcePreviewAccount>.Success(cached);
            inspection = _sourceInspection;
        }
        var prepared = await _proteinPreparation.PreviewSourceModelAsync(inspection, modelIndex, assemblyId,
            WorkDirectory("source-preview", NewId()),
            _catalogue?.MaximumSourceAtoms is > 0 ? _catalogue.MaximumSourceAtoms : int.MaxValue,
            cancellationToken);
        if (prepared.Value is null) return BoundaryOutcome<SourcePreviewAccount>.Unavailable(prepared.Reason);
        lock (_gate)
        {
            if (_sourceInspection != inspection || _selectedSource?.Id != sourceId)
                return BoundaryOutcome<SourcePreviewAccount>.Unavailable("The requested source changed while its view was being prepared.");
            if (_sourcePreviews.TryGetValue(key, out var cached)) return BoundaryOutcome<SourcePreviewAccount>.Success(cached);
            var model = inspection.Models.Single(item => item.Index == modelIndex);
            var chains = assemblyId is null ? model.Chains.Select(item => item.Name).ToImmutableArray() :
                model.Assemblies.Single(item => item.Name == assemblyId).ChainCopies.Select(item => item.CopyId).ToImmutableArray();
            var account = new SourcePreviewAccount(sourceId, modelIndex, assemblyId,
                StructureUrlLocked(prepared.Value.Path), chains);
            _sourcePreviews[key] = account;
            return BoundaryOutcome<SourcePreviewAccount>.Success(account);
        }
    }

    public async Task<BoundaryOutcome<WorkspaceState>> ExecuteAsync(ActorCommand command, CancellationToken cancellationToken)
    {
        await _commands.WaitAsync(cancellationToken);
        try
        {
            lock (_gate)
                if (command.ExpectedRevision != _revision)
                    return BoundaryOutcome<WorkspaceState>.Unavailable("The workspace changed; review its current account before making this decision.");
            string? refusal = command.Kind switch
            {
                ActorActionKind.SearchSource => await SearchSourceAsync(command.Data, cancellationToken),
                ActorActionKind.SelectSource => await SelectSourceAsync(command.Data, cancellationToken),
                ActorActionKind.SelectProteinModel => await SelectProteinModelAsync(command.Data, cancellationToken),
                ActorActionKind.ApprovePreparationChange => await DecidePreparationChangeAsync(command.Data, cancellationToken),
                ActorActionKind.AuthorizePreparationPlan => await AuthorizePreparationPlanAsync(command.Data, cancellationToken),
                ActorActionKind.OverridePreparationPlanChoice => await OverridePreparationPlanChoiceAsync(command.Data, cancellationToken),
                ActorActionKind.RetryPreparationPlan => await RetryPreparationPlanAsync(cancellationToken),
                ActorActionKind.StartProteinPreparation => await StartProteinPreparationAsync(cancellationToken),
                ActorActionKind.RetryProteinPreparation => await RetryProteinPreparationAsync(cancellationToken),
                ActorActionKind.AdoptMembrane => await AdoptMembraneAsync(command.Data, cancellationToken),
                ActorActionKind.RetryMembraneCheck => await RetryMembraneCheckAsync(command.Data, cancellationToken),
                ActorActionKind.ProposePlacement => await ProposePlacementAsync(command.Data, cancellationToken),
                ActorActionKind.RevisePlacement => await RevisePlacementAsync(command.Data, cancellationToken),
                ActorActionKind.AdoptPlacement => AdoptPlacement(command.Data),
                ActorActionKind.BuildAndMinimize => StartPreparation(command.Data),
                ActorActionKind.StopAttempt => StopAttempt(command.Data),
                ActorActionKind.RequestEquilibration => RequestEquilibration(command.Data),
                ActorActionKind.SelectInspectionSubject => SelectInspectionSubject(command.Data),
                ActorActionKind.SetInspectionFocus => SetInspectionFocus(command.Data),
                ActorActionKind.ExportStage => await ExportStageAsync(command.Data, cancellationToken),
                _ => "This actor action has no established route."
            };
            if (refusal is not null) return BoundaryOutcome<WorkspaceState>.Unavailable(refusal);
            lock (_gate) return BoundaryOutcome<WorkspaceState>.Success(SnapshotLocked());
        }
        catch (JsonException)
        {
            return BoundaryOutcome<WorkspaceState>.Unavailable("The action's information could not be interpreted.");
        }
        catch (InvalidOperationException exception)
        {
            return BoundaryOutcome<WorkspaceState>.Unavailable(exception.Message);
        }
        finally
        {
            _commands.Release();
        }
    }

    public async Task<string> UploadAsync(Stream input, string fileName, UploadOriginKind provenance,
        string? provenanceNote, CancellationToken cancellationToken)
    {
        if (!Enum.IsDefined(provenance))
            throw new InvalidDataException("Declare whether this uploaded structure is predicted, experimental, or unknown.");
        provenanceNote = string.IsNullOrWhiteSpace(provenanceNote) ? null : provenanceNote.Trim();
        if (provenanceNote?.Length > 500)
            throw new InvalidDataException("The upload provenance note exceeds 500 characters.");
        var extension = Path.GetExtension(fileName).ToLowerInvariant();
        if (extension is not (".pdb" or ".cif" or ".mmcif"))
            throw new InvalidDataException("Only identified PDB or mmCIF coordinate files can be selected.");
        var directory = WorkDirectory("uploads", NewId());
        var path = Path.Combine(directory, $"source{extension}");
        await using (var output = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None))
        {
            var buffer = new byte[64 * 1024];
            long count = 0;
            int read;
            while ((read = await input.ReadAsync(buffer, cancellationToken)) > 0)
            {
                count += read;
                if (count > 100_000_000) throw new InvalidDataException("The coordinate source exceeds 100 MB.");
                await output.WriteAsync(buffer.AsMemory(0, read), cancellationToken);
            }
            if (count == 0) throw new InvalidDataException("The coordinate source is empty.");
        }
        var token = NewId();
        var source = new StructuralSource($"upload:{token}", SourceRouteKind.Upload,
            $"Local upload: {Path.GetFileName(fileName)}; researcher-declared {provenance.ToString().ToLowerInvariant()} provenance, not independently verified",
            path, Hash(path), null, Path.GetFileName(fileName), null, provenance, provenanceNote);
        lock (_gate) _uploads.Add(token, source);
        return token;
    }

    public string? StructurePath(string token)
    {
        lock (_gate)
            return _structures.TryGetValue(token, out var binding) && IsWorkspaceFile(binding.Path)
                ? binding.Path : null;
    }

    /// <summary>Returns only the bytes selected for this exact structure URL.</summary>
    public (byte[] Bytes, string Extension)? VerifiedStructureContent(string token)
    {
        StructureBinding binding;
        lock (_gate)
            if (!_structures.TryGetValue(token, out binding!)) return null;
        if (!IsWorkspaceFile(binding.Path)) return null;
        try
        {
            var bytes = File.ReadAllBytes(binding.Path);
            var actual = Convert.ToHexString(SHA256.HashData(bytes)).ToLowerInvariant();
            return actual == binding.Sha256 ? (bytes, binding.Extension) : null;
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        {
            return null;
        }
    }

    /// <summary>Reads only a retained, exact-trial diagnostic whose bytes still match its digest.</summary>
    public (byte[] Bytes, string FileName, string Extension)? VerifiedDiagnosticContent(string token)
    {
        DiagnosticBinding binding;
        lock (_gate)
        {
            if (!_diagnostics.TryGetValue(token, out binding!)) return null;
            if (!_trialsByAttempt.TryGetValue(binding.AttemptId, out var trials) ||
                trials.IsDefaultOrEmpty || !trials.Any(trial =>
                    trial.TrialId == binding.TrialId && !trial.DiagnosticArtifacts.IsDefaultOrEmpty &&
                    trial.DiagnosticArtifacts.Any(artifact =>
                        artifact.Role == binding.Role && artifact.Sha256 == binding.Sha256 &&
                        artifact.DownloadUrl == "/api/diagnostics/" + token)))
                return null;
        }
        if (!IsWorkspaceFile(binding.Path)) return null;
        try
        {
            var bytes = File.ReadAllBytes(binding.Path);
            if (bytes.LongLength > 512L * 1024 * 1024 ||
                !string.Equals(Convert.ToHexString(SHA256.HashData(bytes)),
                    binding.Sha256, StringComparison.OrdinalIgnoreCase)) return null;
            return (bytes, binding.FileName, Path.GetExtension(binding.FileName).ToLowerInvariant());
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        {
            return null;
        }
    }

    /// <summary>Projects exact molecular roles for the selected full-system view.</summary>
    public InspectionComponentsAccount? ResolveInspectionComponents(string subjectId, string structureToken)
    {
        if (string.IsNullOrWhiteSpace(subjectId) || string.IsNullOrWhiteSpace(structureToken)) return null;
        lock (_gate)
        {
            var inspection = _inspection.Current;
            if (inspection is null || inspection.SubjectId != subjectId ||
                inspection.StructureUrl is null ||
                !string.Equals(inspection.StructureUrl.Split('?', 2)[0],
                    "/api/structures/" + structureToken, StringComparison.Ordinal) ||
                VerifiedStructureContent(structureToken) is null)
                return null;

            SourceToResultCorrespondence? correspondence;
            int atomCount;
            if (_constructedByAttempt.Values.FirstOrDefault(item => item.Id == subjectId) is
                { } constructed && constructed.Attempt.StudyRevisionId == inspection.StudyRevisionId)
            {
                correspondence = constructed.Correspondence;
                atomCount = constructed.Molecule.AtomCount;
                if (correspondence.ResultId != constructed.Molecule.CoordinateSha256) return null;
            }
            else if (_stages.TryGetValue(subjectId, out var stage) &&
                     stage.Attempt.StudyRevisionId == inspection.StudyRevisionId)
            {
                correspondence = stage.Correspondence;
                atomCount = stage.Molecule.AtomCount;
                if (correspondence.ResultId != stage.Molecule.Id) return null;
            }
            else return null;

            if (correspondence is not { Complete: true } || correspondence.Atoms.IsDefault ||
                correspondence.Atoms.Length != atomCount || atomCount <= 0) return null;
            var atoms = new AtomCorrespondence?[atomCount];
            foreach (var atom in correspondence.Atoms)
            {
                if (atom.ResultAtomIndex < 0 || atom.ResultAtomIndex >= atomCount ||
                    atoms[atom.ResultAtomIndex] is not null || !Enum.IsDefined(atom.MoleculeRole))
                    return null;
                atoms[atom.ResultAtomIndex] = atom;
            }
            var runs = ImmutableArray.CreateBuilder<InspectionComponentRun>();
            var start = 0;
            for (var index = 0; index < atomCount; index++)
            {
                var atom = atoms[index];
                if (atom is null) return null;
                var chain = atom.SourceResidue?.Chain;
                var copyId = atom.SourceResidue?.CopyId;
                if (index == 0) continue;
                var first = atoms[start]!;
                if (atom.MoleculeRole == first.MoleculeRole &&
                    chain == first.SourceResidue?.Chain && copyId == first.SourceResidue?.CopyId)
                    continue;
                runs.Add(new InspectionComponentRun(start, index, first.MoleculeRole,
                    first.SourceResidue?.Chain, first.SourceResidue?.CopyId));
                start = index;
            }
            var last = atoms[start]!;
            runs.Add(new InspectionComponentRun(start, atomCount, last.MoleculeRole,
                last.SourceResidue?.Chain, last.SourceResidue?.CopyId));
            return new InspectionComponentsAccount(subjectId, inspection.StudyRevisionId,
                structureToken, atomCount, runs.ToImmutable());
        }
    }

    /// <summary>Resolves a picked coordinate row only through the selected subject's correspondence.</summary>
    public InspectionAtomAccount? ResolveInspectionAtom(string subjectId, string structureToken,
        int atomSiteIndex)
    {
        if (atomSiteIndex < 0 || string.IsNullOrWhiteSpace(subjectId) ||
            string.IsNullOrWhiteSpace(structureToken)) return null;
        lock (_gate)
        {
            var inspection = _inspection.Current;
            if (inspection is null || inspection.SubjectId != subjectId ||
                inspection.StructureUrl is null ||
                !string.Equals(inspection.StructureUrl.Split('?', 2)[0],
                    "/api/structures/" + structureToken, StringComparison.Ordinal) ||
                VerifiedStructureContent(structureToken) is null)
                return null;

            SourceToResultCorrespondence? correspondence = null;
            int atomCount = 0;
            if (_protein?.Id == subjectId && _protein.StudyRevisionId == inspection.StudyRevisionId)
            {
                correspondence = _protein.Correspondence;
                atomCount = _protein.Molecule.AtomCount;
                if (correspondence.ResultId != _protein.Molecule.CoordinateSha256) return null;
            }
            else if (_constructedByAttempt.Values.FirstOrDefault(item => item.Id == subjectId) is
                     { } constructed && constructed.Attempt.StudyRevisionId == inspection.StudyRevisionId)
            {
                correspondence = constructed.Correspondence;
                atomCount = constructed.Molecule.AtomCount;
                if (correspondence.ResultId != constructed.Molecule.CoordinateSha256) return null;
            }
            else if (_stages.TryGetValue(subjectId, out var stage) &&
                     stage.Attempt.StudyRevisionId == inspection.StudyRevisionId)
            {
                correspondence = stage.Correspondence;
                atomCount = stage.Molecule.AtomCount;
                if (correspondence.ResultId != stage.Molecule.Id) return null;
            }
            if (correspondence is not { Complete: true } || correspondence.Atoms.IsDefault ||
                correspondence.Atoms.Length != atomCount || atomSiteIndex >= atomCount)
                return null;
            var matching = correspondence.Atoms.Where(item => item.ResultAtomIndex == atomSiteIndex)
                .Take(2).ToArray();
            if (matching.Length != 1) return null;
            return new InspectionAtomAccount(subjectId, inspection.StudyRevisionId,
                structureToken, atomSiteIndex, matching[0]);
        }
    }

    /// <summary>Read and hash the exact published bytes before the host sends them.</summary>
    public ExportDeliveryObservation VerifiedExportContent(string stageId, string? expectedSha256)
    {
        ExportDeliveryObservation answer;
        var changed = false;
        lock (_gate)
        {
            if (!_bundles.TryGetValue(stageId, out var bundle) ||
                !_stageAssessments.TryGetValue(stageId, out var assessment) ||
                !assessment.CurrentlyApplicable || assessment.Id != bundle.AssessmentId)
                return new ExportDeliveryObservation(ExportDeliveryStanding.Absent,
                    "No current verified export is available for this completed stage.", null, null, null);
            if (!string.Equals(expectedSha256?.Trim('"'), bundle.Sha256, StringComparison.OrdinalIgnoreCase))
                return new ExportDeliveryObservation(ExportDeliveryStanding.IdentityChanged,
                    "The requested bundle digest does not match this stage's current export.", null, null, null);
            byte[]? bytes = null;
            try
            {
                if (IsWorkspaceFile(bundle.BundlePath)) bytes = File.ReadAllBytes(bundle.BundlePath);
            }
            catch (Exception exception) when (exception is IOException or UnauthorizedAccessException) { }
            if (bytes is not null && bytes.LongLength == bundle.ByteLength &&
                Convert.ToHexString(SHA256.HashData(bytes)).Equals(bundle.Sha256,
                    StringComparison.OrdinalIgnoreCase))
                answer = new ExportDeliveryObservation(ExportDeliveryStanding.Available, null,
                    bytes, bundle.Sha256, bundle.ByteLength);
            else
            {
                _bundles.Remove(stageId);
                _exportAccounts[stageId] = new StageExportAccount(stageId, assessment.Id, "failed",
                    "The published bundle bytes no longer match the verified export. Retry export for this completed stage.",
                    null, null);
                TouchLocked();
                changed = true;
                answer = new ExportDeliveryObservation(ExportDeliveryStanding.Corrupt,
                    _exportAccounts[stageId].Reason, null, null, null);
            }
        }
        if (changed) Changed?.Invoke();
        return answer;
    }

    private async Task<string?> SearchSourceAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var query = Text(data, "query");
        if (query is null) return "Enter a structural identifier or descriptive protein query.";
        var result = await _sources.SearchAsync(query, cancellationToken);
        lock (_gate)
        {
            _candidates = result.Candidates;
            _sourceSearchIssues = result.UnavailableRoutes;
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private async Task<string?> SelectSourceAsync(JsonElement data, CancellationToken cancellationToken)
    {
        StructuralSource source;
        var token = Text(data, "uploadToken");
        var exactIdentifier = Text(data, "exactIdentifier");
        var exactKind = Text(data, "sourceKind");
        var sourceId = Text(data, "sourceId");
        if (new[] { token is not null, exactIdentifier is not null || exactKind is not null,
                sourceId is not null }.Count(present => present) != 1)
            return "Supply exactly one uploaded source, discovered candidate, or exact database reference.";
        if (token is not null)
        {
            lock (_gate)
            {
                if (!_uploads.TryGetValue(token, out var uploaded))
                    return "The uploaded source is not available in this local session.";
                source = uploaded;
            }
        }
        else
        {
            CandidateSourceRecord? candidate;
            if (exactIdentifier is not null || exactKind is not null)
            {
                var kind = exactKind switch
                {
                    "rcsb" => SourceRouteKind.Rcsb,
                    "alphafold" => SourceRouteKind.AlphaFold,
                    _ => (SourceRouteKind?)null
                };
                candidate = kind is null || exactIdentifier is null ? null :
                    ExternalSourceExchange.IdentifyExactReference(kind.Value, exactIdentifier);
                if (candidate is null)
                    return "Supply an exact four-character RCSB PDB entry or full AlphaFold DB model or fragment ID.";
            }
            else
            {
                lock (_gate) candidate = _candidates.FirstOrDefault(item => item.Id == sourceId);
                if (candidate is null) return "Select a source from the current candidate account.";
            }
            try { source = await _sources.RetrieveAsync(candidate, WorkDirectory("intake", NewId()), cancellationToken); }
            catch (Exception exception) when (exception is HttpRequestException or IOException or JsonException or InvalidDataException)
            { return "The selected source could not be retrieved and identified: " + exception.Message; }
        }
        if (!IsWorkspaceFile(source.CoordinatePath) || Hash(source.CoordinatePath) != source.Sha256 ||
            !ValidSourceEvidenceIdentity(source))
            return "The source bytes do not correspond to the selected source identity.";
        var inspected = await _proteinPreparation.InspectSourceAsync(source, WorkDirectory("inspection", NewId()),
            _catalogue?.MaximumSourceAtoms is > 0 ? _catalogue.MaximumSourceAtoms : int.MaxValue, cancellationToken);
        if (inspected.Value is null) return inspected.Reason;
        lock (_gate)
        {
            _selectedSource = source;
            _sourceInspection = inspected.Value;
            _sourcePreviews.Clear();
            AdvanceStudyLocked(null, _study.Membrane, null);
            var carriedMembrane = CarryMembraneLocked(_study);
            _preparationProposals = null;
            _preparationPlan = null;
            _authorizedPreparationPlan = null;
            _planOverrides = ImmutableArray<ResidueVariantChoice>.Empty;
            _recommendationRunning = false;
            _recommendationIssue = null;
            _decisions = ImmutableArray<ResearcherDecision>.Empty;
            _protein = null;
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementMeasurementIssue = null;
            _placementWitness = null;
            _placement = null;
            _membrane = carriedMembrane;
            var sourceSubject = GenericSubject(source.Id, _study.Id, source.CoordinatePath,
                "structuralSource", ImmutableArray<ScientificEvidence>.Empty,
                ImmutableArray<ScientificFinding>.Empty, null);
            var selected = _inspection.Select(_study, sourceSubject);
            if (selected.Value is null)
                throw new InvalidOperationException(selected.Reason ?? "The selected source could not be opened for inspection.");
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private async Task<string?> SelectProteinModelAsync(JsonElement data, CancellationToken cancellationToken)
    {
        StructuralSource source;
        SourceInspectionReport inspection;
        lock (_gate)
        {
            if (_selectedSource is null || _sourceInspection is null) return "Inspect an exact structural source first.";
            source = _selectedSource;
            inspection = _sourceInspection;
        }
        var modelIndex = Integer(data, "modelIndex");
        var model = inspection.Models.FirstOrDefault(item => item.Index == modelIndex);
        if (model is null) return "The selected coordinate model was not observed in this source.";
        var assembly = Text(data, "biologicalAssemblyId");
        if (assembly is not null && !model.Assemblies.Any(item => item.Name == assembly))
            return "The selected biological assembly is not in the inspected source.";
        var chains = ParseArray<ChainSelection>(data, "chains");
        var partners = ParseArray<PartnerSelection>(data, "partners");
        var rawAltlocs = ParseArray<AlternateLocationChoice>(data, "alternateLocations");
        var selectedAssembly = assembly is null ? null : model.Assemblies.First(item => item.Name == assembly);
        var relevantPartners = PartnersForSelection(model, selectedAssembly, chains);
        // An assembly CopyId is the full observed output chain ID (for example A1),
        // not a suffix or an independently invented copy number.
        if (chains.IsDefaultOrEmpty || chains.Any(item => string.IsNullOrWhiteSpace(item.CopyId) ||
                !model.Chains.Any(observed => observed.Name == item.SourceChain) ||
                !model.Residues.Any(residue => residue.ResidueKind == SourceResidueKind.Protein &&
                    residue.Address.Chain == item.SourceChain)) ||
            chains.Select(item => item.CopyId).Distinct(StringComparer.Ordinal).Count() != chains.Length ||
            selectedAssembly is not null && chains.Any(item => !selectedAssembly.ChainCopies.Contains(item)) ||
            selectedAssembly is null && chains.Any(item => item.CopyId != item.SourceChain) ||
            partners.Length != relevantPartners.Length ||
            partners.Select(item => item.SourceId).Distinct(StringComparer.Ordinal).Count() != partners.Length ||
            partners.Any(item => !relevantPartners.Any(observed => observed.SourceId == item.SourceId)) ||
            rawAltlocs.Any(choice => choice.Residue.CopyId != string.Empty ||
                choice.Residue.Model != model.Index ||
                !model.Residues.Any(residue => residue.Address == choice.Residue &&
                    residue.AlternateLocations.Contains(choice.Altloc))) ||
            rawAltlocs.Select(choice => choice.Residue).Distinct().Count() != rawAltlocs.Length)
            return "Choose observed protein chain copies and explicit dispositions for all partners in the selected membership.";
        var altlocs = rawAltlocs.Select(item => item with { DecisionId = NewId() }).ToImmutableArray();
        var intended = new IntendedProteinModel(NewId(), source, modelIndex!.Value, assembly, chains, partners, altlocs,
            model.SourceModelId);
        StudyRevision revision;
        ProteinChemicalStatePolicy? chemicalPolicy;
        ProteinStructuralAssessmentPolicy? structuralPolicy;
        ProteinSelectionIssue? selectionIssue;
        lock (_gate)
        {
            AdvanceStudyLocked(intended, _study.Membrane, null);
            revision = _study;
            var chemicalCompatibility = SelectChemicalPolicyWithCompatibilityLocked(model, intended);
            chemicalPolicy = chemicalCompatibility.Policy;
            structuralPolicy = SelectStructuralPolicyLocked(model, intended);
            selectionIssue = chemicalCompatibility.UnsupportedPartners.Length > 0 &&
                structuralPolicy is not null && chemicalCompatibility.CorePolicyAvailable &&
                chemicalCompatibility.AllPartnersKnown
                    ? new ProteinSelectionIssue(
                        "This preparation method does not yet support the selected molecules: " +
                        string.Join(", ", chemicalCompatibility.UnsupportedPartners.Select(UnsupportedPartnerLabel)) + ".", true)
                    : null;
            var carriedMembrane = CarryMembraneLocked(revision);
            _preparationProposals = null;
            _preparationPlan = null;
            _authorizedPreparationPlan = null;
            _planOverrides = ImmutableArray<ResidueVariantChoice>.Empty;
            _recommendationRunning = false;
            _recommendationIssue = null;
            _decisions = ImmutableArray<ResearcherDecision>.Empty;
            _protein = null;
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementMeasurementIssue = null;
            _placementWitness = null;
            _placement = null;
            _membrane = carriedMembrane;
            _proteinSelectionRunning = true;
            TouchLocked();
        }
        Changed?.Invoke();
        if (chemicalPolicy is null || structuralPolicy is null)
        {
            lock (_gate)
            {
                if (_study.Id == revision.Id)
                {
                    _proteinSelectionRunning = false;
                    _proteinSelectionIssue = selectionIssue ?? new ProteinSelectionIssue(
                        "No qualified protein chemical-state and structural assessment policies cover this selected structure.");
                    TouchLocked();
                }
            }
            Changed?.Invoke();
            SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = source.Id }));
            return null;
        }
        BoundaryOutcome<PreparationProposalReport> proposed;
        try
        {
            proposed = await _proteinPreparation.ProposeChangesAsync(revision, intended, inspection, chemicalPolicy,
                structuralPolicy, WorkDirectory("protein-proposals", revision.Id), cancellationToken);
        }
        catch (Exception exception)
        {
            lock (_gate)
            {
                if (_study.Id != revision.Id) return null;
                _proteinSelectionRunning = false;
                _proteinSelectionIssue = new ProteinSelectionIssue(
                    $"The selected protein could not be assessed: {exception.Message}");
                TouchLocked();
            }
            Changed?.Invoke();
            return null;
        }
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _proteinSelectionRunning = false;
            _preparationProposals = proposed.Value;
            _proteinSelectionIssue = proposed.Value is null && proposed.Reason is { } reason
                ? new ProteinSelectionIssue(reason) : null;
            TouchLocked();
        }
        Changed?.Invoke();
        SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = intended.Id }));
        if (proposed.Value is not null)
            await CalculatePreparationPlanAsync(revision, cancellationToken);
        return null;
    }

    private async Task CalculatePreparationPlanAsync(StudyRevision revision, CancellationToken cancellationToken)
    {
        IntendedProteinModel? intended;
        SourceInspectionReport? inspection;
        PreparationProposalReport? proposals;
        ProteinChemicalStatePolicy? policy;
        ProteinStructuralAssessmentPolicy? structuralPolicy;
        ImmutableArray<ResidueVariantChoice> overrides;
        ImmutableArray<ResearcherDecision> decisions;
        long generation;
        lock (_gate)
        {
            if (_study.Id != revision.Id) return;
            generation = ++_recommendationGeneration;
            intended = revision.IntendedProtein;
            inspection = _sourceInspection;
            proposals = _preparationProposals;
            policy = inspection is null ? null : SelectChemicalPolicyLocked(
                inspection.Models.FirstOrDefault(item => item.Index == intended?.ModelIndex), intended);
            structuralPolicy = inspection is null ? null : SelectStructuralPolicyLocked(
                inspection.Models.FirstOrDefault(item => item.Index == intended?.ModelIndex), intended);
            overrides = _planOverrides;
            decisions = _decisions;
            _preparationPlan = null;
            _recommendationIssue = null;
            var unresolvedExceptions = proposals?.Changes.Any(change =>
                (change.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.Disulfide) &&
                !decisions.Any(decision => decision.SubjectId == change.Id &&
                    (change.Kind == PreparationChangeKind.Disulfide ||
                     decision.ChosenValue == ResearcherDecisionValue.Approved))) == true;
            if (intended is null || inspection is null || proposals is null || policy is null ||
                structuralPolicy is null || !proposals.UnresolvedQuestions.IsDefaultOrEmpty ||
                unresolvedExceptions)
            {
                _recommendationRunning = false;
                _recommendationIssue = "The observed structural exceptions need individual review before a complete recommendation plan is available.";
                TouchLocked();
                return;
            }
            _recommendationRunning = true;
            TouchLocked();
        }
        Changed?.Invoke();
        BoundaryOutcome<RecommendedPreparationPlan> calculated;
        try
        {
            calculated = await _proteinPreparation.RecommendAsync(revision, intended!, inspection!, policy!,
                structuralPolicy!, proposals!, decisions, overrides,
                WorkDirectory("protein-recommendation", NewId()), cancellationToken);
        }
        catch (Exception exception)
        {
            calculated = BoundaryOutcome<RecommendedPreparationPlan>.Unavailable(
                $"The starting-state method could not establish a plan: {exception.Message}");
        }
        lock (_gate)
        {
            if (_study.Id != revision.Id || generation != _recommendationGeneration ||
                !_planOverrides.SequenceEqual(overrides) || !_decisions.SequenceEqual(decisions)) return;
            _recommendationRunning = false;
            _preparationPlan = calculated.Value;
            _recommendationIssue = calculated.Value is null ? calculated.Reason : null;
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private async Task<string?> RetryPreparationPlanAsync(CancellationToken cancellationToken)
    {
        StudyRevision revision;
        lock (_gate)
        {
            if (_study.IntendedProtein is null || _preparationProposals is null ||
                _recommendationRunning || _preparationPlan is not null ||
                _authorizedPreparationPlan is not null)
                return "There is no failed current recommendation calculation to retry.";
            revision = _study;
        }
        await CalculatePreparationPlanAsync(revision, cancellationToken);
        return null;
    }

    private async Task<string?> StartProteinPreparationAsync(CancellationToken cancellationToken)
    {
        StudyRevision revision;
        lock (_gate)
        {
            var review = BuildPreparationReviewLocked();
            if (_protein is not null || _proteinPreparationRunning ||
                _authorizedPreparationPlan is not null || _preparationPlan is not null ||
                _preparationProposals?.StudyRevisionId != _study.Id ||
                review is not { Blockers.Length: 0, RemainingCount: 0 } ||
                !EveryRequiredChoiceSettled(_preparationProposals, _decisions))
                return "Finish the exact manual decisions or authorize the current checked plan before preparation.";
            revision = _study;
        }
        await TryPrepareProteinAsync(revision, cancellationToken);
        return null;
    }

    private async Task<string?> OverridePreparationPlanChoiceAsync(JsonElement data,
        CancellationToken cancellationToken)
    {
        var proposalId = Text(data, "proposalId");
        StudyRevision revision;
        lock (_gate)
        {
            if (_recommendationRunning || _authorizedPreparationPlan is not null || _protein is not null)
                return "The current plan cannot be changed while it is running or after authorization.";
            var proposal = _preparationProposals?.Changes.FirstOrDefault(change =>
                change.Id == proposalId && change.Kind == PreparationChangeKind.ResidueState);
            if (proposal is null || proposal.StudyRevisionId != _study.Id)
                return "Choose a current chemical-state option at its exact residue.";
            _planOverrides = _planOverrides.Where(item => item.Residue != proposal.Residue)
                .Append(new ResidueVariantChoice(proposal.Residue, proposal.ProposedChange, NewId()))
                .ToImmutableArray();
            revision = _study;
            _preparationPlan = null;
            TouchLocked();
        }
        Changed?.Invoke();
        await CalculatePreparationPlanAsync(revision, cancellationToken);
        return null;
    }

    private async Task<string?> AuthorizePreparationPlanAsync(JsonElement data,
        CancellationToken cancellationToken)
    {
        var expectedDigest = Text(data, "planSha256");
        StudyRevision revision;
        lock (_gate)
        {
            var plan = _preparationPlan;
            var intended = _study.IntendedProtein;
            var proposals = _preparationProposals;
            var model = _sourceInspection?.Models.FirstOrDefault(item => item.Index == intended?.ModelIndex);
            var policy = SelectChemicalPolicyLocked(model, intended);
            if (plan is null || proposals is null || intended is null || policy is null ||
                _recommendationRunning || _proteinPreparationRunning || _protein is not null ||
                _authorizedPreparationPlan is not null ||
                expectedDigest != plan.PlanSha256 || plan.StudyRevisionId != _study.Id ||
                plan.IntendedProteinId != intended.Id || plan.SourceSha256 != intended.Source.Sha256 ||
                plan.ChemicalPolicyId != policy.Id || plan.ChemicalPolicyVersion != policy.Version ||
                !IsWorkspaceFile(plan.Candidate.Path) || !File.Exists(plan.Candidate.Path) ||
                !Hash(plan.Candidate.Path).Equals(plan.Candidate.Sha256, StringComparison.OrdinalIgnoreCase) ||
                !proposals.UnresolvedQuestions.IsDefaultOrEmpty ||
                !plan.PrerequisiteDecisionIds.SequenceEqual(_decisions.Select(decision => decision.Id)
                    .Order(StringComparer.Ordinal)) ||
                BuildPreparationReviewLocked() is not { Blockers.Length: 0 })
                return "The exact checked recommendation plan is absent, changed, or already authorized. Recalculate it before preparation.";
            var newDecisions = ImmutableArray.CreateBuilder<ResearcherDecision>();
            newDecisions.AddRange(_decisions);
            foreach (var group in proposals.Changes.GroupBy(DecisionScope))
            {
                var recorded = group.SelectMany(change => _decisions.Where(decision =>
                    decision.SubjectId == change.Id)).ToArray();
                if (recorded.Length > 0)
                {
                    if (group.Key.Kind == PreparationChangeKind.ResidueState &&
                        (recorded.Length != 1 || !group.Any(change => change.Id == recorded[0].SubjectId &&
                            recorded[0].ChosenValue == ResearcherDecisionValue.Approved &&
                            plan.Choices.Any(choice => choice.Residue == change.Residue &&
                                choice.Variant == change.ProposedChange))))
                        return "A recorded state differs from the checked plan; recalculate it.";
                    if (group.Key.Kind == PreparationChangeKind.HeavyAtom &&
                        (recorded.Length != 1 || recorded[0].ChosenValue != ResearcherDecisionValue.Approved))
                        return "A required repair is not approved for this plan.";
                    if (group.Key.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.Disulfide)
                        continue;
                    if (group.Key.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.HeavyAtom)
                        continue;
                }
                var chosen = group.Key.Kind switch
                {
                    PreparationChangeKind.HeavyAtom => group.SingleOrDefault(),
                    PreparationChangeKind.ResidueState => group.SingleOrDefault(item =>
                        plan.Choices.Any(choice => choice.Residue == item.Residue &&
                            choice.Variant == item.ProposedChange)),
                    _ => null
                };
                if (chosen is null)
                    return "A required exact preparation choice is missing from the checked plan.";
                newDecisions.Add(new ResearcherDecision(NewId(), _study.Id, chosen.Id,
                    ResearcherDecisionKind.ApprovePreparationChange, ResearcherDecisionValue.Approved,
                    DateTimeOffset.UtcNow, null));
            }
            _decisions = newDecisions.ToImmutable();
            _allDecisions = _allDecisions.AddRange(newDecisions.Where(decision =>
                !plan.PrerequisiteDecisionIds.Contains(decision.Id)));
            _authorizedPreparationPlan = plan;
            revision = _study;
            TouchLocked();
        }
        Changed?.Invoke();
        await TryPrepareProteinAsync(revision, cancellationToken);
        return null;
    }

    private async Task<string?> DecidePreparationChangeAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var proposalId = Text(data, "proposalId");
        var approved = Boolean(data, "approve");
        // Earlier decisions may carry an authored annotation; new decisions need only
        // the exact choice and applicable scientific checks.
        var rationale = Text(data, "rationale");
        PreparationChangeProposal? proposal;
        StudyRevision revision;
        bool recommendationConfigured;
        lock (_gate)
        {
            proposal = _preparationProposals?.Changes.FirstOrDefault(item => item.Id == proposalId);
            revision = _study;
            if (proposal is null || proposal.StudyRevisionId != revision.Id) return "This preparation proposal is absent or stale.";
            if (approved is null) return "The proposal decision must explicitly approve or decline.";
            if (!approved.Value && proposal.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState)
                return "Choose one alternative for this site; an exclusive state or conformer cannot be declined individually.";
            if (_protein is not null) return "The prepared protein has already been established from the reviewed choices.";
            if (_decisions.Any(item => item.SubjectId == proposal.Id && item.StudyRevisionId == revision.Id))
                return "This exact proposal has already been decided.";
            if (BuildPreparationReviewLocked()?.Decisions.FirstOrDefault(item =>
                    item.Options.Any(option => option.ProposalId == proposal.Id))?.Blocker is { } siteBlocker)
                return siteBlocker;
            if (approved.Value && (proposal.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState) &&
                _preparationProposals!.Changes.Any(other => other.Id != proposal.Id &&
                    other.Kind == proposal.Kind && other.Residue == proposal.Residue &&
                    _decisions.Any(decision => decision.SubjectId == other.Id &&
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved)))
                return "An alternative for this exact residue has already been approved.";
            if (approved.Value && PreparationOptionPrerequisiteLocked(proposal, _preparationProposals!) is { } reason)
                return reason;
            var decision = new ResearcherDecision(NewId(), revision.Id, proposal.Id,
                ResearcherDecisionKind.ApprovePreparationChange,
                approved.Value ? ResearcherDecisionValue.Approved : ResearcherDecisionValue.Declined,
                DateTimeOffset.UtcNow, rationale);
            _decisions = _decisions.Add(decision);
            _allDecisions = _allDecisions.Add(decision);
            if (proposal.Kind == PreparationChangeKind.ResidueState)
                _planOverrides = _planOverrides.Where(item => item.Residue != proposal.Residue)
                    .ToImmutableArray();
            _preparationPlan = null;
            _recommendationIssue = null;
            recommendationConfigured = _sourceInspection is not null &&
                SelectChemicalPolicyLocked(_sourceInspection.Models.FirstOrDefault(item =>
                    item.Index == revision.IntendedProtein?.ModelIndex), revision.IntendedProtein)?
                    .Recommendation is not null;
            TouchLocked();
        }
        Changed?.Invoke();
        if (recommendationConfigured)
            await CalculatePreparationPlanAsync(revision, cancellationToken);
        else
            await TryPrepareProteinAsync(revision, cancellationToken);
        return null;
    }

    private async Task TryPrepareProteinAsync(StudyRevision revision, CancellationToken cancellationToken)
    {
        IntendedProteinModel? intended;
        SourceInspectionReport? inspection;
        PreparationProposalReport? proposals;
        ImmutableArray<ResearcherDecision> decisions;
        ProteinChemicalStatePolicy? policy;
        ProteinStructuralAssessmentPolicy? structuralPolicy;
        lock (_gate)
        {
            intended = revision.IntendedProtein;
            inspection = _sourceInspection;
            proposals = _preparationProposals;
            decisions = _decisions;
            policy = inspection is null ? null : SelectChemicalPolicyLocked(
                inspection.Models.FirstOrDefault(item => item.Index == intended?.ModelIndex), intended);
            structuralPolicy = inspection is null ? null : SelectStructuralPolicyLocked(
                inspection.Models.FirstOrDefault(item => item.Index == intended?.ModelIndex), intended);
            if (_protein is null && !_proteinPreparationRunning &&
                intended is not null && inspection is not null && proposals is not null &&
                policy is not null && structuralPolicy is not null &&
                proposals.UnresolvedQuestions.IsDefaultOrEmpty &&
                EveryRequiredChoiceSettled(proposals, decisions) &&
                BuildPreparationReviewLocked() is { Blockers.Length: 0, RemainingCount: 0 })
            {
                _proteinPreparationRunning = true;
                _proteinPreparationFailure = null;
                _proteinPreparationRetryable = false;
                _proteinPreparationObservationIssue = null;
                TouchLocked();
            }
        }
        if (intended is null || inspection is null || proposals is null || policy is null || structuralPolicy is null ||
            !proposals.UnresolvedQuestions.IsDefaultOrEmpty ||
            !EveryRequiredChoiceSettled(proposals, decisions) || !_proteinPreparationRunning) return;
        Changed?.Invoke();
        BoundaryOutcome<AssessedPreparedProtein> prepared;
        try
        {
            prepared = await _proteinPreparation.PrepareAsync(revision, intended, inspection, policy, structuralPolicy,
                proposals, decisions,
                WorkDirectory("protein-preparation", revision.Id), cancellationToken,
                _authorizedPreparationPlan?.StudyRevisionId == revision.Id,
                _authorizedPreparationPlan?.StudyRevisionId == revision.Id ? _authorizedPreparationPlan : null);
        }
        catch (Exception exception)
        {
            lock (_gate)
            {
                if (_study.Id != revision.Id) return;
                _proteinPreparationRunning = false;
                if (exception is OperationCanceledException)
                    _proteinPreparationObservationIssue = "The preparation outcome could not be observed after the request was interrupted. Refresh the account before any new operation.";
                else
                {
                    _proteinPreparationFailure = $"Protein preparation stopped before an assessed result: {exception.Message}";
                    _proteinPreparationRetryable = exception is IOException or TimeoutException or HttpRequestException;
                }
                TouchLocked();
            }
            Changed?.Invoke();
            return;
        }
        lock (_gate)
        {
            if (_study.Id != revision.Id) return;
            _proteinPreparationRunning = false;
            _proteinPreparationObservationIssue = null;
            _protein = prepared.Value;
            _proteinDiagnostic = prepared.Value is null ? prepared.Diagnostic as ProteinPreparationDiagnostic : null;
            _proteinPreparationFailure = prepared.Value is null ? prepared.Reason : null;
            var exchangeFailure = prepared.Diagnostic as ProteinPreparationExchangeDiagnostic;
            _proteinPreparationRetryable = exchangeFailure?.Standing == WorkerResultStanding.Failed ||
                prepared.Diagnostic is ProteinPreparationInputIntegrityDiagnostic;
            if (exchangeFailure?.Standing == WorkerResultStanding.Unobserved)
            {
                _proteinPreparationObservationIssue = prepared.Reason;
                _proteinPreparationFailure = null;
            }
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private async Task<string?> RetryProteinPreparationAsync(CancellationToken cancellationToken)
    {
        StudyRevision revision;
        lock (_gate)
        {
            if (_proteinPreparationFailure is null || _protein is not null ||
                !_proteinPreparationRetryable || _proteinPreparationRunning ||
                _preparationProposals?.StudyRevisionId != _study.Id)
                return "There is no failed preparation for the current reviewed protein to retry.";
            revision = _study;
        }
        await TryPrepareProteinAsync(revision, cancellationToken);
        return null;
    }

    private static bool EveryRequiredChoiceSettled(PreparationProposalReport proposals, ImmutableArray<ResearcherDecision> decisions)
    {
        if (proposals.Changes.IsDefaultOrEmpty) return true;
        if (proposals.Changes.Any(change => change.Kind == PreparationChangeKind.HeavyAtom &&
            !decisions.Any(decision => decision.SubjectId == change.Id &&
                decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                decision.ChosenValue == ResearcherDecisionValue.Approved))) return false;
        if (proposals.Changes.Any(change => change.Kind == PreparationChangeKind.Disulfide &&
            !decisions.Any(decision => decision.SubjectId == change.Id &&
                decision.Kind == ResearcherDecisionKind.ApprovePreparationChange))) return false;
        return proposals.Changes.Where(change => change.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState)
            .GroupBy(change => (change.Residue, change.Kind))
            .All(group => group.Count(change => decisions.Any(decision => decision.SubjectId == change.Id &&
                decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                decision.ChosenValue == ResearcherDecisionValue.Approved)) == 1);
    }

    private sealed record PreparationDecisionScope(PreparationChangeKind Kind, ResidueAddress Residue,
        ResidueAddress? PartnerResidue, string? AtomName);

    private static ImmutableArray<SourcePartnerObservation> PartnersForSelection(
        SourceModelObservation model, SourceAssemblyObservation? assembly,
        ImmutableArray<ChainSelection> selectedProteinCopies) =>
        model.Partners.Where(partner =>
            (assembly is null || partner.Chain is null ||
                assembly.ChainCopies.Any(copy => copy.SourceChain == partner.Chain)) &&
            (partner.Chain is null ||
                !model.Residues.Any(residue => residue.ResidueKind == SourceResidueKind.Protein &&
                    residue.Address.Chain == partner.Chain) ||
                selectedProteinCopies.Any(copy => copy.SourceChain == partner.Chain)))
            .ToImmutableArray();

    private static PreparationDecisionScope DecisionScope(PreparationChangeProposal proposal)
    {
        var first = proposal.Residue;
        var second = proposal.PartnerResidue;
        if (proposal.Kind == PreparationChangeKind.Disulfide && second is not null &&
            string.CompareOrdinal(AddressKey(first), AddressKey(second)) > 0)
            (first, second) = (second, first);
        return new PreparationDecisionScope(proposal.Kind, first,
            proposal.Kind == PreparationChangeKind.Disulfide ? second : null,
            proposal.Kind == PreparationChangeKind.HeavyAtom ? proposal.ProposedChange : null);
    }

    private static string AddressKey(ResidueAddress address) =>
        $"{address.Model:D10}:{address.Chain.Length}:{address.Chain}:" +
        $"{address.CopyId.Length}:{address.CopyId}:{address.Residue:D10}:" +
        $"{address.InsertionCode.Length}:{address.InsertionCode}";

    private string? PreparationOptionPrerequisiteLocked(PreparationChangeProposal proposal,
        PreparationProposalReport report)
    {
        var evidence = report.Evidence.Where(item => item.SubjectId == proposal.Id).ToArray();
        if (evidence.Length == 0 || evidence.Any(item => string.IsNullOrWhiteSpace(item.Id) ||
                string.IsNullOrWhiteSpace(item.Observation) || string.IsNullOrWhiteSpace(item.Applicability)) ||
            evidence.Select(item => item.Id).Distinct(StringComparer.Ordinal).Count() != evidence.Length)
            return "The exact option has no attributable scientific evidence. Reassess the selected protein.";
        if (report.Preview is not { } preview || !IsWorkspaceFile(preview.Path) ||
            !VerifyHash(preview.Path, preview.Sha256))
            return "The selected-coordinate preview is missing or has changed. Reassess the selected protein.";
        // Build the same exact chain-copy/residue mapping used for inspection, without
        // requiring a browser selection, successful rendering, or an inspection click.
        var subject = _proteinPreparation.ProposalInspectionSubject(report, proposal.Id, "exact-preview");
        return subject.Value is null ? subject.Reason : null;
    }

    private ProteinPreparationReviewAccount? BuildPreparationReviewLocked()
    {
        var report = _preparationProposals;
        var intended = _study.IntendedProtein;
        if (report is null || intended is null || report.StudyRevisionId != _study.Id ||
            report.IntendedProteinId != intended.Id) return null;

        var blockers = ImmutableArray.CreateBuilder<string>();
        blockers.AddRange(report.UnresolvedQuestions);
        if (report.Changes.Select(change => change.Id).Distinct(StringComparer.Ordinal).Count() != report.Changes.Length)
            blockers.Add("The proposal report repeats an exact proposal identity; reassess the selected protein.");
        if (report.Changes.Any(change => change.StudyRevisionId != _study.Id ||
                change.IntendedProteinId != intended.Id))
            blockers.Add("The proposal report contains a change for another study or protein.");
        var model = _sourceInspection?.Models.FirstOrDefault(item => item.Index == intended.ModelIndex);
        if (SelectChemicalPolicyLocked(model, intended) is null ||
            SelectStructuralPolicyLocked(model, intended) is null)
            blockers.Add("An applicable protein preparation policy is unavailable for this selected model.");

        var scopes = report.Changes.GroupBy(DecisionScope).ToArray();
        var decisions = ImmutableArray.CreateBuilder<PreparationDecisionAccount>();
        foreach (var group in scopes)
        {
            var proposals = group.ToArray();
            var recorded = proposals.Select(proposal => (Proposal: proposal,
                Decisions: _decisions.Where(decision => decision.StudyRevisionId == _study.Id &&
                    decision.SubjectId == proposal.Id &&
                    decision.Kind == ResearcherDecisionKind.ApprovePreparationChange).ToArray())).ToArray();
            var approved = recorded.Where(item => item.Decisions.Any(decision =>
                decision.ChosenValue == ResearcherDecisionValue.Approved)).ToArray();
            var explicitlyDeclined = recorded.Where(item => item.Decisions.Any(decision =>
                decision.ChosenValue == ResearcherDecisionValue.Declined)).ToArray();
            string? blocker = null;
            if (group.Key.Kind is PreparationChangeKind.HeavyAtom or PreparationChangeKind.Disulfide &&
                proposals.Length != 1)
                blocker = "This exact repair or bond scope has duplicate proposals; reassess the selected protein.";
            else if (proposals.GroupBy(proposal => proposal.ProposedChange, StringComparer.Ordinal)
                    .Any(variants => variants.Count() != 1))
                blocker = "This exact decision scope has duplicate alternatives; reassess the selected protein.";
            else if (recorded.Any(item => item.Decisions.Length > 1))
                blocker = "This exact proposal has conflicting recorded decisions.";
            else if (group.Key.Kind == PreparationChangeKind.Disulfide && group.Key.PartnerResidue is null)
                blocker = "The possible disulfide has no exact partner residue.";
            else if (approved.Length > 1)
                blocker = "More than one alternative was approved for this exact site.";
            else if (group.Key.Kind == PreparationChangeKind.HeavyAtom && explicitlyDeclined.Length > 0)
                blocker = "A required missing-atom repair was declined. Revise the selected protein to continue.";
            else if (group.Key.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.AlternateLocation &&
                     approved.Length == 0 && explicitlyDeclined.Length == proposals.Length)
                blocker = "Every alternative at this site was declined. Revise the selected protein to continue.";
            if (blocker is not null) blockers.Add(blocker);

            var satisfied = blocker is null && (approved.Length == 1 ||
                group.Key.Kind == PreparationChangeKind.Disulfide && explicitlyDeclined.Length == 1);
            var decisionStanding = blocker is not null ? "blocked" : satisfied ? "confirmed" : "pending";
            var options = ImmutableArray.CreateBuilder<PreparationDecisionOptionAccount>();
            foreach (var item in recorded)
            {
                var choice = item.Decisions.FirstOrDefault();
                var disposition = choice?.ChosenValue == ResearcherDecisionValue.Approved ? "confirmed" :
                    choice?.ChosenValue == ResearcherDecisionValue.Declined ? "declined" :
                    satisfied ? "notChosen" : "available";
                var evidenceReason = PreparationOptionPrerequisiteLocked(item.Proposal, report);
                var confirmationBlocker = disposition != "available" ? "This option is already decided or its site is complete." :
                    blocker ?? evidenceReason;
                options.Add(new PreparationDecisionOptionAccount(item.Proposal.Id, item.Proposal.ProposedChange,
                    disposition, choice?.Id, evidenceReason is null, confirmationBlocker,
                    false, false, report.Evidence.Where(evidence => evidence.SubjectId == item.Proposal.Id).ToImmutableArray(),
                    item.Proposal.Kind == PreparationChangeKind.ResidueState
                        ? PreparationInformationRole.ModelAssumption : PreparationInformationRole.Observed));
            }
            decisions.Add(new PreparationDecisionAccount(proposals[0].Id, group.Key.Kind,
                group.Key.Residue, group.Key.PartnerResidue, group.Key.AtomName, decisionStanding,
                approved.FirstOrDefault().Proposal?.Id ??
                    (group.Key.Kind is PreparationChangeKind.HeavyAtom or PreparationChangeKind.Disulfide
                        ? explicitlyDeclined.FirstOrDefault().Proposal?.Id : null),
                options.ToImmutable(), blocker));
        }
        var result = decisions.ToImmutable();
        // Final-confirmation disclosure depends on the complete obligation set, not iteration order.
        var remaining = result.Count(item => item.Standing != "confirmed");
        var recommendationMethodAvailable = SelectChemicalPolicyLocked(model, intended)?.Recommendation is not null;
        if (remaining != 1 || blockers.Count != 0 || recommendationMethodAvailable)
            result = result.Select(item => item with { Options = item.Options.Select(option =>
                option with { StartsPreparationOnConfirmation = false,
                    StartsPreparationOnDecline = false }).ToImmutableArray() }).ToImmutableArray();
        else
            result = result.Select(item => item with { Options = item.Options.Select(option =>
                option with { StartsPreparationOnConfirmation = item.Standing == "pending" &&
                    option.Disposition == "available",
                    StartsPreparationOnDecline = item.Standing == "pending" &&
                    item.Kind == PreparationChangeKind.Disulfide &&
                    option.Disposition == "available" }).ToImmutableArray() }).ToImmutableArray();
        var standing = _protein is not null ? "assessed" : _proteinPreparationRunning ? "preparing" :
            _proteinPreparationFailure is not null ? "failed" : blockers.Count > 0 ? "blocked" :
            remaining > 0 ? "awaitingDecisions" : "ready";
        return new ProteinPreparationReviewAccount(_study.Id, intended.Id, result,
            result.Length - remaining, remaining, blockers.ToImmutable(), standing,
            standing switch
            {
                "preparing" => "Applying the recorded choices and checking the resulting protein…",
                "failed" => _proteinPreparationFailure,
                "assessed" => "An assessed prepared protein was established from these recorded choices.",
                "blocked" => blockers.FirstOrDefault(),
                "ready" => "All choices are recorded. Protein preparation can start.",
                _ => null
            }, DecisionInspectionRelationLocked(intended));
    }

    private DecisionInspectionRelation DecisionInspectionRelationLocked(IntendedProteinModel intended)
    {
        var inspected = _inspection.Current;
        if (inspected is null) return DecisionInspectionRelation.Unavailable;
        if (inspected.StudyRevisionId != _study.Id) return DecisionInspectionRelation.Historical;
        if (inspected.SubjectId == intended.Id || inspected.SubjectId == intended.Source.Id ||
            _preparationProposals?.Changes.Any(change => change.Id == inspected.SubjectId) == true)
            return DecisionInspectionRelation.BeforePreparation;
        if (_protein is not null && _protein.Intended.Id == intended.Id &&
            inspected.SubjectId == _protein.Id) return DecisionInspectionRelation.PreparedResult;
        if (_proteinDiagnostic is not null && inspected.SubjectId == _proteinDiagnostic.Candidate.Id)
            return DecisionInspectionRelation.UnqualifiedCandidate;
        return DecisionInspectionRelation.OtherSubject;
    }

    private ProteinTaskAccount? BuildProteinTaskLocked(ProteinPreparationReviewAccount? review)
    {
        var intended = _study.IntendedProtein;
        if (intended is null) return null;
        var standing = _protein is not null ? "assessed" :
            _proteinPreparationRunning ? "preparing" :
            _proteinPreparationObservationIssue is not null ? "unavailable" :
            _proteinPreparationFailure is not null ? "failed" :
            _proteinSelectionRunning ? "assessing" :
            _proteinSelectionIssue is not null ? "blocked" :
            review?.Blockers.Length > 0 ? "blocked" :
            _preparationPlan is not null && _authorizedPreparationPlan is null ? "planReady" :
            review?.RemainingCount > 0 ? "awaitingDecisions" :
            review is not null ? "ready" : "unavailable";
        var message = standing switch
        {
            "assessed" => "The selected protein was prepared and passed its applicable preparation checks.",
            "preparing" => "Applying confirmed choices and checking the resulting protein…",
            "unavailable" => _proteinPreparationObservationIssue ??
                "The current protein preparation account is unavailable. Refresh the study before continuing.",
            "failed" => _proteinPreparationFailure,
            "assessing" => "Assessing the selected model and its required preparation choices…",
            "blocked" => _proteinSelectionIssue?.Message ?? review?.Blockers.FirstOrDefault(),
            "planReady" => "A jointly checked starting-state plan is ready for explicit authorization or review.",
            "awaitingDecisions" => $"{review!.RemainingCount} preparation decision{(review.RemainingCount == 1 ? "" : "s")} remain.",
            _ => "The selected protein is ready for preparation once its current prerequisites are confirmed."
        };
        return new ProteinTaskAccount(_study.Id, intended.Id, intended.Source.Id, intended.Chains,
            standing, message, _protein?.Id, standing == "failed" && _proteinPreparationRetryable &&
            _preparationProposals?.StudyRevisionId == _study.Id,
            _proteinSelectionIssue?.ReviewMoleculeSelection == true);
    }

    private PreparationPlanAccount? BuildPreparationPlanAccountLocked()
    {
        if (_study.IntendedProtein is null) return null;
        var plan = _preparationPlan;
        var partial = _preparationProposals is { } proposals &&
            (!proposals.UnresolvedQuestions.IsDefaultOrEmpty ||
             proposals.Changes.Any(change =>
                 (change.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.Disulfide) &&
                 !_decisions.Any(decision => decision.SubjectId == change.Id &&
                     (change.Kind == PreparationChangeKind.Disulfide ||
                      decision.ChosenValue == ResearcherDecisionValue.Approved))));
        var standing = _authorizedPreparationPlan is not null ? _protein is not null ? "applied" :
            _proteinPreparationRunning ? "preparing" : _proteinPreparationFailure is not null ? "failedAfterAuthorization" :
            "authorized" : _recommendationRunning ? "calculating" : plan is not null ? "ready" :
            partial ? "partial" : _recommendationIssue is not null ? "failed" : "unavailable";
        var current = plan ?? _authorizedPreparationPlan;
        var message = standing switch
        {
            "calculating" => "Finding preparation suggestions and checking one complete candidate…",
            "ready" when !_planOverrides.IsDefaultOrEmpty =>
                "Your changed choice and the remaining starting states were jointly checked. Review or prepare with this exact plan.",
            "ready" => "Starting-state suggestions are ready. Review them or explicitly prepare with this checked plan.",
            "partial" => "Review the exact conformer, possible bond, or structural exception before a complete plan can be offered.",
            "preparing" => "Applying the authorized plan and checking the resulting protein…",
            "applied" => "The authorized plan was applied and the resulting protein passed its preparation checks.",
            "failedAfterAuthorization" => _proteinPreparationFailure,
            "failed" => _recommendationIssue,
            _ => "No current starting-state plan is available. The manual review route remains available."
        };
        return new PreparationPlanAccount(standing, current?.Id, current?.PlanSha256,
            current is null ? null : $"{current.Method} · {current.MethodVersion}",
            current?.NominalPh, current?.Choices.Length ?? 0,
            current?.ProposedHeavyAtoms.Length ?? 0,
            current?.RemovedSourceHydrogens.Length ?? 0,
            current?.Choices ?? ImmutableArray<RecommendedStateChoice>.Empty,
            message, standing == "ready" && current?.PrerequisiteDecisionIds.SequenceEqual(
                _decisions.Select(decision => decision.Id).Order(StringComparer.Ordinal)) == true,
            !_planOverrides.IsDefaultOrEmpty);
    }

    private PlacementTaskAccount? BuildPlacementTaskLocked()
    {
        if (_protein is null || _study.Membrane is null) return null;
        var proposalId = _placementProposal?.Id;
        var adopted = proposalId is not null && _study.AdoptedPlacementProposalId == proposalId;
        var routeOutcome = _placementRouteOutcome is { } observed &&
            observed.StudyRevisionId == _study.Id && observed.PreparedProteinId == _protein.Id &&
            observed.MembraneModelId == _study.Membrane.Id ? observed : null;
        var routeDidNotEstablishPosition = routeOutcome?.Standing is
            "noMatch" or "noncorresponding" or "failed" or "unobserved";
        var standing = _placementRunning ? _placementAssessing ? "assessing" : "obtaining" :
            routeDidNotEstablishPosition ? "noProposal" :
            _placementProposal is null ? _placementOperationIssue is null ? "ready" : "noProposal" :
            _placement?.Standing switch
            {
                AssessmentStanding.Supported => "supported",
                AssessmentStanding.Unsupported => "unsupported",
                AssessmentStanding.NotEstablished => "notEstablished",
                _ => "review"
            };
        var message = standing switch
        {
            "obtaining" => "Positioning the complete prepared construct in the chosen membrane frame…",
            "assessing" => "Checking the complete positioned construct and membrane frame…",
            "ready" => _membrane is null ? "Complete the chosen membrane's parameter check before positioning." :
                "Choose a starting position, then move or rotate the complete construct.",
            "noProposal" => routeOutcome?.Message ??
                "No reviewable position was established by the latest attempt: " + _placementOperationIssue,
            "supported" when adopted => "Position selected. The same checked construct can continue to system preparation.",
            "supported" => "Position ready to use. Review it and explicitly choose Use this position.",
            "unsupported" => _placement?.Reason ?? "The current position failed a technical check.",
            "notEstablished" => _placement?.Reason ?? "The current position could not be checked.",
            _ => "Review the exact positioned construct before selecting it."
        };
        return new PlacementTaskAccount(_study.Id, _protein.Id, _study.Membrane.Id,
            proposalId, standing, message, adopted,
            _placementOperationIssue, routeOutcome);
    }

    private async Task<string?> AdoptMembraneAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var upper = ParseArray<LipidFraction>(data, "upper");
        var lower = ParseArray<LipidFraction>(data, "lower");
        if (!Coherent(upper) || !Coherent(lower))
            return "Both leaflets need coherent finite species fractions.";
        MembraneModel proposal;
        StudyRevision revision;
        MembraneSupportPolicy? policy;
        IReadOnlyDictionary<string, MolecularRepresentation> lipids;
        lock (_gate)
        {
            if (_membraneAssessmentRunning) return "This membrane support assessment is already running.";
            if (_study.Membrane is { } chosen &&
                SameComposition(chosen.Upper.Fractions, upper) &&
                SameComposition(chosen.Lower.Fractions, lower))
                return "This membrane composition is already selected. Review its current checks or use the available retry.";
            proposal = new MembraneModel(NewId(), new LeafletComposition(LeafletSide.Upper, upper),
                new LeafletComposition(LeafletSide.Lower, lower), _study.Conditions, null);
            AdvanceStudyLocked(_study.IntendedProtein, proposal, null);
            _allDecisions = _allDecisions.Add(new ResearcherDecision(NewId(), _study.Id,
                proposal.Id, ResearcherDecisionKind.AdoptMembrane, ResearcherDecisionValue.Adopted, DateTimeOffset.UtcNow,
                null));
            _protein = CarryProteinLocked(_study);
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementMeasurementIssue = null;
            _placementWitness = null;
            _placement = null;
            revision = _study;
            _membrane = null;
            _membraneAssessmentReason = null;
            _membraneAssessmentRunning = true;
            _membraneAssessmentUnavailable = false;
            policy = SelectMembranePolicyLocked(proposal);
            lipids = LipidsLocked();
            TouchLocked();
        }
        return await AssessChosenMembraneAsync(proposal, revision, policy, lipids, cancellationToken);
    }

    private async Task<string?> RetryMembraneCheckAsync(JsonElement data, CancellationToken cancellationToken)
    {
        MembraneModel proposal;
        StudyRevision revision;
        MembraneSupportPolicy? policy;
        IReadOnlyDictionary<string, MolecularRepresentation> lipids;
        lock (_gate)
        {
            if (_study.Membrane is not { } selected || selected.Id != Text(data, "modelId") ||
                !_membraneAssessmentUnavailable ||
                _membraneAssessmentRunning)
                return "Only a failed check of the exact currently selected membrane can be retried.";
            proposal = selected;
            revision = _study;
            _membrane = null;
            _membraneAssessmentRunning = true;
            _membraneAssessmentUnavailable = false;
            policy = SelectMembranePolicyLocked(proposal);
            lipids = LipidsLocked();
            TouchLocked();
        }
        return await AssessChosenMembraneAsync(proposal, revision, policy, lipids, cancellationToken);
    }

    private async Task<string?> AssessChosenMembraneAsync(MembraneModel proposal, StudyRevision revision,
        MembraneSupportPolicy? policy, IReadOnlyDictionary<string, MolecularRepresentation> lipids,
        CancellationToken cancellationToken)
    {
        Changed?.Invoke();
        // The adopted intention is available for inspection during its exact check.
        SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = proposal.Id }));
        BoundaryOutcome<AssessedMembraneModel> assessed;
        try
        {
            assessed = await _membraneAssessment.AssessAsync(revision, proposal, lipids, policy,
                WorkDirectory("membrane-assessment", revision.Id), cancellationToken);
        }
        catch (Exception exception)
        {
            lock (_gate)
            {
                if (_study.Id != revision.Id) return null;
                _membraneAssessmentRunning = false;
                _membraneAssessmentUnavailable = true;
                _membraneAssessmentReason = exception is OperationCanceledException
                    ? "The membrane support assessment outcome could not be observed after interruption. Refresh the account before retrying."
                    : $"Membrane support assessment was unavailable: {exception.Message}";
                TouchLocked();
            }
            Changed?.Invoke();
            return null;
        }
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _membraneAssessmentRunning = false;
            _membrane = assessed.Value;
            _membraneAssessmentReason = assessed.Value is null ? assessed.Reason : null;
            TouchLocked();
        }
        Changed?.Invoke();
        SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = proposal.Id }));
        return null;
    }

    private async Task<string?> ProposePlacementAsync(JsonElement data, CancellationToken cancellationToken)
    {
        StudyRevision revision;
        AssessedPreparedProtein protein;
        MembraneModel membrane;
        PlacementSupportPolicy framePolicy;
        IReadOnlyDictionary<string, MolecularRepresentation> frameSpecies;
        bool ppmAvailable;
        lock (_gate)
        {
            if (_protein is null || _study.Membrane is null)
                return "A corresponding assessed protein and chosen membrane model are required before placement.";
            revision = _study;
            protein = _protein;
            membrane = _study.Membrane;
            ppmAvailable = CanRunPpmLocked();
        }
        var route = Text(data, "orientationRoute") ?? "manual";
        if (route is not ("manual" or "ppm" or "opm"))
            return "Choose user-defined positioning or an available orientation method.";
        var startingPosition = Text(data, "startingPosition") switch
        {
            "center" => PlacementStartingPosition.Center,
            "upper" => PlacementStartingPosition.Upper,
            "lower" => PlacementStartingPosition.Lower,
            _ => (PlacementStartingPosition?)null
        };
        var manualValues = new[] { Number(data, "offsetXAngstrom"), Number(data, "offsetYAngstrom"),
            Number(data, "offsetZAngstrom"), Number(data, "rotationXDegrees"),
            Number(data, "rotationYDegrees"), Number(data, "rotationZDegrees") };
        if (route == "manual" && (startingPosition is null ||
            manualValues.Any(value => value is null || !double.IsFinite(value.Value))))
            return "Choose a starting position and finite movement and rotation values.";
        var topologyKind = ParseProteinTopology(Text(data, "topologyKind"));
        var physicalSide = ParsePlacementSide(Text(data, "physicalSide"));
        var ppmNterminalSide = ParsePpmNterminalSide(Text(data, "ppmNterminalSide"));
        if (route != "manual" && (topologyKind is null || physicalSide is null ||
            route == "ppm" && ppmNterminalSide is null)
            ) return "Choose the requested orientation method's topology, physical side and N-terminal assignment.";
        if (route == "opm" && (protein.Intended.Source.Kind != SourceRouteKind.Rcsb ||
            protein.Intended.Source.Accession is null))
            return "OPM lookup requires an identified RCSB entry; another available positioning method may be selected.";
        if (route == "ppm" && !ppmAvailable)
            return "The local PPM installation is not currently identified and hash-verified.";
        var routeRequestId = route == "manual" ? null : NewId();
        lock (_gate)
        {
            if (_placementRunning) return "An exact placement operation is already running.";
            var applicableFrames = ApplicablePlacementPoliciesLocked(
                route == "manual" ? ProteinTopologyKind.Unclassified : topologyKind!.Value, membrane);
            if (applicableFrames.Length != 1 ||
                applicableFrames[0].OuterLeafletEnvelopeAngstrom is not double envelope ||
                !double.IsFinite(envelope) || envelope <= 0)
                return "The chosen membrane has no unique declared placement frame for this position.";
            framePolicy = applicableFrames[0];
            frameSpecies = LipidsLocked();
            _placementRunning = true;
            _placementAssessing = false;
            _placementOperationIssue = null;
            _placementRouteOutcome = routeRequestId is null ? null :
                new PlacementRouteOutcomeAccount(routeRequestId, revision.Id, protein.Id, membrane.Id,
                    route, topologyKind!.Value, physicalSide!.Value, ppmNterminalSide,
                    route == "opm" ? "searching" : "running",
                    route == "opm" ? "Searching OPM for the exact selected source…" :
                        "Running the selected local PPM orientation…");
            TouchLocked();
        }
        Changed?.Invoke();
        try
        {
        BoundaryOutcome<PlacementProposal>? proposed = null;
        if (route == "manual")
            proposed = await _placementAssessment.ProposeManualAsync(revision, protein, membrane,
                frameSpecies, framePolicy,
                startingPosition!.Value, manualValues[0]!.Value, manualValues[1]!.Value,
                manualValues[2]!.Value, manualValues[3]!.Value, manualValues[4]!.Value,
                manualValues[5]!.Value,
                _catalogue?.MaximumSourceAtoms is > 0 ? _catalogue.MaximumSourceAtoms : int.MaxValue,
                WorkDirectory("placement", NewId()), cancellationToken);
        OpmReferenceReview? opm = null;
        string? routeStanding = null;
        string? routeMessage = null;
        if (route == "opm")
        {
            var lookup = await _sources.TryRetrieveOpmReferenceAsync(protein.Intended.Source.Accession!,
                WorkDirectory("opm-reference", revision.Id), cancellationToken);
            if (lookup.Standing != OpmLookupStanding.Found || lookup.Reference is null)
            {
                routeStanding = lookup.Standing switch
                {
                    OpmLookupStanding.NoMatch => "noMatch",
                    OpmLookupStanding.Failed => "failed",
                    _ => "unobserved"
                };
                routeMessage = routeStanding switch
                {
                    "noMatch" => lookup.Reason ?? "No matching OPM reference was found for this source.",
                    "failed" => lookup.Reason ?? "The OPM lookup failed.",
                    _ => lookup.Reason ?? "The OPM lookup did not establish a usable reference."
                };
                routeMessage += " An earlier position remains available if one was checked; retry or choose another method.";
            }
            else
            {
                var reviewed = _placementAssessment.ReviewOpmReference(revision, protein, membrane,
                    lookup.Reference);
                opm = reviewed.Value;
                if (opm is { CorrespondsToSelectedConstruct: false })
                {
                    routeStanding = "noncorresponding";
                    routeMessage = "An OPM reference was found, but it does not correspond to the selected construct. " +
                        string.Join(" ", opm.Limitations) +
                        " An earlier position remains available if one was checked.";
                }
                else if (opm is null)
                {
                    routeStanding = "unobserved";
                    routeMessage = (reviewed.Reason ?? "The OPM reference could not be checked against the selected construct.") +
                        " An earlier position remains available if one was checked.";
                }
                else
                {
                    proposed = _placementAssessment.ProposeFromOpmReference(revision, protein, membrane,
                        framePolicy, lookup.Reference, opm, topologyKind!.Value, physicalSide!.Value,
                        Text(data, "biologicalSidedness"));
                    if (proposed.Value is null)
                    {
                        routeStanding = "unobserved";
                        routeMessage = (proposed.Reason ?? "The OPM reference did not yield a usable mapped position for this construct.") +
                            " An earlier position remains available if one was checked.";
                    }
                }
            }
        }
        if (route == "ppm")
        {
            proposed = await _placementAssessment.ProposeWithPpmAsync(revision, protein, membrane,
                framePolicy, topologyKind!.Value, physicalSide!.Value,
                Text(data, "biologicalSidedness"), ppmNterminalSide!.Value,
                _ppmExecutablePath, _catalogue!.PpmVersion,
                _catalogue.PpmExecutableSha256,
                WorkDirectory("placement", NewId()), cancellationToken,
                _catalogue.PpmResidueLibraryPath, _catalogue.PpmResidueLibrarySha256);
            if (proposed.Value is null)
            {
                routeStanding = "unobserved";
                routeMessage = (proposed.Reason ?? "The local PPM route did not establish a usable position.") +
                    " An earlier position remains available if one was checked.";
            }
        }
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            if (proposed?.Value is not null)
            {
                _placementProposal = proposed.Value;
                _placementAssessing = true;
                _inspection.Clear();
                _opmReview = opm;
                _placementMeasurement = null;
                _placementMeasurementIssue = null;
                _placementWitness = null;
                _placement = null;
                if (routeRequestId is not null && _placementRouteOutcome?.RequestId == routeRequestId)
                    _placementRouteOutcome = _placementRouteOutcome with
                    {
                        Standing = "checking",
                        Message = "A mapped position was obtained; checking it against the selected membrane…",
                        ProposalId = proposed.Value.Id
                    };
            }
            else
            {
                _placementOperationIssue = proposed?.Reason ??
                    (routeStanding is "failed" or "unobserved" ? routeMessage : null);
                if (routeRequestId is not null && _placementRouteOutcome?.RequestId == routeRequestId)
                    _placementRouteOutcome = _placementRouteOutcome with
                    {
                        Standing = routeStanding ?? "unobserved",
                        Message = routeMessage ?? "No usable position was established by this request."
                    };
            }
            TouchLocked();
        }
        Changed?.Invoke();
        if (proposed?.Value is not null)
        {
            await AssessPlacementAsync(revision, proposed.Value, cancellationToken);
            SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = proposed.Value.Id }));
        }
        return null;
        }
        catch (Exception exception)
        {
            lock (_gate)
            {
                if (_study.Id != revision.Id) return null;
                _placementOperationIssue = exception is OperationCanceledException
                    ? "The placement outcome could not be observed after interruption. Refresh before a new operation."
                    : $"The requested position could not be obtained or checked: {exception.Message}";
                if (routeRequestId is not null && _placementRouteOutcome?.RequestId == routeRequestId)
                    _placementRouteOutcome = _placementRouteOutcome with
                    {
                        Standing = exception is OperationCanceledException ? "unobserved" : "failed",
                        Message = _placementOperationIssue
                    };
                TouchLocked();
            }
            Changed?.Invoke();
            return null;
        }
        finally
        {
            lock (_gate)
            {
                if (_study.Id == revision.Id)
                {
                    _placementRunning = false;
                    _placementAssessing = false;
                    TouchLocked();
                }
            }
            Changed?.Invoke();
        }
    }

    private async Task<string?> RevisePlacementAsync(JsonElement data, CancellationToken cancellationToken)
    {
        PlacementProposal source;
        StudyRevision revision;
        lock (_gate)
        {
            if (_placementProposal is null || _placementProposal.Id != Text(data, "proposalId"))
                return "Revise the current identified placement proposal.";
            source = _placementProposal;
            revision = _study;
        }
        var values = new[] { Number(data, "depthShiftAngstrom"), Number(data, "tiltAboutXDegrees"),
            Number(data, "tiltAboutYDegrees"), Number(data, "rotationAboutNormalDegrees") };
        if (values.Any(value => value is null || !double.IsFinite(value.Value)))
            return "A placement revision needs finite depth, tilt and rotation values.";
        lock (_gate)
        {
            if (_placementRunning) return "An exact placement operation is already running.";
            _placementRunning = true;
            _placementAssessing = false;
            _placementOperationIssue = null;
            _placementRouteOutcome = null;
            TouchLocked();
        }
        Changed?.Invoke();
        try
        {
        var revised = await _placementAssessment.ReviseProposalAsync(revision, source,
            values[0]!.Value, values[1]!.Value, values[2]!.Value, values[3]!.Value,
            WorkDirectory("placement", NewId()), cancellationToken);
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            if (revised.Value is not null)
            {
                _placementProposal = revised.Value;
                _placementAssessing = true;
                _inspection.Clear();
                _placementMeasurement = null;
                _placementMeasurementIssue = null;
                _placementWitness = null;
                _placement = null;
            }
            else
            {
                _placementOperationIssue = revised.Reason;
            }
            TouchLocked();
        }
        Changed?.Invoke();
        if (revised.Value is not null)
        {
            await AssessPlacementAsync(revision, revised.Value, cancellationToken);
            SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = revised.Value.Id }));
        }
        return null;
        }
        catch (Exception exception)
        {
            lock (_gate)
            {
                if (_study.Id != revision.Id) return null;
                _placementOperationIssue = exception is OperationCanceledException
                    ? "The corrected placement outcome could not be observed after interruption. Refresh before a new operation."
                    : $"The corrected position could not be obtained or checked: {exception.Message}";
                TouchLocked();
            }
            Changed?.Invoke();
            return null;
        }
        finally
        {
            lock (_gate)
            {
                if (_study.Id == revision.Id)
                {
                    _placementRunning = false;
                    _placementAssessing = false;
                    TouchLocked();
                }
            }
            Changed?.Invoke();
        }
    }

    private async Task AssessPlacementAsync(StudyRevision revision, PlacementProposal proposal, CancellationToken cancellationToken)
    {
        AssessedPreparedProtein? protein;
        MembraneModel? membrane;
        PlacementSupportPolicy? policy;
        PlacementStructuralWitness? witness;
        lock (_gate)
        {
            protein = _protein;
            membrane = _study.Membrane;
            policy = membrane is null ? null : SelectPlacementPolicyLocked(proposal, membrane);
            witness = protein is null || _membrane is null ? null :
                SelectPlacementWitnessLocked(protein, _membrane, proposal);
        }
        if (protein is null) return;
        PlacementMeasurementReport? measurement = null;
        string? measurementIssue = null;
        if (membrane is not null)
        {
            var observed = await _placementAssessment.MeasureAgainstMembraneAsync(revision, protein,
                membrane, proposal, policy, witness,
                WorkDirectory("placement-measurement", proposal.Id), cancellationToken);
            measurement = observed.Value;
            measurementIssue = observed.Value is null ? observed.Reason : null;
        }
        // PPM's candidate is not support for the explicit bilayer. Only the
        // child's policy-evaluated, corresponding measurement can add it.
        var assessed = _placementAssessment.Assess(revision, protein, membrane, proposal, policy,
            measurement, witness,
            measurement?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<ScientificFinding>.Empty);
        lock (_gate)
        {
            if (_study.Id != revision.Id || _placementProposal?.Id != proposal.Id) return;
            _placementMeasurement = measurement;
            _placementMeasurementIssue = measurementIssue;
            _placementWitness = witness;
            _placement = assessed;
            if (_placementRouteOutcome?.ProposalId == proposal.Id)
                _placementRouteOutcome = _placementRouteOutcome with
                {
                    Standing = assessed.Standing switch
                    {
                        AssessmentStanding.Supported => "supported",
                        AssessmentStanding.Unsupported => "unsupported",
                        _ => "notEstablished"
                    },
                    Message = assessed.Standing == AssessmentStanding.Supported
                        ? "The mapped position passed its current membrane checks. Select Use this position to continue."
                        : assessed.Reason
                };
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private string? AdoptPlacement(JsonElement data)
    {
        lock (_gate)
        {
            if (_placementProposal is null || _placement is null ||
                _placementProposal.Id != Text(data, "proposalId") ||
                _study.AdoptedPlacementProposalId == _placementProposal.Id ||
                _placement.Proposal.Id != _placementProposal.Id || _placement.Standing != AssessmentStanding.Supported ||
                _protein is null || _study.Membrane is null || _placementMeasurement is null)
                return "Only the current supported proposal can be adopted.";
            if (_placementMeasurement.Evidence.IsDefaultOrEmpty ||
                _placementMeasurement.Evidence.Any(item => item.SubjectId != _placementProposal.Id ||
                    string.IsNullOrWhiteSpace(item.Observation)))
                return "The exact placement has no attributable assessment evidence. Reassess its position.";
            if (!IsWorkspaceFile(_placementProposal.OrientedProtein.CoordinatePath) ||
                !VerifyHash(_placementProposal.OrientedProtein.CoordinatePath,
                    _placementProposal.OrientedProtein.CoordinateSha256))
                return "The oriented placement coordinates are missing or have changed. Reassess this proposal.";
            var revised = new StudyRevision(NewId(), _study.Number + 1, _study.IntendedProtein,
                _study.Membrane, _placementProposal.Id, _study.Conditions);
            var carriedProtein = CarryProteinLocked(revised);
            var carriedMembrane = CarryMembraneLocked(revised);
            if (carriedProtein is null || revised.Membrane is null)
                return "An affected assessed input cannot be carried into this changed study premise.";
            var policy = SelectPlacementPolicyLocked(_placementProposal, revised.Membrane);
            var reassessed = _placementAssessment.Assess(revised, carriedProtein, revised.Membrane,
                _placementProposal, policy, _placementMeasurement, _placementWitness,
                _placementMeasurement.Evidence,
                _placement.Findings);
            if (reassessed.Standing != AssessmentStanding.Supported)
                return "The adopted premise is not supported for the new study revision: " + reassessed.Reason;
            _study = revised;
            _revisions.Add(revised.Id, revised);
            _workspace.RetainStudy(revised);
            _inspection.Clear();
            _allDecisions = _allDecisions.Add(new ResearcherDecision(NewId(), revised.Id,
                _placementProposal.Id, ResearcherDecisionKind.AdoptPlacement, ResearcherDecisionValue.Adopted, DateTimeOffset.UtcNow,
                null));
            _protein = carriedProtein;
            _membrane = carriedMembrane;
            _placement = reassessed;
            _placementOperationIssue = null;
            TouchLocked();
        }
        Changed?.Invoke();
        SelectInspectionSubject(JsonSerializer.SerializeToElement(new { subjectId = _placementProposal!.Id }));
        return null;
    }

    private string? StartPreparation(JsonElement data)
    {
        string? refusal = null;
        lock (_gate)
        {
            if (_protein is null || _membrane is null || _placement is null ||
                _placement.Standing != AssessmentStanding.Supported ||
                _study.AdoptedPlacementProposalId != _placement.Proposal.Id ||
                _placement.StudyRevisionId != _study.Id || _protein.StudyRevisionId != _study.Id ||
                _membrane.StudyRevisionId != _study.Id)
                return "Current, corresponding supported inputs are required before preparation.";
            if (_attemptTask is { IsCompleted: false }) return "The current attempt is still running.";
            var requestedPolicyId = Text(data, "policyId");
            var policy = SelectPreparationPolicyLocked(_study, _protein, _placement.Proposal, _membrane,
                requestedPolicyId);
            if (policy is null) return "No applicable construction route and exact policy are available for these inputs.";
            if (!ConstructionProviderMatches(policy.Construction, _startupConstructionProvider) ||
                !PreparationAssetsAvailable(policy, _membrane))
                return "The selected construction provider or one of its exact assets is unavailable.";
            var admissionProvider = _constructionProviderProbe();
            _lastObservedConstructionProvider = admissionProvider;
            if (!ConstructionProviderMatches(policy.Construction, admissionProvider) ||
                !PreparationAssetsAvailable(policy, _membrane))
            {
                TouchLocked();
                refusal = "The selected construction provider or one of its exact assets changed before admission.";
            }
            else
            {
                var revision = _study;
                var protein = _protein;
                var membrane = _membrane;
                var placement = _placement;
                var decisionSubjects = protein.Changes.Select(change => change.Id)
                    .Append(membrane.Intended.Id).Append(placement.Proposal.Id)
                    .ToHashSet(StringComparer.Ordinal);
                var decisions = _allDecisions.Where(decision => decisionSubjects.Contains(decision.SubjectId))
                    .ToImmutableArray();
                var attemptStop = new CancellationTokenSource();
                if (_currentAttempt is { } earlier && _execution is { } earlierExecution)
                    _priorAttempts.Insert(0, AttemptAccountLocked(earlier, earlierExecution)!);
                _currentAttempt = null;
                _execution = new StageExecutionState(string.Empty, null, null, StageExecutionStanding.Pending,
                    "The corresponding attempt is being admitted.", null, DateTimeOffset.UtcNow);
                TrackAttemptTaskLocked(Task.Run(() => RunPreparationAsync(revision, protein, membrane, placement, policy,
                    decisions, attemptStop.Token)), attemptStop);
                TouchLocked();
            }
        }
        Changed?.Invoke();
        return refusal;
    }

    private static bool CandidateArtifactsMatch(ConstructedExplicitSystem candidate)
    {
        var molecule = candidate.Molecule;
        return VerifyHash(molecule.CoordinatePath, molecule.CoordinateSha256) &&
            molecule.TopologyPath is { } topologyPath &&
            molecule.TopologySha256 is { } topologySha && VerifyHash(topologyPath, topologySha) &&
            molecule.SystemXmlPath is { } systemPath &&
            molecule.SystemXmlSha256 is { } systemSha && VerifyHash(systemPath, systemSha) &&
            molecule.StateXmlPath is { } statePath &&
            molecule.StateXmlSha256 is { } stateSha && VerifyHash(statePath, stateSha) &&
            molecule.CorrespondencePath is { } mappingPath &&
            molecule.CorrespondenceSha256 is { } mappingSha && VerifyHash(mappingPath, mappingSha);
    }

    private string? StopAttempt(JsonElement data)
    {
        lock (_gate)
        {
            var id = Text(data, "attemptId");
            if (id is null || _currentAttempt?.Id != id || _execution?.AttemptId != id)
                return "No matching unfinished preparation attempt can be stopped.";
            if (_attemptTask is { IsCompleted: false } && _attemptStop is not null &&
                _execution.Standing is StageExecutionStanding.Pending or StageExecutionStanding.Running)
            {
                if (_stopRequestedAttemptId == id)
                    return "A stop request is already pending for this attempt; wait for its observed outcome.";
                _stopRequestedAttemptId = id;
                _attemptStop.Cancel();
            }
            else return "No matching unfinished preparation attempt can be stopped.";
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private string? RequestEquilibration(JsonElement data)
    {
        var id = Text(data, "stageId");
        lock (_gate)
        {
            if (id is null || !_stages.TryGetValue(id, out var minimized) || minimized.Kind != StageKind.Minimization)
                return "Select a completed minimized stage.";
            if (_attemptTask is { IsCompleted: false }) return "Another stage is already running.";
            if (!CanRequestEquilibrationLocked(minimized))
                return "No validated, applicable optional equilibration procedure is available for this stage.";
            var constructed = _constructedByAttempt[minimized.Attempt.Id];
            var policy = _policyByAttempt[minimized.Attempt.Id];
            var stop = new CancellationTokenSource();
            var pendingExecution = new StageExecutionState(minimized.Attempt.Id, null, StageKind.Equilibration,
                StageExecutionStanding.Pending, "Optional equilibration is starting.", null,
                DateTimeOffset.UtcNow);
            if (_workspace.Snapshot().CurrentAttempt?.Id != minimized.Attempt.Id)
                _workspace.RetainAttempt(minimized.Attempt);
            _workspace.RetainExecution(pendingExecution);
            _currentAttempt = minimized.Attempt;
            _execution = pendingExecution;
            TrackAttemptTaskLocked(Task.Run(() => RunEquilibrationAsync(minimized, constructed, policy,
                stop.Token)), stop);
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private void TrackAttemptTaskLocked(Task task, CancellationTokenSource stop)
    {
        _attemptStop = stop;
        _attemptTask = task;
        _stopRequestedAttemptId = null;
        _ = task.ContinueWith(completed =>
        {
            var notify = false;
            lock (_gate)
            {
                if (ReferenceEquals(_attemptTask, completed))
                {
                    _attemptStop = null;
                    _stopRequestedAttemptId = null;
                    TouchLocked();
                    notify = true;
                }
            }
            stop.Dispose();
            if (notify) Changed?.Invoke();
        }, CancellationToken.None, TaskContinuationOptions.None, TaskScheduler.Default);
    }

    private string? SelectInspectionSubject(JsonElement data)
    {
        var id = Text(data, "subjectId");
        if (id is null) return "Select an identified inspection subject.";
        InspectionSubject? subject;
        StudyRevision? revision;
        lock (_gate)
        {
            (subject, revision) = MakeInspectionSubjectLocked(id);
            if (subject is null || revision is null) return "The exact inspection subject is unavailable.";
            var selected = _inspection.Select(revision, subject);
            if (selected.Value is null) return selected.Reason;
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private string? SetInspectionFocus(JsonElement data)
    {
        var annotation = Text(data, "annotationId");
        lock (_gate)
        {
            var current = _inspection.Current;
            if (current is null) return "Select an inspection subject first.";
            var focused = annotation is null
                ? _inspection.ClearFocus(current.SubjectId)
                : _inspection.Focus(current.SubjectId, annotation);
            if (focused.Value is null) return focused.Reason;
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private async Task<string?> ExportStageAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var stageId = Text(data, "stageId");
        CompletedStage stage;
        StudyRevision revision;
        PreparationAssessmentResult assessment;
        ApplicablePreparationPolicy policy;
        AssessedPreparedProtein protein;
        AssessedMembraneModel membrane;
        AssessedProteinMembranePlacement placement;
        ConstructedExplicitSystem constructed;
        ConstructionDerivation derivation;
        ImmutableArray<ResearcherDecision> decisions;
        lock (_gate)
        {
            if (stageId is null || !_stages.TryGetValue(stageId, out var selected))
                return "Select an identified completed stage before export.";
            stage = selected;
            var gateReason = ExportGateReasonLocked(stage);
            if (gateReason is not null) return gateReason;
            assessment = _stageAssessments[stageId];
            revision = _revisions[stage.Attempt.StudyRevisionId];
            policy = _policyByAttempt[stage.Attempt.Id];
            protein = _proteinByAttempt[stage.Attempt.Id];
            membrane = _membraneByAttempt[stage.Attempt.Id];
            placement = _placementByAttempt[stage.Attempt.Id];
            constructed = _constructedByAttempt[stage.Attempt.Id];
            derivation = _derivationByAttempt[stage.Attempt.Id];
            decisions = _decisionsByAttempt[stage.Attempt.Id];
            if (_bundles.TryGetValue(stageId, out var existing) &&
                existing.AssessmentId == assessment.Id && IsWorkspaceFile(existing.BundlePath) &&
                new FileInfo(existing.BundlePath).Length == existing.ByteLength &&
                VerifyHash(existing.BundlePath, existing.Sha256)) return null;
            _bundles.Remove(stageId);
        }
        BoundaryOutcome<CompletedStageBundle> delivered;
        try
        {
            var directory = WorkDirectory("exports", stage.Id + "-" + assessment.Id);
            string? optionalProtocolSha256 = null;
            if (stage.Kind == StageKind.Equilibration && policy.OptionalEquilibration is { } protocol &&
                EquilibrationProtocolFingerprint.TryCompute(protocol, out var digest))
                optionalProtocolSha256 = digest;
            delivered = await _export.ExportAsync(stage, assessment, revision, protein,
                membrane, placement, constructed, derivation, policy, decisions,
                _catalogue?.PpmVersion, _catalogue?.PpmExecutableSha256,
                optionalProtocolSha256,
                directory, cancellationToken);
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException)
        {
            delivered = BoundaryOutcome<CompletedStageBundle>.Unavailable(
                $"The completed-stage export destination could not be used: {exception.Message}");
        }
        string? outcomeReason;
        lock (_gate)
        {
            if (!_stageAssessments.TryGetValue(stage.Id, out var current) || current.Id != assessment.Id ||
                !current.CurrentlyApplicable)
            {
                outcomeReason = "The stage assessment changed during export; review the current assessment before delivery.";
                if (current is not null)
                    _exportAccounts[stage.Id] = new StageExportAccount(stage.Id, current.Id,
                        "failed", outcomeReason, null, null);
            }
            else if (delivered.Value is { } bundle)
            {
                _bundles[stage.Id] = bundle;
                _exportAccounts[stage.Id] = new StageExportAccount(stage.Id, assessment.Id,
                    "verified", null, bundle.Sha256, bundle.ByteLength);
                outcomeReason = null;
            }
            else
            {
                _exportAccounts[stage.Id] = new StageExportAccount(stage.Id, assessment.Id,
                    "failed", delivered.Reason, null, null);
                outcomeReason = delivered.Reason;
            }
            TouchLocked();
        }
        Changed?.Invoke();
        return outcomeReason;
    }

    private string? ExportGateReasonLocked(CompletedStage stage)
    {
        if (stage.Kind is not (StageKind.Minimization or StageKind.Equilibration) ||
            stage.Observation.Kind != stage.Kind ||
            stage.Observation.StageId != stage.Id || stage.Observation.AttemptId != stage.Attempt.Id ||
            stage.Molecule.Id != stage.Id || stage.Correspondence.ResultId != stage.Id ||
            !stage.Correspondence.Complete ||
            stage.Kind == StageKind.Minimization &&
                (stage.SourceStageId is not null || stage.Observation.Termination != StageTermination.Converged) ||
            stage.Kind == StageKind.Equilibration &&
                (stage.Observation.Termination != StageTermination.Completed ||
                 stage.Observation.ObservationAdequacy is null ||
                 stage.Observation.EquilibrationWindows.IsDefaultOrEmpty ||
                 stage.Observation.EquilibrationSamples.IsDefaultOrEmpty ||
                 stage.Observation.EquilibrationAssessments.IsDefaultOrEmpty))
            return "Select a factually completed, observed molecular stage before export.";
        if (!_stageAssessments.TryGetValue(stage.Id, out var assessment) ||
            !assessment.CurrentlyApplicable || assessment.StageId != stage.Id)
            return "The selected completed stage has no currently applicable assessment.";
        if (!_revisions.ContainsKey(stage.Attempt.StudyRevisionId) ||
            !_policyByAttempt.TryGetValue(stage.Attempt.Id, out var policy) ||
            !_proteinByAttempt.ContainsKey(stage.Attempt.Id) ||
            !_membraneByAttempt.ContainsKey(stage.Attempt.Id) ||
            !_placementByAttempt.ContainsKey(stage.Attempt.Id) ||
            !_constructedByAttempt.ContainsKey(stage.Attempt.Id) ||
            !_derivationByAttempt.ContainsKey(stage.Attempt.Id) ||
            !_decisionsByAttempt.ContainsKey(stage.Attempt.Id))
            return "The originating study, attempt and molecular preparation context are unavailable.";
        if (stage.PolicyId != policy.Id || !PreparationPolicyFingerprint.Matches(stage.Attempt, policy))
            return "The completed stage no longer matches its originating preparation policy.";
        if (stage.Kind == StageKind.Equilibration &&
            (string.IsNullOrWhiteSpace(stage.SourceStageId) ||
             !_stages.TryGetValue(stage.SourceStageId, out var source) ||
             source.Kind != StageKind.Minimization ||
             source.Observation.Kind != StageKind.Minimization ||
             source.Observation.StageId != source.Id ||
             source.Observation.Termination != StageTermination.Converged ||
             source.Attempt.Id != stage.Attempt.Id ||
             source.Attempt.StudyRevisionId != stage.Attempt.StudyRevisionId ||
             source.PolicyId != stage.PolicyId ||
             source.Molecule.AtomCount != stage.Molecule.AtomCount ||
             source.Correspondence.SourceId != stage.Correspondence.SourceId ||
             !source.Correspondence.Complete ||
             source.Correspondence.Atoms.IsDefault || stage.Correspondence.Atoms.IsDefault ||
             !source.Correspondence.Atoms.SequenceEqual(stage.Correspondence.Atoms) ||
             policy.OptionalEquilibration is null ||
             !EquilibrationProtocolFingerprint.TryCompute(policy.OptionalEquilibration, out _)))
            return "The equilibrated stage has no corresponding completed minimized source and bound procedure.";
        return null;
    }

    private async Task RunPreparationAsync(StudyRevision revision, AssessedPreparedProtein protein,
        AssessedMembraneModel membrane, AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy, ImmutableArray<ResearcherDecision> decisions,
        CancellationToken cancellationToken)
    {
        var attemptId = NewId();
        var admitted = false;
        try
        {
            var attempt = new PreparationAttempt(attemptId, revision.Id, protein.Id, membrane.Id,
                placement.Id, policy.Id, DateTimeOffset.UtcNow, policy.Version,
                PreparationPolicyFingerprint.Compute(policy), policy.ForceFieldFiles,
                policy.Construction.ProviderVersion, policy.Construction.NativePatchSha256,
                policy.Construction.Route, policy.Construction.SaltConvention,
                policy.Construction.ProviderAssets.IsDefault
                    ? ImmutableArray<ProviderAsset>.Empty : policy.Construction.ProviderAssets);
            var started = await _explicitPreparation.StartAsync(attempt, revision, protein, membrane, placement, policy,
                WorkDirectory("attempts", attempt.Id), accepted =>
                {
                    lock (_gate)
                    {
                        _workspace.RetainAttempt(accepted);
                        _currentAttempt = accepted;
                        admitted = true;
                        _execution = new StageExecutionState(accepted.Id, null, null,
                            StageExecutionStanding.Pending, "The identified attempt was admitted.",
                            null, DateTimeOffset.UtcNow);
                        TouchLocked();
                    }
                    Changed?.Invoke();
                }, Progress(), cancellationToken);
            if (started.Attempt is null)
            {
                SetExecution(started.State);
                return;
            }
            lock (_gate)
            {
                var outcome = PreserveObservedConstructionPhase(_execution, started.State);
                _workspace.RetainExecution(outcome);
                _currentAttempt = started.Attempt;
                _policyByAttempt[started.Attempt.Id] = policy;
                _proteinByAttempt[started.Attempt.Id] = protein;
                _membraneByAttempt[started.Attempt.Id] = membrane;
                _placementByAttempt[started.Attempt.Id] = placement;
                _decisionsByAttempt[started.Attempt.Id] = decisions;
                if (started.Derivation is not null)
                    _derivationByAttempt[started.Attempt.Id] = started.Derivation;
                if (!started.Trials.IsDefault)
                    _trialsByAttempt[started.Attempt.Id] = started.Trials
                        .Select(trial => BindTrialDiagnosticsLocked(started.Attempt.Id, trial))
                        .ToImmutableArray();
                if (started.Constructed is not null)
                    _constructedByAttempt[started.Attempt.Id] = started.Constructed;
                if (_stopRequestedAttemptId == started.Attempt.Id && outcome.Standing is
                    StageExecutionStanding.Completed or StageExecutionStanding.Stopped or
                    StageExecutionStanding.Failed or StageExecutionStanding.ResourceRefused or
                    StageExecutionStanding.Unobserved)
                    _stopRequestedAttemptId = null;
                _execution = outcome;
                TouchLocked();
            }
            Changed?.Invoke();
            if (started.Constructed is null) return;
            if (started.State.Standing is StageExecutionStanding.Failed or StageExecutionStanding.Stopped or
                StageExecutionStanding.ResourceRefused or StageExecutionStanding.Unobserved)
            {
                SetExecution(started.State with
                {
                    Standing = StageExecutionStanding.Unobserved,
                    Message = "The construction handoff had a terminal standing despite supplying a candidate."
                });
                return;
            }
            if (!CandidateArtifactsMatch(started.Constructed))
            {
                SetExecution(new StageExecutionState(attemptId, null, null,
                    StageExecutionStanding.Failed,
                    "The checked construction handoff changed before required minimization could start.",
                    null, DateTimeOffset.UtcNow, PreparationPhase.HandoffChecks,
                    FailureCode: "handoffArtifactChanged"));
                return;
            }
            await RunMinimizationAsync(started.Constructed, protein, membrane, placement, policy,
                cancellationToken);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            SetExecution(new StageExecutionState(admitted ? attemptId : string.Empty, null, null,
                StageExecutionStanding.Stopped,
                admitted ? "The unfinished preparation work was stopped." :
                    "Preparation was stopped before an attempt was admitted.",
                null, DateTimeOffset.UtcNow));
        }
        catch (Exception exception)
        {
            var resourceRefused = !admitted && exception is IOException or UnauthorizedAccessException;
            SetExecution(new StageExecutionState(admitted ? attemptId : string.Empty, null, null,
                resourceRefused ? StageExecutionStanding.ResourceRefused :
                    admitted ? StageExecutionStanding.Unobserved : StageExecutionStanding.Failed,
                (resourceRefused ? "The local workspace could not admit preparation: " :
                    admitted ? "The local worker outcome could not be established: " :
                        "Preparation could not be admitted: ") + exception.Message,
                null, DateTimeOffset.UtcNow,
                FailureCode: resourceRefused ? "resourceRefused" :
                    admitted ? "unobservedWorkerOutcome" : "admissionFailed"));
        }
    }

    private async Task RunMinimizationAsync(ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy,
        CancellationToken cancellationToken)
    {
        try
        {
            var minimized = await _explicitPreparation.MinimizeAsync(constructed, policy,
                WorkDirectory("attempts", constructed.Attempt.Id), Progress(), cancellationToken);
            if (minimized.CompletedStage is not null)
            {
                RetainCompleted(minimized.CompletedStage, constructed, protein, membrane, placement, policy);
                SetExecution(minimized.State);
            }
            else SetExecution(minimized.State.Standing == StageExecutionStanding.Completed
                ? minimized.State with
                {
                    Standing = StageExecutionStanding.Unobserved,
                    Message = "The stage reported completion without a retained completed result."
                }
                : minimized.State);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            SetExecution(new StageExecutionState(constructed.Attempt.Id, null, StageKind.Minimization,
                StageExecutionStanding.Stopped, "Required minimization was stopped; the constructed candidate is retained.",
                null, DateTimeOffset.UtcNow));
        }
        catch (Exception exception)
        {
            SetExecution(new StageExecutionState(constructed.Attempt.Id, null, StageKind.Minimization,
                StageExecutionStanding.Unobserved, "The minimization outcome could not be established: " + exception.Message,
                null, DateTimeOffset.UtcNow));
        }
    }

    private async Task RunEquilibrationAsync(CompletedStage minimized, ConstructedExplicitSystem constructed,
        ApplicablePreparationPolicy policy, CancellationToken cancellationToken)
    {
        try
        {
            var result = await _explicitPreparation.RequestOptionalEquilibrationAsync(minimized, policy,
                WorkDirectory("attempts", minimized.Attempt.Id), Progress(), cancellationToken);
            if (result.CompletedStage is not null)
            {
                AssessedPreparedProtein? protein;
                AssessedMembraneModel? membrane;
                AssessedProteinMembranePlacement? placement;
                lock (_gate)
                {
                    protein = _proteinByAttempt.GetValueOrDefault(minimized.Attempt.Id);
                    membrane = _membraneByAttempt.GetValueOrDefault(minimized.Attempt.Id);
                    placement = _placementByAttempt.GetValueOrDefault(minimized.Attempt.Id);
                }
                if (protein is not null && membrane is not null && placement is not null)
                {
                    RetainCompleted(result.CompletedStage, constructed, protein, membrane, placement, policy);
                    SetExecution(result.State);
                }
                else SetExecution(result.State with
                {
                    Standing = StageExecutionStanding.Unobserved,
                    Message = "The completed optional stage lacks its retained preparation context."
                });
            }
            else SetExecution(result.State.Standing == StageExecutionStanding.Completed
                ? result.State with
                {
                    Standing = StageExecutionStanding.Unobserved,
                    Message = "The optional stage reported completion without a retained completed result."
                }
                : result.State);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            SetExecution(new StageExecutionState(minimized.Attempt.Id, null, StageKind.Equilibration,
                StageExecutionStanding.Stopped, "Optional equilibration stopped; the minimized stage is retained.", null, DateTimeOffset.UtcNow));
        }
        catch (Exception exception)
        {
            SetExecution(new StageExecutionState(minimized.Attempt.Id, null, StageKind.Equilibration,
                StageExecutionStanding.Unobserved, "Optional equilibration outcome is unobserved: " + exception.Message,
                null, DateTimeOffset.UtcNow));
        }
    }

    private IProgress<StageExecutionState> Progress() => new DirectProgress(SetExecution);

    private void ObserveConstructionProgress(string requestId, JsonElement payload)
    {
        if (Text(payload, "operation") != "construct_system" ||
            Text(payload, "attemptId") is not { Length: > 0 } attemptId ||
            Text(payload, "studyRevisionId") is not { Length: > 0 } revisionId ||
            Text(payload, "stage") is not { Length: > 0 } stage)
            return;
        var phase = stage switch
        {
            "providerPopulation" => PreparationPhase.ProviderPopulation,
            "providerPacking" => PreparationPhase.ProviderPacking,
            "providerCleanup" => PreparationPhase.ProviderCleanup,
            "amberParameterization" => PreparationPhase.AmberParameterization,
            "providerConditioningRestrained" => PreparationPhase.ProviderConditioningRestrained,
            "providerConditioningUnrestrained" => PreparationPhase.ProviderConditioningUnrestrained,
            _ => (PreparationPhase?)null
        };
        if (phase is null) return;
        var detail = payload.TryGetProperty("detail", out var value) && value.ValueKind == JsonValueKind.Object
            ? value : default;
        var trialId = detail.ValueKind == JsonValueKind.Object ? Text(detail, "trialId") : null;
        int? trialIndex = detail.ValueKind == JsonValueKind.Object &&
            detail.TryGetProperty("trialIndex", out var indexValue) && indexValue.TryGetInt32(out var index)
            ? index : null;
        double? fraction = detail.ValueKind == JsonValueKind.Object &&
            detail.TryGetProperty("progress", out var progressValue) && progressValue.TryGetDouble(out var progress) &&
            double.IsFinite(progress) && progress is >= 0 and <= 1 ? progress : null;
        var message = detail.ValueKind == JsonValueKind.Object ? Text(detail, "message") : null;
        if (string.IsNullOrWhiteSpace(trialId) || trialIndex is null) return;
        lock (_gate)
        {
            if (_currentAttempt?.Id != attemptId || _currentAttempt.StudyRevisionId != revisionId ||
                _execution?.Standing is not (StageExecutionStanding.Pending or StageExecutionStanding.Running) ||
                _execution.TrialId != trialId || _execution.TrialIndex != trialIndex ||
                !SetExecutionLocked(new StageExecutionState(attemptId, null, null,
                    StageExecutionStanding.Running, message ?? stage, fraction,
                    DateTimeOffset.UtcNow, phase, trialId, trialIndex)))
                return;
        }
        Changed?.Invoke();
    }

    private sealed class DirectProgress(Action<StageExecutionState> receive) : IProgress<StageExecutionState>
    {
        public void Report(StageExecutionState value) => receive(value);
    }

    private static StageExecutionState PreserveObservedConstructionPhase(
        StageExecutionState? previous, StageExecutionState incoming)
    {
        // The construction owner knows that a trial ended without a checked
        // handoff, while only the correlated worker progress identifies its
        // last observed provider phase. A terminal fallback must not move the
        // same trial back to population after packing or conditioning began.
        if (previous is not null && previous.AttemptId == incoming.AttemptId &&
            previous.TrialId is { Length: > 0 } trialId && trialId == incoming.TrialId &&
            previous.Kind is null && incoming.Kind is null &&
            previous.Phase is { } observed &&
            observed >= PreparationPhase.ProviderPacking &&
            observed <= PreparationPhase.HandoffChecks &&
            incoming.Phase is null or PreparationPhase.ProviderPopulation)
            return incoming with { Phase = observed };
        return incoming;
    }

    private void SetExecution(StageExecutionState state)
    {
        lock (_gate)
            if (!SetExecutionLocked(state)) return;
        Changed?.Invoke();
    }

    private bool SetExecutionLocked(StageExecutionState state)
    {
        if (_currentAttempt is { } attempt && state.AttemptId != attempt.Id)
            return false;
        if (_currentAttempt is null && _execution?.AttemptId is { Length: > 0 } pendingId &&
            state.AttemptId != pendingId)
            return false;
        if (_execution is { } previous && previous.AttemptId == state.AttemptId)
        {
            var terminalRegression = previous.Standing is StageExecutionStanding.Completed or
                StageExecutionStanding.Stopped or StageExecutionStanding.Failed or
                StageExecutionStanding.ResourceRefused or StageExecutionStanding.Unobserved &&
                state.Standing is StageExecutionStanding.Pending or StageExecutionStanding.Running;
            if (previous.Kind is not null && previous.Kind != state.Kind ||
                previous.StageId is not null && previous.StageId != state.StageId ||
                terminalRegression)
                return false;
        }
        state = PreserveObservedConstructionPhase(_execution, state);
        // The child reports terminal progress before the root has judged and
        // retained its completed stage. Only the paired stage and assessment
        // may make completion visible to a reattached actor.
        if (state.Standing == StageExecutionStanding.Completed &&
            (state.StageId is null || !_stages.TryGetValue(state.StageId, out var stage) ||
             stage.Attempt.Id != state.AttemptId || !_stageAssessments.ContainsKey(state.StageId)))
            return false;
        if (_currentAttempt?.Id == state.AttemptId)
        {
            if (state.Trial is { } trial)
            {
                trial = BindTrialDiagnosticsLocked(state.AttemptId, trial);
                state = state with { Trial = trial };
                var known = _trialsByAttempt.GetValueOrDefault(state.AttemptId);
                if (known.IsDefault) known = ImmutableArray<ConstructionTrialSummary>.Empty;
                var position = -1;
                for (var index = 0; index < known.Length; index++)
                    if (known[index].TrialId == trial.TrialId) { position = index; break; }
                _trialsByAttempt[state.AttemptId] = position < 0 ? known.Add(trial) :
                    known.SetItem(position, trial);
            }
            _workspace.RetainExecution(state);
        }
        if (_stopRequestedAttemptId == state.AttemptId && state.Standing is
            StageExecutionStanding.Completed or StageExecutionStanding.Stopped or
            StageExecutionStanding.Failed or StageExecutionStanding.ResourceRefused or
            StageExecutionStanding.Unobserved)
            _stopRequestedAttemptId = null;
        _execution = state;
        TouchLocked();
        return true;
    }

    private void RetainCompleted(CompletedStage stage, ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy)
    {
        lock (_gate)
        {
            if (_currentAttempt?.Id != stage.Attempt.Id)
                throw new InvalidOperationException("A completed stage must identify the surviving attempt.");
            var judgedPlacement = placement;
            var newPlacementFindings = stage.Findings.Where(finding => finding.Material &&
                finding.SubjectId == placement.Proposal.Id &&
                (finding.Disposition is FindingDisposition.Challenges or FindingDisposition.Disqualifies))
                .ToImmutableArray();
            var reassessPlacement = !newPlacementFindings.IsEmpty &&
                _study.Id == stage.Attempt.StudyRevisionId &&
                _placement?.Id == placement.Id && _protein?.Id == protein.Id &&
                _membrane?.Id == membrane.Id && _placementMeasurement?.ProposalId == placement.Proposal.Id;
            if (reassessPlacement)
            {
                judgedPlacement = _placementAssessment.Assess(_study, _protein!, _membrane!.Intended,
                    placement.Proposal, SelectPlacementPolicyLocked(placement.Proposal, _membrane!.Intended),
                    _placementMeasurement, _placementWitness, _placementMeasurement!.Evidence,
                    placement.Findings.AddRange(newPlacementFindings));
            }
            // Stage.Findings already belongs to this stage. The additional argument is
            // reserved for genuinely later observations, not a duplicate of its own.
            var assessment = _preparationAssessment.Assess(stage, constructed, protein, membrane,
                judgedPlacement, policy, ImmutableArray<ScientificFinding>.Empty);
            var consequentialLaterFindings = stage.Findings.Where(item => item.Material &&
                (item.Disposition is FindingDisposition.Challenges or FindingDisposition.Disqualifies))
                .ToImmutableArray();
            var renewedEarlier = _stages.Values
                .Where(item => item.Id != stage.Id && item.Attempt.Id == stage.Attempt.Id)
                .SelectMany(earlier =>
                {
                    var relevant = consequentialLaterFindings.Where(finding =>
                        finding.SubjectId == earlier.Id ||
                        finding.SubjectId == earlier.Attempt.Id ||
                        finding.SubjectId == constructed.Id ||
                        finding.SubjectId == protein.Id ||
                        finding.SubjectId == membrane.Id ||
                        finding.SubjectId == membrane.Intended.Id ||
                        finding.SubjectId == placement.Id ||
                        finding.SubjectId == placement.Proposal.Id).ToImmutableArray();
                    if (relevant.IsEmpty)
                        return Array.Empty<(string Id, ImmutableArray<ScientificFinding> Accumulated,
                            PreparationAssessmentResult Renewed)>();
                    var accumulated = _laterFindingsByStage.GetValueOrDefault(earlier.Id);
                    if (accumulated.IsDefault) accumulated = ImmutableArray<ScientificFinding>.Empty;
                    accumulated = accumulated.AddRange(relevant.Where(item =>
                        accumulated.All(previous => previous.Id != item.Id)));
                    var renewed = _preparationAssessment.Assess(earlier, constructed, protein,
                        membrane, judgedPlacement, policy, accumulated);
                    return new[] { (earlier.Id, accumulated, renewed) };
                }).ToArray();

            // Validate and retain the factual result before any root-facing stage,
            // assessment, placement or finding is published. A rejected workspace
            // handoff leaves the previously visible account unchanged.
            _workspace.RetainCompletedStage(stage);
            _workspace.RetainAssessment(assessment);
            foreach (var (_, _, renewed) in renewedEarlier)
                _workspace.RetainAssessment(renewed);
            foreach (var finding in stage.Findings) _workspace.RetainFinding(finding);
            var completedState = new StageExecutionState(stage.Attempt.Id, stage.Id, stage.Kind,
                StageExecutionStanding.Completed, $"{stage.Kind} completed.", 1.0, stage.CompletedAt);
            _workspace.RetainExecution(completedState);

            if (reassessPlacement)
            {
                _placement = judgedPlacement;
            }
            _stages[stage.Id] = stage;
            _stageAssessments[stage.Id] = assessment;
            foreach (var (earlierId, accumulated, renewed) in renewedEarlier)
            {
                _laterFindingsByStage[earlierId] = accumulated;
                _stageAssessments[earlierId] = renewed;
                _bundles.Remove(earlierId);
                _exportAccounts.Remove(earlierId);
            }
            if (_stopRequestedAttemptId == stage.Attempt.Id)
                _stopRequestedAttemptId = null;
            _execution = completedState;
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private (InspectionSubject?, StudyRevision?) MakeInspectionSubjectLocked(string id)
    {
        if (_selectedSource is not null && _sourceInspection is not null && id == _selectedSource.Id &&
            IsWorkspaceFile(_selectedSource.CoordinatePath) &&
            VerifyHash(_selectedSource.CoordinatePath, _selectedSource.Sha256))
            return (GenericSubject(id, _study.Id, _selectedSource.CoordinatePath, "structuralSource",
                ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty, null), _study);
        if (_study.IntendedProtein is not null && id == _study.IntendedProtein.Id &&
            IsWorkspaceFile(_study.IntendedProtein.Source.CoordinatePath) &&
            VerifyHash(_study.IntendedProtein.Source.CoordinatePath,
                _study.IntendedProtein.Source.Sha256))
            return (WithGeometry(GenericSubject(id, _study.Id,
                _preparationProposals?.Preview is { } selectedPreview &&
                IsWorkspaceFile(selectedPreview.Path) && VerifyHash(selectedPreview.Path, selectedPreview.Sha256)
                    ? selectedPreview.Path : _study.IntendedProtein.Source.CoordinatePath,
                "intendedProtein", ImmutableArray<ScientificEvidence>.Empty,
                ImmutableArray<ScientificFinding>.Empty, null), _preparationProposals?.SourceGeometry), _study);
        if (_preparationProposals?.Changes.Any(item => item.Id == id) == true)
        {
            var preview = _preparationProposals.Preview;
            if (preview is null || !IsWorkspaceFile(preview.Path) ||
                !VerifyHash(preview.Path, preview.Sha256)) return (null, null);
            var url = StructureUrlLocked(preview.Path);
            var subject = _proteinPreparation.ProposalInspectionSubject(_preparationProposals, id, url);
            return (subject.Value, _study);
        }
        if (_proteinDiagnostic is { } diagnostic && id == diagnostic.Candidate.Id &&
            IsWorkspaceFile(diagnostic.Candidate.CoordinatePath) &&
            VerifyHash(diagnostic.Candidate.CoordinatePath, diagnostic.Candidate.CoordinateSha256))
            return (WithGeometry(GenericSubject(id, _study.Id, diagnostic.Candidate.CoordinatePath,
                "unqualifiedProteinCandidate", diagnostic.Evidence, diagnostic.Findings, null),
                diagnostic.Geometry), _study);
        if (_protein is not null && id == _protein.Id &&
            IsWorkspaceFile(_protein.Molecule.CoordinatePath) &&
            VerifyHash(_protein.Molecule.CoordinatePath, _protein.Molecule.CoordinateSha256))
            return (WithGeometry(GenericSubject(id, _study.Id, _protein.Molecule.CoordinatePath, "preparedProtein",
                _protein.Evidence.AddRange(_opmReview?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty),
                _protein.Findings, null), _protein.Geometry), _study);
        if (_study.Membrane is { } selectedMembrane && id == selectedMembrane.Id)
        {
            var evidence = _membrane?.Intended.Id == id
                ? _membrane.Evidence
                : _study.Membrane?.Id == id && !string.IsNullOrWhiteSpace(_membraneAssessmentReason)
                    ? ImmutableArray.Create(new ScientificEvidence($"membrane-assessment-{_study.Id}", id,
                        "Membrane Model Assessment",
                        "Membrane-local assessment", _membraneAssessmentReason,
                        $"Chosen membrane model {id}; study revision {_study.Id}",
                        "No assessed membrane model was established for this choice.", EvidenceBearing.Unknown))
                    : ImmutableArray<ScientificEvidence>.Empty;
            return (GenericSubject(id, _study.Id, null, "membraneModel", evidence,
                ImmutableArray<ScientificFinding>.Empty, null), _study);
        }
        if (_placementProposal is not null && id == _placementProposal.Id)
        {
            if (!IsWorkspaceFile(_placementProposal.OrientedProtein.CoordinatePath) ||
                !VerifyHash(_placementProposal.OrientedProtein.CoordinatePath,
                    _placementProposal.OrientedProtein.CoordinateSha256))
                return (null, null);
            var url = StructureUrlLocked(_placementProposal.OrientedProtein.CoordinatePath);
            var subject = _placementAssessment.ProposalInspectionSubject(_study, _placementProposal,
                _placementMeasurement, url, _placementMeasurementIssue);
            return (subject.Value, _study);
        }
        if (_diagnostics.TryGetValue(id, out var diagnosticBinding) &&
            _trialsByAttempt.TryGetValue(diagnosticBinding.AttemptId, out var diagnosticTrials) &&
            !diagnosticTrials.IsDefaultOrEmpty &&
            _revisions.TryGetValue(diagnosticBinding.StudyRevisionId, out var diagnosticOrigin))
        {
            var artifact = diagnosticTrials.Where(trial => trial.TrialId == diagnosticBinding.TrialId)
                .SelectMany(trial => trial.DiagnosticArtifacts.IsDefault
                    ? ImmutableArray<TrialDiagnosticArtifact>.Empty : trial.DiagnosticArtifacts)
                .FirstOrDefault(item => item.SubjectId == id);
            if (artifact is { LocalPath: { } diagnosticPath } &&
                IsWorkspaceFile(diagnosticPath) && VerifyHash(diagnosticPath, artifact.Sha256))
                return (GenericSubject(id, diagnosticOrigin.Id, diagnosticPath,
                    "providerDiagnostic", ImmutableArray<ScientificEvidence>.Empty,
                    ImmutableArray<ScientificFinding>.Empty, null), diagnosticOrigin);
        }
        var constructed = _constructedByAttempt.Values.FirstOrDefault(item => item.Id == id);
        if (constructed is not null &&
            _revisions.TryGetValue(constructed.Attempt.StudyRevisionId, out var constructionOrigin))
        {
            if (!IsWorkspaceFile(constructed.Molecule.CoordinatePath) ||
                !VerifyHash(constructed.Molecule.CoordinatePath, constructed.Molecule.CoordinateSha256))
                return (null, null);
            return (WithConstructionMeasurements(GenericSubject(id, constructionOrigin.Id,
                constructed.Molecule.CoordinatePath, "constructedSystem", constructed.Evidence,
                constructed.Findings, null), constructed), constructionOrigin);
        }
        if (_stages.TryGetValue(id, out var stage) &&
            _revisions.TryGetValue(stage.Attempt.StudyRevisionId, out var origin))
            return (WithGeometry(WithStageMeasurements(GenericSubject(id, origin.Id,
                    IsWorkspaceFile(stage.Molecule.CoordinatePath) &&
                    VerifyHash(stage.Molecule.CoordinatePath, stage.Molecule.CoordinateSha256)
                        ? stage.Molecule.CoordinatePath : null, "completedStage",
                stage.Observation.Evidence, stage.Findings,
                _stageAssessments.GetValueOrDefault(id)), stage,
                _policyByAttempt.GetValueOrDefault(stage.Attempt.Id)),
                stage.Observation.ProteinGeometry), origin);
        return (null, null);
    }

    private InspectionSubject GenericSubject(string id, string revisionId, string? path, string kind,
        ImmutableArray<ScientificEvidence> evidence, ImmutableArray<ScientificFinding> findings,
        PreparationAssessmentResult? assessment)
    {
        var relevantEvidence = evidence.IsDefault ? ImmutableArray<ScientificEvidence>.Empty :
            evidence.Where(item => item.SubjectId == id).ToImmutableArray();
        var relevantFindings = findings.IsDefault ? ImmutableArray<ScientificFinding>.Empty :
            findings.Where(item => item.SubjectId == id).ToImmutableArray();
        var annotations = relevantEvidence.Where(item => path is not null)
            .Select(item => new InspectionAnnotation(NewId(), id, item.Method, item.Observation, item.Id, null))
            .ToImmutableArray();
        var metrics = relevantEvidence.Select(item => new InspectionMetric(item.Method, item.Observation,
            null, id, item.Id)).ToImmutableArray();
        return new InspectionSubject(id, revisionId,
            path is not null && IsWorkspaceFile(path) ? StructureUrlLocked(path) : null, kind,
            ImmutableArray<string>.Empty, relevantEvidence, relevantFindings, annotations, metrics, assessment);
    }

    private static InspectionSubject WithConstructionMeasurements(InspectionSubject subject,
        ConstructedExplicitSystem constructed)
    {
        if (constructed.LocalState?.Standing != ObservationStanding.Observed ||
            constructed.LocalState.Measurements.IsDefault)
            return subject;
        var metrics = subject.Metrics.ToBuilder();
        var evidenceId = constructed.Evidence.FirstOrDefault()?.Id;
        foreach (var measure in constructed.LocalState.Measurements.Where(item => double.IsFinite(item.Value)))
            metrics.Add(new InspectionMetric(measure.Name,
                measure.Value.ToString("G6", CultureInfo.InvariantCulture), measure.Unit,
                measure.Scope, evidenceId));
        return subject with { Metrics = metrics.ToImmutable() };
    }

    // A measured distance is shown with its exact source addresses. Prepared
    // atom serials are not inferred from those addresses, so these metrics do
    // not advertise spatial focus until a verified view mapping exists.
    private static InspectionSubject WithGeometry(InspectionSubject subject, ProteinGeometryObservations? geometry)
    {
        if (geometry is null) return subject;
        var evidence = subject.Evidence.ToBuilder();
        var metrics = subject.Metrics.ToBuilder();
        var geometryEvidenceIds = ImmutableArray.CreateBuilder<string>();
        var kinds = geometry.Kinds.IsDefault
            ? ImmutableArray<ProteinGeometryKindObservation>.Empty : geometry.Kinds;
        for (var index = 0; index < kinds.Length; index++)
        {
            var kind = kinds[index];
            var evidenceId = DerivedInspectionEvidenceId(subject, "geometry-kind", kind.Kind,
                index, evidence);
            geometryEvidenceIds.Add(evidenceId);
            evidence.Add(new ScientificEvidence(evidenceId, subject.Id, "local scientific worker",
                kind.Kind, $"{kind.MeasuredCount} of {kind.EligibleCount} eligible observations; standing {kind.Standing}",
                $"Exact inspected {subject.RepresentationKind} {subject.Id}",
                kind.UnavailableReason ?? "Only the stated measurement scope was examined.",
                EvidenceBearing.Context));
            if (kind.MinimumDistanceAngstrom is double minimum)
                metrics.Add(new InspectionMetric($"{kind.Kind} minimum", minimum.ToString("G6", CultureInfo.InvariantCulture),
                    "Å", subject.Id, evidenceId));
            if (kind.MaximumDistanceAngstrom is double maximum)
                metrics.Add(new InspectionMetric($"{kind.Kind} maximum", maximum.ToString("G6", CultureInfo.InvariantCulture),
                    "Å", subject.Id, evidenceId));
        }
        var located = geometry.LocatedDistances.IsDefault
            ? ImmutableArray<GeometryDistanceObservation>.Empty : geometry.LocatedDistances;
        for (var index = 0; index < located.Length; index++)
        {
            var distance = located[index];
            var first = distance.First.Residue;
            var second = distance.Second.Residue;
            var label = $"{first.Chain}/{first.CopyId}:{first.Residue}{first.InsertionCode}/{distance.First.AtomName} ↔ " +
                $"{second.Chain}/{second.CopyId}:{second.Residue}{second.InsertionCode}/{distance.Second.AtomName}";
            var evidenceId = DerivedInspectionEvidenceId(subject, "located-geometry", label,
                index, evidence);
            geometryEvidenceIds.Add(evidenceId);
            evidence.Add(new ScientificEvidence(evidenceId, subject.Id, "local scientific worker", distance.Kind,
                $"Located atom-pair separation for {label}", $"Exact inspected {subject.RepresentationKind} {subject.Id}",
                "Addressed numerical observation; no unverified spatial focus is inferred.", EvidenceBearing.Context));
            metrics.Add(new InspectionMetric(label, distance.DistanceAngstrom.ToString("G6", CultureInfo.InvariantCulture),
                "Å", subject.Id, evidenceId));
        }
        return subject with { Evidence = evidence.ToImmutable(), Metrics = metrics.ToImmutable(),
            Geometry = new InspectionGeometryAccount(geometry, geometryEvidenceIds.ToImmutable()) };
    }

    private static InspectionSubject WithStageMeasurements(InspectionSubject subject, CompletedStage stage,
        ApplicablePreparationPolicy? policy)
    {
        var evidence = subject.Evidence.ToBuilder();
        var metrics = subject.Metrics.ToBuilder();
        var measured = stage.Observation.Measurements.IsDefault
            ? ImmutableArray<MeasuredValue>.Empty : stage.Observation.Measurements;
        for (var index = 0; index < measured.Length; index++)
        {
            var value = measured[index];
            if (!double.IsFinite(value.Value) || string.IsNullOrWhiteSpace(value.Name) ||
                string.IsNullOrWhiteSpace(value.Unit) || string.IsNullOrWhiteSpace(value.Scope)) continue;
            var evidenceId = DerivedInspectionEvidenceId(subject, "stage-measurement",
                value.Name + " · " + value.Scope, index, evidence);
            evidence.Add(new ScientificEvidence(evidenceId, subject.Id,
                stage.Observation.ProviderVersion, "Observed completed-stage measurement",
                $"{value.Name} in {value.Scope}: {value.Value.ToString("G6", CultureInfo.InvariantCulture)} {value.Unit}",
                $"Stage {stage.Id}; attempt {stage.Attempt.Id}; policy {stage.PolicyId}",
                "This numerical observation does not itself establish scientific qualification or atom-level focus.",
                EvidenceBearing.Context));
            metrics.Add(new InspectionMetric($"{value.Name} · {value.Scope}",
                value.Value.ToString("G6", CultureInfo.InvariantCulture), value.Unit, subject.Id, evidenceId));
        }
        var local = stage.Observation.LocalState;
        if (local is { Standing: ObservationStanding.Observed } && !local.RolePairMeasurements.IsDefault &&
            policy is { LocalStateObservation: { } spec } &&
            double.IsFinite(spec.ContactSearchRadiusAngstrom) && spec.ContactSearchRadiusAngstrom > 0)
        {
            var index = 0;
            foreach (var pair in local.RolePairMeasurements.OrderByDescending(item =>
                         item.FirstMoleculeRole == MoleculeRoleKind.Protein && item.SecondMoleculeRole == MoleculeRoleKind.Lipid))
            {
                var ordinal = index++;
                if (pair.PairsWithinSearchRadius < 0 || !Enum.IsDefined(pair.FirstMoleculeRole) ||
                    !Enum.IsDefined(pair.SecondMoleculeRole)) continue;
                var label = MoleculeRoleTokens.PairKey(pair.FirstMoleculeRole, pair.SecondMoleculeRole);
                var evidenceId = DerivedInspectionEvidenceId(subject, "role-pair-contact",
                    label, ordinal, evidence);
                evidence.Add(new ScientificEvidence(evidenceId, subject.Id,
                    stage.Observation.ProviderVersion, "Observed local role-pair contact",
                    $"{label}: {pair.PairsWithinSearchRadius} atom pairs within " +
                    $"{spec.ContactSearchRadiusAngstrom.ToString("G6", CultureInfo.InvariantCulture)} Å",
                    $"Stage {stage.Id}; attempt {stage.Attempt.Id}; policy {policy.Id}",
                    "Role-pair measurement only; it does not assign atom-level spatial focus or qualify the stage.",
                    EvidenceBearing.Context));
                metrics.Add(new InspectionMetric($"{label} pairs within search radius",
                    pair.PairsWithinSearchRadius.ToString(CultureInfo.InvariantCulture), "pairs", subject.Id,
                    evidenceId));
                if (pair.MinimumDistanceAngstrom is double minimum && double.IsFinite(minimum) && minimum >= 0)
                    metrics.Add(new InspectionMetric($"{label} nearest distance",
                        minimum.ToString("G6", CultureInfo.InvariantCulture), "Å", subject.Id, evidenceId));
            }
        }
        return subject with { Evidence = evidence.ToImmutable(), Metrics = metrics.ToImmutable() };
    }

    private static string DerivedInspectionEvidenceId(InspectionSubject subject, string family,
        string measurement, int ordinal, IEnumerable<ScientificEvidence> existing)
    {
        var identity = JsonSerializer.Serialize(new[] { subject.Id, subject.StudyRevisionId,
            family, measurement, ordinal.ToString(CultureInfo.InvariantCulture) });
        var digest = Convert.ToHexString(SHA256.HashData(
            System.Text.Encoding.UTF8.GetBytes(identity))).ToLowerInvariant();
        var basis = "inspection-derived-" + digest;
        var candidate = basis;
        for (var collision = 1; existing.Any(item => item.Id == candidate); collision++)
            candidate = basis + "-" + collision.ToString(CultureInfo.InvariantCulture);
        return candidate;
    }

    private WorkspaceState SnapshotLocked()
    {
        var retained = _workspace.Snapshot();
        var retainedAssessments = retained.Assessments.ToDictionary(item => item.StageId,
            StringComparer.Ordinal);
        var stages = retained.CompletedStages.Select(stage => new StageAccount(
            stage.Id, stage.Attempt.Id, stage.Attempt.StudyRevisionId, stage.Kind, "completed",
            retainedAssessments.GetValueOrDefault(stage.Id),
            stage.Attempt.StudyRevisionId == retained.CurrentStudy?.Id
                ? $"{stage.Kind} completed for the current study revision {stage.Attempt.StudyRevisionId}."
                : $"Historical {stage.Kind} result from study revision {stage.Attempt.StudyRevisionId}; " +
                  $"current study revision {retained.CurrentStudy?.Id ?? _study.Id}.",
            stage.Observation, _constructedByAttempt.TryGetValue(stage.Attempt.Id, out var stageSource)
                ? ConstructedAccount(stageSource) : null,
            _exportAccounts.TryGetValue(stage.Id, out var export) &&
            export.AssessmentId == retainedAssessments.GetValueOrDefault(stage.Id)?.Id
                ? export : null, stage.SourceStageId,
            StageProteinLabelLocked(stage.Attempt.StudyRevisionId),
            stage.Attempt.Route == ConstructionRouteKind.PackmolMemgen
                ? $"PACKMOL-Memgen {stage.Attempt.ConstructionProviderVersion}"
                : $"OpenMM native construction {stage.Attempt.ConstructionProviderVersion}",
            stage.Attempt.StartedAt)).ToImmutableArray();
        var declinedRequiredChange = _preparationProposals?.StudyRevisionId == _study.Id &&
            _preparationProposals.IntendedProteinId == _study.IntendedProtein?.Id
            ? _preparationProposals.Changes.FirstOrDefault(change =>
                change.Kind == PreparationChangeKind.HeavyAtom && _decisions.Any(decision =>
                    decision.StudyRevisionId == _study.Id && decision.SubjectId == change.Id &&
                    decision.ChosenValue == ResearcherDecisionValue.Declined))
            : null;
        var protein = _protein is not null
            ? new ProteinAccount(_protein.Id, "assessed", "Exact prepared protein; membrane placement remains separate.",
                _protein.Molecule.AtomCount, _preparationProposals?.Changes ?? ImmutableArray<PreparationChangeProposal>.Empty,
                _protein.Findings, PredictionAccount(_protein.Prediction), _protein.Geometry,
                _preparationProposals?.SourceGeometry)
            : _study.IntendedProtein is not null
                ? new ProteinAccount(_study.IntendedProtein.Id,
                    _proteinDiagnostic is not null ? "candidateNotQualified" :
                        declinedRequiredChange is not null ? "declined" :
                        _preparationProposals is null ? "notEstablished" : "review",
                    _proteinDiagnostic is not null
                        ? "The observed candidate geometry does not establish a prepared protein; inspect its located findings."
                        : declinedRequiredChange is not null
                            ? $"Required change {declinedRequiredChange.Id} at " +
                              $"{declinedRequiredChange.Residue.Chain}[{declinedRequiredChange.Residue.CopyId}]:" +
                              $"{declinedRequiredChange.Residue.Residue}{declinedRequiredChange.Residue.InsertionCode} " +
                              $"was declined for study revision {_study.Id}; no assessed prepared protein was established."
                        : _preparationProposals?.UnresolvedQuestions.FirstOrDefault() ??
                          "Review exact protein changes and chemical states.",
                    null, _preparationProposals?.Changes ?? ImmutableArray<PreparationChangeProposal>.Empty,
                    _proteinDiagnostic?.Findings ?? ImmutableArray<ScientificFinding>.Empty,
                    PredictionAccount(_sourceInspection?.Prediction), _proteinDiagnostic?.Geometry,
                    _preparationProposals?.SourceGeometry, _proteinDiagnostic?.Candidate.Id)
                : null;
        var selectedMembrane = _study.Membrane;
        var currentMembrane = selectedMembrane is not null && _membrane?.Intended.Id == selectedMembrane.Id
            ? _membrane : null;
        var membrane = selectedMembrane is null ? null : new MembraneAccount(
            selectedMembrane.Id, currentMembrane is not null ? "assessed" :
                _membraneAssessmentRunning ? "assessing" :
                _membraneAssessmentUnavailable ? "unavailable" : "notEstablished",
            selectedMembrane.ScientificPurpose,
            selectedMembrane.Upper.Fractions, selectedMembrane.Lower.Fractions,
            currentMembrane?.Limitations ?? ImmutableArray<string>.Empty,
            _membraneAssessmentReason,
            currentMembrane?.PolicyId, currentMembrane?.PolicyVersion,
            currentMembrane?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty,
            currentMembrane?.SpeciesRepresentations.Select(item => new MembraneSpeciesSupportAccount(
                item.SpeciesId, item.ChemistryId, item.Category, item.ForceFieldFamily,
                item.ForceFieldVersion, item.CoordinateTemplateSha256, item.TemplateSha256,
                item.Limitations)).ToImmutableArray() ?? ImmutableArray<MembraneSpeciesSupportAccount>.Empty);
        var placementPolicy = _placementProposal is not null && _study.Membrane is not null
            ? SelectPlacementPolicyLocked(_placementProposal, _study.Membrane) : null;
        var placementEvidence = _placementProposal?.Evidence.AddRange(
            _placementMeasurement?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty)
            .AddRange(_opmReview?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty) ??
            ImmutableArray<ScientificEvidence>.Empty;
        var placementLimitations = _placementProposal?.Limitations.AddRange(
            _placementMeasurement?.Limitations ?? ImmutableArray<string>.Empty)
            .AddRange(_opmReview?.Limitations ?? ImmutableArray<string>.Empty) ??
            ImmutableArray<string>.Empty;
        if (placementPolicy is not null && !placementPolicy.Limitations.IsDefaultOrEmpty)
            placementLimitations = placementLimitations.AddRange(placementPolicy.Limitations);
        if (!string.IsNullOrWhiteSpace(_placementMeasurementIssue))
            placementLimitations = placementLimitations.Add(_placementMeasurementIssue);
        var observedContacts = _placementMeasurement?.Residues
            .Where(item => item.AtomsWithinCore > 0)
            .Select(item => $"{item.Address.Chain}[{item.Address.CopyId}]:{item.Address.Residue}{item.Address.InsertionCode} — observed core atoms")
            .ToImmutableArray() ?? ImmutableArray<string>.Empty;
        var placement = _placementProposal is null ? null : new PlacementAccount(
            _placementProposal.Id, _placementAssessing ? "assessing" : _placement?.Standing switch
            {
                AssessmentStanding.Supported => "supported",
                AssessmentStanding.Unsupported => "unsupported",
                AssessmentStanding.NotEstablished => "notEstablished",
                _ => "proposed"
            },
            _placementProposal.TopologyKind,
            _placementProposal.MidplaneAngstrom, _placementProposal.TiltDegrees,
            _placementProposal.BiologicalSidedness,
            _placementAssessing ? "Checking the exact construct, rigid position and membrane frame." :
            _placement is null ? "The positioned candidate still needs technical checks." :
                _placement.Reason +
                    (_placementMeasurementIssue is null ? string.Empty : " " + _placementMeasurementIssue),
            placementEvidence, _placementMeasurement?.Prediction,
            _placementProposal.PreparedProteinId, _placementProposal.MembraneModelId,
            _placementProposal.PhysicalSide, _placementProposal.MidplaneAngstrom,
            _placementProposal.ThicknessAngstrom,
            _placementProposal.ContactingRegions.AddRange(observedContacts),
            placementLimitations, placementPolicy?.Id, placementPolicy?.Version, _placementWitness?.Id,
            _placementProposal.Transform);
        var studyAccount = new StudyAccount(_study.Id, _study.Number,
            _study.IntendedProtein is null ? "Choose an exact protein structure." : "One identified protein–membrane study.",
            _selectedSource?.Id, _study.IntendedProtein?.ModelIndex,
            _study.IntendedProtein?.BiologicalAssemblyId,
            _study.IntendedProtein?.Chains.Select(item => item.CopyId).ToImmutableArray() ?? ImmutableArray<string>.Empty,
            _study.IntendedProtein?.Partners ?? ImmutableArray<PartnerSelection>.Empty,
            _study.IntendedProtein?.AlternateLocations ?? ImmutableArray<AlternateLocationChoice>.Empty,
            _study.Conditions, _selectedSource?.Kind, _selectedSource?.UploadProvenance,
            _selectedSource?.UploadProvenanceNote, _study.AdoptedPlacementProposalId,
            _selectedSource?.Kind == SourceRouteKind.Upload
                ? _selectedSource.SourceModelDescription ?? "Uploaded structure"
                : _selectedSource?.Accession ?? _selectedSource?.SourceModelDescription);
        var shownAttempt = _currentAttempt is null ? null : retained.CurrentAttempt;
        var shownExecution = _currentAttempt is null ? _execution : retained.CurrentExecution;
        var preparationReview = BuildPreparationReviewLocked();
        return new WorkspaceState(_revision, studyAccount,
            _candidates.Select(candidate => new SourceCandidateAccount(candidate.Id, candidate.Label,
                candidate.Kind, candidate.Provenance, candidate.Limitations)).ToImmutableArray(),
            _sourceInspection?.Models ?? ImmutableArray<SourceModelObservation>.Empty,
            LipidsLocked().Values.OrderBy(item => item.SpeciesId).Select(item => new LipidCatalogueAccount(
                item.SpeciesId, item.SpeciesId, item.ChemistryId, item.Limitations)).ToImmutableArray(),
            protein, membrane, placement,
            AttemptAccountLocked(shownAttempt, shownExecution),
            stages, _inspection.Current, ActionsLocked(stages),
            ActiveIssuesLocked(AttemptAccountLocked(shownAttempt, shownExecution), stages, preparationReview),
            PredictionAccount(_sourceInspection?.Prediction), preparationReview,
            BuildProteinTaskLocked(preparationReview), BuildPlacementTaskLocked(),
            ImmutableArray.Create(
                new PlacementMethodAccount("OPM", _selectedSource?.Kind == SourceRouteKind.Rcsb ? "lookupEligible" : "notApplicable",
                    _selectedSource?.Kind == SourceRouteKind.Rcsb ? null : "OPM lookup requires an identified RCSB entry; user positioning remains available."),
                new PlacementMethodAccount("PPM", CanRunPpmLocked() ? "configured" : "unavailable",
                    CanRunPpmLocked() ? null : "The local PPM executable or residue library is not identified and hash-verified.")),
            BuildPreparationPlanAccountLocked(), BuildConstructionRoutesLocked(),
            _priorAttempts.ToImmutableArray());
    }

    private ImmutableArray<WorkspaceNotice> ActiveIssuesLocked(AttemptAccount? attempt,
        ImmutableArray<StageAccount> stages, ProteinPreparationReviewAccount? review)
    {
        var issues = ImmutableArray.CreateBuilder<WorkspaceNotice>();
        var keys = new HashSet<string>(StringComparer.Ordinal);
        void Add(string condition, string severity, string? message, string? subject,
            string correctionArea, params string[] affectedAreas)
        {
            if (string.IsNullOrWhiteSpace(message)) return;
            var key = $"{condition}:{subject ?? "study"}";
            if (!keys.Add(key)) return;
            issues.Add(new WorkspaceNotice(key, severity, message, subject, key,
                affectedAreas.ToImmutableArray(), correctionArea));
        }

        Add("catalogue", "error", _catalogueIssue, null, "preparation", "preparation");
        if (_selectedSource is null)
            for (var index = 0; index < _sourceSearchIssues.Length; index++)
                Add($"source-search-{index}", "warning", _sourceSearchIssues[index], null,
                    "protein", "protein");
        if (_study.IntendedProtein is { } intended && _protein is null)
        {
            Add("protein-selection", "error", _proteinSelectionIssue?.Message, intended.Id, "protein", "protein");
            Add("preparation-outcome", "error", _proteinPreparationObservationIssue ??
                _proteinPreparationFailure, intended.Id, "protein", "protein", "placement", "preparation");
            if (_preparationPlan is null && !_recommendationRunning)
                Add("recommendation", "warning", _recommendationIssue, intended.Id, "protein", "protein");
            if (review is not null)
                for (var index = 0; index < review.Blockers.Length; index++)
                    Add($"preparation-blocker-{index}", "error", review.Blockers[index], intended.Id,
                        "protein", "protein", "placement", "preparation");
        }
        if (_study.Membrane is { } chosen && _membrane is null)
            Add("membrane-check", "error", _membraneAssessmentReason, chosen.Id,
                "membrane", "membrane", "placement", "preparation");
        if (_placementProposal is { } proposed && _study.Membrane is not null)
        {
            Add("placement-measurement", "error", _placementMeasurementIssue, proposed.Id,
                "placement", "placement", "preparation");
            if (_placement is { Standing: not AssessmentStanding.Supported })
                Add("placement-check", "error", _placement.Reason, proposed.Id,
                    "placement", "placement", "preparation");
        }
        if (_placementRouteOutcome?.Standing is "failed" or "unobserved" &&
            _placementRouteOutcome.StudyRevisionId == _study.Id)
            Add("orientation-route", "warning", _placementRouteOutcome.Message,
                _placementRouteOutcome.RequestId, "placement", "placement");
        else if (_placementOperationIssue is not null && _placementRouteOutcome is null)
            Add("placement-operation", "warning", _placementOperationIssue,
                _placementProposal?.Id, "placement", "placement");
        if (attempt is { Status: "failed" or "resourceRefused" or "unobserved" } &&
            attempt.StudyRevisionId == _study.Id)
            Add("attempt-outcome", "error", attempt.Message, attempt.AttemptId,
                "preparation", "preparation");
        foreach (var stage in stages)
        {
            if (stage.Assessment is { CurrentlyApplicable: true,
                    CheckStanding: not PreparationCheckStanding.ChecksPassed } assessment)
                Add("stage-check", "error", assessment.Reason, stage.StageId, "results", "results");
            if (stage.Export is { Status: "failed" } export)
                Add("stage-export", "warning", export.Reason, stage.StageId, "results", "results");
        }
        return issues.ToImmutable();
    }

    private string StageProteinLabelLocked(string studyRevisionId)
    {
        if (!_revisions.TryGetValue(studyRevisionId, out var revision) ||
            revision.IntendedProtein is not { } intended)
            return "Protein source unavailable in the retained account";
        var source = intended.Source;
        var name = !string.IsNullOrWhiteSpace(source.Accession) ? source.Accession :
            !string.IsNullOrWhiteSpace(source.SourceModelDescription) ? source.SourceModelDescription :
            source.Kind == SourceRouteKind.Upload ? "Uploaded structure" : "Selected protein source";
        var assembly = intended.BiologicalAssemblyId is { Length: > 0 } assemblyId
            ? $"assembly {assemblyId}" : "deposited coordinates";
        var chains = intended.Chains.IsDefaultOrEmpty ? "chains unavailable" :
            "chains " + string.Join(", ", intended.Chains.Select(chain =>
                chain.SourceChain == chain.CopyId ? chain.SourceChain :
                    $"{chain.SourceChain} copy {chain.CopyId}"));
        var sourceModelId = string.IsNullOrWhiteSpace(intended.SourceModelId)
            ? (intended.ModelIndex + 1).ToString(System.Globalization.CultureInfo.InvariantCulture)
            : intended.SourceModelId.Trim();
        return $"{name} · coordinate model {sourceModelId} · {assembly} · {chains}";
    }

    private AttemptAccount? AttemptAccountLocked(PreparationAttempt? attempt,
        StageExecutionState? execution)
    {
        if (attempt is null && execution is null) return null;
        return new AttemptAccount(
            attempt?.Id ?? execution?.AttemptId ?? string.Empty,
            execution?.Standing switch
            {
                StageExecutionStanding.Pending => "pending",
                StageExecutionStanding.Running => "running",
                StageExecutionStanding.Completed => "completed",
                StageExecutionStanding.Stopped => "stopped",
                StageExecutionStanding.Failed => "failed",
                StageExecutionStanding.ResourceRefused => "resourceRefused",
                StageExecutionStanding.Unobserved => "unobserved",
                _ => "pending"
            },
            execution?.Kind, execution?.Progress,
            execution?.Message ?? "The local attempt has not started.",
            attempt?.StudyRevisionId, attempt?.PolicyId, attempt?.PolicyVersion,
            execution?.StageId,
            attempt is not null ? _derivationByAttempt.GetValueOrDefault(attempt.Id) : null,
            attempt is not null && _constructedByAttempt.TryGetValue(attempt.Id, out var constructed)
                ? ConstructedAccount(constructed) : null,
            attempt is not null && _stopRequestedAttemptId == attempt.Id,
            attempt is not null && _trialsByAttempt.TryGetValue(attempt.Id, out var trials)
                ? trials : ImmutableArray<ConstructionTrialSummary>.Empty,
            execution?.Phase, execution?.TrialId, execution?.TrialIndex,
            execution?.FailureCode);
    }

    private static PredictionEvidenceAccount? PredictionAccount(PredictionEvidenceObservations? observations) =>
        observations is null ? null : new PredictionEvidenceAccount(observations.RecordId,
            observations.LocalConfidence, observations.PaeStanding, observations.PaeReason,
            observations.PaeAxisResidueCount, observations.Limitations);

    private static ConstructedSystemAccount ConstructedAccount(ConstructedExplicitSystem source) =>
        new(source.Id, source.Attempt.Id, source.Molecule.AtomCount, source.AchievedComposition,
            source.ActualCellAngstrom, source.Derivation.WaterCount, source.Derivation.SodiumCount,
            source.Derivation.ChlorideCount, source.ConditionsTreatment, source.LocalState,
            source.MaximumProteinCoordinateDeviationAngstrom);

    private ImmutableArray<AvailableAction> ActionsLocked(ImmutableArray<StageAccount> stages)
    {
        var actions = ImmutableArray.CreateBuilder<AvailableAction>();
        var preparationReview = BuildPreparationReviewLocked();
        void Add(ActorActionKind kind, string? subject, bool enabled, string reason) =>
            actions.Add(new AvailableAction(kind, subject, enabled, enabled ? null : reason));
        Add(ActorActionKind.SearchSource, null, true, "");
        Add(ActorActionKind.SelectSource, null, true, "");
        Add(ActorActionKind.SelectProteinModel, null, _sourceInspection is not null, "Choose and inspect a structural source first.");
        foreach (var change in _preparationProposals?.Changes ?? ImmutableArray<PreparationChangeProposal>.Empty)
        {
            var siteBlocker = preparationReview?.Decisions.FirstOrDefault(item =>
                item.Options.Any(option => option.ProposalId == change.Id))?.Blocker;
            var alternativeChosen = change.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState &&
                _preparationProposals!.Changes.Any(other => other.Id != change.Id &&
                    other.Kind == change.Kind && other.Residue == change.Residue &&
                    _decisions.Any(decision => decision.SubjectId == other.Id &&
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved));
            var undecided = _protein is null && !alternativeChosen && !_decisions.Any(item => item.SubjectId == change.Id);
            var option = preparationReview?.Decisions.SelectMany(item => item.Options)
                .FirstOrDefault(item => item.ProposalId == change.Id);
            Add(ActorActionKind.ApprovePreparationChange, change.Id,
                siteBlocker is null && undecided && option?.ConfirmationBlocker is null,
                siteBlocker ?? (undecided ? option?.ConfirmationBlocker ?? "This option is not ready for confirmation." : "This proposal has already been decided."));
            Add(ActorActionKind.DeclinePreparationChange, change.Id,
                change.Kind is PreparationChangeKind.HeavyAtom or PreparationChangeKind.Disulfide &&
                siteBlocker is null && undecided,
                change.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.AlternateLocation
                    ? "Choose one exclusive alternative; individual decline is not available."
                    : siteBlocker ?? "This proposal has already been decided.");
        }
        var planAccount = BuildPreparationPlanAccountLocked();
        Add(ActorActionKind.AuthorizePreparationPlan, _preparationPlan?.Id,
            planAccount is { AuthorizationAvailable: true } && !_proteinPreparationRunning,
            planAccount?.Message ?? "Calculate a current, jointly checked preparation plan first.");
        foreach (var change in _preparationProposals?.Changes.Where(item =>
                     item.Kind == PreparationChangeKind.ResidueState) ?? [])
            Add(ActorActionKind.OverridePreparationPlanChoice, change.Id,
                _preparationPlan is not null && _authorizedPreparationPlan is null && !_recommendationRunning,
                "A current plan is needed before changing this suggested state.");
        Add(ActorActionKind.RetryPreparationPlan, null,
            planAccount?.Standing == "failed" && !_recommendationRunning,
            "Retry is available after a failed recommendation calculation.");
        Add(ActorActionKind.StartProteinPreparation, null,
            _protein is null && !_proteinPreparationRunning &&
            _authorizedPreparationPlan is null && _preparationPlan is null &&
            preparationReview is { Blockers.Length: 0, RemainingCount: 0 } &&
            _preparationProposals?.StudyRevisionId == _study.Id,
            _preparationPlan is not null ? "Authorize the checked plan to prepare with its proposed choices." :
            "Resolve the current structural exceptions and exact manual choices before preparation.");
        Add(ActorActionKind.RetryProteinPreparation, null,
            _proteinPreparationFailure is not null && _protein is null && !_proteinPreparationRunning &&
            _proteinPreparationRetryable && _preparationProposals?.StudyRevisionId == _study.Id,
            "Retry is available only after a recoverable failure of the current reviewed protein.");
        Add(ActorActionKind.AdoptMembrane, null, !_membraneAssessmentRunning,
            "This membrane support assessment is already running.");
        Add(ActorActionKind.RetryMembraneCheck, null,
            _study.Membrane is not null &&
            _membraneAssessmentUnavailable && !_membraneAssessmentRunning,
            "Retry is available only after an unfinished check of the current selected membrane.");
        Add(ActorActionKind.ProposePlacement, null, !_placementRunning && _protein is not null && _study.Membrane is not null,
            _placementRunning ? "An exact placement operation is running." :
            _protein is null ? "Prepare a protein first." :
            _study.Membrane is null ? "Choose a membrane first." :
            "The prepared protein and chosen membrane must belong to the current study.");
        Add(ActorActionKind.RevisePlacement, null, !_placementRunning && _placementProposal is not null,
            _placementRunning ? "An exact placement operation is running." : "A positioned candidate is needed first.");
        Add(ActorActionKind.AdoptPlacement, null, !_placementRunning && _placement?.Standing == AssessmentStanding.Supported &&
            _placementProposal is { } currentPlacement &&
            _placementMeasurement is { Evidence.IsDefaultOrEmpty: false } &&
            IsWorkspaceFile(currentPlacement.OrientedProtein.CoordinatePath) &&
            VerifyHash(currentPlacement.OrientedProtein.CoordinatePath,
                currentPlacement.OrientedProtein.CoordinateSha256) &&
            _study.AdoptedPlacementProposalId != _placementProposal?.Id,
            "A current supported proposal, its exact assessment evidence and intact oriented coordinates are required.");
        var constructionRoutes = BuildConstructionRoutesLocked();
        var constructionReady = _placement is { Standing: AssessmentStanding.Supported } &&
            _study.AdoptedPlacementProposalId == _placement.Proposal.Id &&
            _attemptTask is not { IsCompleted: false } &&
            _protein is not null && _membrane is not null &&
            _protein.StudyRevisionId == _study.Id && _membrane.StudyRevisionId == _study.Id &&
            _placement.StudyRevisionId == _study.Id;
        var constructionBlocker = _protein is null ? "Prepare a protein first." :
            _membrane is null ? "Choose and check a membrane composition first." :
            _placement is null or { Standing: not AssessmentStanding.Supported } ?
                "Check a position for this exact protein and membrane first." :
            _study.AdoptedPlacementProposalId != _placement.Proposal.Id ?
                "Use the checked position before constructing the system." :
            _attemptTask is { IsCompleted: false } ? "The current system operation is still running." :
            constructionRoutes.IsDefaultOrEmpty ?
                $"Construction is unavailable for the chosen {string.Join("/", _membrane.SpeciesRepresentations.Select(item => item.SpeciesId))} composition; no verified recipe is enabled." :
                "The exact protein, membrane and position must belong to the current study revision.";
        if (constructionRoutes.IsDefaultOrEmpty)
            Add(ActorActionKind.BuildAndMinimize, null, false, constructionBlocker);
        else
            foreach (var route in constructionRoutes)
                Add(ActorActionKind.BuildAndMinimize, route.PolicyId,
                    constructionReady && route.Available,
                    !route.Available ? route.Reason ?? "The selected construction provider is unavailable." :
                        constructionBlocker);
        Add(ActorActionKind.StopAttempt, null, _currentAttempt is not null &&
            _execution?.AttemptId == _currentAttempt.Id &&
            (_attemptTask is { IsCompleted: false } && _attemptStop is not null &&
                _stopRequestedAttemptId != _currentAttempt.Id &&
                _execution.Standing is StageExecutionStanding.Pending or StageExecutionStanding.Running),
            "No identified unfinished attempt can be stopped.");
        Add(ActorActionKind.SelectInspectionSubject, null, _sourceInspection is not null || _study.Membrane is not null ||
            _placementProposal is not null || _proteinDiagnostic is not null ||
            _constructedByAttempt.Count > 0 || stages.Length > 0,
            "There is no identified subject to inspect.");
        Add(ActorActionKind.SetInspectionFocus, null, _inspection.Current is not null, "Select an inspection subject first.");
        foreach (var stage in stages)
        {
            var completed = _stages[stage.StageId];
            var exportReason = ExportGateReasonLocked(completed);
            Add(ActorActionKind.ExportStage, stage.StageId, exportReason is null,
                exportReason ?? "Export this completed stage with its currently applicable assessment.");
            if (completed.Kind == StageKind.Minimization)
                Add(ActorActionKind.RequestEquilibration, stage.StageId,
                    _attemptTask is not { IsCompleted: false } &&
                    CanRequestEquilibrationLocked(completed),
                    "No validated applicable optional procedure is available for this completed minimized stage.");
        }
        return actions.ToImmutable();
    }

    private void AdvanceStudyLocked(IntendedProteinModel? intended, MembraneModel? membrane, string? placementId)
    {
        _study = new StudyRevision(NewId(), _study.Number + 1, intended, membrane, placementId, _study.Conditions);
        _revisions.Add(_study.Id, _study);
        _workspace.RetainStudy(_study);
        _inspection.Clear();
        _proteinDiagnostic = null;
        _proteinSelectionRunning = false;
        _proteinSelectionIssue = null;
        _proteinPreparationRunning = false;
        _proteinPreparationFailure = null;
        _proteinPreparationRetryable = false;
        _proteinPreparationObservationIssue = null;
        _membraneAssessmentRunning = false;
        _membraneAssessmentUnavailable = false;
        _placementRunning = false;
        _placementAssessing = false;
        _placementOperationIssue = null;
        _placementRouteOutcome = null;
    }

    private AssessedPreparedProtein? CarryProteinLocked(StudyRevision revision)
    {
        var chemicalPolicy = revision.IntendedProtein is null ? null : SelectChemicalPolicyLocked(
            _sourceInspection?.Models.FirstOrDefault(item => item.Index == revision.IntendedProtein.ModelIndex),
            revision.IntendedProtein);
        var structuralPolicy = revision.IntendedProtein is null ? null : SelectStructuralPolicyLocked(
            _sourceInspection?.Models.FirstOrDefault(item => item.Index == revision.IntendedProtein.ModelIndex),
            revision.IntendedProtein);
        if (_protein is null || revision.IntendedProtein?.Id != _protein.Intended.Id ||
            _protein.Intended.Source.Sha256 != revision.IntendedProtein.Source.Sha256 ||
            chemicalPolicy is null || structuralPolicy is null ||
            _protein.ChemicalStatePolicyId != chemicalPolicy.Id ||
            _protein.ChemicalStatePolicyVersion != chemicalPolicy.Version ||
            _protein.StructuralAssessmentPolicyId != structuralPolicy.Id ||
            _protein.StructuralAssessmentPolicyVersion != structuralPolicy.Version) return null;
        return _protein with { StudyRevisionId = revision.Id };
    }

    private AssessedMembraneModel? CarryMembraneLocked(StudyRevision revision)
    {
        var policy = _membrane is null ? null : SelectMembranePolicyLocked(_membrane.Intended);
        if (_membrane is null || revision.Membrane?.Id != _membrane.Intended.Id ||
            _membrane.Intended.Conditions != revision.Conditions ||
            policy is null || policy.Id != _membrane.PolicyId ||
            policy.Version != _membrane.PolicyVersion) return null;
        return _membrane with { StudyRevisionId = revision.Id };
    }

    private ProteinChemicalStatePolicy? SelectChemicalPolicyLocked(SourceModelObservation? model,
        IntendedProteinModel? intended) => SelectChemicalPolicyWithCompatibilityLocked(model, intended).Policy;

    private (ProteinChemicalStatePolicy? Policy, ImmutableArray<SourcePartnerObservation> UnsupportedPartners,
        bool CorePolicyAvailable, bool AllPartnersKnown) SelectChemicalPolicyWithCompatibilityLocked(
        SourceModelObservation? model, IntendedProteinModel? intended)
    {
        if (model is null || intended is null || _catalogue is null)
            return (null, ImmutableArray<SourcePartnerObservation>.Empty, false, false);
        var corePolicy = SelectChemicalPolicyForProteinCoreLocked(model, intended);
        var unsupported = ImmutableArray.CreateBuilder<SourcePartnerObservation>();
        var allPartnersKnown = true;
        foreach (var partner in intended.Partners.Where(partner => partner.Retain))
        {
            var source = model.Partners.FirstOrDefault(item => item.SourceId == partner.SourceId);
            if (source is null) allPartnersKnown = false;
            else if (!SupportedRetainedPartner(source)) unsupported.Add(source);
        }
        return (allPartnersKnown && unsupported.Count == 0 ? corePolicy : null,
            unsupported.ToImmutable(), corePolicy is not null, allPartnersKnown);
    }

    private static bool SupportedRetainedPartner(SourcePartnerObservation source) => source.Kind switch
    {
        "water" => source.Label is ("HOH" or "WAT" or "H2O") && source.AtomCount is (1 or 3),
        "ion" => source.Label is ("NA" or "CL") && source.AtomCount == 1,
        _ => false
    };

    private static string UnsupportedPartnerLabel(SourcePartnerObservation source)
    {
        var name = source.DisplayName?.Trim();
        var component = string.IsNullOrEmpty(name) || name.All(char.IsDigit) ||
            string.Equals(name, source.Label, StringComparison.OrdinalIgnoreCase)
                ? source.Label : $"{name} ({source.Label})";
        var chain = string.IsNullOrWhiteSpace(source.Chain) ? "" : $" · Chain {source.Chain}";
        var residue = source.Residue is null ? "" : $" · Residue {source.Residue}{source.InsertionCode}";
        return component + chain + residue;
    }

    private ProteinChemicalStatePolicy? SelectChemicalPolicyForProteinCoreLocked(
        SourceModelObservation model, IntendedProteinModel intended)
    {
        if (_catalogue is null) return null;
        var canonical = new HashSet<string>("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
            .Split(' '), StringComparer.Ordinal);
        var selectedChainResidues = model.Residues.Where(item =>
            intended.Chains.Any(chain => chain.SourceChain == item.Address.Chain)).ToArray();
        if (selectedChainResidues.Any(item => item.ResidueKind is not (SourceResidueKind.Protein or SourceResidueKind.Solvent or SourceResidueKind.Heterogen)) ||
            selectedChainResidues.Any(item => item.ResidueKind == SourceResidueKind.Heterogen &&
                !model.Partners.Any(partner => partner.Chain == item.Address.Chain &&
                    partner.Residue == item.Address.Residue))) return null;
        // Source solvent is not a retained protein residue. Other non-protein
        // entities must have an explicit partner disposition, and noncanonical
        // residues classified as protein still fail the supported core.
        var selectedResidues = selectedChainResidues.Where(item => item.ResidueKind == SourceResidueKind.Protein).ToArray();
        if (selectedResidues.Length == 0 || selectedResidues.Any(item => !canonical.Contains(item.Name)))
            return null;
        var applicable = _catalogue.ProteinChemicalStates.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty &&
            policy.ApplicableChemistry == "canonical-amino-acid-assembly" &&
            double.IsFinite(policy.DisulfideCandidateMaxSgDistanceAngstrom) &&
            policy.DisulfideCandidateMaxSgDistanceAngstrom > 0 &&
            policy.DefaultVariants is not null && !policy.PermittedVariants.IsDefaultOrEmpty &&
            !policy.ForceFieldFiles.IsDefaultOrEmpty && policy.ForceFieldFiles.All(VerifiedAsset)).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private ProteinStructuralAssessmentPolicy? SelectStructuralPolicyLocked(SourceModelObservation? model,
        IntendedProteinModel? intended)
    {
        if (model is null || intended is null || _catalogue is null) return null;
        var selected = model.Residues.Where(residue => residue.ResidueKind == SourceResidueKind.Protein &&
            intended.Chains.Any(chain => chain.SourceChain == residue.Address.Chain)).ToArray();
        if (selected.Length == 0 || selected.Any(residue => !residue.BackboneHeavyAtomsComplete)) return null;
        var applicable = _catalogue.ProteinStructuralPolicies.Where(policy =>
            policy is not null && !string.IsNullOrWhiteSpace(policy.Id) &&
            !string.IsNullOrWhiteSpace(policy.Version) && !policy.EvidenceReferences.IsDefaultOrEmpty &&
            policy.ApplicableProteinClass == "canonical-amino-acid-assembly" &&
            policy.Measurement is { } measurement && !measurement.RequiredKinds.IsDefaultOrEmpty &&
            measurement.RequiredKinds.Distinct(StringComparer.Ordinal).Count() == measurement.RequiredKinds.Length &&
            measurement.AtomRadiusByElementAngstrom is { Count: > 0 } &&
            measurement.AtomRadiusByElementAngstrom.All(item =>
                !string.IsNullOrWhiteSpace(item.Key) && double.IsFinite(item.Value) && item.Value > 0) &&
            double.IsFinite(measurement.NeighborSearchRadiusAngstrom) &&
            measurement.NeighborSearchRadiusAngstrom > 0 && measurement.ExcludedBondHops >= 0 &&
            measurement.MaximumReportedPairs > 0 && !policy.Criteria.IsDefaultOrEmpty &&
            measurement.RequiredKinds.All(kind => policy.Criteria.Count(criterion => criterion.Kind == kind &&
                (criterion.MinimumObservedAngstrom is null ||
                 double.IsFinite(criterion.MinimumObservedAngstrom.Value)) &&
                (criterion.MaximumObservedAngstrom is null ||
                 double.IsFinite(criterion.MaximumObservedAngstrom.Value))) == 1)).Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private MembraneSupportPolicy? SelectMembranePolicyLocked(MembraneModel membrane)
    {
        if (_catalogue is null) return null;
        var applicable = _catalogue.MembranePolicies.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty && !policy.CoveredSpeciesIds.IsDefaultOrEmpty &&
            membrane.Upper.Fractions.Concat(membrane.Lower.Fractions).Where(item => item.Fraction > 0).All(item =>
                policy.CoveredSpeciesIds.Contains(item.SpeciesId))).Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private PlacementSupportPolicy? SelectPlacementPolicyLocked(PlacementProposal proposal, MembraneModel membrane)
    {
        var applicable = ApplicablePlacementPoliciesLocked(proposal, membrane);
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private string PlacementPolicyIssueLocked(PlacementProposal proposal, MembraneModel membrane) =>
        ApplicablePlacementPoliciesLocked(proposal, membrane).Length > 1
            ? "Multiple applicable placement support policies are present; ambiguous placement support policies cannot be selected."
            : "No placement support policy covers this exact topology and membrane.";

    private PlacementSupportPolicy[] ApplicablePlacementPoliciesLocked(PlacementProposal proposal,
        MembraneModel membrane) => ApplicablePlacementPoliciesLocked(proposal.TopologyKind, membrane);

    private PlacementSupportPolicy[] ApplicablePlacementPoliciesLocked(ProteinTopologyKind topologyKind,
        MembraneModel membrane)
    {
        if (_catalogue is null) return [];
        return _catalogue.PlacementPolicies.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty &&
            !policy.CoveredTopologyKinds.IsDefaultOrEmpty &&
            policy.CoveredTopologyKinds.All(kind => Enum.IsDefined(kind)) &&
            !policy.CoveredSpeciesIds.IsDefaultOrEmpty &&
            double.IsFinite(policy.InterfaceBandAngstrom) && policy.InterfaceBandAngstrom > 0 &&
            policy.CoveredTopologyKinds.Contains(topologyKind) &&
            membrane.Upper.Fractions.Concat(membrane.Lower.Fractions)
                .Where(item => item.Fraction > 0)
                .All(item => policy.CoveredSpeciesIds.Contains(item.SpeciesId)))
            .Take(2).ToArray();
    }

    private PlacementStructuralWitness? SelectPlacementWitnessLocked(AssessedPreparedProtein protein,
        AssessedMembraneModel membrane, PlacementProposal proposal)
    {
        if (_catalogue is null || string.IsNullOrWhiteSpace(proposal.BiologicalSidedness)) return null;
        var intended = protein.Intended;
        var applicable = _catalogue.PlacementWitnesses.Where(witness =>
            witness is not null && !string.IsNullOrWhiteSpace(witness.Id) && !string.IsNullOrWhiteSpace(witness.Version) &&
            !string.IsNullOrWhiteSpace(witness.Source) && !witness.EvidenceReferences.IsDefaultOrEmpty &&
            !witness.Residues.IsDefaultOrEmpty && witness.Limitations.IsDefaultOrEmpty &&
            witness.SourceCoordinateSha256 == intended.Source.Sha256 &&
            witness.SourceModelIndex == intended.ModelIndex &&
            witness.BiologicalAssemblyId == intended.BiologicalAssemblyId &&
            !witness.ChainCopies.IsDefault && witness.ChainCopies.ToHashSet().SetEquals(intended.Chains) &&
            !witness.RetainedPartnerSourceIds.IsDefault &&
            witness.RetainedPartnerSourceIds.ToHashSet(StringComparer.Ordinal).SetEquals(
                intended.Partners.Where(item => item.Retain).Select(item => item.SourceId)) &&
            (witness.PreparedCoordinateSha256 is null ||
                witness.PreparedCoordinateSha256 == protein.Molecule.CoordinateSha256) &&
            !string.IsNullOrWhiteSpace(witness.PreparedBondGraphSha256) &&
            witness.PreparedBondGraphSha256 == protein.Molecule.TopologySha256 &&
            witness.Upper is not null && SameLeaflet(witness.Upper, membrane.Intended.Upper) &&
            witness.Lower is not null && SameLeaflet(witness.Lower, membrane.Intended.Lower) &&
            witness.Conditions == membrane.Intended.Conditions &&
            witness.TopologyKind == proposal.TopologyKind &&
            witness.BiologicalSidedness == proposal.BiologicalSidedness).Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private static bool SameLeaflet(LeafletComposition left, LeafletComposition right) =>
        !left.Fractions.IsDefaultOrEmpty && !right.Fractions.IsDefaultOrEmpty &&
        left.PhysicalSide == right.PhysicalSide && left.Fractions.Length == right.Fractions.Length &&
        left.Fractions.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() == left.Fractions.Length &&
        right.Fractions.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() == right.Fractions.Length &&
        left.Fractions.OrderBy(item => item.SpeciesId, StringComparer.Ordinal)
            .Zip(right.Fractions.OrderBy(item => item.SpeciesId, StringComparer.Ordinal))
            .All(pair => pair.First.SpeciesId == pair.Second.SpeciesId &&
                double.IsFinite(pair.First.Fraction) && double.IsFinite(pair.Second.Fraction) &&
                Math.Abs(pair.First.Fraction - pair.Second.Fraction) <= 1e-9);

    private bool CanRunPpmLocked() => _catalogue is not null &&
        !string.IsNullOrWhiteSpace(_catalogue.PpmVersion) &&
        VerifyHash(_ppmExecutablePath, _catalogue.PpmExecutableSha256) &&
        VerifyHash(_catalogue.PpmResidueLibraryPath, _catalogue.PpmResidueLibrarySha256);

    private static bool ConstructionProviderMatches(ConstructionPolicy policy,
        ConstructionProviderInstallation? installed)
    {
        if (installed is null)
            return false;
        if (policy.Route == ConstructionRouteKind.PackmolMemgen)
            return policy.ProviderName == "PACKMOL-Memgen" &&
                string.Equals(installed.MemgenVersion, policy.ProviderVersion, StringComparison.Ordinal) &&
                !string.IsNullOrWhiteSpace(installed.MemgenHome) &&
                !policy.ProviderAssets.IsDefaultOrEmpty &&
                policy.ProviderAssets.All(asset =>
                    Path.GetFullPath(asset.Path).StartsWith(
                        Path.GetFullPath(installed.MemgenHome) + Path.DirectorySeparatorChar,
                        StringComparison.Ordinal) && VerifyHash(asset.Path, asset.Sha256));
        if (policy.Route != ConstructionRouteKind.NativeOpenMm ||
            !string.Equals(installed.FullVersion, policy.ProviderVersion, StringComparison.Ordinal))
            return false;
        var mode = policy.NativePatchMode ?? "installed";
        if (mode == "installed")
        {
            if (policy.LipidTypeArgument == "DMPC")
                return string.Equals(installed.NativePatchPath, policy.NativePatchPath,
                           StringComparison.Ordinal) &&
                       string.Equals(installed.NativePatchSha256, policy.NativePatchSha256,
                           StringComparison.OrdinalIgnoreCase);
            return !installed.AdditionalPatches.IsDefault &&
                installed.AdditionalPatches.Count(patch => patch.SpeciesId == policy.LipidTypeArgument &&
                    string.Equals(patch.Path, policy.NativePatchPath, StringComparison.Ordinal) &&
                    string.Equals(patch.Sha256, policy.NativePatchSha256,
                        StringComparison.OrdinalIgnoreCase)) == 1;
        }
        if (mode == "mapped-lipid21-zenodo-popc")
            return policy.LipidTypeArgument == "POPC" &&
                VerifyHash(policy.NativePatchPath ?? string.Empty,
                    policy.NativePatchSha256 ?? string.Empty) &&
                VerifyHash(policy.NativeSourcePatchPath ?? string.Empty,
                    policy.NativeSourcePatchSha256 ?? string.Empty);
        return false;
    }

    private ApplicablePreparationPolicy? SelectPreparationPolicyLocked(StudyRevision revision,
        AssessedPreparedProtein protein, PlacementProposal proposal, AssessedMembraneModel membrane,
        string? requestedPolicyId = null)
    {
        var applicable = ApplicablePreparationPoliciesLocked(revision, protein, proposal, membrane);
        if (requestedPolicyId is not null)
            return applicable.SingleOrDefault(policy => policy.Id == requestedPolicyId);
        return applicable.SingleOrDefault(policy => policy.Construction.Route == ConstructionRouteKind.PackmolMemgen);
    }

    private ApplicablePreparationPolicy[] ApplicablePreparationPoliciesLocked(StudyRevision revision,
        AssessedPreparedProtein protein, PlacementProposal proposal, AssessedMembraneModel membrane)
    {
        if (_catalogue is null) return Array.Empty<ApplicablePreparationPolicy>();
        var applicable = _catalogue.PreparationPolicies.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty &&
            policy.ApplicableMolecularClass == "canonical-amino-acid-assembly" &&
            policy.MaximumMinimizationIterations == 20000 &&
            double.IsFinite(policy.FinalUnrestrainedRmsForceTargetKjMolNm) &&
            policy.FinalUnrestrainedRmsForceTargetKjMolNm == 10.0 &&
            policy.Construction is { } construction &&
            ExplicitPreparationBoundary.ValidConstructionPolicy(construction) &&
            !string.IsNullOrWhiteSpace(construction.Id) &&
            !string.IsNullOrWhiteSpace(construction.Version) &&
            !construction.EvidenceReferences.IsDefaultOrEmpty &&
            (construction.Route == ConstructionRouteKind.NativeOpenMm
                ? construction.ProviderName == "OpenMM Modeller.addMembrane" &&
                  !string.IsNullOrWhiteSpace(construction.NativePatchPath) &&
                  !string.IsNullOrWhiteSpace(construction.NativePatchSha256) &&
                  construction.LipidTypeArgument is "DLPC" or "DLPE" or "DMPC" or "DOPC" or
                      "DPPC" or "POPC" or "POPE" &&
                  construction.WaterMolarityForIonRounding == 55.4
                : construction.Route == ConstructionRouteKind.PackmolMemgen &&
                  construction.ProviderName == "PACKMOL-Memgen" &&
                  construction.SaltConvention == SaltConventionKind.MemgenChargeCompensated &&
                  construction.NativePatchPath is null &&
                  construction.NativePatchSha256 is null &&
                  construction.LipidTypeArgument is null &&
                  construction.Memgen is not null &&
                  !construction.ProviderAssets.IsDefaultOrEmpty) &&
            !string.IsNullOrWhiteSpace(construction.ProviderVersion) &&
            construction.PositiveIonArgument == "Na+" &&
            construction.NegativeIonArgument == "Cl-" &&
            double.IsFinite(construction.MinimumPaddingNanometers) &&
            construction.MinimumPaddingNanometers > 0 &&
            construction.MaximumAtomCount > 0 &&
            double.IsFinite(construction.MaximumCellDimensionAngstrom) &&
            construction.MaximumCellDimensionAngstrom > 0 &&
            construction.MaximumConstructionSeconds is > 0 and <= 86400 &&
            !string.IsNullOrWhiteSpace(construction.ApproximationStatement) &&
            !construction.CoveredTopologyKinds.IsDefaultOrEmpty &&
            construction.CoveredTopologyKinds.All(kind => Enum.IsDefined(kind)) &&
            !construction.CoveredSpeciesIds.IsDefaultOrEmpty &&
            policy.SystemSettings is { NonbondedMethod: "PME", Constraints: "HBonds" } &&
            policy.SystemSettings.NonbondedCutoffNanometers == 1.0 &&
            policy.SystemSettings.RigidWater &&
            policy.SystemSettings.EwaldErrorTolerance == 0.0005 &&
            policy.SystemSettings.SwitchDistanceNanometers is null &&
            policy.SystemSettings.UseDispersionCorrection &&
            policy.SystemSettings.RemoveCMMotion &&
            policy.SystemSettings.HydrogenMassDaltons is null &&
            construction.CoveredTopologyKinds.Contains(proposal.TopologyKind) &&
            membrane.SpeciesRepresentations.All(item => construction.CoveredSpeciesIds.Contains(item.SpeciesId)) &&
            ExplicitPreparationBoundary.PolicyScopeMatches(policy.Scope, revision, protein, membrane, proposal) &&
            ValidLocalStatePolicy(policy) &&
            ValidStageProteinGeometryPolicy(policy) &&
            !policy.ForceFieldFiles.IsDefaultOrEmpty)
            .ToArray();
        return applicable;
    }

    private ImmutableArray<ConstructionRouteAccount> BuildConstructionRoutesLocked()
    {
        if (_protein is null || _membrane is null || _placement is null)
            return ImmutableArray<ConstructionRouteAccount>.Empty;
        return ApplicablePreparationPoliciesLocked(_study, _protein, _placement.Proposal, _membrane)
            .OrderBy(policy => policy.Construction.Route == ConstructionRouteKind.PackmolMemgen ? 0 : 1)
            .ThenBy(policy => policy.Id, StringComparer.Ordinal)
            .Select(policy =>
            {
                var providerAvailable = ConstructionProviderMatches(policy.Construction,
                    _lastObservedConstructionProvider);
                var assetsAvailable = PreparationAssetsAvailable(policy, _membrane);
                var available = providerAvailable && assetsAvailable;
                return new ConstructionRouteAccount(policy.Id,
                    policy.Construction.Route == ConstructionRouteKind.PackmolMemgen
                        ? "PACKMOL-Memgen: general construction" :
                          $"OpenMM native: {policy.Construction.LipidTypeArgument}",
                    policy.Construction.Route, policy.Construction.SaltConvention, available,
                    available ? null : !providerAvailable
                        ? "The exact local provider is unavailable or changed. Restore it and restart the local app to recheck this route."
                        : "One or more exact construction assets are unavailable or changed.",
                    policy.Construction.ProviderVersion,
                    policy.Construction.MaximumConstructionSeconds,
                    policy.Construction.MaximumAtomCount,
                    policy.Construction.MaximumCellDimensionAngstrom);
            }).ToImmutableArray();
    }

    private bool PreparationAssetsAvailable(ApplicablePreparationPolicy policy,
        AssessedMembraneModel membrane) =>
        policy.ForceFieldFiles.All(VerifiedAsset) &&
        new[] { policy.Water, policy.Sodium, policy.Chloride }.All(VerifiedRepresentation) &&
        !membrane.SpeciesRepresentations.IsDefaultOrEmpty &&
        membrane.SpeciesRepresentations.All(VerifiedRepresentation);

    private static bool ValidLocalStatePolicy(ApplicablePreparationPolicy policy)
    {
        var observation = policy.LocalStateObservation;
        var criteria = policy.ConstructionCriteria;
        var contacts = policy.ContactCriteria;
        if (observation is null || observation.ContactRolePairs.IsDefaultOrEmpty ||
            observation.RequiredMetricNames.IsDefaultOrEmpty || criteria.IsDefaultOrEmpty ||
            contacts.IsDefault ||
            observation.AtomRadiusByElementAngstrom is not { Count: > 0 } ||
            observation.AtomRadiusByElementAngstrom.Any(item =>
                string.IsNullOrWhiteSpace(item.Key) || !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(observation.ContactSearchRadiusAngstrom) ||
            observation.ContactSearchRadiusAngstrom <= 0 ||
            observation.MaximumReportedPairs <= 0 || !observation.UsePeriodicBoundary ||
            (observation.ReferenceMidplaneZAngstrom is double referenceMidplane &&
                !double.IsFinite(referenceMidplane)) ||
            observation.ContactRolePairs.Any(pair => !Enum.IsDefined(pair.FirstMoleculeRole) ||
                !Enum.IsDefined(pair.SecondMoleculeRole)) ||
            observation.ContactRolePairs.Distinct().Count() != observation.ContactRolePairs.Length ||
            observation.RequiredMetricNames.Any(string.IsNullOrWhiteSpace) ||
            observation.RequiredMetricNames.Distinct(StringComparer.Ordinal).Count() !=
                observation.RequiredMetricNames.Length ||
            observation.RequiredMetricNames.Any(name => name is not (
                "minimumIntermolecularDistanceAngstrom" or
                "minimumIntermolecularHeavyAtomDistanceAngstrom" or "upperLipidHeadMeanZAngstrom" or
                "lowerLipidHeadMeanZAngstrom" or "leafletHeadSeparationAngstrom" or
                "proteinBilayerMidplaneOffsetAngstrom")) ||
            criteria.Select(criterion => criterion.MeasurementName).Distinct(StringComparer.Ordinal).Count() !=
                criteria.Length ||
            criteria.Any(criterion => criterion.MeasurementName is
                "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom") ||
            (!policy.AssessmentCriteria.IsDefault && policy.AssessmentCriteria.Any(criterion =>
                criterion.MeasurementName is
                    "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom")) ||
            !criteria.Any(criterion => criterion.MeasurementName == "minimumIntermolecularHeavyAtomDistanceAngstrom" &&
                criterion.Minimum is double distance && double.IsFinite(distance) && distance == 1.5) ||
            contacts.Select(item => (item.StageKind, item.FirstMoleculeRole, item.SecondMoleculeRole))
                .Distinct().Count() != contacts.Length ||
            contacts.Any(item => item.StageKind is not
                    (null or StageKind.Minimization or StageKind.Equilibration) ||
                item.MinimumPairsWithinSearchRadius < 0 ||
                !observation.ContactRolePairs.Any(pair =>
                    pair.FirstMoleculeRole == item.FirstMoleculeRole &&
                    pair.SecondMoleculeRole == item.SecondMoleculeRole) ||
                (item.MinimumNearestDistanceAngstrom is double lower &&
                 (!double.IsFinite(lower) || lower <= 0)) ||
                (item.MaximumNearestDistanceAngstrom is double upper &&
                 (!double.IsFinite(upper) || upper <= 0 || upper > observation.ContactSearchRadiusAngstrom)) ||
                (item.MinimumNearestDistanceAngstrom is double minimumContact &&
                 item.MaximumNearestDistanceAngstrom is double maximumContact &&
                 minimumContact > maximumContact)) ||
            contacts.Any(item => item.StageKind == StageKind.Equilibration &&
                policy.OptionalEquilibration is null) ||
            criteria.Any(criterion => string.IsNullOrWhiteSpace(criterion.MeasurementName) ||
                !observation.RequiredMetricNames.Contains(criterion.MeasurementName) ||
                string.IsNullOrWhiteSpace(criterion.Unit) || string.IsNullOrWhiteSpace(criterion.Scope) ||
                (criterion.Minimum is double minimum && !double.IsFinite(minimum)) ||
                (criterion.Maximum is double maximum && !double.IsFinite(maximum)) ||
                (criterion.Minimum is double low && criterion.Maximum is double high && low > high))) return false;

        if (!observation.RequiredMetricNames.Contains("leafletHeadSeparationAngstrom") ||
            !observation.RequiredMetricNames.Contains("proteinBilayerMidplaneOffsetAngstrom"))
            return false;
        // Explicitly empty assessment criteria keep the result observed and
        // exportable without asserting biological suitability.
        return !policy.AssessmentCriteria.IsDefault;
    }

    private static bool ValidStageProteinGeometryPolicy(ApplicablePreparationPolicy policy)
    {
        var measurement = policy.StageProteinGeometryMeasurement;
        var criteria = policy.StageProteinGeometryCriteria;
        if (measurement is null || measurement.RequiredKinds.IsDefaultOrEmpty || criteria.IsDefaultOrEmpty ||
            measurement.RequiredKinds.Any(string.IsNullOrWhiteSpace) ||
            measurement.RequiredKinds.Distinct(StringComparer.Ordinal).Count() !=
                measurement.RequiredKinds.Length ||
            !measurement.RequiredKinds.ToHashSet(StringComparer.Ordinal).SetEquals(
                new[] { "covalentBond", "chainContinuity", "nonbondedDistance" }) ||
            measurement.AtomRadiusByElementAngstrom is not { Count: > 0 } ||
            measurement.AtomRadiusByElementAngstrom.Any(item =>
                string.IsNullOrWhiteSpace(item.Key) || !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(measurement.NeighborSearchRadiusAngstrom) ||
            measurement.NeighborSearchRadiusAngstrom <= 0 || measurement.ExcludedBondHops < 0 ||
            measurement.MaximumReportedPairs <= 0 ||
            criteria.Select(item => (item.StageKind, item.Criterion?.Kind)).Distinct().Count() != criteria.Length ||
            criteria.Any(item => item.Criterion is null ||
                item.StageKind is not (StageKind.Minimization or StageKind.Equilibration) ||
                !measurement.RequiredKinds.Contains(item.Criterion.Kind) ||
                (item.Criterion.Kind == "covalentBond" && item.Criterion.AllowNotApplicable) ||
                item.Criterion.MinimumObservedAngstrom is null &&
                    item.Criterion.MaximumObservedAngstrom is null ||
                (item.Criterion.MinimumObservedAngstrom is double minimum && !double.IsFinite(minimum)) ||
                (item.Criterion.MaximumObservedAngstrom is double maximum && !double.IsFinite(maximum)) ||
                (item.Criterion.MinimumObservedAngstrom is double lower &&
                 item.Criterion.MaximumObservedAngstrom is double upper && lower > upper))) return false;
        return new[] { StageKind.Minimization, StageKind.Equilibration }
            .Where(kind => kind != StageKind.Equilibration || policy.OptionalEquilibration is not null)
            .All(kind => measurement.RequiredKinds.All(required =>
                criteria.Count(item => item.StageKind == kind && item.Criterion.Kind == required) == 1));
    }

    private bool CanRequestEquilibrationLocked(CompletedStage minimized)
    {
        var attempt = minimized.Attempt;
        return minimized.Kind == StageKind.Minimization &&
            minimized.Observation.Kind == StageKind.Minimization &&
            minimized.Observation.StageId == minimized.Id &&
            minimized.Observation.AttemptId == attempt.Id &&
            minimized.Molecule.Id == minimized.Id &&
            minimized.Correspondence.ResultId == minimized.Id &&
            minimized.Correspondence.Complete &&
            _stageAssessments.TryGetValue(minimized.Id, out var sourceAssessment) &&
            sourceAssessment.StageId == minimized.Id &&
            sourceAssessment.CurrentlyApplicable &&
            sourceAssessment.CheckStanding != PreparationCheckStanding.IssuesFound &&
            _constructedByAttempt.TryGetValue(attempt.Id, out var constructed) &&
            constructed.Attempt.Id == attempt.Id &&
            constructed.Attempt.StudyRevisionId == attempt.StudyRevisionId &&
            constructed.Attempt.PolicyFingerprintSha256 == attempt.PolicyFingerprintSha256 &&
            _revisions.TryGetValue(attempt.StudyRevisionId, out var revision) &&
            _proteinByAttempt.TryGetValue(attempt.Id, out var protein) &&
            protein.Id == attempt.ProteinId &&
            _membraneByAttempt.TryGetValue(attempt.Id, out var membrane) &&
            membrane.Id == attempt.MembraneId &&
            _placementByAttempt.TryGetValue(attempt.Id, out var placement) &&
            placement.Id == attempt.PlacementId &&
            _policyByAttempt.TryGetValue(attempt.Id, out var policy) &&
            minimized.PolicyId == policy.Id &&
            PreparationPolicyFingerprint.Matches(attempt, policy) &&
            ExplicitPreparationBoundary.PolicyScopeMatches(policy.Scope, revision, protein,
                membrane, placement.Proposal) &&
            QualifiedEquilibrationLocked(policy, attempt.Id);
    }

    private bool QualifiedEquilibrationLocked(ApplicablePreparationPolicy policy, string attemptId)
    {
        if (policy.OptionalEquilibration is not { } protocol || _catalogue is null ||
            string.IsNullOrWhiteSpace(protocol.Id) || protocol.Stages.IsDefaultOrEmpty ||
            protocol.ExtensionWindow is null || protocol.MaximumExtensions < 0 ||
            protocol.MaximumSampleCount <= 0 ||
            protocol.MaximumFrameBytes <= 0 || protocol.MaximumFrameBytes > 16L * 1024 * 1024 * 1024 ||
            protocol.RequiredObservations.IsDefaultOrEmpty ||
            protocol.Observables.IsDefaultOrEmpty || protocol.SufficiencyRules.IsDefaultOrEmpty ||
            string.IsNullOrWhiteSpace(protocol.ComparisonBasis) ||
            !_membraneByAttempt.TryGetValue(attemptId, out var membrane) ||
            !_placementByAttempt.TryGetValue(attemptId, out var placement)) return false;
        var chosenSpecies = membrane.Intended.Upper.Fractions.Concat(membrane.Intended.Lower.Fractions)
            .Where(fraction => fraction.Fraction > 0).Select(fraction => fraction.SpeciesId)
            .Distinct(StringComparer.Ordinal).ToArray();
        var mixed = membrane.Intended.Upper.Fractions.Count(fraction => fraction.Fraction > 0) > 1 ||
            membrane.Intended.Lower.Fractions.Count(fraction => fraction.Fraction > 0) > 1;
        var asymmetric = !membrane.Intended.Upper.Fractions.SequenceEqual(membrane.Intended.Lower.Fractions);
        if (!EquilibrationProtocolFingerprint.TryCompute(protocol, out var protocolSha256)) return false;
        var matches = _catalogue.EquilibrationQualifications.Where(record =>
            record.PolicyId == policy.Id && record.PolicyVersion == policy.Version &&
            record.ProtocolId == protocol.Id &&
            string.Equals(record.ProtocolSha256, protocolSha256, StringComparison.OrdinalIgnoreCase) &&
            !string.IsNullOrWhiteSpace(record.Version) && !record.EvidenceReferences.IsDefaultOrEmpty &&
            !record.CoveredTopologyKinds.IsDefaultOrEmpty &&
            record.CoveredTopologyKinds.All(kind => Enum.IsDefined(kind)) &&
            record.CoveredTopologyKinds.Contains(placement.Proposal.TopologyKind) &&
            chosenSpecies.All(record.CoveredSpeciesIds.Contains) &&
            (!mixed || record.AllowsMixtures) && (!asymmetric || record.AllowsAsymmetry)).Take(2).ToArray();
        return matches.Length == 1;
    }

    private IReadOnlyDictionary<string, MolecularRepresentation> LipidsLocked() =>
        (_catalogue?.Lipids ?? ImmutableArray<MolecularRepresentation>.Empty)
            .Where(item => VerifiedRepresentation(item))
            .GroupBy(item => item.SpeciesId, StringComparer.Ordinal)
            .Where(group => group.Count() == 1)
            .ToDictionary(group => group.Key, group => group.Single(), StringComparer.Ordinal);

    private bool VerifiedRepresentation(MolecularRepresentation item) =>
        !string.IsNullOrWhiteSpace(item.SpeciesId) && item.AtomCount > 0 &&
        !string.IsNullOrWhiteSpace(item.ForceFieldFamily) &&
        !string.IsNullOrWhiteSpace(item.ForceFieldVersion) &&
        VerifyHash(item.TemplatePath, item.TemplateSha256) &&
        VerifyHash(item.CoordinateTemplatePath, item.CoordinateTemplateSha256);

    private static bool VerifiedAsset(ForceFieldAsset asset) =>
        !string.IsNullOrWhiteSpace(asset.Id) && !string.IsNullOrWhiteSpace(asset.Version) &&
        !string.IsNullOrWhiteSpace(asset.Family) && VerifyHash(asset.Path, asset.Sha256);

    private bool ValidSourceEvidenceIdentity(StructuralSource source)
    {
        if (source.Kind == SourceRouteKind.Upload)
            return source.Prediction is null && source.UploadProvenance is { } origin &&
                Enum.IsDefined(origin);
        if (source.UploadProvenance is not null || source.UploadProvenanceNote is not null) return false;
        if (source.Kind != SourceRouteKind.AlphaFold) return source.Kind == SourceRouteKind.Rcsb && source.Prediction is null;
        var prediction = source.Prediction;
        if (prediction is null || string.IsNullOrWhiteSpace(prediction.RecordId) ||
            source.Id != $"alphafold:{prediction.RecordId}" ||
            string.IsNullOrWhiteSpace(prediction.CoordinateUrl) ||
            prediction.CoordinateSha256 != source.Sha256) return false;
        return prediction.PaeStanding switch
        {
            PaeAcquisitionStanding.Available => prediction.PaePath is { } path && IsWorkspaceFile(path) &&
                VerifyHash(path, prediction.PaeSha256 ?? string.Empty) &&
                !string.IsNullOrWhiteSpace(prediction.PaeUrl),
            PaeAcquisitionStanding.Unavailable => prediction.PaePath is null && prediction.PaeSha256 is null &&
                !string.IsNullOrWhiteSpace(prediction.PaeReason),
            _ => false
        };
    }

    private static bool VerifyHash(string path, string expected)
    {
        if (string.IsNullOrWhiteSpace(path) || string.IsNullOrWhiteSpace(expected) || !File.Exists(path)) return false;
        try { return Hash(path).Equals(expected, StringComparison.OrdinalIgnoreCase); }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException) { return false; }
    }

    private static (CatalogueDocument?, string?) ReadCatalogue(string path)
    {
        if (string.IsNullOrWhiteSpace(path) || !File.Exists(path))
            return (null, "No qualified local policy catalogue is installed. Scientific routes that require one remain unavailable.");
        try
        {
            var options = new JsonSerializerOptions(JsonSerializerDefaults.Web);
            options.Converters.Add(new JsonStringEnumConverter());
            var document = JsonSerializer.Deserialize<CatalogueDocument>(File.ReadAllText(path), options);
            if (document is null || string.IsNullOrWhiteSpace(document.Version) ||
                document.EvidenceReferences.IsDefaultOrEmpty || document.Lipids.IsDefault ||
                document.ProteinChemicalStates.IsDefault || document.ProteinStructuralPolicies.IsDefault ||
                document.MembranePolicies.IsDefault ||
                document.PlacementPolicies.IsDefault || document.PreparationPolicies.IsDefault)
                return (null, "The local policy catalogue lacks its version, evidence or asset identities.");
            var baseDirectory = Path.GetDirectoryName(Path.GetFullPath(path))!;
            string Resolve(string source) => Path.GetFullPath(Path.IsPathRooted(source) ? source : Path.Combine(baseDirectory, source));
            ForceFieldAsset MapAsset(ForceFieldAsset asset) => asset with { Path = Resolve(asset.Path) };
            MolecularRepresentation Map(MolecularRepresentation item) => item with
            {
                TemplatePath = Resolve(item.TemplatePath),
                CoordinateTemplatePath = Resolve(item.CoordinateTemplatePath),
                StereoChecks = item.StereoChecks.IsDefault
                    ? ImmutableArray<MolecularStereoCheck>.Empty : item.StereoChecks
            };
            ProteinChemicalStatePolicy MapChemical(ProteinChemicalStatePolicy item) => item with
            { ForceFieldFiles = item.ForceFieldFiles.Select(MapAsset).ToImmutableArray() };
            ApplicablePreparationPolicy MapPreparation(ApplicablePreparationPolicy item) => item with
            {
                ForceFieldFiles = item.ForceFieldFiles.Select(MapAsset).ToImmutableArray(),
                Construction = item.Construction with
                {
                    NativePatchPath = item.Construction.NativePatchPath is null ?
                        null : Resolve(item.Construction.NativePatchPath),
                    NativeSourcePatchPath = item.Construction.NativeSourcePatchPath is null ?
                        null : Resolve(item.Construction.NativeSourcePatchPath),
                    ProviderAssets = item.Construction.ProviderAssets.IsDefault
                        ? ImmutableArray<ProviderAsset>.Empty
                        : item.Construction.ProviderAssets.Select(asset => asset with
                            { Path = Resolve(asset.Path) }).ToImmutableArray()
                },
                Water = Map(item.Water), Sodium = Map(item.Sodium), Chloride = Map(item.Chloride)
            };
            return (document with
            {
                Lipids = document.Lipids.Select(Map).ToImmutableArray(),
                ProteinChemicalStates = document.ProteinChemicalStates.Select(MapChemical).ToImmutableArray(),
                PreparationPolicies = document.PreparationPolicies.Select(MapPreparation).ToImmutableArray(),
                PlacementWitnesses = document.PlacementWitnesses.IsDefault
                    ? ImmutableArray<PlacementStructuralWitness>.Empty
                    : document.PlacementWitnesses,
                PpmResidueLibraryPath = string.IsNullOrWhiteSpace(document.PpmResidueLibraryPath)
                    ? string.Empty : Resolve(document.PpmResidueLibraryPath),
                EquilibrationQualifications = document.EquilibrationQualifications.IsDefault
                    ? ImmutableArray<EquilibrationQualification>.Empty
                    : document.EquilibrationQualifications
            }, null);
        }
        catch (Exception exception) when (exception is IOException or JsonException or ArgumentException or
                                         InvalidOperationException or NullReferenceException or UnauthorizedAccessException)
        {
            return (null, "The local policy catalogue could not be read coherently: " + exception.Message);
        }
    }

    private sealed record CatalogueDocument(
        string Version,
        ImmutableArray<string> EvidenceReferences,
        ImmutableArray<MolecularRepresentation> Lipids,
        ImmutableArray<ProteinChemicalStatePolicy> ProteinChemicalStates,
        ImmutableArray<ProteinStructuralAssessmentPolicy> ProteinStructuralPolicies,
        ImmutableArray<MembraneSupportPolicy> MembranePolicies,
        ImmutableArray<PlacementSupportPolicy> PlacementPolicies,
        ImmutableArray<PlacementStructuralWitness> PlacementWitnesses,
        ImmutableArray<ApplicablePreparationPolicy> PreparationPolicies,
        ImmutableArray<EquilibrationQualification> EquilibrationQualifications,
        string PpmVersion,
        string PpmExecutableSha256,
        int MaximumSourceAtoms,
        string PpmResidueLibraryPath = "",
        string PpmResidueLibrarySha256 = "");

    private sealed record EquilibrationQualification(
        string PolicyId,
        string PolicyVersion,
        string ProtocolId,
        string ProtocolSha256,
        string Version,
        ImmutableArray<string> EvidenceReferences,
        ImmutableArray<string> CoveredSpeciesIds,
        ImmutableArray<ProteinTopologyKind> CoveredTopologyKinds,
        bool AllowsMixtures,
        bool AllowsAsymmetry);

    private string WorkDirectory(string region, string id)
    {
        var directory = Path.GetFullPath(Path.Combine(_workspaceRoot, region, id));
        if (!directory.StartsWith(_workspaceRoot + Path.DirectorySeparatorChar, StringComparison.Ordinal))
            throw new InvalidOperationException("The requested work directory is outside the local workspace.");
        Directory.CreateDirectory(directory);
        return directory;
    }

    private bool IsWorkspaceFile(string path)
    {
        if (string.IsNullOrWhiteSpace(path)) return false;
        var resolved = Path.GetFullPath(path);
        if (!resolved.StartsWith(_workspaceRoot + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase) ||
            !File.Exists(resolved)) return false;
        try
        {
            var current = _workspaceRoot;
            foreach (var segment in Path.GetRelativePath(_workspaceRoot, resolved)
                         .Split(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar))
            {
                current = Path.Combine(current, segment);
                if ((File.GetAttributes(current) & FileAttributes.ReparsePoint) != 0) return false;
            }
            return true;
        }
        catch (Exception exception) when (exception is IOException or UnauthorizedAccessException) { return false; }
    }

    private string StructureUrlLocked(string path)
    {
        if (!IsWorkspaceFile(path)) throw new InvalidOperationException("The selected structure is not a workspace artifact.");
        var sha256 = Hash(path);
        // Different proposal subjects can inspect one immutable selected-coordinate
        // preview. Reuse its verified token so changing evidence does not reload
        // the same molecular geometry or lose the researcher's camera.
        var token = _structures.FirstOrDefault(item =>
            item.Value.Path == path && item.Value.Sha256 == sha256).Key;
        if (token is null)
        {
            token = NewId();
            _structures[token] = new StructureBinding(path, sha256, Path.GetExtension(path).ToLowerInvariant());
        }
        var format = Path.GetExtension(path).ToLowerInvariant() is ".cif" or ".mmcif" ? "mmcif" : "pdb";
        return "/api/structures/" + token + "?format=" + format;
    }

    private ConstructionTrialSummary BindTrialDiagnosticsLocked(string attemptId,
        ConstructionTrialSummary trial)
    {
        if (trial.DiagnosticArtifacts.IsDefaultOrEmpty) return trial with
        {
            DiagnosticArtifacts = ImmutableArray<TrialDiagnosticArtifact>.Empty
        };
        var bound = ImmutableArray.CreateBuilder<TrialDiagnosticArtifact>();
        foreach (var artifact in trial.DiagnosticArtifacts)
        {
            if (artifact.LocalPath is not { } path || !IsWorkspaceFile(path) ||
                !VerifyHash(path, artifact.Sha256) ||
                Path.GetFileName(path) != artifact.FileName) continue;
            var existing = _diagnostics.FirstOrDefault(item =>
                item.Value.AttemptId == attemptId && item.Value.TrialId == trial.TrialId &&
                item.Value.Role == artifact.Role && item.Value.Sha256 == artifact.Sha256 &&
                item.Value.Path == path);
            var token = existing.Key ?? NewId();
            _diagnostics[token] = new DiagnosticBinding(attemptId,
                _currentAttempt?.StudyRevisionId ?? _study.Id, trial.TrialId, artifact.Role,
                path, artifact.Sha256, artifact.FileName);
            var structural = Path.GetExtension(path).ToLowerInvariant() is ".pdb" or ".cif" or ".mmcif";
            bound.Add(artifact with
            {
                DownloadUrl = "/api/diagnostics/" + token,
                StructureUrl = structural ? StructureUrlLocked(path) : null,
                SubjectId = structural ? token : null
            });
        }
        return trial with { DiagnosticArtifacts = bound.ToImmutable() };
    }

    private bool SelectedStructureIsVerifiedLocked(string subjectId, string studyRevisionId,
        string expectedPath, string expectedSha256)
    {
        var selected = _inspection.Current;
        if (selected?.SubjectId != subjectId || selected.StudyRevisionId != studyRevisionId ||
            selected.StructureUrl is null) return false;
        const string prefix = "/api/structures/";
        var structurePath = selected.StructureUrl.Split('?', 2)[0];
        if (!structurePath.StartsWith(prefix, StringComparison.Ordinal)) return false;
        var token = structurePath[prefix.Length..];
        return _structures.TryGetValue(token, out var binding) &&
            string.Equals(binding.Path, expectedPath, StringComparison.Ordinal) &&
            string.Equals(binding.Sha256, expectedSha256, StringComparison.OrdinalIgnoreCase) &&
            IsWorkspaceFile(binding.Path) && VerifyHash(binding.Path, binding.Sha256);
    }

    private sealed record StructureBinding(string Path, string Sha256, string Extension);
    private sealed record DiagnosticBinding(string AttemptId, string StudyRevisionId,
        string TrialId, string Role, string Path, string Sha256, string FileName);

    private static string Hash(string path)
    {
        using var input = File.OpenRead(path);
        return Convert.ToHexString(SHA256.HashData(input)).ToLowerInvariant();
    }

    private static string NewId() => Guid.NewGuid().ToString("N");
    private void TouchLocked() => _revision++;

    private static string? Text(JsonElement data, string name) =>
        data.TryGetProperty(name, out var value) && value.ValueKind == JsonValueKind.String &&
        !string.IsNullOrWhiteSpace(value.GetString()) ? value.GetString()!.Trim() : null;
    private static ProteinTopologyKind? ParseProteinTopology(string? value) => value switch
    {
        "membrane-spanning" => ProteinTopologyKind.MembraneSpanning,
        "one-surface-associated" => ProteinTopologyKind.OneSurfaceAssociated,
        _ => null
    };
    private static PlacementPhysicalSide? ParsePlacementSide(string? value) => value switch
    {
        "upper" => PlacementPhysicalSide.Upper,
        "lower" => PlacementPhysicalSide.Lower,
        "both" => PlacementPhysicalSide.Both,
        _ => null
    };
    private static PpmNterminalSide? ParsePpmNterminalSide(string? value) => value switch
    {
        "in" => PpmNterminalSide.In,
        "out" => PpmNterminalSide.Out,
        _ => null
    };
    private static int? Integer(JsonElement data, string name) =>
        data.TryGetProperty(name, out var value) && value.TryGetInt32(out var number) ? number : null;
    private static double? Number(JsonElement data, string name) =>
        data.TryGetProperty(name, out var value) && value.TryGetDouble(out var number) ? number : null;
    private static bool? Boolean(JsonElement data, string name) =>
        data.TryGetProperty(name, out var value) && value.ValueKind is JsonValueKind.True or JsonValueKind.False
            ? value.GetBoolean() : null;
    private static ImmutableArray<T> ParseArray<T>(JsonElement data, string name) =>
        data.TryGetProperty(name, out var value) && value.ValueKind == JsonValueKind.Array
            ? JsonSerializer.Deserialize<ImmutableArray<T>>(value.GetRawText(), new JsonSerializerOptions(JsonSerializerDefaults.Web))
            : ImmutableArray<T>.Empty;
    private static bool Coherent(ImmutableArray<LipidFraction> fractions) =>
        !fractions.IsDefaultOrEmpty && fractions.All(item => !string.IsNullOrWhiteSpace(item.SpeciesId) &&
            double.IsFinite(item.Fraction) && item.Fraction >= 0) &&
        fractions.Select(item => item.SpeciesId).Distinct(StringComparer.Ordinal).Count() == fractions.Length &&
        Math.Abs(fractions.Sum(item => item.Fraction) - 1.0) <= 1e-9;

    private static bool SameComposition(ImmutableArray<LipidFraction> first,
        ImmutableArray<LipidFraction> second) =>
        first.Length == second.Length && first.OrderBy(item => item.SpeciesId, StringComparer.Ordinal)
            .Zip(second.OrderBy(item => item.SpeciesId, StringComparer.Ordinal),
                (left, right) => left.SpeciesId == right.SpeciesId &&
                    Math.Abs(left.Fraction - right.Fraction) <= 1e-12).All(equal => equal);
}
