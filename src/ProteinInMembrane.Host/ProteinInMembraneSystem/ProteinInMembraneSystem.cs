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

/// <summary>
/// Owns the study, decisions, and handoffs among direct children. Provider
/// observations are never promoted here merely because a process succeeded.
/// </summary>
public sealed class ProteinInMembraneSystem
{
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
    private readonly string _packmolExecutablePath;
    private readonly CatalogueDocument? _catalogue;
    private readonly string? _catalogueIssue;
    private readonly Dictionary<string, string> _structures = new(StringComparer.Ordinal);
    private readonly Dictionary<string, StructuralSource> _uploads = new(StringComparer.Ordinal);
    private readonly Dictionary<string, StudyRevision> _revisions = new(StringComparer.Ordinal);
    private readonly Dictionary<string, CompletedStage> _stages = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ConstructedExplicitSystem> _constructedByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ApplicablePreparationPolicy> _policyByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedPreparedProtein> _proteinByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedMembraneModel> _membraneByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, AssessedProteinMembranePlacement> _placementByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ImmutableArray<ResearcherDecision>> _decisionsByAttempt = new(StringComparer.Ordinal);
    private readonly Dictionary<string, PreparationAssessmentResult> _stageAssessments = new(StringComparer.Ordinal);
    private readonly Dictionary<string, ImmutableArray<ScientificFinding>> _laterFindingsByStage = new(StringComparer.Ordinal);
    private readonly Dictionary<string, CompletedStageBundle> _bundles = new(StringComparer.Ordinal);
    private readonly List<WorkspaceNotice> _notices = new();
    private ImmutableArray<CandidateSourceRecord> _candidates = ImmutableArray<CandidateSourceRecord>.Empty;
    private StudyRevision _study;
    private StructuralSource? _selectedSource;
    private SourceInspectionReport? _sourceInspection;
    private PreparationProposalReport? _preparationProposals;
    private ImmutableArray<ResearcherDecision> _decisions = ImmutableArray<ResearcherDecision>.Empty;
    private ImmutableArray<ResearcherDecision> _allDecisions = ImmutableArray<ResearcherDecision>.Empty;
    private AssessedPreparedProtein? _protein;
    private ProteinPreparationDiagnostic? _proteinDiagnostic;
    private MembraneModel? _membraneProposal;
    private AssessedMembraneModel? _membrane;
    private PlacementProposal? _placementProposal;
    private OpmReferenceReview? _opmReview;
    private PlacementMeasurementReport? _placementMeasurement;
    private PlacementStructuralWitness? _placementWitness;
    private AssessedProteinMembranePlacement? _placement;
    private PreparationAttempt? _currentAttempt;
    private StageExecutionState? _execution;
    private CancellationTokenSource? _attemptStop;
    private Task? _attemptTask;
    private long _revision;

    public event Action? Changed;

    public ProteinInMembraneSystem(
        IScientificWorkerExchange worker,
        ExternalSourceExchange sources,
        string workspaceRoot,
        string policyCataloguePath,
        string ppmExecutablePath,
        string packmolExecutablePath)
    {
        _sources = sources;
        _proteinPreparation = new ProteinPreparationBoundary(worker);
        _membraneAssessment = new MembraneAssessmentBoundary(worker);
        _placementAssessment = new PlacementAssessmentBoundary(worker);
        _explicitPreparation = new ExplicitPreparationBoundary(worker, worker, worker);
        _export = new ExportBoundary(worker);
        _workspaceRoot = Path.GetFullPath(workspaceRoot);
        Directory.CreateDirectory(_workspaceRoot);
        _ppmExecutablePath = ppmExecutablePath;
        _packmolExecutablePath = packmolExecutablePath;
        (_catalogue, _catalogueIssue) = ReadCatalogue(policyCataloguePath);
        _study = new StudyRevision(NewId(), 1, null, null, null, FixedStudyConditions.Initial);
        _revisions.Add(_study.Id, _study);
        _workspace.RetainStudy(_study);
        if (_catalogueIssue is not null) Notice("warning", _catalogueIssue, null);
    }

    public WorkspaceState Snapshot()
    {
        lock (_gate) return SnapshotLocked();
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
                ActorActionKind.ProposeMembrane => ProposeMembrane(command.Data),
                ActorActionKind.AdoptMembrane => await AdoptMembraneAsync(command.Data, cancellationToken),
                ActorActionKind.ProposePlacement => await ProposePlacementAsync(command.Data, cancellationToken),
                ActorActionKind.RevisePlacement => await RevisePlacementAsync(command.Data, cancellationToken),
                ActorActionKind.AdoptPlacement => AdoptPlacement(command.Data),
                ActorActionKind.StartPreparation => StartPreparation(),
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
            return _structures.TryGetValue(token, out var path) && IsWorkspaceFile(path) ? path : null;
    }

    public string? ExportPath(string stageId)
    {
        lock (_gate)
            return _bundles.TryGetValue(stageId, out var bundle) &&
                _stageAssessments.TryGetValue(stageId, out var assessment) && assessment.CurrentlyApplicable &&
                IsWorkspaceFile(bundle.BundlePath)
                ? bundle.BundlePath : null;
    }

    private async Task<string?> SearchSourceAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var query = Text(data, "query");
        if (query is null) return "Enter a structural identifier or descriptive protein query.";
        var result = await _sources.SearchAsync(query, cancellationToken);
        lock (_gate)
        {
            _candidates = result.Candidates;
            foreach (var issue in result.UnavailableRoutes) Notice("warning", issue, null);
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
            AdvanceStudyLocked(null, _study.Membrane, null);
            var carriedMembrane = CarryMembraneLocked(_study);
            _preparationProposals = null;
            _decisions = ImmutableArray<ResearcherDecision>.Empty;
            _protein = null;
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementWitness = null;
            _placement = null;
            _membrane = carriedMembrane;
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
        // An assembly CopyId is the full observed output chain ID (for example A1),
        // not a suffix or an independently invented copy number.
        if (chains.IsDefaultOrEmpty || chains.Any(item => string.IsNullOrWhiteSpace(item.CopyId) ||
                !model.Chains.Any(observed => observed.Name == item.SourceChain)) ||
            chains.Select(item => item.CopyId).Distinct(StringComparer.Ordinal).Count() != chains.Length ||
            selectedAssembly is not null && chains.Any(item => !selectedAssembly.ChainCopies.Contains(item)) ||
            selectedAssembly is null && chains.Any(item => item.CopyId != item.SourceChain) ||
            partners.Length != model.Partners.Length ||
            partners.Select(item => item.SourceId).Distinct(StringComparer.Ordinal).Count() != partners.Length ||
            partners.Any(item => string.IsNullOrWhiteSpace(item.Reason) ||
                !model.Partners.Any(observed => observed.SourceId == item.SourceId)) ||
            rawAltlocs.Any(choice => choice.Residue.CopyId != string.Empty ||
                choice.Residue.Model != model.Index ||
                !model.Residues.Any(residue => residue.Address == choice.Residue &&
                    residue.AlternateLocations.Contains(choice.Altloc))) ||
            rawAltlocs.Select(choice => choice.Residue).Distinct().Count() != rawAltlocs.Length)
            return "Choose observed chains and explicit dispositions for the selected source's partners.";
        var altlocs = rawAltlocs.Select(item => item with { DecisionId = NewId() }).ToImmutableArray();
        var intended = new IntendedProteinModel(NewId(), source, modelIndex!.Value, assembly, chains, partners, altlocs);
        StudyRevision revision;
        ProteinChemicalStatePolicy? chemicalPolicy;
        ProteinStructuralAssessmentPolicy? structuralPolicy;
        lock (_gate)
        {
            AdvanceStudyLocked(intended, _study.Membrane, null);
            revision = _study;
            chemicalPolicy = SelectChemicalPolicyLocked(model, intended);
            structuralPolicy = SelectStructuralPolicyLocked(model, intended);
            var carriedMembrane = CarryMembraneLocked(revision);
            _preparationProposals = null;
            _decisions = ImmutableArray<ResearcherDecision>.Empty;
            _protein = null;
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementWitness = null;
            _placement = null;
            _membrane = carriedMembrane;
            TouchLocked();
        }
        Changed?.Invoke();
        if (chemicalPolicy is null || structuralPolicy is null)
        {
            AddNotice("warning", "No qualified protein chemical-state and structural assessment policies cover this selected structure; preparation remains unavailable.", intended.Id);
            return null;
        }
        var proposed = await _proteinPreparation.ProposeChangesAsync(revision, intended, inspection, chemicalPolicy,
            structuralPolicy,
            WorkDirectory("protein-proposals", revision.Id), cancellationToken);
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _preparationProposals = proposed.Value;
            if (proposed.Value is null) Notice("warning", proposed.Reason, intended.Id);
            TouchLocked();
        }
        Changed?.Invoke();
        if (proposed.Value is { Changes.Length: 0, UnresolvedQuestions.Length: 0 })
            await TryPrepareProteinAsync(revision, cancellationToken);
        return null;
    }

    private async Task<string?> DecidePreparationChangeAsync(JsonElement data, CancellationToken cancellationToken)
    {
        var proposalId = Text(data, "proposalId");
        var approved = Boolean(data, "approve");
        var rationale = Text(data, "rationale") ?? string.Empty;
        PreparationChangeProposal? proposal;
        StudyRevision revision;
        lock (_gate)
        {
            proposal = _preparationProposals?.Changes.FirstOrDefault(item => item.Id == proposalId);
            revision = _study;
            if (proposal is null || proposal.StudyRevisionId != revision.Id) return "This preparation proposal is absent or stale.";
            if (approved is null) return "The proposal decision must explicitly approve or decline.";
            if (_protein is not null) return "The prepared protein has already been established from the reviewed choices.";
            if (_decisions.Any(item => item.SubjectId == proposal.Id && item.StudyRevisionId == revision.Id))
                return "This exact proposal has already been decided.";
            if (approved.Value && (proposal.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState) &&
                _preparationProposals!.Changes.Any(other => other.Id != proposal.Id &&
                    other.Kind == proposal.Kind && other.Residue == proposal.Residue &&
                    _decisions.Any(decision => decision.SubjectId == other.Id &&
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved)))
                return "An alternative for this exact residue has already been approved.";
            if (approved.Value && (proposal.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.Disulfide) &&
                string.IsNullOrWhiteSpace(rationale))
                return "A site-specific reason is required for this chemical-state assumption.";
            if (approved.Value && !_inspection.HasRequiredEvidenceForApproval(proposal.Id,
                    _preparationProposals!.Evidence.Where(item => item.SubjectId == proposal.Id)
                        .Select(item => item.Id).ToImmutableArray(), out var reason)) return reason;
            var decision = new ResearcherDecision(NewId(), revision.Id, proposal.Id,
                ResearcherDecisionKind.ApprovePreparationChange,
                approved.Value ? ResearcherDecisionValue.Approved : ResearcherDecisionValue.Declined,
                DateTimeOffset.UtcNow, rationale);
            _decisions = _decisions.Add(decision);
            _allDecisions = _allDecisions.Add(decision);
            if (decision.ChosenValue == ResearcherDecisionValue.Declined &&
                proposal.Kind == PreparationChangeKind.HeavyAtom)
                Notice("warning", $"Required heavy-atom change {proposal.ProposedChange} at " +
                    $"{proposal.Residue.Chain}[{proposal.Residue.CopyId}]:{proposal.Residue.Residue}" +
                    $"{proposal.Residue.InsertionCode} was declined for study revision {revision.Id}; " +
                    "no assessed prepared protein can be established from this selection.", proposal.Id);
            TouchLocked();
        }
        Changed?.Invoke();
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
        }
        if (intended is null || inspection is null || proposals is null || policy is null || structuralPolicy is null ||
            !proposals.UnresolvedQuestions.IsDefaultOrEmpty ||
            !EveryRequiredChoiceSettled(proposals, decisions)) return;
        var prepared = await _proteinPreparation.PrepareAsync(revision, intended, inspection, policy, structuralPolicy,
            proposals, decisions,
            WorkDirectory("protein-preparation", revision.Id), cancellationToken);
        lock (_gate)
        {
            if (_study.Id != revision.Id) return;
            _protein = prepared.Value;
            _proteinDiagnostic = prepared.Value is null ? prepared.Diagnostic as ProteinPreparationDiagnostic : null;
            if (prepared.Value is null) Notice("warning", prepared.Reason, intended.Id);
            TouchLocked();
        }
        Changed?.Invoke();
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

    private string? ProposeMembrane(JsonElement data)
    {
        var upper = ParseArray<LipidFraction>(data, "upper");
        var lower = ParseArray<LipidFraction>(data, "lower");
        var purpose = Text(data, "scientificPurpose");
        if (!Coherent(upper) || !Coherent(lower) || purpose is null)
            return "Both leaflets need coherent finite species fractions and an explicit scientific purpose.";
        lock (_gate)
        {
            _membraneProposal = new MembraneModel(NewId(), new LeafletComposition(LeafletSide.Upper, upper),
                new LeafletComposition(LeafletSide.Lower, lower), _study.Conditions, purpose);
            _inspection.Clear();
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private async Task<string?> AdoptMembraneAsync(JsonElement data, CancellationToken cancellationToken)
    {
        MembraneModel proposal;
        StudyRevision revision;
        MembraneSupportPolicy? policy;
        IReadOnlyDictionary<string, MolecularRepresentation> lipids;
        lock (_gate)
        {
            if (_membraneProposal is null || _membraneProposal.Id != Text(data, "modelId"))
                return "Choose the currently proposed membrane model.";
            proposal = _membraneProposal;
            AdvanceStudyLocked(_study.IntendedProtein, proposal, null);
            revision = _study;
            _allDecisions = _allDecisions.Add(new ResearcherDecision(NewId(), revision.Id,
                proposal.Id, ResearcherDecisionKind.AdoptMembrane, ResearcherDecisionValue.Adopted, DateTimeOffset.UtcNow,
                proposal.ScientificPurpose));
            _protein = CarryProteinLocked(revision);
            _placementProposal = null;
            _opmReview = null;
            _placementMeasurement = null;
            _placementWitness = null;
            _placement = null;
            _membrane = null;
            policy = SelectMembranePolicyLocked(proposal);
            lipids = LipidsLocked();
            TouchLocked();
        }
        Changed?.Invoke();
        var assessed = await _membraneAssessment.AssessAsync(revision, proposal, lipids, policy,
            WorkDirectory("membrane-assessment", revision.Id), cancellationToken);
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _membrane = assessed.Value;
            if (assessed.Value is null) Notice("warning", assessed.Reason, proposal.Id);
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private async Task<string?> ProposePlacementAsync(JsonElement data, CancellationToken cancellationToken)
    {
        StudyRevision revision;
        AssessedPreparedProtein protein;
        MembraneModel membrane;
        string? ppmVersion;
        lock (_gate)
        {
            if (_protein is null || _study.Membrane is null || _membrane is null)
                return "A corresponding assessed protein and membrane model are required before placement.";
            revision = _study;
            protein = _protein;
            membrane = _study.Membrane;
            ppmVersion = _catalogue?.PpmVersion;
            if (string.IsNullOrWhiteSpace(ppmVersion) ||
                !VerifyHash(_ppmExecutablePath, _catalogue?.PpmExecutableSha256 ?? string.Empty))
                return "The local PPM executable does not match an identified, versioned installation.";
        }
        var topologyKind = ParseProteinTopology(Text(data, "topologyKind"));
        var physicalSide = ParsePlacementSide(Text(data, "physicalSide"));
        var ppmNterminalSide = ParsePpmNterminalSide(Text(data, "ppmNterminalSide"));
        if (topologyKind is null || physicalSide is null || ppmNterminalSide is null)
            return "Choose an established protein topology, physical bilayer side and PPM N-terminal assignment.";
        var opmTask = protein.Intended.Source.Kind == SourceRouteKind.Rcsb &&
            protein.Intended.Source.Accession is { } accession
            ? TryOptionalOpmReferenceAsync(accession,
                WorkDirectory("opm-reference", revision.Id), cancellationToken)
            : Task.FromResult<OpmReferenceRecord?>(null);
        var ppmTask = _placementAssessment.ProposeWithPpmAsync(revision, protein, membrane,
            topologyKind.Value, physicalSide.Value,
            Text(data, "biologicalSidedness"), ppmNterminalSide.Value,
            _ppmExecutablePath, ppmVersion ?? string.Empty,
            _catalogue?.PpmExecutableSha256 ?? string.Empty,
            WorkDirectory("placement", NewId()), cancellationToken);
        await Task.WhenAll(opmTask, ppmTask);
        var proposed = ppmTask.Result;
        var reference = opmTask.Result;
        var opm = reference is null ? null :
            _placementAssessment.ReviewOpmReference(revision, protein, membrane, reference).Value;
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _placementProposal = proposed.Value;
            _inspection.Clear();
            _opmReview = opm;
            _placementMeasurement = null;
            _placementWitness = null;
            _placement = null;
            if (proposed.Value is null) Notice("warning", proposed.Reason, protein.Id);
            if (protein.Intended.Source.Kind == SourceRouteKind.Rcsb && reference is null)
                Notice("information", "No identified OPM reference was available for this source; local PPM remains an independent placement route.", protein.Id);
            if (opm is { CorrespondsToSelectedConstruct: false })
                Notice("information", "An OPM reference was found but its exact construct applicability is not established.", protein.Id);
            TouchLocked();
        }
        Changed?.Invoke();
        if (proposed.Value is not null) await AssessPlacementAsync(revision, proposed.Value, cancellationToken);
        return null;
    }

    private async Task<OpmReferenceRecord?> TryOptionalOpmReferenceAsync(string accession,
        string directory, CancellationToken cancellationToken)
    {
        try { return await _sources.TryRetrieveOpmReferenceAsync(accession, directory, cancellationToken); }
        catch (Exception exception) when (exception is HttpRequestException or IOException or JsonException or
                                         InvalidDataException || exception is OperationCanceledException &&
                                         !cancellationToken.IsCancellationRequested)
        {
            // A reference is independent contextual evidence, not a prerequisite
            // for invoking local PPM or judging its explicit-membrane candidate.
            return null;
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
        var revised = await _placementAssessment.ReviseProposalAsync(revision, source,
            values[0]!.Value, values[1]!.Value, values[2]!.Value, values[3]!.Value,
            Text(data, "rationale") ?? string.Empty, WorkDirectory("placement", NewId()), cancellationToken);
        lock (_gate)
        {
            if (_study.Id != revision.Id) return null;
            _placementProposal = revised.Value;
            _inspection.Clear();
            _placementMeasurement = null;
            _placementWitness = null;
            _placement = null;
            if (revised.Value is null) Notice("warning", revised.Reason, source.Id);
            TouchLocked();
        }
        Changed?.Invoke();
        if (revised.Value is not null) await AssessPlacementAsync(revision, revised.Value, cancellationToken);
        return null;
    }

    private async Task AssessPlacementAsync(StudyRevision revision, PlacementProposal proposal, CancellationToken cancellationToken)
    {
        AssessedPreparedProtein? protein;
        AssessedMembraneModel? membrane;
        PlacementSupportPolicy? policy;
        PlacementStructuralWitness? witness;
        lock (_gate)
        {
            protein = _protein;
            membrane = _membrane;
            policy = membrane is null ? null : SelectPlacementPolicyLocked(proposal, membrane);
            witness = protein is null || membrane is null ? null :
                SelectPlacementWitnessLocked(protein, membrane, proposal);
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
            _placementWitness = witness;
            _placement = assessed;
            if (measurementIssue is not null) Notice("warning", measurementIssue, proposal.Id);
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
                _placement.Proposal.Id != _placementProposal.Id || _placement.Standing != AssessmentStanding.Supported ||
                _protein is null || _membrane is null || _placementMeasurement is null)
                return "Only the current supported proposal can be adopted.";
            if (!_inspection.HasRequiredEvidenceForApproval(_placementProposal.Id,
                    _placementMeasurement.Evidence.Select(item => item.Id).ToImmutableArray(), out var reason)) return reason;
            var revised = new StudyRevision(NewId(), _study.Number + 1, _study.IntendedProtein,
                _study.Membrane, _placementProposal.Id, _study.Conditions);
            var carriedProtein = CarryProteinLocked(revised);
            var carriedMembrane = CarryMembraneLocked(revised);
            if (carriedProtein is null || carriedMembrane is null)
                return "An affected assessed input cannot be carried into this changed study premise.";
            var policy = SelectPlacementPolicyLocked(_placementProposal, carriedMembrane);
            var reassessed = _placementAssessment.Assess(revised, carriedProtein, carriedMembrane,
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
                "The inspected proposal met the corresponding placement policy."));
            _protein = carriedProtein;
            _membrane = carriedMembrane;
            _placement = reassessed;
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private string? StartPreparation()
    {
        lock (_gate)
        {
            if (_protein is null || _membrane is null || _placement is null ||
                _placement.Standing != AssessmentStanding.Supported ||
                _study.AdoptedPlacementProposalId != _placement.Proposal.Id ||
                _placement.StudyRevisionId != _study.Id || _protein.StudyRevisionId != _study.Id ||
                _membrane.StudyRevisionId != _study.Id)
                return "Current, corresponding supported inputs are required before preparation.";
            if (_attemptTask is { IsCompleted: false }) return "The current attempt is still running.";
            var policy = SelectPreparationPolicyLocked(_placement.Proposal, _membrane);
            if (policy is null) return "No qualified, applicable preparation policy and exact assets are available.";
            if (string.IsNullOrWhiteSpace(_catalogue?.PackmolVersion) ||
                !VerifyHash(_packmolExecutablePath, _catalogue.PackmolExecutableSha256))
                return "The local Packmol executable does not match an identified, versioned installation.";
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
            _attemptStop = attemptStop;
            _currentAttempt = null;
            _execution = new StageExecutionState(string.Empty, null, null, StageExecutionStanding.Pending,
                "The corresponding attempt is being admitted.", null, DateTimeOffset.UtcNow);
            _attemptTask = Task.Run(() => RunPreparationAsync(revision, protein, membrane, placement, policy, decisions,
                attemptStop.Token));
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
    }

    private string? StopAttempt(JsonElement data)
    {
        lock (_gate)
        {
            var id = Text(data, "attemptId");
            if (_attemptTask is not { IsCompleted: false } || _attemptStop is null ||
                _currentAttempt?.Id != id)
                return "No matching unfinished preparation attempt can be stopped.";
            _attemptStop.Cancel();
            Notice("information", "Stop requested; the actual stage standing will be reported when the worker stops.", id);
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
            if (!_constructedByAttempt.TryGetValue(minimized.Attempt.Id, out var constructed) ||
                !_revisions.ContainsKey(minimized.Attempt.StudyRevisionId) ||
                !_policyByAttempt.TryGetValue(minimized.Attempt.Id, out var policy))
                return "The selected stage has no surviving, corresponding attempt and policy context.";
            if (policy.OptionalEquilibration is null || !QualifiedEquilibrationLocked(policy, minimized.Attempt.Id))
                return "No validated, applicable optional equilibration procedure is available for this stage.";
            var stop = new CancellationTokenSource();
            _attemptStop = stop;
            _currentAttempt = minimized.Attempt;
            _attemptTask = Task.Run(() => RunEquilibrationAsync(minimized, constructed, policy, stop.Token));
            TouchLocked();
        }
        Changed?.Invoke();
        return null;
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
        ImmutableArray<ResearcherDecision> decisions;
        lock (_gate)
        {
            if (stageId is null || !_stages.TryGetValue(stageId, out var selected) ||
                !_stageAssessments.TryGetValue(stageId, out var currentAssessment) ||
                !_revisions.TryGetValue(selected.Attempt.StudyRevisionId, out var originating) ||
                !_policyByAttempt.TryGetValue(selected.Attempt.Id, out var boundPolicy) ||
                !_proteinByAttempt.TryGetValue(selected.Attempt.Id, out var originatingProtein) ||
                !_decisionsByAttempt.TryGetValue(selected.Attempt.Id, out decisions))
                return "The selected completed stage and its identified assessment are unavailable.";
            stage = selected;
            assessment = currentAssessment;
            revision = originating;
            policy = boundPolicy;
            protein = originatingProtein;
            if (_bundles.ContainsKey(stageId)) return null;
        }
        var delivered = await _export.ExportAsync(stage, assessment, revision, protein, policy, decisions,
            assessment.Findings, WorkDirectory("exports", stage.Id + "-" + assessment.Id), cancellationToken);
        lock (_gate)
        {
            if (!_stageAssessments.TryGetValue(stage.Id, out var current) || current.Id != assessment.Id ||
                !current.CurrentlyApplicable)
                return "The stage assessment changed during export; review the current assessment before delivery.";
            if (delivered.Value is not null) _bundles[stage.Id] = delivered.Value;
            else Notice("warning", delivered.Reason, stage.Id);
            TouchLocked();
        }
        Changed?.Invoke();
        return delivered.Value is null ? delivered.Reason : null;
    }

    private async Task RunPreparationAsync(StudyRevision revision, AssessedPreparedProtein protein,
        AssessedMembraneModel membrane, AssessedProteinMembranePlacement placement,
        ApplicablePreparationPolicy policy, ImmutableArray<ResearcherDecision> decisions,
        CancellationToken cancellationToken)
    {
        try
        {
            var attempt = new PreparationAttempt(NewId(), revision.Id, protein.Id, membrane.Id,
                placement.Id, policy.Id, DateTimeOffset.UtcNow);
            var started = await _explicitPreparation.StartAsync(attempt, revision, protein, membrane, placement, policy,
                _packmolExecutablePath, _catalogue!.PackmolVersion, _catalogue.PackmolExecutableSha256,
                WorkDirectory("attempts", attempt.Id), accepted =>
                {
                    lock (_gate)
                    {
                        _currentAttempt = accepted;
                        if (_study.Id == revision.Id) _workspace.RetainAttempt(accepted);
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
                _currentAttempt = started.Attempt;
                _policyByAttempt[started.Attempt.Id] = policy;
                _proteinByAttempt[started.Attempt.Id] = protein;
                _membraneByAttempt[started.Attempt.Id] = membrane;
                _placementByAttempt[started.Attempt.Id] = placement;
                _decisionsByAttempt[started.Attempt.Id] = decisions;
                if (_study.Id == revision.Id)
                {
                    _workspace.RetainExecution(started.State);
                }
                _execution = started.State;
                TouchLocked();
            }
            Changed?.Invoke();
            if (started.Constructed is null) return;
            lock (_gate) _constructedByAttempt[started.Attempt.Id] = started.Constructed;
            var minimized = await _explicitPreparation.MinimizeAsync(started.Constructed, policy,
                WorkDirectory("attempts", started.Attempt.Id), Progress(), cancellationToken);
            SetExecution(minimized.State);
            if (minimized.CompletedStage is not null)
                RetainCompleted(minimized.CompletedStage, started.Constructed, protein, membrane, placement, policy);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
        {
            SetExecution(new StageExecutionState(_currentAttempt?.Id ?? string.Empty, null, null,
                StageExecutionStanding.Stopped, "The unfinished preparation work was stopped.", null, DateTimeOffset.UtcNow));
        }
        catch (Exception exception)
        {
            SetExecution(new StageExecutionState(_currentAttempt?.Id ?? string.Empty, null, null,
                StageExecutionStanding.Unobserved, "The local worker outcome could not be established: " + exception.Message,
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
            SetExecution(result.State);
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
                    RetainCompleted(result.CompletedStage, constructed, protein, membrane, placement, policy);
            }
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

    private sealed class DirectProgress(Action<StageExecutionState> receive) : IProgress<StageExecutionState>
    {
        public void Report(StageExecutionState value) => receive(value);
    }

    private void SetExecution(StageExecutionState state)
    {
        lock (_gate)
        {
            _execution = state;
            if (_currentAttempt?.Id == state.AttemptId && _study.Id == _currentAttempt.StudyRevisionId)
                _workspace.RetainExecution(state);
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private void RetainCompleted(CompletedStage stage, ConstructedExplicitSystem constructed,
        AssessedPreparedProtein protein, AssessedMembraneModel membrane,
        AssessedProteinMembranePlacement placement, ApplicablePreparationPolicy policy)
    {
        lock (_gate)
        {
            var judgedPlacement = placement;
            var newPlacementFindings = stage.Findings.Where(finding => finding.Material &&
                finding.SubjectId == placement.Proposal.Id &&
                (finding.Disposition is FindingDisposition.Challenges or FindingDisposition.Disqualifies))
                .ToImmutableArray();
            if (!newPlacementFindings.IsEmpty && _study.Id == stage.Attempt.StudyRevisionId &&
                _placement?.Id == placement.Id && _protein?.Id == protein.Id &&
                _membrane?.Id == membrane.Id && _placementMeasurement?.ProposalId == placement.Proposal.Id)
            {
                judgedPlacement = _placementAssessment.Assess(_study, _protein, _membrane,
                    placement.Proposal, SelectPlacementPolicyLocked(placement.Proposal, _membrane),
                    _placementMeasurement, _placementWitness, _placementMeasurement.Evidence,
                    placement.Findings.AddRange(newPlacementFindings));
                _placement = judgedPlacement;
                Notice("warning", "A completed-stage observation challenges the previously adopted placement; its support has been reassessed.",
                    placement.Proposal.Id);
            }
            // Stage.Findings already belongs to this stage. The additional argument is
            // reserved for genuinely later observations, not a duplicate of its own.
            var assessment = _preparationAssessment.Assess(stage, constructed, protein, membrane,
                judgedPlacement, policy, ImmutableArray<ScientificFinding>.Empty);
            _stages[stage.Id] = stage;
            _stageAssessments[stage.Id] = assessment;
            foreach (var finding in stage.Findings) _workspace.RetainFinding(finding);
            var consequentialLaterFindings = stage.Findings.Where(item => item.Material &&
                (item.Disposition is FindingDisposition.Challenges or FindingDisposition.Disqualifies))
                .ToImmutableArray();
            foreach (var earlier in _stages.Values.Where(item => !consequentialLaterFindings.IsEmpty &&
                         item.Id != stage.Id && item.Attempt.Id == stage.Attempt.Id).ToArray())
            {
                var accumulated = _laterFindingsByStage.GetValueOrDefault(earlier.Id);
                if (accumulated.IsDefault) accumulated = ImmutableArray<ScientificFinding>.Empty;
                var newFindings = consequentialLaterFindings.Where(item =>
                    accumulated.All(previous => previous.Id != item.Id));
                accumulated = accumulated.AddRange(newFindings);
                _laterFindingsByStage[earlier.Id] = accumulated;
                var renewed = _preparationAssessment.Assess(earlier, constructed, protein, membrane,
                    judgedPlacement, policy, accumulated);
                _stageAssessments[earlier.Id] = renewed;
                _bundles.Remove(earlier.Id);
                if (_currentAttempt?.Id == earlier.Attempt.Id && _study.Id == earlier.Attempt.StudyRevisionId)
                    _workspace.RetainAssessment(renewed);
            }
            if (_currentAttempt?.Id == stage.Attempt.Id && _study.Id == stage.Attempt.StudyRevisionId)
            {
                _workspace.RetainCompletedStage(stage);
                _workspace.RetainAssessment(assessment);
            }
            TouchLocked();
        }
        Changed?.Invoke();
    }

    private (InspectionSubject?, StudyRevision?) MakeInspectionSubjectLocked(string id)
    {
        if (_selectedSource is not null && _sourceInspection is not null && id == _selectedSource.Id)
            return (GenericSubject(id, _study.Id, _selectedSource.CoordinatePath, "structuralSource",
                ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty, null), _study);
        if (_study.IntendedProtein is not null && id == _study.IntendedProtein.Id)
            return (WithGeometry(GenericSubject(id, _study.Id, _study.IntendedProtein.Source.CoordinatePath,
                "intendedProtein", ImmutableArray<ScientificEvidence>.Empty,
                ImmutableArray<ScientificFinding>.Empty, null), _preparationProposals?.SourceGeometry), _study);
        if (_preparationProposals?.Changes.Any(item => item.Id == id) == true)
        {
            var preview = _preparationProposals.Preview;
            if (preview is null || !IsWorkspaceFile(preview.Path)) return (null, null);
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
        if (_protein is not null && id == _protein.Id)
            return (WithGeometry(GenericSubject(id, _study.Id, _protein.Molecule.CoordinatePath, "preparedProtein",
                _protein.Evidence.AddRange(_opmReview?.Evidence ?? ImmutableArray<ScientificEvidence>.Empty),
                _protein.Findings, null), _protein.Geometry), _study);
        if (_membraneProposal is not null && id == _membraneProposal.Id)
            return (GenericSubject(id, _study.Id, null, "membraneModel", ImmutableArray<ScientificEvidence>.Empty,
                ImmutableArray<ScientificFinding>.Empty, null), _study);
        if (_placementProposal is not null && id == _placementProposal.Id)
        {
            if (_placementMeasurement is null || !IsWorkspaceFile(_placementProposal.OrientedProtein.CoordinatePath))
                return (null, null);
            var url = StructureUrlLocked(_placementProposal.OrientedProtein.CoordinatePath);
            var subject = _placementAssessment.ProposalInspectionSubject(_study, _placementProposal,
                _placementMeasurement, url);
            return (subject.Value, _study);
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

    // A measured distance is shown with its exact source addresses. Prepared
    // atom serials are not inferred from those addresses, so these metrics do
    // not advertise spatial focus until a verified view mapping exists.
    private static InspectionSubject WithGeometry(InspectionSubject subject, ProteinGeometryObservations? geometry)
    {
        if (geometry is null || geometry.Kinds.IsDefault) return subject;
        var evidence = subject.Evidence.ToBuilder();
        var metrics = subject.Metrics.ToBuilder();
        foreach (var kind in geometry.Kinds)
        {
            var evidenceId = NewId();
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
        foreach (var distance in geometry.LocatedDistances.IsDefault
            ? ImmutableArray<GeometryDistanceObservation>.Empty : geometry.LocatedDistances)
        {
            var first = distance.First.Residue;
            var second = distance.Second.Residue;
            var label = $"{first.Chain}/{first.CopyId}:{first.Residue}{first.InsertionCode}/{distance.First.AtomName} ↔ " +
                $"{second.Chain}/{second.CopyId}:{second.Residue}{second.InsertionCode}/{distance.Second.AtomName}";
            var evidenceId = NewId();
            evidence.Add(new ScientificEvidence(evidenceId, subject.Id, "local scientific worker", distance.Kind,
                $"Located atom-pair separation for {label}", $"Exact inspected {subject.RepresentationKind} {subject.Id}",
                "Addressed numerical observation; no unverified spatial focus is inferred.", EvidenceBearing.Context));
            metrics.Add(new InspectionMetric(label, distance.DistanceAngstrom.ToString("G6", CultureInfo.InvariantCulture),
                "Å", subject.Id, evidenceId));
        }
        return subject with { Evidence = evidence.ToImmutable(), Metrics = metrics.ToImmutable() };
    }

    private static InspectionSubject WithStageMeasurements(InspectionSubject subject, CompletedStage stage,
        ApplicablePreparationPolicy? policy)
    {
        var evidence = subject.Evidence.ToBuilder();
        var metrics = subject.Metrics.ToBuilder();
        foreach (var value in stage.Observation.Measurements.IsDefault
            ? ImmutableArray<MeasuredValue>.Empty : stage.Observation.Measurements)
        {
            if (!double.IsFinite(value.Value) || string.IsNullOrWhiteSpace(value.Name) ||
                string.IsNullOrWhiteSpace(value.Unit) || string.IsNullOrWhiteSpace(value.Scope)) continue;
            var evidenceId = NewId();
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
            foreach (var pair in local.RolePairMeasurements.OrderByDescending(item =>
                         item.FirstMoleculeRole == MoleculeRoleKind.Protein && item.SecondMoleculeRole == MoleculeRoleKind.Lipid))
            {
                if (pair.PairsWithinSearchRadius < 0 || !Enum.IsDefined(pair.FirstMoleculeRole) ||
                    !Enum.IsDefined(pair.SecondMoleculeRole)) continue;
                var label = MoleculeRoleTokens.PairKey(pair.FirstMoleculeRole, pair.SecondMoleculeRole);
                var evidenceId = NewId();
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
        return subject with { Evidence = evidence.ToImmutable(), Metrics = metrics.ToImmutable() };
    }

    private WorkspaceState SnapshotLocked()
    {
        var stages = _stages.Values.OrderBy(item => item.CompletedAt).Select(stage => new StageAccount(
            stage.Id, stage.Attempt.Id, stage.Attempt.StudyRevisionId, stage.Kind, "completed",
            _stageAssessments.GetValueOrDefault(stage.Id),
            $"{stage.Kind} completed for study revision {stage.Attempt.StudyRevisionId}.")).ToImmutableArray();
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
        var membrane = _membraneProposal is null ? null : new MembraneAccount(
            _membraneProposal.Id, _membrane?.Intended.Id == _membraneProposal.Id ? "assessed" :
            _study.Membrane?.Id == _membraneProposal.Id ? "notEstablished" : "proposed",
            _membraneProposal.ScientificPurpose,
            _membraneProposal.Upper.Fractions, _membraneProposal.Lower.Fractions,
            _membrane?.Limitations ?? ImmutableArray<string>.Empty);
        var placement = _placementProposal is null ? null : new PlacementAccount(
            _placementProposal.Id, _placement?.Standing.ToString() ?? "proposed",
            _placementProposal.TopologyKind,
            _placementProposal.MidplaneAngstrom, _placementProposal.TiltDegrees,
            _placementProposal.BiologicalSidedness,
            _placement?.Reason ?? "The positioned candidate still needs support assessment.",
            _placementProposal.Evidence, _placementMeasurement?.Prediction);
        var studyAccount = new StudyAccount(_study.Id, _study.Number,
            _study.IntendedProtein is null ? "Choose an exact protein structure." : "One identified protein–membrane study.",
            _selectedSource?.Id, _study.IntendedProtein?.ModelIndex,
            _study.IntendedProtein?.BiologicalAssemblyId,
            _study.IntendedProtein?.Chains.Select(item => item.CopyId).ToImmutableArray() ?? ImmutableArray<string>.Empty,
            _study.IntendedProtein?.Partners ?? ImmutableArray<PartnerSelection>.Empty,
            _study.IntendedProtein?.AlternateLocations ?? ImmutableArray<AlternateLocationChoice>.Empty,
            _study.Conditions, _selectedSource?.Kind, _selectedSource?.UploadProvenance,
            _selectedSource?.UploadProvenanceNote);
        return new WorkspaceState(_revision, studyAccount,
            _candidates.Select(candidate => new SourceCandidateAccount(candidate.Id, candidate.Label,
                candidate.Kind, candidate.Provenance, candidate.Limitations)).ToImmutableArray(),
            _sourceInspection?.Models ?? ImmutableArray<SourceModelObservation>.Empty,
            LipidsLocked().Values.OrderBy(item => item.SpeciesId).Select(item => new LipidCatalogueAccount(
                item.SpeciesId, item.SpeciesId, item.ChemistryId, item.Limitations)).ToImmutableArray(),
            protein, membrane, placement,
            _currentAttempt is null && _execution is null ? null : new AttemptAccount(
                _currentAttempt?.Id ?? _execution?.AttemptId ?? string.Empty,
                _execution?.Standing.ToString() ?? "pending",
                _execution?.Kind, _execution?.Progress,
                _execution?.Message ?? "The local attempt has not started."),
            stages, _inspection.Current, ActionsLocked(stages), _notices.ToImmutableArray(),
            PredictionAccount(_sourceInspection?.Prediction));
    }

    private static PredictionEvidenceAccount? PredictionAccount(PredictionEvidenceObservations? observations) =>
        observations is null ? null : new PredictionEvidenceAccount(observations.RecordId,
            observations.LocalConfidence, observations.PaeStanding, observations.PaeReason,
            observations.PaeAxisResidueCount, observations.Limitations);

    private ImmutableArray<AvailableAction> ActionsLocked(ImmutableArray<StageAccount> stages)
    {
        var actions = ImmutableArray.CreateBuilder<AvailableAction>();
        void Add(ActorActionKind kind, string? subject, bool enabled, string reason) =>
            actions.Add(new AvailableAction(kind, subject, enabled, enabled ? null : reason));
        Add(ActorActionKind.SearchSource, null, true, "");
        Add(ActorActionKind.SelectSource, null, true, "");
        Add(ActorActionKind.SelectProteinModel, null, _sourceInspection is not null, "Choose and inspect a structural source first.");
        foreach (var change in _preparationProposals?.Changes ?? ImmutableArray<PreparationChangeProposal>.Empty)
        {
            var alternativeChosen = change.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState &&
                _preparationProposals!.Changes.Any(other => other.Id != change.Id &&
                    other.Kind == change.Kind && other.Residue == change.Residue &&
                    _decisions.Any(decision => decision.SubjectId == other.Id &&
                        decision.Kind == ResearcherDecisionKind.ApprovePreparationChange &&
                        decision.ChosenValue == ResearcherDecisionValue.Approved));
            var undecided = _protein is null && !alternativeChosen && !_decisions.Any(item => item.SubjectId == change.Id);
            var visible = _inspection.Current?.SubjectId == change.Id;
            Add(ActorActionKind.ApprovePreparationChange, change.Id, undecided && visible,
                undecided ? "Inspect this exact proposal and its evidence before approval." : "This proposal has already been decided.");
            Add(ActorActionKind.DeclinePreparationChange, change.Id, undecided, "This proposal has already been decided.");
        }
        Add(ActorActionKind.ProposeMembrane, null, LipidsLocked().Count > 0, "No qualified lipid representation is available.");
        Add(ActorActionKind.AdoptMembrane, null, _membraneProposal is not null, "Propose a complete membrane model first.");
        Add(ActorActionKind.ProposePlacement, null, _protein is not null && _membrane is not null &&
            !string.IsNullOrWhiteSpace(_catalogue?.PpmVersion) &&
            VerifyHash(_ppmExecutablePath, _catalogue?.PpmExecutableSha256 ?? string.Empty),
            "A corresponding prepared protein, assessed membrane and identified PPM executable are required.");
        Add(ActorActionKind.RevisePlacement, null, _placementProposal is not null, "A positioned candidate is needed first.");
        Add(ActorActionKind.AdoptPlacement, null, _placement?.Standing == AssessmentStanding.Supported &&
            _inspection.Current?.SubjectId == _placementProposal?.Id,
            "Review a current, supported proposal before adopting it.");
        Add(ActorActionKind.StartPreparation, null, _placement?.Standing == AssessmentStanding.Supported &&
            _study.AdoptedPlacementProposalId == _placement?.Proposal.Id && _attemptTask is not { IsCompleted: false } &&
            _membrane is not null && SelectPreparationPolicyLocked(_placement!.Proposal, _membrane) is not null &&
            !string.IsNullOrWhiteSpace(_catalogue?.PackmolVersion) &&
            VerifyHash(_packmolExecutablePath, _catalogue?.PackmolExecutableSha256 ?? string.Empty),
            "A corresponding adopted placement, qualified preparation policy and identified Packmol executable are required.");
        Add(ActorActionKind.StopAttempt, null, _attemptTask is { IsCompleted: false } && _currentAttempt is not null,
            "No identified unfinished attempt can be stopped.");
        Add(ActorActionKind.SelectInspectionSubject, null, _sourceInspection is not null || _membraneProposal is not null ||
            _placementProposal is not null || _proteinDiagnostic is not null || stages.Length > 0,
            "There is no identified subject to inspect.");
        Add(ActorActionKind.SetInspectionFocus, null, _inspection.Current is not null, "Select an inspection subject first.");
        foreach (var stage in stages)
        {
            var completed = _stages[stage.StageId];
            Add(ActorActionKind.ExportStage, stage.StageId, _stageAssessments.TryGetValue(stage.StageId, out var assessment) &&
                assessment.CurrentlyApplicable, "This stage needs a current corresponding assessment before export.");
            if (completed.Kind == StageKind.Minimization)
                Add(ActorActionKind.RequestEquilibration, stage.StageId,
                    _attemptTask is not { IsCompleted: false } &&
                    _policyByAttempt.TryGetValue(completed.Attempt.Id, out var policy) &&
                    QualifiedEquilibrationLocked(policy, completed.Attempt.Id),
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
        if (_membrane is null || revision.Membrane?.Id != _membrane.Intended.Id ||
            _membrane.Intended.Conditions != revision.Conditions ||
            SelectMembranePolicyLocked(_membrane.Intended) is null) return null;
        return _membrane with { StudyRevisionId = revision.Id };
    }

    private ProteinChemicalStatePolicy? SelectChemicalPolicyLocked(SourceModelObservation? model,
        IntendedProteinModel? intended)
    {
        if (model is null || intended is null || _catalogue is null ||
            intended.Partners.Any(partner => partner.Retain)) return null;
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
            membrane.Upper.Fractions.Concat(membrane.Lower.Fractions).All(item =>
                policy.CoveredSpeciesIds.Contains(item.SpeciesId))).Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private PlacementSupportPolicy? SelectPlacementPolicyLocked(PlacementProposal proposal, AssessedMembraneModel membrane)
    {
        if (_catalogue is null) return null;
        var applicable = _catalogue.PlacementPolicies.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty && !policy.GeometryCriteria.IsDefaultOrEmpty &&
            !policy.CoveredTopologyKinds.IsDefaultOrEmpty &&
            policy.CoveredTopologyKinds.All(kind => Enum.IsDefined(kind)) &&
            !policy.CoveredSpeciesIds.IsDefaultOrEmpty &&
            double.IsFinite(policy.InterfaceBandAngstrom) && policy.InterfaceBandAngstrom > 0 &&
            policy.CoveredTopologyKinds.Contains(proposal.TopologyKind) &&
            membrane.SpeciesRepresentations.All(item => policy.CoveredSpeciesIds.Contains(item.SpeciesId)))
            .Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
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
            witness.PreparedCoordinateSha256 == protein.Molecule.CoordinateSha256 &&
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
        left.Fractions.All(item => right.Fractions.Any(other => other.SpeciesId == item.SpeciesId &&
            Math.Abs(other.Fraction - item.Fraction) <= 1e-9));

    private ApplicablePreparationPolicy? SelectPreparationPolicyLocked(PlacementProposal proposal, AssessedMembraneModel membrane)
    {
        if (_catalogue is null) return null;
        var applicable = _catalogue.PreparationPolicies.Where(policy =>
            !string.IsNullOrWhiteSpace(policy.Id) && !string.IsNullOrWhiteSpace(policy.Version) &&
            !policy.EvidenceReferences.IsDefaultOrEmpty &&
            policy.ApplicableMolecularClass == "canonical-amino-acid-assembly" &&
            policy.Construction is { } construction &&
            !string.IsNullOrWhiteSpace(construction.Id) &&
            !string.IsNullOrWhiteSpace(construction.Version) &&
            !construction.EvidenceReferences.IsDefaultOrEmpty &&
            !construction.CoveredTopologyKinds.IsDefaultOrEmpty &&
            construction.CoveredTopologyKinds.All(kind => Enum.IsDefined(kind)) &&
            !construction.CoveredSpeciesIds.IsDefaultOrEmpty &&
            policy.SystemSettings is { NonbondedMethod: "PME", Constraints: "HBonds" } &&
            double.IsFinite(policy.SystemSettings.NonbondedCutoffNanometers) &&
            policy.SystemSettings.NonbondedCutoffNanometers > 0 &&
            double.IsFinite(policy.SystemSettings.EwaldErrorTolerance) &&
            policy.SystemSettings.EwaldErrorTolerance is > 0 and < 1 &&
            (policy.SystemSettings.SwitchDistanceNanometers is null ||
             policy.SystemSettings.SwitchDistanceNanometers is double switchDistance &&
             double.IsFinite(switchDistance) && switchDistance > 0 &&
             switchDistance < policy.SystemSettings.NonbondedCutoffNanometers) &&
            (policy.SystemSettings.HydrogenMassDaltons is null ||
             policy.SystemSettings.HydrogenMassDaltons is double hydrogenMass &&
             double.IsFinite(hydrogenMass) && hydrogenMass > 0) &&
            construction.CoveredTopologyKinds.Contains(proposal.TopologyKind) &&
            membrane.SpeciesRepresentations.All(item => construction.CoveredSpeciesIds.Contains(item.SpeciesId)) &&
            ValidLocalStatePolicy(policy) &&
            ValidStageProteinGeometryPolicy(policy) &&
            !policy.ForceFieldFiles.IsDefaultOrEmpty && policy.ForceFieldFiles.All(VerifiedAsset) &&
            new[] { policy.Water, policy.Sodium, policy.Chloride }.All(VerifiedRepresentation))
            .Take(2).ToArray();
        return applicable.Length == 1 ? applicable[0] : null;
    }

    private static bool ValidLocalStatePolicy(ApplicablePreparationPolicy policy)
    {
        var observation = policy.LocalStateObservation;
        var criteria = policy.ConstructionCriteria;
        var contacts = policy.ContactCriteria;
        if (observation is null || observation.ContactRolePairs.IsDefaultOrEmpty ||
            observation.RequiredMetricNames.IsDefaultOrEmpty || criteria.IsDefaultOrEmpty ||
            contacts.IsDefaultOrEmpty ||
            observation.AtomRadiusByElementAngstrom is not { Count: > 0 } ||
            observation.AtomRadiusByElementAngstrom.Any(item =>
                string.IsNullOrWhiteSpace(item.Key) || !double.IsFinite(item.Value) || item.Value <= 0) ||
            !double.IsFinite(observation.ContactSearchRadiusAngstrom) ||
            observation.ContactSearchRadiusAngstrom <= 0 ||
            observation.MaximumReportedPairs <= 0 || !observation.UsePeriodicBoundary ||
            observation.ContactRolePairs.Any(pair => !Enum.IsDefined(pair.FirstMoleculeRole) ||
                !Enum.IsDefined(pair.SecondMoleculeRole)) ||
            observation.ContactRolePairs.Distinct().Count() != observation.ContactRolePairs.Length ||
            observation.RequiredMetricNames.Any(string.IsNullOrWhiteSpace) ||
            observation.RequiredMetricNames.Distinct(StringComparer.Ordinal).Count() !=
                observation.RequiredMetricNames.Length ||
            observation.RequiredMetricNames.Any(name => name is not (
                "minimumIntermolecularDistanceAngstrom" or "upperLipidHeadMeanZAngstrom" or
                "lowerLipidHeadMeanZAngstrom" or "leafletHeadSeparationAngstrom" or
                "proteinBilayerMidplaneOffsetAngstrom")) ||
            criteria.Select(criterion => criterion.MeasurementName).Distinct(StringComparer.Ordinal).Count() !=
                criteria.Length ||
            criteria.Any(criterion => criterion.MeasurementName is
                "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom") ||
            (!policy.AssessmentCriteria.IsDefault && policy.AssessmentCriteria.Any(criterion =>
                criterion.MeasurementName is
                    "upperLipidHeadMeanZAngstrom" or "lowerLipidHeadMeanZAngstrom")) ||
            !criteria.Any(criterion => criterion.MeasurementName == "minimumIntermolecularDistanceAngstrom" &&
                criterion.Minimum is double distance && double.IsFinite(distance) && distance > 0) ||
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
            !new StageKind?[] { null, StageKind.Minimization, StageKind.Equilibration }
                .Where(kind => kind is null || kind != StageKind.Equilibration ||
                    policy.OptionalEquilibration is not null)
                .All(kind => contacts.Any(item => item.StageKind == kind &&
                    item.FirstMoleculeRole == MoleculeRoleKind.Protein && item.SecondMoleculeRole == MoleculeRoleKind.Lipid &&
                    item.MinimumPairsWithinSearchRadius > 0 &&
                    item.MaximumNearestDistanceAngstrom is double maximum &&
                    double.IsFinite(maximum) && maximum > 0)) ||
            criteria.Any(criterion => string.IsNullOrWhiteSpace(criterion.MeasurementName) ||
                !observation.RequiredMetricNames.Contains(criterion.MeasurementName) ||
                string.IsNullOrWhiteSpace(criterion.Unit) || string.IsNullOrWhiteSpace(criterion.Scope) ||
                (criterion.Minimum is double minimum && !double.IsFinite(minimum)) ||
                (criterion.Maximum is double maximum && !double.IsFinite(maximum)) ||
                (criterion.Minimum is double low && criterion.Maximum is double high && low > high))) return false;

        static bool Bounded(string unit, string scope, double? minimum, double? maximum,
            string requiredScope, bool requirePositiveMinimum) =>
            unit == "angstrom" && scope == requiredScope &&
            minimum is double lower && maximum is double upper &&
            double.IsFinite(lower) && double.IsFinite(upper) && lower <= upper &&
            (!requirePositiveMinimum || lower > 0);
        bool ConstructionBounded(string name, string scope, bool positive)
        {
            var matching = criteria.Where(item => item.MeasurementName == name).Take(2).ToArray();
            return matching.Length == 1 && Bounded(matching[0].Unit, matching[0].Scope,
                matching[0].Minimum, matching[0].Maximum, scope, positive);
        }
        bool StageBounded(StageKind kind, string name, string scope, bool positive)
        {
            if (policy.AssessmentCriteria.IsDefaultOrEmpty) return false;
            var matching = policy.AssessmentCriteria.Where(item => item.StageKind == kind &&
                item.MeasurementName == name).Take(2).ToArray();
            return matching.Length == 1 && Bounded(matching[0].Unit, matching[0].Scope,
                matching[0].Minimum, matching[0].Maximum, scope, positive);
        }
        if (!observation.RequiredMetricNames.Contains("leafletHeadSeparationAngstrom") ||
            !observation.RequiredMetricNames.Contains("proteinBilayerMidplaneOffsetAngstrom") ||
            !ConstructionBounded("leafletHeadSeparationAngstrom", "bilayer", true) ||
            !ConstructionBounded("proteinBilayerMidplaneOffsetAngstrom", "proteinVsBilayer", false))
            return false;
        return new[] { StageKind.Minimization, StageKind.Equilibration }
            .Where(kind => kind != StageKind.Equilibration || policy.OptionalEquilibration is not null)
            .All(kind => StageBounded(kind, "leafletHeadSeparationAngstrom", "bilayer", true) &&
                StageBounded(kind, "proteinBilayerMidplaneOffsetAngstrom", "proteinVsBilayer", false));
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

    private bool QualifiedEquilibrationLocked(ApplicablePreparationPolicy policy, string attemptId)
    {
        if (policy.OptionalEquilibration is not { } protocol || _catalogue is null ||
            string.IsNullOrWhiteSpace(protocol.Id) || protocol.Stages.IsDefaultOrEmpty ||
            protocol.ExtensionWindow is null || protocol.MaximumExtensions < 0 ||
            protocol.MaximumSampleCount <= 0 || protocol.RequiredObservations.IsDefaultOrEmpty ||
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
                CoordinateTemplatePath = Resolve(item.CoordinateTemplatePath)
            };
            ProteinChemicalStatePolicy MapChemical(ProteinChemicalStatePolicy item) => item with
            { ForceFieldFiles = item.ForceFieldFiles.Select(MapAsset).ToImmutableArray() };
            ApplicablePreparationPolicy MapPreparation(ApplicablePreparationPolicy item) => item with
            {
                ForceFieldFiles = item.ForceFieldFiles.Select(MapAsset).ToImmutableArray(),
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
        string PackmolVersion,
        string PackmolExecutableSha256,
        int MaximumSourceAtoms);

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
        var token = NewId();
        _structures[token] = path;
        var format = Path.GetExtension(path).ToLowerInvariant() is ".cif" or ".mmcif" ? "mmcif" : "pdb";
        return "/api/structures/" + token + "?format=" + format;
    }

    private static string Hash(string path)
    {
        using var input = File.OpenRead(path);
        return Convert.ToHexString(SHA256.HashData(input)).ToLowerInvariant();
    }

    private static string NewId() => Guid.NewGuid().ToString("N");
    private void TouchLocked() => _revision++;
    private void Notice(string severity, string message, string? subjectId)
    {
        _notices.Add(new WorkspaceNotice(NewId(), severity, message, subjectId));
        if (_notices.Count > 16) _notices.RemoveAt(0);
    }
    private void AddNotice(string severity, string message, string? subjectId)
    {
        lock (_gate) { Notice(severity, message, subjectId); TouchLocked(); }
        Changed?.Invoke();
    }

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
}
