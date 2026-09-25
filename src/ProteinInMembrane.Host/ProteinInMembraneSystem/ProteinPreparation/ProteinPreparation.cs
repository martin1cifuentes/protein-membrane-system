using System.Collections.Immutable;
using System.Text.Json;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinPreparation;

/// <summary>Establishes the exact prepared protein, never placement or whole-system support.</summary>
public sealed class ProteinPreparation
{
    private readonly IProteinPreparationWork _worker;

    public ProteinPreparation(IProteinPreparationWork worker) => _worker = worker;

    public async Task<BoundaryOutcome<SourceInspectionReport>> InspectSourceAsync(
        StructuralSource source,
        string workingDirectory,
        int maxAtoms,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(source);
        if (string.IsNullOrWhiteSpace(source.CoordinatePath) || string.IsNullOrWhiteSpace(source.Sha256) || maxAtoms <= 0)
            return BoundaryOutcome<SourceInspectionReport>.Unavailable("A readable, identified source and positive intake limit are required.");
        if (source.Kind == SourceRouteKind.Upload &&
                (source.UploadProvenance is not { } origin || !Enum.IsDefined(origin)) ||
            source.Kind != SourceRouteKind.Upload && source.UploadProvenance is not null)
            return BoundaryOutcome<SourceInspectionReport>.Unavailable(
                "An uploaded structure needs an explicit, separate predicted, experimental, or unknown provenance declaration.");

        var request = new ScientificWorkRequest<SourceInspectionPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new SourceInspectionPayload(source.CoordinatePath, source.Sha256, maxAtoms,
                source.Kind, source.Prediction));
        var result = await _worker.InspectSourceAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<SourceInspectionReport>.Unavailable(result.FailureMessage ?? "Source inspection was not observed.");

        var observation = result.Observations;
        if (observation.Models.IsDefaultOrEmpty || observation.Models.Any(model => model.AtomCount <= 0) ||
            observation.Models.GroupBy(model => model.Index).Any(group => group.Count() > 1))
            return BoundaryOutcome<SourceInspectionReport>.Unavailable("The structural source has no unambiguous, nonempty model account.");
        if (source.Kind == SourceRouteKind.AlphaFold &&
            (source.Prediction is null || source.Prediction.CoordinateSha256 != source.Sha256 ||
             observation.Prediction is null ||
             observation.Prediction.RecordId != source.Prediction.RecordId ||
             observation.Prediction.CoordinateSha256 != source.Sha256))
            return BoundaryOutcome<SourceInspectionReport>.Unavailable(
                "The predicted coordinates and their source-attributed confidence account do not identify the same prediction.");
        if (source.Kind == SourceRouteKind.AlphaFold && !ValidPredictedSourceAccount(observation))
            return BoundaryOutcome<SourceInspectionReport>.Unavailable(
                "The prediction's residue-local confidence or relative-position evidence cannot be mapped or accounted for without ambiguity.");
        if (source.Kind != SourceRouteKind.AlphaFold && observation.Prediction is not null)
            return BoundaryOutcome<SourceInspectionReport>.Unavailable(
                "Prediction-confidence values cannot be attributed to an uploaded or experimental coordinate source.");

        return BoundaryOutcome<SourceInspectionReport>.Success(new SourceInspectionReport(
            source, observation.SourceFormat, observation.Models, ImmutableArray<string>.Empty,
            observation.Prediction));
    }

    public async Task<BoundaryOutcome<PreparationProposalReport>> ProposeChangesAsync(
        StudyRevision revision,
        IntendedProteinModel intended,
        SourceInspectionReport inspection,
        ProteinChemicalStatePolicy policy,
        ProteinStructuralAssessmentPolicy structuralPolicy,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        if (revision.IntendedProtein?.Id != intended.Id || inspection.Source.Id != intended.Source.Id ||
            inspection.Source.Sha256 != intended.Source.Sha256 || string.IsNullOrWhiteSpace(policy.Id) ||
            string.IsNullOrWhiteSpace(policy.Version) || policy.EvidenceReferences.IsDefaultOrEmpty ||
            !double.IsFinite(policy.DisulfideCandidateMaxSgDistanceAngstrom) ||
            policy.DisulfideCandidateMaxSgDistanceAngstrom <= 0)
            return BoundaryOutcome<PreparationProposalReport>.Unavailable("The intended protein, inspected source, and chemical-state policy must belong to the current study.");
        if (!ValidStructuralPolicy(structuralPolicy))
            return BoundaryOutcome<PreparationProposalReport>.Unavailable(
                "No identified, applicable protein structural assessment policy is available.");
        var model = inspection.Models.FirstOrDefault(item => item.Index == intended.ModelIndex);
        if (model is null || intended.Chains.IsDefaultOrEmpty ||
            intended.Chains.Any(chain => !model.Chains.Any(observed => observed.Name == chain.SourceChain)) ||
            model.Residues.Any(residue => intended.Chains.Any(chain => chain.SourceChain == residue.Address.Chain) &&
                residue.ResidueKind is not (SourceResidueKind.Protein or SourceResidueKind.Solvent or SourceResidueKind.Heterogen)) ||
            intended.Partners.Any(partner => !model.Partners.Any(observed => observed.SourceId == partner.SourceId)) ||
            intended.BiologicalAssemblyId is not null &&
            !model.Assemblies.Any(assembly => assembly.Name == intended.BiologicalAssemblyId &&
                intended.Chains.All(chain => assembly.ChainCopies.Contains(chain))))
            return BoundaryOutcome<PreparationProposalReport>.Unavailable("The selected model and protein chains have not been inspected.");

        var selectedResidues = SelectedCopyResidues(model, intended);
        if (selectedResidues.IsDefaultOrEmpty)
            return BoundaryOutcome<PreparationProposalReport>.Unavailable("The selected model has no inspectable retained residues.");
        if (intended.AlternateLocations.Any(choice => !string.IsNullOrEmpty(choice.Residue.CopyId) ||
                !selectedResidues.Any(residue => SourceAddress(residue.Address) == choice.Residue)))
            return BoundaryOutcome<PreparationProposalReport>.Unavailable("A source-level alternate location does not belong to the selected model and chain copies.");
        foreach (var residue in selectedResidues.Where(item => item.AlternateLocations.Length > 0))
            if (intended.AlternateLocations.Count(choice => choice.Residue == SourceAddress(residue.Address) &&
                    residue.AlternateLocations.Contains(choice.Altloc)) != 1)
                return BoundaryOutcome<PreparationProposalReport>.Unavailable("One observed alternate location must be provisionally selected for each affected residue before repair assessment.");

        var expandedAltlocs = selectedResidues.Where(item => item.AlternateLocations.Length > 0)
            .Select(residue => intended.AlternateLocations.Single(choice =>
                choice.Residue == SourceAddress(residue.Address)) with { Residue = residue.Address })
            .ToImmutableArray();

        var request = new ScientificWorkRequest<PreparationChangeInspectionPayload>(Guid.NewGuid().ToString("N"),
            workingDirectory, new PreparationChangeInspectionPayload(revision.Id,
                policy.DisulfideCandidateMaxSgDistanceAngstrom,
                intended.Source.CoordinatePath, intended.Source.Sha256, intended.ModelIndex,
                intended.BiologicalAssemblyId,
                intended.Chains, expandedAltlocs,
                intended.Partners.Where(partner => partner.Retain)
                    .Select(partner => model.Partners.Single(source => source.SourceId == partner.SourceId)).ToImmutableArray(),
                structuralPolicy.Measurement));
        var result = await _worker.InspectPreparationChangesAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<PreparationProposalReport>.Unavailable(result.FailureMessage ?? "Exact preparation-change assessment was not observed.");
        var observed = result.Observations;
        if (observed.PreviewChains.IsDefault ||
            observed.PreviewChains.Length != intended.Chains.Length ||
            observed.PreviewChains.Any(mapping =>
                string.IsNullOrWhiteSpace(mapping.SourceChain) ||
                string.IsNullOrWhiteSpace(mapping.CopyId) ||
                string.IsNullOrWhiteSpace(mapping.PreviewChain) ||
                mapping.PreviewChain.Length != 1 ||
                !intended.Chains.Any(chain => chain.SourceChain == mapping.SourceChain &&
                    chain.CopyId == mapping.CopyId)) ||
            observed.PreviewChains.Select(mapping => (mapping.SourceChain, mapping.CopyId))
                .Distinct().Count() != intended.Chains.Length ||
            observed.PreviewChains.Select(mapping => mapping.PreviewChain)
                .Distinct(StringComparer.Ordinal).Count() != intended.Chains.Length)
            return BoundaryOutcome<PreparationProposalReport>.Unavailable(
                "The exact selected chain copies cannot be mapped to the preview coordinates.");

        var changes = ImmutableArray.CreateBuilder<PreparationChangeProposal>();
        var unresolved = ImmutableArray.CreateBuilder<string>();
        if (observed.AssessmentStanding != ObservationStanding.Observed || observed.SelectedAtomCount <= 0)
            unresolved.Add("Missing-atom and disulfide assessment of the selected coordinates was not established.");
        foreach (var residue in selectedResidues)
        {
            if (!residue.BackboneHeavyAtomsComplete)
                unresolved.Add($"Backbone heavy-atom coordinates are missing at {residue.Address.Chain}:{residue.Address.Residue}{residue.Address.InsertionCode}.");

            if (residue.AlternateLocations.Length > 0)
            {
                var choice = intended.AlternateLocations.Single(candidate =>
                    candidate.Residue == SourceAddress(residue.Address));
                changes.Add(new PreparationChangeProposal(Guid.NewGuid().ToString("N"), revision.Id, intended.Id,
                    residue.Address, PreparationChangeKind.AlternateLocation, choice.Altloc, "Retain this provisionally selected observed coordinate location.",
                    residue.Limitations, true));
            }

            if (residue.Name == "HIS")
                foreach (var variant in new[] { "HID", "HIE", "HIP" }.Where(policy.PermittedVariants.Contains))
                    changes.Add(new PreparationChangeProposal(Guid.NewGuid().ToString("N"), revision.Id, intended.Id,
                        residue.Address, PreparationChangeKind.ResidueState, variant, "Choose a site-specific histidine state at the fixed study pH.",
                        residue.Limitations, true));
        }
        if (observed.AssessmentStanding == ObservationStanding.Observed)
            foreach (var atom in observed.MissingNonbackboneHeavyAtoms.Distinct())
            {
                if (!selectedResidues.Any(residue => residue.Address == atom.Residue))
                    unresolved.Add($"Missing atom {atom.AtomName} does not belong to a selected residue.");
                else
                    changes.Add(new PreparationChangeProposal(Guid.NewGuid().ToString("N"), revision.Id, intended.Id,
                        atom.Residue, PreparationChangeKind.HeavyAtom, atom.AtomName, "Add this exact missing nonbackbone heavy atom from the supported residue template.",
                        observed.Limitations, true));
            }
        foreach (var link in observed.PossibleDisulfides)
        {
            if (!selectedResidues.Any(residue => residue.Address == link.First) ||
                !selectedResidues.Any(residue => residue.Address == link.Second))
            {
                unresolved.Add("A candidate disulfide refers to a residue outside the selected construct.");
                continue;
            }
            changes.Add(new PreparationChangeProposal(Guid.NewGuid().ToString("N"), revision.Id,
                intended.Id, link.First, PreparationChangeKind.Disulfide,
                $"{link.First.Chain}:{link.First.Residue}–{link.Second.Chain}:{link.Second.Residue}",
                $"Review whether this exact cysteine pair forms a bond; observed SG separation {link.DistanceAngstrom:G6} Å.",
                observed.Limitations, true, link.Second));
        }
        var preview = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "selectedProteinPreview");
        if (preview is null || !File.Exists(preview.Path))
            unresolved.Add("The exact selected coordinates are unavailable for connected review of preparation changes.");
        var evidence = changes.Select(change => new ScientificEvidence(Guid.NewGuid().ToString("N"), change.Id,
            result.Provider?.Name ?? inspection.Source.Provenance,
            change.Kind == PreparationChangeKind.HeavyAtom ? "Selected-residue missing-atom assessment" :
                change.Kind == PreparationChangeKind.AlternateLocation ? "Observed coordinate alternative" :
                change.Kind == PreparationChangeKind.Disulfide ? "Observed cysteine proximity" : "Fixed-pH chemical-state proposal",
            $"{change.Kind}: {change.ProposedChange} at {change.Residue.Chain}:{change.Residue.Residue}{change.Residue.InsertionCode}",
            $"Study revision {revision.Id}; model {intended.ModelIndex}; policy {policy.Id}",
            change.Kind == PreparationChangeKind.ResidueState ? "Site-specific histidine state is a researcher-reviewed model assumption, not an automated optimal-protonation result." :
                change.Kind == PreparationChangeKind.Disulfide ? "Proximity suggests a possible link but does not establish a covalent bond; actual prepared bonds must match the decision." :
                "Selected-source observation; a proposal is not a prepared protein or placement.",
            change.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.Disulfide ? EvidenceBearing.Unknown : EvidenceBearing.Supports)).ToImmutableArray();
        return BoundaryOutcome<PreparationProposalReport>.Success(new PreparationProposalReport(
            revision.Id, intended.Id, changes.ToImmutable(), preview, observed.PreviewChains,
            evidence, unresolved.ToImmutable(), observed.Geometry));
    }

    public BoundaryOutcome<InspectionSubject> ProposalInspectionSubject(
        PreparationProposalReport report, string proposalId, string previewUrl)
    {
        var proposal = report.Changes.FirstOrDefault(change => change.Id == proposalId);
        var evidence = report.Evidence.FirstOrDefault(item => item.SubjectId == proposalId);
        if (proposal is null || evidence is null || report.Preview is null ||
            string.IsNullOrWhiteSpace(previewUrl))
            return BoundaryOutcome<InspectionSubject>.Unavailable("The exact preparation proposal has no attributable selected-coordinate preview and evidence.");
        var residue = proposal.Residue;
        var previewChain = PreviewChainFor(report, residue);
        var partnerPreviewChain = proposal.PartnerResidue is null ? null :
            PreviewChainFor(report, proposal.PartnerResidue);
        if (previewChain is null || (proposal.PartnerResidue is not null && partnerPreviewChain is null))
            return BoundaryOutcome<InspectionSubject>.Unavailable(
                "The proposal's exact selected chain copy is not identifiable in the preview.");
        var focus = new StructureFocus(previewChain, residue.Residue,
            string.IsNullOrWhiteSpace(residue.InsertionCode) ? null : residue.InsertionCode, null);
        var annotation = new InspectionAnnotation(Guid.NewGuid().ToString("N"),
            $"{residue.Chain}[{residue.CopyId}]:{residue.Residue}{residue.InsertionCode}",
            proposal.Kind.ToString(), proposal.Rationale, evidence.Id, focus);
        var metric = new InspectionMetric(proposal.Kind.ToString(), proposal.ProposedChange, null,
            annotation.SubjectPartId, evidence.Id);
        var annotations = ImmutableArray.CreateBuilder<InspectionAnnotation>();
        annotations.Add(annotation);
        if (proposal.PartnerResidue is not null)
        {
            var partner = proposal.PartnerResidue;
            annotations.Add(new InspectionAnnotation(Guid.NewGuid().ToString("N"),
                $"{partner.Chain}[{partner.CopyId}]:{partner.Residue}{partner.InsertionCode}", "possible disulfide partner",
                "Second observed cysteine in the proposed exact pair.", evidence.Id,
                new StructureFocus(partnerPreviewChain!, partner.Residue,
                    string.IsNullOrWhiteSpace(partner.InsertionCode) ? null : partner.InsertionCode, null)));
        }
        return BoundaryOutcome<InspectionSubject>.Success(new InspectionSubject(
            proposal.Id, report.StudyRevisionId, previewUrl, "selected-protein-before-repair",
            ImmutableArray<string>.Empty, ImmutableArray.Create(evidence), ImmutableArray<ScientificFinding>.Empty,
            annotations.ToImmutable(), ImmutableArray.Create(metric), null));
    }

    private static string? PreviewChainFor(PreparationProposalReport report, ResidueAddress address)
    {
        if (report.PreviewChains.IsDefaultOrEmpty || string.IsNullOrWhiteSpace(address.CopyId))
            return null;
        var matches = report.PreviewChains.Where(mapping =>
            mapping.SourceChain == address.Chain && mapping.CopyId == address.CopyId).Take(2).ToArray();
        return matches.Length == 1 && !string.IsNullOrWhiteSpace(matches[0].PreviewChain)
            ? matches[0].PreviewChain : null;
    }

    public async Task<BoundaryOutcome<AssessedPreparedProtein>> PrepareAsync(
        StudyRevision revision,
        IntendedProteinModel intended,
        SourceInspectionReport sourceInspection,
        ProteinChemicalStatePolicy policy,
        ProteinStructuralAssessmentPolicy structuralPolicy,
        PreparationProposalReport proposals,
        ImmutableArray<ResearcherDecision> decisions,
        string workingDirectory,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(revision);
        ArgumentNullException.ThrowIfNull(intended);
        ArgumentNullException.ThrowIfNull(sourceInspection);
        ArgumentNullException.ThrowIfNull(policy);
        ArgumentNullException.ThrowIfNull(proposals);

        if (revision.IntendedProtein?.Id != intended.Id || sourceInspection.Source.Id != intended.Source.Id ||
            sourceInspection.Source.Sha256 != intended.Source.Sha256)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The selected structure no longer corresponds to this study revision and inspected source.");
        if (proposals.StudyRevisionId != revision.Id || proposals.IntendedProteinId != intended.Id ||
            !proposals.UnresolvedQuestions.IsDefaultOrEmpty)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The exact preparation proposals are stale or contain unresolved structural questions.");
        if (string.IsNullOrWhiteSpace(policy.Id) || string.IsNullOrWhiteSpace(policy.Version) ||
            policy.EvidenceReferences.IsDefaultOrEmpty || policy.ForceFieldFiles.IsDefaultOrEmpty ||
            string.IsNullOrWhiteSpace(policy.ApplicableChemistry) ||
            !double.IsFinite(policy.DisulfideCandidateMaxSgDistanceAngstrom) ||
            policy.DisulfideCandidateMaxSgDistanceAngstrom <= 0 ||
            policy.ForceFieldFiles.GroupBy(asset => asset.Path, StringComparer.Ordinal)
                .Any(group => group.Select(asset => asset.Sha256).Distinct(StringComparer.OrdinalIgnoreCase).Count() != 1) ||
            policy.ForceFieldFiles.GroupBy(asset => asset.Sha256, StringComparer.OrdinalIgnoreCase)
                .Any(group => group.Select(asset => (asset.Family, asset.Version)).Distinct().Count() != 1))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("No identified, applicable protein chemical-state policy is available.");
        if (!ValidStructuralPolicy(structuralPolicy))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(
                "No identified, applicable protein structural assessment policy is available.");
        var forceFieldFiles = policy.ForceFieldFiles.DistinctBy(asset => asset.Sha256,
            StringComparer.OrdinalIgnoreCase).ToImmutableArray();
        if (intended.Chains.IsDefaultOrEmpty || intended.Chains.Any(chain => string.IsNullOrWhiteSpace(chain.CopyId)) ||
            intended.Chains.Select(chain => chain.CopyId).Distinct(StringComparer.Ordinal).Count() != intended.Chains.Length)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Selected protein chain copies are absent or ambiguous.");

        var model = sourceInspection.Models.FirstOrDefault(candidate => candidate.Index == intended.ModelIndex);
        if (model is null || intended.Chains.Any(chain => !model.Chains.Any(source => source.Name == chain.SourceChain)) ||
            model.Residues.Any(residue => intended.Chains.Any(chain => chain.SourceChain == residue.Address.Chain) &&
                residue.ResidueKind is not (SourceResidueKind.Protein or SourceResidueKind.Solvent or SourceResidueKind.Heterogen)) ||
            (intended.BiologicalAssemblyId is not null && !model.Assemblies.Any(assembly =>
                assembly.Name == intended.BiologicalAssemblyId &&
                intended.Chains.All(chain => assembly.ChainCopies.Contains(chain)))) ||
            intended.Partners.Any(partner => !model.Partners.Any(source => source.SourceId == partner.SourceId)))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The selected model, assembly, chain, or partner does not occur in the inspected source.");
        var selectedResidues = SelectedCopyResidues(model, intended);
        if (selectedResidues.IsDefaultOrEmpty || selectedResidues.Any(residue => !residue.BackboneHeavyAtomsComplete))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Retained protein residues require their observed backbone heavy-atom coordinates.");
        var approvedProposals = new List<(PreparationChangeProposal Proposal, ResearcherDecision Decision)>();
        foreach (var proposal in proposals.Changes)
        {
            if (proposal.StudyRevisionId != revision.Id || proposal.IntendedProteinId != intended.Id)
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("A preparation proposal is not bound to the selected study and protein.");
            var matching = decisions.Where(decision => decision.StudyRevisionId == revision.Id &&
                decision.SubjectId == proposal.Id).ToArray();
            if (matching.Length > 1 || matching.Any(decision => decision.Kind != ResearcherDecisionKind.ApprovePreparationChange ||
                    decision.ChosenValue is not (ResearcherDecisionValue.Approved or ResearcherDecisionValue.Declined)))
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("A preparation change has conflicting or invalid decisions.");
            if ((proposal.Kind is PreparationChangeKind.ResidueState or PreparationChangeKind.Disulfide) && matching.Any(decision =>
                    decision.ChosenValue == ResearcherDecisionValue.Approved && string.IsNullOrWhiteSpace(decision.Rationale)))
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("A site-specific chemical-state or disulfide choice requires the researcher's recorded rationale.");
            if (proposal.Kind == PreparationChangeKind.Disulfide && matching.Length != 1)
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Each possible disulfide must be explicitly approved or declined before preparation.");
            if (matching.Length == 1 && matching[0].ChosenValue == ResearcherDecisionValue.Approved)
                approvedProposals.Add((proposal, matching[0]));
        }
        if (proposals.Changes.Where(proposal => proposal.Kind == PreparationChangeKind.HeavyAtom).Any(proposal =>
                !approvedProposals.Any(approved => approved.Proposal.Id == proposal.Id)))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Every proposed heavy-atom addition requires exact approval before preparation.");
        foreach (var group in proposals.Changes.Where(proposal => proposal.Kind is PreparationChangeKind.AlternateLocation or PreparationChangeKind.ResidueState)
                     .GroupBy(proposal => (proposal.Residue, proposal.Kind)))
            if (approvedProposals.Count(approved => group.Any(proposal => proposal.Id == approved.Proposal.Id)) != 1)
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Exactly one reviewed alternate location or residue state is required per affected residue.");

        var heavyAtomApprovals = approvedProposals.Where(item => item.Proposal.Kind == PreparationChangeKind.HeavyAtom)
            .Select(item => new HeavyAtomApproval(item.Proposal.Residue, item.Proposal.ProposedChange, item.Decision.Id))
            .ToImmutableArray();
        var altlocChoices = approvedProposals.Where(item => item.Proposal.Kind == PreparationChangeKind.AlternateLocation)
            .Select(item => new AlternateLocationChoice(item.Proposal.Residue, item.Proposal.ProposedChange, item.Decision.Id))
            .ToImmutableArray();
        var residueVariants = approvedProposals.Where(item => item.Proposal.Kind == PreparationChangeKind.ResidueState)
            .Select(item => new ResidueVariantChoice(item.Proposal.Residue, item.Proposal.ProposedChange, item.Decision.Id))
            .ToImmutableArray();
        var approvedDisulfides = approvedProposals.Where(item => item.Proposal.Kind == PreparationChangeKind.Disulfide &&
                item.Proposal.PartnerResidue is not null)
            .Select(item => new DisulfideChoice(item.Proposal.Residue, item.Proposal.PartnerResidue!, item.Decision.Id))
            .ToImmutableArray();
        if (approvedDisulfides.SelectMany(choice => new[] { choice.First, choice.Second }).Distinct().Count() !=
            approvedDisulfides.Length * 2)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("A cysteine cannot belong to more than one approved disulfide pair.");
        if (residueVariants.Any(choice => !policy.PermittedVariants.Contains(choice.Variant)) ||
            intended.AlternateLocations.Any(choice => selectedResidues
                .Where(residue => SourceAddress(residue.Address) == choice.Residue)
                .Any(residue => !altlocChoices.Any(approved =>
                    approved.Residue == residue.Address && approved.Altloc == choice.Altloc))))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("An intended alternate location or residue state lacks a valid reviewed choice.");

        var effectiveVariants = ImmutableArray.CreateBuilder<ResidueVariantChoice>();
        foreach (var residue in selectedResidues)
        {
            var disulfide = approvedDisulfides.FirstOrDefault(choice =>
                choice.First == residue.Address || choice.Second == residue.Address);
            if (disulfide is not null)
            {
                if (!policy.PermittedVariants.Contains("CYX"))
                    return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The identified chemical-state policy cannot represent an approved disulfide cysteine.");
                effectiveVariants.Add(new ResidueVariantChoice(residue.Address, "CYX", disulfide.DecisionId));
                continue;
            }
            var chosen = residueVariants.FirstOrDefault(choice => choice.Residue == residue.Address);
            if (chosen is not null)
            {
                effectiveVariants.Add(chosen);
                continue;
            }
            if (policy.DefaultVariants.TryGetValue(residue.Name, out var ordinary))
            {
                effectiveVariants.Add(new ResidueVariantChoice(residue.Address, ordinary, null));
                continue;
            }
            if (residue.Name == "HIS")
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Each histidine needs an explicit, reviewed chemical-state choice.");
        }

        var request = new ScientificWorkRequest<ProteinPreparationPayload>(
            Guid.NewGuid().ToString("N"), workingDirectory,
            new ProteinPreparationPayload(
                revision.Id, revision.Conditions.NominalPh,
                policy.DisulfideCandidateMaxSgDistanceAngstrom,
                intended.Source.CoordinatePath, intended.Source.Sha256,
                intended.ModelIndex, intended.BiologicalAssemblyId,
                intended.Chains,
                altlocChoices,
                intended.Partners.Where(partner => partner.Retain).Select(partner =>
                    model.Partners.Single(source => source.SourceId == partner.SourceId)).ToImmutableArray(),
                heavyAtomApprovals, approvedDisulfides, effectiveVariants.ToImmutable(), forceFieldFiles,
                structuralPolicy.Measurement));
        var result = await _worker.PrepareProteinAsync(request, cancellationToken);
        if (result.RequestId != request.RequestId || result.StudyRevisionId != revision.Id ||
            result.Standing != WorkerResultStanding.Observed || result.Observations is null)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(result.FailureMessage ?? "Protein preparation was not observed.");

        var observed = result.Observations;
        if (!observed.MissingBackboneResidues.IsDefaultOrEmpty || !observed.NoncanonicalResidues.IsDefaultOrEmpty ||
            !observed.UnresolvedAlternateLocations.IsDefaultOrEmpty || !observed.UnparameterizedResidues.IsDefaultOrEmpty)
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The exact retained backbone, residue chemistry, alternate selection, or parameter route is unsupported.");
        var approved = heavyAtomApprovals.Select(approval => new AtomAddress(approval.Residue, approval.AtomName)).ToHashSet();
        if (heavyAtomApprovals.Length != approved.Count || observed.AddedHeavyAtoms.IsDefault ||
            observed.AddedHeavyAtoms.Length != approved.Count ||
            observed.AddedHeavyAtoms.Distinct().Count() != approved.Count ||
            observed.AddedHeavyAtoms.Any(atom => !approved.Contains(atom)))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(
                "The actual non-hydrogen heavy-atom additions do not match the exact approvals.");
        var expectedVariants = effectiveVariants.Select(choice => (choice.Residue, choice.Variant)).ToHashSet();
        if (effectiveVariants.Count != expectedVariants.Count || observed.ActualResidueVariants.IsDefault ||
            observed.ActualResidueVariants.Length != expectedVariants.Count ||
            observed.ActualResidueVariants.Select(actual => (actual.Residue, actual.Variant)).Distinct().Count() !=
                expectedVariants.Count ||
            observed.ActualResidueVariants.Any(actual => !expectedVariants.Contains((actual.Residue, actual.Variant))))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(
                "The actual residue states do not match the selected explicit variants.");
        var approvedPairs = approvedDisulfides.Select(choice => new DisulfideBond(choice.First, choice.Second))
            .ToImmutableArray();
        if (observed.ActualDisulfides.Length != approvedPairs.Length ||
            observed.ActualDisulfides.Any(actual => !approvedPairs.Any(approved =>
                approved == actual || approved.First == actual.Second && approved.Second == actual.First)))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("Actual prepared disulfide bonds do not match the exact approved pairs.");

        var prepared = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "preparedPdb");
        var bondGraph = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "preparedBondGraph");
        var correspondenceArtifact = result.Artifacts.FirstOrDefault(artifact => artifact.Role == "correspondenceJson");
        if (prepared is null || bondGraph is null || correspondenceArtifact is null ||
            !File.Exists(bondGraph.Path) || !File.Exists(correspondenceArtifact.Path))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The corresponding prepared structure, explicit bonds, or atom mapping was not observed.");

        SourceToResultCorrespondence? correspondence;
        try
        {
            correspondence = JsonSerializer.Deserialize<SourceToResultCorrespondence>(
                await File.ReadAllTextAsync(correspondenceArtifact.Path, cancellationToken),
                new JsonSerializerOptions(JsonSerializerDefaults.Web));
        }
        catch (JsonException)
        {
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The prepared atom correspondence has unrecognized values.");
        }
        if (correspondence is null || !correspondence.Complete ||
            correspondence.SourceId != intended.Source.Sha256 || correspondence.ResultId != prepared.Sha256 ||
            correspondence.Atoms.Length != observed.PreparedAtomCount ||
            observed.CorrespondedResultAtomCount != observed.PreparedAtomCount ||
            correspondence.Atoms.Select(atom => atom.ResultAtomId).Distinct(StringComparer.Ordinal).Count() != observed.PreparedAtomCount ||
            correspondence.Atoms.Any(atom => !Enum.IsDefined(atom.Role) || !Enum.IsDefined(atom.MoleculeRole) ||
                !Enum.IsDefined(atom.AtomRole) ||
                atom.MoleculeRole is not (MoleculeRoleKind.Protein or MoleculeRoleKind.RetainedPartner) ||
                (atom.MoleculeRole == MoleculeRoleKind.Protein &&
                    atom.AtomRole is not (AtomRoleKind.Backbone or AtomRoleKind.Sidechain)) ||
                (atom.MoleculeRole == MoleculeRoleKind.RetainedPartner &&
                    atom.AtomRole != AtomRoleKind.PartnerAtom) ||
                atom.GeneratedComponentRole is not null || atom.PhysicalSide is not null))
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable("The prepared structure cannot be traced atom-for-atom to its source and approved additions.");

        var molecule = new MolecularArtifact(
            correspondence.ResultId, prepared.Path, prepared.Sha256,
            bondGraph.Path, null, null, observed.PreparedAtomCount, null, bondGraph.Sha256);
        var evidence = ImmutableArray.Create(new ScientificEvidence(
            Guid.NewGuid().ToString("N"), molecule.Id, result.Provider?.Name ?? "local scientific worker",
            "Observed structural preparation", $"{observed.PreparedAtomCount} atoms and {observed.RetainedResidueCount} retained residues",
            $"Selected source {intended.Source.Id}, model {intended.ModelIndex} and policy {policy.Id}",
            "Protein-local preparation only; placement and assembled-system compatibility remain unassessed.",
            EvidenceBearing.Context));
        if (!GeometrySupports(observed.Geometry, structuralPolicy) ||
            !observed.GeometryWarnings.IsDefaultOrEmpty)
        {
            if (observed.Geometry is null)
                return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(
                    "The actual prepared candidate has no structural-geometry observation; support cannot be established.");
            var kinds = observed.Geometry.Kinds.IsDefault
                ? ImmutableArray<ProteinGeometryKindObservation>.Empty : observed.Geometry.Kinds;
            var geometryEvidence = new ScientificEvidence(Guid.NewGuid().ToString("N"), molecule.Id,
                result.Provider?.Name ?? "local scientific worker", "Measured candidate protein geometry",
                string.Join("; ", kinds.Select(kind =>
                    $"{kind.Kind}: {kind.Standing}, {kind.MeasuredCount}/{kind.EligibleCount} measured")),
                $"Prepared candidate {molecule.Id}; structural policy {structuralPolicy.Id} version {structuralPolicy.Version}",
                "Actual measurements require policy interpretation; an unpromoted candidate is not an assessed protein.",
                EvidenceBearing.Context);
            var geometryFindings = kinds
                .Where(kind => structuralPolicy.Measurement.RequiredKinds.Contains(kind.Kind) &&
                    !GeometryKindSupports(kind, structuralPolicy.Criteria.Single(item => item.Kind == kind.Kind)))
                .Select(kind => new ScientificFinding(Guid.NewGuid().ToString("N"), molecule.Id,
                    geometryEvidence.Id,
                    $"{kind.Kind}: standing {kind.Standing}; observed {kind.MeasuredCount}/{kind.EligibleCount}; " +
                    $"distance range {kind.MinimumDistanceAngstrom?.ToString("G6") ?? "unavailable"}–" +
                    $"{kind.MaximumDistanceAngstrom?.ToString("G6") ?? "unavailable"} Å. " +
                    (kind.UnavailableReason ?? string.Empty),
                    "The actual candidate's geometry does not establish this product's bounded preparation support.",
                    kind.Standing == GeometryKindStanding.Observed ? FindingDisposition.Disqualifies : FindingDisposition.Challenges,
                    true, DateTimeOffset.UtcNow)).ToImmutableArray();
            if (!observed.GeometryWarnings.IsDefaultOrEmpty)
                geometryFindings = geometryFindings.AddRange(observed.GeometryWarnings.Select(warning =>
                    new ScientificFinding(Guid.NewGuid().ToString("N"), molecule.Id, geometryEvidence.Id,
                        warning, "Required candidate geometry observation is not established.",
                        FindingDisposition.Challenges, true, DateTimeOffset.UtcNow)));
            var geometryReason = !observed.GeometryWarnings.IsDefaultOrEmpty
                ? observed.GeometryWarnings[0]
                : kinds.Select(kind => kind.UnavailableReason).FirstOrDefault(reason => !string.IsNullOrWhiteSpace(reason));
            var diagnostic = new ProteinPreparationDiagnostic(molecule, observed.Geometry,
                correspondence, evidence.Add(geometryEvidence), geometryFindings);
            return BoundaryOutcome<AssessedPreparedProtein>.Unavailable(
                "Measured geometry of the actual prepared candidate is unavailable, incomplete, or outside the applicable structural policy." +
                    (geometryReason is null ? string.Empty : $" {geometryReason}"),
                geometryFindings, diagnostic);
        }
        return BoundaryOutcome<AssessedPreparedProtein>.Success(new AssessedPreparedProtein(
            molecule.Id, revision.Id, intended, molecule, policy.Id, effectiveVariants.ToImmutable(),
            approvedProposals.Select(item => item.Proposal).ToImmutableArray(), correspondence,
            evidence, ImmutableArray<ScientificFinding>.Empty,
            ImmutableArray.Create("Protein-local assessment does not establish membrane placement."),
            sourceInspection.Prediction, observed.Geometry,
            policy.Version, structuralPolicy.Id, structuralPolicy.Version));
    }

    private static bool ValidStructuralPolicy(ProteinStructuralAssessmentPolicy? policy)
    {
        if (policy is null || string.IsNullOrWhiteSpace(policy.Id) ||
            string.IsNullOrWhiteSpace(policy.Version) ||
            string.IsNullOrWhiteSpace(policy.ApplicableProteinClass) ||
            policy.EvidenceReferences.IsDefaultOrEmpty ||
            policy.Measurement is null || policy.Criteria.IsDefaultOrEmpty ||
            policy.Measurement.RequiredKinds.IsDefaultOrEmpty ||
            policy.Measurement.RequiredKinds.Distinct(StringComparer.Ordinal).Count() !=
                policy.Measurement.RequiredKinds.Length ||
            !double.IsFinite(policy.Measurement.NeighborSearchRadiusAngstrom) ||
            policy.Measurement.NeighborSearchRadiusAngstrom <= 0 ||
            policy.Measurement.ExcludedBondHops < 0 ||
            policy.Measurement.MaximumReportedPairs <= 0 ||
            policy.Measurement.AtomRadiusByElementAngstrom.IsEmpty ||
            policy.Measurement.AtomRadiusByElementAngstrom.Any(item =>
                string.IsNullOrWhiteSpace(item.Key) || !double.IsFinite(item.Value) || item.Value <= 0))
            return false;
        return policy.Measurement.RequiredKinds.All(kind =>
            policy.Criteria.Count(criterion => criterion.Kind == kind &&
                (criterion.MinimumObservedAngstrom is not null ||
                 criterion.MaximumObservedAngstrom is not null) &&
                (criterion.MinimumObservedAngstrom is null ||
                 double.IsFinite(criterion.MinimumObservedAngstrom.Value)) &&
                (criterion.MaximumObservedAngstrom is null ||
                 double.IsFinite(criterion.MaximumObservedAngstrom.Value)) &&
                (criterion.MinimumObservedAngstrom is null || criterion.MaximumObservedAngstrom is null ||
                 criterion.MinimumObservedAngstrom <= criterion.MaximumObservedAngstrom)) == 1);
    }

    private static bool ValidPredictedSourceAccount(SourceInspectionObservations observed)
    {
        var prediction = observed.Prediction;
        if (prediction is null || prediction.LocalConfidence.IsDefault ||
            prediction.PaeStanding is not (PredictionObservationStanding.Observed or PredictionObservationStanding.Unavailable or PredictionObservationStanding.Unmapped) ||
            prediction.LocalConfidence.Select(item => item.Residue).Distinct().Count() !=
                prediction.LocalConfidence.Length)
            return false;
        var proteinResidues = observed.Models.SelectMany(model => model.Residues)
            .Where(residue => residue.ResidueKind == SourceResidueKind.Protein)
            .Select(residue => residue.Address).ToHashSet();
        if (proteinResidues.Count == 0 || prediction.LocalConfidence.Length != proteinResidues.Count ||
            prediction.LocalConfidence.Any(item => !proteinResidues.Contains(item.Residue) ||
                item.Standing is not (PredictionObservationStanding.Observed or PredictionObservationStanding.Unavailable or PredictionObservationStanding.Unmapped) ||
                item.Standing == PredictionObservationStanding.Observed &&
                (item.PLddt is not double score || !double.IsFinite(score) || score < 0 || score > 100) ||
                item.Standing != PredictionObservationStanding.Observed && item.PLddt is not null))
            return false;
        return prediction.PaeStanding != PredictionObservationStanding.Observed ||
            prediction.PaeAxisResidueCount is > 0 &&
            !string.IsNullOrWhiteSpace(prediction.PaeMappingPath) &&
            !string.IsNullOrWhiteSpace(prediction.PaeMappingSha256);
    }

    private static bool GeometrySupports(ProteinGeometryObservations? observed,
        ProteinStructuralAssessmentPolicy policy)
    {
        if (observed is null || observed.Standing != ObservationStanding.Observed || observed.Kinds.IsDefault ||
            observed.Kinds.Select(kind => kind.Kind).Distinct(StringComparer.Ordinal).Count() !=
                observed.Kinds.Length)
            return false;
        foreach (var required in policy.Measurement.RequiredKinds)
        {
            var kind = observed.Kinds.FirstOrDefault(item => item.Kind == required);
            var criterion = policy.Criteria.Single(item => item.Kind == required);
            if (kind is null || !GeometryKindSupports(kind, criterion))
                return false;
        }
        return true;
    }

    private static bool GeometryKindSupports(ProteinGeometryKindObservation kind,
        ProteinGeometryCriterion criterion)
    {
        if (kind.Standing == GeometryKindStanding.NotApplicable && criterion.AllowNotApplicable &&
            kind.EligibleCount == 0 && kind.MeasuredCount == 0 &&
            kind.MinimumDistanceAngstrom is null && kind.MaximumDistanceAngstrom is null)
            return true;
        if (kind.Standing != GeometryKindStanding.Observed || kind.EligibleCount <= 0 ||
            kind.MeasuredCount < 0 || kind.MeasuredCount > kind.EligibleCount ||
            kind.MinimumDistanceAngstrom is not double minimum ||
            kind.MaximumDistanceAngstrom is not double maximum ||
            !double.IsFinite(minimum) || !double.IsFinite(maximum) || minimum > maximum ||
            kind.MeasuredCount != kind.EligibleCount ||
            criterion.MinimumObservedAngstrom is double lower && minimum < lower ||
            criterion.MaximumObservedAngstrom is double upper && maximum > upper)
            return false;
        return true;
    }

    private static ImmutableArray<SourceResidueObservation> SelectedCopyResidues(
        SourceModelObservation model, IntendedProteinModel intended)
    {
        var selected = ImmutableArray.CreateBuilder<SourceResidueObservation>();
        foreach (var chain in intended.Chains)
            foreach (var residue in model.Residues.Where(item => item.Address.Chain == chain.SourceChain &&
                         item.ResidueKind == SourceResidueKind.Protein))
                selected.Add(residue with { Address = residue.Address with { CopyId = chain.CopyId } });
        return selected.ToImmutable();
    }

    private static ResidueAddress SourceAddress(ResidueAddress address) => address with { CopyId = string.Empty };
}
