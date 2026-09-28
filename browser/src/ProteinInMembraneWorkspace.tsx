import { useEffect, useRef, useState, type ReactNode } from 'react';
import { ConnectedStructuralInspection, viewerSubjectHeading, type InspectionAccount, type StructureLoadStatus } from './ConnectedStructuralInspection';

type SourceRouteKind = 'rcsb' | 'alphafold' | 'upload';
type UploadOriginKind = 'predicted' | 'experimental' | 'unknown';
type PreparationChangeKind = 'alternateLocation' | 'residueState' | 'heavyAtom' | 'disulfide';
type ProteinTopologyKind = 'membrane-spanning' | 'one-surface-associated';
type LeafletSide = 'upper' | 'lower';
type PlacementPhysicalSide = LeafletSide | 'both';
type PpmNterminalSide = 'in' | 'out';
type PlacementStartingPosition = 'center' | 'upper' | 'lower';
type SourceResidueKind = 'protein' | 'solvent' | 'heterogen';
type PredictionObservationStanding = 'Observed' | 'Unavailable' | 'Unmapped';
type ObservationStanding = 'Observed' | 'Unavailable';
type GeometryKindStanding = ObservationStanding | 'NotApplicable';
type StageKind = 'Minimization' | 'Equilibration';
type WorkArea = 'protein' | 'membrane' | 'placement' | 'preparation' | 'results';
type PendingAction = ActorActionKind | 'uploadSource';
interface ActionNotice { kind: PendingAction; tone: 'success' | 'warning' | 'error'; message: string; }
interface SourceRequest { serial: number; key: string; label: string; data: object; }
interface SourceIntent extends SourceRequest { phase: 'retrieving' | 'rendering' | 'failed'; reason?: string; }
interface SourcePreviewAccount { sourceId: string; modelIndex: number; assemblyId: string | null;
  structureUrl: string; chainIds: string[]; }
interface SourcePreviewState { key: string; phase: 'loading' | 'ready' | 'failed';
  account?: SourcePreviewAccount; reason?: string; }
const workAreas: { id: WorkArea; label: string; description: string }[] = [
  { id: 'protein', label: 'Protein', description: 'Discover, inspect, and prepare a structure' },
  { id: 'membrane', label: 'Membrane', description: 'Define and assess the intended bilayer' },
  { id: 'placement', label: 'Placement', description: 'Assess and review the protein position' },
  { id: 'preparation', label: 'Preparation', description: 'Construct and minimize an explicit system' },
  { id: 'results', label: 'Results', description: 'Inspect and export completed stages' },
];
const activityText: Partial<Record<PendingAction, string>> = {
  searchSource: 'Searching structural databases…',
  selectSource: 'Retrieving source…',
  uploadSource: 'Uploading structure…',
  selectProteinModel: 'Assessing the selected protein model…',
  approvePreparationChange: 'Recording the exact choice; checking the protein if review is complete…',
  retryProteinPreparation: 'Retrying protein preparation and assessment…',
  authorizePreparationPlan: 'Applying the checked preparation plan and assessing the protein…',
  overridePreparationPlanChoice: 'Checking the updated preparation plan…',
  retryPreparationPlan: 'Finding preparation suggestions…',
  startProteinPreparation: 'Preparing and checking the selected protein…',
  proposeMembrane: 'Recording the intended leaflet composition…',
  adoptMembrane: 'Assessing the chosen membrane model…',
  proposePlacement: 'Assessing protein–membrane placement…',
  revisePlacement: 'Reassessing the adjusted placement…',
  adoptPlacement: 'Adopting the supported placement…',
  startPreparation: 'Starting system construction…',
  continueMinimization: 'Starting minimization…',
  stopAttempt: 'Submitting the stop request…',
  exportStage: 'Verifying and transferring the completed stage…',
  requestEquilibration: 'Requesting optional equilibration…',
  selectInspectionSubject: 'Loading the selected subject…',
};
export type ActorActionKind =
  | 'searchSource' | 'selectSource' | 'selectProteinModel'
  | 'approvePreparationChange' | 'declinePreparationChange'
  | 'retryProteinPreparation'
  | 'authorizePreparationPlan' | 'overridePreparationPlanChoice' | 'retryPreparationPlan'
  | 'startProteinPreparation'
  | 'proposeMembrane' | 'adoptMembrane' | 'proposePlacement'
  | 'revisePlacement' | 'adoptPlacement' | 'startPreparation'
  | 'continueMinimization' | 'stopAttempt' | 'requestEquilibration' | 'selectInspectionSubject'
  | 'setInspectionFocus' | 'exportStage';
type ActorCommandKind = Exclude<ActorActionKind, 'declinePreparationChange'>;

interface FixedStudyConditions {
  nominalPh: number;
  targetNaClMolar: number;
  optionalTemperatureKelvin: number;
}

interface StudyAccount {
  id: string;
  number: number;
  summary: string;
  selectedSourceId: string | null;
  selectedSourceLabel?: string | null;
  selectedSourceKind?: SourceRouteKind | null;
  uploadProvenance?: UploadOriginKind | null;
  uploadProvenanceNote?: string | null;
  modelIndex: number | null;
  adoptedPlacementProposalId: string | null;
  biologicalAssemblyId: string | null;
  chainIds: string[];
  partners: PartnerSelection[];
  alternateLocations: AlternateLocationChoice[];
  conditions: FixedStudyConditions;
}

interface SourceCandidateAccount {
  id: string;
  label: string;
  kind: SourceRouteKind;
  provenance: string;
  limitations: string[];
}

interface ChainSelection { sourceChain: string; copyId: string; }
interface PartnerSelection { sourceId: string; retain: boolean; reason?: string | null; }
interface ResidueAddress { model: number; chain: string; residue: number; insertionCode: string; copyId: string; }
interface AlternateLocationChoice { residue: ResidueAddress; altloc: string; decisionId: string; }
interface SourceResidueObservation { address: ResidueAddress; name: string; residueKind: SourceResidueKind; backboneHeavyAtomsComplete: boolean; alternateLocations: string[]; }
interface SourceModelObservation {
  index: number;
  sourceModelId?: string | null;
  chains: { name: string; residueCount: number; atomCount: number }[];
  assemblies: { name: string; chainCopies: ChainSelection[] }[];
  partners: SourcePartnerObservation[];
  residues: SourceResidueObservation[];
  atomCount: number;
}
interface SourcePartnerObservation {
  sourceId: string; label: string; kind: string; atomCount: number;
  chain: string | null; residue: number | null;
  insertionCode?: string | null; displayName?: string | null;
}

interface PredictedResidueConfidence {
  residue: ResidueAddress;
  pLddt: number | null;
  standing: PredictionObservationStanding;
  reason: string | null;
}
interface PredictionEvidenceAccount {
  recordId: string;
  localConfidence: PredictedResidueConfidence[];
  paeStanding: PredictionObservationStanding;
  paeReason: string | null;
  paeAxisResidueCount: number | null;
  limitations: string[];
}
interface ProteinGeometryKindObservation {
  kind: string;
  standing: GeometryKindStanding;
  eligibleCount: number;
  measuredCount: number;
  minimumDistanceAngstrom: number | null;
  maximumDistanceAngstrom: number | null;
  unavailableReason: string | null;
}
interface ProteinGeometryObservations {
  standing: ObservationStanding;
  kinds: ProteinGeometryKindObservation[];
  limitations: string[];
}
interface DirectionalPredictionSummary {
  possiblePairCount: number;
  validPairCount: number;
  minimumAngstrom: number | null;
  maximumAngstrom: number | null;
  meanAngstrom: number | null;
}
interface PredictionRegionSummaryObservations {
  recordId: string;
  standing: PredictionObservationStanding;
  reason: string | null;
  firstAlignedOnSecond: DirectionalPredictionSummary | null;
  secondAlignedOnFirst: DirectionalPredictionSummary | null;
}

interface LipidCatalogueAccount {
  speciesId: string;
  displayName: string;
  chemistryId: string;
  limitations: string[];
}

interface PreparationChangeProposal {
  id: string;
  studyRevisionId: string;
  intendedProteinId: string;
  kind: PreparationChangeKind;
  proposedChange: string;
  rationale: string;
  limitations: string[];
  approvalRequired: boolean;
  residue: ResidueAddress;
  partnerResidue?: ResidueAddress | null;
}

interface PreparationDecisionOptionAccount {
  proposalId: string;
  proposedChange: string;
  disposition: 'available' | 'confirmed' | 'declined' | 'notChosen';
  decisionId: string | null;
  evidenceCurrent: boolean;
  confirmationBlocker: string | null;
  startsPreparationOnConfirmation: boolean;
  startsPreparationOnDecline: boolean;
  evidence: ScientificEvidence[];
  informationRole?: 'modelAssumption' | 'observed';
}
interface PreparationDecisionAccount {
  id: string;
  kind: PreparationChangeKind;
  residue: ResidueAddress;
  partnerResidue: ResidueAddress | null;
  atomName: string | null;
  standing: 'pending' | 'confirmed' | 'blocked';
  chosenProposalId: string | null;
  options: PreparationDecisionOptionAccount[];
  blocker: string | null;
}
interface ProteinPreparationReviewAccount {
  studyRevisionId: string;
  intendedProteinId: string;
  decisions: PreparationDecisionAccount[];
  confirmedCount: number;
  remainingCount: number;
  blockers: string[];
  preparationStanding: 'awaitingDecisions' | 'blocked' | 'ready' | 'preparing' | 'failed' | 'assessed';
  preparationMessage: string | null;
  inspectionRelation?: 'beforePreparation' | 'preparedResult' | 'unqualifiedCandidate' | 'historical' | 'otherSubject' | 'unavailable';
}

interface PreparationPlanAccount {
  standing: 'calculating' | 'ready' | 'partial' | 'failed' | 'authorized' | 'preparing' | 'applied' | 'failedAfterAuthorization' | 'unavailable';
  planId: string | null;
  planSha256: string | null;
  method: string | null;
  nominalPh: number | null;
  stateChoiceCount: number;
  heavyAtomCount: number;
  removedSourceHydrogenCount: number;
  choices: { residue: ResidueAddress; variant: string; overridden: boolean }[];
  message: string | null;
  authorizationAvailable: boolean;
  hasOverrides: boolean;
}

interface ProteinTaskAccount {
  studyRevisionId: string;
  intendedProteinId: string;
  sourceId: string;
  chains: ChainSelection[];
  standing: 'assessing' | 'awaitingDecisions' | 'planReady' | 'ready' | 'blocked' | 'preparing' | 'assessed' | 'failed' | 'unavailable';
  message: string | null;
  preparedProteinId: string | null;
  retryAvailable: boolean;
}
interface PlacementTaskAccount {
  studyRevisionId: string;
  preparedProteinId: string;
  membraneModelId: string;
  proposalId: string | null;
  standing: 'ready' | 'obtaining' | 'assessing' | 'noProposal' | 'review' | 'supported' | 'unsupported' | 'notEstablished';
  message: string;
  adopted: boolean;
  latestAttemptIssue: string | null;
}
interface PlacementMethodAccount { method: 'OPM' | 'PPM'; standing: 'lookupEligible' | 'notApplicable' | 'configured' | 'unavailable'; reason: string | null; }

export interface ScientificFinding { id: string; subjectId: string; evidenceId: string; meaning: string; consequence: string; disposition: string; material: boolean; }
export interface ScientificEvidence { id: string; subjectId: string; source: string; method: string; observation: string; applicability: string; uncertainty: string; bearing: string; }

interface ProteinAccount {
  subjectId: string;
  status: string;
  summary: string;
  atomCount: number | null;
  changes: PreparationChangeProposal[];
  findings: ScientificFinding[];
  candidateId?: string | null;
  prediction: PredictionEvidenceAccount | null;
  geometry: ProteinGeometryObservations | null;
  sourceGeometry: ProteinGeometryObservations | null;
}

interface LipidFraction { speciesId: string; fraction: number; }
interface MembraneSpeciesSupportAccount {
  speciesId: string;
  chemistryId: string;
  category: string;
  forceFieldFamily: string;
  forceFieldVersion: string;
  coordinateSha256: string;
  parameterSha256: string;
  limitations: string[];
}
interface MembraneAccount {
  modelId: string;
  status: string;
  upper: LipidFraction[];
  lower: LipidFraction[];
  scientificPurpose: string | null;
  limitations: string[];
  reason: string | null;
  policyId: string | null;
  policyVersion: string | null;
  evidence: ScientificEvidence[];
  speciesSupport: MembraneSpeciesSupportAccount[];
}

interface PlacementAccount {
  proposalId: string;
  status: string;
  preparedProteinId: string | null;
  membraneModelId: string | null;
  topologyKind: ProteinTopologyKind;
  physicalSide: PlacementPhysicalSide;
  midplaneAngstrom: number | null;
  thicknessAngstrom: number | null;
  depthAngstrom: number | null;
  tiltDegrees: number | null;
  sidedness: string | null;
  contactingRegions: string[];
  limitations: string[];
  policyId: string | null;
  policyVersion: string | null;
  witnessId: string | null;
  reason: string;
  evidence: ScientificEvidence[];
  prediction: PredictionRegionSummaryObservations | null;
  transform?: { startingPosition: PlacementStartingPosition;
    offsetXAngstrom: number; offsetYAngstrom: number; offsetZAngstrom: number;
    rotationXDegrees: number; rotationYDegrees: number; rotationZDegrees: number;
    appliedTranslationXAngstrom: number; appliedTranslationYAngstrom: number;
    appliedTranslationZAngstrom: number; headgroupBoundaryAngstrom: number } | null;
}

interface SpeciesCountAccount {
  physicalSide: LeafletSide;
  speciesId: string;
  count: number;
  intendedFraction: number;
}

interface MeasuredValueAccount { name: string; value: number; unit: string; scope: string; }
interface LocalStateAccount {
  standing: string;
  unavailableReason: string | null;
  measurements: MeasuredValueAccount[];
  limitations: string[];
  rolePairMeasurements: {
    firstMoleculeRole: string;
    secondMoleculeRole: string;
    pairsWithinSearchRadius: number;
    minimumDistanceAngstrom: number | null;
  }[];
}
interface ConstructionDerivationAccount {
  attemptId: string;
  lipidCounts: SpeciesCountAccount[];
  cellAngstrom: number[];
  waterCount: number;
  sodiumCount: number;
  chlorideCount: number;
  proteinNetChargeElementary: number;
  intendedNaClMolar: number;
  estimatedNaClMolar: number;
  estimatedAqueousVolumeAngstromCubed: number;
  approximations: string[];
  limitations: string[];
}
interface ConstructedSystemAccount {
  subjectId: string;
  attemptId: string;
  atomCount: number;
  achievedComposition: SpeciesCountAccount[];
  cellAngstrom: number[];
  waterCount: number;
  sodiumCount: number;
  chlorideCount: number;
  conditionsTreatment: string;
  localState: LocalStateAccount | null;
}

interface AttemptAccount {
  attemptId: string;
  status: string;
  stageKind: StageKind | null;
  progress: number | null;
  message: string;
  studyRevisionId: string | null;
  policyId: string | null;
  policyVersion: string | null;
  currentStageId: string | null;
  derivation: ConstructionDerivationAccount | null;
  constructed: ConstructedSystemAccount | null;
  stopRequested: boolean;
}

interface StageObservationAccount {
  stageId: string;
  attemptId: string;
  kind: StageKind;
  measurements: MeasuredValueAccount[];
  evidence: ScientificEvidence[];
  termination: string;
  providerVersion: string;
  observedAt: string;
  localState: LocalStateAccount | null;
  proteinGeometry: ProteinGeometryObservations | null;
}

interface ExportAccount {
  stageId: string;
  assessmentId: string;
  status: 'verified' | 'failed';
  reason: string | null;
  sha256: string | null;
  byteLength: number | null;
}

export interface PreparationAssessmentResult {
  id: string;
  stageId: string;
  qualification: string;
  reason: string;
  evidence: ScientificEvidence[];
  findings: ScientificFinding[];
  limitations: string[];
  currentlyApplicable: boolean;
}

interface StageAccount {
  stageId: string;
  attemptId: string;
  studyRevisionId: string;
  kind: StageKind;
  sourceStageId?: string | null;
  status: string;
  assessment: PreparationAssessmentResult | null;
  summary: string;
  observation: StageObservationAccount | null;
  constructed?: ConstructedSystemAccount | null;
  export: ExportAccount | null;
}

interface AvailableAction { kind: ActorActionKind; subjectId: string | null; enabled: boolean; reason: string | null; }
interface WorkspaceNotice { id: string; severity: string; message: string; subjectId: string | null; }

export interface WorkspaceState {
  revision: number;
  study: StudyAccount | null;
  sourceCandidates: SourceCandidateAccount[];
  sourceModels: SourceModelObservation[];
  sourcePrediction: PredictionEvidenceAccount | null;
  availableLipids: LipidCatalogueAccount[];
  protein: ProteinAccount | null;
  preparationReview?: ProteinPreparationReviewAccount | null;
  preparationPlan?: PreparationPlanAccount | null;
  proteinTask?: ProteinTaskAccount | null;
  placementTask?: PlacementTaskAccount | null;
  placementMethods?: PlacementMethodAccount[];
  membrane: MembraneAccount | null;
  placement: PlacementAccount | null;
  attempt: AttemptAccount | null;
  stages: StageAccount[];
  inspection: InspectionAccount | null;
  actions: AvailableAction[];
  notices: WorkspaceNotice[];
}

interface EditableFraction { speciesId: string; percent: string; manual: boolean; }
const blankFraction = (): EditableFraction => ({ speciesId: '', percent: '100', manual: false });

function action(state: WorkspaceState, kind: ActorActionKind, subjectId: string | null = null): AvailableAction | undefined {
  return state.actions.find(item => item.kind === kind && item.subjectId === subjectId);
}

function sameResidue(first: ResidueAddress, second: ResidueAddress): boolean {
  return first.model === second.model && first.chain === second.chain && first.copyId === second.copyId &&
    first.residue === second.residue && first.insertionCode === second.insertionCode;
}

function readable(value: string | null | undefined): string {
  if (!value) return 'Not established';
  if (value === 'readyForMinimization') return 'Ready for minimization';
  return value.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/[-_]/g, ' ');
}

function decisionSiteLabel(decision: PreparationDecisionAccount, models: SourceModelObservation[]): string {
  const residue = decision.residue;
  const observedName = models.find(model => model.index === residue.model)?.residues.find(item =>
    item.address.chain === residue.chain && item.address.residue === residue.residue &&
    item.address.insertionCode === residue.insertionCode)?.name;
  const name: Record<string, string> = { HIS: 'Histidine', ASP: 'Aspartate', GLU: 'Glutamate',
    CYS: 'Cysteine', LYS: 'Lysine', SER: 'Serine', THR: 'Threonine', TYR: 'Tyrosine' };
  const model = models.find(item => item.index === residue.model);
  const modelScope = models.length > 1 && model ? ` · ${sourceModelLabel(model, models.length)}` : '';
  const copyScope = residue.copyId && residue.copyId !== residue.chain ? ` · Copy ${residue.copyId}` : '';
  return `${name[observedName ?? ''] ?? observedName ?? 'Residue'} ${residue.residue}${residue.insertionCode} · Chain ${residue.chain}${copyScope}${modelScope}`;
}

function decisionKindLabel(kind: PreparationChangeKind): string {
  return kind === 'residueState' ? 'Chemical-state choice' :
    kind === 'alternateLocation' ? 'Alternate-location choice' :
      kind === 'heavyAtom' ? 'Missing-atom repair' : 'Possible disulfide';
}

function optionMeaning(variant: string): string | null {
  switch (variant.trim().toUpperCase()) {
    case 'HID': return 'Neutral histidine with its ring proton on ND1.';
    case 'HIE': return 'Neutral histidine with its ring proton on NE2.';
    case 'HIP': return 'Positively charged histidine with ring protons on ND1 and NE2.';
    case 'ASP': return 'Negatively charged, deprotonated aspartate.';
    case 'ASH': return 'Neutral, protonated aspartate.';
    default: return null;
  }
}

function PredictionSummary({ prediction, label }: { prediction: PredictionEvidenceAccount | null | undefined; label: string }) {
  if (!prediction) return null;
  const numeric = prediction.localConfidence.filter(item => item.pLddt !== null);
  return <div className="hint-box">
    <strong>{label} · predicted-structure evidence</strong><br />
    <span className="tabular">Record {prediction.recordId}</span><br />
    Local confidence available for {numeric.length} of {prediction.localConfidence.length} mapped residues; PAE {readable(prediction.paeStanding)}
    {prediction.paeAxisResidueCount !== null && <> · {prediction.paeAxisResidueCount} PAE axes</>}
    {prediction.paeReason && <><br />{prediction.paeReason}</>}
    {prediction.limitations.length > 0 && <><br />{prediction.limitations.join('; ')}</>}
    <p className="help-text">These values describe prediction uncertainty; they do not by themselves establish that the selected protein or its placement is suitable.</p>
  </div>;
}

function GeometrySummary({ geometry, label }: { geometry: ProteinGeometryObservations | null | undefined; label: string }) {
  if (!geometry) return null;
  return <div className="hint-box">
    <strong>{label} · geometry measurements</strong>
    <p className="help-text">{label.toLowerCase().includes('source') ?
      'These source measurements alone do not establish structural suitability.' :
      'The measurements are observations; the protein outcome states the applicable assessment.'}</p>
    <details><summary>Measurement details</summary>
      {geometry.kinds.map(item => <div key={item.kind} className="help-text">
        <strong>{item.kind === 'covalentBond' ? 'Measured bond lengths' : item.kind === 'chainContinuity' ?
          'Distances between connected residues' : item.kind === 'nonbondedDistance' ? 'Nonbonded distances' : readable(item.kind)}</strong> · {readable(item.standing)}
        <br />{item.measuredCount} of {item.eligibleCount} applicable distances measured
        {item.minimumDistanceAngstrom !== null && <> · shortest {item.minimumDistanceAngstrom.toFixed(2)} Å</>}
        {item.maximumDistanceAngstrom !== null && <> · longest {item.maximumDistanceAngstrom.toFixed(2)} Å</>}
        {item.unavailableReason && <> · {item.unavailableReason}</>}
      </div>)}
      {geometry.limitations.length > 0 && <p className="help-text">{geometry.limitations.join('; ')}</p>}
    </details>
  </div>;
}

function PlacementPredictionSummary({ prediction }: { prediction: PredictionRegionSummaryObservations | null | undefined }) {
  if (!prediction) return null;
  return <div className="hint-box">
    <strong>Predicted relative-position evidence · {readable(prediction.standing)}</strong><br />
    <span className="tabular">Record {prediction.recordId}</span>
    {prediction.reason && <><br />{prediction.reason}</>}
    {([
      ['First contact region aligned on second', prediction.firstAlignedOnSecond],
      ['Second contact region aligned on first', prediction.secondAlignedOnFirst],
    ] as const).map(([label, summary]) => summary && <div className="help-text" key={label}>
      {label}: {summary.validPairCount}/{summary.possiblePairCount} residue pairs
      {summary.meanAngstrom !== null && <> · mean predicted aligned error {summary.meanAngstrom.toFixed(2)} Å</>}
    </div>)}
    <p className="help-text">Both directions are reported separately; this prediction is not evidence of membrane insertion on its own.</p>
  </div>;
}

function ActionButton({ state, kind, children, onClick, busy, variant = '', subjectId = null }: {
  state: WorkspaceState;
  kind: ActorActionKind;
  children: ReactNode;
  onClick: () => void;
  busy: boolean;
  variant?: string;
  subjectId?: string | null;
}) {
  const available = action(state, kind, subjectId);
  return <button
    type="button"
    className={`button ${variant}`}
    disabled={busy || available?.enabled !== true}
    title={available?.reason ?? (!available ? 'Unavailable for the current subject.' : undefined)}
    onClick={onClick}
  >{children}</button>;
}

function ActionFeedback({ kind, pending, notice }: {
  kind: PendingAction;
  pending: PendingAction | null;
  notice: ActionNotice | null;
}) {
  if (pending !== kind && notice?.kind !== kind) return null;
  const running = pending === kind;
  return <p className={`action-feedback ${running ? 'pending' : notice?.tone ?? ''}`}
    role={running || notice?.tone !== 'error' ? 'status' : 'alert'}>
    {running && <span className="activity-spinner" aria-hidden="true" />}
    {running ? activityText[kind] ?? 'Working…' : notice?.message}
  </p>;
}

function outcomeMessage(kind: ActorCommandKind, updated: WorkspaceState, previous: WorkspaceState): ActionNotice | null {
  const success = (message: string): ActionNotice => ({ kind, tone: 'success', message });
  switch (kind) {
    case 'searchSource': {
      if (updated.sourceCandidates.length > 0)
        return success(`${updated.sourceCandidates.length} matching structural source${updated.sourceCandidates.length === 1 ? '' : 's'} found. Select one to inspect.`);
      const newWarnings = updated.notices.filter(item => !previous.notices.some(before => before.id === item.id) &&
        item.severity.toLowerCase() === 'warning');
      return newWarnings.length > 0
        ? { kind, tone: 'warning', message: 'A structural search service is unavailable. Try an exact database reference or upload a structure.' }
        : { kind, tone: 'warning', message: 'No matching structures found. Try another name or accession, or enter an exact reference.' };
    }
    case 'selectSource': return success('Structure loaded. Review its coordinate models, chains, and partners.');
    case 'selectProteinModel': return success(updated.protein?.status === 'assessed'
      ? 'Protein preparation assessed for this revision. Review its evidence.'
      : 'Protein assessment returned findings. Review any proposed changes and their evidence.');
    case 'approvePreparationChange': return updated.preparationReview?.preparationStanding === 'failed'
      ? { kind, tone: 'warning', message: updated.preparationReview.preparationMessage ?? 'Protein preparation did not establish an assessed result. The choices remain recorded.' }
      : updated.preparationReview?.preparationStanding === 'blocked'
        ? { kind, tone: 'warning', message: updated.preparationReview.preparationMessage ?? 'The current selection has a preparation blocker.' }
      : success(updated.preparationReview?.preparationStanding === 'assessed'
        ? 'Choice recorded. Protein preparation and assessment completed.'
        : 'Choice recorded. Continue with the remaining review decisions.');
    case 'retryProteinPreparation': return updated.preparationReview?.preparationStanding === 'assessed'
      ? success('Protein preparation and assessment completed from the recorded choices.')
      : { kind, tone: 'warning', message: updated.preparationReview?.preparationMessage ?? 'Protein preparation did not establish an assessed result.' };
    case 'authorizePreparationPlan': return updated.proteinTask?.standing === 'assessed'
      ? success('The authorized plan was applied and the resulting protein passed its preparation checks.')
      : { kind, tone: 'warning', message: updated.proteinTask?.message ?? 'Preparation did not establish a checked protein.' };
    case 'overridePreparationPlanChoice': return updated.preparationPlan?.standing === 'ready'
      ? success('The updated plan was jointly checked. Review it before preparation.')
      : { kind, tone: 'warning', message: updated.preparationPlan?.message ?? 'The updated plan is not ready.' };
    case 'retryPreparationPlan': return updated.preparationPlan?.standing === 'ready'
      ? success('Recommendations are ready for review or authorization.')
      : { kind, tone: 'warning', message: updated.preparationPlan?.message ?? 'Recommendations remain unavailable.' };
    case 'startProteinPreparation': return updated.proteinTask?.standing === 'assessed'
      ? success('The reviewed choices were applied and the resulting protein passed its preparation checks.')
      : { kind, tone: 'warning', message: updated.proteinTask?.message ?? 'Manual preparation did not establish a checked protein.' };
    case 'proposeMembrane': return success('Membrane proposal ready. Review the intended leaflets before adoption.');
    case 'adoptMembrane': return success('Membrane model assessed for this revision. Review its support and limits.');
    case 'proposePlacement':
    case 'revisePlacement': return success('Placement assessed. Review its support and oriented structure before adoption.');
    case 'adoptPlacement': return success('Supported placement adopted for this study revision.');
    case 'startPreparation': return success('Construction request accepted. Follow the identified attempt below.');
    case 'continueMinimization': return success('Minimization request accepted. Follow the identified attempt below.');
    case 'stopAttempt': return updated.attempt?.status.toLowerCase() === 'stopped'
      ? success('The attempt has stopped. No unfinished result was promoted to a completed stage.')
      : { kind, tone: 'warning', message: 'Stop requested. The attempt remains unfinished until its actual outcome is reported.' };
    case 'requestEquilibration': return success('Optional equilibration request accepted. Follow its identified attempt.');
    default: return null;
  }
}

function FractionEditor({ label, rows, onChange, options }: {
  label: string;
  rows: EditableFraction[];
  onChange: (rows: EditableFraction[]) => void;
  options: LipidCatalogueAccount[];
}) {
  function update(index: number, change: Partial<EditableFraction>) {
    onChange(rows.map((row, current) => current === index ? { ...row, ...change } : row));
  }
  const sum = rows.reduce((total, row) => total + (Number(row.percent) || 0), 0);
  return <div>
    <span className="field-label">{label} leaflet · <span className="tabular">{sum.toFixed(1)}%</span></span>
    {rows.map((row, index) => {
      const selectedLipid = row.manual ? undefined : options.find(lipid => lipid.speciesId === row.speciesId);
      return <div key={index}>
      <div className="fraction-row">
        <select className="select-input" aria-label={`${label} leaflet lipid ${index + 1}`} value={row.manual ? '__other__' : row.speciesId} onChange={event => {
          const value = event.target.value;
          update(index, { speciesId: value === '__other__' ? '' : value, manual: value === '__other__' });
        }}>
          <option value="">Choose lipid</option>
          {options.map(lipid => <option value={lipid.speciesId} key={lipid.speciesId}>{lipid.displayName} ({lipid.speciesId})</option>)}
          <option value="__other__">Other exact species ID…</option>
        </select>
        <input className="number-input tabular" type="number" min="0" max="100" step="0.1" aria-label={`${label} leaflet percentage ${index + 1}`} value={row.percent} onChange={event => update(index, { percent: event.target.value })} />
        <button className="button compact" type="button" aria-label={`Remove ${label} leaflet lipid ${index + 1}`} disabled={rows.length === 1} onClick={() => onChange(rows.filter((_, current) => current !== index))}>×</button>
      </div>
      {row.manual && <>
        <input className="text-input tabular" aria-label={`${label} leaflet exact species ID ${index + 1}`}
          value={row.speciesId} onChange={event => update(index, { speciesId: event.target.value })}
          placeholder="Exact chemical species ID" autoCapitalize="off" spellCheck={false} />
        <p className="help-text">This exact choice will be assessed as entered. An unsupported species is reported; no catalogue species is substituted.</p>
      </>}
      {selectedLipid && selectedLipid.limitations.length > 0 && <details className="fraction-species-basis">
        <summary>{selectedLipid.speciesId} basis and limits · {selectedLipid.limitations.length} notes</summary>
        {selectedLipid.limitations.map(limit => <p className="help-text" key={limit}>{limit}</p>)}
      </details>}
    </div>;
    })}
    <button className="button compact" type="button" onClick={() => onChange([...rows, { speciesId: '', percent: '0', manual: false }])}>Add lipid</button>
  </div>;
}

function validFractions(rows: EditableFraction[]): boolean {
  if (rows.length === 0 || rows.some(row => row.percent.trim() === '' ||
      !Number.isFinite(Number(row.percent)) || Number(row.percent) < 0)) return false;
  const present = rows.filter(row => Number(row.percent) > 0);
  if (present.length === 0 || present.some(row => !row.speciesId.trim())) return false;
  if (new Set(present.map(row => row.speciesId.trim())).size !== present.length) return false;
  return Math.abs(rows.reduce((sum, row) => sum + Number(row.percent), 0) - 100) < 0.01;
}

function membraneInputIssue(upper: EditableFraction[], lower: EditableFraction[]): string | null {
  for (const [label, rows] of [['Upper', upper], ['Lower', lower]] as const) {
    if (rows.some(row => !row.percent.trim() || !Number.isFinite(Number(row.percent)) || Number(row.percent) < 0))
      return `Enter valid nonnegative percentages for the ${label.toLowerCase()} leaflet.`;
    const positive = rows.filter(row => Number(row.percent) > 0);
    if (!positive.length || positive.some(row => !row.speciesId.trim()))
      return `Choose an identified lipid for the ${label.toLowerCase()} leaflet.`;
    if (new Set(positive.map(row => row.speciesId.trim())).size !== positive.length)
      return `List each ${label.toLowerCase()} leaflet lipid only once.`;
    const total = rows.reduce((sum, row) => sum + Number(row.percent), 0);
    if (Math.abs(total - 100) >= 0.01)
      return `The ${label.toLowerCase()} leaflet totals ${total.toFixed(1)}%; adjust it to 100%.`;
  }
  return null;
}

function sameMembraneFractions(rows: EditableFraction[], proposed: LipidFraction[]): boolean {
  const actual = toFractions(rows).sort((a, b) => a.speciesId.localeCompare(b.speciesId));
  const expected = proposed.filter(item => item.fraction > 0)
    .sort((a, b) => a.speciesId.localeCompare(b.speciesId));
  return actual.length === expected.length && actual.every((item, index) =>
    item.speciesId === expected[index].speciesId && Math.abs(item.fraction - expected[index].fraction) < 0.0001);
}

function fractionSummary(fractions: LipidFraction[]): string {
  return fractions.filter(item => item.fraction > 0)
    .map(item => `${item.speciesId} ${(item.fraction * 100).toFixed(1)}%`).join(' · ') || 'none';
}

interface Guidance { message: string; target?: WorkArea; }
function actionGuidance(state: WorkspaceState, kind: ActorActionKind, subjectId: string | null = null): Guidance | null {
  if (action(state, kind, subjectId)?.enabled !== false) return null;
  const proteinReady = state.protein?.status === 'assessed';
  const membraneReady = state.membrane?.status === 'assessed';
  if (kind === 'proposePlacement' || kind === 'startPreparation') {
    if (!proteinReady) return state.protein?.status === 'declined'
      ? { message: 'The proposed protein change was declined. Review the protein findings or choose another model.', target: 'protein' }
      : state.protein
        ? { message: 'Protein preparation is not established. Review its findings and resolve any required change.', target: 'protein' }
        : { message: 'Prepare and assess a protein before continuing.', target: 'protein' };
    if (!membraneReady) return state.membrane?.status === 'notEstablished'
      ? { message: 'This membrane composition is unsupported. Review its reason or choose another composition.', target: 'membrane' }
      : state.membrane?.status === 'proposed'
        ? { message: 'Review and adopt the proposed membrane before continuing.', target: 'membrane' }
        : { message: 'Choose and assess a membrane composition first.', target: 'membrane' };
    if (kind === 'proposePlacement') {
      if (state.placementTask?.standing === 'obtaining' || state.placementTask?.standing === 'assessing')
        return { message: 'The current position is being calculated or checked. You can review other work while it runs.' };
      return { message: action(state, kind)?.reason ?? 'The current pair cannot be positioned yet.' };
    }
    if (!state.placement) return { message: 'Assess the protein–membrane placement first.', target: 'placement' };
    if (state.placement.status !== 'supported')
      return { message: 'This exact position has not passed its technical checks. Review the issue and try another position.', target: 'placement' };
    if (state.study?.adoptedPlacementProposalId !== state.placement.proposalId)
      return { message: 'Review and adopt the supported placement before constructing the system.', target: 'placement' };
    if (state.attempt?.status.toLowerCase() === 'readyforminimization')
      return { message: 'A constructed candidate is ready. Review its actual counts and continue minimization below.' };
    if (state.attempt && ['pending', 'running'].includes(state.attempt.status.toLowerCase()))
      return { message: 'The identified preparation attempt is already in progress.' };
    if (state.notices.some(item => /OpenMM|native lipid patch|Python installation/i.test(item.message)))
      return { message: 'The native construction tool or its identified assets are unavailable. Check the local installation.' };
    return { message: 'Construction is unavailable for this protein and membrane combination. No qualified construction policy and exact assets apply.' };
  }
  if (kind === 'adoptMembrane')
    return state.membrane?.status === 'notEstablished'
      ? { message: 'This membrane model could not be established. Review the reason and revise its composition.', target: 'membrane' }
      : state.membrane?.status === 'assessed'
        ? { message: 'This membrane model is already assessed for the current study revision.' }
        : { message: 'Propose a complete membrane model first.', target: 'membrane' };
  if (kind === 'adoptPlacement') {
    if (state.placement?.status === 'unsupported')
      return { message: 'This placement is unsupported. Review its evidence and revise the proposal if justified.', target: 'placement' };
    if (!state.placement || state.placement.status !== 'supported')
      return { message: 'Placement support is not established. Review the proposal and its evidence.', target: 'placement' };
    if (state.study?.adoptedPlacementProposalId === state.placement.proposalId)
      return { message: 'This placement is already adopted for the current study revision.' };
    return { message: 'The current placement needs its exact assessment evidence and intact oriented coordinates before adoption.', target: 'placement' };
  }
  if (kind === 'requestEquilibration')
    return { message: 'No validated optional equilibration procedure applies to this completed stage. It remains available for inspection and export.' };
  if (kind === 'exportStage')
    return { message: action(state, kind, subjectId)?.reason ?? 'Select a completed stage with a current assessment before exporting.', target: 'results' };
  return { message: action(state, kind, subjectId)?.reason ?? 'This action is unavailable for the current subject.' };
}

function Prerequisite({ state, kind, subjectId = null, onOpen }: {
  state: WorkspaceState;
  kind: ActorActionKind;
  subjectId?: string | null;
  onOpen: (area: WorkArea) => void;
}) {
  const guidance = actionGuidance(state, kind, subjectId);
  if (!guidance) return null;
  const technical = action(state, kind, subjectId)?.reason;
  return <div className="prerequisite" role="note">
    <span>{guidance.message}</span>
    {guidance.target && <button className="text-link" type="button" onClick={() => onOpen(guidance.target!)}>Open {workAreas.find(area => area.id === guidance.target)?.label}</button>}
    {technical && technical !== guidance.message && <details><summary>Technical detail</summary><p>{technical}</p></details>}
  </div>;
}

function toFractions(rows: EditableFraction[]): LipidFraction[] {
  return rows.filter(row => row.speciesId.trim() && Number(row.percent) > 0)
    .map(row => ({ speciesId: row.speciesId.trim(), fraction: Number(row.percent) / 100 }));
}

function selectedSourceModel(state: WorkspaceState, choice: string): SourceModelObservation | undefined {
  if (choice === '') return state.sourceModels.length === 1 ? state.sourceModels[0] : undefined;
  return state.sourceModels.find(model => model.index === Number(choice));
}

function sourceModelLabel(model: SourceModelObservation, count: number): string {
  const sourceId = model.sourceModelId?.trim();
  if (count === 1) return sourceId ? `Only coordinate model · source model ${sourceId}` : 'Only coordinate model';
  return sourceId ? `Source model ${sourceId}` : `Coordinate model ${model.index + 1} of ${count}`;
}

function partnerLabel(partner: SourcePartnerObservation): string {
  const name = partner.displayName?.trim() || partner.label;
  const component = name === partner.label ? partner.label : `${name} (${partner.label})`;
  const address = partner.chain ? ` · Chain ${partner.chain}` : '';
  const residue = partner.residue != null ? ` · Residue ${partner.residue}${partner.insertionCode ?? ''}` : '';
  return `${component}${address}${residue}`;
}

function partnersForAssembly(model: SourceModelObservation, assembly: SourceModelObservation['assemblies'][number] | undefined,
  deposited: boolean, selectedProteinChains: Set<string>, proteinSourceChains: Set<string>): SourcePartnerObservation[] {
  if (!deposited && !assembly) return [];
  return model.partners.filter(partner =>
    (deposited || !partner.chain || assembly?.chainCopies.some(copy => copy.sourceChain === partner.chain)) &&
    (!partner.chain || !proteinSourceChains.has(partner.chain) || selectedProteinChains.has(partner.chain)));
}

function selectedAssembly(model: SourceModelObservation | undefined, name: string) {
  return model?.assemblies.find(item => item.name === name);
}

function areaContainsSubject(area: WorkArea, account: WorkspaceState, subjectId: string): boolean {
  if (area === 'protein') return subjectId === account.study?.selectedSourceId ||
    subjectId === account.protein?.subjectId || subjectId === account.protein?.candidateId ||
    account.protein?.changes.some(item => item.id === subjectId) === true;
  if (area === 'membrane') return subjectId === account.membrane?.modelId;
  if (area === 'placement') return subjectId === account.placement?.proposalId ||
    subjectId === account.protein?.subjectId;
  if (area === 'preparation') return subjectId === account.attempt?.constructed?.subjectId ||
    account.stages.some(item => item.stageId === subjectId && item.attemptId === account.attempt?.attemptId);
  return account.stages.some(item => item.stageId === subjectId);
}

export function ProteinInMembraneWorkspace() {
  const [state, setState] = useState<WorkspaceState | null>(null);
  const [sourceExplorerOpen, setSourceExplorerOpen] = useState(true);
  const latestState = useRef<WorkspaceState | null>(null);
  const [communication, setCommunication] = useState<string | null>(null);
  const [streamInterrupted, setStreamInterrupted] = useState(false);
  const [feedback, setFeedback] = useState<string | null>(null);
  const [pendingAction, setPendingAction] = useState<PendingAction | null>(null);
  const [actionNotice, setActionNotice] = useState<ActionNotice | null>(null);
  const [inspectionFailure, setInspectionFailure] = useState<{ subjectId: string; reason: string } | null>(null);
  const [sourceIntent, setSourceIntent] = useState<SourceIntent | null>(null);
  const [structureLoad, setStructureLoad] = useState<StructureLoadStatus | null>(null);
  const [structureReload, setStructureReload] = useState(0);
  const [sourceChainColors, setSourceChainColors] = useState<{ subjectId: string; structureUrl: string;
    colors: Record<string, string> } | null>(null);
  const [sourcePreview, setSourcePreview] = useState<SourcePreviewState | null>(null);
  const [sourcePreviewRetry, setSourcePreviewRetry] = useState(0);
  const [chainFocus, setChainFocus] = useState<{ chainId: string; serial: number } | null>(null);
  const [busy, setBusy] = useState(false);
  const [exportBusy, setExportBusy] = useState(false);
  const [exportTargetId, setExportTargetId] = useState<string | null>(null);
  const [localExportFault, setLocalExportFault] = useState<{ stageId: string; reason: string } | null>(null);
  const [localExportNotice, setLocalExportNotice] = useState<{ stageId: string; message: string } | null>(null);
  const [workflowOpen, setWorkflowOpen] = useState(true);
  const [activeArea, setActiveArea] = useState<WorkArea>('protein');
  const activeAreaRef = useRef<WorkArea>('protein');
  const areaSubject = useRef<Record<WorkArea, string | null>>({ protein: null, membrane: null,
    placement: null, preparation: null, results: null });
  const previousSourceIdentity = useRef<string | null | undefined>(undefined);
  const [viewTarget, setViewTarget] = useState<string | null>(null);
  const [viewRestoreFailure, setViewRestoreFailure] = useState<{ subjectId: string; reason: string } | null>(null);
  const viewRestoreRunning = useRef(false);
  const [attemptReviewRequested, setAttemptReviewRequested] = useState(false);
  const [selectedDecisionId, setSelectedDecisionId] = useState<string | null>(null);
  const [selectedOptionId, setSelectedOptionId] = useState<string | null>(null);
  const [decisionCommandTarget, setDecisionCommandTarget] = useState<string | null>(null);
  const [decisionFilter, setDecisionFilter] = useState<'all' | 'remaining' | 'confirmed'>('all');
  const [reviewOtherOpen, setReviewOtherOpen] = useState(false);
  const [selectedStageId, setSelectedStageId] = useState<string | null>(null);
  const busyRef = useRef(false);
  const sourceRequestSerial = useRef(0);
  const newestSourceRequest = useRef<SourceRequest | null>(null);
  const sourceRequestRunning = useRef(false);
  const focusRequestGeneration = useRef(0);
  const refreshIndex = useRef(0);
  const initialAccountPresented = useRef(false);
  const automaticallyInspectedConstruction = useRef<string | null>(null);
  const manualPositionAttempt = useRef<string | null>(null);
  const manualDraftInitialized = useRef(false);

  const [query, setQuery] = useState('');
  const [exactSourceKind, setExactSourceKind] = useState<'rcsb' | 'alphafold'>('rcsb');
  const [exactIdentifier, setExactIdentifier] = useState('');
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [uploadProvenance, setUploadProvenance] = useState<UploadOriginKind>('unknown');
  const [uploadProvenanceNote, setUploadProvenanceNote] = useState('');
  const uploadInput = useRef<HTMLInputElement>(null);
  const [modelChoice, setModelChoice] = useState('');
  const [assemblyChoice, setAssemblyChoice] = useState('');
  const [selectionDetailsOpen, setSelectionDetailsOpen] = useState(false);
  const [reviewConfirmedOpen, setReviewConfirmedOpen] = useState(false);
  const [reviewPlanOpen, setReviewPlanOpen] = useState(false);
  const [chainChoices, setChainChoices] = useState<string[]>([]);
  const [partnerChoices, setPartnerChoices] = useState<Record<string, PartnerSelection>>({});
  const [altlocChoices, setAltlocChoices] = useState<Record<string, number>>({});
  const [upper, setUpper] = useState<EditableFraction[]>([blankFraction()]);
  const [lower, setLower] = useState<EditableFraction[]>([blankFraction()]);
  const [topologyKind, setTopologyKind] = useState<ProteinTopologyKind | ''>('');
  const [orientationRoute, setOrientationRoute] = useState<'manual' | 'ppm' | 'opm'>('manual');
  const [startingPosition, setStartingPosition] = useState<PlacementStartingPosition>('center');
  const [offsetX, setOffsetX] = useState('0');
  const [offsetY, setOffsetY] = useState('0');
  const [offsetZ, setOffsetZ] = useState('0');
  const [rotationX, setRotationX] = useState('0');
  const [rotationY, setRotationY] = useState('0');
  const [rotationZ, setRotationZ] = useState('0');
  const [manualDraftTouched, setManualDraftTouched] = useState(false);
  const [physicalSide, setPhysicalSide] = useState<PlacementPhysicalSide | ''>('');
  const [ppmNterminalSide, setPpmNterminalSide] = useState<PpmNterminalSide | ''>('');
  const [depthShift, setDepthShift] = useState('0');
  const [tiltX, setTiltX] = useState('0');
  const [tiltY, setTiltY] = useState('0');
  const [rotationNormal, setRotationNormal] = useState('0');

  async function refresh() {
    const index = ++refreshIndex.current;
    try {
      const response = await fetch('/api/state', { cache: 'no-store' });
      if (!response.ok) throw new Error(`The local workspace could not be read (${response.status}).`);
      const account = await response.json() as WorkspaceState;
      if (index === refreshIndex.current) {
        if (!latestState.current || latestState.current.revision <= account.revision) {
          latestState.current = account;
          // An earlier source request may complete while a newer selection is
          // queued. Keep the preceding account until the newest request resolves;
          // the inspection pane masks its molecular view during retrieval.
          if (!sourceRequestRunning.current || !newestSourceRequest.current)
            setState(account);
          if (!initialAccountPresented.current) {
            initialAccountPresented.current = true;
            const inspected = account.inspection?.subjectId;
            const inferredArea: WorkArea = account.stages.some(stage => stage.stageId === inspected) ? 'results' :
              account.attempt && (!inspected || inspected === account.attempt.constructed?.subjectId) ? 'preparation' :
              inspected && inspected === account.placement?.proposalId ? 'placement' :
              inspected && inspected === account.membrane?.modelId ? 'membrane' : 'protein';
            activeAreaRef.current = inferredArea;
            setActiveArea(inferredArea);
            areaSubject.current[inferredArea] = inspected ?? null;
          }
        }
        setCommunication(null);
      }
    } catch (cause) {
      if (index === refreshIndex.current) setCommunication(cause instanceof Error ? cause.message : 'The local workspace is unavailable.');
    }
  }

  useEffect(() => {
    void refresh();
    const events = new EventSource('/api/events');
    events.onopen = () => { setStreamInterrupted(false); void refresh(); };
    events.onmessage = () => { setStreamInterrupted(false); void refresh(); };
    events.onerror = () => setStreamInterrupted(true);
    return () => events.close();
  }, []);

  useEffect(() => {
    if (!state?.membrane) return;
    const toEditable = (fractions: LipidFraction[]) => fractions.map(item => ({
      speciesId: item.speciesId,
      percent: String(item.fraction * 100),
      manual: !state.availableLipids.some(candidate => candidate.speciesId === item.speciesId),
    }));
    setUpper(toEditable(state.membrane.upper));
    setLower(toEditable(state.membrane.lower));
  }, [state?.membrane?.modelId]);

  useEffect(() => {
    const account = state?.study;
    const currentSource = account?.selectedSourceId ?? null;
    if (previousSourceIdentity.current !== undefined && previousSourceIdentity.current !== currentSource) {
      areaSubject.current.protein = currentSource;
      setViewTarget(null);
      setViewRestoreFailure(null);
    }
    previousSourceIdentity.current = currentSource;
    const adopted = state?.sourceModels.find(item => item.index === account?.modelIndex);
    const draft = adopted ?? (state?.sourceModels.length === 1 ? state.sourceModels[0] : undefined);
    setModelChoice(draft ? String(draft.index) : '');
    setAssemblyChoice(adopted ? account?.biologicalAssemblyId ?? 'deposited' :
      draft?.assemblies.length ? '' : draft ? 'deposited' : '');
    if (adopted && account) {
      const copies = account.biologicalAssemblyId
        ? adopted.assemblies.find(item => item.name === account.biologicalAssemblyId)?.chainCopies ?? []
        : adopted.chains.map(item => ({ sourceChain: item.name, copyId: item.name }));
      setChainChoices(copies.filter(copy => account.chainIds.includes(copy.copyId))
        .map(copy => `${copy.sourceChain}:${copy.copyId}`));
      setPartnerChoices(Object.fromEntries(account.partners.map(item => [item.sourceId, item])));
      setAltlocChoices(Object.fromEntries(account.alternateLocations.flatMap(choice => {
        const source = adopted.residues.find(item => item.address.model === choice.residue.model &&
          item.address.chain === choice.residue.chain && item.address.residue === choice.residue.residue &&
          item.address.insertionCode === choice.residue.insertionCode);
        const index = source?.alternateLocations.indexOf(choice.altloc) ?? -1;
        return index < 0 ? [] : [[`${choice.residue.model}:${choice.residue.chain}:${choice.residue.residue}:${choice.residue.insertionCode}:${choice.residue.copyId}`, index]];
      })));
    } else {
      setChainChoices([]);
      setPartnerChoices({});
      setAltlocChoices({});
    }
    setSelectionDetailsOpen(false);
    setReviewConfirmedOpen(false);
    setReviewPlanOpen(false);
    setSourceExplorerOpen(!state?.study?.selectedSourceId);
  }, [state?.study?.selectedSourceId]);

  const previewSourceId = state?.study?.selectedSourceId ?? null;
  const previewModel = state?.sourceModels.find(item => String(item.index) === modelChoice);
  const previewAssemblyId = assemblyChoice && assemblyChoice !== 'deposited' ? assemblyChoice : null;
  const previewRequired = !!previewSourceId && state?.inspection?.subjectId === previewSourceId &&
    !!previewModel && (state.sourceModels.length > 1 || previewAssemblyId !== null);
  const previewKey = previewRequired ? `${previewSourceId}|${previewModel!.index}|${previewAssemblyId ?? ''}` : null;

  useEffect(() => {
    if (!previewKey || !previewSourceId || !previewModel) return;
    const controller = new AbortController();
    setSourceChainColors(null);
    setSourcePreview({ key: previewKey, phase: 'loading' });
    const query = new URLSearchParams({ sourceId: previewSourceId, modelIndex: String(previewModel.index) });
    if (previewAssemblyId) query.set('assemblyId', previewAssemblyId);
    void fetch(`/api/source-preview?${query}`, { signal: controller.signal, cache: 'no-store' })
      .then(async response => {
        const body = await response.json();
        if (!response.ok) throw new Error(body?.reason ?? `Source preview failed (${response.status}).`);
        return body as SourcePreviewAccount;
      }).then(account => {
        if (!controller.signal.aborted && account.sourceId === previewSourceId &&
            account.modelIndex === previewModel.index && account.assemblyId === previewAssemblyId)
          setSourcePreview({ key: previewKey, phase: 'ready', account });
      }).catch(cause => {
        if (!controller.signal.aborted) setSourcePreview({ key: previewKey, phase: 'failed',
          reason: cause instanceof Error ? cause.message : 'This exact model could not be displayed.' });
      });
    return () => controller.abort();
  }, [previewKey, sourcePreviewRetry]);

  useEffect(() => {
    setSelectedDecisionId(null);
    setSelectedOptionId(null);
    setDecisionCommandTarget(null);
    setDecisionFilter('all');
    setReviewOtherOpen(false);
  }, [state?.preparationReview?.studyRevisionId]);

  useEffect(() => {
    if (['preparing', 'assessed', 'failed', 'unavailable'].includes(state?.proteinTask?.standing ?? ''))
      setReviewConfirmedOpen(false);
  }, [state?.proteinTask?.standing, state?.proteinTask?.studyRevisionId]);

  useEffect(() => {
    setDepthShift('0');
    setTiltX('0');
    setTiltY('0');
    setRotationNormal('0');
  }, [state?.placement?.proposalId]);

  useEffect(() => {
    manualDraftInitialized.current = false;
    manualPositionAttempt.current = null;
    setManualDraftTouched(false);
    setStartingPosition('center');
    setOffsetX('0'); setOffsetY('0'); setOffsetZ('0');
    setRotationX('0'); setRotationY('0'); setRotationZ('0');
  }, [state?.protein?.subjectId, state?.membrane?.modelId]);

  useEffect(() => {
    const transform = state?.placement?.transform;
    if (!transform || manualDraftInitialized.current || manualDraftTouched) return;
    manualDraftInitialized.current = true;
    setStartingPosition(transform.startingPosition);
    setOffsetX(String(transform.offsetXAngstrom));
    setOffsetY(String(transform.offsetYAngstrom));
    setOffsetZ(String(transform.offsetZAngstrom));
    setRotationX(String(transform.rotationXDegrees));
    setRotationY(String(transform.rotationYDegrees));
    setRotationZ(String(transform.rotationZDegrees));
  }, [state?.placement?.proposalId, state?.protein?.subjectId,
    state?.membrane?.modelId, manualDraftTouched]);

  useEffect(() => {
    if (!state || activeArea !== 'placement' || orientationRoute !== 'manual' || busy || busyRef.current ||
        action(state, 'proposePlacement')?.enabled !== true ||
        state.placement && !state.placement.transform && !manualDraftTouched) return;
    const raw = [offsetX, offsetY, offsetZ, rotationX, rotationY, rotationZ];
    if (raw.some(value => value.trim() === '' || !Number.isFinite(Number(value)))) return;
    const values = raw.map(Number);
    const current = state.placement?.transform;
    if (current && startingPosition === current.startingPosition &&
        values.every((value, index) => value === [current.offsetXAngstrom,
          current.offsetYAngstrom, current.offsetZAngstrom, current.rotationXDegrees,
          current.rotationYDegrees, current.rotationZDegrees][index])) return;
    const key = [state.study?.id, state.protein?.subjectId, state.membrane?.modelId,
      startingPosition, ...values].join('|');
    if (manualPositionAttempt.current === key) return;
    const timer = window.setTimeout(() => {
      if (busyRef.current) return;
      manualPositionAttempt.current = key;
      void command('proposePlacement', { orientationRoute: 'manual', startingPosition,
        offsetXAngstrom: values[0], offsetYAngstrom: values[1], offsetZAngstrom: values[2],
        rotationXDegrees: values[3], rotationYDegrees: values[4], rotationZDegrees: values[5] });
    }, 450);
    return () => window.clearTimeout(timer);
  }, [state?.revision, activeArea, orientationRoute, busy, manualDraftTouched,
    startingPosition, offsetX, offsetY, offsetZ, rotationX, rotationY, rotationZ]);

  useEffect(() => {
    if (!state || activeArea === 'membrane') return;
    const desired = areaSubject.current[activeArea];
    if (!desired) return;
    if (!areaContainsSubject(activeArea, state, desired)) {
      areaSubject.current[activeArea] = null;
      setViewTarget(null);
      setViewRestoreFailure(null);
      return;
    }
    if (state.inspection?.subjectId === desired) {
      if (viewTarget) setViewTarget(null);
      if (viewRestoreFailure) setViewRestoreFailure(null);
    } else if (!busy && viewTarget !== desired && viewRestoreFailure?.subjectId !== desired) {
      setViewTarget(desired);
    }
  }, [activeArea, state?.inspection?.subjectId, state?.revision, busy, viewTarget, viewRestoreFailure]);

  useEffect(() => {
    if (!state || !viewTarget || busy || viewRestoreRunning.current || activeArea === 'membrane' ||
        state.inspection?.subjectId === viewTarget || !areaContainsSubject(activeArea, state, viewTarget)) return;
    const target = viewTarget;
    viewRestoreRunning.current = true;
    void command('selectInspectionSubject', { subjectId: target }, reason =>
      setViewRestoreFailure({ subjectId: target, reason })).then(ok => {
        viewRestoreRunning.current = false;
        if (!ok) setViewTarget(current => current === target ? null : current);
      });
  }, [viewTarget, busy, activeArea, state?.inspection?.subjectId]);

  async function command(kind: ActorCommandKind, data: object,
                         onFailure?: (reason: string) => void): Promise<boolean> {
    const current = latestState.current;
    if (!current || busyRef.current) return false;
    const originatingArea = activeAreaRef.current;
    const inspectionOnly = kind === 'selectInspectionSubject' || kind === 'setInspectionFocus';
    busyRef.current = true;
    setBusy(true);
    setPendingAction(kind);
    if (!inspectionOnly) setActionNotice(null);
    setFeedback(null);
    setInspectionFailure(null);
    try {
      const response = await fetch('/api/commands', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ kind, data, expectedRevision: current.revision }),
      });
      const result: unknown = await response.json();
      if (!response.ok) {
        const reason = typeof result === 'object' && result !== null && 'reason' in result
          ? String(result.reason) : `The request was refused (${response.status}).`;
        if (inspectionOnly) {
          setFeedback(reason);
          if (kind === 'selectInspectionSubject' && 'subjectId' in data && typeof data.subjectId === 'string')
            setInspectionFailure({ subjectId: data.subjectId, reason });
        }
        else setActionNotice({ kind, tone: 'error', message: response.status === 409
          ? `The study changed while this request was pending. Review the current account and retry. ${reason}`
          : response.status >= 500 ? `The service is unavailable. Retry when it is restored. ${reason}`
            : reason });
        onFailure?.(reason);
        await refresh();
        return false;
      }
      const updated = result as WorkspaceState;
      if (updated.inspection && activeAreaRef.current === originatingArea &&
          kind !== 'setInspectionFocus' && kind !== 'exportStage')
        areaSubject.current[originatingArea] = updated.inspection.subjectId;
      if (!latestState.current || latestState.current.revision <= updated.revision) {
        latestState.current = updated;
        setState(updated);
      }
      setCommunication(null);
      setFeedback(null);
      if (!inspectionOnly && kind !== 'exportStage') setActionNotice(outcomeMessage(kind, updated, current));
      return true;
    } catch (cause) {
      const reason = cause instanceof Error ? cause.message : 'The local request did not complete.';
      if (inspectionOnly) {
        setFeedback(reason);
        if (kind === 'selectInspectionSubject' && 'subjectId' in data && typeof data.subjectId === 'string')
          setInspectionFailure({ subjectId: data.subjectId, reason });
      }
      else setActionNotice({ kind, tone: 'error', message: `The request did not complete. Check the local service and retry. ${reason}` });
      onFailure?.(reason);
      await refresh();
      return false;
    } finally {
      busyRef.current = false;
      setBusy(false);
      setPendingAction(null);
    }
  }

  function chooseSource(data: object, key: string, label: string) {
    if (busyRef.current && !sourceRequestRunning.current) return;
    if (sourceRequestRunning.current && newestSourceRequest.current?.key === key) return;
    const request = { serial: ++sourceRequestSerial.current, key, label, data };
    focusRequestGeneration.current += 1;
    newestSourceRequest.current = request;
    setSourceIntent({ ...request, phase: 'retrieving' });
    setActionNotice(null);
    setFeedback(null);
    if (!sourceRequestRunning.current) void drainSourceRequests();
  }

  async function drainSourceRequests() {
    if (sourceRequestRunning.current) return;
    sourceRequestRunning.current = true;
    busyRef.current = true;
    setBusy(true);
    setPendingAction('selectSource');
    try {
      while (newestSourceRequest.current) {
        const request = newestSourceRequest.current;
        const current = latestState.current;
        if (!current) break;
        try {
          const response = await fetch('/api/commands', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ kind: 'selectSource', data: request.data, expectedRevision: current.revision }),
          });
          const result: unknown = await response.json();
          if (!response.ok) {
            if (response.status === 409) {
              await refresh();
              if (newestSourceRequest.current?.serial === request.serial &&
                  latestState.current?.revision !== current.revision) continue;
            }
            if (newestSourceRequest.current?.serial === request.serial) {
              const reason = typeof result === 'object' && result !== null && 'reason' in result
                ? String(result.reason) : `The source could not be retrieved (${response.status}).`;
              newestSourceRequest.current = null;
              setSourceIntent({ ...request, phase: 'failed', reason });
              // A superseded earlier request may have established a source in
              // the host while this newer request failed. Show that actual
              // account under its own identity with the failed request nearby.
              if (latestState.current) setState(latestState.current);
            }
            continue;
          }
          const updated = result as WorkspaceState;
          if (!latestState.current || latestState.current.revision <= updated.revision) {
            latestState.current = updated;
            if (newestSourceRequest.current?.serial === request.serial)
              setState(updated);
          }
          if (newestSourceRequest.current?.serial === request.serial && updated.inspection)
            areaSubject.current.protein = updated.inspection.subjectId;
          setCommunication(null);
          if (newestSourceRequest.current?.serial === request.serial) {
            newestSourceRequest.current = null;
            setSourceIntent({ ...request, key: updated.study?.selectedSourceId ?? request.key,
              phase: 'rendering' });
            if (request.key.startsWith('upload:')) {
              setUploadFile(null);
              if (uploadInput.current) uploadInput.current.value = '';
            }
          }
        } catch (cause) {
          if (newestSourceRequest.current?.serial === request.serial) {
            newestSourceRequest.current = null;
            setSourceIntent({ ...request, phase: 'failed', reason: cause instanceof Error
              ? cause.message : 'The source request did not complete.' });
            await refresh();
          }
        }
      }
    } finally {
      sourceRequestRunning.current = false;
      busyRef.current = false;
      setBusy(false);
      setPendingAction(null);
      if (newestSourceRequest.current) void drainSourceRequests();
    }
  }

  function onStructureLoad(status: StructureLoadStatus) {
    const current = latestState.current?.inspection;
    const exactPreview = sourcePreview?.account;
    if (current?.subjectId !== status.subjectId ||
        current.structureUrl !== status.structureUrl &&
        (exactPreview?.sourceId !== status.subjectId || exactPreview.structureUrl !== status.structureUrl)) return;
    setStructureLoad(status);
    if (status.phase !== 'loading')
      setSourceIntent(intent => intent?.phase === 'rendering' && intent.key === status.subjectId ? null : intent);
  }

  useEffect(() => {
    const subjectId = state?.attempt?.constructed?.subjectId;
    if (!attemptReviewRequested || !subjectId || busy || busyRef.current ||
        state?.inspection?.subjectId === subjectId ||
        automaticallyInspectedConstruction.current === subjectId ||
        action(state, 'selectInspectionSubject')?.enabled !== true) return;
    automaticallyInspectedConstruction.current = subjectId;
    void command('selectInspectionSubject', { subjectId });
  }, [attemptReviewRequested, busy, state?.attempt?.constructed?.subjectId,
      state?.inspection?.subjectId, state?.revision]);

  async function inspectSubject(subjectId: string): Promise<boolean> {
    const area = activeAreaRef.current;
    const previous = areaSubject.current[area];
    areaSubject.current[area] = subjectId;
    if (await command('selectInspectionSubject', { subjectId })) {
      setAttemptReviewRequested(false);
      return true;
    }
    areaSubject.current[area] = previous;
    return false;
  }

  async function focusDecisionOption(proposalId: string) {
    const generation = ++focusRequestGeneration.current;
    const current = latestState.current?.inspection;
    if (current?.subjectId !== proposalId && !await inspectSubject(proposalId)) return;
    if (generation !== focusRequestGeneration.current) return;
    const exact = latestState.current?.inspection;
    if (exact?.subjectId !== proposalId) return;
    const annotation = exact.annotations.find(item => item.geometryFocus);
    if (!annotation) {
      setInspectionFailure({ subjectId: proposalId,
        reason: 'No mapped residue is available for local focus. The option remains reviewable from its evidence.' });
      return;
    }
    await command('setInspectionFocus', { annotationId: annotation.id });
  }

  async function confirmDecisionOption(proposalId: string, approve: boolean) {
    setDecisionCommandTarget(proposalId);
    const accepted = await command('approvePreparationChange', { proposalId, approve });
    if (accepted && activeAreaRef.current === 'protein' &&
        ['preparing', 'assessed', 'failed', 'unavailable'].includes(latestState.current?.proteinTask?.standing ?? ''))
      window.requestAnimationFrame(() => document.getElementById('protein-task-status')?.focus());
  }

  async function startPreparation() {
    if (await command('startPreparation', {})) {
      setAttemptReviewRequested(true);
      showWorkArea('preparation');
    }
  }

  async function continueMinimization(attempt: AttemptAccount) {
    if (!attempt.constructed) return;
    if (await command('continueMinimization', {
      attemptId: attempt.attemptId,
      constructedSubjectId: attempt.constructed.subjectId,
    })) {
      setAttemptReviewRequested(true);
      showWorkArea('preparation');
    }
  }

  function showWorkArea(area: WorkArea) {
    focusRequestGeneration.current += 1;
    activeAreaRef.current = area;
    setActiveArea(area);
    setViewRestoreFailure(null);
    const current = latestState.current;
    if (current && areaSubject.current[area] &&
        !areaContainsSubject(area, current, areaSubject.current[area]!))
      areaSubject.current[area] = null;
    const preferred = areaSubject.current[area] ?? (area === 'protein' ? current?.study?.selectedSourceId ?? current?.protein?.subjectId :
      area === 'placement' ? current?.placement?.proposalId ?? current?.protein?.subjectId :
      area === 'preparation' ? current?.attempt?.constructed?.subjectId ?? null :
      area === 'results' ? selectedStageId ?? current?.stages[0]?.stageId ?? null : null);
    if (preferred) areaSubject.current[area] = preferred;
    setViewTarget(area === 'membrane' || !preferred || preferred === current?.inspection?.subjectId ? null : preferred);
    setWorkflowOpen(true);
    window.requestAnimationFrame(() => {
      const rail = document.getElementById('researcher-workflow');
      if (rail) rail.scrollTop = 0;
    });
  }

  async function uploadSource() {
    if (!state || !uploadFile || busyRef.current || action(state, 'selectSource')?.enabled !== true) return;
    busyRef.current = true;
    setBusy(true);
    setPendingAction('uploadSource');
    setActionNotice(null);
    setFeedback(null);
    let token: string | null = null;
    try {
      const form = new FormData();
      form.append('file', uploadFile, uploadFile.name);
      form.append('provenance', uploadProvenance);
      if (uploadProvenanceNote.trim()) form.append('provenanceNote', uploadProvenanceNote.trim());
      const response = await fetch('/api/uploads', { method: 'POST', body: form });
      const result: unknown = await response.json();
      if (!response.ok) {
        const reason = typeof result === 'object' && result !== null && 'reason' in result
          ? String(result.reason) : `The source could not be uploaded (${response.status}).`;
        setActionNotice({ kind: 'uploadSource', tone: 'error', message: reason });
        return;
      }
      if (typeof result === 'object' && result !== null && 'uploadToken' in result && typeof result.uploadToken === 'string') {
        token = result.uploadToken;
      }
      if (!token) throw new Error('The local host did not identify the uploaded source.');
    } catch (cause) {
      setActionNotice({ kind: 'uploadSource', tone: 'error', message: cause instanceof Error
        ? cause.message : 'The source upload did not complete.' });
      return;
    } finally {
      busyRef.current = false;
      setBusy(false);
      setPendingAction(null);
    }
    chooseSource({ uploadToken: token }, `upload:${token}`, uploadFile.name);
  }

  async function exportStage(stageId: string) {
    const current = latestState.current;
    if (exportBusy || !current ||
        action(current, 'exportStage', stageId)?.enabled !== true) return;
    setExportTargetId(stageId);
    setExportBusy(true);
    setActionNotice(null);
    setLocalExportNotice(null);
    try {
      if (!await command('exportStage', { stageId }, reason => {
        setLocalExportFault({ stageId, reason });
        setActionNotice({ kind: 'exportStage', tone: 'error', message: `Export not delivered. ${reason}` });
      })) return;
      const selected = latestState.current?.stages.find(stage => stage.stageId === stageId);
      const delivery = selected?.export;
      if (delivery?.stageId !== stageId || delivery.status !== 'verified' ||
          !selected?.assessment?.currentlyApplicable ||
          delivery.assessmentId !== selected.assessment.id ||
          !delivery.sha256 || !/^[0-9a-f]{64}$/.test(delivery.sha256) ||
          typeof delivery.byteLength !== 'number' ||
          !Number.isSafeInteger(delivery.byteLength) || delivery.byteLength <= 0)
        throw new Error('The host did not identify a verified bundle for this completed stage.');

      const response = await fetch(`/api/export/${encodeURIComponent(stageId)}`, {
        cache: 'no-store', headers: { 'If-Match': `"${delivery.sha256}"` },
      });
      if (!response.ok) {
        let reason = `The verified bundle could not be transferred (${response.status}).`;
        try {
          const body: unknown = await response.json();
          if (typeof body === 'object' && body !== null && 'reason' in body &&
              typeof body.reason === 'string' && body.reason.trim()) reason = body.reason;
        } catch { /* Preserve the HTTP status when no JSON reason is available. */ }
        throw new Error(reason);
      }
      if (!response.headers.get('content-type')?.toLowerCase().includes('application/zip') ||
          response.headers.get('etag') !== `"${delivery.sha256}"` ||
          response.headers.get('x-content-sha256') !== delivery.sha256)
        throw new Error('The downloaded bundle identity was not confirmed by the local host.');
      const bytes = await response.arrayBuffer();
      if (bytes.byteLength !== delivery.byteLength)
        throw new Error('The downloaded bytes did not match the verified bundle length.');
      const digest = Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', bytes)),
        value => value.toString(16).padStart(2, '0')).join('');
      if (digest !== delivery.sha256)
        throw new Error('The downloaded bytes did not match the verified bundle digest.');

      const url = URL.createObjectURL(new Blob([bytes], { type: 'application/zip' }));
      try {
        const link = document.createElement('a');
        link.href = url;
        link.download = `protein-membrane-${stageId}.zip`;
        document.body.appendChild(link);
        link.click();
        link.remove();
      } finally {
        window.setTimeout(() => URL.revokeObjectURL(url), 60_000);
      }
      setLocalExportFault(null);
      setLocalExportNotice({ stageId, message: 'Verified bundle received; browser download requested.' });
      setActionNotice({ kind: 'exportStage', tone: 'success', message: 'Verified bundle received; browser download requested.' });
    } catch (cause) {
      const reason = cause instanceof Error ? cause.message : 'The verified bundle could not be transferred.';
      setLocalExportFault({ stageId, reason });
      setActionNotice({ kind: 'exportStage', tone: 'error', message: `Export not delivered. ${reason}` });
      await refresh();
    } finally {
      setExportBusy(false);
    }
  }

  if (!state) return <main className="workspace"><header className="workspace-header"><div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation</div></div></div></header><div className="scene-empty" role="status"><strong>Connecting to the local workspace</strong><span>{communication ?? (streamInterrupted ? 'The live connection is interrupted. Reading established local standing…' : 'Reading the current scientific account…')}</span></div></main>;

  const model = selectedSourceModel(state, modelChoice);
  const assembly = selectedAssembly(model, assemblyChoice);
  const assemblyResolved = assemblyChoice === 'deposited' || !!assembly;
  const proteinSourceChains = new Set(model?.residues.filter(residue => residue.residueKind === 'protein')
    .map(residue => residue.address.chain) ?? []);
  const availableChains = (assembly?.chainCopies ?? (assemblyChoice === 'deposited'
    ? model?.chains.map(chain => ({ sourceChain: chain.name, copyId: chain.name })) : []) ?? [])
    .filter(chain => proteinSourceChains.has(chain.sourceChain));
  const chosenChains = availableChains.filter(chain => chainChoices.includes(`${chain.sourceChain}:${chain.copyId}`));
  const selectedSourceChains = new Set(chosenChains.map(chain => chain.sourceChain));
  const relevantPartners = model ? partnersForAssembly(model, assembly, assemblyChoice === 'deposited',
    selectedSourceChains, proteinSourceChains) : [];
  const omittedPartners = model?.partners.filter(partner => !relevantPartners.some(item => item.sourceId === partner.sourceId)) ?? [];
  const allPartnersDecided = !!model && assemblyResolved && relevantPartners.every(partner => !!partnerChoices[partner.sourceId]);
  const chosenPartners = relevantPartners.map(partner => partnerChoices[partner.sourceId]).filter((item): item is PartnerSelection => !!item);
  const selectedProteinBackboneIncomplete = model?.residues.some(residue =>
    residue.residueKind === 'protein' && selectedSourceChains.has(residue.address.chain) &&
    !residue.backboneHeavyAtomsComplete) ?? false;
  const ambiguousResidues = model?.residues.filter(residue =>
    residue.residueKind === 'protein' && selectedSourceChains.has(residue.address.chain) &&
    residue.alternateLocations.length > 0) ?? [];
  const residueKey = (address: ResidueAddress) => `${address.model}:${address.chain}:${address.residue}:${address.insertionCode}:${address.copyId}`;
  const allAltlocsChosen = ambiguousResidues.every(residue => Number.isInteger(altlocChoices[residueKey(residue.address)]));
  const chosenAltlocs = ambiguousResidues.map(residue => ({ residue: residue.address, altloc: residue.alternateLocations[altlocChoices[residueKey(residue.address)]] }));
  const membraneReady = validFractions(upper) && validFractions(lower);
  const membraneDraftMatchesProposal = !!state.membrane &&
    sameMembraneFractions(upper, state.membrane.upper) && sameMembraneFractions(lower, state.membrane.lower);
  const correctionValues = [depthShift, tiltX, tiltY, rotationNormal].map(Number);
  const correctionReady = state.placement !== null && correctionValues.every(Number.isFinite)
    && correctionValues.some(value => value !== 0);
  const manualFields = [offsetX, offsetY, offsetZ, rotationX, rotationY, rotationZ];
  const manualValues = manualFields.map(Number);
  const manualInputReady = manualFields.every(value => value.trim() !== '') && manualValues.every(Number.isFinite);
  const currentTransform = state.placement?.transform;
  const manualDraftMatchesProposal = !!currentTransform && startingPosition === currentTransform.startingPosition &&
    manualValues.every((value, index) => Math.abs(value - [currentTransform.offsetXAngstrom,
      currentTransform.offsetYAngstrom, currentTransform.offsetZAngstrom,
      currentTransform.rotationXDegrees, currentTransform.rotationYDegrees,
      currentTransform.rotationZDegrees][index]) <= 1e-9);
  const selectedStage = state.stages.find(stage => stage.stageId === selectedStageId)
    ?? state.stages.find(stage => stage.stageId === state.inspection?.subjectId);
  const selectedExportFault = selectedStage ? (localExportFault?.stageId === selectedStage.stageId
    ? localExportFault.reason : selectedStage.export?.status === 'failed'
      ? selectedStage.export.reason ?? 'The verified bundle could not be delivered.' : null) : null;
  const selectedExportNotice = selectedStage && localExportNotice?.stageId === selectedStage.stageId
    ? localExportNotice.message : null;
  const selectedSourceStage = selectedStage?.kind === 'Equilibration' && selectedStage.sourceStageId
    ? state.stages.find(stage => stage.stageId === selectedStage.sourceStageId &&
      stage.attemptId === selectedStage.attemptId && stage.kind === 'Minimization') : undefined;
  const completedStageSummary = selectedStage
    ? (['Minimization', 'Equilibration'] as const)
      .filter(kind => state.stages.some(stage => stage.attemptId === selectedStage.attemptId &&
        stage.kind === kind && stage.status.toLowerCase() === 'completed'))
      .map(kind => kind === 'Minimization' ? 'Minimized' : 'Equilibrated').join(' · ')
    : state.stages.length === 0 ? 'None' : String(state.stages.length);
  const preparationReview = state.preparationReview;
  const preparationPlan = state.preparationPlan ?? null;
  const partialPlan = preparationPlan?.standing === 'partial';
  const primaryDecisions = preparationReview?.decisions.filter(item => !partialPlan ||
    ((item.kind === 'alternateLocation' || item.kind === 'disulfide') && item.standing !== 'confirmed')) ?? [];
  const otherPartialDecisions = partialPlan ? preparationReview?.decisions.filter(item =>
    !primaryDecisions.some(primary => primary.id === item.id)) ?? [] : [];
  const proteinTask = state.proteinTask ?? null;
  const placementTask = state.placementTask ?? null;
  const proteinTaskOutcome = !!proteinTask && ['preparing', 'assessed', 'failed', 'blocked', 'unavailable'].includes(proteinTask.standing);
  const selectedDecision = preparationReview?.decisions.find(item => item.id === selectedDecisionId)
    ?? preparationReview?.decisions.find(item => item.options.some(option => option.proposalId === state.inspection?.subjectId))
    ?? primaryDecisions.find(item => item.standing === 'pending')
    ?? primaryDecisions[0] ?? preparationReview?.decisions[0];
  const selectedOption = selectedDecision?.options.find(option => option.proposalId === selectedOptionId)
    ?? selectedDecision?.options.find(option => option.proposalId === selectedDecision.chosenProposalId)
    ?? selectedDecision?.options.find(option => preparationPlan?.choices.some(choice =>
      sameResidue(choice.residue, selectedDecision.residue) && choice.variant === option.proposedChange))
    ?? (selectedDecision?.options.length === 1 ? selectedDecision.options[0] : undefined);
  const reviewPanelVisible = activeArea === 'protein' && !!preparationReview && !!selectedDecision &&
    !(preparationPlan?.standing === 'ready' && !reviewPlanOpen) &&
    (reviewConfirmedOpen || !proteinTaskOutcome ||
      (proteinTask?.standing === 'blocked' && preparationReview.remainingCount > 0));
  const proteinOutcomeVisible = activeArea === 'protein' && !!proteinTask && proteinTaskOutcome && !reviewPanelVisible;
  const reviewChoices = primaryDecisions.filter(item =>
    item.kind === 'residueState' || item.kind === 'alternateLocation') ?? [];
  const reviewOther = (partialPlan ? [] : primaryDecisions).filter(item =>
    item.kind === 'heavyAtom' || item.kind === 'disulfide') ?? [];
  const visibleReviewChoices = partialPlan ? primaryDecisions : reviewChoices.filter(item =>
    decisionFilter === 'all' || (decisionFilter === 'confirmed' ? item.standing === 'confirmed' :
      item.standing !== 'confirmed'));
  const nextUnresolved = selectedDecision && primaryDecisions.length
    ? primaryDecisions.slice(primaryDecisions.findIndex(item => item.id === selectedDecision.id) + 1)
      .concat(primaryDecisions.slice(0, primaryDecisions.findIndex(item => item.id === selectedDecision.id)))
      .find(item => item.standing === 'pending')
    : primaryDecisions.find(item => item.standing === 'pending');
  const chooseDecision = (id: string) => {
    focusRequestGeneration.current += 1;
    setSelectedDecisionId(id);
    setSelectedOptionId(null);
    if (reviewOther.some(item => item.id === id)) setReviewOtherOpen(true);
  };
  const goToNextUnresolved = () => {
    if (!nextUnresolved) return;
    setDecisionFilter('all');
    chooseDecision(nextUnresolved.id);
  };
  const sourceId = state.study?.selectedSourceId ?? null;
  const sourceIsViewed = !!sourceId && state.inspection?.subjectId === sourceId;
  const readyPreview = previewRequired && sourcePreview?.key === previewKey && sourcePreview.phase === 'ready'
    ? sourcePreview.account : null;
  const sourcePreviewLabel = previewRequired && model ? `${sourceModelLabel(model, state.sourceModels.length)} · ${
    previewAssemblyId ? `biological assembly ${previewAssemblyId}` : 'deposited coordinates'}` : null;
  const previewPending = previewRequired && !readyPreview;
  const previewFailure = sourcePreview?.key === previewKey && sourcePreview.phase === 'failed'
    ? sourcePreview.reason : null;
  const displayedState: WorkspaceState = readyPreview && sourceIsViewed && state.inspection
    ? { ...state, inspection: { ...state.inspection, structureUrl: readyPreview.structureUrl } } : state;
  const sourceLoad = sourceIsViewed && structureLoad?.subjectId === sourceId &&
    structureLoad.structureUrl === displayedState.inspection?.structureUrl ? structureLoad : null;
  const sourceLabel = sourceIntent?.label ?? state.study?.selectedSourceLabel ??
    state.sourceCandidates.find(candidate => candidate.id === sourceId)?.label ?? sourceId;
  const sourceStatus = sourceIntent?.phase === 'retrieving' ? 'Retrieving source…'
    : sourceIntent?.phase === 'failed' ? 'Retrieval failed'
      : previewFailure ? 'Selected model view unavailable'
        : previewPending ? 'Loading selected model…'
      : sourceIntent?.phase === 'rendering' || sourceLoad?.phase === 'loading' ? 'Loading structure…'
        : sourceLoad?.phase === 'failed' ? 'Visualization unavailable'
          : sourceLoad?.phase === 'displayed' ? 'Structure displayed'
            : sourceIsViewed ? 'Loading structure…' : 'Selected source is not currently displayed';
  const selectedMembrane = state.membrane?.modelId === state.inspection?.subjectId ? state.membrane : null;
  const selectedPlacement = state.placement?.proposalId === state.inspection?.subjectId ? state.placement : null;
  const selectedPlacementAdopted = selectedPlacement !== null &&
    state.study?.adoptedPlacementProposalId === selectedPlacement.proposalId;
  const activeAttempt = !!state.attempt && ['pending', 'running'].includes(state.attempt.status.toLowerCase());
  const minimizationRunning = activeAttempt && state.attempt?.stageKind === 'Minimization';
  const awaitingMinimization = state.attempt?.status.toLowerCase() === 'readyforminimization';
  const selectedConstructed = !!state.attempt?.constructed &&
    state.inspection?.subjectId === state.attempt.constructed.subjectId;
  const reviewAttempt = !selectedStage && !!state.attempt &&
    (attemptReviewRequested || selectedConstructed || (activeAttempt && !state.inspection));
  const executionReview = reviewAttempt || !!selectedStage;
  const currentAttemptStage = state.stages.find(stage => stage.attemptId === state.attempt?.attemptId);
  const pendingAreaSubject = activeArea === 'membrane' ? null : areaSubject.current[activeArea] &&
    areaContainsSubject(activeArea, state, areaSubject.current[activeArea]!) &&
    areaSubject.current[activeArea] !== state.inspection?.subjectId ? areaSubject.current[activeArea] : null;
  return <main className={`workspace ${workflowOpen ? 'workflow-open' : 'workflow-closed'}${selectedMembrane ? ' membrane-review' : ''}${selectedPlacement && !executionReview ? ' placement-review' : ''}${executionReview ? ' execution-review' : ''}`}>
    <header className="workspace-header">
      <div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation · one planar bilayer</div></div></div>
      <nav className="workspace-context-tabs" aria-label="Research work areas">
        {workAreas.map(area => <button key={area.id} type="button"
          className={activeArea === area.id ? 'active' : ''}
          aria-current={activeArea === area.id ? 'page' : undefined}
          title={area.description} onClick={() => showWorkArea(area.id)}>{area.label}</button>)}
      </nav>
      <div className="header-context">
        <button className="button workflow-toggle" type="button" aria-expanded={workflowOpen} aria-controls="researcher-workflow"
          onClick={() => setWorkflowOpen(value => !value)}>{workflowOpen ? 'Collapse inputs' : 'Show inputs'}</button>
        <span className="context-chip" title={state.study?.id ?? 'No study established'}>Study <strong className="id">revision {state.study?.number ?? 'not established'}</strong></span>
        <span className="context-chip" title={activeArea === 'membrane' ? 'Composition diagram' : pendingAreaSubject ?? state.inspection?.subjectId ?? 'No inspected subject'}>Viewing <strong>{pendingAreaSubject
          ? `Restoring ${workAreas.find(item => item.id === activeArea)?.label} view`
          : activeArea === 'results' && state.stages.length === 0 ? 'No completed stage'
          : activeArea === 'preparation' && !state.attempt && state.stages.length === 0 ? 'No system attempt'
          : activeArea === 'placement' && !state.protein && !state.placement ? 'No prepared protein'
          : activeArea === 'membrane'
          ? !membraneDraftMatchesProposal || !state.membrane ? 'Membrane composition draft' : state.membrane.status === 'proposed' ? 'Membrane proposal' : 'Chosen membrane'
          : viewerSubjectHeading(state.inspection, state.stages.find(stage => stage.stageId === state.inspection?.subjectId), state.study?.selectedSourceLabel)}</strong></span>
        {state.study?.modelIndex !== null && state.study?.modelIndex !== undefined && <span className="context-chip compact-model-context" title={`${state.study.biologicalAssemblyId ? `Biological assembly ${state.study.biologicalAssemblyId}` : 'Deposited coordinates'}; chains ${state.study.chainIds.join(', ')}`}>
          <strong>{state.study.biologicalAssemblyId ? `Assembly ${state.study.biologicalAssemblyId}` : 'Deposited coordinates'} · {state.study.chainIds.join(', ')}</strong>
        </span>}
        {state.attempt && <span className="context-chip" title={`Attempt ${state.attempt.attemptId}`}>Attempt <strong>{readable(state.attempt.status)}</strong></span>}
        {selectedStage && <span className="context-chip" title={`Stage ${selectedStage.stageId}`}>Stage <strong>{readable(selectedStage.kind)}</strong></span>}
      </div>
    </header>

    <div className="workspace-body">
      <aside id="researcher-workflow" className="rail" aria-label="Researcher choices and available actions" hidden={!workflowOpen}>
        <h2 className="panel-heading">{workAreas.find(area => area.id === activeArea)?.label}</h2>
        <p className="work-area-description">{workAreas.find(area => area.id === activeArea)?.description}</p>
        <div className="work-context" aria-label="Current work status">
          <small className="work-summary">Protein {state.protein ? readable(state.protein.status) : 'not prepared'} · Membrane {state.membrane ? readable(state.membrane.status) : 'not chosen'} · Stage {selectedStage ? readable(selectedStage.kind) : 'none selected'}</small>
        </div>
        {pendingAction && <div className="work-pending" role="status"><span className="activity-spinner" aria-hidden="true" />{activityText[pendingAction] ?? 'Working…'}</div>}
        {feedback && <div className="notice warning" role="alert">{feedback}</div>}
        {state.notices.length > 0 && <details className="notice-list" aria-label="Scientific and workflow notices"
          open={state.notices.some(notice => notice.severity.toLowerCase() === 'error')}>
          <summary>Study notices ({state.notices.length})</summary>
          <div className="notice-list-items">{state.notices.map(notice => <div key={notice.id} className={`notice ${notice.severity.toLowerCase()}`}><strong>{readable(notice.severity)}</strong> · {notice.message}{notice.subjectId && <small className="tabular">Subject {notice.subjectId}</small>}</div>)}</div>
        </details>}

        <section className="rail-section" id="protein-workflow" hidden={activeArea !== 'protein'}>
          {proteinTask && (proteinTaskOutcome || proteinTask.standing === 'assessing') &&
            <section className={`protein-task-account ${proteinTask.standing}`} aria-label="Protein task outcome">
              <h3 id="protein-task-status" tabIndex={-1} role="status">{proteinTask.standing === 'assessed' ? 'Protein prepared' :
                proteinTask.standing === 'preparing' ? 'Preparing and checking protein…' :
                proteinTask.standing === 'assessing' ? 'Assessing selected protein…' :
                proteinTask.standing === 'failed' ? 'Protein could not be prepared' :
                proteinTask.standing === 'unavailable' ? 'Preparation status unavailable' :
                preparationReview?.remainingCount === 0 ? 'Choices complete — preparation blocked' : 'Protein preparation blocked'}</h3>
              {!proteinOutcomeVisible && <p>{proteinTask.message}</p>}
              <p className="task-subject">{state.study?.selectedSourceLabel ?? proteinTask.sourceId} · {proteinTask.chains.map(chain =>
                chain.sourceChain === chain.copyId ? `Chain ${chain.sourceChain}` : `Chain ${chain.sourceChain}, copy ${chain.copyId}`).join(' and ')}</p>
              {proteinTask.standing === 'preparing' && <p className="review-activity" role="status"><span className="activity-spinner" aria-hidden="true" />The confirmed choices remain recorded while the result is checked.</p>}
              {proteinTask.standing === 'assessed' && <div className="button-row"><button className="button primary" type="button"
                onClick={() => showWorkArea(state.membrane?.status === 'assessed' ? 'placement' : 'membrane')}>
                {state.membrane?.status === 'assessed' ? 'Review placement' : 'Choose membrane'}</button></div>}
              {proteinTask.standing === 'failed' && proteinTask.retryAvailable && <div className="button-row"><ActionButton
                state={state} kind="retryProteinPreparation" busy={busy} onClick={() => void command('retryProteinPreparation', {})}>Retry protein preparation</ActionButton></div>}
              {proteinTask.standing === 'unavailable' && <button className="button" type="button" onClick={() => void refresh()}>Refresh preparation status</button>}
              {preparationReview && preparationReview.confirmedCount > 0 && <button className="button compact" type="button" onClick={() => {
                setReviewConfirmedOpen(value => !value);
                if (!reviewConfirmedOpen) {
                  setDecisionFilter('confirmed');
                  setSelectedDecisionId(preparationReview.decisions.find(item => item.standing === 'confirmed')?.id ?? null);
                }
              }}>{reviewConfirmedOpen ? 'Back to protein result' : `Review confirmed choices (${preparationReview.confirmedCount})`}</button>}
              <details className="task-technical"><summary>Exact task identity</summary><p className="tabular">Study revision {proteinTask.studyRevisionId}<br />Intended protein {proteinTask.intendedProteinId}{proteinTask.preparedProteinId && <><br />Prepared result {proteinTask.preparedProteinId}</>}</p></details>
            </section>}
          {preparationPlan && ['calculating', 'ready', 'partial', 'failed'].includes(preparationPlan.standing) &&
            <section className="preparation-plan-card" aria-label="Preparation plan">
              <h3>{preparationPlan.standing === 'calculating' ? 'Finding preparation suggestions…' :
                preparationPlan.standing === 'ready' ? preparationPlan.hasOverrides ? 'Updated plan ready' : 'Recommendations ready' :
                  preparationPlan.standing === 'partial' ? 'Choices need your input' : 'Recommendations unavailable'}</h3>
              {preparationPlan.standing === 'calculating' ?
                <p className="review-activity" role="status"><span className="activity-spinner" aria-hidden="true" />{preparationPlan.message}</p> :
                <p>{preparationPlan.message}</p>}
              {preparationPlan.standing === 'ready' && <>
                <p><strong>{preparationPlan.stateChoiceCount} state choices</strong> · {preparationPlan.heavyAtomCount} atom repairs
                  {preparationPlan.removedSourceHydrogenCount > 0 && <> · {preparationPlan.removedSourceHydrogenCount} source hydrogens to replace</>}</p>
                <p className="help-text">Starting-state suggestions at pH {preparationPlan.nominalPh}; you can change any choice. These changes have not been applied.</p>
                <p className="review-consequence">Preparing will apply the exact checked plan and assess the resulting protein.</p>
                <div className="button-row">
                  <ActionButton state={state} kind="authorizePreparationPlan" subjectId={preparationPlan.planId}
                    busy={busy || !preparationPlan.authorizationAvailable} variant="primary"
                    onClick={() => void command('authorizePreparationPlan', { planSha256: preparationPlan.planSha256 })}>
                    {preparationPlan.hasOverrides ? 'Prepare with these choices' : 'Prepare with recommendations'}
                  </ActionButton>
                  {(preparationPlan.stateChoiceCount > 0 || preparationPlan.heavyAtomCount > 0) &&
                    <button className="button" type="button" onClick={() => {
                      setDecisionFilter('all'); setReviewPlanOpen(value => !value);
                    }}>
                      {reviewPlanOpen ? 'Hide choice review' : 'Review or change choices'}
                    </button>}
                </div>
                <details className="task-technical"><summary>Method and exact plan</summary>
                  <p>{preparationPlan.method} · plan {preparationPlan.planSha256}</p>
                  <p>Suggested states are starting choices. Joint parameter matching is a technical check, not a claim of optimal chemistry.</p>
                </details>
              </>}
              {preparationPlan.standing === 'failed' && <div className="button-row">
                <ActionButton state={state} kind="retryPreparationPlan" busy={busy}
                  onClick={() => void command('retryPreparationPlan', {})}>Retry suggestions</ActionButton>
                <button className="button" type="button" onClick={() => setReviewPlanOpen(true)}>Review choices manually</button>
              </div>}
              <ActionFeedback kind="authorizePreparationPlan" pending={pendingAction} notice={actionNotice} />
              <ActionFeedback kind="retryPreparationPlan" pending={pendingAction} notice={actionNotice} />
            </section>}
          {proteinTask?.standing === 'ready' && preparationPlan?.standing !== 'ready' &&
            <section className="preparation-plan-card" aria-label="Manual protein preparation">
              <h3>Ready to prepare</h3>
              <p>The selected protein has no remaining manual decisions. Preparing will apply its recorded choices and check the actual result.</p>
              <ActionButton state={state} kind="startProteinPreparation" busy={busy} variant="primary"
                onClick={() => void command('startProteinPreparation', {})}>Prepare selected protein</ActionButton>
              <ActionFeedback kind="startProteinPreparation" pending={pendingAction} notice={actionNotice} />
            </section>}
          {preparationReview && (!proteinTaskOutcome || reviewConfirmedOpen ||
            (proteinTask?.standing === 'blocked' && preparationReview.remainingCount > 0)) &&
            (preparationPlan?.standing !== 'ready' || reviewPlanOpen) &&
            (preparationPlan?.standing !== 'failed' || reviewPlanOpen) &&
            (preparationReview.decisions.length > 0 ||
            ['preparing', 'failed', 'blocked'].includes(preparationReview.preparationStanding)) &&
            <section className="review-queue" aria-label="Protein preparation decisions">
            {preparationReview.decisions.length > 0 && <div className="review-progress" role="status"><strong>{partialPlan ?
              `${primaryDecisions.length} structural exception${primaryDecisions.length === 1 ? '' : 's'}` :
              `${reviewChoices.length} site choice${reviewChoices.length === 1 ? '' : 's'}`}</strong>
              <span>{partialPlan ? `${preparationReview.decisions.filter(item =>
                item.kind === 'alternateLocation' || item.kind === 'disulfide').length - primaryDecisions.length} of ${preparationReview.decisions.filter(item =>
                item.kind === 'alternateLocation' || item.kind === 'disulfide').length} structural exceptions resolved` :
                preparationPlan?.standing === 'ready' ?
                `${preparationPlan.stateChoiceCount} method-resolved states · ${preparationPlan.heavyAtomCount} proposed repairs` :
                `${preparationReview.confirmedCount} of ${preparationReview.decisions.length} obligations settled · ${preparationReview.remainingCount} remain`}</span>
              <small>{partialPlan ? 'Resolve these exact conformer or bond choices. Other decisions remain available below.' :
                preparationPlan?.standing === 'ready' ? 'Inspect or override any suggested state. Preparation starts only after whole-plan authorization.' :
                `${reviewOther.length > 0 ? `${reviewOther.length} separate repair or bond decisions. ` : ''}Alternatives at one site count as one decision.`}</small></div>}
            {preparationReview.preparationStanding === 'preparing' && <p className="review-activity" role="status"><span className="activity-spinner" aria-hidden="true" />{preparationReview.preparationMessage}</p>}
            {preparationReview.preparationStanding === 'failed' && !proteinTask && <div className="review-blocker" role="alert"><strong>Preparation did not establish a protein</strong><p>{preparationReview.preparationMessage}</p></div>}
            {preparationReview.blockers.length > 0 && <div className="review-blocker" role="alert"><strong>Preparation blocked</strong>
              {preparationReview.blockers.map((blocker, index) => <p key={`${index}:${blocker}`}>{blocker}</p>)}</div>}
            {reviewOther.length > 0 && <p className="help-text">Repairs and possible bonds have their own disposition. Hydrogens are added during preparation for confirmed states.</p>}
            {(state.protein?.findings.length ?? 0) > 0 && <details className="review-other"><summary>Findings to read, not approve ({state.protein!.findings.length})</summary>
              {state.protein!.findings.map(finding => <p key={finding.id}><strong>{finding.meaning}</strong><br />{finding.consequence}</p>)}</details>}
            {preparationReview.decisions.length === 0 && <p className="help-text">No site or repair choice is pending for this selection.</p>}
            {preparationReview.decisions.length > 0 && <>
            {preparationPlan?.standing !== 'ready' && !partialPlan && <div className="review-filters" role="group" aria-label="Filter site choices">
              {(['all', 'remaining', 'confirmed'] as const).map(filter => <button key={filter} type="button"
                className={`button compact${decisionFilter === filter ? ' active' : ''}`}
                aria-pressed={decisionFilter === filter} onClick={() => setDecisionFilter(filter)}>{filter === 'all' ? 'All sites' : filter === 'remaining' ? 'Open sites' : 'Reviewed sites'} ({filter === 'all' ? reviewChoices.length :
                  filter === 'remaining' ? reviewChoices.filter(item => item.standing !== 'confirmed').length :
                    reviewChoices.filter(item => item.standing === 'confirmed').length})</button>)}
            </div>}
            <div className="review-site-list" aria-label="Site decision queue">
              {visibleReviewChoices.map(item => <button className={`review-site${selectedDecision?.id === item.id ? ' selected' : ''}`}
                type="button" key={item.id} aria-current={selectedDecision?.id === item.id ? 'true' : undefined}
                onClick={() => chooseDecision(item.id)}>
                <span className={`review-site-mark ${item.standing}`} aria-hidden="true">{item.standing === 'confirmed' ? '✓' : item.standing === 'blocked' ? '!' : '○'}</span>
                <span><strong>{decisionSiteLabel(item, state.sourceModels)}</strong><small>{partialPlan ?
                  `${decisionKindLabel(item.kind)} · ${item.standing === 'blocked' ? 'Review reason' : 'Needs choice'}` :
                  preparationPlan?.standing === 'ready' ?
                  `${preparationPlan.choices.find(choice => sameResidue(choice.residue, item.residue))?.overridden ? 'Your choice' : 'Suggested'} · ${preparationPlan.choices.find(choice => sameResidue(choice.residue, item.residue))?.variant ?? 'review required'}` : item.standing === 'confirmed'
                  ? `Confirmed · ${item.options.find(option => option.proposalId === item.chosenProposalId)?.proposedChange ?? 'declined'}`
                  : item.standing === 'blocked' ? 'Blocked · review reason' : 'Needs choice'}</small></span>
                <span aria-hidden="true">›</span>
              </button>)}
              {visibleReviewChoices.length === 0 && <p className="help-text">No sites in this view. Choose another filter to inspect them.</p>}
            </div>
            {reviewOther.length > 0 && <details className="review-other" open={reviewOtherOpen || reviewChoices.length === 0}
              onToggle={event => setReviewOtherOpen(event.currentTarget.open)}><summary>Repairs and other findings ({reviewOther.length})</summary>
              {reviewOther.map(item => <button type="button" key={item.id} className={`review-site${selectedDecision?.id === item.id ? ' selected' : ''}`}
                onClick={() => chooseDecision(item.id)}><span className={`review-site-mark ${item.standing}`}>{item.standing === 'confirmed' ? '✓' : item.standing === 'blocked' ? '!' : '○'}</span>
                <span><strong>{decisionKindLabel(item.kind)} · {decisionSiteLabel(item, state.sourceModels)}</strong><small>{readable(item.standing)}</small></span><span aria-hidden="true">›</span></button>)}
              </details>}
            {otherPartialDecisions.length > 0 && <details className="review-other"><summary>Other choices available for review ({otherPartialDecisions.length})</summary>
              {otherPartialDecisions.map(item => <button type="button" key={item.id}
                className={`review-site${selectedDecision?.id === item.id ? ' selected' : ''}`}
                onClick={() => chooseDecision(item.id)}><span className={`review-site-mark ${item.standing}`} aria-hidden="true">{item.standing === 'confirmed' ? '✓' : item.standing === 'blocked' ? '!' : '○'}</span>
                <span><strong>{decisionKindLabel(item.kind)} · {decisionSiteLabel(item, state.sourceModels)}</strong><small>{readable(item.standing)}</small></span><span aria-hidden="true">›</span></button>)}
            </details>}
            </>}
          </section>}
          {proteinTask && <button className="button compact selection-details-toggle" type="button"
            aria-expanded={selectionDetailsOpen} aria-controls="protein-selection-inputs protein-model-inputs"
            onClick={() => setSelectionDetailsOpen(value => !value)}>{selectionDetailsOpen ? 'Hide selection details' : 'Selection details and source'}</button>}
          <div id="protein-selection-inputs" className="selection-inputs" hidden={!!proteinTask && !selectionDetailsOpen}>
              <h3 className="section-title">Structural source</h3>
          {(state.placement || state.stages.length > 0) && <p className="change-impact">Changing the chosen protein starts a new study revision. Placement must be reassessed; completed stages remain attached to their original revision.</p>}
          {(sourceId || sourceIntent) && <div className="source-context" role="status" aria-label="Selected source and loading state">
            <strong className="tabular">{sourceLabel}</strong>
            <span>{sourceStatus}</span>
            {sourceIntent?.phase === 'failed' && <><span className="source-failure">{sourceIntent.reason}</span>
              <small>View and evidence remain on {state.inspection?.subjectId ?? 'no previously displayed subject'}.</small>
              <button className="button compact" type="button" onClick={() => chooseSource(sourceIntent.data, sourceIntent.key, sourceIntent.label)}>Retry retrieval</button></>}
            {sourceLoad?.phase === 'failed' && !sourceIntent && <><span className="source-failure">{sourceLoad.reason}</span>
              <button className="button compact" type="button" onClick={() => setStructureReload(value => value + 1)}>Retry visualization</button></>}
            {sourceId && !sourceIsViewed && !sourceIntent && <button className="button compact" type="button" disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true}
              onClick={() => void inspectSubject(sourceId)}>View selected source</button>}
            {!sourceIntent && sourceLoad?.phase === 'displayed' && !model && <small>Source coordinates are displayed for inspection; no preparation model is chosen by viewing them.</small>}
          </div>}
          <details className="source-explorer" open={sourceExplorerOpen}
            onToggle={event => setSourceExplorerOpen(event.currentTarget.open)}>
            <summary>{sourceId ? 'Find or load another source' : 'Find or load a source'}</summary>
          <label className="field-label" htmlFor="source-query">Discover structural sources</label>
          <input id="source-query" className="text-input" value={query} onChange={event => setQuery(event.target.value)} placeholder="Protein name or accession" />
          <div className="button-row"><ActionButton state={state} kind="searchSource" busy={busy || !query.trim()} onClick={() => void command('searchSource', { query: query.trim() })}>Find sources</ActionButton></div>
          {!query.trim() && <p className="input-guidance">Enter a protein name or accession to search.</p>}
          <ActionFeedback kind="searchSource" pending={pendingAction} notice={actionNotice} />
          {action(state, 'searchSource')?.reason && <p className="action-reason">{action(state, 'searchSource')?.reason}</p>}
          {state.sourceCandidates.length > 0 && <div className="source-list">
            {state.sourceCandidates.map(candidate => <button type="button"
              className={`source-item ${sourceId === candidate.id ? 'selected' : ''} ${sourceIntent?.key === candidate.id ? 'requested' : ''}`}
              key={candidate.id} disabled={(busy && pendingAction !== 'selectSource') || action(state, 'selectSource')?.enabled !== true}
              onClick={() => chooseSource({ sourceId: candidate.id }, candidate.id, candidate.label)}>
              <span className="item-title">{candidate.label}</span>
              <span className="item-detail">{candidate.kind} · {candidate.provenance}</span>
              {(sourceIntent?.key === candidate.id || sourceId === candidate.id) && <span className="item-detail source-item-state">
                {sourceIntent?.key === candidate.id && sourceIntent.phase === 'retrieving' ? 'Requested · Retrieving source…'
                  : sourceIntent?.key === candidate.id && sourceIntent.phase === 'failed' ? 'Requested · Retrieval failed'
                  : sourceId === candidate.id ? sourceIntent && sourceIntent.key !== candidate.id
                    ? 'Current study source' : `Selected · ${sourceStatus}` : 'Requested'}</span>}
            </button>)}
          </div>}
          <details className="source-input-routes" open={!sourceId && state.sourceCandidates.length === 0 ? true : undefined}>
            <summary>Exact reference or upload</summary>
            <div className="source-routes-body">
              <section className="source-route" aria-labelledby="source-reference-heading">
                <h4 id="source-reference-heading" className="source-route-heading">Database reference</h4>
                <p className="help-text">Use an exact reference when search is unavailable. Loading displays its coordinates; choose the model and chains separately.</p>
                <label className="field-label" htmlFor="exact-source-kind">Source archive</label>
                <select id="exact-source-kind" className="select-input" value={exactSourceKind} onChange={event => {
                  const value = event.target.value;
                  if (value === 'rcsb' || value === 'alphafold') setExactSourceKind(value);
                }}>
                  <option value="rcsb">RCSB PDB</option>
                  <option value="alphafold">AlphaFold DB</option>
                </select>
                <label className="field-label" htmlFor="exact-identifier">{exactSourceKind === 'rcsb' ? 'PDB entry ID' : 'Exact prediction model or fragment ID'}</label>
                <input id="exact-identifier" className="text-input tabular" value={exactIdentifier} onChange={event => setExactIdentifier(event.target.value)}
                  placeholder={exactSourceKind === 'rcsb' ? 'For example, 1BL8' : 'For example, AF-Q9Y2J2-F1'} autoCapitalize="off" spellCheck={false} />
                <div className="button-row"><ActionButton state={state} kind="selectSource" busy={busy || !exactIdentifier.trim()}
                  onClick={() => chooseSource({ sourceKind: exactSourceKind, exactIdentifier: exactIdentifier.trim() },
                    `exact:${exactSourceKind}:${exactIdentifier.trim()}`, `${exactSourceKind === 'rcsb' ? 'RCSB PDB' : 'AlphaFold DB'} ${exactIdentifier.trim()}`)}>Load source</ActionButton></div>
                {!exactIdentifier.trim() && <p className="input-guidance">Enter an exact database ID to inspect this route.</p>}
              </section>
              <section className="source-route" aria-labelledby="source-upload-heading">
                <h4 id="source-upload-heading" className="source-route-heading">File upload</h4>
                <label className="field-label" htmlFor="source-upload">PDB/mmCIF file</label>
                <input id="source-upload" ref={uploadInput} className="text-input" type="file" accept=".pdb,.cif,.mmcif" onChange={event => setUploadFile(event.target.files?.[0] ?? null)} />
                <label className="field-label" htmlFor="upload-provenance">Researcher-declared structure origin</label>
                <select id="upload-provenance" className="select-input" value={uploadProvenance} onChange={event => {
                  const value = event.target.value;
                  if (value === 'unknown' || value === 'experimental' || value === 'predicted') setUploadProvenance(value);
                }}>
                  <option value="unknown">Unknown</option>
                  <option value="experimental">Experimental</option>
                  <option value="predicted">Predicted</option>
                </select>
                <label className="field-label" htmlFor="upload-provenance-note">Source or method note, if available</label>
                <input id="upload-provenance-note" className="text-input" maxLength={500} value={uploadProvenanceNote} onChange={event => setUploadProvenanceNote(event.target.value)} placeholder="Archive, experiment, or prediction method" />
                <p className="help-text">This declaration identifies the source; it does not verify its method or turn uploaded B-factors into prediction confidence. A predicted or unknown upload needs an attributable qualification route before placement can be supported.</p>
                <div className="button-row"><ActionButton state={state} kind="selectSource" busy={busy || !uploadFile} onClick={() => void uploadSource()}>Upload source</ActionButton></div>
                {uploadFile && <p className="input-guidance tabular">Selected file: {uploadFile.name}</p>}
                <ActionFeedback kind="uploadSource" pending={pendingAction} notice={actionNotice} />
                {action(state, 'selectSource')?.reason && <p className="action-reason">{action(state, 'selectSource')?.reason}</p>}
              </section>
            </div>
          </details>
          </details>
          </div>
        </section>

        <section id="protein-model-inputs" className="rail-section" hidden={activeArea !== 'protein' || (!!proteinTask && !selectionDetailsOpen)}>
              <h3 className="section-title">Protein model</h3>
          {state.study?.modelIndex !== null && state.study?.modelIndex !== undefined && <p className="help-text">The current assessed selection uses {sourceModelLabel(state.sourceModels.find(item => item.index === state.study!.modelIndex) ??
            { index: state.study.modelIndex, sourceModelId: null } as SourceModelObservation, state.sourceModels.length)}. Changes below are drafts until you assess them.</p>}
          {state.sourceModels.length === 0 ? <p className="help-text">Inspect a structural source to see its observed models, assemblies, chains and partners.</p> : <>
            <span className="field-label">Coordinate model</span>
            {state.sourceModels.length === 1 ? <div className="sole-model" role="status">{sourceModelLabel(state.sourceModels[0], 1)} — selected for this draft</div> :
              <select id="model-index" className="select-input" aria-label="Coordinate model" value={modelChoice} onChange={event => {
                const next = state.sourceModels.find(item => item.index === Number(event.target.value));
                setModelChoice(event.target.value); setAssemblyChoice(next && next.assemblies.length === 0 ? 'deposited' : '');
                setChainChoices([]); setPartnerChoices({}); setAltlocChoices({});
              }}>
                <option value="">Choose a coordinate model</option>
                {state.sourceModels.map(item => <option key={item.index} value={item.index}>{sourceModelLabel(item, state.sourceModels.length)} · {item.atomCount.toLocaleString()} atoms</option>)}
              </select>}
            {model && <>
              {state.sourceModels.length > 1 && <p className="help-text">This structure contains several sets of coordinates. Choose the model to prepare; the displayed source model is only an inspection view.</p>}
              <label className="field-label" htmlFor="assembly-choice">Assembly</label>
              <p className="help-text">The set of protein chains arranged together in the model.</p>
              <select id="assembly-choice" className="select-input" value={assemblyChoice} onChange={event => { setAssemblyChoice(event.target.value); setChainChoices([]); }}>
                {model.assemblies.length > 0 && <option value="" disabled>Choose an assembly</option>}
                <option value="deposited">Deposited coordinates</option>
                {model.assemblies.map(item => <option key={item.name} value={item.name}>Biological assembly {item.name}</option>)}
              </select>
              {assemblyChoice === 'deposited' && <p className="help-text">Chains as supplied in this file. This does not establish the biological unit.</p>}
              {assembly && <p className="help-text">Archive-described biological assembly {assembly.name}. Its protein copies include {availableChains.map(copy => copy.sourceChain === copy.copyId
                ? `chain ${copy.sourceChain}` : `chain ${copy.sourceChain}, copy ${copy.copyId}`).join(', ')}. Choose the copies to retain.</p>}
              {previewPending && sourcePreviewLabel && <p className={`help-text ${previewFailure ? 'source-failure' : ''}`} role="status">
                {previewFailure ? <>The exact {sourcePreviewLabel.toLowerCase()} view could not be loaded: {previewFailure}{' '}
                  <button className="button compact" type="button" onClick={() => setSourcePreviewRetry(value => value + 1)}>Retry model view</button></>
                  : <>Loading {sourcePreviewLabel.toLowerCase()} for inspection…</>}
              </p>}
              {model.assemblies.length === 0 && <p className="help-text">No usable biological assembly was supplied with this source. Deposited coordinates remain available for assessment.</p>}
              {assemblyResolved && <><span className="field-label">Protein chains to retain</span>
              {state.inspection?.subjectId === sourceId && <p className="help-text">Showing all source chains. The checkboxes choose preparation membership; viewer highlighting only changes inspection.</p>}
              {availableChains.map(chain => {
                const key = `${chain.sourceChain}:${chain.copyId}`;
                const modelNumber = model.sourceModelId?.trim() || String(model.index + 1);
                const visibleColor = state.inspection?.subjectId === sourceId && !previewPending &&
                  sourceChainColors?.subjectId === sourceId &&
                  sourceChainColors.structureUrl === displayedState.inspection?.structureUrl &&
                  (!previewRequired || readyPreview?.chainIds.includes(chain.copyId))
                    ? sourceChainColors.colors[chain.copyId] ?? sourceChainColors.colors[`${modelNumber}:${chain.copyId}`] : undefined;
                return <div className="chain-choice" key={key}>
                  <label className="checkbox-line"><input type="checkbox" checked={chainChoices.includes(key)} onChange={event => setChainChoices(previous => event.target.checked ? [...previous, key] : previous.filter(value => value !== key))} />
                    {visibleColor && <i className="chain-color-swatch" style={{ backgroundColor: visibleColor }} aria-label={`Viewer colour ${visibleColor}`} title={`Actual viewer colour ${visibleColor}`} />}
                    <span>{chain.sourceChain === chain.copyId ? `Chain ${chain.sourceChain}` : `Chain ${chain.sourceChain} · copy ${chain.copyId}`}</span></label>
                  {visibleColor && <button className="text-link" type="button" onClick={() => setChainFocus(previous => ({ chainId: chain.copyId, serial: (previous?.serial ?? 0) + 1 }))}>Focus in viewer</button>}
                  {chain.copyId !== chain.sourceChain && previewFailure && <small>This assembly copy cannot be focused until its exact view loads.</small>}
                </div>;
              })}
              {relevantPartners.length > 0 && <>
                <span className="field-label">Other molecules in this selection</span>
                {relevantPartners.map(partner => {
                  const copies = assembly?.chainCopies.filter(copy => copy.sourceChain === partner.chain &&
                    (!partner.chain || !proteinSourceChains.has(partner.chain) ||
                      chosenChains.some(selected => selected.copyId === copy.copyId))) ?? [];
                  return <div className="partner-choice" key={partner.sourceId}>
                  <strong>{partnerLabel(partner)}</strong>
                  {copies.length > 1 && <small>This choice applies to every assembly copy ({copies.map(copy => copy.copyId).join(', ')}). Independently different occupancy needs a separately identified structure.</small>}
                  {copies.length === 1 && copies[0].copyId !== partner.chain && <small>Applies to copy {copies[0].copyId}.</small>}
                  <div className="inline-fields">
                    <label className="checkbox-line"><input type="radio" name={`partner-${partner.sourceId}`} checked={partnerChoices[partner.sourceId]?.retain === true} onChange={() => setPartnerChoices(previous => ({ ...previous, [partner.sourceId]: { sourceId: partner.sourceId, retain: true } }))} />Keep</label>
                    <label className="checkbox-line"><input type="radio" name={`partner-${partner.sourceId}`} checked={partnerChoices[partner.sourceId]?.retain === false} onChange={() => setPartnerChoices(previous => ({ ...previous, [partner.sourceId]: { sourceId: partner.sourceId, retain: false } }))} />Exclude</label>
                  </div>
                </div>;})}
              </>}
              {omittedPartners.length > 0 && <details className="source-omissions"><summary>Other source contents ({omittedPartners.length})</summary>
                <p>These components are outside the currently selected protein copies or chosen assembly membership; they are not preparation dispositions.</p>
                {omittedPartners.map(partner => <p key={partner.sourceId}>{partnerLabel(partner)}</p>)}
              </details>}
              {selectedProteinBackboneIncomplete && <div className="notice warning">The selected protein chains contain observed residues with incomplete backbone heavy-atom coordinates. Assessment may refuse them; no missing backbone is silently reconstructed.</div>}
              {ambiguousResidues.length > 0 && <>
                <span className="field-label">Observed alternate locations · choose each explicitly</span>
                <div className="altloc-list">
                  {ambiguousResidues.map(residue => <label key={residueKey(residue.address)} className="altloc-row">
                    <span className="tabular">{residue.name} · {residue.address.chain}{residue.address.residue}{residue.address.insertionCode} · copy {residue.address.copyId}</span>
                    <select className="select-input" value={altlocChoices[residueKey(residue.address)] ?? ''} onChange={event => setAltlocChoices(previous => {
                      const next = { ...previous };
                      if (event.target.value === '') delete next[residueKey(residue.address)];
                      else next[residueKey(residue.address)] = Number(event.target.value);
                      return next;
                    })}>
                      <option value="">Choose location</option>
                      {residue.alternateLocations.map((altloc, index) => <option value={index} key={`${index}:${altloc}`}>{altloc || '(unlabelled)'}</option>)}
                    </select>
                  </label>)}
                </div>
              </>}
              {chosenChains.length > 0 && <div className="hint-box" aria-label="Protein membership to assess">
                <strong>Protein to assess</strong><br />{sourceModelLabel(model, state.sourceModels.length)} · {assemblyChoice === 'deposited' ? 'deposited coordinates' : `biological assembly ${assemblyChoice}`}<br />
                Selected chains: {chosenChains.map(chain => chain.sourceChain === chain.copyId ? `Chain ${chain.sourceChain}` : `Chain ${chain.sourceChain}, copy ${chain.copyId}`).join(', ')}
                {chosenPartners.some(partner => partner.retain) && <><br />Kept: {chosenPartners.filter(partner => partner.retain).map(partner => partnerLabel(relevantPartners.find(item => item.sourceId === partner.sourceId)!)).join('; ')}</>}
                {chosenPartners.some(partner => !partner.retain) && <><br />Excluded: {chosenPartners.filter(partner => !partner.retain).map(partner => partnerLabel(relevantPartners.find(item => item.sourceId === partner.sourceId)!)).join('; ')}</>}
                {chosenAltlocs.length > 0 && <><br />Observed conformers: {chosenAltlocs.map(choice => `${choice.residue.chain}[${choice.residue.copyId}]:${choice.residue.residue}${choice.residue.insertionCode} → ${choice.altloc || '(unlabelled)'}`).join('; ')}</>}
              </div>}
              <div className="button-row"><ActionButton state={state} kind="selectProteinModel" busy={busy || chosenChains.length === 0 || !allPartnersDecided || !allAltlocsChosen} onClick={() => void command('selectProteinModel', { modelIndex: model.index, biologicalAssemblyId: assemblyChoice === 'deposited' ? null : assemblyChoice, chains: chosenChains, partners: chosenPartners, alternateLocations: chosenAltlocs })}>Assess selected protein</ActionButton></div>
              <ActionFeedback kind="selectProteinModel" pending={pendingAction} notice={actionNotice} />
              {action(state, 'selectProteinModel')?.reason && <p className="action-reason">{action(state, 'selectProteinModel')?.reason}</p>}
              {chosenChains.length === 0 && <p className="help-text">Choose at least one observed chain copy.</p>}
              {!allAltlocsChosen && <p className="help-text">Choose each observed alternate location before assessment.</p>}
              </>}
            </>}
          </>}
          {state.protein && <><div className={`protein-standing ${state.protein.status === 'assessed' ? 'established' : 'unestablished'}`} role="status">
            <strong>{state.protein.status === 'assessed' ? 'Assessed prepared protein' :
              preparationReview?.preparationStanding === 'preparing' ? 'Preparing and checking protein…' :
                preparationReview?.preparationStanding === 'failed' ? 'Protein preparation did not complete' :
                  preparationReview?.preparationStanding === 'blocked' || state.protein.status === 'declined' ? 'Prepared protein not established' :
                  preparationReview ? 'Protein review in progress' : 'Prepared protein not established'}</strong>
            <span>{state.protein.summary}</span>
            <small className="tabular">Study revision {state.study?.number ?? 'unknown'} · subject {state.protein.subjectId}</small>
          </div><div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.protein!.subjectId)}>Inspect protein</ActionButton>
              {state.protein.candidateId && <ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.protein!.candidateId!)}>Inspect unqualified candidate</ActionButton>}
            </div>
            {(state.protein.prediction || state.protein.sourceGeometry || state.protein.geometry) &&
              <details className="review-context"><summary>Source confidence and geometry observations</summary>
                <PredictionSummary prediction={state.protein.prediction} label={state.protein.status === 'assessed' ? 'Assessed protein' : 'Selected protein'} />
                <GeometrySummary geometry={state.protein.sourceGeometry} label="Selected source geometry" />
                <GeometrySummary geometry={state.protein.geometry} label="Prepared protein geometry" />
              </details>}
          </>}
        </section>

        <section className="rail-section" id="membrane-workflow" hidden={activeArea !== 'membrane'}>
          {state.membrane && state.membrane.status !== 'proposed' && <section className={`membrane-task-account ${state.membrane.status}`} aria-label="Membrane task outcome">
            <h3>{state.membrane.status === 'assessed' ? 'Membrane ready for placement' :
              state.membrane.status === 'assessing' ? 'Checking membrane…' :
              state.membrane.status === 'unavailable' ? 'Membrane check could not finish' : 'Membrane support not established'}</h3>
            <p>Upper: {fractionSummary(state.membrane.upper)}<br />Lower: {fractionSummary(state.membrane.lower)}</p>
            {state.membrane.status === 'assessing' && <p role="status"><span className="activity-spinner" aria-hidden="true" />Checking the chosen composition…</p>}
            {state.membrane.reason && state.membrane.status !== 'assessed' && <p>{state.membrane.reason}</p>}
            {state.membrane.status === 'assessed' && <button className="button primary" type="button" onClick={() => showWorkArea(state.protein?.status === 'assessed' ? 'placement' : 'protein')}>
              {state.protein?.status === 'assessed' ? 'Position protein' : 'Prepare protein'}</button>}
            {state.membrane.status === 'unavailable' && <ActionButton state={state} kind="adoptMembrane" busy={busy} onClick={() => void command('adoptMembrane', { modelId: state.membrane!.modelId })}>Retry membrane check</ActionButton>}
          </section>}
              <h3 className="section-title">Membrane composition draft</h3>
          {(state.placement || state.stages.length > 0) && <p className="change-impact">Changing the membrane starts a new study revision. Placement must be reassessed; completed stages remain attached to their original revision.</p>}
          <p className="help-text">Specify each physical leaflet. These percentages describe an intended model, not achieved molecule counts.</p>
          {state.availableLipids.length === 0 && <div className="hint-box">No lipid catalogue is currently qualified for selection. Enter an exact species ID below to have its support assessed without substitution.</div>}
          <FractionEditor label="Upper" rows={upper} onChange={setUpper} options={state.availableLipids} />
          <FractionEditor label="Lower" rows={lower} onChange={setLower} options={state.availableLipids} />
          <div className="button-row"><ActionButton state={state} kind="proposeMembrane" busy={busy || !membraneReady} onClick={() => void command('proposeMembrane', { upper: toFractions(upper), lower: toFractions(lower) })}>Propose membrane model</ActionButton></div>
          <ActionFeedback kind="proposeMembrane" pending={pendingAction} notice={actionNotice} />
          {action(state, 'proposeMembrane')?.reason && <p className="action-reason">{action(state, 'proposeMembrane')?.reason}</p>}
          {!membraneReady && <p className="input-guidance">{membraneInputIssue(upper, lower)}</p>}
          <p className="help-text">Zero-percent rows are absent from the proposal; support is assessed after adoption.</p>
          {state.membrane?.status === 'proposed' && <><div className="hint-box membrane-proposal-context" aria-label="Identified membrane proposal and support">
            <strong>Proposal awaiting adoption</strong>
            <span>Upper: {fractionSummary(state.membrane.upper)}</span>
            <span>Lower: {fractionSummary(state.membrane.lower)}</span>
            {!membraneDraftMatchesProposal && <p className="change-impact">The editable draft differs from this proposal. Adopting acts on the displayed proposal; propose the edited draft first if that is your intended composition.</p>}
            {state.membrane.reason && <p className="membrane-support-reason">{state.membrane.reason}</p>}
            {state.membrane.limitations.length > 0 && <p>{state.membrane.limitations.join('; ')}</p>}
            <details><summary>Proposal identity and provenance</summary><p className="tabular">{state.membrane.modelId}</p>
              {state.membrane.scientificPurpose?.trim() && <p>Legacy recorded purpose: {state.membrane.scientificPurpose}</p>}
              {state.membrane.policyId && <p>Policy {state.membrane.policyId} · {state.membrane.policyVersion}</p>}
            </details>
          </div><div className="button-row"><ActionButton state={state} kind="adoptMembrane" busy={busy} onClick={() => void command('adoptMembrane', { modelId: state.membrane!.modelId })}>Adopt and assess displayed proposal</ActionButton></div></>}
          <ActionFeedback kind="adoptMembrane" pending={pendingAction} notice={actionNotice} />
          {state.membrane?.status === 'proposed' && <Prerequisite state={state} kind="adoptMembrane" onOpen={showWorkArea} />}
          {state.study && <p className="help-text tabular">Fixed conditions: pH {state.study.conditions.nominalPh}; NaCl {state.study.conditions.targetNaClMolar} M; temperature {state.study.conditions.optionalTemperatureKelvin} K.</p>}
        </section>

        <section className="rail-section" id="placement-workflow" hidden={activeArea !== 'placement'}>
              <h3 className="section-title">Position the prepared protein</h3>
          {placementTask && placementTask.standing !== 'ready' && <section className={`placement-task-account ${placementTask.standing}`} aria-label="Placement task outcome">
            <strong>{placementTask.standing === 'obtaining' ? 'Obtaining a position…' :
              placementTask.standing === 'assessing' ? 'Assessing the exact pair…' :
              placementTask.standing === 'noProposal' ? 'No position established' :
              placementTask.standing === 'supported' ? placementTask.adopted ? 'Position selected' : 'Position ready to use' :
              placementTask.standing === 'unsupported' ? 'Position failed a technical check' :
              placementTask.standing === 'notEstablished' ? 'Could not check this position' : 'Review current position'}</strong>
            <p>{placementTask.message}</p>
            {['obtaining', 'assessing'].includes(placementTask.standing) && <p role="status"><span className="activity-spinner" aria-hidden="true" />The current pair remains identified while this operation runs.</p>}
            {placementTask.latestAttemptIssue && <p className="notice warning" role="alert">The latest new-position attempt did not replace the current proposal: {placementTask.latestAttemptIssue}</p>}
            {placementTask.standing === 'supported' && placementTask.adopted && <button className="button primary" type="button" onClick={() => showWorkArea('preparation')}>Continue to system preparation</button>}
            <details><summary>Exact pair and proposal identity</summary><p className="tabular">Prepared protein {placementTask.preparedProteinId}<br />Membrane {placementTask.membraneModelId}{placementTask.proposalId && <><br />Proposal {placementTask.proposalId}</>}</p></details>
          </section>}
          <p className="help-text">Position the complete prepared construct in the chosen membrane frame. Molecular rotation changes coordinates; rotating the viewer does not.</p>
          <label className="field-label" htmlFor="orientation-route">Starting position method</label>
          <select id="orientation-route" className="select-input" value={orientationRoute}
            onChange={event => { setOrientationRoute(event.target.value as 'manual' | 'ppm' | 'opm');
              if (event.target.value === 'manual') setManualDraftTouched(true); }}>
            <option value="manual">Set position directly</option>
            <option value="ppm">Use local PPM orientation estimate</option>
            <option value="opm">Use matching OPM reference</option>
          </select>
          {orientationRoute === 'manual' ? <div className="manual-position-inputs" role="group" aria-label="User-defined protein position">
            <label className="field-label" htmlFor="starting-position">Starting position</label>
            <select id="starting-position" className="select-input" value={startingPosition}
              onChange={event => { setStartingPosition(event.target.value as PlacementStartingPosition);
                setManualDraftTouched(true); }}>
              <option value="center">Centered in membrane frame</option>
              <option value="upper">Beside upper headgroups · 2 Å starting gap</option>
              <option value="lower">Beside lower headgroups · 2 Å starting gap</option>
            </select>
            <div className="position-coordinate-grid">
              {([["Move X (Å)", offsetX, setOffsetX], ["Move Y (Å)", offsetY, setOffsetY], ["Move Z (Å)", offsetZ, setOffsetZ],
                ["Rotate X (°)", rotationX, setRotationX], ["Rotate Y (°)", rotationY, setRotationY],
                ["Rotate Z (°)", rotationZ, setRotationZ]] as const).map(([label, value, setValue]) =>
                <label key={label}><span className="field-label">{label}</span><input className="number-input tabular" type="number" step="0.1"
                  value={value} onChange={event => { setValue(event.target.value); setManualDraftTouched(true); }} /></label>)}
            </div>
            <p className="help-text">X and Y run across the membrane; Z follows its normal. Rotations act around the complete construct’s heavy-atom center in X, then Y, then Z order.</p>
            <button className="button compact" type="button" onClick={() => {
              setOffsetX('0'); setOffsetY('0'); setOffsetZ('0');
              setRotationX('0'); setRotationY('0'); setRotationZ('0'); setManualDraftTouched(true);
            }}>Reset movement</button>
          </div> : <div className="optional-orientation-inputs">
            <p className="help-text">This method offers an orientation estimate. You can return to direct positioning at any time.</p>
            <label className="field-label" htmlFor="topology-kind">Method topology input</label>
            <select id="topology-kind" className="select-input" value={topologyKind} onChange={event => {
              const value = event.target.value;
              if (value === '' || value === 'membrane-spanning' || value === 'one-surface-associated') {
                setTopologyKind(value); setPhysicalSide(value === 'membrane-spanning' ? 'both' : '');
              }
            }}><option value="">Choose for this method</option>
              <option value="membrane-spanning">Membrane-spanning</option>
              <option value="one-surface-associated">Associated with one surface</option></select>
            {topologyKind === 'one-surface-associated' && <><label className="field-label" htmlFor="physical-side">Physical side</label>
              <select id="physical-side" className="select-input" value={physicalSide} onChange={event =>
                setPhysicalSide(event.target.value as PlacementPhysicalSide | '')}>
                <option value="">Choose physical side</option><option value="upper">Upper</option><option value="lower">Lower</option>
              </select></>}
            {orientationRoute === 'ppm' && <><label className="field-label" htmlFor="ppm-nterminal-side">N-terminus in the PPM convention</label>
              <select id="ppm-nterminal-side" className="select-input" value={ppmNterminalSide} onChange={event =>
                setPpmNterminalSide(event.target.value as PpmNterminalSide | '')}>
                <option value="">Choose N-terminal assignment</option><option value="in">Inside</option><option value="out">Outside</option>
              </select><p className="help-text">PPM inside/outside is a biological N-terminal assignment; it does not mean upper/lower in this viewer.</p></>}
          </div>}
          <div className="button-row">{(orientationRoute !== 'manual' || !manualDraftMatchesProposal) && <ActionButton state={state} kind="proposePlacement" busy={busy ||
            (orientationRoute === 'manual' ? !manualInputReady : !topologyKind || !physicalSide ||
              orientationRoute === 'ppm' && (!ppmNterminalSide || state.placementMethods?.find(item => item.method === 'PPM')?.standing !== 'configured') ||
              orientationRoute === 'opm' && state.placementMethods?.find(item => item.method === 'OPM')?.standing !== 'lookupEligible')}
            onClick={() => { manualPositionAttempt.current = null; void command('proposePlacement', orientationRoute === 'manual'
              ? { orientationRoute: 'manual', startingPosition,
                offsetXAngstrom: Number(offsetX), offsetYAngstrom: Number(offsetY), offsetZAngstrom: Number(offsetZ),
                rotationXDegrees: Number(rotationX), rotationYDegrees: Number(rotationY), rotationZDegrees: Number(rotationZ) }
              : { orientationRoute, topologyKind, physicalSide, ppmNterminalSide: orientationRoute === 'ppm' ? ppmNterminalSide : null }); }}>
              {orientationRoute === 'manual' ? state.placementTask?.latestAttemptIssue ? 'Retry position check' : 'Check position now' : 'Calculate orientation estimate'}</ActionButton>}</div>
          <ActionFeedback kind="proposePlacement" pending={pendingAction} notice={actionNotice} />
          <Prerequisite state={state} kind="proposePlacement" onOpen={showWorkArea} />
          <details className="method-availability"><summary>Orientation method availability</summary>
            {state.placementMethods?.map(method => <p key={method.method}><strong>{method.method}</strong> · {method.standing === 'lookupEligible' ? 'RCSB lookup can be attempted; a matching record is not guaranteed' :
              method.standing === 'configured' ? 'Local installation verified' : method.reason ?? 'Unavailable'}</p>)}
          </details>
          {action(state, 'proposePlacement')?.enabled && orientationRoute === 'manual' && !manualInputReady &&
            <p className="input-guidance">Enter finite movement and rotation values.</p>}
          {orientationRoute !== 'manual' && state.placementMethods?.find(item => item.method.toLowerCase() === orientationRoute)?.standing === 'unavailable' &&
            <p className="input-guidance">{state.placementMethods.find(item => item.method.toLowerCase() === orientationRoute)?.reason}</p>}
          {state.placement && <><div className="hint-box membrane-proposal-context" aria-label="Position proposal and scientific support">
            <strong>{state.placement.status === 'supported' ? 'Position ready to use' :
              state.placement.status === 'unsupported' ? 'Position failed a technical check' :
                state.placement.status === 'assessing' ? 'Checking this position…' :
                'Position check incomplete'}</strong>
            <span>Complete prepared construct · chosen membrane frame.</span>
            <p>{state.placement.status === 'supported'
              ? 'The exact positioned construct and frame passed their technical checks. Construction checks the assembled system separately.'
              : state.placement.reason}</p>
            {state.placement.transform && <p className="tabular">{state.placement.transform.startingPosition} start · applied translation ({
              state.placement.transform.appliedTranslationXAngstrom.toFixed(2)}, {
              state.placement.transform.appliedTranslationYAngstrom.toFixed(2)}, {
              state.placement.transform.appliedTranslationZAngstrom.toFixed(2)}) Å · rotations ({
              state.placement.transform.rotationXDegrees}, {state.placement.transform.rotationYDegrees}, {
              state.placement.transform.rotationZDegrees})°</p>}
            {state.study?.adoptedPlacementProposalId === state.placement.proposalId && <span>Selected in the current study revision.</span>}
            {state.placement.transform && !manualDraftMatchesProposal && <p className="change-impact">The movement draft differs from this checked position. Its technical check will update after editing stops.</p>}
            <details><summary>Proposal identity and basis</summary><p className="tabular">Proposal {state.placement.proposalId}</p>
              <p className="tabular">Prepared protein {state.placement.preparedProteinId ?? 'not identified'}<br />Membrane model {state.placement.membraneModelId ?? 'not identified'}</p>
              {state.placement.policyId && <p>Policy {state.placement.policyId} · {state.placement.policyVersion}</p>}
              {state.placement.limitations.map(limit => <p key={limit}>{limit}</p>)}</details>
          </div><div className="button-row">{state.study?.adoptedPlacementProposalId !== state.placement.proposalId && <ActionButton state={state} kind="adoptPlacement" busy={busy || (!!state.placement.transform && !manualDraftMatchesProposal)} onClick={() => void command('adoptPlacement', { proposalId: state.placement!.proposalId })}>Use this position</ActionButton>}<ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.placement!.proposalId)}>View positioned protein</ActionButton></div>
            <ActionFeedback kind="adoptPlacement" pending={pendingAction} notice={actionNotice} />
            {state.study?.adoptedPlacementProposalId !== state.placement.proposalId && <Prerequisite state={state} kind="adoptPlacement" onOpen={showWorkArea} />}
            <PlacementPredictionSummary prediction={state.placement.prediction} />
            {!state.placement.transform && <div className="change-card">
              <strong>Bounded placement correction</strong>
              <p>Adjust this method-derived proposal, then check the complete positioned construct again.</p>
              <div className="inline-fields">
                <label><span className="field-label">Depth shift (Å)</span><input className="number-input tabular" type="number" step="0.1" value={depthShift} onChange={event => setDepthShift(event.target.value)} /></label>
                <label><span className="field-label">Tilt about X (°)</span><input className="number-input tabular" type="number" step="0.1" value={tiltX} onChange={event => setTiltX(event.target.value)} /></label>
                <label><span className="field-label">Tilt about Y (°)</span><input className="number-input tabular" type="number" step="0.1" value={tiltY} onChange={event => setTiltY(event.target.value)} /></label>
                <label><span className="field-label">Rotation about normal (°)</span><input className="number-input tabular" type="number" step="0.1" value={rotationNormal} onChange={event => setRotationNormal(event.target.value)} /></label>
              </div>
              <div className="button-row"><ActionButton state={state} kind="revisePlacement" busy={busy || !correctionReady} onClick={() => void command('revisePlacement', { proposalId: state.placement!.proposalId, depthShiftAngstrom: Number(depthShift), tiltAboutXDegrees: Number(tiltX), tiltAboutYDegrees: Number(tiltY), rotationAboutNormalDegrees: Number(rotationNormal) })}>Reassess corrected placement</ActionButton></div>
              <ActionFeedback kind="revisePlacement" pending={pendingAction} notice={actionNotice} />
              {action(state, 'revisePlacement')?.reason && <p className="action-reason">{action(state, 'revisePlacement')?.reason}</p>}
            </div>}
          </>}
        </section>

        <section className="rail-section" id="preparation-workflow" hidden={activeArea !== 'preparation'}>
          <h3 className="section-title">Construct and minimize</h3>
          <div className="button-row"><ActionButton state={state} kind="startPreparation" busy={busy} variant="primary" onClick={() => void startPreparation()}>Construct system</ActionButton></div>
          <ActionFeedback kind="startPreparation" pending={pendingAction} notice={actionNotice} />
          <Prerequisite state={state} kind="startPreparation" onOpen={showWorkArea} />
          {state.attempt && <div className="attempt-account">
            <strong>Current system attempt</strong>
            <span>{readable(state.attempt.status)}{state.attempt.stageKind && ` · ${readable(state.attempt.stageKind)}`}</span>
            {activeAttempt && <p className="attempt-running" role="status"><span className="activity-spinner" aria-hidden="true" />
              {state.attempt.stageKind === 'Minimization' ? 'Minimization in progress' : 'Construction in progress'}</p>}
            {state.attempt.stopRequested && <p className="attempt-running" role="status">Stop requested · waiting for the worker’s observed outcome.</p>}
            <p>{state.attempt.message}</p>
            <details><summary>Exact attempt identity</summary><p className="tabular">Attempt {state.attempt.attemptId}<br />Study revision {state.attempt.studyRevisionId ?? 'not established'}</p></details>
            {state.attempt.progress !== null && <progress max="1" value={state.attempt.progress} aria-label="Observed attempt progress" />}
            {state.attempt.constructed && <div className="candidate-summary" aria-label="Checked constructed candidate">
              <strong>Checked candidate · {state.attempt.constructed.atomCount.toLocaleString()} atoms</strong>
              <span>Achieved lipids: {state.attempt.constructed.achievedComposition.map(item =>
                `${item.physicalSide} ${item.speciesId} ${item.count}`).join(' · ') || 'none reported'}</span>
              <span>Cell: {state.attempt.constructed.cellAngstrom.map(value => `${value.toFixed(1)} Å`).join(' × ')}</span>
              <span>Water {state.attempt.constructed.waterCount.toLocaleString()} · Na⁺ {state.attempt.constructed.sodiumCount} · Cl⁻ {state.attempt.constructed.chlorideCount}</span>
              {awaitingMinimization && <p>Continuing explicitly authorizes minimization of this same identified candidate. Construction alone is not a completed minimized stage.</p>}
              <details><summary>Exact candidate identity and conditions</summary><p className="tabular">Subject {state.attempt.constructed.subjectId}</p>
                <p>{state.attempt.constructed.conditionsTreatment}</p></details>
            </div>}
            <div className="button-row">
              <button className="button" type="button" onClick={() => setAttemptReviewRequested(true)}>Review attempt</button>
              {awaitingMinimization && state.attempt.constructed && <ActionButton state={state} kind="continueMinimization" subjectId={state.attempt.attemptId}
                busy={busy} variant="primary" onClick={() => void continueMinimization(state.attempt!)}>Authorize minimization of this candidate</ActionButton>}
              <ActionButton state={state} kind="stopAttempt" busy={busy} variant="danger" onClick={() => void command('stopAttempt', { attemptId: state.attempt!.attemptId })}>{awaitingMinimization ? 'Decline candidate' : 'Stop unfinished work'}</ActionButton>
            </div>
            <ActionFeedback kind="continueMinimization" pending={pendingAction} notice={actionNotice} />
            <ActionFeedback kind="stopAttempt" pending={pendingAction} notice={actionNotice} />
            {awaitingMinimization && <Prerequisite state={state} kind="continueMinimization" subjectId={state.attempt.attemptId} onOpen={showWorkArea} />}
          </div>}
        </section>

        <section className="rail-section" id="results-workflow" hidden={activeArea !== 'results'}>
          <h3 className="section-title">Completed stages</h3>
          {state.stages.length === 0 ? <><p className="help-text">No completed stage is established. You can inspect this area before its prerequisites are met.</p>
            <button className="button compact" type="button" onClick={() => showWorkArea('preparation')}>Open Preparation</button></> :
            <div className="result-stage-list" aria-label="Completed stage choices">{state.stages.map(stage =>
              <button key={stage.stageId} type="button" className={`source-item ${selectedStage?.stageId === stage.stageId ? 'selected' : ''}`}
                disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true}
                onClick={() => { setSelectedStageId(stage.stageId); void inspectSubject(stage.stageId); }}>
                <span className="item-title">{readable(stage.kind)} · {readable(stage.status)}</span>
                <span className="item-detail">{stage.studyRevisionId === state.study?.id ? 'Current study' : 'Historical study'} · checked result</span>
                <span className="item-detail">{stage.assessment ? `${stage.assessment.currentlyApplicable ? '' : 'Historical '}${readable(stage.assessment.qualification)}` : 'Assessment not established'}</span>
              </button>)}</div>}
          {state.stages.length > 0 && !selectedStage && <p className="help-text">Select a completed stage to inspect its structure, assessment, and export options.</p>}
          {selectedStage && <div className="stage-export">
            <strong>{readable(selectedStage.kind)} · {readable(selectedStage.status)}</strong>
            <details><summary>Exact stage identity</summary><p className="tabular">Stage {selectedStage.stageId}<br />Attempt {selectedStage.attemptId}<br />Study revision {selectedStage.studyRevisionId}</p></details>
            {state.inspection?.subjectId !== selectedStage.stageId && <p className="help-text">The selected {readable(selectedStage.kind).toLowerCase()} result differs from the molecular view ({viewerSubjectHeading(state.inspection,
              state.stages.find(stage => stage.stageId === state.inspection?.subjectId), state.study?.selectedSourceLabel)}). Assessment and export below belong to the selected result.</p>}
            <p>{selectedStage.assessment ? `${selectedStage.assessment.currentlyApplicable ? readable(selectedStage.assessment.qualification) : `Historical ${readable(selectedStage.assessment.qualification)}; assessment not current`} — ${selectedStage.assessment.reason}` : 'No scientific assessment established for this stage.'}</p>
            <div className="button-row"><ActionButton state={state} kind="exportStage" subjectId={selectedStage.stageId} busy={busy || exportBusy} onClick={() => void exportStage(selectedStage.stageId)}>Export this completed stage</ActionButton></div>
            {exportTargetId === selectedStage.stageId && (exportBusy && pendingAction !== 'exportStage' ?
              <p className="action-feedback pending" role="status"><span className="activity-spinner" aria-hidden="true" />Transferring and checking stage {selectedStage.stageId}…</p> :
              <ActionFeedback kind="exportStage" pending={pendingAction} notice={actionNotice} />)}
            <Prerequisite state={state} kind="exportStage" subjectId={selectedStage.stageId} onOpen={showWorkArea} />
            {selectedExportFault && <button className="button compact" type="button" onClick={() => document.querySelector('.export-validation')?.scrollIntoView({ block: 'start' })}>Inspect export issue</button>}
            {selectedStage.kind === 'Minimization' && <><div className="button-row"><ActionButton state={state} kind="requestEquilibration" subjectId={selectedStage.stageId} busy={busy} onClick={() => void command('requestEquilibration', { stageId: selectedStage.stageId })}>Request optional equilibration from this stage</ActionButton></div><ActionFeedback kind="requestEquilibration" pending={pendingAction} notice={actionNotice} /><Prerequisite state={state} kind="requestEquilibration" subjectId={selectedStage.stageId} onOpen={showWorkArea} /></>}
          </div>}
        </section>
      </aside>

      <ConnectedStructuralInspection state={displayedState} reviewAttempt={reviewAttempt} onSelectFocus={annotationId => void command('setInspectionFocus', { annotationId })}
        onInspectSubject={subjectId => { void inspectSubject(subjectId); }}
        onChainColors={(subjectId, structureUrl, colors) => setSourceChainColors({ subjectId, structureUrl, colors })} chainFocus={chainFocus}
        sourcePreviewLabel={sourcePreviewLabel}
        activeArea={activeArea} membraneDraft={{ upper: toFractions(upper), lower: toFractions(lower) }}
        membraneDraftMatchesProposal={membraneDraftMatchesProposal}
        requestedAreaView={pendingAreaSubject ? { area: workAreas.find(item => item.id === activeArea)!.label,
          subjectId: pendingAreaSubject, failure: viewRestoreFailure?.subjectId === pendingAreaSubject ? viewRestoreFailure.reason : null } : null}
        structureReload={structureReload} structureStatus={structureLoad} onStructureLoad={onStructureLoad}
        requestedSource={sourceIntent && sourceIntent.phase !== 'failed'
          ? { label: sourceIntent.label, phase: sourceIntent.phase } : previewPending && sourcePreviewLabel
            ? { label: sourcePreviewLabel, phase: previewFailure ? 'previewFailed' : 'previewing', reason: previewFailure ?? undefined } : null}
        exportFault={selectedExportFault}
        connectionMessage={communication ?? (streamInterrupted ? 'The live connection is interrupted. Showing the last account read from the host; progress may change until reconnection.' : null)}
        onRefreshAccount={() => void refresh()}
        placementOutcome={activeArea === 'placement' && placementTask && !state.placement &&
          ['obtaining', 'noProposal'].includes(placementTask.standing) ?
          <section className={`account-card placement-outcome-panel ${placementTask.standing}`} aria-label="Current placement outcome">
            <strong>{placementTask.standing === 'obtaining' ? 'Position being obtained' : 'No reviewable position established'}</strong>
            <p>{placementTask.message}</p>
            {placementTask.standing === 'obtaining' && <p role="status"><span className="activity-spinner" aria-hidden="true" />Awaiting the identified orientation route.</p>}
            <p className="help-text">The molecular viewer still shows {viewerSubjectHeading(state.inspection,
              state.stages.find(stage => stage.stageId === state.inspection?.subjectId), state.study?.selectedSourceLabel)}. This is not a placement proposal.</p>
            <details><summary>Exact pair identity</summary><p className="tabular">Prepared protein {placementTask.preparedProteinId}<br />Intended membrane {placementTask.membraneModelId}<br />Study revision {placementTask.studyRevisionId}</p></details>
          </section> : null}
        proteinDraft={activeArea === 'protein' && sourceId && !sourceIntent && !proteinTask ?
          <section className="account-card protein-draft-account" aria-label="Protein selection draft">
            <dl className="detail-grid">
              <dt>Coordinate model</dt><dd>{model ? sourceModelLabel(model, state.sourceModels.length) : 'Choose a coordinate model'}</dd>
              <dt>Assembly</dt><dd>{!model ? 'Choose a coordinate model first' : !assemblyResolved ?
                'Choose deposited coordinates or a biological assembly' : assemblyChoice === 'deposited' ?
                  'Deposited coordinates' : `Biological assembly ${assemblyChoice}`}</dd>
              <dt>Protein chains</dt><dd>{chosenChains.length ? chosenChains.map(chain => chain.sourceChain === chain.copyId ?
                `Chain ${chain.sourceChain}` : `Chain ${chain.sourceChain}, copy ${chain.copyId}`).join(' and ') :
                'Choose chain copies to retain'}</dd>
              {assemblyResolved && chosenChains.length > 0 && <><dt>Kept partners</dt><dd>{chosenPartners.filter(partner => partner.retain).length ?
                chosenPartners.filter(partner => partner.retain).map(partner => partnerLabel(relevantPartners.find(item => item.sourceId === partner.sourceId)!)).join('; ') : 'None'}</dd>
                <dt>Excluded partners</dt><dd>{chosenPartners.filter(partner => !partner.retain).length ?
                  chosenPartners.filter(partner => !partner.retain).map(partner => partnerLabel(relevantPartners.find(item => item.sourceId === partner.sourceId)!)).join('; ') : 'None'}</dd>
                {relevantPartners.length > chosenPartners.length && <><dt>Still to decide</dt><dd>{relevantPartners.filter(partner => !partnerChoices[partner.sourceId]).map(partner => partnerLabel(partner)).join('; ')}</dd></>}
              </>}
            </dl>
            <p className="help-text">This is a draft selection. Assessment checks whether its exact members and observed structure can be prepared.</p>
          </section> : null}
        proteinOutcome={proteinOutcomeVisible && proteinTask ? <section className={`protein-outcome-panel ${proteinTask.standing}`} aria-label="Current protein result">
          <div className="protein-outcome-summary"><strong>{proteinTask.standing === 'assessed' ? 'Preparation complete' :
            proteinTask.standing === 'preparing' ? 'Preparing and checking protein…' :
            proteinTask.standing === 'failed' ? 'No prepared protein established' :
            proteinTask.standing === 'unavailable' ? 'Outcome not verified' : 'Preparation blocked'}</strong>
            <p>{proteinTask.message}</p></div>
          <p><strong>Selected protein</strong><br />{state.study?.selectedSourceLabel ?? proteinTask.sourceId} · {proteinTask.chains.map(chain =>
            chain.sourceChain === chain.copyId ? `Chain ${chain.sourceChain}` : `Chain ${chain.sourceChain}, copy ${chain.copyId}`).join(' and ')}</p>
          {proteinTask.standing === 'assessed' && proteinTask.preparedProteinId && <>
            <button className="button primary" type="button" disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true}
              onClick={() => void inspectSubject(proteinTask.preparedProteinId!)}>{state.inspection?.subjectId === proteinTask.preparedProteinId ? 'Prepared protein displayed' : 'View prepared protein'}</button>
            <p className="help-text">{state.inspection?.subjectId === proteinTask.preparedProteinId ?
              'The viewer shows the assessed prepared artifact.' :
              'The viewer retains its current subject. Open the prepared result to see the checked coordinates.'}</p>
            {inspectionFailure?.subjectId === proteinTask.preparedProteinId && <p className="review-blocker-inline" role="alert">Prepared protein visualization unavailable: {inspectionFailure.reason}. The assessed result remains established; retry viewing it when available.</p>}
            <p className="help-text">Membrane placement has not been established by protein preparation.</p>
          </>}
          {proteinTask.standing === 'failed' && <p>Confirmed choices remain available under Review confirmed choices. A failed operation has not created an assessed protein.</p>}
          {proteinTask.standing === 'preparing' && <p role="status"><span className="activity-spinner" aria-hidden="true" />The result is being checked. You may work in another area.</p>}
          {state.protein?.findings.length ? <div className="protein-outcome-findings"><strong>Structure findings</strong>
            {state.protein.findings.map(finding => <p key={finding.id}>{finding.consequence} · {finding.meaning}</p>)}</div> : null}
          {(state.protein?.sourceGeometry || state.protein?.geometry || state.protein?.prediction) && <details className="review-evidence-details"><summary>Structure checks and limitations</summary>
            <PredictionSummary prediction={state.protein?.prediction} label="Selected protein" />
            <GeometrySummary geometry={state.protein?.sourceGeometry} label="Source coordinates" />
            <GeometrySummary geometry={state.protein?.geometry} label={proteinTask.standing === 'assessed' ? 'Prepared result' : 'Unqualified candidate'} />
          </details>}
          <details className="review-evidence-details"><summary>Technical details</summary><p className="tabular">Source {proteinTask.sourceId}<br />Intended protein {proteinTask.intendedProteinId}<br />Study revision {proteinTask.studyRevisionId}{proteinTask.preparedProteinId && <><br />Prepared artifact {proteinTask.preparedProteinId}</>}</p></details>
        </section> : null}
        proteinReview={reviewPanelVisible && preparationReview && selectedDecision ? <section className="protein-review-panel" aria-label="Current protein decision">
          <p className="review-site-title">{decisionSiteLabel(selectedDecision, state.sourceModels)}</p>
          <p className="review-site-scope">{decisionKindLabel(selectedDecision.kind)}
            {selectedDecision.partnerResidue && ` · Bond partner: Cysteine ${selectedDecision.partnerResidue.residue}${selectedDecision.partnerResidue.insertionCode}, chain ${selectedDecision.partnerResidue.chain}${selectedDecision.partnerResidue.copyId !== selectedDecision.partnerResidue.chain ? `, copy ${selectedDecision.partnerResidue.copyId}` : ''}`}</p>
          {preparationReview.preparationStanding === 'preparing' && <div className="review-preparation-standing" role="status"><span className="activity-spinner" aria-hidden="true" /><span>{preparationReview.preparationMessage}</span></div>}
          {preparationReview.preparationStanding === 'failed' && <div className="review-preparation-standing failed" role="alert"><strong>Preparation unsuccessful</strong><span>{preparationReview.preparationMessage}</span><span>The choices remain recorded; use Retry protein preparation in the input pane if the cause is transient.</span></div>}
          <p className="review-viewing-context">{preparationReview.inspectionRelation === 'beforePreparation'
            ? state.inspection?.representationKind === 'structuralSource'
              ? 'Viewing the full source structure. This decision concerns the selected protein within it; the reviewed change is not yet applied.'
              : 'Viewing corresponding protein coordinates before preparation. The reviewed change is not yet applied.'
            : preparationReview.inspectionRelation === 'preparedResult' ? 'Viewing the corresponding prepared result.'
              : preparationReview.inspectionRelation === 'unqualifiedCandidate' ? 'Viewing the corresponding candidate that did not establish a prepared protein.'
                : preparationReview.inspectionRelation === 'historical' ? 'Viewing an earlier revision. This decision concerns the current selected protein.'
                  : preparationReview.inspectionRelation === 'otherSubject' ? 'Viewing another molecular subject. This decision concerns the selected protein.'
                    : 'No corresponding molecular view is available. The decision and its scientific information remain identified.'}</p>
          {selectedDecision.chosenProposalId && selectedDecision.options.some(option =>
            option.proposalId === selectedDecision.chosenProposalId && option.disposition !== 'available') &&
            <div className="review-receipt" role="status"><strong>
            {selectedDecision.options.find(option => option.proposalId === selectedDecision.chosenProposalId)?.disposition === 'declined'
              ? selectedDecision.kind === 'heavyAtom' ? 'Required repair rejected for this site' : 'No bond confirmed for this pair' :
                `${selectedDecision.options.find(option => option.proposalId === selectedDecision.chosenProposalId)?.proposedChange ?? 'Choice'} confirmed for this site`}</strong>
            <span>Recorded for this selected protein and study revision.</span></div>}
          {selectedDecision.blocker && <p className="review-blocker-inline" role="alert">{selectedDecision.blocker}</p>}
          {selectedDecision.kind === 'residueState' && <p className="help-text">{preparationPlan?.standing === 'ready'
            ? preparationPlan.choices.find(choice => sameResidue(choice.residue, selectedDecision.residue))?.overridden
              ? 'Your selected state was checked with the rest of the plan. No state has been applied yet.'
              : `The identified method suggested a starting state at pH ${preparationPlan.nominalPh}. Review or override it; no state has been applied.`
            : `Choose this site’s state. Nominal pH ${state.study?.conditions.nominalPh ?? 'unknown'} alone does not establish it.`}</p>}
          <div className="review-options" role="group" aria-label="Available alternatives">
            {selectedDecision.options.map(option => {
              const change = state.protein?.changes.find(item => item.id === option.proposalId);
                const planChoice = preparationPlan?.standing === 'ready' ? preparationPlan.choices.find(choice =>
                  sameResidue(choice.residue, selectedDecision.residue) && choice.variant === option.proposedChange) : undefined;
                const picked = selectedDecision.standing !== 'pending' && selectedDecision.chosenProposalId
                ? selectedDecision.chosenProposalId === option.proposalId : selectedOption?.proposalId === option.proposalId;
              return <div className={`review-option${picked ? ' picked' : ''}`} key={option.proposalId}>
                <label><input type="radio" name={`decision-${selectedDecision.id}`} checked={picked}
                  disabled={selectedDecision.standing === 'confirmed' || option.disposition !== 'available'}
                  onChange={() => { focusRequestGeneration.current += 1; setSelectedOptionId(option.proposalId); }} />
                  <span><strong>{option.proposedChange}</strong>{optionMeaning(option.proposedChange) && <small className="review-option-meaning">{optionMeaning(option.proposedChange)}</small>}<small>{option.disposition === 'notChosen' ? 'Not chosen' :
                    option.disposition === 'confirmed' ? 'Confirmed' : option.disposition === 'declined' ? 'Declined' :
                      planChoice ? planChoice.overridden ? 'Your choice · checked with the method' : 'Suggested by the identified method' :
                        change?.kind === 'disulfide' ? 'Possible bond' : 'Available for review'}</small></span></label>
                {picked && <button className="button compact" type="button"
                  disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true}
                  onClick={() => void focusDecisionOption(option.proposalId)}>Focus on this {selectedDecision.kind === 'disulfide' ? 'pair' : 'residue'}</button>}
              </div>;
            })}
          </div>
          {inspectionFailure && selectedDecision.options.some(option => option.proposalId === inspectionFailure.subjectId) &&
            <p className="review-blocker-inline" role="alert">Local view unavailable: {inspectionFailure.reason} You can retry Focus on this residue. The recorded scientific evidence remains separate.</p>}
          {selectedOption && <>
            <div className="review-evidence-current" aria-label="Selected option explanation and evidence">
              <strong>{selectedOption.informationRole === 'modelAssumption' ? 'Model choice and its limits' : 'Why this change is proposed'}</strong>
              {selectedOption.informationRole !== 'modelAssumption' &&
                state.protein?.changes.find(item => item.id === selectedOption.proposalId)?.rationale &&
                <p>{state.protein.changes.find(item => item.id === selectedOption.proposalId)!.rationale}</p>}
              {selectedOption.informationRole === 'modelAssumption' ? <p>This supported state is a researcher-reviewed model assumption. The app has not measured which state is best at this site.</p> :
                selectedOption.evidence.map(item => <p key={item.id}>{item.observation}<br /><span>Uncertainty: {item.uncertainty}</span></p>)}
              {selectedOption.evidence.length === 0 && <p>No attributable observation is available for this option.</p>}
            </div>
            {preparationPlan?.standing === 'ready' && selectedDecision.kind === 'residueState' ?
              <div className="button-row review-confirm-actions">
                {preparationPlan.choices.some(choice => sameResidue(choice.residue, selectedDecision.residue) &&
                  choice.variant === selectedOption.proposedChange) ?
                  <p className="help-text">This is the current plan choice. The plan action in Protein applies all checked choices together.</p> :
                  <ActionButton state={state} kind="overridePreparationPlanChoice" subjectId={selectedOption.proposalId}
                    busy={busy} variant="primary" onClick={() => void command('overridePreparationPlanChoice',
                      { proposalId: selectedOption.proposalId })}>Use this choice</ActionButton>}
              </div> : selectedDecision.standing !== 'confirmed' && selectedOption.disposition === 'available' && <>
            {selectedOption.startsPreparationOnConfirmation && <p className="review-consequence">Your confirmed choices will be applied and the resulting protein checked.</p>}
            {selectedOption.startsPreparationOnDecline && <p className="review-consequence">This final disposition starts protein preparation and checking.</p>}
            {selectedDecision.kind === 'heavyAtom' && <p className="review-consequence">Declining this required repair blocks preparation of the current selected protein.</p>}
            <div className="button-row review-confirm-actions">
              <button className="button primary" type="button" disabled={busy || !selectedOption.evidenceCurrent ||
                action(state, 'approvePreparationChange', selectedOption.proposalId)?.enabled !== true}
                onClick={() => void confirmDecisionOption(selectedOption.proposalId, true)}>{selectedOption.startsPreparationOnConfirmation
                  ? selectedDecision.kind === 'heavyAtom' ? 'Approve repair and prepare protein' : selectedDecision.kind === 'disulfide' ? 'Confirm bond and prepare protein' :
                    selectedDecision.kind === 'residueState' ? 'Confirm state and prepare protein' : 'Confirm conformer and prepare protein'
                  : selectedDecision.kind === 'heavyAtom' ? 'Approve repair' : selectedDecision.kind === 'disulfide' ? 'Confirm possible bond' : selectedDecision.kind === 'residueState' ? 'Confirm state' : 'Confirm conformer'}</button>
              {(selectedDecision.kind === 'heavyAtom' || selectedDecision.kind === 'disulfide') && <button className="button" type="button" disabled={busy ||
                action(state, 'declinePreparationChange', selectedOption.proposalId)?.enabled !== true}
                onClick={() => void confirmDecisionOption(selectedOption.proposalId, false)}>{selectedDecision.kind === 'heavyAtom' ? 'Reject required repair' :
                  selectedOption.startsPreparationOnDecline ? 'Reject bond and prepare protein' : 'Reject possible bond'}</button>}
            </div>
            {selectedOption.confirmationBlocker && <p className="action-reason">{selectedOption.confirmationBlocker}</p>}
            </>}
          </>}
          <ActionFeedback kind="overridePreparationPlanChoice" pending={pendingAction} notice={actionNotice} />
          <ActionFeedback kind="selectInspectionSubject" pending={pendingAction} notice={actionNotice} />
          {decisionCommandTarget && selectedDecision.options.some(option => option.proposalId === decisionCommandTarget) &&
            pendingAction === 'approvePreparationChange' && <p className="action-feedback pending" role="status"><span className="activity-spinner" aria-hidden="true" />Recording {selectedDecision.options.find(option => option.proposalId === decisionCommandTarget)?.proposedChange} for {decisionSiteLabel(selectedDecision, state.sourceModels)}…</p>}
          {decisionCommandTarget && selectedDecision.options.some(option => option.proposalId === decisionCommandTarget) &&
            actionNotice?.kind === 'approvePreparationChange' && actionNotice.tone === 'error' &&
            <p className="action-feedback error" role="alert">{selectedDecision.options.find(option => option.proposalId === decisionCommandTarget)?.proposedChange} at {decisionSiteLabel(selectedDecision, state.sourceModels)}: {actionNotice.message}</p>}
          {nextUnresolved && <div className="review-next"><button className="button primary" type="button" onClick={goToNextUnresolved}>{preparationPlan?.standing === 'ready' ? 'Next site to review' :
            `Next unresolved ${nextUnresolved.kind === 'heavyAtom' ? 'repair' : nextUnresolved.kind === 'disulfide' ? 'bond decision' : 'site'}`} ›</button>
            {preparationPlan?.standing !== 'ready' && <small>{preparationReview.remainingCount} total obligation{preparationReview.remainingCount === 1 ? '' : 's'} remain</small>}</div>}
          {selectedDecision.standing === 'pending' && !nextUnresolved && preparationPlan?.standing !== 'ready' && <p className="help-text">Review this site’s available alternatives to continue.</p>}
          <details className="review-evidence-details"><summary>Evidence and limitations</summary>
            {selectedOption ? <>
              {selectedOption.evidence.map(item => <div key={item.id} className="review-evidence-item"><strong>{item.source} · {item.method}</strong><p>{item.observation}</p><p>Applicability: {item.applicability}</p><p>Uncertainty: {item.uncertainty}</p></div>)}
              {state.protein?.changes.find(item => item.id === selectedOption.proposalId)?.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
            </> : <p>Choose an option to read its explanation and attributable evidence. Viewing coordinates is optional.</p>}
          </details>
          <details className="review-evidence-details"><summary>Technical provenance</summary><p className="tabular">Study revision {preparationReview.studyRevisionId} · intended protein {preparationReview.intendedProteinId}</p>
            {selectedOption && <p className="tabular">Decision option {selectedOption.proposalId} · recorded decision {selectedOption.decisionId ?? 'none'}</p>}
            {state.inspection && <p className="tabular">Viewed subject {state.inspection.subjectId} · revision {state.inspection.studyRevisionId}</p>}
            <p>Visualization does not approve a change or alter coordinates. The preparation outcome remains separate from this recorded decision.</p></details>
        </section> : null}
        proposalDecision={executionReview && !workflowOpen && !reviewPanelVisible ? <div className="proposal-decision execution-decision" aria-label={selectedStage ? 'Completed stage next steps' : 'Attempt decision'}>
          {selectedStage ? <>
            {selectedExportFault ? <>
              <strong className="export-next-title">Next steps</strong>
              <div className="button-row decision-actions export-decision-actions">
                <button className="button" type="button" onClick={() =>
                  document.querySelector('.export-validation')?.scrollIntoView({ block: 'start' })}>Inspect export issue</button>
                <button className="button primary" type="button" disabled={busy || exportBusy ||
                  action(state, 'exportStage', selectedStage.stageId)?.enabled !== true}
                  onClick={() => void exportStage(selectedStage.stageId)}>Retry export</button>
              </div>
            </> : <>
            <div className="decision-standing"><strong>{selectedStage.kind === 'Minimization' ? 'Completed minimized stage' : 'Completed equilibrated stage'}</strong>
              <span>{selectedStage.assessment ? `${selectedStage.assessment.currentlyApplicable ? '' : 'Historical '}${readable(selectedStage.assessment.qualification)}` : 'Scientific assessment unavailable'}</span></div>
            {selectedStage.assessment && <p className="decision-reason">{selectedStage.assessment.reason}</p>}
            <div className="button-row decision-actions">
              {selectedSourceStage ? <button className="button" type="button"
                disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true}
                onClick={() => void inspectSubject(selectedSourceStage.stageId)}>Select minimized stage</button> :
                <button className="button" type="button" onClick={() => document.querySelector('.execution-findings-account, .execution-assessment-account')?.scrollIntoView({ block: 'start' })}>Inspect findings</button>}
              {action(state, 'exportStage', selectedStage.stageId)?.enabled && <button className="button primary" type="button" disabled={busy || exportBusy} onClick={() => void exportStage(selectedStage.stageId)}>Export with status</button>}
            </div>
            {selectedExportNotice && <p className="export-transfer-notice" role="status">{selectedExportNotice}</p>}
            {selectedStage.assessment?.qualification.toLowerCase() === 'indeterminate' && <div className="execution-decision-alert" role="status">
              <strong>{selectedStage.assessment.currentlyApplicable ? 'Assessment is indeterminate.' : 'Historical assessment was indeterminate.'}</strong>
              <span>{selectedStage.assessment.reason}</span>
            </div>}
            {selectedStage.assessment?.qualification.toLowerCase() === 'notqualified' && <div className="execution-decision-alert danger" role="status">
              <strong>{selectedStage.assessment.currentlyApplicable ? 'Stage is not qualified.' : 'Historical stage assessment was not qualified.'}</strong>
              <span>{selectedStage.assessment.reason}</span>
            </div>}
            </>}
          </> : <>
            <div className="decision-standing"><strong>{minimizationRunning ? 'Minimization in progress' : awaitingMinimization ? 'Constructed candidate ready' : activeAttempt ? 'Construction in progress' : 'Preparation attempt'}</strong>
              <span>{currentAttemptStage ? 'Completed stage available for inspection' : awaitingMinimization ? 'Review the actual system before required minimization' : activeAttempt ? 'No completed stage yet' : 'No completed stage from this attempt'}</span></div>
            {awaitingMinimization && state.attempt?.constructed && <div className="execution-candidate-brief" aria-label="Actual constructed system before minimization">
              <span>Upper / lower: {(['upper', 'lower'] as const).map(side =>
                state.attempt!.constructed!.achievedComposition.filter(item => item.physicalSide === side)
                  .map(item => `${item.speciesId} ${item.count.toLocaleString()}`).join(', ') || 'None').join(' / ')}</span>
              <span>Cell: {state.attempt.constructed.cellAngstrom.length === 3
                ? `${state.attempt.constructed.cellAngstrom.map(value => value.toFixed(1)).join(' × ')} Å` : 'Not established'}</span>
              <span>Water {state.attempt.constructed.waterCount.toLocaleString()} · Na⁺ {state.attempt.constructed.sodiumCount.toLocaleString()} · Cl⁻ {state.attempt.constructed.chlorideCount.toLocaleString()}</span>
            </div>}
            <div className="button-row decision-actions">
              <button className="button" type="button" onClick={() => document.querySelector(awaitingMinimization ? '.execution-basis-account' : '.execution-progress-account, .execution-identity-account')?.scrollIntoView({ block: 'start' })}>{awaitingMinimization ? 'Inspect counts and evidence' : 'Inspect'}</button>
              {awaitingMinimization && state.attempt?.constructed && <ActionButton state={state} kind="continueMinimization"
                subjectId={state.attempt.attemptId} busy={busy} variant="primary"
                onClick={() => void continueMinimization(state.attempt!)}>Continue minimization</ActionButton>}
              {currentAttemptStage && <button className="button primary" type="button" disabled={busy} onClick={() => void inspectSubject(currentAttemptStage.stageId)}>Inspect completed stage</button>}
              {state.attempt && action(state, 'stopAttempt')?.enabled && <button className="button danger" type="button" disabled={busy}
                onClick={() => void command('stopAttempt', { attemptId: state.attempt!.attemptId })}>{awaitingMinimization ? 'Decline candidate' : 'Stop unfinished work'}</button>}
            </div>
            {awaitingMinimization && state.attempt && action(state, 'continueMinimization', state.attempt.attemptId)?.reason &&
              <p className="action-reason">{action(state, 'continueMinimization', state.attempt.attemptId)?.reason}</p>}
          </>}
        </div> : selectedMembrane && !workflowOpen ? <div className="proposal-decision membrane-decision" aria-label="Membrane model decision">
          <div className="decision-standing">
            <strong>{selectedMembrane.status === 'proposed' ? 'Choice awaits adoption' : selectedMembrane.status === 'assessed' ? 'Assessed membrane model' : 'Membrane model not established'}</strong>
            <span>{selectedMembrane.status === 'proposed' ? 'Proposed intention only' : selectedMembrane.status === 'assessed' ? 'Membrane-local support established' : 'Support not established'}</span>
          </div>
          {selectedMembrane.reason && <p className="decision-reason">{selectedMembrane.reason}</p>}
          <div className="button-row decision-actions">
            {selectedMembrane.status === 'proposed' && <ActionButton state={state} kind="adoptMembrane" busy={busy} variant="primary" onClick={() => void command('adoptMembrane', { modelId: selectedMembrane.modelId })}>Adopt and assess model</ActionButton>}
            <button className="button" type="button" onClick={() => setWorkflowOpen(true)}>{selectedMembrane.status === 'proposed' ? 'Revise proposal' : 'Revise membrane choice'}</button>
          </div>
        </div> : selectedPlacement && !workflowOpen ? <div className="proposal-decision placement-decision" aria-label="Placement decision">
          <div className="decision-standing">
            <strong>{selectedPlacementAdopted ? 'Adopted placement' : selectedPlacement.status === 'supported' ? 'Supported placement' :
              selectedPlacement.status === 'unsupported' ? 'Unsupported placement' : 'Support not established'}</strong>
            <span>{selectedPlacementAdopted ? 'Current study revision' :
              selectedPlacement.status === 'supported' ? 'Adoption available for this exact proposal' : 'Placement cannot be adopted'}</span>
          </div>
          <p className="decision-reason">{selectedPlacement.reason}</p>
          <div className="button-row decision-actions">
            {selectedPlacement.status === 'supported' && !selectedPlacementAdopted && <ActionButton state={state} kind="adoptPlacement" busy={busy} variant="primary"
              onClick={() => void command('adoptPlacement', { proposalId: selectedPlacement.proposalId })}>Adopt supported placement</ActionButton>}
            <button className="button" type="button" onClick={() =>
              document.querySelector('.placement-support-account')?.scrollIntoView({ block: 'start' })}>Inspect evidence</button>
            <button className="button" type="button" onClick={() => showWorkArea('placement')}>Revise proposal</button>
          </div>
        </div> : undefined} />
    </div>

    <nav className="stage-strip" aria-label="System stages and attempts">
      <div className={`stage-summary${selectedStage ? '' : ' single'}`}>
        <div><span className="stage-strip-heading">System stages</span><strong>{state.stages.length === 0 ? 'None yet' : completedStageSummary}</strong>{state.stages.length === 0 && <small>No minimized or equilibrated system stage has completed.</small>}
          {state.attempt && <small title={`Attempt ${state.attempt.attemptId}`}>Current attempt · {readable(state.attempt.stageKind ?? 'preparation')} {readable(state.attempt.status)}</small>}
        </div>
        {selectedStage && <div><span className="stage-strip-heading">Stage assessment</span><strong>{selectedStage.assessment ? `${selectedStage.assessment.currentlyApplicable ? '' : 'Historical '}${readable(selectedStage.assessment.qualification)}` : 'Not established'}</strong></div>}
      </div>
      {state.stages.map(stage => <button type="button" title={`Stage ${stage.stageId} · Attempt ${stage.attemptId}`} className={`stage-card ${state.inspection?.subjectId === stage.stageId ? 'selected' : ''}`} key={stage.stageId} disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true} onClick={() => { showWorkArea('results'); void inspectSubject(stage.stageId); }}>
        <span className="stage-title">{readable(stage.kind)}</span>
        <span className="stage-detail">{stage.studyRevisionId === state.study?.id ? 'Current study' : 'Historical study'}</span>
        <span className="stage-status">{readable(stage.status)}{stage.assessment && ` · ${stage.assessment.currentlyApplicable ? readable(stage.assessment.qualification) : 'Assessment not current'}`}</span>
      </button>)}
    </nav>
  </main>;
}
