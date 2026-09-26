import { useEffect, useRef, useState, type ReactNode } from 'react';
import { ConnectedStructuralInspection, type InspectionAccount } from './ConnectedStructuralInspection';

type SourceRouteKind = 'rcsb' | 'alphafold' | 'upload';
type UploadOriginKind = 'predicted' | 'experimental' | 'unknown';
type PreparationChangeKind = 'alternateLocation' | 'residueState' | 'heavyAtom' | 'disulfide';
type ProteinTopologyKind = 'membrane-spanning' | 'one-surface-associated';
type LeafletSide = 'upper' | 'lower';
type PlacementPhysicalSide = LeafletSide | 'both';
type PpmNterminalSide = 'in' | 'out';
type SourceResidueKind = 'protein' | 'solvent' | 'heterogen';
type PredictionObservationStanding = 'Observed' | 'Unavailable' | 'Unmapped';
type ObservationStanding = 'Observed' | 'Unavailable';
type GeometryKindStanding = ObservationStanding | 'NotApplicable';
type StageKind = 'Minimization' | 'Equilibration';
export type ActorActionKind =
  | 'searchSource' | 'selectSource' | 'selectProteinModel'
  | 'approvePreparationChange' | 'declinePreparationChange'
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
interface PartnerSelection { sourceId: string; retain: boolean; reason: string; }
interface ResidueAddress { model: number; chain: string; residue: number; insertionCode: string; copyId: string; }
interface AlternateLocationChoice { residue: ResidueAddress; altloc: string; decisionId: string; }
interface SourceResidueObservation { address: ResidueAddress; name: string; residueKind: SourceResidueKind; backboneHeavyAtomsComplete: boolean; alternateLocations: string[]; }
interface SourceModelObservation {
  index: number;
  chains: { name: string; residueCount: number; atomCount: number }[];
  assemblies: { name: string; chainCopies: ChainSelection[] }[];
  partners: { sourceId: string; label: string; kind: string; atomCount: number; chain: string | null; residue: number | null }[];
  residues: SourceResidueObservation[];
  atomCount: number;
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
  residue: { model: number; chain: string; residue: number; insertionCode: string; copyId: string };
}

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
  scientificPurpose: string;
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

function readable(value: string | null | undefined): string {
  if (!value) return 'Not established';
  if (value === 'readyForMinimization') return 'Ready for minimization';
  return value.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/[-_]/g, ' ');
}

const preparationChangeLabel: Record<PreparationChangeKind, string> = {
  alternateLocation: 'Alternate location',
  residueState: 'Residue state',
  heavyAtom: 'Heavy atom',
  disulfide: 'Disulfide',
};

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
    <strong>{label} · {readable(geometry.standing)}</strong>
    {geometry.kinds.map(item => <div key={item.kind} className="help-text">
      {readable(item.kind)}: {readable(item.standing)} · {item.measuredCount}/{item.eligibleCount} eligible measured
      {item.minimumDistanceAngstrom !== null && <> · shortest {item.minimumDistanceAngstrom.toFixed(2)} Å</>}
      {item.maximumDistanceAngstrom !== null && <> · longest {item.maximumDistanceAngstrom.toFixed(2)} Å</>}
      {item.unavailableReason && <> · {item.unavailableReason}</>}
    </div>)}
    {geometry.limitations.length > 0 && <p className="help-text">{geometry.limitations.join('; ')}</p>}
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
    {rows.map((row, index) => <div key={index}>
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
      {!row.manual && options.find(lipid => lipid.speciesId === row.speciesId)?.limitations.map(limit => <p className="help-text" key={limit}>{limit}</p>)}
    </div>)}
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

function toFractions(rows: EditableFraction[]): LipidFraction[] {
  return rows.filter(row => Number(row.percent) > 0)
    .map(row => ({ speciesId: row.speciesId.trim(), fraction: Number(row.percent) / 100 }));
}

function selectedSourceModel(state: WorkspaceState, choice: string): SourceModelObservation | undefined {
  if (choice === '') return undefined;
  return state.sourceModels.find(model => model.index === Number(choice));
}

function selectedAssembly(model: SourceModelObservation | undefined, name: string) {
  return model?.assemblies.find(item => item.name === name);
}

export function ProteinInMembraneWorkspace() {
  const [state, setState] = useState<WorkspaceState | null>(null);
  const latestState = useRef<WorkspaceState | null>(null);
  const [communication, setCommunication] = useState<string | null>(null);
  const [streamInterrupted, setStreamInterrupted] = useState(false);
  const [feedback, setFeedback] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [exportBusy, setExportBusy] = useState(false);
  const [localExportFault, setLocalExportFault] = useState<{ stageId: string; reason: string } | null>(null);
  const [localExportNotice, setLocalExportNotice] = useState<{ stageId: string; message: string } | null>(null);
  const [requestedStageTab, setRequestedStageTab] = useState<{
    stageId: string; tab: 'stage' | 'export';
  } | null>(null);
  const [workflowOpen, setWorkflowOpen] = useState(true);
  const [attemptReviewRequested, setAttemptReviewRequested] = useState(false);
  const busyRef = useRef(false);
  const refreshIndex = useRef(0);
  const initialAccountPresented = useRef(false);
  const automaticallyInspectedConstruction = useRef<string | null>(null);

  const [query, setQuery] = useState('');
  const [exactSourceKind, setExactSourceKind] = useState<'rcsb' | 'alphafold'>('rcsb');
  const [exactIdentifier, setExactIdentifier] = useState('');
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [uploadProvenance, setUploadProvenance] = useState<UploadOriginKind>('unknown');
  const [uploadProvenanceNote, setUploadProvenanceNote] = useState('');
  const uploadInput = useRef<HTMLInputElement>(null);
  const [modelChoice, setModelChoice] = useState('');
  const [assemblyChoice, setAssemblyChoice] = useState('');
  const [chainChoices, setChainChoices] = useState<string[]>([]);
  const [partnerChoices, setPartnerChoices] = useState<Record<string, PartnerSelection>>({});
  const [altlocChoices, setAltlocChoices] = useState<Record<string, number>>({});
  const [decisionRationales, setDecisionRationales] = useState<Record<string, string>>({});
  const [upper, setUpper] = useState<EditableFraction[]>([blankFraction()]);
  const [lower, setLower] = useState<EditableFraction[]>([blankFraction()]);
  const [membranePurpose, setMembranePurpose] = useState('');
  const [topologyKind, setTopologyKind] = useState<ProteinTopologyKind | ''>('');
  const [physicalSide, setPhysicalSide] = useState<PlacementPhysicalSide | ''>('');
  const [ppmNterminalSide, setPpmNterminalSide] = useState<PpmNterminalSide | ''>('');
  const [biologicalSidedness, setBiologicalSidedness] = useState('');
  const [depthShift, setDepthShift] = useState('0');
  const [tiltX, setTiltX] = useState('0');
  const [tiltY, setTiltY] = useState('0');
  const [rotationNormal, setRotationNormal] = useState('0');
  const [placementRationale, setPlacementRationale] = useState('');

  async function refresh() {
    const index = ++refreshIndex.current;
    try {
      const response = await fetch('/api/state', { cache: 'no-store' });
      if (!response.ok) throw new Error(`The local workspace could not be read (${response.status}).`);
      const account = await response.json() as WorkspaceState;
      if (index === refreshIndex.current) {
        if (!latestState.current || latestState.current.revision <= account.revision) {
          latestState.current = account;
          setState(account);
          if (!initialAccountPresented.current) {
            initialAccountPresented.current = true;
            if (account.attempt && (!account.inspection ||
                account.inspection.subjectId === account.attempt.constructed?.subjectId ||
                account.stages.some(stage => stage.stageId === account.inspection?.subjectId)))
              setWorkflowOpen(false);
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
    setModelChoice('');
    setAssemblyChoice('');
    setChainChoices([]);
    setPartnerChoices({});
    setAltlocChoices({});
  }, [state?.study?.selectedSourceId]);

  useEffect(() => {
    setDepthShift('0');
    setTiltX('0');
    setTiltY('0');
    setRotationNormal('0');
    setPlacementRationale('');
  }, [state?.placement?.proposalId]);

  async function command(kind: ActorCommandKind, data: object,
                         onFailure?: (reason: string) => void): Promise<boolean> {
    const current = latestState.current;
    if (!current || busyRef.current) return false;
    busyRef.current = true;
    setBusy(true);
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
        setFeedback(reason);
        onFailure?.(reason);
        await refresh();
        return false;
      }
      const updated = result as WorkspaceState;
      if (!latestState.current || latestState.current.revision <= updated.revision) {
        latestState.current = updated;
        setState(updated);
      }
      setCommunication(null);
      setFeedback(null);
      return true;
    } catch (cause) {
      const reason = cause instanceof Error ? cause.message : 'The local request did not complete.';
      setFeedback(reason);
      onFailure?.(reason);
      await refresh();
      return false;
    } finally {
      busyRef.current = false;
      setBusy(false);
    }
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

  async function inspectSubject(subjectId: string) {
    if (await command('selectInspectionSubject', { subjectId })) {
      setAttemptReviewRequested(false);
      setWorkflowOpen(false);
    }
  }

  async function startPreparation() {
    if (await command('startPreparation', {})) {
      setAttemptReviewRequested(true);
      setWorkflowOpen(false);
    }
  }

  async function continueMinimization(attempt: AttemptAccount) {
    if (!attempt.constructed) return;
    if (await command('continueMinimization', {
      attemptId: attempt.attemptId,
      constructedSubjectId: attempt.constructed.subjectId,
    })) {
      setAttemptReviewRequested(true);
      setWorkflowOpen(false);
    }
  }

  function showWorkflowSection(id: string) {
    setWorkflowOpen(true);
    window.requestAnimationFrame(() => document.getElementById(id)?.scrollIntoView({ block: 'start' }));
  }

  async function uploadSource() {
    if (!state || !uploadFile || busyRef.current || action(state, 'selectSource')?.enabled !== true) return;
    busyRef.current = true;
    setBusy(true);
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
        setFeedback(reason);
        return;
      }
      if (typeof result === 'object' && result !== null && 'uploadToken' in result && typeof result.uploadToken === 'string') {
        token = result.uploadToken;
      }
      if (!token) throw new Error('The local host did not identify the uploaded source.');
      setUploadFile(null);
      if (uploadInput.current) uploadInput.current.value = '';
    } catch (cause) {
      setFeedback(cause instanceof Error ? cause.message : 'The source upload did not complete.');
      return;
    } finally {
      busyRef.current = false;
      setBusy(false);
    }
    await command('selectSource', { uploadToken: token });
  }

  async function exportStage(stageId: string) {
    const current = latestState.current;
    if (exportBusy || !current || current.inspection?.subjectId !== stageId ||
        action(current, 'exportStage', stageId)?.enabled !== true) return;
    setExportBusy(true);
    setLocalExportNotice(null);
    try {
      if (!await command('exportStage', { stageId }, reason =>
          setLocalExportFault({ stageId, reason }))) return;
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
    } catch (cause) {
      setLocalExportFault({ stageId, reason: cause instanceof Error
        ? cause.message : 'The verified bundle could not be transferred.' });
      await refresh();
    } finally {
      setExportBusy(false);
    }
  }

  if (!state) return <main className="workspace"><header className="workspace-header"><div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation</div></div></div></header><div className="scene-empty" role="status"><strong>Connecting to the local workspace</strong><span>{communication ?? (streamInterrupted ? 'The live connection is interrupted. Reading established local standing…' : 'Reading the current scientific account…')}</span></div></main>;

  const model = selectedSourceModel(state, modelChoice);
  const assembly = selectedAssembly(model, assemblyChoice);
  const availableChains = assembly?.chainCopies ?? model?.chains.map(chain => ({ sourceChain: chain.name, copyId: chain.name })) ?? [];
  const chosenChains = availableChains.filter(chain => chainChoices.includes(`${chain.sourceChain}:${chain.copyId}`));
  const allPartnersDecided = model?.partners.every(partner => !!partnerChoices[partner.sourceId] && !!partnerChoices[partner.sourceId].reason.trim()) ?? false;
  const chosenPartners = model?.partners.map(partner => partnerChoices[partner.sourceId]).filter((item): item is PartnerSelection => !!item) ?? [];
  const selectedSourceChains = new Set(chosenChains.map(chain => chain.sourceChain));
  const ambiguousResidues = model?.residues.filter(residue =>
    residue.residueKind === 'protein' && selectedSourceChains.has(residue.address.chain) &&
    residue.alternateLocations.length > 0) ?? [];
  const residueKey = (address: ResidueAddress) => `${address.model}:${address.chain}:${address.residue}:${address.insertionCode}:${address.copyId}`;
  const allAltlocsChosen = ambiguousResidues.every(residue => Number.isInteger(altlocChoices[residueKey(residue.address)]));
  const chosenAltlocs = ambiguousResidues.map(residue => ({ residue: residue.address, altloc: residue.alternateLocations[altlocChoices[residueKey(residue.address)]] }));
  const membraneReady = validFractions(upper) && validFractions(lower) && membranePurpose.trim().length > 0;
  const correctionValues = [depthShift, tiltX, tiltY, rotationNormal].map(Number);
  const correctionReady = state.placement !== null && placementRationale.trim().length > 0
    && correctionValues.every(Number.isFinite) && correctionValues.some(value => value !== 0);
  const selectedStage = state.stages.find(stage => stage.stageId === state.inspection?.subjectId);
  const selectedExportFault = selectedStage ? (localExportFault?.stageId === selectedStage.stageId
    ? localExportFault.reason : selectedStage.export?.status === 'failed'
      ? selectedStage.export.reason ?? 'The verified bundle could not be delivered.' : null) : null;
  const selectedExportNotice = selectedStage && localExportNotice?.stageId === selectedStage.stageId
    ? localExportNotice.message : null;
  const selectedSourceStage = selectedStage?.kind === 'Equilibration' && selectedStage.sourceStageId
    ? state.stages.find(stage => stage.stageId === selectedStage.sourceStageId &&
      stage.attemptId === selectedStage.attemptId && stage.kind === 'Minimization') : undefined;
  const pairedMinimizedStage = selectedStage?.kind === 'Minimization'
    ? selectedStage : selectedSourceStage;
  const equilibratedDescendants = selectedStage?.kind === 'Minimization'
    ? state.stages.filter(stage => stage.kind === 'Equilibration' &&
      stage.sourceStageId === selectedStage.stageId &&
      stage.attemptId === selectedStage.attemptId) : [];
  const pairedEquilibratedStage = selectedStage?.kind === 'Equilibration'
    ? selectedStage : equilibratedDescendants.length === 1 ? equilibratedDescendants[0] : undefined;
  const completedStageSummary = selectedStage
    ? (['Minimization', 'Equilibration'] as const)
      .filter(kind => state.stages.some(stage => stage.attemptId === selectedStage.attemptId &&
        stage.kind === kind && stage.status.toLowerCase() === 'completed'))
      .map(kind => kind === 'Minimization' ? 'Minimized' : 'Equilibrated').join(' · ')
    : state.stages.length === 0 ? 'None' : String(state.stages.length);
  const stageViewTab = selectedStage && selectedExportFault &&
    requestedStageTab?.stageId === selectedStage.stageId
    ? requestedStageTab.tab : selectedExportFault ? 'export' : 'stage';
  const selectedChange = state.protein?.changes.find(change => change.id === state.inspection?.subjectId);
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
  const activeContext = state.placement ? 'placement' : state.membrane ? 'membrane' : 'protein';

  return <main className={`workspace ${workflowOpen ? 'workflow-open' : 'workflow-closed'}${selectedMembrane ? ' membrane-review' : ''}${selectedPlacement && !executionReview ? ' placement-review' : ''}${executionReview ? ' execution-review' : ''}`}>
    <header className="workspace-header">
      <div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation · one planar bilayer</div></div></div>
      <nav className="workspace-context-tabs" aria-label="Jump to researcher workflow">
        {executionReview ? <>
          <button type="button" onClick={() => showWorkflowSection('protein-workflow')}>Study revision</button>
          {selectedStage && selectedExportFault && <button type="button"
            className={stageViewTab === 'export' ? 'active' : ''}
            aria-current={stageViewTab === 'export' ? 'page' : undefined}
            onClick={() => {
              setWorkflowOpen(false);
              setRequestedStageTab({ stageId: selectedStage.stageId, tab: 'export' });
              document.querySelector('.export-validation')?.scrollIntoView({ block: 'start' });
            }}>Export validation</button>}
          {pairedMinimizedStage && pairedEquilibratedStage ? <>
            <button type="button" className={selectedStage?.stageId === pairedMinimizedStage.stageId && stageViewTab === 'stage' ? 'active' : ''}
              aria-current={selectedStage?.stageId === pairedMinimizedStage.stageId && stageViewTab === 'stage' ? 'page' : undefined}
              onClick={() => { setWorkflowOpen(false); void inspectSubject(pairedMinimizedStage.stageId); }}>
              Minimized stage review</button>
            <button type="button" className={selectedStage?.stageId === pairedEquilibratedStage.stageId && stageViewTab === 'stage' ? 'active' : ''}
              aria-current={selectedStage?.stageId === pairedEquilibratedStage.stageId && stageViewTab === 'stage' ? 'page' : undefined}
              onClick={() => { setWorkflowOpen(false); void inspectSubject(pairedEquilibratedStage.stageId); }}>
              Equilibrated stage review</button>
          </> : <button type="button" className={stageViewTab === 'stage' ? 'active' : ''}
            aria-current={stageViewTab === 'stage' ? 'page' : undefined}
            onClick={() => {
              setWorkflowOpen(false);
              if (selectedStage) {
                setRequestedStageTab({ stageId: selectedStage.stageId, tab: 'stage' });
                document.querySelector('.execution-assessment-account')?.scrollIntoView({ block: 'start' });
              } else if (state.attempt) setAttemptReviewRequested(true);
            }}>
            {selectedStage ? `${selectedStage.kind === 'Minimization' ? 'Minimized' : 'Equilibrated'} stage review` : minimizationRunning ? 'Minimization running' : awaitingMinimization ? 'Construction review' : activeAttempt ? 'Preparation running' : 'Preparation attempt review'}
          </button>}
        </> : <>
        <button type="button" className={activeContext === 'protein' ? 'active' : ''} aria-current={activeContext === 'protein' ? 'page' : undefined} onClick={() => showWorkflowSection('protein-workflow')}>Protein preparation</button>
        <button type="button" className={activeContext === 'membrane' ? 'active' : ''} aria-current={activeContext === 'membrane' ? 'page' : undefined} onClick={() => showWorkflowSection('membrane-workflow')}>Study revision</button>
        {state.attempt && <button type="button" onClick={() => { setAttemptReviewRequested(true); setWorkflowOpen(false); }}>
          {state.attempt.stageKind === 'Minimization' && activeAttempt ? 'Minimization running' : 'Current preparation attempt'}
        </button>}
        <button type="button" className={activeContext === 'placement' ? 'active' : ''} aria-current={activeContext === 'placement' ? 'page' : undefined} onClick={() => showWorkflowSection('placement-workflow')}>Placement review</button>
        </>}
      </nav>
      <div className="header-context">
        <button className="button workflow-toggle" type="button" aria-expanded={workflowOpen} aria-controls="researcher-workflow"
          aria-label={workflowOpen ? 'Hide workflow' : 'Show workflow'}
          title={workflowOpen ? 'Hide researcher workflow' : 'Show researcher workflow'}
          onClick={() => setWorkflowOpen(value => !value)}>Workspace <span aria-hidden="true">⌄</span></button>
        <span className="context-chip">Study <strong className="id">{state.study?.id ?? 'Not established'}</strong></span>
        <span className="context-chip">Study revision <strong className="id">{state.study?.number ?? 'Not established'}</strong></span>
        {state.study?.selectedSourceId && <span className="context-chip" title={state.study.selectedSourceId}>Source <strong className="id">{state.study.selectedSourceId}</strong></span>}
        {state.study?.modelIndex !== null && state.study?.modelIndex !== undefined && <span className="context-chip compact-model-context" title={`Model ${state.study.modelIndex}; ${state.study.biologicalAssemblyId ? `assembly ${state.study.biologicalAssemblyId}` : 'no assembly transformation'}; chains ${state.study.chainIds.join(', ')}`}>
          <strong>Model {state.study.modelIndex} · {state.study.biologicalAssemblyId ? `assembly ${state.study.biologicalAssemblyId}` : 'no assembly transform'} · {state.study.chainIds.join(', ')}</strong>
        </span>}
        {state.attempt && <span className="context-chip">Attempt <strong className="id">{state.attempt.attemptId}</strong></span>}
        {selectedStage && <span className="context-chip">Stage <strong className="id">{selectedStage.stageId}</strong></span>}
      </div>
    </header>

    <div className="workspace-body">
      <aside id="researcher-workflow" className="rail" aria-label="Researcher choices and available actions" hidden={!workflowOpen}>
        <p className="eyebrow">Researcher workflow</p>
        <h2 className="panel-heading">Prepare the system</h2>
        {feedback && <div className="notice warning" role="alert">{feedback}</div>}
        {state.notices.length > 0 && <div className="notice-list" aria-label="Scientific and workflow notices">
          {state.notices.map(notice => <div key={notice.id} className={`notice ${notice.severity.toLowerCase()}`}><strong>{readable(notice.severity)}</strong> · {notice.message}{notice.subjectId && <small className="tabular">Subject {notice.subjectId}</small>}</div>)}
        </div>}

        <section className="rail-section" id="protein-workflow">
          <h3 className="section-title">1 · Structural source</h3>
          <label className="field-label" htmlFor="source-query">Discover structural sources</label>
          <input id="source-query" className="text-input" value={query} onChange={event => setQuery(event.target.value)} placeholder="Protein name or accession" />
          <div className="button-row"><ActionButton state={state} kind="searchSource" busy={busy} onClick={() => void command('searchSource', { query: query.trim() })}>Find sources</ActionButton></div>
          {action(state, 'searchSource')?.reason && <p className="action-reason">{action(state, 'searchSource')?.reason}</p>}
          {state.sourceCandidates.length > 0 && <div className="source-list">
            {state.sourceCandidates.map(candidate => <button type="button" className={`source-item ${state.study?.selectedSourceId === candidate.id ? 'selected' : ''}`} key={candidate.id} disabled={busy || action(state, 'selectSource')?.enabled !== true} onClick={() => void command('selectSource', { sourceId: candidate.id })}>
              <span className="item-title">{candidate.label}</span>
              <span className="item-detail">{candidate.kind} · {candidate.provenance}</span>
              {candidate.limitations.length > 0 && <span className="item-detail">{candidate.limitations.join('; ')}</span>}
            </button>)}
          </div>}
          <fieldset className="source-route">
            <legend>Or load an exact database reference</legend>
            <p className="help-text">An exact reference can be used when candidate search is unavailable. The coordinates and model still need inspection.</p>
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
              onClick={() => void command('selectSource', { sourceKind: exactSourceKind, exactIdentifier: exactIdentifier.trim() })}>Inspect exact source</ActionButton></div>
          </fieldset>
          <label className="field-label" htmlFor="source-upload">Or supply a PDB/mmCIF file</label>
          <input id="source-upload" ref={uploadInput} className="text-input" type="file" accept=".pdb,.cif,.mmcif" onChange={event => setUploadFile(event.target.files?.[0] ?? null)} />
          <label className="field-label" htmlFor="upload-provenance">Researcher-declared structure origin</label>
          <select id="upload-provenance" className="select-input" value={uploadProvenance} onChange={event => {
            const value = event.target.value;
            if (value === 'unknown' || value === 'experimental' || value === 'predicted') setUploadProvenance(value);
          }}>
            <option value="unknown">Unknown or not established</option>
            <option value="experimental">Experimental coordinates</option>
            <option value="predicted">Predicted model</option>
          </select>
          <label className="field-label" htmlFor="upload-provenance-note">Source or method note, if available</label>
          <input id="upload-provenance-note" className="text-input" maxLength={500} value={uploadProvenanceNote} onChange={event => setUploadProvenanceNote(event.target.value)} placeholder="Archive, experiment, or prediction method" />
          <p className="help-text">This declaration identifies the source; it does not verify its method or turn uploaded B-factors into prediction confidence. A predicted or unknown upload needs an attributable qualification route before placement can be supported.</p>
          <div className="button-row"><ActionButton state={state} kind="selectSource" busy={busy || !uploadFile} onClick={() => void uploadSource()}>Inspect uploaded source</ActionButton></div>
          {action(state, 'selectSource')?.reason && <p className="action-reason">{action(state, 'selectSource')?.reason}</p>}
          {state.study?.selectedSourceId && <><div className="source-context">
            <strong>Selected structural source</strong>
            <span className="tabular">{state.study.selectedSourceId}</span>
            {state.study.selectedSourceKind === 'rcsb' && <span>Retrieved RCSB PDB coordinate entry</span>}
            {state.study.selectedSourceKind === 'alphafold' && <span>AlphaFold DB predicted coordinate record</span>}
            {state.study.uploadProvenance && <span>Researcher-declared {readable(state.study.uploadProvenance)} origin</span>}
            {state.study.uploadProvenanceNote && <span>{state.study.uploadProvenanceNote}</span>}
            <span className="help-text">Source coordinates remain distinct from any prepared protein.</span>
          </div>
            <div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy}
              onClick={() => void inspectSubject(state.study!.selectedSourceId!)}>Inspect source structure</ActionButton></div></>}
          <PredictionSummary prediction={state.sourcePrediction} label="Selected source" />
        </section>

        <section className="rail-section">
          <h3 className="section-title">2 · Protein model</h3>
          {state.study?.modelIndex !== null && state.study?.modelIndex !== undefined && <div className="hint-box">
            <strong>Current adopted model {state.study.modelIndex}</strong><br />
            Assembly {state.study.biologicalAssemblyId ?? 'none'} · chains {state.study.chainIds.join(', ') || 'not identified'}
            {state.study.partners.length > 0 && <><br />Partners: {state.study.partners.map(partner => `${partner.sourceId} ${partner.retain ? 'retained' : 'excluded'}`).join('; ')}</>}
            {state.study.alternateLocations.length > 0 && <><br />Alternate locations: {state.study.alternateLocations.map(choice => `${choice.residue.chain}${choice.residue.residue}${choice.residue.insertionCode} → ${choice.altloc || '(unlabelled)'}`).join('; ')}</>}
          </div>}
          {state.sourceModels.length === 0 ? <p className="help-text">Inspect a structural source to see its observed models, assemblies, chains and partners.</p> : <>
            <label className="field-label" htmlFor="model-index">Coordinate model</label>
            <select id="model-index" className="select-input" value={modelChoice} onChange={event => { setModelChoice(event.target.value); setAssemblyChoice(''); setChainChoices([]); setPartnerChoices({}); setAltlocChoices({}); }}>
              <option value="">Choose an observed model</option>
              {state.sourceModels.map(item => <option key={item.index} value={item.index}>Model {item.index} · {item.atomCount.toLocaleString()} atoms</option>)}
            </select>
            {model && <>
              <label className="field-label" htmlFor="assembly-choice">Biological assembly</label>
              <select id="assembly-choice" className="select-input" value={assemblyChoice} onChange={event => { setAssemblyChoice(event.target.value); setChainChoices([]); }}>
                <option value="">No assembly transformation</option>
                {model.assemblies.map(item => <option key={item.name} value={item.name}>{item.name}</option>)}
              </select>
              <span className="field-label">Retained chain copies</span>
              {availableChains.map(chain => {
                const key = `${chain.sourceChain}:${chain.copyId}`;
                return <label className="checkbox-line" key={key}><input type="checkbox" checked={chainChoices.includes(key)} onChange={event => setChainChoices(previous => event.target.checked ? [...previous, key] : previous.filter(value => value !== key))} /><span>{chain.sourceChain} · copy {chain.copyId}</span></label>;
              })}
              {model.partners.length > 0 && <>
                <span className="field-label">Bound partners · decide each explicitly</span>
                {assemblyChoice && <p className="supporting-text">A partner decision applies to every corresponding selected assembly copy. Different occupancy between copies requires a source structure that identifies those partners separately.</p>}
                {model.partners.map(partner => <div className="partner-choice" key={partner.sourceId}>
                  <strong>{partner.label}</strong><small>{partner.kind} · {partner.atomCount} atoms</small>
                  <div className="inline-fields">
                    <label className="checkbox-line"><input type="radio" name={`partner-${partner.sourceId}`} checked={partnerChoices[partner.sourceId]?.retain === true} onChange={() => setPartnerChoices(previous => ({ ...previous, [partner.sourceId]: { sourceId: partner.sourceId, retain: true, reason: previous[partner.sourceId]?.reason ?? '' } }))} />Retain</label>
                    <label className="checkbox-line"><input type="radio" name={`partner-${partner.sourceId}`} checked={partnerChoices[partner.sourceId]?.retain === false} onChange={() => setPartnerChoices(previous => ({ ...previous, [partner.sourceId]: { sourceId: partner.sourceId, retain: false, reason: previous[partner.sourceId]?.reason ?? '' } }))} />Exclude</label>
                  </div>
                  {partnerChoices[partner.sourceId] && <input className="text-input" aria-label={`Reason for ${partner.label}`} placeholder="Reason for this decision" value={partnerChoices[partner.sourceId].reason} onChange={event => setPartnerChoices(previous => ({ ...previous, [partner.sourceId]: { ...previous[partner.sourceId], reason: event.target.value } }))} />}
                </div>)}
              </>}
              {model.residues.some(residue => !residue.backboneHeavyAtomsComplete) && <div className="notice warning">The selected model contains observed residues with incomplete backbone heavy-atom coordinates. Assessment may refuse them; no missing backbone is silently reconstructed.</div>}
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
              <div className="button-row"><ActionButton state={state} kind="selectProteinModel" busy={busy || chosenChains.length === 0 || !allPartnersDecided || !allAltlocsChosen} onClick={() => void command('selectProteinModel', { modelIndex: model.index, biologicalAssemblyId: assemblyChoice || null, chains: chosenChains, partners: chosenPartners, alternateLocations: chosenAltlocs })}>Assess selected protein</ActionButton></div>
              {action(state, 'selectProteinModel')?.reason && <p className="action-reason">{action(state, 'selectProteinModel')?.reason}</p>}
              {chosenChains.length === 0 && <p className="help-text">Choose at least one observed chain copy.</p>}
              {!allAltlocsChosen && <p className="help-text">Choose each observed alternate location before assessment.</p>}
            </>}
          </>}
          {state.protein && <><div className={`protein-standing ${state.protein.status === 'assessed' ? 'established' : 'unestablished'}`} role="status">
            <strong>{state.protein.status === 'assessed' ? 'Assessed prepared protein' : 'Prepared protein not established'}</strong>
            <span>{state.protein.summary}</span>
            <small className="tabular">Study revision {state.study?.number ?? 'unknown'} · subject {state.protein.subjectId}</small>
          </div><div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.protein!.subjectId)}>Inspect protein</ActionButton>
              {state.protein.candidateId && <ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.protein!.candidateId!)}>Inspect unqualified candidate</ActionButton>}
            </div>
            <PredictionSummary prediction={state.protein.prediction} label={state.protein.status === 'assessed' ? 'Assessed protein' : 'Selected protein'} />
            <GeometrySummary geometry={state.protein.sourceGeometry} label="Selected source geometry" />
            <GeometrySummary geometry={state.protein.geometry} label="Prepared protein geometry" />
          </>}
          {state.protein?.changes.map(change => <div className={`change-card ${state.inspection?.subjectId === change.id ? 'selected' : ''}`} key={change.id}>
            <strong>{preparationChangeLabel[change.kind]} · {change.residue.chain}{change.residue.residue}{change.residue.insertionCode} · copy {change.residue.copyId}</strong>
            <small className="tabular">Proposal {change.id} · study {change.studyRevisionId}</small>
            {state.inspection?.subjectId === change.id && <span className="selected-context">Selected proposal · linked structure and evidence shown beside this decision</span>}
            <p>{change.proposedChange}</p><p>{change.rationale}</p>
            {change.limitations.length > 0 && <small>{change.limitations.join('; ')}</small>}
            <div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(change.id)}>Inspect change and evidence</ActionButton></div>
          </div>)}
        </section>

        <section className="rail-section" id="membrane-workflow">
          <h3 className="section-title">3 · Planar membrane</h3>
          <p className="help-text">Specify each leaflet. Percentages express the intended model, not achieved molecule counts.</p>
          {state.availableLipids.length === 0 && <div className="hint-box">No lipid catalogue is currently qualified for selection. Enter an exact species ID below to have its support assessed without substitution.</div>}
          <FractionEditor label="Upper" rows={upper} onChange={setUpper} options={state.availableLipids} />
          <FractionEditor label="Lower" rows={lower} onChange={setLower} options={state.availableLipids} />
          <label className="field-label" htmlFor="scientific-purpose">Scientific purpose</label>
          <input id="scientific-purpose" className="text-input" value={membranePurpose} onChange={event => setMembranePurpose(event.target.value)} placeholder="Why this membrane model?" />
          <div className="button-row"><ActionButton state={state} kind="proposeMembrane" busy={busy || !membraneReady} onClick={() => void command('proposeMembrane', { upper: toFractions(upper), lower: toFractions(lower), scientificPurpose: membranePurpose.trim() })}>Propose membrane model</ActionButton></div>
          {action(state, 'proposeMembrane')?.reason && <p className="action-reason">{action(state, 'proposeMembrane')?.reason}</p>}
          {!membraneReady && <p className="help-text">Each leaflet needs distinct identified positive-fraction species totaling 100%, plus a study purpose. Zero-percent rows are absent from the proposed model; support is assessed after adoption.</p>}
          {state.membrane && <><div className="hint-box"><strong>{readable(state.membrane.status)}</strong><br /><span className="tabular">{state.membrane.modelId}</span><br />Purpose: {state.membrane.scientificPurpose}{state.membrane.limitations.length > 0 && <><br />{state.membrane.limitations.join('; ')}</>}</div><div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.membrane!.modelId)}>Inspect membrane</ActionButton><ActionButton state={state} kind="adoptMembrane" busy={busy} onClick={() => void command('adoptMembrane', { modelId: state.membrane!.modelId })}>Adopt and assess model</ActionButton></div></>}
          {state.study && <p className="help-text tabular">Fixed conditions: pH {state.study.conditions.nominalPh}; NaCl {state.study.conditions.targetNaClMolar} M; temperature {state.study.conditions.optionalTemperatureKelvin} K.</p>}
        </section>

        <section className="rail-section" id="placement-workflow">
          <h3 className="section-title">4 · Protein placement</h3>
          <label className="field-label" htmlFor="topology-kind">Membrane relationship</label>
          <select id="topology-kind" className="select-input" value={topologyKind} onChange={event => {
            const value = event.target.value;
            if (value === '' || value === 'membrane-spanning' || value === 'one-surface-associated') {
              setTopologyKind(value);
              setPhysicalSide(value === 'membrane-spanning' ? 'both' : '');
            }
          }}>
            <option value="">Choose relationship</option>
            <option value="membrane-spanning">Membrane-spanning</option>
            <option value="one-surface-associated">Associated with one surface</option>
          </select>
          {topologyKind === 'one-surface-associated' && <><label className="field-label" htmlFor="physical-side">Contact side</label><select id="physical-side" className="select-input" value={physicalSide} onChange={event => {
            const value = event.target.value;
            if (value === '' || value === 'upper' || value === 'lower') setPhysicalSide(value);
          }}><option value="">Choose physical side</option><option value="upper">Upper</option><option value="lower">Lower</option></select></>}
          <label className="field-label" htmlFor="ppm-nterminal-side">PPM N-terminal topology side</label>
          <select id="ppm-nterminal-side" className="select-input" value={ppmNterminalSide} onChange={event => {
            const value = event.target.value;
            if (value === '' || value === 'in' || value === 'out') setPpmNterminalSide(value);
          }}>
            <option value="">Choose N-terminal side</option>
            <option value="in">In</option>
            <option value="out">Out</option>
          </select>
          <p className="help-text">PPM’s in/out topology is a biological N-terminal assignment; it is not the upper/lower side of the displayed bilayer.</p>
          <label className="field-label" htmlFor="biological-sidedness">Biological sidedness premise, if known</label>
          <input id="biological-sidedness" className="text-input" value={biologicalSidedness} onChange={event => setBiologicalSidedness(event.target.value)} placeholder="Evidence-backed description or leave blank" />
          <p className="help-text">You can inspect a geometric proposal without this premise. Supported placement additionally requires independent, construct-specific topology, contact and sidedness evidence; geometry or this description alone cannot establish it.</p>
          <div className="button-row"><ActionButton state={state} kind="proposePlacement" busy={busy || !topologyKind || !physicalSide || !ppmNterminalSide} onClick={() => void command('proposePlacement', { topologyKind, physicalSide, ppmNterminalSide, biologicalSidedness: biologicalSidedness.trim() || null })}>Assess placement</ActionButton></div>
          {action(state, 'proposePlacement')?.reason && <p className="action-reason">{action(state, 'proposePlacement')?.reason}</p>}
          {state.placement && <><div className="hint-box"><strong>{readable(state.placement.status)}</strong><br />{state.placement.reason}</div><div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.placement!.proposalId)}>Inspect placement</ActionButton><ActionButton state={state} kind="adoptPlacement" busy={busy} onClick={() => void command('adoptPlacement', { proposalId: state.placement!.proposalId })}>Adopt supported placement</ActionButton></div>
            <PlacementPredictionSummary prediction={state.placement.prediction} />
            <div className="change-card">
              <strong>Bounded placement correction</strong>
              <p>Adjust the identified proposal, then reassess it. A geometric change cannot override contradictory scientific evidence.</p>
              <div className="inline-fields">
                <label><span className="field-label">Depth shift (Å)</span><input className="number-input tabular" type="number" step="0.1" value={depthShift} onChange={event => setDepthShift(event.target.value)} /></label>
                <label><span className="field-label">Tilt about X (°)</span><input className="number-input tabular" type="number" step="0.1" value={tiltX} onChange={event => setTiltX(event.target.value)} /></label>
                <label><span className="field-label">Tilt about Y (°)</span><input className="number-input tabular" type="number" step="0.1" value={tiltY} onChange={event => setTiltY(event.target.value)} /></label>
                <label><span className="field-label">Rotation about normal (°)</span><input className="number-input tabular" type="number" step="0.1" value={rotationNormal} onChange={event => setRotationNormal(event.target.value)} /></label>
              </div>
              <label className="field-label" htmlFor="placement-rationale">Scientific rationale</label>
              <input id="placement-rationale" className="text-input" value={placementRationale} onChange={event => setPlacementRationale(event.target.value)} placeholder="Evidence for this adjustment" />
              <div className="button-row"><ActionButton state={state} kind="revisePlacement" busy={busy || !correctionReady} onClick={() => void command('revisePlacement', { proposalId: state.placement!.proposalId, depthShiftAngstrom: Number(depthShift), tiltAboutXDegrees: Number(tiltX), tiltAboutYDegrees: Number(tiltY), rotationAboutNormalDegrees: Number(rotationNormal), rationale: placementRationale.trim() })}>Reassess corrected placement</ActionButton></div>
              {action(state, 'revisePlacement')?.reason && <p className="action-reason">{action(state, 'revisePlacement')?.reason}</p>}
            </div>
          </>}
        </section>

        <section className="rail-section">
          <h3 className="section-title">5 · Prepare and inspect</h3>
          <div className="button-row"><ActionButton state={state} kind="startPreparation" busy={busy} variant="primary" onClick={() => void startPreparation()}>Construct system</ActionButton></div>
          {action(state, 'startPreparation')?.reason && <p className="action-reason">{action(state, 'startPreparation')?.reason}</p>}
          {state.attempt && <div className="attempt-account">
            <strong className="tabular">Attempt {state.attempt.attemptId}</strong>
            <span>{readable(state.attempt.status)}{state.attempt.stageKind && ` · ${readable(state.attempt.stageKind)}`}</span>
            <p>{state.attempt.message}</p>
            {state.attempt.progress !== null && <progress max="1" value={state.attempt.progress} aria-label="Observed attempt progress" />}
            <div className="button-row">
              <button className="button" type="button" onClick={() => { setAttemptReviewRequested(true); setWorkflowOpen(false); }}>Review attempt</button>
              {awaitingMinimization && state.attempt.constructed && <ActionButton state={state} kind="continueMinimization" subjectId={state.attempt.attemptId}
                busy={busy} variant="primary" onClick={() => void continueMinimization(state.attempt!)}>Continue minimization</ActionButton>}
              <ActionButton state={state} kind="stopAttempt" busy={busy} variant="danger" onClick={() => void command('stopAttempt', { attemptId: state.attempt!.attemptId })}>{awaitingMinimization ? 'Decline candidate' : 'Stop unfinished work'}</ActionButton>
            </div>
          </div>}
          {selectedStage && <div className="stage-export">
            <strong>{readable(selectedStage.kind)} · {readable(selectedStage.status)}</strong>
            <p className="tabular">Attempt {selectedStage.attemptId} · study {selectedStage.studyRevisionId}</p>
            <p>{selectedStage.assessment ? `${selectedStage.assessment.currentlyApplicable ? readable(selectedStage.assessment.qualification) : `Historical ${readable(selectedStage.assessment.qualification)}; assessment not current`} — ${selectedStage.assessment.reason}` : 'No scientific assessment established for this stage.'}</p>
            <div className="button-row"><ActionButton state={state} kind="exportStage" subjectId={selectedStage.stageId} busy={busy || exportBusy} onClick={() => void exportStage(selectedStage.stageId)}>Export this completed stage</ActionButton></div>
            {action(state, 'exportStage', selectedStage.stageId)?.reason && <p className="action-reason">{action(state, 'exportStage', selectedStage.stageId)?.reason}</p>}
            {selectedExportFault && <p className="action-reason">Export not delivered: {selectedExportFault}</p>}
            {selectedExportNotice && <p className="help-text" role="status">{selectedExportNotice}</p>}
            {selectedStage.kind === 'Minimization' && <><div className="button-row"><ActionButton state={state} kind="requestEquilibration" subjectId={selectedStage.stageId} busy={busy} onClick={() => void command('requestEquilibration', { stageId: selectedStage.stageId })}>Request optional equilibration from this stage</ActionButton></div>{action(state, 'requestEquilibration', selectedStage.stageId)?.reason && <p className="action-reason">{action(state, 'requestEquilibration', selectedStage.stageId)?.reason}</p>}</>}
          </div>}
        </section>
      </aside>

      <ConnectedStructuralInspection state={state} reviewAttempt={reviewAttempt} onSelectFocus={annotationId => void command('setInspectionFocus', { annotationId })}
        onInspectSubject={subjectId => { void inspectSubject(subjectId); }}
        exportFault={selectedExportFault}
        connectionMessage={communication ?? (streamInterrupted ? 'The live connection is interrupted. Showing the last account read from the host; progress may change until reconnection.' : null)}
        onRefreshAccount={() => void refresh()}
        proposalDecision={selectedChange && !executionReview ? <div className="proposal-decision" aria-label="Preparation change decision">
          <div className="decision-standing">
            <strong>{state.protein?.status === 'declined' ? 'Change declined' : selectedChange.approvalRequired ? 'Approval required' : 'Proposal under review'}</strong>
            <span>{state.protein?.status === 'assessed' ? 'Assessed prepared protein' : 'Prepared protein not established'}</span>
          </div>
          {state.protein?.status === 'declined' && <p className="decision-reason">{state.protein.summary}</p>}
          {(selectedChange.kind === 'residueState' || selectedChange.kind === 'disulfide') && <label className="field-label" htmlFor={`decision-${selectedChange.id}`}>
            Scientific rationale for this {selectedChange.kind === 'disulfide' ? 'disulfide choice' : 'site-specific state'}
            <input id={`decision-${selectedChange.id}`} className="text-input" value={decisionRationales[selectedChange.id] ?? ''}
              onChange={event => setDecisionRationales(current => ({ ...current, [selectedChange.id]: event.target.value }))}
              placeholder="Explain why this is the intended model assumption" />
          </label>}
          {selectedChange.approvalRequired && <div className="button-row decision-actions">
            <ActionButton state={state} kind="approvePreparationChange" subjectId={selectedChange.id} variant="primary"
              busy={busy || ((selectedChange.kind === 'residueState' || selectedChange.kind === 'disulfide') && !(decisionRationales[selectedChange.id] ?? '').trim())}
              onClick={() => void command('approvePreparationChange', { proposalId: selectedChange.id, approve: true, rationale: (decisionRationales[selectedChange.id] ?? '').trim() })}>Approve change</ActionButton>
            <ActionButton state={state} kind="declinePreparationChange" subjectId={selectedChange.id} busy={busy}
              onClick={() => void command('approvePreparationChange', { proposalId: selectedChange.id, approve: false, rationale: (decisionRationales[selectedChange.id] ?? '').trim() })}>Decline</ActionButton>
          </div>}
          {action(state, 'approvePreparationChange', selectedChange.id)?.reason && <p className="action-reason">{action(state, 'approvePreparationChange', selectedChange.id)?.reason}</p>}
        </div> : executionReview && !workflowOpen ? <div className="proposal-decision execution-decision" aria-label={selectedStage ? 'Completed stage next steps' : 'Attempt decision'}>
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
            <button className="button" type="button" onClick={() => showWorkflowSection('placement-workflow')}>Revise proposal</button>
          </div>
        </div> : undefined} />
    </div>

    <nav className="stage-strip" aria-label="Completed and attempted stages">
      <div className="stage-summary">
        <div><span className="stage-strip-heading">Completed stages</span><strong>{completedStageSummary}</strong>{state.stages.length === 0 && <small>No completed stage is currently established.</small>}
          {state.attempt && <small className="tabular">Attempt {state.attempt.attemptId} · {readable(state.attempt.stageKind ?? 'preparation')} {readable(state.attempt.status)} · origin {state.attempt.studyRevisionId ?? 'pending'}</small>}
        </div>
        <div><span className="stage-strip-heading">Stage assessment</span><strong>{selectedStage?.assessment ? `${selectedStage.assessment.currentlyApplicable ? '' : 'Historical '}${readable(selectedStage.assessment.qualification)}` : selectedStage ? 'Not established' : 'Not applicable'}</strong></div>
      </div>
      {state.stages.map(stage => <button type="button" className={`stage-card ${state.inspection?.subjectId === stage.stageId ? 'selected' : ''}`} key={stage.stageId} disabled={busy || action(state, 'selectInspectionSubject')?.enabled !== true} onClick={() => void inspectSubject(stage.stageId)}>
        <span className="stage-title">{readable(stage.kind)}</span>
        <span className="stage-detail tabular">{stage.stageId}</span>
        <span className="stage-detail tabular">Attempt {stage.attemptId}</span>
        <span className="stage-status">{readable(stage.status)}{stage.assessment && ` · ${stage.assessment.currentlyApplicable ? readable(stage.assessment.qualification) : 'Assessment not current'}`}</span>
      </button>)}
    </nav>
  </main>;
}
