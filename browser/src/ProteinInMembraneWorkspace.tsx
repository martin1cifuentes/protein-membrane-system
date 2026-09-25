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
  | 'stopAttempt' | 'requestEquilibration' | 'selectInspectionSubject'
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

export interface ScientificFinding { id: string; meaning: string; consequence: string; disposition: string; material: boolean; }
export interface ScientificEvidence { id: string; source: string; method: string; observation: string; applicability: string; uncertainty: string; bearing: string; }

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
interface MembraneAccount {
  modelId: string;
  status: string;
  upper: LipidFraction[];
  lower: LipidFraction[];
  scientificPurpose: string;
  limitations: string[];
}

interface PlacementAccount {
  proposalId: string;
  status: string;
  topologyKind: ProteinTopologyKind;
  depthAngstrom: number | null;
  tiltDegrees: number | null;
  sidedness: string | null;
  reason: string;
  evidence: ScientificEvidence[];
  prediction: PredictionRegionSummaryObservations | null;
}

interface AttemptAccount {
  attemptId: string;
  status: string;
  stageKind: StageKind | null;
  progress: number | null;
  message: string;
}

interface PreparationAssessmentResult {
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
  status: string;
  assessment: PreparationAssessmentResult | null;
  summary: string;
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

interface EditableFraction { speciesId: string; percent: string; }
const blankFraction = (): EditableFraction => ({ speciesId: '', percent: '100' });

function action(state: WorkspaceState, kind: ActorActionKind, subjectId: string | null = null): AvailableAction | undefined {
  return state.actions.find(item => item.kind === kind && item.subjectId === subjectId);
}

function readable(value: string | null | undefined): string {
  if (!value) return 'Not established';
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
        <select className="select-input" aria-label={`${label} leaflet lipid ${index + 1}`} value={row.speciesId} onChange={event => update(index, { speciesId: event.target.value })}>
          <option value="">Choose lipid</option>
          {options.map(lipid => <option value={lipid.speciesId} key={lipid.speciesId}>{lipid.displayName} ({lipid.speciesId})</option>)}
        </select>
        <input className="number-input tabular" type="number" min="0" max="100" step="0.1" aria-label={`${label} leaflet percentage ${index + 1}`} value={row.percent} onChange={event => update(index, { percent: event.target.value })} />
        <button className="button compact" type="button" aria-label={`Remove ${label} leaflet lipid ${index + 1}`} disabled={rows.length === 1} onClick={() => onChange(rows.filter((_, current) => current !== index))}>×</button>
      </div>
      {options.find(lipid => lipid.speciesId === row.speciesId)?.limitations.map(limit => <p className="help-text" key={limit}>{limit}</p>)}
    </div>)}
    <button className="button compact" type="button" onClick={() => onChange([...rows, { speciesId: '', percent: '0' }])}>Add lipid</button>
  </div>;
}

function validFractions(rows: EditableFraction[]): boolean {
  if (rows.length === 0 || rows.some(row => !row.speciesId || !Number.isFinite(Number(row.percent)) || Number(row.percent) <= 0)) return false;
  if (new Set(rows.map(row => row.speciesId)).size !== rows.length) return false;
  return Math.abs(rows.reduce((sum, row) => sum + Number(row.percent), 0) - 100) < 0.01;
}

function toFractions(rows: EditableFraction[]): LipidFraction[] {
  return rows.map(row => ({ speciesId: row.speciesId, fraction: Number(row.percent) / 100 }));
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
  const [feedback, setFeedback] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [workflowOpen, setWorkflowOpen] = useState(true);
  const busyRef = useRef(false);
  const refreshIndex = useRef(0);

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
    events.onmessage = () => { void refresh(); };
    events.onerror = () => setCommunication('The live connection is interrupted. Reconnecting to the local workspace…');
    return () => events.close();
  }, []);

  useEffect(() => {
    if (!state?.membrane) return;
    const toEditable = (fractions: LipidFraction[]) => fractions.map(item => ({ speciesId: item.speciesId, percent: String(item.fraction * 100) }));
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

  async function command(kind: ActorCommandKind, data: object): Promise<boolean> {
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
      setFeedback(cause instanceof Error ? cause.message : 'The local request did not complete.');
      await refresh();
      return false;
    } finally {
      busyRef.current = false;
      setBusy(false);
    }
  }

  async function inspectSubject(subjectId: string) {
    if (await command('selectInspectionSubject', { subjectId })) setWorkflowOpen(false);
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
    if (!await command('exportStage', { stageId })) return;
    // Let the browser stream the verified bundle instead of copying it into JavaScript memory.
    const link = document.createElement('a');
    link.href = `/api/export/${encodeURIComponent(stageId)}`;
    link.download = `protein-membrane-${stageId}.zip`;
    link.click();
  }

  if (!state) return <main className="workspace"><header className="workspace-header"><div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation</div></div></div></header><div className="scene-empty" role="status"><strong>Connecting to the local workspace</strong><span>{communication ?? 'Reading the current scientific account…'}</span></div></main>;

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
  const selectedChange = state.protein?.changes.find(change => change.id === state.inspection?.subjectId);

  return <main className={`workspace ${workflowOpen ? 'workflow-open' : 'workflow-closed'}`}>
    <header className="workspace-header">
      <div className="brand-lockup"><div className="brand-mark">PM</div><div><div className="brand-name">Protein–Membrane Workspace</div><div className="brand-subtitle">Local scientific preparation · one planar bilayer</div></div></div>
      <div className="header-context">
        <button className="button workflow-toggle" type="button" aria-expanded={workflowOpen} aria-controls="researcher-workflow"
          onClick={() => setWorkflowOpen(value => !value)}>{workflowOpen ? 'Hide workflow' : 'Show workflow'}</button>
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
        {communication && <div className="notice warning" role="alert">{communication}<div className="button-row"><button className="button compact" type="button" onClick={() => void refresh()}>Refresh account</button></div></div>}
        {feedback && <div className="notice warning" role="alert">{feedback}</div>}
        {state.notices.length > 0 && <div className="notice-list" aria-label="Scientific and workflow notices">
          {state.notices.map(notice => <div key={notice.id} className={`notice ${notice.severity.toLowerCase()}`}><strong>{readable(notice.severity)}</strong> · {notice.message}{notice.subjectId && <small className="tabular">Subject {notice.subjectId}</small>}</div>)}
        </div>}

        <section className="rail-section">
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

        <section className="rail-section">
          <h3 className="section-title">3 · Planar membrane</h3>
          <p className="help-text">Specify each leaflet. Percentages express the intended model, not achieved molecule counts.</p>
          {state.availableLipids.length === 0 && <div className="hint-box">No lipid catalogue is currently qualified for selection.</div>}
          <FractionEditor label="Upper" rows={upper} onChange={setUpper} options={state.availableLipids} />
          <FractionEditor label="Lower" rows={lower} onChange={setLower} options={state.availableLipids} />
          <label className="field-label" htmlFor="scientific-purpose">Scientific purpose</label>
          <input id="scientific-purpose" className="text-input" value={membranePurpose} onChange={event => setMembranePurpose(event.target.value)} placeholder="Why this membrane model?" />
          <div className="button-row"><ActionButton state={state} kind="proposeMembrane" busy={busy || !membraneReady} onClick={() => void command('proposeMembrane', { upper: toFractions(upper), lower: toFractions(lower), scientificPurpose: membranePurpose.trim() })}>Propose membrane model</ActionButton></div>
          {action(state, 'proposeMembrane')?.reason && <p className="action-reason">{action(state, 'proposeMembrane')?.reason}</p>}
          {!membraneReady && <p className="help-text">Each leaflet needs distinct qualified lipids totaling 100%, plus a study purpose.</p>}
          {state.membrane && <><div className="hint-box"><strong>{readable(state.membrane.status)}</strong><br /><span className="tabular">{state.membrane.modelId}</span><br />Purpose: {state.membrane.scientificPurpose}{state.membrane.limitations.length > 0 && <><br />{state.membrane.limitations.join('; ')}</>}</div><div className="button-row"><ActionButton state={state} kind="selectInspectionSubject" busy={busy} onClick={() => void inspectSubject(state.membrane!.modelId)}>Inspect membrane</ActionButton><ActionButton state={state} kind="adoptMembrane" busy={busy} onClick={() => void command('adoptMembrane', { modelId: state.membrane!.modelId })}>Adopt and assess model</ActionButton></div></>}
          {state.study && <p className="help-text tabular">Fixed conditions: pH {state.study.conditions.nominalPh}; NaCl {state.study.conditions.targetNaClMolar} M; temperature {state.study.conditions.optionalTemperatureKelvin} K.</p>}
        </section>

        <section className="rail-section">
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
          <div className="button-row"><ActionButton state={state} kind="startPreparation" busy={busy} variant="primary" onClick={() => void command('startPreparation', {})}>Construct &amp; minimize</ActionButton></div>
          {action(state, 'startPreparation')?.reason && <p className="action-reason">{action(state, 'startPreparation')?.reason}</p>}
          {state.attempt && <div className="attempt-account">
            <strong className="tabular">Attempt {state.attempt.attemptId}</strong>
            <span>{readable(state.attempt.status)}{state.attempt.stageKind && ` · ${readable(state.attempt.stageKind)}`}</span>
            <p>{state.attempt.message}</p>
            {state.attempt.progress !== null && <progress max="1" value={state.attempt.progress} aria-label="Observed attempt progress" />}
            <div className="button-row"><ActionButton state={state} kind="stopAttempt" busy={busy} variant="danger" onClick={() => void command('stopAttempt', { attemptId: state.attempt!.attemptId })}>Stop unfinished work</ActionButton></div>
          </div>}
          {selectedStage && <div className="stage-export">
            <strong>{readable(selectedStage.kind)} · {readable(selectedStage.status)}</strong>
            <p className="tabular">Attempt {selectedStage.attemptId} · study {selectedStage.studyRevisionId}</p>
            <p>{selectedStage.assessment ? `${selectedStage.assessment.currentlyApplicable ? readable(selectedStage.assessment.qualification) : `Historical ${readable(selectedStage.assessment.qualification)}; assessment not current`} — ${selectedStage.assessment.reason}` : 'No scientific assessment established for this stage.'}</p>
            <div className="button-row"><ActionButton state={state} kind="exportStage" subjectId={selectedStage.stageId} busy={busy} onClick={() => void exportStage(selectedStage.stageId)}>Export this completed stage</ActionButton></div>
            {action(state, 'exportStage', selectedStage.stageId)?.reason && <p className="action-reason">{action(state, 'exportStage', selectedStage.stageId)?.reason}</p>}
            {selectedStage.kind === 'Minimization' && <><div className="button-row"><ActionButton state={state} kind="requestEquilibration" subjectId={selectedStage.stageId} busy={busy} onClick={() => void command('requestEquilibration', { stageId: selectedStage.stageId })}>Request optional equilibration from this stage</ActionButton></div>{action(state, 'requestEquilibration', selectedStage.stageId)?.reason && <p className="action-reason">{action(state, 'requestEquilibration', selectedStage.stageId)?.reason}</p>}</>}
          </div>}
        </section>
      </aside>

      <ConnectedStructuralInspection state={state} onSelectFocus={annotationId => void command('setInspectionFocus', { annotationId })}
        proposalDecision={selectedChange && <div className="proposal-decision" aria-label="Preparation change decision">
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
        </div>} />
    </div>

    <nav className="stage-strip" aria-label="Completed and attempted stages">
      <div className="stage-summary">
        <div><span className="stage-strip-heading">Completed stages</span><strong>{state.stages.length === 0 ? 'None' : state.stages.length}</strong>{state.stages.length === 0 && <small>No completed stage is currently established.</small>}</div>
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
