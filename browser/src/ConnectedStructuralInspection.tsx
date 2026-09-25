import { useEffect, useRef, useState, type ReactNode } from 'react';
import { Viewer } from 'molstar/lib/apps/viewer/app';
import type { WorkspaceState } from './ProteinInMembraneWorkspace';

export interface StructureFocus {
  authAsymId: string;
  authSeqId: number;
  insertionCode: string | null;
  authAtomId: number | null;
}

export interface InspectionAnnotation {
  id: string;
  subjectPartId: string;
  label: string;
  meaning: string;
  evidenceId: string | null;
  geometryFocus: StructureFocus | null;
}

export interface InspectionMetric {
  name: string;
  value: string;
  unit: string | null;
  subjectPartId: string;
  evidenceId: string | null;
}

export interface InspectionAccount {
  subjectId: string;
  structureUrl: string | null;
  representationKind: string;
  omittedMolecules: string[];
  focusId: string | null;
  focus: StructureFocus | null;
  annotations: InspectionAnnotation[];
  metrics: InspectionMetric[];
}

function localStructureUrl(path: string): string {
  const url = new URL(path, window.location.origin);
  if (url.origin !== window.location.origin || !url.pathname.startsWith('/api/structures/')) {
    throw new Error('The selected structure is not available from this local workspace.');
  }
  return url.toString();
}

function MolecularScene({ inspection }: { inspection: InspectionAccount | null }) {
  const mount = useRef<HTMLDivElement>(null);
  const [viewer, setViewer] = useState<Viewer | null>(null);
  const [error, setError] = useState<string | null>(null);
  const structureUrl = inspection?.structureUrl ?? null;
  const hasSelectedFocus = useRef(false);
  hasSelectedFocus.current = !!inspection?.focus;

  useEffect(() => {
    if (!structureUrl || !mount.current) return;
    let closed = false;
    const active = { viewer: null as Viewer | null };
    const target = mount.current;
    setError(null);
    setViewer(null);

    async function load() {
      try {
        const url = localStructureUrl(structureUrl!);
        const created = await Viewer.create(target, {
          extensions: [],
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowRemoteState: false,
          layoutShowSequence: false,
          layoutShowLog: false,
          layoutShowLeftPanel: false,
          viewportShowExpand: false,
          viewportShowControls: false,
          viewportShowSettings: false,
          viewportShowSelectionMode: true,
          viewportShowAnimation: false,
          viewportBackgroundColor: '#edf3f6',
        });
        if (closed) { created.dispose(); return; }
        active.viewer = created;
        const format = new URL(url).searchParams.get('format');
        if (format !== 'pdb' && format !== 'mmcif')
          throw new Error('The selected structure has no recognized PDB or mmCIF representation.');
        await created.loadStructureFromUrl(url, format, false, { label: inspection?.subjectId });
        if (!closed) setViewer(created);
      } catch (cause) {
        if (!closed) setError(cause instanceof Error ? cause.message : 'The structure could not be rendered.');
      }
    }

    void load();
    return () => {
      closed = true;
      active.viewer?.dispose();
      target.replaceChildren();
    };
  }, [structureUrl]);

  useEffect(() => {
    if (!viewer) return;
    const focus = inspection?.focus;
    if (!focus || !focus.authAsymId || !Number.isInteger(focus.authSeqId)) {
      viewer.structureInteractivity({ action: 'select' });
      return;
    }
    viewer.structureInteractivity({
      action: ['select', 'focus'],
      elements: {
        auth_asym_id: focus.authAsymId,
        auth_seq_id: focus.authSeqId,
        ...(focus.insertionCode ? { pdbx_PDB_ins_code: focus.insertionCode } : {}),
        ...(focus.authAtomId !== null ? { atom_id: focus.authAtomId } : {}),
      },
    });
  }, [viewer, inspection?.focusId, inspection?.focus?.authAsymId, inspection?.focus?.authSeqId, inspection?.focus?.insertionCode, inspection?.focus?.authAtomId]);

  useEffect(() => {
    if (!viewer || !mount.current) return;
    let resetTimer: number | undefined;
    const observer = new ResizeObserver(() => {
      viewer.handleResize();
      window.clearTimeout(resetTimer);
      resetTimer = window.setTimeout(() => {
        if (!hasSelectedFocus.current) viewer.plugin.canvas3d?.requestCameraReset({ durationMs: 0 });
      }, 80);
    });
    observer.observe(mount.current);
    return () => { observer.disconnect(); window.clearTimeout(resetTimer); };
  }, [viewer]);

  if (!structureUrl) {
    return <div className="scene-empty"><strong>{inspection ? 'No renderable structure is established' : 'No molecular subject selected'}</strong><span>{inspection ? 'The selected subject’s established values and limitations remain available beside this view.' : 'Choose a source or completed subject to inspect its available structure.'}</span></div>;
  }

  return <>
    <div className="viewer-mount" ref={mount} aria-label={`3D structure for ${inspection?.subjectId ?? 'selected subject'}`} />
    {!viewer && !error && <div className="scene-loading" role="status">Loading the identified structure…</div>}
    {error && <div className="scene-error" role="alert"><strong>Structure unavailable</strong><span>{error}</span></div>}
  </>;
}

function statusClass(status: string): string {
  const value = status.toLowerCase().replace(/[^a-z]/g, '');
  if (value.includes('fail') || value.includes('declined') || value.includes('unsupported') || value.includes('notsupported') || value.includes('notqualified') || value.includes('unqualified') || value.includes('disqualified') || value.includes('refused')) return 'danger';
  if (value.includes('notestablished') || value.includes('notcurrent') || value.includes('indeterminate') || value.includes('unavailable')) return 'warning';
  if (value.includes('qualified') || value.includes('supported') || value.includes('completed')) return 'success';
  return '';
}

function representationLabel(kind: string): string {
  switch (kind) {
    case 'structuralSource': return 'Structural source';
    case 'preparedProtein': return 'Prepared protein';
    case 'selected-protein-before-repair': return 'Selected protein before repair';
    default: return kind;
  }
}

export function ConnectedStructuralInspection({
  state,
  onSelectFocus,
  proposalDecision,
}: {
  state: WorkspaceState;
  onSelectFocus: (annotationId: string | null) => void;
  proposalDecision?: ReactNode;
}) {
  const inspection = state.inspection;
  const subject = inspection?.subjectId ?? state.study?.id ?? null;
  const stage = state.stages.find(item => item.stageId === subject);
  const preparationChange = state.protein?.changes.find(change => change.id === subject);
  const subjectStatus = stage?.assessment && !stage.assessment.currentlyApplicable ? 'Assessment not current'
    : stage?.assessment?.qualification ?? stage?.status
    ?? (state.placement?.proposalId === subject ? state.placement.status : null)
    ?? (state.protein?.subjectId === subject ? state.protein.status : null)
    ?? (preparationChange ? state.protein?.status === 'declined' ? 'Declined proposal' : 'Preparation proposal' : null)
    ?? (state.membrane?.modelId === subject ? state.membrane.status : null);
  const selectedAction = state.actions.find(item => item.kind === 'setInspectionFocus' && item.subjectId === null);
  const canSelectFocus = selectedAction?.enabled === true;
  const affectedAnnotation = preparationChange && inspection?.annotations.find(item => item.geometryFocus);
  const evidenceContent = useRef<HTMLDivElement>(null);
  const compactReview = useRef<boolean | null>(null);
  const currentSourceContext = inspection && state.study?.selectedSourceId && (
    subject === state.study.selectedSourceId || subject === state.protein?.subjectId || preparationChange);

  useEffect(() => {
    if (!preparationChange) { compactReview.current = null; return; }
    let frame: number | undefined;
    const orientReview = () => {
      const compact = window.innerWidth <= 880;
      if (compactReview.current === compact) return;
      compactReview.current = compact;
      window.cancelAnimationFrame(frame ?? 0);
      frame = window.requestAnimationFrame(() => {
        const content = evidenceContent.current;
        const proposal = content?.querySelector('.proposal-account');
        if (!content || !proposal) return;
        content.scrollTop = compact
          ? proposal.getBoundingClientRect().top - content.getBoundingClientRect().top + content.scrollTop - 4
          : 0;
      });
    };
    window.addEventListener('resize', orientReview);
    orientReview();
    return () => { window.removeEventListener('resize', orientReview); window.cancelAnimationFrame(frame ?? 0); };
  }, [preparationChange?.id]);

  return <>
    <section className="scene-panel" aria-label="Molecular structure">
      <div className="scene-heading">
        <div>
          <p className="eyebrow">Connected structural inspection</p>
          <h1>{stage?.summary ?? (inspection && representationLabel(inspection.representationKind)) ?? 'Scientific subject'}</h1>
          <div className="scene-subtitle tabular">{subject ?? 'No subject selected'}{inspection && ` · ${representationLabel(inspection.representationKind)}`}</div>
        </div>
        {subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      <div className="scene-surface"><MolecularScene inspection={inspection} /></div>
      <div className="scene-caption">
        <strong>Spatial view is evidence, not assessment.</strong>{' '}
        {inspection?.omittedMolecules.length
          ? `Not shown: ${inspection.omittedMolecules.join(', ')}.`
          : 'Any available omissions and exact measurements appear with this subject’s account.'}
      </div>
    </section>

    <aside className="evidence-panel" aria-label="Evidence and assessment">
      <div className="evidence-content" ref={evidenceContent}>
      <p className="eyebrow">Exact subject account</p>
      <div className="status-line">
        <h2 className="panel-heading">Evidence &amp; standing</h2>
        {subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      {!inspection && <div className="hint-box">Select a scientific subject to connect its structure with the applicable values and findings.</div>}
      {inspection && <>
        {!preparationChange && <section className="account-card" aria-label="Selected structure">
          <h2>Selected structure</h2>
          <dl className="detail-grid">
          <dt>Subject</dt><dd>{inspection.subjectId}</dd>
          <dt>Representation</dt><dd>{representationLabel(inspection.representationKind)}</dd>
          {inspection.focusId && <><dt>Selected part</dt><dd>{inspection.focusId}</dd></>}
          </dl>
        </section>}
        {currentSourceContext && <section className="account-card source-account" aria-label="Source and assembly">
          <h2>Source and assembly</h2>
          <dl className="detail-grid">
            <dt>Source</dt><dd>{state.study!.selectedSourceId}</dd>
            {state.study!.selectedSourceKind && <><dt>Route</dt><dd>{state.study!.selectedSourceKind === 'rcsb' ? 'RCSB PDB' : state.study!.selectedSourceKind === 'alphafold' ? 'AlphaFold DB' : 'Researcher upload'}</dd></>}
            {state.study!.uploadProvenance && <><dt>Declared origin</dt><dd>{state.study!.uploadProvenance}</dd></>}
            {state.study!.modelIndex !== null && <><dt>Coordinate model</dt><dd>{state.study!.modelIndex}</dd></>}
            {state.study!.modelIndex !== null && <><dt>Assembly</dt><dd>{state.study!.biologicalAssemblyId ?? 'No assembly transformation'}</dd></>}
            {state.study!.chainIds.length > 0 && <><dt>Retained chains</dt><dd>{state.study!.chainIds.join(', ')}</dd></>}
          </dl>
          {state.study!.uploadProvenance && !preparationChange && <p className="help-text">Origin is researcher declared; it is not independently verified by this label.</p>}
        </section>}
        {preparationChange && <>
          <section className="account-card proposal-account" aria-label="Proposed change">
            <h2>Proposed change</h2>
            <dl className="detail-grid">
              <dt>Type</dt><dd>{preparationChange.kind === 'heavyAtom' ? 'Complete missing heavy atom' : preparationChange.kind === 'alternateLocation' ? 'Alternate location' : preparationChange.kind === 'residueState' ? 'Residue state' : 'Disulfide'}</dd>
              <dt>Proposal</dt><dd>{preparationChange.id}</dd>
              <dt>Study revision</dt><dd>{preparationChange.studyRevisionId}</dd>
            </dl>
            <p>{preparationChange.proposedChange}</p>
            <p className="help-text">Rationale: {preparationChange.rationale}</p>
          </section>
          <section className="account-card affected-account" aria-label="Affected region">
            <h2>Affected region</h2>
            <dl className="detail-grid">
              <dt>Selection</dt><dd>One residue</dd>
              <dt>Location</dt><dd>Model {preparationChange.residue.model} · chain {preparationChange.residue.chain} · residue {preparationChange.residue.residue}{preparationChange.residue.insertionCode} · copy {preparationChange.residue.copyId}</dd>
            </dl>
            {affectedAnnotation && <div className="button-row"><button className="button compact" type="button" disabled={!canSelectFocus}
              title={!canSelectFocus ? selectedAction?.reason ?? 'Selection is not available.' : undefined}
              onClick={() => onSelectFocus(affectedAnnotation.id)}>{inspection.focusId === affectedAnnotation.subjectPartId ? 'Affected region selected' : 'Focus affected region in structure'}</button></div>}
          </section>
        </>}
        {inspection.focusId && <div className="button-row"><button className="button compact" type="button" disabled={!canSelectFocus} onClick={() => onSelectFocus(null)}>Clear selected part</button></div>}
        {inspection.annotations.length > 0 && <section className="evidence-section">
          <h2>Located evidence and findings</h2>
          <div className="annotation-list">
            {inspection.annotations.map(annotation => <button
              className={`annotation-item ${inspection.focusId === annotation.subjectPartId ? 'selected' : ''}`}
              type="button"
              key={annotation.id}
              disabled={!canSelectFocus}
              title={!canSelectFocus ? selectedAction?.reason ?? 'Selection is not available.' : undefined}
              onClick={() => onSelectFocus(annotation.id)}
            >
              <span className="item-title">{annotation.label}</span>
              <span className="item-detail">{annotation.meaning}</span>
              <span className="item-detail tabular">{annotation.subjectPartId}{annotation.evidenceId && ` · Evidence ${annotation.evidenceId}`}</span>
            </button>)}
          </div>
        </section>}
        {inspection.metrics.length > 0 && <section className="evidence-section">
          <h2>Measured and derived values</h2>
          <dl className="detail-grid">
            {inspection.metrics.map(metric => <FragmentMetric key={`${metric.name}:${metric.subjectPartId}`} metric={metric} />)}
          </dl>
        </section>}
        {preparationChange && state.protein?.sourceGeometry && <section className="evidence-section" aria-label="Source geometry observations">
          <h2>Source geometry observations</h2>
          <p>Standing: {state.protein.sourceGeometry.standing}</p>
          {state.protein.sourceGeometry.kinds.map(item => <p key={item.kind}>{item.kind}: {item.standing} · {item.measuredCount}/{item.eligibleCount} eligible measured{item.minimumDistanceAngstrom !== null && ` · shortest ${item.minimumDistanceAngstrom.toFixed(2)} Å`}{item.maximumDistanceAngstrom !== null && ` · longest ${item.maximumDistanceAngstrom.toFixed(2)} Å`}{item.unavailableReason && ` · ${item.unavailableReason}`}</p>)}
          {state.protein.sourceGeometry.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
        </section>}
        {preparationChange && preparationChange.limitations.length > 0 && <section className="evidence-section">
          <h2>Proposal limitations</h2>
          {preparationChange.limitations.map(limit => <p key={limit}>{limit}</p>)}
        </section>}
        {inspection.omittedMolecules.length > 0 && <section className="evidence-section">
          <h2>Visual omissions</h2>
          <p>{inspection.omittedMolecules.join(', ')}</p>
        </section>}
      </>}
      {stage?.assessment && <section className="evidence-section">
        <h2>Scientific assessment · {stage.kind}</h2>
        <span className={`status-badge ${stage.assessment.currentlyApplicable ? statusClass(stage.assessment.qualification) : 'warning'}`}>
          {stage.assessment.currentlyApplicable ? stage.assessment.qualification : `Historical conclusion: ${stage.assessment.qualification}`}
        </span>
        <p>{stage.assessment.reason}</p>
        {!stage.assessment.currentlyApplicable && <div className="notice warning">This assessment is not currently applicable. A positive conclusion must be reassessed before it is presented as current.</div>}
        {stage.assessment.limitations.length > 0 && <p>Limits: {stage.assessment.limitations.join('; ')}</p>}
        {stage.assessment.findings.length > 0 && <div className="finding-list">
          {stage.assessment.findings.map(finding => <div className="finding-item" key={finding.id}><strong>{finding.disposition} · {finding.consequence}</strong><span>{finding.meaning}</span></div>)}
        </div>}
        {stage.assessment.evidence.length > 0 && <div className="finding-list">
          {stage.assessment.evidence.map(evidence => <div className="finding-item" key={evidence.id}><strong>{evidence.source} · {evidence.method}</strong><span>{evidence.observation} · Applicability: {evidence.applicability} · Uncertainty: {evidence.uncertainty}</span></div>)}
        </div>}
      </section>}
      {state.placement && state.placement.proposalId === subject && <section className="evidence-section">
        <h2>Placement proposal</h2>
        <p>{state.placement.reason}</p>
        <dl className="detail-grid">
          <dt>Topology</dt><dd>{state.placement.topologyKind}</dd>
          {state.placement.depthAngstrom !== null && <><dt>Depth</dt><dd>{state.placement.depthAngstrom.toFixed(2)} Å</dd></>}
          {state.placement.tiltDegrees !== null && <><dt>Tilt</dt><dd>{state.placement.tiltDegrees.toFixed(1)}°</dd></>}
          {state.placement.sidedness && <><dt>Sidedness</dt><dd>{state.placement.sidedness}</dd></>}
        </dl>
        {state.placement.evidence.map(evidence => <div className="finding-item" key={evidence.id}>
          <strong>{evidence.bearing} · {evidence.source}</strong>
          <span>{evidence.observation} · {evidence.uncertainty}</span>
        </div>)}
      </section>}
      {state.protein && state.protein.subjectId === subject && <section className="evidence-section">
        <h2>Protein preparation</h2>
        <p>{state.protein.summary}</p>
        {state.protein.atomCount !== null && <p className="tabular">{state.protein.atomCount.toLocaleString()} atoms</p>}
        {state.protein.findings.map(finding => <div className="finding-item" key={finding.id}>
          <strong>{finding.disposition} · {finding.consequence}</strong><span>{finding.meaning}</span>
        </div>)}
      </section>}
      </div>
      {proposalDecision}
    </aside>
  </>;
}

function FragmentMetric({ metric }: { metric: InspectionMetric }) {
  return <>
    <dt title={metric.subjectPartId}>{metric.name}</dt>
    <dd title={metric.evidenceId ?? undefined}>{metric.value}{metric.unit && ` ${metric.unit}`}</dd>
  </>;
}
