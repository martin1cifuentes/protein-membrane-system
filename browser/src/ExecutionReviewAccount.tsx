import type { WorkspaceState } from './ProteinInMembraneWorkspace';

type AttemptAccount = NonNullable<WorkspaceState['attempt']>;
type StageAccount = WorkspaceState['stages'][number];
type SpeciesCount = NonNullable<AttemptAccount['derivation']>['lipidCounts'][number];

function label(value: string | null | undefined): string {
  if (!value) return 'Not established';
  if (value === 'readyForMinimization') return 'Ready for minimization';
  return value.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/_/g, ' ');
}

function currentAttemptStages(state: WorkspaceState, attempt: AttemptAccount): StageAccount[] {
  return state.stages.filter(stage => stage.attemptId === attempt.attemptId);
}

function counts(items: SpeciesCount[], side: 'upper' | 'lower'): string {
  const selected = items.filter(item => item.physicalSide === side && item.count > 0);
  return selected.length ? selected.map(item => `${item.speciesId} ${item.count.toLocaleString()}`).join(', ') : 'None established';
}

function cell(value: number[]): string {
  return value.length === 3 && value.every(Number.isFinite)
    ? `${value.map(item => item.toFixed(1)).join(' × ')} Å` : 'Not established';
}

function measured(value: {value: number; unit: string}): string {
  return `${Number.isFinite(value.value) ? value.value.toLocaleString(undefined, { maximumSignificantDigits: 5 }) : 'Unavailable'} ${value.unit}`;
}

export function AttemptReviewAccount({ state, attempt, onInspectConstructed }: {
  state: WorkspaceState;
  attempt: AttemptAccount;
  onInspectConstructed: (subjectId: string) => void;
}) {
  const stages = currentAttemptStages(state, attempt);
  const hasIdentity = attempt.attemptId.length > 0;
  const status = attempt.status.toLowerCase();
  const running = status === 'running' || status === 'pending';
  const progress = attempt.progress !== null && Number.isFinite(attempt.progress)
    ? Math.min(100, Math.max(0, attempt.progress * 100)) : null;
  const actual = attempt.constructed;
  const planned = attempt.derivation;
  const originCurrent = attempt.studyRevisionId === state.study?.id;
  const currentSource = state.study?.id === attempt.studyRevisionId
    ? state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId)?.label : null;
  const membraneName = state.study?.id === attempt.studyRevisionId && state.membrane?.status === 'assessed'
    ? [...new Set([...state.membrane.upper, ...state.membrane.lower]
      .filter(item => item.fraction > 0).map(item => item.speciesId))].join(' / ') : null;
  const inspectAction = state.actions.find(item => item.kind === 'selectInspectionSubject' && item.subjectId === null);
  const awaitingMinimization = status === 'readyforminimization';

  return <>
    <section className="account-card execution-identity-account" aria-label="Current preparation attempt">
      <h2>{attempt.stageKind === 'Minimization' && running ? 'Current minimization' : awaitingMinimization ? 'Constructed system review' : 'Current preparation attempt'}</h2>
      <dl className="detail-grid">
        <dt>Attempt</dt><dd>{hasIdentity ? `Accepted · ${label(attempt.status)}` : status === 'pending' ? 'Admission pending' : `Not admitted · ${label(attempt.status)}`}</dd>
        {attempt.studyRevisionId && <><dt>Study context</dt><dd>{originCurrent ?
          `Current study revision ${state.study?.number ?? 'unknown'}` :
          `Surviving attempt from earlier study revision ${attempt.studyRevisionId}; current revision ${state.study?.number ?? 'unknown'}`}</dd></>}
        <dt>Constructed system</dt><dd>{actual ? 'Validated protein + bilayer + water and ions' : 'Not yet validated'}</dd>
        <dt>Protein</dt><dd>{currentSource ?? 'Bound protein source'}</dd>
        <dt>Bilayer</dt><dd>{actual ? `${membraneName ?? 'Selected species'} · constructed planar bilayer` : membraneName ? `${membraneName} target` : 'No constructed bilayer claim'}</dd>
        <dt>Water and ions</dt><dd>{actual ? `${actual.waterCount.toLocaleString()} water · ${actual.sodiumCount.toLocaleString()} Na⁺ · ${actual.chlorideCount.toLocaleString()} Cl⁻` : 'Actual counts not established'}</dd>
        <dt>Status</dt><dd className={status === 'running' ? 'execution-running' : ''}>{attempt.stageKind ? `${label(attempt.stageKind)} · ${label(attempt.status)}` : label(attempt.status)}</dd>
      </dl>
      {actual && <div className="button-row"><button className="button compact" type="button"
        disabled={inspectAction?.enabled !== true} title={inspectAction?.reason ?? undefined}
        onClick={() => onInspectConstructed(actual.subjectId)}>Inspect verified constructed system</button></div>}
      <details className="execution-exact-identity"><summary>Exact bound attempt and policy</summary>
        <dl className="detail-grid">
          <dt>Attempt ID</dt><dd className="tabular">{hasIdentity ? attempt.attemptId : 'Not yet accepted'}</dd>
          <dt>Study revision</dt><dd className="tabular">{attempt.studyRevisionId ?? 'Not yet bound'}</dd>
          <dt>Policy</dt><dd className="tabular">{attempt.policyId ?? 'Not yet bound'}{attempt.policyVersion && ` · ${attempt.policyVersion}`}</dd>
          {actual && <><dt>Constructed subject</dt><dd className="tabular">{actual.subjectId}</dd><dt>Atoms</dt><dd>{actual.atomCount.toLocaleString()}</dd></>}
        </dl>
      </details>
    </section>

    <section className="account-card execution-progress-account" aria-label="Observed preparation progress">
      <h2>Observed progress</h2>
      <div className="execution-progress-body">
        {running && <span className="execution-progress-spinner" aria-hidden="true" />}
        <div>
          <p>{attempt.message}</p>
          <strong>{stages.length === 0 ? 'No completed stage yet.' : `${stages.length} completed stage${stages.length === 1 ? '' : 's'} retained.`}</strong>
        </div>
      </div>
      {progress !== null && <><progress max="100" value={progress} aria-label="Observed attempt progress" />
        <p className="help-text">{progress.toFixed(0)}% of the reported current operation; this is not a completion or qualification claim.</p></>}
    </section>

    <section className="account-card execution-basis-account" aria-label="Preparation basis">
      <h2>Preparation basis</h2>
      <dl className="detail-grid">
        <dt>Placement</dt><dd>{attempt.studyRevisionId && state.study?.id === attempt.studyRevisionId && state.placement?.status === 'supported'
          ? 'Supported placement at admitted revision' : attempt.studyRevisionId ? 'Bound to earlier study revision' : 'Admission pending'}</dd>
        <dt>Chemistry</dt><dd>{actual ? 'Combined system represented and validated' : planned ? 'Provider-selected values recorded; validation pending' : 'Not yet established'}</dd>
        <dt>Stage</dt><dd>{attempt.stageKind ? `${label(attempt.stageKind)} · ${label(attempt.status)}` : label(attempt.status)}</dd>
      </dl>
      {planned && <details className="execution-candidate-values" aria-label="Provider-selected construction values" open={awaitingMinimization}>
        <summary>Provider-selected candidate</summary>
        <p className="help-text">These finite counts and cell come from the same construction invocation as the molecular candidate.</p>
        <dl className="detail-grid">
          <dt>Upper physical leaflet</dt><dd>{counts(planned.lipidCounts, 'upper')}</dd>
          <dt>Lower physical leaflet</dt><dd>{counts(planned.lipidCounts, 'lower')}</dd>
          <dt>Candidate cell</dt><dd>{cell(planned.cellAngstrom)}</dd>
          <dt>Water</dt><dd>{planned.waterCount.toLocaleString()}</dd>
          <dt>Na⁺ / Cl⁻</dt><dd>{planned.sodiumCount.toLocaleString()} / {planned.chlorideCount.toLocaleString()}</dd>
          <dt>Target / estimated NaCl</dt><dd>{planned.intendedNaClMolar.toFixed(3)} / {planned.estimatedNaClMolar.toFixed(3)} M</dd>
        </dl>
        {planned.approximations.length > 0 && <div className="execution-approximations">
          <strong>Approximations</strong>
          {planned.approximations.map((item, index) => <p key={`${item}:${index}`}>{item}</p>)}
        </div>}
        {planned.limitations.map((item, index) => <p className="help-text" key={`${item}:${index}`}>Limit: {item}</p>)}
      </details>}
      {actual && <details className="execution-candidate-values execution-actual-values" aria-label="Validated actual constructed system" open={awaitingMinimization}>
        <summary>Validated actual system</summary>
        <dl className="detail-grid">
          <dt>Upper physical leaflet</dt><dd>{counts(actual.achievedComposition, 'upper')}</dd>
          <dt>Lower physical leaflet</dt><dd>{counts(actual.achievedComposition, 'lower')}</dd>
          <dt>Actual cell</dt><dd>{cell(actual.cellAngstrom)}</dd>
          <dt>Actual water</dt><dd>{actual.waterCount.toLocaleString()}</dd>
          <dt>Actual Na⁺ / Cl⁻</dt><dd>{actual.sodiumCount.toLocaleString()} / {actual.chlorideCount.toLocaleString()}</dd>
          <dt>Local observation</dt><dd>{label(actual.localState?.standing)}</dd>
        </dl>
        <p className="help-text">{actual.conditionsTreatment}</p>
        <details className="execution-exact-identity"><summary>Inspect local construction observations</summary>
          {actual.localState?.measurements?.map(item => <p key={`${item.name}:${item.scope}`}>{label(item.name)}: {measured(item)} · {item.scope}</p>)}
          {actual.localState?.limitations?.map((item, index) => <p key={`${item}:${index}`}>Limit: {item}</p>)}
        </details>
      </details>}
      {awaitingMinimization && <p className="help-text">Candidate values and the validated actual system are separate accounts. Required minimization starts only after the researcher continues this exact attempt.</p>}
    </section>
  </>;
}

export function StageReviewAccount({ state, stage, exportFault }: {
  state: WorkspaceState;
  stage: StageAccount;
  exportFault: string | null;
}) {
  const actual = stage.constructed;
  const assessment = stage.assessment;
  const observation = stage.observation;
  const qualification = assessment ? label(assessment.qualification) : 'Not established';
  const originCurrent = stage.studyRevisionId === state.study?.id;
  const sourceStage = stage.kind === 'Equilibration' && stage.sourceStageId
    ? state.stages.find(item => item.stageId === stage.sourceStageId &&
      item.attemptId === stage.attemptId && item.kind === 'Minimization') : undefined;
  const finalForce = observation?.measurements.find(item => item.name === 'finalRmsForce' &&
    item.scope === 'final unrestrained constraint tangent; per particle');
  const completed = stage.status.toLowerCase() === 'completed';

  return <>
    <section className="account-card execution-identity-account" aria-label="Completed stage information">
      <h2>Stage information</h2>
      <dl className="detail-grid">
        <dt>Attempt</dt><dd>{stage.attemptId === state.attempt?.attemptId ? 'Current attempt (completed)' : 'Earlier accepted attempt'}</dd>
        {!originCurrent && <><dt>Study context</dt><dd>
          Historical stage from study revision {stage.studyRevisionId}; current revision {state.study?.number ?? 'unknown'}
        </dd></>}
        <dt>Constructed system</dt><dd>{actual ? 'Verified protein + bilayer + water and ions' : 'Historical construction account unavailable here'}</dd>
        <dt>Protein</dt><dd>{stage.studyRevisionId === state.study?.id
          ? state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId)?.label ?? 'Selected prepared protein'
          : 'Bound protein at source revision'}</dd>
        <dt>Bilayer</dt><dd>{actual ? `${counts(actual.achievedComposition, 'upper')} / ${counts(actual.achievedComposition, 'lower')} physical leaflets` : 'Bound constructed bilayer'}</dd>
        <dt>Water and ions</dt><dd>{actual ? `${actual.waterCount.toLocaleString()} water · ${actual.sodiumCount.toLocaleString()} Na⁺ · ${actual.chlorideCount.toLocaleString()} Cl⁻` : 'Actual counts unavailable here'}</dd>
      </dl>
      <details className="execution-exact-identity"><summary>Exact stage, attempt, revision and constructed identity</summary>
        <dl className="detail-grid">
          <dt>Stage ID</dt><dd className="tabular">{stage.stageId}</dd>
          <dt>Attempt ID</dt><dd className="tabular">{stage.attemptId}</dd>
          <dt>Study revision</dt><dd className="tabular">{stage.studyRevisionId}</dd>
          {actual && <><dt>Constructed subject</dt><dd className="tabular">{actual.subjectId}</dd><dt>Actual atoms</dt><dd>{actual.atomCount.toLocaleString()}</dd><dt>Construction cell</dt><dd>{cell(actual.cellAngstrom)}</dd></>}
        </dl>
      </details>
    </section>

    <section className={`account-card execution-assessment-account${exportFault ? ' export-validation' : ''}`}
      aria-label={exportFault ? 'Export validation and unchanged stage review' :
        stage.kind === 'Minimization' ? 'Minimized stage review and distinct scientific assessment' :
          'Equilibration stage review and distinct scientific assessment'}>
      {exportFault ? <>
        <h2>Export validation</h2>
        <div className="export-failure" role="alert">
          <span className="export-failure-icon" aria-hidden="true">!</span>
          <div><strong>Export not delivered</strong><span>{exportFault}</span></div>
          <p>The {stage.kind === 'Minimization' ? 'minimized' : 'equilibrated'} stage remains completed, but the export could not be delivered. The molecular stage and scientific assessment are unchanged.</p>
        </div>
      </> : <h2>{stage.kind === 'Minimization' ? 'Minimized stage review' : 'Equilibration stage'}</h2>}
      <dl className="detail-grid">
        <dt>Selected stage</dt><dd>{completed ? stage.kind === 'Minimization' ? 'Minimized' : 'Equilibrated' : `${label(stage.kind)} · ${label(stage.status)}`}</dd>
        {stage.kind === 'Equilibration' && <>
          <dt>Earlier stage</dt><dd>{sourceStage ? 'Minimized — available' : 'Source minimized stage unavailable'}</dd>
          <dt>Optional procedure</dt><dd>{completed ? 'Completed' : label(stage.status)}</dd>
        </>}
        <dt>Construction correspondence</dt><dd>{completed && observation && actual ? 'Bound to verified constructed system' : 'Not established here'}</dd>
        {stage.kind === 'Minimization' && <><dt>Unrestrained convergence</dt><dd>{observation?.termination === 'converged' ? 'Observed' : 'Not established here'}</dd></>}
        {finalForce && <><dt>Final RMS force</dt><dd>{measured(finalForce)}</dd></>}
        <dt>Scientific assessment</dt><dd>{assessment?.currentlyApplicable === false ? `Historical ${qualification}` : qualification}</dd>
        {exportFault && <><dt>Reason</dt><dd>{assessment?.reason ?? 'No stage-specific scientific assessment is established.'}</dd></>}
      </dl>
      {!exportFault && <><p>{assessment?.reason ?? 'No stage-specific scientific assessment is established.'}</p>
        <p className="help-text">{stage.kind === 'Minimization'
          ? 'Factual completion and scientific qualification are separate. This stage is not automatically equilibrated or production ready.'
          : 'Procedure completion and scientific qualification are separate. This stage does not establish global thermodynamic equilibrium or production study validity.'}</p></>}
    </section>

    {(observation || assessment) &&
      <section className="account-card execution-findings-account" aria-label="Stage-specific findings and evidence">
        <h2>Stage-specific findings and evidence</h2>
        {observation && <details><summary>Observed final operation and measurements</summary>
          <p>Termination: {label(observation.termination)} · provider {observation.providerVersion}</p>
          {observation.measurements.map(item => <p key={`${item.name}:${item.scope}`}>{label(item.name)}: {measured(item)} · {item.scope}</p>)}
          <p>Fresh local-state observation: {observation.localState ? label(observation.localState.standing) : 'Unavailable'}{observation.localState?.unavailableReason && ` · ${observation.localState.unavailableReason}`}</p>
          {observation.localState?.rolePairMeasurements?.map((pair, index) => <p key={`${pair.firstMoleculeRole}:${pair.secondMoleculeRole}:${index}`}>
            {label(pair.firstMoleculeRole)}–{label(pair.secondMoleculeRole)} contact: {pair.pairsWithinSearchRadius.toLocaleString()} pairs{pair.minimumDistanceAngstrom !== null && ` · nearest ${pair.minimumDistanceAngstrom.toFixed(2)} Å`}
          </p>)}
          {observation.localState?.measurements?.map(item => <p key={`local:${item.name}:${item.scope}`}>
            {label(item.name)}: {measured(item)} · {item.scope}
          </p>)}
          <p>Fresh protein geometry: {observation.proteinGeometry ? label(observation.proteinGeometry.standing) : 'Unavailable'}</p>
          {observation.proteinGeometry?.kinds?.map(item => <p key={item.kind}>
            {label(item.kind)}: {label(item.standing)} · {item.measuredCount}/{item.eligibleCount} measured{item.minimumDistanceAngstrom !== null && ` · shortest ${item.minimumDistanceAngstrom.toFixed(2)} Å`}{item.maximumDistanceAngstrom !== null && ` · longest ${item.maximumDistanceAngstrom.toFixed(2)} Å`}{item.unavailableReason && ` · ${item.unavailableReason}`}
          </p>)}
        </details>}
        {assessment?.findings.map(finding => <div className="finding-item" key={finding.id}>
          <strong>{label(finding.disposition)} · {finding.consequence}</strong><span>{finding.meaning}</span>
        </div>)}
        {assessment?.evidence.map(evidence => <div className="finding-item" key={evidence.id}>
          <strong>{evidence.source} · {evidence.method}</strong><span>{evidence.observation}</span>
          <span>Applicability: {evidence.applicability}</span><span>Uncertainty: {evidence.uncertainty}</span>
        </div>)}
        {assessment?.limitations.map((limitation, index) => <p key={`${limitation}:${index}`}>Limit: {limitation}</p>)}
      </section>}
  </>;
}
