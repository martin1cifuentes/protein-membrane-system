import type { WorkspaceState } from './ProteinInMembraneWorkspace';
import { measurementLabel, operationLabel, runStartedLabel } from './executionDisplay';

type AttemptAccount = NonNullable<WorkspaceState['attempt']>;
type StageAccount = WorkspaceState['stages'][number];
type SpeciesCount = NonNullable<AttemptAccount['derivation']>['lipidCounts'][number];

function label(value: string | null | undefined): string {
  if (!value) return 'Not established';
  if (value === 'checksPassed') return 'Checks passed';
  if (value === 'issuesFound') return 'Issues found';
  if (value === 'checksIncomplete') return 'Checks incomplete';
  return value.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/_/g, ' ');
}

function diagnosticLabel(value: string): string {
  const words = label(value);
  return words.charAt(0).toUpperCase() + words.slice(1);
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
  return `${Number.isFinite(value.value) ? value.value.toLocaleString(undefined, { maximumSignificantDigits: 12 }) : 'Unavailable'} ${value.unit}`;
}

function MeasurementGroups({ values }: { values: { name: string; value: number; unit: string; scope: string }[] }) {
  const groups = new Map<string, typeof values>();
  values.forEach(value => groups.set(value.name, [...(groups.get(value.name) ?? []), value]));
  return <div className="measurement-groups">{[...groups].map(([name, items]) =>
    <details key={name}><summary>{measurementLabel(name)} · {items.length.toLocaleString()} measured</summary>
      <ul>{items.map((item, index) => <li key={`${item.scope}:${index}`}>
        <strong>{measured(item)}</strong><span>{item.scope}</span>
      </li>)}</ul>
    </details>)}</div>;
}

function atomAddress(value: { residue: { model?: number; chain: string; copyId: string; residue: number; insertionCode: string };
  atomName: string }): string {
  const residue = value.residue;
  return `${Number.isInteger(residue.model) ? `Model ${residue.model! + 1} · ` : ''}chain ${residue.chain}${residue.copyId && residue.copyId !== residue.chain ? `, copy ${residue.copyId}` : ''} · residue ${residue.residue}${residue.insertionCode} · atom ${value.atomName}`;
}

export function AttemptReviewAccount({ state, attempt, historical = false, onInspectSubject }: {
  state: WorkspaceState;
  attempt: AttemptAccount;
  historical?: boolean;
  onInspectSubject: (subjectId: string) => void;
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

  return <>
    <section className="account-card execution-identity-account" aria-label={historical ? 'Earlier preparation attempt' : 'Current preparation attempt'}>
      <h2>{historical ? 'Earlier preparation attempt' : attempt.stageKind === 'Minimization' && running ? 'Current minimization' : 'Current preparation attempt'}</h2>
      {historical && <p className="help-text">Read-only account from a prior attempt. Current Build and Stop controls apply only to the current attempt.</p>}
      <dl className="detail-grid">
        <dt>Run</dt><dd>{hasIdentity ? label(attempt.status) : status === 'pending' ? 'Waiting for admission' : `Request ${label(attempt.status)}`}</dd>
        {!originCurrent && attempt.studyRevisionId && <><dt>Inputs</dt><dd>This result uses earlier inputs.</dd></>}
        <dt>Constructed system</dt><dd>{actual ? 'Validated protein + bilayer + water and ions' : 'Not yet validated'}</dd>
        <dt>Protein</dt><dd>{currentSource ?? 'Bound protein source'}</dd>
        <dt>Bilayer</dt><dd>{actual ? `${membraneName ?? 'Selected species'} · constructed planar bilayer` : membraneName ? `${membraneName} target` : 'No constructed bilayer claim'}</dd>
        <dt>Water and ions</dt><dd>{actual ? `${actual.waterCount.toLocaleString()} water · ${actual.sodiumCount.toLocaleString()} Na⁺ · ${actual.chlorideCount.toLocaleString()} Cl⁻` : 'Actual counts not established'}</dd>
        <dt>Completed stage</dt><dd>{stages.length ? `${stages.length} retained` : 'None yet'}</dd>
      </dl>
      {actual && state.inspection?.subjectId !== actual.subjectId && <div className="button-row"><button className="button compact" type="button"
        disabled={inspectAction?.enabled !== true} title={inspectAction?.reason ?? undefined}
        onClick={() => onInspectSubject(actual.subjectId)}>Inspect verified constructed system</button></div>}
      {actual && <p className="help-text">{actual.atomCount.toLocaleString()} atoms in the checked constructed system. Final minimization is a separate observed stage.</p>}
    </section>

    <section className="account-card execution-progress-account" aria-label="Observed preparation progress">
      <h2>Observed progress</h2>
      <div className={`execution-progress-body${running ? ' is-running' : ''}`}>
        {running && <span className="execution-progress-spinner" aria-hidden="true" />}
        <div>
          <strong>{operationLabel(attempt.phase, status, attempt.stopRequested)}</strong>
          {!running && <p>{attempt.message}</p>}
          {attempt.stopRequested && running && <p>Stop was requested. The worker’s observed outcome is pending.</p>}
        </div>
      </div>
      {running && progress !== null ? <><progress max="100" value={progress} aria-label="Observed operation progress" />
        <p className="help-text">{progress.toFixed(0)}% reported for the current operation.</p></> :
        running && <progress aria-label="Operation in progress; completion fraction unknown" />}
      <details className="execution-diagnostics"><summary>Run details</summary>
        <p>{attempt.trials?.length ? `${attempt.trials.length} construction trial${attempt.trials.length === 1 ? '' : 's'} recorded.` : 'No construction trial recorded yet.'}</p>
        {attempt.failureCode && <p>Failure reported: {label(attempt.failureCode)}.</p>}
      </details>
      {!!attempt.trials?.length && <details className="execution-exact-identity"><summary>Construction trials and actual causes</summary>
        {attempt.trials.map(trial => <div className="finding-item" key={trial.trialId}>
          <strong>Trial {trial.trialIndex + 1} · {label(trial.standing)}</strong>
          <span>{trial.message ?? (trial.failureCode ? label(trial.failureCode) : 'Provider account recorded.')}</span>
          <span>Geometry: {trial.lateralPaddingAngstrom} Å lateral, {trial.aqueousPaddingAngstrom} Å aqueous · cell {cell(trial.actualCellAngstrom)}</span>
          {!!trial.proposedLipidCounts?.length && <span>Proposed upper / lower: {counts(trial.proposedLipidCounts, 'upper')} / {counts(trial.proposedLipidCounts, 'lower')}</span>}
          {!!trial.achievedLipidCounts?.length && <span>Actual upper / lower: {counts(trial.achievedLipidCounts, 'upper')} / {counts(trial.achievedLipidCounts, 'lower')}</span>}
          {!!trial.cleanupRemovedLipidCounts?.length && <span>Provider cleanup removed: {counts(trial.cleanupRemovedLipidCounts, 'upper')} / {counts(trial.cleanupRemovedLipidCounts, 'lower')}</span>}
          {!!trial.diagnosticArtifacts?.length && <details className="execution-diagnostics">
            <summary>Method-internal diagnostics ({trial.diagnosticArtifacts.length})</summary>
            <p className="help-text">These files belong to trial {trial.trialIndex + 1} ({label(trial.standing)}). They are provider working states, not a checked constructed system or completed stage.</p>
            <ul>{trial.diagnosticArtifacts.map(artifact => <li key={`${trial.trialId}:${artifact.role}`}>
              <strong>{diagnosticLabel(artifact.role)}</strong> · {diagnosticLabel(artifact.phase)}
              <div className="button-row">
                {artifact.subjectId && artifact.structureUrl && state.inspection?.subjectId !== artifact.subjectId && <button className="button compact"
                  type="button" disabled={inspectAction?.enabled !== true}
                  title={inspectAction?.reason ?? undefined}
                  onClick={() => onInspectSubject(artifact.subjectId!)}>Inspect diagnostic structure</button>}
                {artifact.downloadUrl && <a className="button compact" href={artifact.downloadUrl}
                  download={artifact.fileName}>Download diagnostic</a>}
              </div>
              <small>{artifact.fileName}</small>
            </li>)}</ul>
          </details>}
        </div>)}
      </details>}
    </section>

    <section className="account-card execution-basis-account" aria-label="Preparation basis">
      <h2>Preparation basis</h2>
      <dl className="detail-grid">
        <dt>Placement</dt><dd>{attempt.studyRevisionId && state.study?.id === attempt.studyRevisionId && state.placement?.status === 'supported'
          ? 'Supported position used for this run' : attempt.studyRevisionId ? 'Earlier accepted inputs' : 'Admission pending'}</dd>
        <dt>Chemistry</dt><dd>{actual ? 'Combined system represented and validated' : planned ? 'Provider-selected values recorded; validation pending' : 'Not yet established'}</dd>
      </dl>
      {planned && <details className="execution-candidate-values" aria-label="Provider-selected construction values">
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
        {planned.conditions && <details className="execution-exact-identity"><summary>Actual hydration and salt treatment</summary><dl className="detail-grid">
          <dt>Provider salt branch</dt><dd>{label(planned.conditions.saltBranch)}</dd>
          <dt>Nominal input</dt><dd>{planned.intendedNaClMolar.toFixed(3)} M</dd>
          <dt>Final Na⁺ / Cl⁻</dt><dd>{planned.conditions.finalSodiumCount} / {planned.conditions.finalChlorideCount}</dd>
          <dt>Retained Na⁺ / Cl⁻</dt><dd>{planned.conditions.retainedSodiumCount} / {planned.conditions.retainedChlorideCount}</dd>
          <dt>Provider Na⁺ / Cl⁻</dt><dd>{planned.conditions.providerGeneratedSodiumCount} / {planned.conditions.providerGeneratedChlorideCount}</dd>
          <dt>LEaP additions Na⁺ / Cl⁻</dt><dd>{planned.conditions.leapAddedSodiumCount} / {planned.conditions.leapAddedChlorideCount}</dd>
          <dt>LEaP removed generated Na⁺ / Cl⁻</dt><dd>{planned.conditions.leapRemovedGeneratedSodiumCount} / {planned.conditions.leapRemovedGeneratedChlorideCount}</dd>
          <dt>Water retained / generated / final</dt><dd>{planned.conditions.retainedWaterCount} / {planned.conditions.providerGeneratedWaterCount} / {planned.conditions.finalWaterCount}</dd>
          <dt>LEaP removed generated water</dt><dd>{planned.conditions.leapRemovedGeneratedWaterCount}</dd>
          <dt>Estimated aqueous volume</dt><dd>{planned.conditions.estimatedAqueousVolumeAngstromCubed.toLocaleString()} Å³</dd>
          <dt>Estimated aqueous Na⁺ / Cl⁻</dt><dd>{planned.conditions.sodiumAqueousMolar?.toFixed(4) ?? 'Unavailable'} / {planned.conditions.chlorideAqueousMolar?.toFixed(4) ?? 'Unavailable'} M</dd>
          <dt>Finite-water Na⁺ / Cl⁻</dt><dd>{planned.conditions.sodiumFiniteWaterMolar?.toFixed(4) ?? 'Unavailable'} / {planned.conditions.chlorideFiniteWaterMolar?.toFixed(4) ?? 'Unavailable'} M</dd>
        </dl><p className="help-text">Concentrations use estimated aqueous volume or the declared finite-water convention; neither is a measured bulk concentration.</p></details>}
        {planned.approximations.length > 0 && <div className="execution-approximations">
          <strong>Approximations</strong>
          {planned.approximations.map((item, index) => <p key={`${item}:${index}`}>{item}</p>)}
        </div>}
        {planned.limitations.map((item, index) => <p className="help-text" key={`${item}:${index}`}>Limit: {item}</p>)}
      </details>}
      {actual && <details className="execution-candidate-values execution-actual-values" aria-label="Validated actual constructed system">
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
          {actual.localState?.measurements && <MeasurementGroups values={actual.localState.measurements} />}
          {actual.localState?.limitations?.map((item, index) => <p key={`${item}:${index}`}>Limit: {item}</p>)}
        </details>
      </details>}
      {actual && <p className="help-text">The checked constructed system is an intermediate in the authorized Build and minimize operation. It is not a completed minimized result.</p>}
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
  const checkStanding = assessment ? label(assessment.checkStanding) : 'Not established';
  const originCurrent = stage.studyRevisionId === state.study?.id;
  const sourceStage = stage.kind === 'Equilibration' && stage.sourceStageId
    ? state.stages.find(item => item.stageId === stage.sourceStageId &&
      item.attemptId === stage.attemptId && item.kind === 'Minimization') : undefined;
  const finalForce = observation?.measurements.find(item => item.name === 'finalRmsForce' &&
    item.scope === 'final unrestrained constraint tangent; per particle');
  const completed = stage.status.toLowerCase() === 'completed';
  const addressedKinds = [...new Set((observation?.proteinGeometry?.locatedDistances ?? []).map(item => item.kind))];

  return <>
    <section className="account-card execution-identity-account" aria-label="Completed stage information">
      <h2>Stage information</h2>
      <dl className="detail-grid">
        <dt>Attempt</dt><dd>{stage.attemptId === state.attempt?.attemptId ? 'Current attempt (completed)' : 'Earlier accepted attempt'}</dd>
        {!originCurrent && <><dt>Inputs</dt><dd>This result uses earlier inputs.</dd></>}
        {runStartedLabel(stage.runStartedAt) && <><dt>Run</dt><dd>{runStartedLabel(stage.runStartedAt)}</dd></>}
        {stage.originMethodLabel && <><dt>Construction method</dt><dd>{stage.originMethodLabel}</dd></>}
        <dt>Constructed system</dt><dd>{actual ? 'Verified protein + bilayer + water and ions' : 'Historical construction account unavailable here'}</dd>
        <dt>Protein</dt><dd>{stage.originProteinLabel ?? (stage.studyRevisionId === state.study?.id
          ? state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId)?.label ?? 'Selected prepared protein'
          : 'Protein source unavailable in retained account')}</dd>
        <dt>Bilayer</dt><dd>{actual ? `${counts(actual.achievedComposition, 'upper')} / ${counts(actual.achievedComposition, 'lower')} physical leaflets` : 'Bound constructed bilayer'}</dd>
        <dt>Water and ions</dt><dd>{actual ? `${actual.waterCount.toLocaleString()} water · ${actual.sodiumCount.toLocaleString()} Na⁺ · ${actual.chlorideCount.toLocaleString()} Cl⁻` : 'Actual counts unavailable here'}</dd>
      </dl>
      {actual && <p className="help-text">{actual.atomCount.toLocaleString()} atoms · construction cell {cell(actual.cellAngstrom)}</p>}
    </section>

    <section className={`account-card execution-assessment-account${exportFault ? ' export-validation' : ''}`}
      aria-label={exportFault ? 'Export validation and unchanged stage review' :
        stage.kind === 'Minimization' ? 'Minimized stage review and technical checks' :
          'Equilibration stage review and technical checks'}>
      {exportFault ? <>
        <h2>Export validation</h2>
        <div className="export-failure" role="alert">
          <span className="export-failure-icon" aria-hidden="true">!</span>
          <div><strong>Export not delivered</strong><span>{exportFault}</span></div>
          <p>The {stage.kind === 'Minimization' ? 'minimized' : 'equilibrated'} stage remains completed, but the export could not be delivered. The stage and its technical checks are unchanged.</p>
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
        <dt>Technical checks</dt><dd>{assessment?.currentlyApplicable === false ? `Historical ${checkStanding}` : checkStanding}</dd>
        {exportFault && <><dt>Reason</dt><dd>{assessment?.reason ?? 'No stage-specific technical assessment is established.'}</dd></>}
      </dl>
      {!exportFault && <><p>{assessment?.reason ?? 'No stage-specific technical assessment is established.'}</p>
        <p className="help-text">{stage.kind === 'Minimization'
          ? 'Factual completion and the technical checks are recorded separately. This stage is not automatically equilibrated.'
          : 'Procedure completion and the technical checks are recorded separately. This stage does not establish global thermodynamic equilibrium.'}</p></>}
    </section>

    {(observation || assessment) &&
      <section className="account-card execution-findings-account" aria-label="Stage-specific findings and evidence">
        <h2>Stage-specific findings and evidence</h2>
        {observation && <details><summary>Observed final operation and measurements</summary>
          <p>Termination: {label(observation.termination)} · method version {observation.providerVersion}</p>
          <MeasurementGroups values={observation.measurements} />
          <p>Local-state observation: {observation.localState ? label(observation.localState.standing) : 'Unavailable'}{observation.localState?.unavailableReason && ` · ${observation.localState.unavailableReason}`}</p>
          {observation.localState?.rolePairMeasurements && observation.localState.rolePairMeasurements.length > 0 &&
            <details><summary>Contacts between molecular components ({observation.localState.rolePairMeasurements.length})</summary>
              <ul>{observation.localState.rolePairMeasurements.map((pair, index) => <li key={`${pair.firstMoleculeRole}:${pair.secondMoleculeRole}:${index}`}>
                {label(pair.firstMoleculeRole)}–{label(pair.secondMoleculeRole)}: {pair.pairsWithinSearchRadius.toLocaleString()} pairs within the search radius
                {pair.minimumDistanceAngstrom !== null && ` · nearest ${measured({ value: pair.minimumDistanceAngstrom, unit: 'Å' })}`}
              </li>)}</ul>
            </details>}
          {observation.localState?.measurements && <MeasurementGroups values={observation.localState.measurements} />}
          <p>Protein geometry coverage: {observation.proteinGeometry ? label(observation.proteinGeometry.standing) : 'Unavailable'}.</p>
          {observation.proteinGeometry?.kinds?.map(item => <details key={item.kind}>
            <summary>{measurementLabel(item.kind)} · {item.measuredCount.toLocaleString()} of {item.eligibleCount.toLocaleString()} eligible measurements</summary>
            <p>{label(item.standing)}{item.minimumDistanceAngstrom !== null && ` · shortest ${measured({ value: item.minimumDistanceAngstrom, unit: 'Å' })}`}{item.maximumDistanceAngstrom !== null && ` · longest ${measured({ value: item.maximumDistanceAngstrom, unit: 'Å' })}`}{item.unavailableReason && ` · ${item.unavailableReason}`}</p>
            <p className="help-text">Measurement coverage does not itself establish a passing technical check.</p>
          </details>)}
          {addressedKinds.map(kind => {
            const items = observation.proteinGeometry?.locatedDistances?.filter(item => item.kind === kind) ?? [];
            return <details className="geometry-address-list" key={kind}>
              <summary>{measurementLabel(kind)} · {items.length.toLocaleString()} addressed measurements</summary>
              <p className="help-text">Method: distance between identified atoms in the final stage coordinates. Scope: this completed stage.</p>
              <ul>{items.map((item, index) => <li key={`${kind}:${index}`}>
                <strong>{measured({ value: item.distanceAngstrom, unit: 'Å' })}</strong>
                <span>{atomAddress(item.first)} ↔ {atomAddress(item.second)}</span>
                {item.radiusSumAngstrom !== null && <span>Reference radii sum {measured({ value: item.radiusSumAngstrom, unit: 'Å' })}</span>}
              </li>)}</ul>
            </details>;
          })}
          {observation.proteinGeometry?.limitations.map((item, index) => <p key={`${item}:${index}`}>Limit: {item}</p>)}
        </details>}
        {assessment?.findings.map(finding => {
          const evidence = assessment.evidence.find(item => item.id === finding.evidenceId);
          return <div className="finding-item" key={finding.id}>
            <strong>{finding.consequence}</strong><span>{finding.meaning}</span>
            {evidence && <details><summary>Observation and method</summary>
              <p>{evidence.observation}</p><p>Method: {evidence.method}. Scope: {evidence.applicability}.</p>
              {evidence.uncertainty && <p>Limitation: {evidence.uncertainty}</p>}
            </details>}
          </div>;
        })}
        {assessment?.evidence.filter(evidence => !assessment.findings.some(finding => finding.evidenceId === evidence.id))
          .map(evidence => <div className="finding-item" key={evidence.id}>
            <strong>{evidence.observation}</strong>
            <details><summary>Method and scope</summary><p>Method: {evidence.method}. Scope: {evidence.applicability}.</p>
              {evidence.uncertainty && <p>Limitation: {evidence.uncertainty}</p>}
            </details>
          </div>)}
        {assessment?.limitations.map((limitation, index) => <p key={`${limitation}:${index}`}>Limit: {limitation}</p>)}
      </section>}
  </>;
}
