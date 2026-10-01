import { measurementLabel } from './executionDisplay';

export interface ProteinGeometryObservations {
  standing: string;
  kinds: {
    kind: string;
    standing: string;
    eligibleCount: number;
    measuredCount: number;
    minimumDistanceAngstrom: number | null;
    maximumDistanceAngstrom: number | null;
    unavailableReason: string | null;
  }[];
  locatedDistances?: {
    kind: string;
    first: { residue: { model: number; chain: string; copyId: string; residue: number; insertionCode: string }; atomName: string };
    second: { residue: { model: number; chain: string; copyId: string; residue: number; insertionCode: string }; atomName: string };
    distanceAngstrom: number;
    radiusSumAngstrom: number | null;
  }[];
  limitations: string[];
}

interface GeometryEvidence {
  id: string;
  source: string;
  uncertainty: string;
}

export function measuredDistance(value: number, addressed = false): string {
  const rounded = Number(value.toPrecision(4));
  const nearDisplayedBoundary = value !== rounded && Math.abs(value - rounded) < 0.0001;
  return `${value.toLocaleString(undefined, { maximumSignificantDigits: addressed || nearDisplayedBoundary ? 12 : 4 })} Å`;
}

function atomAddress(item: NonNullable<ProteinGeometryObservations['locatedDistances']>[number]['first'],
                     models: { index: number; sourceModelId?: string | null }[]): string {
  const residue = item.residue;
  const model = models.find(value => value.index === residue.model);
  const modelLabel = model?.sourceModelId ?? residue.model + 1;
  return `Model ${modelLabel} · chain ${residue.chain}${residue.copyId && residue.copyId !== residue.chain ? `, copy ${residue.copyId}` : ''} · residue ${residue.residue}${residue.insertionCode} · atom ${item.atomName}`;
}

function readable(value: string): string {
  const spaced = value.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/_/g, ' ');
  return spaced.charAt(0).toUpperCase() + spaced.slice(1);
}

function sourceLabel(source: string): string {
  return source === 'local scientific worker' ? 'Local structural measurement' : source;
}

/** One presentation of the exact geometry observations, with linked evidence details when supplied. */
export function GeometrySummary({ geometry, label, models = [], evidence = [], evidenceIds = [], scope }: {
  geometry: ProteinGeometryObservations | null | undefined;
  label: string;
  models?: { index: number; sourceModelId?: string | null }[];
  evidence?: GeometryEvidence[];
  evidenceIds?: string[];
  scope?: string;
}) {
  if (!geometry) return null;
  const addressed = new Map<string, NonNullable<ProteinGeometryObservations['locatedDistances']>>();
  for (const item of geometry.locatedDistances ?? [])
    addressed.set(item.kind, [...(addressed.get(item.kind) ?? []), item]);
  const linkedIds = new Set(evidenceIds);
  const linkedEvidence = evidence.filter(item => linkedIds.has(item.id));
  const sources = [...new Set(linkedEvidence.map(item => sourceLabel(item.source)))];
  const uncertainties = [...new Set(linkedEvidence.map(item => item.uncertainty).filter(Boolean))];
  const measurementScope = scope ?? label.toLowerCase();
  return <div className="hint-box">
    <strong>{label} · geometry measurements</strong>
    <p className="help-text">{label.toLowerCase().includes('source') ?
      'These source measurements alone do not establish structural suitability.' :
      'The measurements are observations; the protein outcome states the applicable assessment.'}</p>
    <details><summary>Measurement details</summary>
      <p className="help-text">Overall observation: {readable(geometry.standing)}.</p>
      <p className="help-text">Method: distance between the identified atom coordinates. Scope: {measurementScope}.</p>
      {sources.length > 0 && <p className="help-text">Source: {sources.join('; ')}.</p>}
      {geometry.kinds.map(item => <div key={item.kind} className="help-text">
        <strong>{measurementLabel(item.kind)}</strong> · {readable(item.standing)}
        <br />{item.measuredCount} of {item.eligibleCount} applicable distances measured
        {item.minimumDistanceAngstrom !== null && <> · shortest {measuredDistance(item.minimumDistanceAngstrom)}</>}
        {item.maximumDistanceAngstrom !== null && <> · longest {measuredDistance(item.maximumDistanceAngstrom)}</>}
        {item.unavailableReason && <> · {item.unavailableReason}</>}
      </div>)}
      {geometry.limitations.length > 0 && <p className="help-text">{geometry.limitations.join('; ')}</p>}
      {[...addressed].map(([kind, items]) => <details className="geometry-address-list" key={kind}>
        <summary>{measurementLabel(kind)} · {items.length.toLocaleString()} addressed measurements</summary>
        <ul>{items.map((item, index) => <li key={`${index}:${item.first.residue.model}:${item.first.residue.chain}:${item.first.residue.residue}:${item.first.atomName}`}>
          <strong>{measuredDistance(item.distanceAngstrom, true)}</strong>
          <span>{atomAddress(item.first, models)} ↔ {atomAddress(item.second, models)}</span>
          {item.radiusSumAngstrom !== null && <span>Reference radii sum {measuredDistance(item.radiusSumAngstrom, true)}</span>}
        </li>)}</ul>
      </details>)}
      {uncertainties.map(item => <p key={item} className="help-text">Limitation: {item}</p>)}
    </details>
  </div>;
}
