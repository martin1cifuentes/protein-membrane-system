/** Actor-facing wording for observed operation phases; exact phase codes remain in the account. */
export function operationLabel(phase: string | null | undefined, status?: string | null,
                               stopRequested = false): string {
  if (stopRequested && (status === 'running' || status === 'pending')) return 'Stop requested — awaiting observed stop';
  if (status === 'stopped') return 'Work stopped';
  if (status === 'unobserved') return 'Worker outcome unavailable';
  if (status?.toLowerCase() === 'resourcerefused') return 'Build stopped at the resource limit';
  if (status === 'failed') return 'Build and minimize failed';
  if (status === 'completed') return 'Build and minimize completed';
  switch (phase) {
    case 'admission': case 'providerPopulation': return 'Preparing inputs';
    case 'nativeConstruction': return 'Building the starting system';
    case 'providerPacking': return 'Packing — arranging molecules around the protein';
    case 'providerCleanup': return 'Checking and cleaning the packed system';
    case 'amberParameterization': return 'Preparing molecular parameters';
    case 'providerConditioningRestrained': case 'providerConditioningUnrestrained':
      return 'Conditioning the starting system';
    case 'handoffChecks': return 'Checking the constructed system';
    case 'finalMinimization': return 'Minimizing — adjusting atomic positions to reduce energy';
    case 'stageObservation': case 'assessment': return 'Checking the result';
    default: return 'Preparing the system';
  }
}

export function measurementLabel(name: string): string {
  const names: Record<string, string> = {
    covalentBond: 'Measured bond lengths', chainContinuity: 'Backbone connections',
    nonbondedDistance: 'Nonbonded distances',
    bondLength: 'Measured bond lengths', bondLengths: 'Measured bond lengths',
    backboneConnection: 'Backbone connections', backboneConnections: 'Backbone connections',
    atomsWithinMembraneGuide: 'Atoms within membrane guide',
    finalRmsForce: 'Final RMS force', maximumConstraintError: 'Maximum constraint error',
    minimumIntermolecularHeavyAtomDistance: 'Nearest intercomponent heavy atoms',
    lipidHeadSeparation: 'Leaflet head separation',
  };
  return names[name] ?? name.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/_/g, ' ');
}

export function runStartedLabel(value: string | null | undefined): string | null {
  if (!value) return null;
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? null :
    `Started ${new Intl.DateTimeFormat(undefined, { dateStyle: 'medium', timeStyle: 'short' }).format(date)}`;
}
