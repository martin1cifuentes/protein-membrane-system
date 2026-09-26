import { useCallback, useEffect, useRef, useState, type ReactNode } from 'react';
import { Viewer } from 'molstar/lib/apps/viewer/app';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { DefaultTrackballBindings } from 'molstar/lib/mol-canvas3d/controls/trackball';
import { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'molstar/lib/mol-geo/geometry/mesh/mesh-builder';
import { Box } from 'molstar/lib/mol-geo/primitive/box';
import { Mat4, Vec3 } from 'molstar/lib/mol-math/linear-algebra';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Shape } from 'molstar/lib/mol-model/shape';
import { ShapeRepresentation } from 'molstar/lib/mol-repr/shape/representation';
import { Color } from 'molstar/lib/mol-util/color';
import { Binding } from 'molstar/lib/mol-util/binding';
import { ButtonsType, ModifiersKeys } from 'molstar/lib/mol-util/input/input-observer';
import { StructureElement, Unit } from 'molstar/lib/mol-model/structure';
import type { PreparationAssessmentResult, ScientificEvidence, ScientificFinding, WorkspaceState } from './ProteinInMembraneWorkspace';
import { AttemptReviewAccount, StageReviewAccount } from './ExecutionReviewAccount';

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
  studyRevisionId: string;
  studyRevisionNumber: number;
  structureUrl: string | null;
  representationKind: string;
  omittedMolecules: string[];
  evidence: ScientificEvidence[];
  findings: ScientificFinding[];
  assessment: PreparationAssessmentResult | null;
  focusId: string | null;
  focus: StructureFocus | null;
  annotations: InspectionAnnotation[];
  metrics: InspectionMetric[];
}

interface AtomPickAccount {
  subjectId: string;
  studyRevisionId: string;
  structureToken: string;
  atomSiteIndex: number;
  atom: {
    resultAtomIndex: number;
    resultAtomId: string;
    sourceAtomId: string | null;
    role: string;
    moleculeRole: string;
    atomRole: string;
    element: string;
    sourceResidue: { model: number; chain: string; residue: number; insertionCode: string; copyId: string } | null;
    approvedChangeId: string | null;
    physicalSide: string | null;
    generatedSpeciesId: string | null;
    generatedComponentRole: string | null;
  };
}

function localStructureUrl(path: string): string {
  const url = new URL(path, window.location.origin);
  if (url.origin !== window.location.origin || !url.pathname.startsWith('/api/structures/')) {
    throw new Error('The selected structure is not available from this local workspace.');
  }
  return url.toString();
}

type PlacementAccount = NonNullable<WorkspaceState['placement']>;

// PPM's positioned protein and the measured membrane planes share one coordinate
// frame. These boxes are translucent extent guides; they contain no lipid atoms.
function focusPlacement(viewer: Viewer, placement: PlacementAccount) {
  const canvas = viewer.plugin.canvas3d;
  const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  const midplane = placement.midplaneAngstrom;
  const thickness = placement.thicknessAngstrom;
  if (!canvas || !structure || midplane === null || thickness === null || !Number.isFinite(midplane) ||
      !Number.isFinite(thickness) || thickness <= 0) {
    canvas?.requestCameraReset({ durationMs: 0 });
    return;
  }
  const bounds = structure.boundary.box;
  const cx = (bounds.min[0] + bounds.max[0]) / 2;
  const cy = (bounds.min[1] + bounds.max[1]) / 2;
  const width = Math.max(bounds.max[0] - bounds.min[0] + 22, 54);
  const compactScene = canvas.camera.viewport.width / canvas.camera.viewport.height < 0.9;
  const radius = Math.max(width * 0.72, thickness * 1.15, bounds.max[2] - bounds.min[2])
    * 0.58 * (compactScene ? 0.67 : 1);
  canvas.camera.setState(canvas.camera.getInvariantFocus(Vec3.create(cx, cy, midplane), radius,
    Vec3.create(0, 0, 1), Vec3.create(0, -1, 0)), 0);
}

// The membrane normal is the prepared protein's z axis. Mol*'s generic reset
// looks down that axis and turns a spanning protein into a small top-down blob.
function focusExplicitSystem(viewer: Viewer) {
  const canvas = viewer.plugin.canvas3d;
  const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!canvas || !structure) return;
  const bounds = structure.boundary.box;
  const center = Vec3.create(
    (bounds.min[0] + bounds.max[0]) / 2,
    (bounds.min[1] + bounds.max[1]) / 2,
    (bounds.min[2] + bounds.max[2]) / 2,
  );
  const aspect = Math.max(0.7, canvas.camera.viewport.width / canvas.camera.viewport.height);
  const xSpan = bounds.max[0] - bounds.min[0];
  const zSpan = bounds.max[2] - bounds.min[2];
  const heightRadius = Math.max(zSpan / 2, 20);
  // At reduced desktop widths, fitting every peripheral lipid across the
  // narrow canvas makes the protein and membrane core too small to inspect.
  // Keep the vertical system prominent and allow the patch edges to clip.
  const widthRadius = Math.min(xSpan / (2 * aspect), heightRadius * 1.08);
  const radius = Math.max(heightRadius, widthRadius) * 0.93;
  canvas.camera.setState(canvas.camera.getInvariantFocus(center, radius,
    Vec3.create(0, 0, 1), Vec3.create(0, -1, 0)), 0);
}

const sampledWaterTag = 'structure-component-explicit-water-sample';

// Keep complete water coordinates in the inspected structure. The sampled
// component only changes the initial depiction and has a disclosed rule.
async function addSampledWater(viewer: Viewer): Promise<boolean> {
  const root = viewer.plugin.managers.structure.hierarchy.current.structures[0];
  if (!root) return false;
  const water = await viewer.plugin.builders.structure.tryCreateComponentFromExpression(root.cell,
    MS.struct.generator.atomGroups({
      'residue-test': MS.core.logic.and([
        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'HOH']),
        MS.core.rel.eq([MS.core.math.mod([MS.struct.atomProperty.macromolecular.auth_seq_id(), 32]), 0]),
      ]),
    }), 'explicit-water-sample', { label: 'Displayed water sample' });
  if (!water) return false;
  await viewer.plugin.builders.structure.representation.addRepresentation(water, {
    type: 'ball-and-stick', typeParams: { sizeFactor: 0.28 }, color: 'element-symbol',
  });
  return true;
}

async function emphasizeIons(viewer: Viewer) {
  const root = viewer.plugin.managers.structure.hierarchy.current.structures[0];
  const ions = root?.components.find(component =>
    component.cell.transform.tags?.includes('structure-component-static-ion'));
  if (!ions) return;
  const original = [...ions.representations];
  await viewer.plugin.builders.structure.representation.addRepresentation(ions.cell, {
    type: 'spacefill', typeParams: { sizeFactor: 0.45 }, color: 'element-symbol',
  });
  viewer.plugin.managers.structure.hierarchy.toggleVisibility(original, 'hide');
}

async function emphasizeLipids(viewer: Viewer) {
  const root = viewer.plugin.managers.structure.hierarchy.current.structures[0];
  if (!root) return;
  // OpenMM's DMPC patch writes residue DMP. Mol* 5.11 recognizes DMPC as a
  // lipid name but classifies DMP as a generic ligand, whose chain color is
  // orange. Replace that exact DMP-only preset component in the view.
  const lipids = await viewer.plugin.builders.structure.tryCreateComponentFromExpression(root.cell,
    MS.struct.generator.atomGroups({
      'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'DMP']),
    }), 'explicit-dmp-lipids', { label: 'Displayed DMPC lipids' });
  if (!lipids) return;
  const count = lipids.cell?.obj?.data.elementCount;
  const preset = root.components.find(component =>
    (component.cell.transform.tags?.includes('structure-component-static-ligand') ||
      component.cell.transform.tags?.includes('structure-component-static-lipid')) &&
    component.cell.obj?.data.elementCount === count);
  if (!preset) return;
  await viewer.plugin.builders.structure.representation.addRepresentation(lipids, {
    type: 'ball-and-stick', typeParams: { sizeFactor: 0.15, alpha: 0.78 },
    color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } },
  });
  viewer.plugin.managers.structure.hierarchy.toggleVisibility([preset], 'hide');
}

interface BilayerPresentation { setVisible(visible: boolean): void; dispose(): void; }

async function addIntendedBilayer(viewer: Viewer, placement: PlacementAccount): Promise<BilayerPresentation | null> {
  const canvas = viewer.plugin.canvas3d;
  const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  const midplane = placement.midplaneAngstrom;
  const thickness = placement.thicknessAngstrom;
  if (!canvas || !structure || midplane === null || thickness === null || !Number.isFinite(midplane) ||
      !Number.isFinite(thickness) || thickness <= 0) return null;

  const bounds = structure.boundary.box;
  const cx = (bounds.min[0] + bounds.max[0]) / 2;
  const cy = (bounds.min[1] + bounds.max[1]) / 2;
  const width = Math.max(bounds.max[0] - bounds.min[0] + 22, 54) * 1.55;
  const depth = Math.max(bounds.max[1] - bounds.min[1] + 22, 48);
  const layerHeight = Math.max(1.4, Math.min(2.6, thickness * 0.08));
  const meshState = MeshBuilder.createState(192, 96);
  const box = Box();
  const addLayer = (group: number, z: number, height: number) => {
    meshState.currentGroup = group;
    const scale = Mat4.fromScaling(Mat4(), Vec3.create(width, depth, height));
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(cx, cy, z));
    MeshBuilder.addPrimitive(meshState, Mat4.mul(Mat4(), move, scale), box);
  };
  addLayer(0, midplane + thickness / 2, layerHeight);
  addLayer(1, midplane - thickness / 2, layerHeight);
  addLayer(2, midplane, Math.max(0.2, thickness - layerHeight));
  const mesh = MeshBuilder.getMesh(meshState);
  const colors = [Color(0x3d94b8), Color(0x688ebf), Color(0xb7d9e5)];
  const shape = Shape.create('Measured intended-bilayer extent', placement, mesh,
    group => colors[group] ?? colors[2], () => 1,
    group => group === 0 ? 'Upper physical leaflet target plane' :
      group === 1 ? 'Lower physical leaflet target plane' : 'Intended bilayer core; no packed lipids');
  const representation = ShapeRepresentation(() => shape, Mesh.Utils);
  await representation.createOrUpdate({ alpha: 0.42, doubleSided: true, ignoreLight: true }, placement).run();
  canvas.add(representation);
  focusPlacement(viewer, placement);
  let visible = true;
  return {
    setVisible(next) {
      if (next === visible) return;
      visible = next;
      if (next) canvas.add(representation);
      else canvas.remove(representation);
    },
    dispose() { if (visible) canvas.remove(representation); representation.destroy(); },
  };
}

function MolecularScene({ inspection, placement, executionReview = false, constructed, onPickAtom, onPickUnavailable }: {
  inspection: InspectionAccount | null;
  placement?: PlacementAccount | null;
  executionReview?: boolean;
  constructed: NonNullable<NonNullable<WorkspaceState['attempt']>['constructed']> | null;
  onPickAtom: (atomSiteIndex: number) => void;
  onPickUnavailable: (reason: string) => void;
}) {
  const mount = useRef<HTMLDivElement>(null);
  const [viewer, setViewer] = useState<Viewer | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [bilayerRegistered, setBilayerRegistered] = useState(false);
  const [toolMode, setToolMode] = useState<'select' | 'rotate' | 'pan' | 'zoom'>('select');
  const [showBilayer, setShowBilayer] = useState(true);
  const [displayOpen, setDisplayOpen] = useState(false);
  const [waterDisplay, setWaterDisplay] = useState<'sample' | 'all' | 'hidden'>('sample');
  const [showIons, setShowIons] = useState(true);
  const [displayAvailable, setDisplayAvailable] = useState({ protein: false, lipids: false, water: false, sampledWater: false, ions: false });
  const bilayer = useRef<BilayerPresentation | null>(null);
  const structureUrl = inspection?.structureUrl ?? null;
  const fullSystemReview = executionReview && (inspection?.representationKind === 'constructedSystem'
    || inspection?.representationKind === 'completedStage');
  const hasSelectedFocus = useRef(false);
  hasSelectedFocus.current = !!inspection?.focus;

  useEffect(() => {
    if (!structureUrl || !mount.current) return;
    let closed = false;
    const active = { viewer: null as Viewer | null };
    let releaseBilayer: BilayerPresentation | null = null;
    const target = mount.current;
    setError(null);
    setViewer(null);
    setBilayerRegistered(false);
    setDisplayAvailable({ protein: false, lipids: false, water: false, sampledWater: false, ions: false });

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
          viewportShowSelectionMode: false,
          viewportShowAnimation: false,
          viewportBackgroundColor: '#edf3f6',
        });
        if (closed) { created.dispose(); return; }
        active.viewer = created;
        // Keep the ordinary Mol* instance attached to its mount for local
        // inspection integrations; it confers no scientific authority.
        Object.defineProperty(target, Symbol.for('molstar.viewer'),
          { value: created, configurable: true });
        if (executionReview) created.plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset,
          'preset-structure-representation-polymer-and-ligand');
        const format = new URL(url).searchParams.get('format');
        if (format !== 'pdb' && format !== 'mmcif')
          throw new Error('The selected structure has no recognized PDB or mmCIF representation.');
        await created.loadStructureFromUrl(url, format, false, { label: inspection?.subjectId });
        if (closed) return;
        // Mol* may report a failed URL load in its own log without rejecting
        // the promise. A renderer without an actual structure is unavailable.
        const loaded = created.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
        if (!loaded || loaded.units.length === 0)
          throw new Error('The identified structure could not be loaded.');
        if (fullSystemReview) {
          const components = created.plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [];
          const hasComponent = (kind: string) => components.some(component =>
            component.cell.transform.tags?.includes(`structure-component-static-${kind}`));
          const hasWater = hasComponent('water');
          const fullWater = components.filter(component =>
            component.cell.transform.tags?.includes('structure-component-static-water'));
          created.plugin.managers.structure.hierarchy.toggleVisibility(fullWater, 'hide');
          let sampledWater = false;
          if (hasWater) {
            try { sampledWater = await addSampledWater(created); }
            catch { /* The complete water renderer remains available. */ }
          }
          try { await emphasizeLipids(created); }
          catch { /* The native lipid representation remains visible. */ }
          try { await emphasizeIons(created); }
          catch { /* The native ion representation remains visible. */ }
          if (closed) return;
          setWaterDisplay(sampledWater ? 'sample' : 'hidden');
          setDisplayAvailable({ protein: hasComponent('polymer'),
            lipids: hasComponent('lipid') || hasComponent('ligand'),
            water: hasWater, sampledWater, ions: hasComponent('ion') });
          focusExplicitSystem(created);
        }
        if (placement) {
          const added = await addIntendedBilayer(created, placement);
          if (closed) { added?.dispose(); return; }
          releaseBilayer = added;
          bilayer.current = releaseBilayer;
          setBilayerRegistered(releaseBilayer !== null);
        }
        if (!closed) setViewer(created);
      } catch (cause) {
        if (!closed) {
          if (Reflect.get(target, Symbol.for('molstar.viewer')) === active.viewer)
            Reflect.deleteProperty(target, Symbol.for('molstar.viewer'));
          active.viewer?.dispose();
          active.viewer = null;
          setError(cause instanceof Error ? cause.message : 'The structure could not be rendered.');
        }
      }
    }

    void load();
    return () => {
      closed = true;
      releaseBilayer?.dispose();
      bilayer.current = null;
      if (Reflect.get(target, Symbol.for('molstar.viewer')) === active.viewer)
        Reflect.deleteProperty(target, Symbol.for('molstar.viewer'));
      active.viewer?.dispose();
    };
  }, [structureUrl, placement?.proposalId, placement?.midplaneAngstrom, placement?.thicknessAngstrom,
    executionReview, fullSystemReview]);

  useEffect(() => { bilayer.current?.setVisible(showBilayer); }, [showBilayer, viewer]);

  useEffect(() => {
    if (!viewer || !fullSystemReview) return;
    const components = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [];
    const setComponentVisibility = (tag: string, visible: boolean) => {
      const matching = components.filter(component =>
        component.cell.transform.tags?.includes(tag));
      viewer.plugin.managers.structure.hierarchy.toggleVisibility(matching, visible ? 'show' : 'hide');
    };
    setComponentVisibility('structure-component-static-water', waterDisplay === 'all');
    setComponentVisibility(sampledWaterTag, waterDisplay === 'sample');
    setComponentVisibility('structure-component-static-ion', showIons);
  }, [viewer, fullSystemReview, waterDisplay, showIons]);

  useEffect(() => {
    if (!viewer || (!placement && !executionReview)) return;
    const canvas = viewer.plugin.canvas3d;
    if (!canvas) return;
    const primary = Binding([Binding.Trigger(ButtonsType.Flag.Primary, ModifiersKeys.create())]);
    canvas.setAttribs({ trackball: { bindings: {
      ...DefaultTrackballBindings,
      dragRotate: toolMode === 'pan' || toolMode === 'zoom' ? Binding.Empty : DefaultTrackballBindings.dragRotate,
      dragPan: toolMode === 'pan' ? primary : DefaultTrackballBindings.dragPan,
      dragZoom: toolMode === 'zoom' ? primary : DefaultTrackballBindings.dragZoom,
    } } });
    viewer.plugin.selectionMode = toolMode === 'select';
  }, [viewer, placement, executionReview, toolMode]);

  useEffect(() => {
    if (!viewer || toolMode !== 'select' || !inspection?.structureUrl) return;
    const subscription = viewer.plugin.behaviors.interaction.click.subscribe(({ current }) => {
      const loci = current.loci;
      // Mol* publishes its initial empty BehaviorSubject value on subscription.
      // That value is not a researcher selection and must not create a card.
      if (loci.kind === 'empty-loci' ||
          StructureElement.Loci.is(loci) && StructureElement.Loci.isEmpty(loci)) return;
      if (!StructureElement.Loci.is(loci) || StructureElement.Loci.size(loci) !== 1) {
        onPickUnavailable('Select one visible atom to inspect its verified molecular identity.');
        return;
      }
      const location = StructureElement.Loci.getFirstLocation(loci);
      if (!location || !Unit.isAtomic(location.unit)) {
        onPickUnavailable('An atom-level identity is unavailable for this part of the view.');
        return;
      }
      const atomSiteIndex = location.unit.model.atomicHierarchy.atomSourceIndex.value(location.element);
      if (!Number.isInteger(atomSiteIndex) || atomSiteIndex < 0) {
        onPickUnavailable('The selected atom has no verified coordinate-row identity.');
        return;
      }
      onPickAtom(atomSiteIndex);
    });
    return () => subscription.unsubscribe();
  }, [viewer, toolMode, inspection?.structureUrl, onPickAtom, onPickUnavailable]);

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
    let frame = 0;
    let generation = 0;
    let disposed = false;
    const fitWhenViewportMatches = (current: number, attempts: number) => {
      if (disposed || current !== generation) return;
      const bounds = mount.current?.getBoundingClientRect();
      const viewport = viewer.plugin.canvas3d?.camera.viewport;
      // Mol* updates its camera viewport after handleResize. Fitting against
      // the preceding narrow viewport leaves a newly widened system tiny.
      if (bounds && viewport &&
          Math.abs(viewport.width - bounds.width) <= 4 &&
          Math.abs(viewport.height - bounds.height) <= 4) {
        if (!hasSelectedFocus.current) {
          if (placement) focusPlacement(viewer, placement);
          else if (fullSystemReview) focusExplicitSystem(viewer);
          else viewer.plugin.canvas3d?.requestCameraReset({ durationMs: 0 });
        }
        mount.current?.setAttribute('data-camera-viewport-width', String(viewport.width));
        mount.current?.setAttribute('data-camera-viewport-height', String(viewport.height));
        mount.current?.setAttribute('data-camera-ready', 'true');
      } else if (attempts < 30) {
        frame = window.requestAnimationFrame(() => fitWhenViewportMatches(current, attempts + 1));
      }
    };
    const observer = new ResizeObserver(() => {
      generation += 1;
      mount.current?.setAttribute('data-camera-ready', 'false');
      window.cancelAnimationFrame(frame);
      viewer.handleResize();
      const current = generation;
      frame = window.requestAnimationFrame(() => fitWhenViewportMatches(current, 0));
    });
    observer.observe(mount.current);
    return () => { disposed = true; window.cancelAnimationFrame(frame); observer.disconnect(); };
  }, [viewer, placement, fullSystemReview]);

  if (!structureUrl) {
    return <div className="scene-empty"><strong>{inspection ? 'No renderable structure is established' : 'No molecular subject selected'}</strong><span>{inspection ? 'The selected subject’s established values and limitations remain available beside this view.' : 'Choose a source or completed subject to inspect its available structure.'}</span></div>;
  }

  return <>
    <div className="viewer-mount" ref={mount} aria-label={`3D structure for ${inspection?.subjectId ?? 'selected subject'}`} />
    {(placement || executionReview) && <nav className="placement-view-tools" aria-label="Molecular view tools">
      {(['select', 'rotate', 'pan', 'zoom'] as const).map(mode =>
        <button type="button" key={mode} className={toolMode === mode ? 'active' : ''}
          aria-pressed={toolMode === mode} onClick={() => setToolMode(mode)}
          title={`${mode[0].toUpperCase() + mode.slice(1)} in the molecular view`}>
          <span aria-hidden="true">{{ select: '↖', rotate: '↻', pan: '✣', zoom: '⌕' }[mode]}</span>
          {mode[0].toUpperCase() + mode.slice(1)}
        </button>)}
      <button type="button" disabled={!inspection?.metrics.length}
        onClick={() => document.querySelector('.placement-metrics, .execution-metrics')?.scrollIntoView({ block: 'start' })}
        title="Show measured and derived values"><span aria-hidden="true">⌁</span>Measure</button>
      <button type="button" onClick={() => viewer && (placement ? focusPlacement(viewer, placement) : fullSystemReview ? focusExplicitSystem(viewer) : viewer.plugin.canvas3d?.requestCameraReset({ durationMs: 0 }))}
        title={placement ? 'Focus the positioned protein and intended bilayer' : 'Fit the identified system in the molecular view'}><span aria-hidden="true">◎</span>Focus</button>
      {placement && <button type="button" className={showBilayer ? 'active' : ''} aria-pressed={showBilayer}
        onClick={() => setShowBilayer(value => !value)} title="Show or hide the intended bilayer bounds">
        <span aria-hidden="true">◉</span>Display
      </button>}
      {fullSystemReview && <div className="execution-display-tool">
        <button type="button" className={displayOpen ? 'active' : ''} aria-expanded={displayOpen}
          aria-controls="execution-display-options" onClick={() => setDisplayOpen(value => !value)}
          title="Choose which explicit molecule classes are rendered">
          <span aria-hidden="true">◉</span>Display
        </button>
        {displayOpen && <div id="execution-display-options" className="execution-display-options"
          aria-label="Rendered explicit molecule classes">
          <strong>Displayed in the molecular view</strong>
          <p>Protein: {viewer ? displayAvailable.protein ? 'shown' : 'renderer unavailable' : 'loading'} · lipids: {viewer ? displayAvailable.lipids ? 'shown' : 'renderer unavailable' : 'loading'}.</p>
          <fieldset className="execution-water-options" disabled={!displayAvailable.water}>
            <legend>Water molecules</legend>
            <label><input type="radio" name="execution-water" checked={waterDisplay === 'sample'}
              disabled={!displayAvailable.sampledWater}
              onChange={() => setWaterDisplay('sample')} />Representative sample</label>
            <label><input type="radio" name="execution-water" checked={waterDisplay === 'all'}
              onChange={() => setWaterDisplay('all')} />All</label>
            <label><input type="radio" name="execution-water" checked={waterDisplay === 'hidden'}
              onChange={() => setWaterDisplay('hidden')} />Hidden</label>
          </fieldset>
          <label><input type="checkbox" checked={showIons} disabled={!displayAvailable.ions}
            onChange={event => setShowIons(event.target.checked)} />Na⁺ / Cl⁻ ions</label>
          {!displayAvailable.sampledWater && displayAvailable.water &&
            <small>Representative water display is unavailable; all waters can still be shown.</small>}
          {constructed && <small>Full model: {constructed.atomCount.toLocaleString()} atoms; {constructed.waterCount.toLocaleString()} waters; {constructed.sodiumCount.toLocaleString()} Na⁺ and {constructed.chlorideCount.toLocaleString()} Cl⁻.</small>}
          {(!displayAvailable.protein || !displayAvailable.lipids ||
            (!!constructed && constructed.waterCount > 0 && !displayAvailable.water) ||
            (!!constructed && constructed.sodiumCount + constructed.chlorideCount > 0 && !displayAvailable.ions)) && viewer &&
            <small>Any molecule class without a recognized renderer is unavailable in this view; its coordinates remain in the identified structure.</small>}
        </div>}
      </div>}
    </nav>}
    {fullSystemReview && viewer &&
      <div className="execution-display-disclosure" role="status">
        Protein {displayAvailable.protein ? 'shown' : 'not rendered'} · lipids {displayAvailable.lipids ? 'shown' : 'not rendered'} · water {constructed?.waterCount === 0 ? 'none' :
          !displayAvailable.water ? 'not rendered' : waterDisplay === 'sample' && displayAvailable.sampledWater
            ? 'sampled (residue IDs divisible by 32)' : waterDisplay === 'all' ? 'all shown' : 'hidden'} · ions {
              constructed && constructed.sodiumCount + constructed.chlorideCount === 0 ? 'none' :
                !displayAvailable.ions ? 'not rendered' : showIons ? 'shown' : 'hidden'}.
        {' '}The identified all-atom coordinates and counts are unchanged.
      </div>}
    {placement && <div className="placement-scene-legend" role="img"
      aria-label={bilayerRegistered && placement.midplaneAngstrom !== null && placement.thicknessAngstrom !== null
        ? `Intended bilayer around positioned protein. Upper and lower physical leaflet planes registered to measured midplane ${placement.midplaneAngstrom.toFixed(2)} angstrom and thickness ${placement.thicknessAngstrom.toFixed(2)} angstrom. No packed lipid positions or achieved membrane are shown.`
        : 'Intended bilayer and positioned protein. No measured membrane plane is available to register in this view; no packed lipid positions or achieved membrane are shown.'}>
      <span><i className="placement-upper-swatch" /> Upper physical leaflet</span>
      <span><i className="placement-lower-swatch" /> Lower physical leaflet</span>
      <strong>{bilayerRegistered && placement.midplaneAngstrom !== null && placement.thicknessAngstrom !== null
        ? 'Measured target planes · no packed lipids' : 'Bilayer registration unavailable'}</strong>
    </div>}
    {!viewer && !error && <div className="scene-loading" role="status">Loading the identified structure…</div>}
    {error && <div className="scene-error" role="alert"><strong>Structure unavailable</strong><span>{error}</span></div>}
  </>;
}

function statusClass(status: string): string {
  const value = status.toLowerCase().replace(/[^a-z]/g, '');
  if (value.includes('fail') || value.includes('declined') || value.includes('unsupported') || value.includes('notsupported') || value.includes('notqualified') || value.includes('unqualified') || value.includes('disqualified') || value.includes('refused')) return 'danger';
  if (value.includes('notestablished') || value.includes('notcurrent') || value.includes('indeterminate') || value.includes('unavailable')) return 'warning';
  if (value.includes('assessed') || value.includes('qualified') || value.includes('supported') || value.includes('completed')) return 'success';
  return '';
}

function representationLabel(kind: string): string {
  switch (kind) {
    case 'structuralSource': return 'Structural source';
    case 'preparedProtein': return 'Prepared protein';
    case 'selected-protein-before-repair': return 'Selected protein before repair';
    case 'membraneModel': return 'Intended planar bilayer';
    case 'oriented-protein-with-proposed-membrane-bounds': return 'Oriented protein with intended bilayer';
    case 'constructedSystem': return 'Verified constructed explicit system';
    case 'completedStage': return 'Completed molecular stage';
    default: return kind;
  }
}

type MembraneAccount = NonNullable<WorkspaceState['membrane']>;
type LipidFraction = MembraneAccount['upper'][number];

function presentFractions(fractions: LipidFraction[]) {
  return fractions.filter(item => item.fraction > 0);
}

function targetPercent(fraction: number): string {
  return `${String(fraction * 100)}%`;
}

function IntendedBilayerPreview({ membrane }: { membrane: MembraneAccount }) {
  const upper = presentFractions(membrane.upper);
  const lower = presentFractions(membrane.lower);
  const speciesIds = [...new Set([...upper, ...lower].map(item => item.speciesId))];
  const color = (speciesId: string) => {
    const hue = (206 + speciesIds.indexOf(speciesId) * 137.5) % 360;
    return `hsl(${hue} 48% 60%)`;
  };
  const describe = (items: LipidFraction[]) => items
    .map(item => `${item.speciesId} ${targetPercent(item.fraction)}`).join(', ');
  const band = (items: LipidFraction[], side: 'upper' | 'lower') => <div className={`bilayer-leaflet-band ${side}`} aria-hidden="true">
    {items.map(item => <div key={item.speciesId} className="bilayer-species-segment"
      style={{ width: `${item.fraction * 100}%`, backgroundColor: color(item.speciesId) }} />)}
  </div>;

  return <figure className="bilayer-preview">
    <div className="bilayer-preview-title">
      <span>Intended planar bilayer</span>
      <strong>Composition preview</strong>
    </div>
    <div className="bilayer-preview-graphic" role="img"
      aria-label={`Intended bilayer preview. Upper physical leaflet target fractions: ${describe(upper)}. Lower physical leaflet target fractions: ${describe(lower)}. This is a composition diagram, not placed molecular coordinates or achieved packing.`}>
      <div className="bilayer-side-label">Upper aqueous side</div>
      <div className="bilayer-physical-label">Upper physical leaflet</div>
      {band(upper, 'upper')}
      <div className="bilayer-core"><span>Bilayer core · no molecule positions shown</span></div>
      {band(lower, 'lower')}
      <div className="bilayer-physical-label">Lower physical leaflet</div>
      <div className="bilayer-side-label">Lower aqueous side</div>
    </div>
    <div className="bilayer-preview-legend" aria-label="Selected lipid and sterol species">
      {speciesIds.map(speciesId => <span className="bilayer-legend-item" key={speciesId}>
        <span className="bilayer-legend-swatch" style={{ backgroundColor: color(speciesId) }} aria-hidden="true" />
        {speciesId}
      </span>)}
    </div>
    <figcaption>Intended fractions only. No lipid coordinates, finite molecule counts, packing result, or biological sidedness is established by this diagram.</figcaption>
  </figure>;
}

function MembraneLeafletAccount({ side, fractions, state }: {
  side: 'Upper' | 'Lower';
  fractions: LipidFraction[];
  state: WorkspaceState;
}) {
  const present = presentFractions(fractions);
  return <section className="account-card membrane-leaflet-account" aria-label={`${side} physical leaflet`}>
    <h2>{side} physical leaflet</h2>
    <dl className="detail-grid">
      {present.map(item => {
        const species = state.availableLipids.find(candidate => candidate.speciesId === item.speciesId);
        return <div className="membrane-fraction-row" key={item.speciesId}>
          <dt>{species?.displayName ?? item.speciesId} <span className="tabular">({item.speciesId})</span></dt>
          <dd>{targetPercent(item.fraction)} target</dd>
        </div>;
      })}
    </dl>
    <p className="help-text">Physical {side.toLowerCase()} is a coordinate label, not an established biological side. Zero-fraction entries are absent.</p>
  </section>;
}

function PlacementReviewAccount({ placement, state, inspection, canSelectFocus, onSelectFocus }: {
  placement: PlacementAccount;
  state: WorkspaceState;
  inspection: InspectionAccount;
  canSelectFocus: boolean;
  onSelectFocus: (annotationId: string | null) => void;
}) {
  const measured = placement.midplaneAngstrom !== null && placement.thicknessAngstrom !== null;
  const standing = placement.status === 'supported' ? 'Supported placement' :
    placement.status === 'unsupported' ? 'Unsupported placement' :
      placement.status === 'proposed' ? 'Proposal under assessment' : 'Support not established';
  const candidateBasis = [...placement.evidence].reverse().find(item =>
    /candidate|adjustment/i.test(item.method)) ?? placement.evidence[0];
  const selectedSource = state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId);
  const sourceLabel = selectedSource?.label ?? (state.study?.selectedSourceKind === 'upload'
    ? 'Researcher-supplied structural source' : 'Selected structural source');
  const chainSummary = state.study?.chainIds.length ? state.study.chainIds.join(', ') : 'Not established';
  const upperTarget = state.membrane?.modelId === placement.membraneModelId
    ? presentFractions(state.membrane.upper).map(item => `${item.speciesId} ${targetPercent(item.fraction)}`).join(', ')
    : 'Not established';
  const lowerTarget = state.membrane?.modelId === placement.membraneModelId
    ? presentFractions(state.membrane.lower).map(item => `${item.speciesId} ${targetPercent(item.fraction)}`).join(', ')
    : 'Not established';
  return <>
    <section className="account-card placement-identity-account" aria-label="Selected protein and membrane model">
      <h2>Selected protein / model</h2>
      <dl className="detail-grid">
        <dt>Type</dt><dd>{placement.topologyKind === 'membrane-spanning' ? 'Membrane-spanning protein' : 'Surface-associated protein'}</dd>
        <dt>Source</dt><dd>{sourceLabel}</dd>
        <dt>Model</dt><dd>{state.study?.modelIndex !== null && state.study?.modelIndex !== undefined
          ? `Source model ${state.study.modelIndex + 1} · chain ${chainSummary}` : `Chain ${chainSummary}`}</dd>
        <dt>Representation</dt><dd>Positioned cartoon + intended bounds</dd>
        <dt>Environment</dt><dd>{upperTarget === lowerTarget ? `${upperTarget} planar bilayer target` : 'Intended planar bilayer target'}</dd>
      </dl>
      <details className="placement-identity-details"><summary>Exact subject identifiers and leaflet targets</summary>
        <dl className="detail-grid">
          <dt>Prepared protein ID</dt><dd className="tabular">{placement.preparedProteinId ?? 'Not established'}</dd>
          <dt>Membrane model ID</dt><dd className="tabular">{placement.membraneModelId ?? 'Not established'}</dd>
          {state.study?.selectedSourceId && <><dt>Structural source ID</dt><dd className="tabular">{state.study.selectedSourceId}</dd></>}
          {state.study?.uploadProvenance && <><dt>Declared origin</dt><dd>{state.study.uploadProvenance} (researcher declared)</dd></>}
          <dt>Upper physical leaflet</dt><dd>{upperTarget}</dd>
          <dt>Lower physical leaflet</dt><dd>{lowerTarget}</dd>
        </dl>
      </details>
    </section>
    <section className="account-card placement-proposal-account" aria-label="Placement proposal">
      <h2>Placement proposal</h2>
      <dl className="detail-grid">
        <dt>Relationship</dt><dd>{placement.topologyKind === 'membrane-spanning' ? 'Membrane-spanning' : 'One-surface associated'}</dd>
        <dt>Physical contact side</dt><dd>{placement.physicalSide === 'both' ? 'Both physical leaflets' : `${placement.physicalSide} physical side`}</dd>
        {candidateBasis && <><dt>Candidate basis</dt><dd>{candidateBasis.method}</dd></>}
        {placement.depthAngstrom !== null && <><dt>Depth</dt><dd>{placement.depthAngstrom.toFixed(2)} Å</dd></>}
        {placement.tiltDegrees !== null && <><dt>Tilt</dt><dd>{placement.tiltDegrees.toFixed(1)}°</dd></>}
        {measured && <><dt>Midplane</dt><dd>{placement.midplaneAngstrom!.toFixed(2)} Å</dd>
          <dt>Bilayer thickness</dt><dd>{placement.thicknessAngstrom!.toFixed(2)} Å</dd></>}
      </dl>
      <details className="placement-identity-details"><summary>Exact proposal and candidate evidence</summary>
        <dl className="detail-grid">
          <dt>Proposal ID</dt><dd className="tabular">{placement.proposalId}</dd>
          {candidateBasis && <><dt>Candidate source</dt><dd>{candidateBasis.source}</dd>
            <dt>Method</dt><dd>{candidateBasis.method}</dd>
            <dt>Observed basis</dt><dd>{candidateBasis.observation}</dd></>}
        </dl>
        <p>Upper and lower are physical coordinate sides. The translucent planes mark intended bounds; they show no lipid positions or achieved packing.</p>
      </details>
    </section>
    <section className="account-card placement-regions-account" aria-label="Orientation evidence and contacting regions">
      <h2>Orientation evidence</h2>
      <div className="placement-orientation-list">
        {(['upper side', 'membrane core', 'lower side'] as const).map(label => {
          const anchor = inspection.annotations.find(item => item.label === label);
          return <button type="button" key={label} disabled={!anchor?.geometryFocus || !inspection.structureUrl || !canSelectFocus}
            className={anchor && inspection.focusId === anchor.subjectPartId ? 'selected' : ''}
            onClick={() => anchor && onSelectFocus(anchor.id)}>
            <span>{label === 'membrane core' ? 'Membrane core' : `${label === 'upper side' ? 'Upper' : 'Lower'} physical side`}</span>
            <strong>{anchor?.geometryFocus && inspection.structureUrl ? `Observed · ${anchor.subjectPartId}` : 'No verified spatial anchor'}</strong>
            {anchor?.evidenceId && <small className="tabular">Evidence {anchor.evidenceId}</small>}
          </button>;
        })}
      </div>
      {placement.contactingRegions.length > 0
        ? <details className="placement-contact-details"><summary>{placement.contactingRegions.length} exact contacting region{placement.contactingRegions.length === 1 ? '' : 's'} · inspect list</summary>
          <ul>{placement.contactingRegions.map((region, index) => <li key={`${region}:${index}`}>{region}</li>)}</ul>
        </details>
        : <p>Contacting regions have not been independently established for this proposal.</p>}
      <dl className="detail-grid">
        <dt>Biological premise</dt><dd>{placement.sidedness || 'Not supplied'}</dd>
      </dl>
      <p className="help-text">The biological premise and PPM in/out assignment are distinct from physical upper and lower coordinates.</p>
    </section>
    <section className={`account-card placement-support-account ${statusClass(placement.status)}`} aria-label="Placement standing and uncertainty">
      <h2>Evidence and uncertainty</h2>
      <p className="placement-standing">{standing}</p>
      <p>{placement.reason}</p>
      {(placement.policyId || placement.witnessId) && <dl className="detail-grid">
        {placement.policyId && <><dt>Policy</dt><dd>{placement.policyId}{placement.policyVersion && ` · version ${placement.policyVersion}`}</dd></>}
        {placement.witnessId && <><dt>Independent witness</dt><dd>{placement.witnessId}</dd></>}
      </dl>}
      {placement.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
    </section>
    {placement.evidence.length > 0 && <section className="account-card placement-evidence-account" aria-label="Attributed placement evidence">
      <h2>Attributed evidence</h2>
      {placement.evidence.map(item => <div className="finding-item" key={item.id}>
        <strong>{item.bearing} · {item.source} · {item.method}</strong>
        <span>{item.observation}</span>
        <span>Applicability: {item.applicability}</span>
        <span>Uncertainty: {item.uncertainty}</span>
      </div>)}
    </section>}
  </>;
}

export function ConnectedStructuralInspection({
  state,
  onSelectFocus,
  onInspectSubject,
  proposalDecision,
  exportFault = null,
  reviewAttempt = false,
  connectionMessage = null,
  onRefreshAccount,
}: {
  state: WorkspaceState;
  onSelectFocus: (annotationId: string | null) => void;
  onInspectSubject: (subjectId: string) => void;
  proposalDecision?: ReactNode;
  exportFault?: string | null;
  reviewAttempt?: boolean;
  connectionMessage?: string | null;
  onRefreshAccount?: () => void;
}) {
  const inspection = state.inspection;
  const subject = inspection?.subjectId ?? state.study?.id ?? null;
  const stage = state.stages.find(item => item.stageId === subject);
  const preparationChange = state.protein?.changes.find(change => change.id === subject);
  const membrane = state.membrane?.modelId === subject ? state.membrane : null;
  const placement = state.placement?.proposalId === subject ? state.placement : null;
  const executionReview = reviewAttempt || !!stage;
  const showSubjectEvidence = !reviewAttempt || inspection?.representationKind === 'constructedSystem';
  const subjectStatus = stage?.assessment && !stage.assessment.currentlyApplicable ? 'Assessment not current'
    : stage?.assessment?.qualification ?? stage?.status
    ?? placement?.status
    ?? (state.protein?.subjectId === subject ? state.protein.status : null)
    ?? (preparationChange ? state.protein?.status === 'declined' ? 'Declined proposal' : 'Preparation proposal' : null)
    ?? (membrane ? membrane.status === 'notEstablished' ? 'Not established' : membrane.status === 'assessed' ? 'Assessed' : 'Proposed' : null);
  const selectedAction = state.actions.find(item => item.kind === 'setInspectionFocus' && item.subjectId === null);
  const canSelectFocus = selectedAction?.enabled === true;
  const affectedAnnotation = preparationChange && inspection?.annotations.find(item => item.geometryFocus);
  const [pickedAtom, setPickedAtom] = useState<AtomPickAccount | null>(null);
  const [pickMessage, setPickMessage] = useState<string | null>(null);
  const pickSequence = useRef(0);
  const selectedStructureKey = `${inspection?.subjectId ?? ''}|${inspection?.structureUrl ?? ''}`;
  const selectedStructureKeyRef = useRef(selectedStructureKey);
  selectedStructureKeyRef.current = selectedStructureKey;
  const visiblePickedAtom = pickedAtom && inspection?.subjectId === pickedAtom.subjectId &&
    inspection.structureUrl?.split('?')[0] === `/api/structures/${pickedAtom.structureToken}` ? pickedAtom : null;
  const constructed = stage?.constructed ?? state.attempt?.constructed ?? null;
  const evidenceContent = useRef<HTMLDivElement>(null);
  const [placementTab, setPlacementTab] = useState<'evidence' | 'models' | 'notes'>('evidence');
  const compactReview = useRef<boolean | null>(null);
  const currentSourceContext = inspection && state.study?.selectedSourceId && (
    subject === state.study.selectedSourceId || subject === state.protein?.subjectId || preparationChange);

  useEffect(() => {
    pickSequence.current += 1;
    setPickedAtom(null);
    setPickMessage(null);
  }, [selectedStructureKey]);

  const onPickUnavailable = useCallback((reason: string) => {
    pickSequence.current += 1;
    setPickedAtom(null);
    setPickMessage(reason);
  }, []);

  const onPickAtom = useCallback((atomSiteIndex: number) => {
    if (!inspection?.structureUrl) {
      onPickUnavailable('The selected subject has no verified structure for atom inspection.');
      return;
    }
    let token: string;
    try {
      const pathname = new URL(localStructureUrl(inspection.structureUrl)).pathname;
      const match = /^\/api\/structures\/([A-Za-z0-9._-]+)$/.exec(pathname);
      if (!match) throw new Error('No exact structure token');
      token = match[1];
    } catch {
      onPickUnavailable('The selected structure has no verified local identity.');
      return;
    }
    const requestIndex = ++pickSequence.current;
    const subjectId = inspection.subjectId;
    const revisionId = inspection.studyRevisionId;
    const selectedKey = selectedStructureKey;
    setPickedAtom(null);
    setPickMessage('Resolving selected atom against the host correspondence…');
    void (async () => {
      try {
        const response = await fetch(`/api/inspection/atoms/${encodeURIComponent(subjectId)}/${encodeURIComponent(token)}/${atomSiteIndex}`, { cache: 'no-store' });
        if (requestIndex !== pickSequence.current || selectedStructureKeyRef.current !== selectedKey) return;
        if (!response.ok) {
          setPickMessage(response.status === 404
            ? 'Verified atom identity is unavailable for this selected subject and structure.'
            : `The selected atom could not be checked against the local account (${response.status}).`);
          return;
        }
        const account = await response.json() as AtomPickAccount;
        if (requestIndex !== pickSequence.current || selectedStructureKeyRef.current !== selectedKey) return;
        if (account.subjectId !== subjectId || account.studyRevisionId !== revisionId ||
            account.structureToken !== token || account.atomSiteIndex !== atomSiteIndex ||
            account.atom?.resultAtomIndex !== atomSiteIndex || !account.atom.resultAtomId) {
          setPickMessage('The selected atom response does not match this exact structure and revision.');
          return;
        }
        setPickedAtom(account);
        setPickMessage(null);
      } catch {
        if (requestIndex === pickSequence.current && selectedStructureKeyRef.current === selectedKey)
          setPickMessage('The selected atom identity is unavailable while the local connection is interrupted.');
      }
    })();
  }, [inspection?.subjectId, inspection?.structureUrl, inspection?.studyRevisionId, selectedStructureKey, onPickUnavailable]);

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
    <section className="scene-panel" aria-label={stage ? 'Selected completed molecular stage' : reviewAttempt ? 'Current attempt inspection' : membrane ? 'Intended membrane model' : placement ? 'Placed protein and intended bilayer' : 'Molecular structure'}>
      <div className="scene-heading">
        <div>
          <p className="eyebrow">Connected structural inspection</p>
          <h1>{stage?.summary ?? (inspection && representationLabel(inspection.representationKind)) ?? 'Scientific subject'}</h1>
          <div className="scene-subtitle tabular">{subject ?? 'No subject selected'}{inspection &&
            ` · ${representationLabel(inspection.representationKind)} · study revision ${inspection.studyRevisionNumber} (${inspection.studyRevisionId}) · ${inspection.studyRevisionId === state.study?.id ? 'current' : 'historical'}`}</div>
        </div>
        {subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      <div className={`scene-surface${placement && !reviewAttempt ? ' placement-scene-surface' : ''}${executionReview ? ' execution-scene-surface' : ''}`}>
        {membrane && !reviewAttempt ? <IntendedBilayerPreview membrane={membrane} /> : <MolecularScene inspection={inspection} placement={reviewAttempt ? null : placement}
          executionReview={executionReview} constructed={constructed} onPickAtom={onPickAtom} onPickUnavailable={onPickUnavailable} />}
        {executionReview && <div className="execution-scene-label">{stage
          ? `Selected ${stage.kind === 'Minimization' ? 'minimized' : 'equilibrated'} stage${stage.studyRevisionId === state.study?.id ? '' : ' · historical revision'}`
          : inspection?.representationKind === 'constructedSystem' ? 'Current constructed system · verified'
            : 'Current inspected subject · constructed system view not established'}</div>}
      </div>
      <div className="scene-caption">
        {executionReview ? <><strong>{stage ? 'Selected completed molecular stage.' : inspection?.representationKind === 'constructedSystem' ? 'Verified constructed starting system.' : 'Earlier inspected subject.'}</strong>{' '}
          {inspection?.omittedMolecules.length ? `Not rendered: ${inspection.omittedMolecules.join(', ')}.` : 'See the exact subject account for representation limits and measurements.'}
          {' '}The scene itself does not establish a scientific qualification.</>
          : placement ? <><strong>Oriented protein with intended bilayer bounds.</strong> The upper and lower translucent planes are registered to the measured midplane and thickness where available. PPM uses an implicit symmetric DOPC membrane; the selected membrane is assessed separately. No packed lipids or achieved system are shown.</>
          : membrane ? <><strong>Intended model, not achieved packing.</strong> The two physical leaflets show target fractions only; membrane-local support is stated in the evidence account.</> : <><strong>Spatial view is evidence, not assessment.</strong>{' '}
        {inspection?.omittedMolecules.length
          ? `Not shown: ${inspection.omittedMolecules.join(', ')}.`
          : 'Any available omissions and exact measurements appear with this subject’s account.'}</>}
      </div>
    </section>

    <aside className="evidence-panel" aria-label="Evidence and assessment">
      {connectionMessage && <div className="notice warning" role="alert">{connectionMessage}
        {onRefreshAccount && <div className="button-row"><button className="button compact" type="button" onClick={onRefreshAccount}>Refresh account</button></div>}
      </div>}
      {(placement || executionReview) && <nav className="placement-account-tabs" aria-label="Review account sections">
        {(['evidence', 'models', 'notes'] as const).map(tab => <button type="button" key={tab}
          className={placementTab === tab ? 'active' : ''} aria-current={placementTab === tab ? 'page' : undefined}
          onClick={() => {
            setPlacementTab(tab);
            const content = evidenceContent.current;
            const target = tab === 'models' ? content?.querySelector('.execution-identity-account, .placement-identity-account') :
              tab === 'notes' ? content?.querySelector('.execution-assessment-account, .execution-basis-account, .placement-support-account') : null;
            if (!content) return;
            content.scrollTop = target
              ? target.getBoundingClientRect().top - content.getBoundingClientRect().top + content.scrollTop - 4
              : 0;
          }}>{tab[0].toUpperCase() + tab.slice(1)}</button>)}
      </nav>}
      <div className="evidence-content" ref={evidenceContent}>
      <p className="eyebrow">Exact subject account</p>
      <div className="status-line">
        <h2 className="panel-heading">Evidence &amp; standing</h2>
        {subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      {!inspection && !reviewAttempt && <div className="hint-box">Select a scientific subject to connect its structure with the applicable values and findings.</div>}
      {reviewAttempt && state.attempt && <AttemptReviewAccount state={state} attempt={state.attempt} onInspectConstructed={onInspectSubject} />}
      {inspection && !stage && <section className="account-card inspection-origin-account" aria-label="Selected subject and study revision">
        <h2>Selected subject and study revision</h2>
        <dl className="detail-grid">
          <dt>Subject</dt><dd className="tabular">{inspection.subjectId}</dd>
          <dt>Origin</dt><dd>Study revision {inspection.studyRevisionNumber} <span className="tabular">({inspection.studyRevisionId})</span></dd>
          <dt>Context</dt><dd>{inspection.studyRevisionId === state.study?.id ? 'Current study revision' :
            `Historical subject; current study is revision ${state.study?.number ?? 'unknown'} (${state.study?.id ?? 'unavailable'})`}</dd>
          {inspection.assessment && <><dt>Stage assessment</dt><dd>{inspection.assessment.currentlyApplicable ? '' : 'Historical assessment · '}{inspection.assessment.qualification}: {inspection.assessment.reason}</dd></>}
        </dl>
      </section>}
      {inspection && showSubjectEvidence && <>
        {!preparationChange && !placement && !executionReview && <section className="account-card" aria-label={membrane ? 'Selected membrane model' : 'Selected structure'}>
          <h2>{membrane ? 'Selected membrane model' : 'Selected structure'}</h2>
          <dl className="detail-grid">
          <dt>Subject</dt><dd>{inspection.subjectId}</dd>
          <dt>Representation</dt><dd>{representationLabel(inspection.representationKind)}</dd>
          {membrane && <><dt>Standing</dt><dd>{membrane.status === 'notEstablished' ? 'Not established' : membrane.status === 'assessed' ? 'Assessed membrane model' : 'Proposed intention'}</dd></>}
          {inspection.focusId && <><dt>Selected part</dt><dd>{inspection.focusId}</dd></>}
          </dl>
        </section>}
        {membrane && !executionReview && <>
          <MembraneLeafletAccount side="Upper" fractions={membrane.upper} state={state} />
          <MembraneLeafletAccount side="Lower" fractions={membrane.lower} state={state} />
          <section className="account-card membrane-purpose-account" aria-label="Membrane purpose and study conditions">
            <h2>Purpose and fixed study conditions</h2>
            <p>{membrane.scientificPurpose}</p>
            <dl className="detail-grid">
              <dt>Target salt and neutrality</dt><dd>{state.study ? `${state.study.conditions.targetNaClMolar} M NaCl with charge-neutralizing compatible monovalent counterions` : 'Not established'}</dd>
              <dt>Optional equilibration target</dt><dd>{state.study ? `${state.study.conditions.optionalTemperatureKelvin} K` : 'Not established'}</dd>
              {state.study && <><dt>Nominal pH</dt><dd>{state.study.conditions.nominalPh}</dd></>}
            </dl>
            <p className="help-text">These disclosed targets are not achieved conditions or controls for this intended model.</p>
          </section>
          <section className={`account-card membrane-support-account ${membrane.status === 'assessed' ? 'established' : membrane.status === 'notEstablished' ? 'unestablished' : ''}`} aria-label="Membrane-local support and uncertainty">
            <h2>Membrane-local support and uncertainty</h2>
            <p className="membrane-support-standing">{membrane.status === 'assessed' ? 'Assessed membrane model' : membrane.status === 'notEstablished' ? 'Membrane model not established' : 'Proposal not yet assessed'}</p>
            {membrane.reason && <p className="membrane-support-reason">{membrane.reason}</p>}
            {membrane.policyId && <dl className="detail-grid"><dt>Policy</dt><dd>{membrane.policyId}{membrane.policyVersion && ` · version ${membrane.policyVersion}`}</dd></dl>}
            {membrane.status === 'proposed' && <p>The selected fractions are an editable proposal until deliberately adopted and assessed.</p>}
            {membrane.status === 'assessed' && <p>Support applies to this exact model and the identified molecular representations. Protein placement and an assembled membrane remain separate.</p>}
            {membrane.limitations.map(limit => <p className="membrane-limit" key={limit}>Limit: {limit}</p>)}
          </section>
          {membrane.speciesSupport?.length > 0 && <section className="account-card membrane-assets-account" aria-label="Selected species representation support">
            <h2>Selected species and representations</h2>
            {membrane.speciesSupport.map(species => <div className="membrane-species-support" key={species.speciesId}>
              <strong>{species.speciesId} · {species.chemistryId}</strong>
              <span>{species.category} · {species.forceFieldFamily} {species.forceFieldVersion}</span>
              <details><summary>Coordinate and parameter asset identities</summary>
                <dl className="detail-grid">
                  <dt>Starting-coordinate SHA-256</dt><dd className="tabular">{species.coordinateSha256}</dd>
                  <dt>Parameter SHA-256</dt><dd className="tabular">{species.parameterSha256}</dd>
                </dl>
              </details>
              {species.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
            </div>)}
          </section>}
          {membrane.evidence?.length > 0 && <section className="account-card membrane-evidence-account" aria-label="Membrane assessment evidence">
            <h2>Assessment evidence</h2>
            {membrane.evidence.map(item => <div className="finding-item" key={item.id}>
              <strong>{item.source} · {item.method}</strong>
              <span>{item.observation}</span>
              <span>Applicability: {item.applicability}</span>
              <span>Uncertainty: {item.uncertainty}</span>
            </div>)}
          </section>}
        </>}
        {placement && !reviewAttempt && <PlacementReviewAccount placement={placement} state={state}
          inspection={inspection} canSelectFocus={canSelectFocus} onSelectFocus={onSelectFocus} />}
        {stage && <StageReviewAccount state={state} stage={stage} exportFault={exportFault} />}
        {(visiblePickedAtom || pickMessage) && <section className="account-card inspection-atom-account" aria-label="Selected atom correspondence">
          <h2>Selected atom correspondence</h2>
          {visiblePickedAtom ? <>
            <dl className="detail-grid">
              <dt>Result atom</dt><dd className="tabular">{visiblePickedAtom.atom.resultAtomId}</dd>
              <dt>Role</dt><dd>{visiblePickedAtom.atom.moleculeRole} · {visiblePickedAtom.atom.atomRole} · {visiblePickedAtom.atom.element}</dd>
              <dt>Origin</dt><dd>{visiblePickedAtom.atom.sourceAtomId ? `Source atom ${visiblePickedAtom.atom.sourceAtomId}` :
                `Generated ${visiblePickedAtom.atom.generatedComponentRole ?? visiblePickedAtom.atom.role}`}</dd>
              {visiblePickedAtom.atom.sourceResidue && <><dt>Source residue</dt><dd>Model {visiblePickedAtom.atom.sourceResidue.model} · chain {visiblePickedAtom.atom.sourceResidue.chain} · residue {visiblePickedAtom.atom.sourceResidue.residue}{visiblePickedAtom.atom.sourceResidue.insertionCode} · copy {visiblePickedAtom.atom.sourceResidue.copyId}</dd></>}
              {visiblePickedAtom.atom.generatedSpeciesId && <><dt>Species</dt><dd>{visiblePickedAtom.atom.generatedSpeciesId}</dd></>}
              {visiblePickedAtom.atom.physicalSide && <><dt>Physical leaflet</dt><dd>{visiblePickedAtom.atom.physicalSide}</dd></>}
              {visiblePickedAtom.atom.approvedChangeId && <><dt>Approved change</dt><dd>{visiblePickedAtom.atom.approvedChangeId}</dd></>}
            </dl>
            <p className="help-text">Identity is resolved against this exact host-bound coordinate row and correspondence. Selecting it does not change the molecular model or assessment.</p>
          </> : <p className="help-text">{pickMessage}</p>}
        </section>}
        {(inspection.findings.length > 0 || inspection.evidence.length > 0) &&
          <section className="evidence-section inspection-connected-account" aria-label="Selected subject findings and evidence">
            <h2>Selected subject findings and evidence</h2>
            {inspection.findings.map(finding => {
              const basis = inspection.evidence.find(item => item.id === finding.evidenceId);
              return <div className="finding-item" key={finding.id}>
                <strong>{finding.disposition} · {finding.consequence}</strong>
                <span>{finding.meaning}</span>
                <span className="tabular">Finding {finding.id} · Subject {finding.subjectId}</span>
                {basis ? <><span>Evidence: {basis.source} · {basis.method}</span><span>{basis.observation}</span>
                  <span>Applicability: {basis.applicability}</span><span>Uncertainty: {basis.uncertainty}</span></>
                  : <span>Linked evidence unavailable for this finding.</span>}
              </div>;
            })}
            {inspection.evidence.filter(item => !inspection.findings.some(finding => finding.evidenceId === item.id))
              .map(item => <div className="finding-item" key={item.id}>
                <strong>{item.source} · {item.method}</strong><span>{item.observation}</span>
                <span>Applicability: {item.applicability}</span><span>Uncertainty: {item.uncertainty}</span>
                <span className="tabular">Evidence {item.id} · Subject {item.subjectId}</span>
              </div>)}
          </section>}
        {currentSourceContext && !executionReview && <section className="account-card source-account" aria-label="Source and assembly">
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
            {affectedAnnotation && <div className="button-row"><button className="button compact" type="button" disabled={!canSelectFocus || !inspection.structureUrl}
              title={!canSelectFocus ? selectedAction?.reason ?? 'Selection is not available.' : !inspection.structureUrl ? 'The exact structure is unavailable for spatial focus.' : undefined}
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
              disabled={!canSelectFocus || !inspection.structureUrl || !annotation.geometryFocus}
              title={!canSelectFocus ? selectedAction?.reason ?? 'Selection is not available.' :
                !inspection.structureUrl || !annotation.geometryFocus ? 'No verified atom-level spatial focus is available for this evidence.' : undefined}
              onClick={() => onSelectFocus(annotation.id)}
            >
              <span className="item-title">{annotation.label}</span>
              <span className="item-detail">{annotation.meaning}</span>
              <span className="item-detail tabular">{annotation.subjectPartId}{annotation.evidenceId && ` · Evidence ${annotation.evidenceId}`}</span>
              {(!annotation.geometryFocus || !inspection.structureUrl) && <span className="item-detail">Spatial focus unavailable; evidence remains inspectable.</span>}
            </button>)}
          </div>
        </section>}
        {inspection.metrics.length > 0 && <section className={`evidence-section${placement ? ' placement-metrics' : ''}${executionReview ? ' execution-metrics' : ''}`}>
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
