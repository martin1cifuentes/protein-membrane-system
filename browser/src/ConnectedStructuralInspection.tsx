import { useCallback, useEffect, useRef, useState, type ReactNode } from 'react';
import { Viewer } from 'molstar/lib/apps/viewer/app';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { DefaultTrackballBindings } from 'molstar/lib/mol-canvas3d/controls/trackball';
import type { Camera } from 'molstar/lib/mol-canvas3d/camera';
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
import { StructureElement, StructureProperties, Unit } from 'molstar/lib/mol-model/structure';
import { StructureFocusRepresentation } from 'molstar/lib/mol-plugin/behavior/dynamic/selection/structure-focus-representation';
import { clearStructureTransparency, setStructureTransparency } from 'molstar/lib/mol-plugin-state/helpers/structure-transparency';
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

export interface StructureLoadStatus {
  subjectId: string;
  structureUrl: string;
  phase: 'loading' | 'displayed' | 'failed';
  reason?: string;
}

interface PickedLocation {
  atomSiteIndex: number;
  residueName: string | null;
  chain: string;
  residueNumber: number | null;
  insertionCode: string;
  modelNumber: number;
  entityKind: string;
  atomName: string;
  element: string;
  atomLevel: boolean;
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

async function emphasizeLipids(viewer: Viewer): Promise<boolean> {
  const root = viewer.plugin.managers.structure.hierarchy.current.structures[0];
  if (!root) return false;
  // OpenMM's DMPC patch writes residue DMP. Mol* 5.11 recognizes DMPC as a
  // lipid name but classifies DMP as a generic ligand, whose chain color is
  // orange. Replace that exact DMP-only preset component in the view.
  const lipids = await viewer.plugin.builders.structure.tryCreateComponentFromExpression(root.cell,
    MS.struct.generator.atomGroups({
      'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'DMP']),
    }), 'explicit-dmp-lipids', { label: 'Displayed DMPC lipids' });
  if (!lipids) return false;
  const count = lipids.cell?.obj?.data.elementCount;
  const preset = root.components.find(component =>
    (component.cell.transform.tags?.includes('structure-component-static-ligand') ||
      component.cell.transform.tags?.includes('structure-component-static-lipid')) &&
    component.cell.obj?.data.elementCount === count);
  if (!preset) return false;
  await viewer.plugin.builders.structure.representation.addRepresentation(lipids, {
    type: 'ball-and-stick', typeParams: { sizeFactor: 0.15, alpha: 0.78 },
    color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } },
  });
  viewer.plugin.managers.structure.hierarchy.toggleVisibility([preset], 'hide');
  return true;
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

type MolecularRepresentation = 'cartoon' | 'sticks' | 'surface' | 'spacefill';
const rememberedCamera = new Map<string, Camera.Snapshot>();
const rememberedDisplay = new Map<string, {
  representation: MolecularRepresentation; protein: boolean; lipids: boolean; ligands: boolean;
  water: 'sample' | 'all' | 'hidden'; ions: boolean;
}>();
const representationType: Record<MolecularRepresentation, 'cartoon' | 'ball-and-stick' | 'molecular-surface' | 'spacefill'> = {
  cartoon: 'cartoon', sticks: 'ball-and-stick', surface: 'molecular-surface', spacefill: 'spacefill',
};

function MolecularScene({ inspection, structureLabel, placement, executionReview = false, constructed, onPickAtom, onPickUnavailable,
  onStructureLoad, reloadToken, localInspection, hasSelection, onClearSelection, onShowWholeSystem, clearSelectionSerial,
  wholeSystemSerial, onChainColors, chainFocus }: {
  inspection: InspectionAccount | null;
  structureLabel: string;
  placement?: PlacementAccount | null;
  executionReview?: boolean;
  constructed: NonNullable<NonNullable<WorkspaceState['attempt']>['constructed']> | null;
  onPickAtom: (location: PickedLocation) => void;
  onPickUnavailable: (reason: string) => void;
  onStructureLoad: (status: StructureLoadStatus) => void;
  reloadToken: number;
  localInspection: boolean;
  hasSelection: boolean;
  onClearSelection: () => void;
  onShowWholeSystem: () => void;
  clearSelectionSerial: number;
  wholeSystemSerial: number;
  onChainColors?: (subjectId: string, structureUrl: string, colors: Record<string, string>) => void;
  chainFocus?: { chainId: string; serial: number } | null;
}) {
  const mount = useRef<HTMLDivElement>(null);
  const [viewer, setViewer] = useState<Viewer | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [bilayerRegistered, setBilayerRegistered] = useState(false);
  const [toolMode, setToolMode] = useState<'select' | 'rotate' | 'pan' | 'zoom'>('select');
  const [showBilayer, setShowBilayer] = useState(true);
  const [displayOpen, setDisplayOpen] = useState(false);
  const [representation, setRepresentation] = useState<MolecularRepresentation>('cartoon');
  const [showProtein, setShowProtein] = useState(true);
  const [showLipids, setShowLipids] = useState(true);
  const [showLigands, setShowLigands] = useState(true);
  const [waterDisplay, setWaterDisplay] = useState<'sample' | 'all' | 'hidden'>('sample');
  const [showIons, setShowIons] = useState(true);
  const [showContacts, setShowContacts] = useState(false);
  const [contextCartoon, setContextCartoon] = useState<'faded' | 'shown' | 'hidden'>('faded');
  const [displayAvailable, setDisplayAvailable] = useState({ protein: false, chains: [] as string[], lipids: false,
    ligands: false, water: false, sampledWater: false, ions: false });
  const [displayError, setDisplayError] = useState<string | null>(null);
  const [pickedLoci, setPickedLoci] = useState<StructureElement.Loci | null>(null);
  const [annotationLoci, setAnnotationLoci] = useState<StructureElement.Loci | null>(null);
  const representationGeneration = useRef(0);
  const bilayer = useRef<BilayerPresentation | null>(null);
  const structureUrl = inspection?.structureUrl ?? null;
  const displayMemory = useRef({ representation, protein: showProtein, lipids: showLipids,
    ligands: showLigands, water: waterDisplay, ions: showIons });
  displayMemory.current = { representation, protein: showProtein, lipids: showLipids,
    ligands: showLigands, water: waterDisplay, ions: showIons };
  const fullSystemReview = executionReview && (inspection?.representationKind === 'constructedSystem'
    || inspection?.representationKind === 'completedStage');

  useEffect(() => {
    const saved = structureUrl ? rememberedDisplay.get(structureUrl) : undefined;
    setRepresentation(saved?.representation ?? 'cartoon');
    setShowProtein(saved?.protein ?? true);
    setShowLipids(saved?.lipids ?? true);
    setShowLigands(saved?.ligands ?? true);
    setShowIons(saved?.ions ?? true);
    setWaterDisplay(saved?.water ?? 'hidden');
    setShowContacts(false);
    setContextCartoon('faded');
    setPickedLoci(null);
    setAnnotationLoci(null);
  }, [structureUrl]);

  useEffect(() => () => {
    if (structureUrl) rememberedDisplay.set(structureUrl, displayMemory.current);
  }, [structureUrl]);

  useEffect(() => {
    if (!structureUrl || !mount.current) return;
    const loadingUrl = structureUrl;
    let closed = false;
    const active = { viewer: null as Viewer | null };
    let releaseBilayer: BilayerPresentation | null = null;
    const target = mount.current;
    setError(null);
    setViewer(null);
    setBilayerRegistered(false);
    setDisplayAvailable({ protein: false, chains: [], lipids: false, ligands: false,
      water: false, sampledWater: false, ions: false });
    onStructureLoad({ subjectId: inspection!.subjectId, structureUrl, phase: 'loading' });

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
        created.plugin.selectionMode = true;
        // Keep the ordinary Mol* instance attached to its mount for local
        // inspection integrations; it confers no scientific authority.
        Object.defineProperty(target, Symbol.for('molstar.viewer'),
          { value: created, configurable: true });
        // The automatic Mol* preset switches small proteins to an all-atom
        // component, leaving no polymer component for the documented
        // representation controls. Keep one explicit component structure.
        created.plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset,
          'preset-structure-representation-polymer-and-ligand');
        const format = new URL(url).searchParams.get('format');
        if (format !== 'pdb' && format !== 'mmcif')
          throw new Error('The selected structure has no recognized PDB or mmCIF representation.');
        await created.loadStructureFromUrl(url, format, false, { label: structureLabel });
        if (closed) return;
        // Mol* may report a failed URL load in its own log without rejecting
        // the promise. A renderer without an actual structure is unavailable.
        const loaded = created.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
        if (!loaded || loaded.units.length === 0)
          throw new Error('The identified structure could not be loaded.');
        let sampledWater = false;
        let customLipids = false;
        if (fullSystemReview) {
          const components = created.plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [];
          const hasComponent = (kind: string) => components.some(component =>
            component.cell.transform.tags?.includes(`structure-component-static-${kind}`));
          const hasWater = hasComponent('water');
          const fullWater = components.filter(component =>
            component.cell.transform.tags?.includes('structure-component-static-water'));
          created.plugin.managers.structure.hierarchy.toggleVisibility(fullWater, 'hide');
          if (hasWater) {
            try { sampledWater = await addSampledWater(created); }
            catch { /* The complete water renderer remains available. */ }
          }
          try { customLipids = await emphasizeLipids(created); }
          catch { /* The native lipid representation remains visible. */ }
          try { await emphasizeIons(created); }
          catch { /* The native ion representation remains visible. */ }
          if (closed) return;
          setWaterDisplay(rememberedDisplay.get(structureUrl!)?.water ?? (sampledWater ? 'sample' : 'hidden'));
          focusExplicitSystem(created);
        } else {
          setWaterDisplay(rememberedDisplay.get(structureUrl!)?.water ?? 'hidden');
        }
        const components = created.plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [];
        const hasComponent = (kind: string) => components.some(component =>
          component.cell.transform.tags?.includes(`structure-component-static-${kind}`));
        const polymer = components.find(component =>
          component.cell.transform.tags?.includes('structure-component-static-polymer'));
        const chains = [...new Set((polymer?.cell.obj?.data.units ?? []).filter(Unit.isAtomic).map(unit =>
          StructureProperties.chain.auth_asym_id(StructureElement.Location.create(polymer!.cell.obj!.data, unit, unit.elements[0]))))];
        const cartoon = polymer?.representations.find(item => item.cell.params?.values?.type?.name === 'cartoon');
        const colorTheme = cartoon?.cell.obj?.data.repr.theme.color;
        const displayedColors: Record<string, string> = {};
        if (polymer?.cell.obj?.data && colorTheme && 'color' in colorTheme) {
          for (const unit of polymer.cell.obj.data.units.filter(Unit.isAtomic)) {
            const location = StructureElement.Location.create(polymer.cell.obj.data, unit, unit.elements[0]);
            const chain = StructureProperties.chain.auth_asym_id(location);
            const tint = Color.toHexString(colorTheme.color(location, false));
            displayedColors[`${StructureProperties.unit.model_num(location)}:${chain}`] = tint;
            if (loaded.models.length === 1) displayedColors[chain] = tint;
          }
        }
        onChainColors?.(inspection!.subjectId, loadingUrl, displayedColors);
        setDisplayAvailable({ protein: !!polymer, chains, lipids: customLipids || hasComponent('lipid'),
          ligands: components.some(component => component.cell.transform.tags?.includes('structure-component-static-ligand') &&
            !component.cell.state.isHidden), water: hasComponent('water'), sampledWater,
          ions: hasComponent('ion') });
        if (placement) {
          const added = await addIntendedBilayer(created, placement);
          if (closed) { added?.dispose(); return; }
          releaseBilayer = added;
          bilayer.current = releaseBilayer;
          setBilayerRegistered(releaseBilayer !== null);
        }
        if (!closed) {
          setViewer(created);
          onStructureLoad({ subjectId: inspection!.subjectId, structureUrl: loadingUrl, phase: 'displayed' });
        }
      } catch (cause) {
        if (!closed) {
          if (Reflect.get(target, Symbol.for('molstar.viewer')) === active.viewer)
            Reflect.deleteProperty(target, Symbol.for('molstar.viewer'));
          active.viewer?.dispose();
          active.viewer = null;
          const reason = cause instanceof Error ? cause.message : 'The structure could not be rendered.';
          setError(reason);
          onStructureLoad({ subjectId: inspection!.subjectId, structureUrl: loadingUrl, phase: 'failed', reason });
        }
      }
    }

    void load();
    return () => {
      closed = true;
      if (active.viewer?.plugin.canvas3d && structureUrl)
        rememberedCamera.set(structureUrl, active.viewer.plugin.canvas3d.camera.getSnapshot());
      releaseBilayer?.dispose();
      bilayer.current = null;
      if (Reflect.get(target, Symbol.for('molstar.viewer')) === active.viewer)
        Reflect.deleteProperty(target, Symbol.for('molstar.viewer'));
      active.viewer?.dispose();
    };
  }, [structureUrl, placement?.proposalId, placement?.midplaneAngstrom, placement?.thicknessAngstrom,
    executionReview, fullSystemReview, reloadToken]);

  useEffect(() => {
    if (!viewer || !structureUrl) return;
    const saved = rememberedCamera.get(structureUrl);
    if (!saved) return;
    let inner = 0;
    const outer = window.requestAnimationFrame(() => {
      inner = window.requestAnimationFrame(() => {
        viewer.handleResize();
        viewer.plugin.canvas3d?.camera.setState(saved, 0);
      });
    });
    return () => { window.cancelAnimationFrame(outer); window.cancelAnimationFrame(inner); };
  }, [viewer, structureUrl]);

  useEffect(() => {
    if (!viewer || !chainFocus) return;
    const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
    if (!structure) return;
    const loci = StructureElement.Loci.fromSchema(structure, { auth_asym_id: chainFocus.chainId });
    if (StructureElement.Loci.isEmpty(loci)) return;
    viewer.plugin.managers.interactivity.lociSelects.selectOnly({ loci }, false);
    viewer.plugin.managers.camera.focusLoci(loci, { extraRadius: 8, durationMs: 0 });
  }, [viewer, chainFocus?.serial]);

  useEffect(() => { bilayer.current?.setVisible(showBilayer); }, [showBilayer, viewer]);

  useEffect(() => {
    if (!viewer) return;
    const generation = ++representationGeneration.current;
    const apply = async () => {
      const root = viewer.plugin.managers.structure.hierarchy.current.structures[0];
      const polymer = root?.components.find(component =>
        component.cell.transform.tags?.includes('structure-component-static-polymer'));
      if (!polymer) return;
      const type = representationType[representation];
      let selected = polymer.representations.find(item => item.cell.params?.values?.type?.name === type);
      try {
        if (!selected) {
          await viewer.plugin.builders.structure.representation.addRepresentation(polymer.cell, {
            type, color: representation === 'sticks' || representation === 'spacefill' ? 'element-symbol' : 'chain-id',
          });
          selected = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.components.find(component =>
            component.cell.transform.tags?.includes('structure-component-static-polymer'))?.representations
            .find(item => item.cell.params?.values?.type?.name === type);
        }
        if (generation !== representationGeneration.current || !selected) return;
        const all = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.components.find(component =>
          component.cell.transform.tags?.includes('structure-component-static-polymer'))?.representations ?? [];
        viewer.plugin.managers.structure.hierarchy.toggleVisibility(all.filter(item => item !== selected), 'hide');
        viewer.plugin.managers.structure.hierarchy.toggleVisibility([selected], 'show');
        setDisplayError(null);
      } catch (cause) {
        if (generation === representationGeneration.current)
          setDisplayError(cause instanceof Error ? cause.message : 'The representation could not be changed.');
      }
    };
    void apply();
    return () => { representationGeneration.current += 1; };
  }, [viewer, representation]);

  useEffect(() => {
    if (!viewer) return;
    const components = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [];
    const matching = (tag: string) => components.filter(component => component.cell.transform.tags?.includes(tag));
    const dmp = matching('structure-component-explicit-dmp-lipids')[0];
    const ligands = matching('structure-component-static-ligand').filter(component =>
      !dmp || component.cell.obj?.data.elementCount !== dmp.cell.obj?.data.elementCount);
    const setVisible = (items: typeof components, visible: boolean) => {
      if (items.length) viewer.plugin.managers.structure.hierarchy.toggleVisibility(items, visible ? 'show' : 'hide');
    };
    setVisible(matching('structure-component-static-polymer'), showProtein &&
      !(localInspection && contextCartoon === 'hidden'));
    setVisible([...matching('structure-component-static-lipid'), ...matching('structure-component-explicit-dmp-lipids')], showLipids);
    setVisible(ligands, showLigands);
    setVisible(matching('structure-component-static-water'), waterDisplay === 'all');
    setVisible(matching(sampledWaterTag), waterDisplay === 'sample');
    setVisible(matching('structure-component-static-ion'), showIons);
  }, [viewer, showProtein, showLipids, showLigands, waterDisplay, showIons, localInspection, contextCartoon]);

  useEffect(() => {
    if (!viewer) return;
    const canvas = viewer.plugin.canvas3d;
    if (!canvas) return;
    const primary = Binding([Binding.Trigger(ButtonsType.Flag.Primary, ModifiersKeys.create())]);
    canvas.setAttribs({ trackball: { bindings: {
      ...DefaultTrackballBindings,
      dragRotate: toolMode === 'pan' || toolMode === 'zoom' ? Binding.Empty : DefaultTrackballBindings.dragRotate,
      dragPan: toolMode === 'pan' ? primary : DefaultTrackballBindings.dragPan,
      dragZoom: toolMode === 'zoom' ? primary : DefaultTrackballBindings.dragZoom,
    } } });
    viewer.plugin.selectionMode = true;
  }, [viewer, placement, executionReview, toolMode]);

  useEffect(() => {
    if (!viewer || toolMode !== 'select' || !inspection?.structureUrl) return;
    const subscription = viewer.plugin.behaviors.interaction.click.subscribe(({ current }) => {
      const loci = current.loci;
      // Mol* publishes its initial empty BehaviorSubject value on subscription.
      // That value is not a researcher selection and must not create a card.
      if (loci.kind === 'empty-loci' ||
          StructureElement.Loci.is(loci) && StructureElement.Loci.isEmpty(loci)) return;
      if (!StructureElement.Loci.is(loci)) {
        onPickUnavailable('This displayed item has no mapped molecular identity.');
        return;
      }
      const residueLoci = StructureElement.Loci.firstResidue(loci);
      if (StructureElement.Loci.isEmpty(residueLoci) || !StructureElement.Loci.isSubset(residueLoci, loci)) {
        onPickUnavailable('Choose one visible residue or molecule to identify it.');
        return;
      }
      const location = StructureElement.Loci.getFirstLocation(loci);
      if (!location || !Unit.isAtomic(location.unit)) {
        onPickUnavailable('A mapped molecular identity is unavailable for this part of the view.');
        return;
      }
      const atomSiteIndex = location.unit.model.atomicHierarchy.atomSourceIndex.value(location.element);
      if (!Number.isInteger(atomSiteIndex) || atomSiteIndex < 0) {
        onPickUnavailable('The selected atom has no verified coordinate-row identity.');
        return;
      }
      viewer.plugin.managers.structure.focus.clear();
      viewer.plugin.managers.interactivity.lociSelects.selectOnly({ loci: residueLoci }, false);
      setPickedLoci(residueLoci);
      const residueNumber = StructureProperties.residue.auth_seq_id(location);
      const insertionCode = StructureProperties.residue.pdbx_PDB_ins_code(location);
      onPickAtom({ atomSiteIndex, residueName: StructureProperties.residue.auth_comp_id(location) || null,
        chain: StructureProperties.chain.auth_asym_id(location) || StructureProperties.chain.label_asym_id(location),
        residueNumber: Number.isFinite(residueNumber) && residueNumber > 0 ? residueNumber : null,
        insertionCode: insertionCode === '.' || insertionCode === '?' ? '' : insertionCode,
        modelNumber: StructureProperties.unit.model_num(location),
        entityKind: StructureProperties.entity.type(location),
        atomName: StructureProperties.atom.auth_atom_id(location),
        element: StructureProperties.atom.type_symbol(location),
        atomLevel: StructureElement.Loci.size(loci) === 1 &&
          (representation === 'sticks' || representation === 'spacefill') });
    });
    return () => subscription.unsubscribe();
  }, [viewer, toolMode, inspection?.structureUrl, representation, onPickAtom, onPickUnavailable]);

  useEffect(() => {
    if (!viewer) return;
    const focus = inspection?.focus;
    if (!focus || !focus.authAsymId || !Number.isInteger(focus.authSeqId)) {
      setAnnotationLoci(null);
      return;
    }
    const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
    if (!structure) return;
    const loci = StructureElement.Loci.fromSchema(structure, {
      auth_asym_id: focus.authAsymId,
      auth_seq_id: focus.authSeqId,
      ...(focus.insertionCode ? { pdbx_PDB_ins_code: focus.insertionCode } : {}),
      ...(focus.authAtomId !== null ? { atom_id: focus.authAtomId } : {}),
    });
    if (StructureElement.Loci.isEmpty(loci)) { setAnnotationLoci(null); return; }
    setPickedLoci(null);
    setAnnotationLoci(StructureElement.Loci.firstResidue(loci));
    viewer.plugin.managers.interactivity.lociSelects.selectOnly({ loci }, false);
  }, [viewer, inspection?.focusId, inspection?.focus?.authAsymId, inspection?.focus?.authSeqId, inspection?.focus?.insertionCode, inspection?.focus?.authAtomId]);

  useEffect(() => {
    if (!viewer) return;
    const target = pickedLoci ?? annotationLoci;
    if (localInspection && target) {
      viewer.plugin.managers.structure.focus.setFromLoci(target);
      viewer.plugin.managers.camera.focusLoci(target, { extraRadius: 5, durationMs: 0 });
    } else {
      viewer.plugin.managers.structure.focus.clear();
    }
  }, [viewer, localInspection, pickedLoci, annotationLoci]);

  useEffect(() => {
    if (!viewer) return;
    let current = true;
    void viewer.plugin.state.updateBehavior(StructureFocusRepresentation, params => {
      params.expandRadius = 5;
      params.components = localInspection && showContacts
        ? ['target', 'surroundings', 'interactions'] : ['target', 'surroundings'];
    }).then(() => {
      // Mol* creates a newly enabled interaction representation when the
      // focused loci are next applied; updating the option alone cannot add it.
      const target = pickedLoci ?? annotationLoci;
      if (current && localInspection && showContacts && target)
        viewer.plugin.managers.structure.focus.setFromLoci(target);
    }).catch(cause => { if (current) setDisplayError(cause instanceof Error
      ? cause.message : 'Contact display is unavailable.'); });
    return () => { current = false; };
  }, [viewer, localInspection, showContacts, pickedLoci, annotationLoci]);

  useEffect(() => {
    if (!viewer) return;
    const polymer = viewer.plugin.managers.structure.hierarchy.current.structures[0]?.components.filter(component =>
      component.cell.transform.tags?.includes('structure-component-static-polymer')) ?? [];
    if (!polymer.length) return;
    let current = true;
    const update = async () => {
      await clearStructureTransparency(viewer.plugin, polymer);
      if (current && localInspection && showProtein && contextCartoon === 'faded')
        await setStructureTransparency(viewer.plugin, polymer, 0.78,
          async structure => StructureElement.Loci.all(structure));
    };
    void update().catch(cause => { if (current) setDisplayError(cause instanceof Error
      ? cause.message : 'Context transparency is unavailable.'); });
    return () => { current = false; };
  }, [viewer, localInspection, showProtein, contextCartoon, representation]);

  useEffect(() => {
    if (!viewer || clearSelectionSerial === 0) return;
    setPickedLoci(null);
    setAnnotationLoci(null);
    setShowContacts(false);
    viewer.plugin.managers.interactivity.lociSelects.deselectAll();
    viewer.plugin.managers.structure.focus.clear();
  }, [viewer, clearSelectionSerial]);

  useEffect(() => {
    if (!viewer || wholeSystemSerial === 0) return;
    let frame = window.requestAnimationFrame(() => {
      // The input pane may have changed the canvas width while the researcher
      // was inspecting locally. Resize before either overview camera reads it.
      viewer.plugin.canvas3d?.handleResize();
      if (placement) focusPlacement(viewer, placement);
      else if (fullSystemReview) focusExplicitSystem(viewer);
      else viewer.plugin.canvas3d?.requestCameraReset({ durationMs: 0 });
    });
    return () => window.cancelAnimationFrame(frame);
  }, [viewer, wholeSystemSerial, placement, fullSystemReview]);

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
    <div className="scene-inspection-toolbar" aria-label="Molecular inspection display controls">
      <label className="scene-representation-label">Protein view
        <select className="select-input" aria-label="Protein representation" value={representation}
          disabled={!viewer || !displayAvailable.protein} onChange={event => setRepresentation(event.target.value as MolecularRepresentation)}>
          <option value="cartoon">Cartoon</option>
          <option value="sticks">Sticks</option>
          <option value="surface">Surface</option>
          <option value="spacefill">Space-filling</option>
        </select>
      </label>
      <div className="scene-component-tool">
        <button className="button compact" type="button" aria-expanded={displayOpen} aria-controls="scene-component-options"
          onClick={() => setDisplayOpen(value => !value)}>Components <span aria-hidden="true">⌄</span></button>
        {displayOpen && <div id="scene-component-options" className="scene-component-options" aria-label="Visible molecular components">
          <strong>Visible components</strong>
          <label><input type="checkbox" checked={showProtein} disabled={!displayAvailable.protein}
            onChange={event => setShowProtein(event.target.checked)} />Protein chains{displayAvailable.chains.length > 0 && ` (${displayAvailable.chains.join(', ')})`}</label>
          <label><input type="checkbox" checked={showLipids} disabled={!displayAvailable.lipids}
            onChange={event => setShowLipids(event.target.checked)} />Lipids{!displayAvailable.lipids && ' · absent'}</label>
          <label><input type="checkbox" checked={showLigands} disabled={!displayAvailable.ligands}
            onChange={event => setShowLigands(event.target.checked)} />Ligands / retained partners{!displayAvailable.ligands && ' · absent'}</label>
          <fieldset className="execution-water-options" disabled={!displayAvailable.water}>
            <legend>Water molecules{!displayAvailable.water && ' · absent'}</legend>
            {displayAvailable.sampledWater && <label><input type="radio" name="inspection-water" checked={waterDisplay === 'sample'}
              onChange={() => setWaterDisplay('sample')} />Representative sample</label>}
            <label><input type="radio" name="inspection-water" checked={waterDisplay === 'all'}
              onChange={() => setWaterDisplay('all')} />All</label>
            <label><input type="radio" name="inspection-water" checked={waterDisplay === 'hidden'}
              onChange={() => setWaterDisplay('hidden')} />Hidden</label>
          </fieldset>
          <label><input type="checkbox" checked={showIons} disabled={!displayAvailable.ions}
            onChange={event => setShowIons(event.target.checked)} />Ions{!displayAvailable.ions && ' · absent'}</label>
          {placement && <label><input type="checkbox" checked={showBilayer}
            onChange={event => setShowBilayer(event.target.checked)} />Membrane preview · not constructed lipids</label>}
          {localInspection && <>
            <label>Surrounding protein
              <select value={contextCartoon} onChange={event => setContextCartoon(event.target.value as typeof contextCartoon)}>
                <option value="faded">Faded</option><option value="shown">Shown</option><option value="hidden">Hidden</option>
              </select>
            </label>
            <label><input type="checkbox" checked={showContacts} onChange={event => setShowContacts(event.target.checked)} />Geometric contact aids</label>
            <small>When shown, Mol* draws computed noncovalent candidates as dashed lines by interaction type. These are display aids, not validated bonds or scientific assessment.</small>
          </>}
          {constructed && <small>Full model: {constructed.atomCount.toLocaleString()} atoms; {constructed.waterCount.toLocaleString()} waters; {constructed.sodiumCount.toLocaleString()} Na⁺ and {constructed.chlorideCount.toLocaleString()} Cl⁻.</small>}
          {localInspection && <small>The 5 Å neighborhood uses the selected structure's nearby atoms, including waters or ions even when those categories are hidden in the whole-system view.</small>}
          <small>Visibility and rendering never change coordinates, retained chemistry or qualification.</small>
        </div>}
      </div>
      {hasSelection && <button className="button compact" type="button" onClick={onClearSelection}>Clear selection</button>}
      <button className="button compact" type="button" disabled={!viewer} onClick={onShowWholeSystem}>Show whole system</button>
      {localInspection && <span className="scene-local-state">Local inspection · surroundings 5 Å</span>}
    </div>
    <div className="viewer-mount" ref={mount} aria-label={`3D structure for ${structureLabel}`} />
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
    </nav>}
    {displayError && <div className="scene-display-error" role="alert">Display change unavailable: {displayError}</div>}
    {fullSystemReview && viewer &&
      <div className="execution-display-disclosure" role="status">
        Protein {displayAvailable.protein ? showProtein ? 'shown' : 'hidden' : 'not rendered'} · lipids {displayAvailable.lipids ? showLipids ? 'shown' : 'hidden' : 'not rendered'} · water {constructed?.waterCount === 0 ? 'none' :
          !displayAvailable.water ? 'not rendered' : waterDisplay === 'sample' && displayAvailable.sampledWater
            ? 'sampled (residue IDs divisible by 32)' : waterDisplay === 'all' ? 'all shown' : 'hidden'} · ions {
              constructed && constructed.sodiumCount + constructed.chlorideCount === 0 ? 'none' :
                !displayAvailable.ions ? 'not rendered' : showIons ? 'shown' : 'hidden'}.
        {' '}The identified all-atom coordinates and counts are unchanged.
      </div>}
    {!fullSystemReview && !placement && viewer && <div className="scene-view-disclosure" role="status">
      Protein {displayAvailable.protein ? showProtein ? 'shown' : 'hidden' : 'absent'} · ligands {displayAvailable.ligands ? showLigands ? 'shown' : 'hidden' : 'absent'} · water {displayAvailable.water ? waterDisplay === 'hidden' ? 'hidden in overview' : 'shown' : 'absent'} · ions {displayAvailable.ions ? showIons ? 'shown' : 'hidden' : 'absent'}.
      {' '}{localInspection ? 'The local neighborhood may show components hidden in the overview.' : 'Cartoon and surface use chain colours; atomic views use element colours.'} Display only; coordinates are unchanged.
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

function rcsbEntryId(subjectId: string): string | null {
  return /^rcsb:([a-z0-9]{4})$/i.exec(subjectId)?.[1].toUpperCase() ?? null;
}

export function viewerSubjectHeading(inspection: InspectionAccount | null,
  stage: WorkspaceState['stages'][number] | undefined, sourceLabel?: string | null): string {
  if (!inspection) return 'Molecular viewer';
  const identifiedSource = sourceLabel?.trim() || rcsbEntryId(inspection.subjectId);
  switch (inspection.representationKind) {
    case 'structuralSource': {
      return identifiedSource ? `Source structure · ${identifiedSource}` : 'Source structure';
    }
    case 'intendedProtein':
    case 'selected-protein-before-repair': return `Protein before preparation${identifiedSource ? ` · ${identifiedSource}` : ''}`;
    case 'unqualifiedProteinCandidate': return 'Unqualified protein candidate';
    case 'preparedProtein': return `Prepared protein${identifiedSource ? ` · ${identifiedSource}` : ''}`;
    case 'membraneModel': return 'Intended membrane model';
    case 'oriented-protein-with-proposed-membrane-bounds': return 'Protein placement proposal';
    case 'constructedSystem': return 'Protein–membrane system';
    case 'completedStage': return stage?.kind === 'Minimization' ? 'Minimized system'
      : stage?.kind === 'Equilibration' ? 'Equilibrated system' : 'Completed molecular stage';
    default: return 'Molecular viewer';
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

function EditableBilayerPreview({ upper, lower, standing }: {
  upper: LipidFraction[]; lower: LipidFraction[]; standing: 'draft' | 'proposal' | 'chosen';
}) {
  const species = [...new Set([...upper, ...lower].map(item => item.speciesId))];
  const color = (id: string) => `hsl(${(206 + species.indexOf(id) * 137.5) % 360} 48% 60%)`;
  const leaflet = (side: 'upper' | 'lower', fractions: LipidFraction[]) => {
    const total = fractions.reduce((sum, item) => sum + item.fraction, 0);
    return <div className={`bilayer-leaflet-band ${side} ${fractions.length === 0 ? 'unspecified' : ''}`}>
      {fractions.map(item => <div className="bilayer-species-segment" key={item.speciesId}
        style={{ width: `${item.fraction * 100}%`, backgroundColor: color(item.speciesId) }} title={`${item.speciesId} ${(item.fraction * 100).toFixed(1)}%`} />)}
      {total < 0.9999 && <span className="bilayer-unassigned" style={{ width: `${Math.max(0, 100 - total * 100)}%` }}>
        {fractions.length ? `${Math.max(0, 100 - total * 100).toFixed(1)}% unassigned` : 'Not specified'}
      </span>}
    </div>;
  };
  return <figure className="bilayer-preview membrane-area-preview">
    <div className="bilayer-preview-title"><span>Physical leaflets</span><strong>{standing === 'draft' ? 'Composition draft' : standing === 'proposal' ? 'Membrane proposal' : 'Chosen membrane'}</strong></div>
    <div className="bilayer-preview-graphic" role="img" aria-label={`Membrane ${standing}. Upper leaflet ${upper.length ? upper.map(item => `${item.speciesId} ${(item.fraction * 100).toFixed(1)} percent`).join(', ') : 'not specified'}. Lower leaflet ${lower.length ? lower.map(item => `${item.speciesId} ${(item.fraction * 100).toFixed(1)} percent`).join(', ') : 'not specified'}. No lipid coordinates are shown.`}>
      <div className="bilayer-side-label">Upper aqueous side</div>
      <div className="bilayer-physical-label">Upper physical leaflet</div>{leaflet('upper', upper)}
      <div className="bilayer-core">Bilayer core · no molecule positions shown</div>
      {leaflet('lower', lower)}<div className="bilayer-physical-label">Lower physical leaflet</div>
      <div className="bilayer-side-label">Lower aqueous side</div>
    </div>
    {species.length > 0 && <div className="bilayer-preview-legend" aria-label="Selected lipid species">{species.map(id =>
      <span className="bilayer-legend-item" key={id}><i className="bilayer-legend-swatch" style={{ backgroundColor: color(id) }} />{id}</span>)}</div>}
    <figcaption>Composition preview — membrane not yet built.</figcaption>
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
  const standing = placement.status === 'supported' ? 'Position ready to use' :
    placement.status === 'unsupported' ? 'Position failed a technical check' :
      placement.status === 'assessing' ? 'Checking position…' :
      placement.status === 'proposed' ? 'Position under assessment' : 'Position check incomplete';
  const candidateBasis = [...placement.evidence].reverse().find(item =>
    /candidate|adjustment/i.test(item.method)) ?? placement.evidence[0];
  const selectedSource = state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId);
  const sourceLabel = state.study?.selectedSourceLabel ?? selectedSource?.label ?? (state.study?.selectedSourceKind === 'upload'
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
        <dt>Type</dt><dd>{placement.transform ? 'Complete prepared construct' :
          placement.topologyKind === 'membrane-spanning' ? 'Membrane-spanning protein' : 'Surface-associated protein'}</dd>
        <dt>Source</dt><dd>{sourceLabel}</dd>
        <dt>Model</dt><dd>{state.study?.modelIndex !== null && state.study?.modelIndex !== undefined
          ? `Source model ${state.sourceModels.find(model => model.index === state.study?.modelIndex)?.sourceModelId ?? state.study.modelIndex + 1} · chain ${chainSummary}` : `Chain ${chainSummary}`}</dd>
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
        <dt>Starting position</dt><dd>{placement.transform ? `${placement.transform.startingPosition} in membrane frame` :
          placement.topologyKind === 'membrane-spanning' ? 'Membrane-spanning method estimate' : 'Surface-associated method estimate'}</dd>
        <dt>Physical region</dt><dd>{placement.physicalSide === 'both' ? 'Membrane region' : `${placement.physicalSide} physical side`}</dd>
        {candidateBasis && <><dt>Candidate basis</dt><dd>{candidateBasis.method}</dd></>}
        {placement.transform && <><dt>Applied translation</dt><dd className="tabular">({placement.transform.appliedTranslationXAngstrom.toFixed(2)}, {placement.transform.appliedTranslationYAngstrom.toFixed(2)}, {placement.transform.appliedTranslationZAngstrom.toFixed(2)}) Å</dd>
          <dt>Rotations X/Y/Z</dt><dd className="tabular">{placement.transform.rotationXDegrees}° / {placement.transform.rotationYDegrees}° / {placement.transform.rotationZDegrees}°</dd></>}
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
      <h2>Spatial observations</h2>
      <div className="placement-orientation-list">
        {(['upper side', 'membrane core', 'lower side'] as const).map(label => {
          const anchor = inspection.annotations.find(item => item.label === label);
          return <button type="button" key={label} disabled={!anchor?.geometryFocus || !inspection.structureUrl || !canSelectFocus}
            className={anchor && inspection.focusId === anchor.subjectPartId ? 'selected' : ''}
            onClick={() => anchor && onSelectFocus(anchor.id)}>
            <span>{label === 'membrane core' ? 'Membrane core' : `${label === 'upper side' ? 'Upper' : 'Lower'} physical side`}</span>
            <strong>{anchor?.geometryFocus && inspection.structureUrl ? 'Located in displayed coordinates' : 'No verified spatial anchor'}</strong>
            {anchor?.evidenceId && <small>See attributed evidence below</small>}
          </button>;
        })}
      </div>
      {placement.contactingRegions.length > 0
        ? <details className="placement-contact-details"><summary>{placement.contactingRegions.length} exact contacting region{placement.contactingRegions.length === 1 ? '' : 's'} · inspect list</summary>
          <ul>{placement.contactingRegions.map((region, index) => <li key={`${region}:${index}`}>{region}</li>)}</ul>
        </details>
        : <p>No located protein region has yet been measured within the intended membrane core.</p>}
      {!placement.transform && placement.sidedness && <p className="help-text">Method context: {placement.sidedness}</p>}
    </section>
    <section className={`account-card placement-support-account ${statusClass(placement.status)}`} aria-label="Placement standing and uncertainty">
      <h2>Technical position checks</h2>
      <p className="placement-standing">{standing}</p>
      <p>{placement.reason}</p>
      {(placement.policyId || placement.witnessId) && <details className="source-provenance"><summary>Method details</summary><dl className="detail-grid">
        {placement.policyId && <><dt>Policy</dt><dd>{placement.policyId}{placement.policyVersion && ` · version ${placement.policyVersion}`}</dd></>}
        {placement.witnessId && <><dt>Independent witness</dt><dd>{placement.witnessId}</dd></>}
      </dl></details>}
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
  proteinReview,
  proteinOutcome,
  proteinDraft,
  placementOutcome,
  proposalDecision,
  exportFault = null,
  reviewAttempt = false,
  connectionMessage = null,
  onRefreshAccount,
  structureReload = 0,
  structureStatus = null,
  onStructureLoad,
  requestedSource = null,
  sourcePreviewLabel = null,
  onChainColors,
  chainFocus = null,
  activeArea,
  membraneDraft,
  membraneDraftMatchesProposal = false,
  requestedAreaView = null,
}: {
  state: WorkspaceState;
  onSelectFocus: (annotationId: string | null) => void;
  onInspectSubject: (subjectId: string) => void;
  proteinReview?: ReactNode;
  proteinOutcome?: ReactNode;
  proteinDraft?: ReactNode;
  placementOutcome?: ReactNode;
  proposalDecision?: ReactNode;
  exportFault?: string | null;
  reviewAttempt?: boolean;
  connectionMessage?: string | null;
  onRefreshAccount?: () => void;
  structureReload?: number;
  structureStatus?: StructureLoadStatus | null;
  onStructureLoad?: (status: StructureLoadStatus) => void;
  requestedSource?: { label: string; phase: 'retrieving' | 'rendering' | 'previewing' | 'previewFailed'; reason?: string } | null;
  sourcePreviewLabel?: string | null;
  onChainColors?: (subjectId: string, structureUrl: string, colors: Record<string, string>) => void;
  chainFocus?: { chainId: string; serial: number } | null;
  activeArea?: 'protein' | 'membrane' | 'placement' | 'preparation' | 'results';
  membraneDraft?: { upper: LipidFraction[]; lower: LipidFraction[] };
  membraneDraftMatchesProposal?: boolean;
  requestedAreaView?: { area: string; subjectId: string; failure: string | null } | null;
}) {
  const inspection = state.inspection;
  // An option switch can retain the same verified coordinate artifact and Mol*
  // canvas. Rendering standing follows those bytes; decision evidence still
  // follows the exact inspection subject and study revision below.
  const currentStructureStatus = structureStatus &&
    structureStatus.structureUrl === inspection?.structureUrl ? structureStatus : null;
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
    ?? (membrane ? membrane.status === 'notEstablished' ? 'Not established' : membrane.status === 'assessing' ? 'Assessing' :
      membrane.status === 'unavailable' ? 'Unavailable' : membrane.status === 'assessed' ? 'Assessed' : 'Proposed' : null);
  const selectedAction = state.actions.find(item => item.kind === 'setInspectionFocus' && item.subjectId === null);
  const canSelectFocus = selectedAction?.enabled === true;
  const affectedAnnotation = preparationChange && inspection?.annotations.find(item => item.geometryFocus);
  const [pickedAtom, setPickedAtom] = useState<AtomPickAccount | null>(null);
  const [pickedLocation, setPickedLocation] = useState<PickedLocation | null>(null);
  const [pickMessage, setPickMessage] = useState<string | null>(null);
  const [localInspection, setLocalInspection] = useState(false);
  const [clearSelectionSerial, setClearSelectionSerial] = useState(0);
  const [wholeSystemSerial, setWholeSystemSerial] = useState(0);
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
  const isSourceInspection = inspection?.representationKind === 'structuralSource' &&
    subject === state.study?.selectedSourceId;
  const subjectHeading = viewerSubjectHeading(inspection, stage, state.study?.selectedSourceLabel);
  const sourceCandidate = state.sourceCandidates.find(candidate => candidate.id === inspection?.subjectId);
  const sourceRoute = (isSourceInspection ? state.study?.selectedSourceKind : sourceCandidate?.kind)
    ?? (subject?.startsWith('rcsb:') ? 'rcsb' : subject?.startsWith('alphafold:') ? 'alphafold'
      : subject?.startsWith('upload:') ? 'upload' : null);
  const sourceOrigin = sourceRoute === 'alphafold' ? 'AlphaFold DB prediction'
    : sourceRoute === 'upload' ? `Researcher upload · ${isSourceInspection && state.study?.uploadProvenance &&
      state.study.uploadProvenance !== 'unknown'
      ? `researcher-declared ${state.study.uploadProvenance} origin` : 'origin not established'}`
      : sourceRoute === 'rcsb' ? sourceCandidate?.provenance.toLowerCase().includes('experimental entry search')
        ? 'RCSB PDB · experimental entry' : 'RCSB PDB' : 'Source origin unavailable';
  const revisionContext = inspection ? `Study revision ${inspection.studyRevisionNumber} · ${
    inspection.studyRevisionId === state.study?.id ? 'current study' : 'historical revision'}` : '';
  const subjectDetail = !inspection ? 'No subject selected'
    : inspection.representationKind === 'structuralSource'
      ? sourceOrigin
      : stage && inspection.representationKind === 'completedStage'
        ? 'Completed molecular result'
        : inspection.representationKind === 'constructedSystem'
          ? 'Checked constructed candidate; minimization not complete'
          : inspection.representationKind === 'oriented-protein-with-proposed-membrane-bounds'
            ? 'Positioned protein and proposed membrane bounds'
            : inspection.representationKind === 'membraneModel'
              ? 'Intended fractions; no lipid coordinates' : inspection.representationKind === 'intendedProtein' ||
                inspection.representationKind === 'selected-protein-before-repair'
                  ? 'Selected protein coordinates; reviewed changes are not yet applied' :
                    inspection.representationKind === 'preparedProtein'
                      ? 'Assessed prepared coordinates' : inspection.subjectId;
  const viewerSubtitle = requestedSource
    ? `${requestedSource.label} · ${requestedSource.phase === 'retrieving' ? 'retrieving source' :
      requestedSource.phase === 'previewFailed' ? 'visualization unavailable' : 'loading coordinates'}`
    : inspection ? `${subjectDetail}${isSourceInspection && sourcePreviewLabel ? ` · ${sourcePreviewLabel}` : ''} · ${revisionContext}` : subjectDetail;

  useEffect(() => {
    pickSequence.current += 1;
    setPickedAtom(null);
    setPickedLocation(null);
    setPickMessage(null);
    setLocalInspection(false);
  }, [selectedStructureKey]);

  useEffect(() => {
    if (inspection?.focusId && inspection.focus) setLocalInspection(true);
  }, [inspection?.focusId, selectedStructureKey]);

  const onPickUnavailable = useCallback((reason: string) => {
    pickSequence.current += 1;
    setPickedAtom(null);
    setPickedLocation(null);
    setPickMessage(reason);
    setLocalInspection(false);
    setClearSelectionSerial(value => value + 1);
  }, []);

  const onPickAtom = useCallback((location: PickedLocation) => {
    if (!inspection?.structureUrl) {
      onPickUnavailable('The selected subject has no verified structure for atom inspection.');
      return;
    }
    // Mol* reads only the structure bytes bound to this exact inspection URL.
    // A source has no source-to-result correspondence yet: report its actual
    // coordinate identity without inventing an adopted preparation atom.
    if (inspection.representationKind === 'structuralSource') {
      pickSequence.current += 1;
      setPickedAtom(null);
      setPickedLocation(location);
      setPickMessage(null);
      setLocalInspection(false);
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
    const atomSiteIndex = location.atomSiteIndex;
    const subjectId = inspection.subjectId;
    const revisionId = inspection.studyRevisionId;
    const selectedKey = selectedStructureKey;
    setPickedAtom(null);
    setPickedLocation(location);
    setLocalInspection(false);
    setPickMessage('Resolving selection against the host correspondence…');
    void (async () => {
      try {
        const response = await fetch(`/api/inspection/atoms/${encodeURIComponent(subjectId)}/${encodeURIComponent(token)}/${atomSiteIndex}`, { cache: 'no-store' });
        if (requestIndex !== pickSequence.current || selectedStructureKeyRef.current !== selectedKey) return;
        if (!response.ok) {
          setPickMessage(response.status === 404
            ? 'Source-to-result correspondence is unavailable for this selected coordinate; its displayed coordinate identity remains inspectable.'
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
  }, [inspection?.subjectId, inspection?.structureUrl, inspection?.studyRevisionId,
    inspection?.representationKind, selectedStructureKey, onPickUnavailable]);

  function clearSelection(restoreOverview: boolean) {
    pickSequence.current += 1;
    setPickedAtom(null);
    setPickedLocation(null);
    setPickMessage(null);
    setLocalInspection(false);
    setClearSelectionSerial(value => value + 1);
    if (restoreOverview) setWholeSystemSerial(value => value + 1);
    if (inspection?.focusId) onSelectFocus(null);
  }

  function showWholeSystem() {
    clearSelection(true);
  }

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

  if (requestedAreaView && activeArea !== 'membrane') return <>
    <section className="scene-panel" aria-label={`${requestedAreaView.area} view`}>
      <div className="scene-heading"><div><h1>{requestedAreaView.area} view</h1><div className="scene-subtitle">Restoring the area’s selected subject</div></div></div>
      <div className="scene-surface"><div className="scene-empty" role={requestedAreaView.failure ? 'alert' : 'status'}>
        <strong>{requestedAreaView.failure ? 'View unavailable' : `Restoring ${requestedAreaView.area.toLowerCase()} view…`}</strong>
        <span>{requestedAreaView.failure ?? 'Loading its structure and corresponding details without changing the scientific selection.'}</span>
      </div></div>
      <div className="scene-caption">The previous area’s molecule is not shown under this subject.</div>
    </section>
    <aside className="evidence-panel" aria-label={`${requestedAreaView.area} details`}><div className="evidence-content">
      <p className="eyebrow">{requestedAreaView.area} task</p><h2 className="panel-heading">Current view</h2>
      <div className={`hint-box${requestedAreaView.failure ? ' warning' : ''}`}>{requestedAreaView.failure ?? 'The exact subject and its details are being restored.'}</div>
    </div></aside>
  </>;

  if ((activeArea === 'results' && state.stages.length === 0) ||
      (activeArea === 'preparation' && !state.attempt && state.stages.length === 0) ||
      (activeArea === 'placement' && !state.protein && !state.placement)) {
    const label = activeArea === 'results' ? 'No completed system stage yet' :
      activeArea === 'preparation' ? 'No system preparation attempt yet' : 'Prepare a protein to position it';
    return <>
      <section className="scene-panel" aria-label={`${activeArea} view`}><div className="scene-heading"><div><h1>{activeArea === 'results' ? 'System results' : activeArea === 'preparation' ? 'System preparation' : 'Protein placement'}</h1>
        <div className="scene-subtitle">{label}</div></div></div><div className="scene-surface"><div className="scene-empty"><strong>{label}</strong><span>Return here when the corresponding molecular subject is available.</span></div></div>
        <div className="scene-caption">No completed structure is implied by opening this work area.</div></section>
      <aside className="evidence-panel"><div className="evidence-content"><p className="eyebrow">{activeArea} task</p><h2 className="panel-heading">Current standing</h2><div className="hint-box">{label}</div></div></aside>
    </>;
  }

  if (activeArea === 'membrane') {
    const selected = state.membrane;
    const draft = membraneDraft ?? { upper: [], lower: [] };
    const showDraft = !selected || !membraneDraftMatchesProposal;
    const shown = showDraft ? draft : selected;
    const standing = showDraft ? 'draft' : selected?.status === 'proposed' ? 'proposal' : 'chosen';
    const heading = standing === 'draft' ? 'Membrane composition draft' : standing === 'proposal' ? 'Membrane proposal' : 'Chosen membrane';
    const outcome = selected?.status === 'assessed' ? 'Membrane ready for placement' :
      selected?.status === 'assessing' ? 'Checking membrane…' :
      selected?.status === 'notEstablished' ? 'Membrane support not established' :
      selected?.status === 'unavailable' ? 'Membrane check could not finish' :
      selected ? 'Proposal awaiting adoption' : 'No membrane chosen yet';
    return <>
      <section className="scene-panel" aria-label="Membrane composition view">
        <div className="scene-heading"><div><h1>{heading}</h1><div className="scene-subtitle">Upper and lower physical leaflets · {showDraft ? 'editable composition' : selected?.status === 'proposed' ? 'proposed intention' : 'adopted intention'}</div></div></div>
        <div className="scene-surface"><EditableBilayerPreview upper={shown.upper} lower={shown.lower} standing={standing} /></div>
        <div className="scene-caption">Composition preview — membrane not yet built.</div>
      </section>
      <aside className="evidence-panel" aria-label="Membrane details">
        <div className="evidence-content"><p className="eyebrow">Membrane task</p><h2 className="panel-heading">Membrane checks</h2>
          <section className="account-card" aria-label="Membrane outcome"><h3>{outcome}</h3>
            {selected?.status === 'assessing' && <p role="status"><span className="activity-spinner" aria-hidden="true" />Checking the chosen composition and its molecular representations…</p>}
            {selected?.reason && selected.status !== 'assessed' && <p>{selected.reason}</p>}
            {selected?.status === 'assessed' && <p>Parameters and a construction route are available for this chosen membrane. Positioning and packing remain separate steps.</p>}
            {!selected && <p>Choose species and target fractions for each leaflet in the input pane.</p>}
            {showDraft && selected && <p>The displayed draft differs from the {selected.status === 'proposed' ? 'proposal' : 'chosen membrane'} below.</p>}
          </section>
          {selected && <section className="account-card" aria-label="Recorded membrane composition"><h3>{selected.status === 'proposed' ? 'Proposed composition' : 'Chosen composition'}</h3>
            <dl className="detail-grid"><dt>Upper leaflet</dt><dd>{presentFractions(selected.upper).map(item => `${item.speciesId} ${targetPercent(item.fraction)}`).join(' · ')}</dd>
              <dt>Lower leaflet</dt><dd>{presentFractions(selected.lower).map(item => `${item.speciesId} ${targetPercent(item.fraction)}`).join(' · ')}</dd></dl>
            {selected.limitations.length > 0 && <details><summary>Technical details and limitations</summary>{selected.limitations.map(limit => <p key={limit}>{limit}</p>)}
              {selected.policyId && <p>{selected.policyId} · {selected.policyVersion}</p>}</details>}
          </section>}
          {state.study && <details><summary>Fixed study conditions</summary><p>Nominal pH {state.study.conditions.nominalPh}; target NaCl {state.study.conditions.targetNaClMolar} M; optional equilibration target {state.study.conditions.optionalTemperatureKelvin} K.</p></details>}
        </div>
      </aside>
    </>;
  }

  return <>
    <section className="scene-panel" aria-label={stage ? 'Selected completed molecular stage' : reviewAttempt ? 'Current attempt inspection' : membrane ? 'Intended membrane model' : placement ? 'Placed protein and intended bilayer' : 'Molecular structure'}>
      <div className="scene-heading">
        <div>
          <h1>{requestedSource ? 'Loading source structure' : subjectHeading}</h1>
          <div className="scene-subtitle tabular" title={requestedSource ? requestedSource.label : inspection?.studyRevisionId}>
            {viewerSubtitle}
          </div>
        </div>
        {!requestedSource && subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      <div className={`scene-surface${inspection?.structureUrl && requestedSource?.phase !== 'retrieving' ? ' has-inspection-toolbar' : ''}${placement && !reviewAttempt ? ' placement-scene-surface' : ''}${executionReview ? ' execution-scene-surface' : ''}`}>
        {requestedSource?.phase === 'retrieving' || requestedSource?.phase === 'previewing' || requestedSource?.phase === 'previewFailed' ? null : membrane && !reviewAttempt ? <IntendedBilayerPreview membrane={membrane} /> : <MolecularScene inspection={inspection} structureLabel={subjectHeading} placement={reviewAttempt ? null : placement}
          executionReview={executionReview} constructed={constructed} onPickAtom={onPickAtom} onPickUnavailable={onPickUnavailable}
          onStructureLoad={status => onStructureLoad?.(status)} reloadToken={structureReload}
          localInspection={localInspection} hasSelection={!!(pickedLocation || inspection?.focusId || localInspection)}
          onClearSelection={() => clearSelection(localInspection)} onShowWholeSystem={showWholeSystem}
          clearSelectionSerial={clearSelectionSerial} wholeSystemSerial={wholeSystemSerial}
          onChainColors={onChainColors} chainFocus={chainFocus} />}
        {requestedSource && <div className="scene-source-loading" role="status" aria-label={`Loading source ${requestedSource.label}`}>
          {requestedSource.phase !== 'previewFailed' && <span className="activity-spinner" aria-hidden="true" />}
          <strong>{requestedSource.phase === 'retrieving' ? 'Retrieving source…' :
            requestedSource.phase === 'previewFailed' ? 'Selected model view unavailable' : 'Loading structure…'}</strong>
          <span className="tabular">{requestedSource.label}</span>
          {requestedSource.reason && <span>{requestedSource.reason}</span>}
        </div>}
        {!requestedSource && executionReview && <div className="execution-scene-label">{inspection
          ? `${subjectHeading} · ${
            inspection.studyRevisionId === state.study?.id ? 'current study revision' : 'historical revision'}`
          : 'No molecular subject selected'}</div>}
      </div>
      <div className="scene-caption">
        {requestedSource ? <><strong>{requestedSource.phase === 'previewFailed' ? 'The selected model cannot be displayed.' : 'Loading the selected source view.'}</strong> {requestedSource.phase === 'previewFailed' ? 'Review the issue beside the model choice and retry.' : 'Its molecular view will appear when the exact coordinates finish loading.'}</>
          : executionReview ? <><strong>{stage ? 'Selected completed molecular stage.' : inspection?.representationKind === 'constructedSystem' ? 'Verified constructed starting system.' : 'Earlier inspected subject.'}</strong>{' '}
          {inspection?.omittedMolecules.length ? `Not rendered: ${inspection.omittedMolecules.join(', ')}.` : 'See the exact subject account for representation limits and measurements.'}
          {' '}The scene itself does not establish a scientific qualification.</>
          : placement ? <><strong>Positioned protein with intended bilayer bounds.</strong> The upper and lower translucent planes mark this proposal’s membrane frame. {placement.transform ? 'The full construct was moved and rotated as one rigid body.' : 'The optional orientation method supplies a starting estimate.'} No packed lipids or achieved system are shown.</>
          : membrane ? <><strong>Intended model, not achieved packing.</strong> The two physical leaflets show target fractions only; membrane-local support is stated in the evidence account.</> : <><strong>Spatial view is evidence, not assessment.</strong>{' '}
        {inspection?.omittedMolecules.length
          ? `Not shown: ${inspection.omittedMolecules.join(', ')}.`
          : 'Any available omissions and exact measurements appear with this subject’s account.'}</>}
      </div>
    </section>

    <aside className="evidence-panel" aria-label="Evidence and assessment">
      {requestedSource && inspection && <div className="notice source-evidence-pending" role="note">
        {requestedSource.phase === 'retrieving'
          ? `Loading ${requestedSource.label}. The evidence below still belongs to ${subjectHeading}.`
          : `Loading the coordinates for ${subjectHeading}. Its source account is available below.`}
      </div>}
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
      <p className="eyebrow">{proteinReview || proteinOutcome ? 'Protein task' : placementOutcome ? 'Placement task' : proteinDraft ? 'Protein selection' : 'Exact subject account'}</p>
      <div className="status-line">
        <h2 className="panel-heading">{proteinOutcome ? 'Protein result' : proteinReview ? 'Review site' : placementOutcome ? 'Placement outcome' : proteinDraft ? 'Protein to assess' : 'Evidence & standing'}</h2>
        {!proteinReview && !proteinOutcome && !placementOutcome && !proteinDraft && subjectStatus && <span className={`status-badge ${statusClass(subjectStatus)}`}>{subjectStatus}</span>}
      </div>
      {proteinReview}
      {proteinOutcome}
      {placementOutcome}
      {proteinDraft}
      {!proteinReview && !proteinOutcome && !placementOutcome && <>
      {!inspection && !reviewAttempt && <div className="hint-box">Select a scientific subject to connect its structure with the applicable values and findings.</div>}
      {reviewAttempt && state.attempt && <AttemptReviewAccount state={state} attempt={state.attempt} onInspectConstructed={onInspectSubject} />}
      {inspection && !stage && !isSourceInspection && <section className="account-card inspection-origin-account" aria-label="Selected subject and study revision">
        <h2>Viewing context</h2>
        <dl className="detail-grid">
          <dt>Subject</dt><dd>{subjectHeading}</dd>
          <dt>Origin</dt><dd>Study revision {inspection.studyRevisionNumber}</dd>
          <dt>Context</dt><dd>{inspection.studyRevisionId === state.study?.id ? 'Current study revision' :
            `Historical subject; current study is revision ${state.study?.number ?? 'unknown'} (${state.study?.id ?? 'unavailable'})`}</dd>
          {membrane && <><dt>Standing</dt><dd>{membrane.status === 'notEstablished' ? 'Support not established' :
            membrane.status === 'assessing' ? 'Assessing support for adopted intention' :
            membrane.status === 'unavailable' ? 'Assessment unavailable' :
            membrane.status === 'assessed' ? 'Assessed membrane model' : 'Proposed intention'}</dd></>}
          {inspection.assessment && <><dt>Stage assessment</dt><dd>{inspection.assessment.currentlyApplicable ? '' : 'Historical assessment · '}{inspection.assessment.qualification}: {inspection.assessment.reason}</dd></>}
        </dl>
        <details className="source-provenance"><summary>Exact subject and revision IDs</summary>
          <p className="tabular">Subject {inspection.subjectId}<br />Study revision {inspection.studyRevisionId}</p>
        </details>
      </section>}
      {isSourceInspection && state.study && <details className="source-details" open={!proteinDraft || undefined}><summary>Source details</summary><section className="account-card selected-source-account" aria-label="Selected source">
        <h2>Selected source</h2>
        <dl className="detail-grid">
          <dt>Structure</dt><dd>{state.study.selectedSourceLabel ?? 'Identified source'}</dd>
          <dt>Route</dt><dd>{state.study.selectedSourceKind === 'rcsb' ? 'RCSB PDB' : state.study.selectedSourceKind === 'alphafold' ? 'AlphaFold DB' : 'Researcher upload'}</dd>
          <dt>Revision</dt><dd>{inspection!.studyRevisionNumber} · current source</dd>
          <dt>Preparation model</dt><dd>{state.study.modelIndex === null ? state.sourceModels.length === 1
            ? 'Only coordinate model is a draft; protein selection not assessed' : 'Not chosen' :
              `Source model ${state.sourceModels.find(model => model.index === state.study?.modelIndex)?.sourceModelId ?? state.study.modelIndex + 1} · selected for assessment`}</dd>
        </dl>
        {state.study.selectedSourceKind === 'alphafold' && <p className="source-qualification-note">Prediction coordinates; not experimental validation.</p>}
        {state.study.selectedSourceKind === 'upload' && <p className="source-qualification-note">{state.study.uploadProvenance ? `Researcher-declared ${state.study.uploadProvenance} origin` : 'Origin not established'}; this label is not independently verified.</p>}
        {state.sourcePrediction && <p className="source-qualification-note">Prediction confidence: {state.sourcePrediction.localConfidence.filter(item => item.pLddt !== null).length}/{state.sourcePrediction.localConfidence.length} mapped residues with local values; PAE {state.sourcePrediction.paeStanding.toLowerCase()}. Prediction uncertainty does not establish preparation or placement support.</p>}
        <details className="source-provenance"><summary>Provenance and limitations</summary>
          <p className="tabular">Exact source {state.study.selectedSourceId}<br />Study revision {inspection!.studyRevisionNumber} · {inspection!.studyRevisionId}</p>
          {state.study.uploadProvenanceNote && <p>{state.study.uploadProvenanceNote}</p>}
          {state.sourceCandidates.find(candidate => candidate.id === state.study?.selectedSourceId)?.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
          {state.sourcePrediction?.paeReason && <p>PAE: {state.sourcePrediction.paeReason}</p>}
          {state.sourcePrediction?.limitations.map(limit => <p key={limit}>Limit: {limit}</p>)}
          <p>Viewing source coordinates does not select chains, retain partners or establish a prepared protein.</p>
        </details>
      </section></details>}
      {inspection && showSubjectEvidence && <>
        {!isSourceInspection && !preparationChange && !placement && !executionReview && !membrane && <section className="account-card" aria-label="Selected structure">
          <h2>Selected structure</h2>
          <dl className="detail-grid">
          <dt>Representation</dt><dd>{representationLabel(inspection.representationKind)}</dd>
          {inspection.focusId && <><dt>Selected part</dt><dd>{inspection.focusId}</dd></>}
          </dl>
        </section>}
        {membrane && !executionReview && <>
          <MembraneLeafletAccount side="Upper" fractions={membrane.upper} state={state} />
          <MembraneLeafletAccount side="Lower" fractions={membrane.lower} state={state} />
          <section className="account-card membrane-conditions-account" aria-label="Fixed study conditions">
            <h2>Fixed study conditions</h2>
            {membrane.scientificPurpose?.trim() && <p>Recorded purpose: {membrane.scientificPurpose}</p>}
            <dl className="detail-grid">
              <dt>Target salt and neutrality</dt><dd>{state.study ? `${state.study.conditions.targetNaClMolar} M NaCl with charge-neutralizing compatible monovalent counterions` : 'Not established'}</dd>
              <dt>Optional equilibration target</dt><dd>{state.study ? `${state.study.conditions.optionalTemperatureKelvin} K` : 'Not established'}</dd>
              {state.study && <><dt>Nominal pH</dt><dd>{state.study.conditions.nominalPh}</dd></>}
            </dl>
            <p className="help-text">These disclosed targets are not achieved conditions or controls for this intended model.</p>
          </section>
          <section className={`account-card membrane-support-account ${membrane.status === 'assessed' ? 'established' : membrane.status === 'notEstablished' ? 'unestablished' : ''}`} aria-label="Membrane-local support and uncertainty">
            <h2>Membrane-local support and uncertainty</h2>
            <p className="membrane-support-standing">{membrane.status === 'assessed' ? 'Assessed membrane model' :
              membrane.status === 'assessing' ? 'Assessing adopted membrane model…' :
              membrane.status === 'unavailable' ? 'Assessment unavailable' :
              membrane.status === 'notEstablished' ? 'Membrane model not established' : 'Proposal not yet assessed'}</p>
            {membrane.reason && <p className="membrane-support-reason">{membrane.reason}</p>}
            {membrane.status === 'assessing' && <p role="status"><span className="activity-spinner" aria-hidden="true" />Checking the exact leaflet composition and selected representations…</p>}
            {membrane.status === 'unavailable' && <p>Return to Membrane to retry the same adopted intention after the service is available.</p>}
            {membrane.policyId && <details className="source-provenance"><summary>Support method</summary>
              <p>{membrane.policyId}{membrane.policyVersion && ` · version ${membrane.policyVersion}`}</p>
            </details>}
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
        {inspection.structureUrl && <section className="account-card inspection-atom-account" aria-label="Inspection selection">
          <h2>Inspection selection</h2>
          {visiblePickedAtom ? <>
            <p className="inspection-selection-name">{visiblePickedAtom.atom.sourceResidue
              ? `${pickedLocation?.residueName ?? 'Residue'} ${visiblePickedAtom.atom.sourceResidue.residue}${visiblePickedAtom.atom.sourceResidue.insertionCode} · chain ${visiblePickedAtom.atom.sourceResidue.chain}`
              : `${pickedLocation?.residueName ?? visiblePickedAtom.atom.generatedSpeciesId ?? visiblePickedAtom.atom.moleculeRole} · ${visiblePickedAtom.atom.moleculeRole}`}</p>
            {pickedLocation?.atomLevel && <p className="inspection-selection-atom tabular">Atom {visiblePickedAtom.atom.resultAtomId} · {visiblePickedAtom.atom.element}</p>}
            <p className="help-text">{localInspection ? 'Local atomic surroundings · 5 Å display extent. Contact aids are off until shown in Components.'
              : 'Selected in the whole-system view. Use Inspect locally to reveal atomic surroundings.'}</p>
            <div className="button-row inspection-selection-actions">
              <button className="button primary compact" type="button" disabled={localInspection} onClick={() => setLocalInspection(true)}>Inspect locally</button>
              <button className="button compact" type="button" onClick={() => clearSelection(localInspection)}>Clear selection</button>
            </div>
            <details className="inspection-correspondence"><summary>Exact coordinate correspondence</summary>
              <dl className="detail-grid">
                <dt>Coordinate atom</dt><dd className="tabular">{visiblePickedAtom.atom.resultAtomId}</dd>
                <dt>Role</dt><dd>{visiblePickedAtom.atom.moleculeRole} · {visiblePickedAtom.atom.atomRole} · {visiblePickedAtom.atom.element}</dd>
                <dt>Origin</dt><dd>{visiblePickedAtom.atom.sourceAtomId ? `Source atom ${visiblePickedAtom.atom.sourceAtomId}` :
                  `Generated ${visiblePickedAtom.atom.generatedComponentRole ?? visiblePickedAtom.atom.role}`}</dd>
                {visiblePickedAtom.atom.sourceResidue && <><dt>Source residue</dt><dd>Model {visiblePickedAtom.atom.sourceResidue.model} · chain {visiblePickedAtom.atom.sourceResidue.chain} · residue {visiblePickedAtom.atom.sourceResidue.residue}{visiblePickedAtom.atom.sourceResidue.insertionCode} · copy {visiblePickedAtom.atom.sourceResidue.copyId}</dd></>}
                {visiblePickedAtom.atom.generatedSpeciesId && <><dt>Species</dt><dd>{visiblePickedAtom.atom.generatedSpeciesId}</dd></>}
                {visiblePickedAtom.atom.physicalSide && <><dt>Physical leaflet</dt><dd>{visiblePickedAtom.atom.physicalSide}</dd></>}
                {visiblePickedAtom.atom.approvedChangeId && <><dt>Approved change</dt><dd>{visiblePickedAtom.atom.approvedChangeId}</dd></>}
              </dl>
            </details>
          </> : pickedLocation ? <>
            <p className="inspection-selection-name">{pickedLocation.residueName ?? pickedLocation.entityKind}{pickedLocation.residueNumber !== null
              ? ` ${pickedLocation.residueNumber}${pickedLocation.insertionCode}` : ''} · chain {pickedLocation.chain}</p>
            <p className="inspection-selection-atom">Model {pickedLocation.modelNumber} · {pickedLocation.entityKind}{pickedLocation.atomLevel
              ? ` · atom ${pickedLocation.atomName} (${pickedLocation.element})` : ''}</p>
            <p className="help-text">{localInspection ? 'Local atomic surroundings · 5 Å display extent.'
              : 'Selected from the exact displayed coordinates. Use Inspect locally to reveal atomic surroundings.'}</p>
            {pickMessage && <p className="inspection-mapping-note">{pickMessage}</p>}
            <div className="button-row inspection-selection-actions">
              <button className="button primary compact" type="button" disabled={localInspection} onClick={() => setLocalInspection(true)}>Inspect locally</button>
              <button className="button compact" type="button" onClick={() => clearSelection(localInspection)}>Clear selection</button>
            </div>
          </> : pickMessage ? <>
            <p className="help-text">{pickMessage}</p>
            <button className="button compact" type="button" onClick={() => clearSelection(localInspection)}>Clear selection</button>
          </> : inspection.focusId && inspection.focus ? <>
            <p className="inspection-selection-name">{inspection.annotations.find(item => item.subjectPartId === inspection.focusId)?.label
              ?? `Chain ${inspection.focus.authAsymId} · residue ${inspection.focus.authSeqId}`}</p>
            <p className="help-text">{localInspection ? 'Local atomic surroundings · 5 Å display extent.' : 'Selected finding on this subject.'}</p>
            <div className="button-row inspection-selection-actions">
              <button className="button primary compact" type="button" disabled={localInspection} onClick={() => setLocalInspection(true)}>Inspect locally</button>
              <button className="button compact" type="button" onClick={() => clearSelection(localInspection)}>Clear selection</button>
            </div>
          </> : currentStructureStatus?.phase === 'failed' ? <>
            <p className="inspection-selection-name">Structure unavailable for selection</p>
            <p className="help-text">{isSourceInspection
              ? 'Use Retry visualization beside the selected source. Its identity and evidence remain available.'
              : 'The identified subject remains available in its evidence account.'}</p>
          </> : currentStructureStatus?.phase === 'loading' ? <>
            <p className="inspection-selection-name">Loading structure…</p>
            <p className="help-text">Selection becomes available when the coordinates display.</p>
          </> : <><p className="inspection-selection-name">No residue selected</p>
            <p className="help-text">Click a visible molecule to identify it.</p></>}
          <p className="help-text">Selection and display do not alter coordinates, retained chemistry or scientific standing.</p>
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
        {currentSourceContext && !isSourceInspection && !executionReview && <section className="account-card source-account" aria-label="Source and assembly">
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
        {inspection.focusId && <div className="button-row"><button className="button compact" type="button" disabled={!canSelectFocus} onClick={() => clearSelection(localInspection)}>Clear selected part</button></div>}
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
      </>}
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
