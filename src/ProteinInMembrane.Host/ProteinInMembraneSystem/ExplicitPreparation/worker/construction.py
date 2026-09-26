"""Owner-local native OpenMM explicit construction and completed-stage observations.

The worker returns exact provider mechanics and observations. The C# owner decides
whether a candidate or completed stage meets an identified product policy.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (
    WorkError, artifact, require_integer, require_mapping, require_number,
    require_text, sha256, verify_sha256, work_path,
)
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _force_field_files, _provider, _coordinate_file, _topology_data,
    _read_topology_data, _nonbonded_charge, _full_atom_sequence,
    _bond_indices, _system_bonds_match, _load_stage, _stage_measurements,
    _state_context, _final_state,
)


def _coordinate_extent(points: list[tuple[float, float, float]]) -> tuple[float, float, float, float, float, float]:
    if not points:
        raise WorkError("invalidStructure", "No atoms can define the requested geometry")
    return (min(p[0] for p in points), max(p[0] for p in points),
            min(p[1] for p in points), max(p[1] for p in points),
            min(p[2] for p in points), max(p[2] for p in points))


def _signed_volume_at(points: list[tuple[float, float, float]], indices: tuple[int, int, int, int]) -> float:
    a, b, c, d = (points[index] for index in indices)
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    w = tuple(d[i] - a[i] for i in range(3))
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) -
            u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0]))


def _verify_stage_readback(pdb: Any, system: Any, state: Any) -> None:
    """Bind the rounded mmCIF and bonded sidecar to this exact serialized State.

    OpenMM writes atom coordinates and cell lengths to four decimal places in
    mmCIF. The small fixed tolerances cover that
    serialization only; a newly rehashed but different artifact remains a
    mismatch. Compare vectors, rather than just lengths, so cell orientation
    cannot change unnoticed.
    """
    import numpy as np
    from openmm import unit

    coordinates = np.asarray(pdb.positions.value_in_unit(unit.angstrom), dtype=float)
    state_coordinates = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
                                   dtype=float)
    if (coordinates.shape != state_coordinates.shape or not len(coordinates) or
            not np.isfinite(coordinates).all() or
            not np.isfinite(state_coordinates).all() or
            np.max(np.linalg.norm(coordinates - state_coordinates, axis=1)) > 0.0002):
        raise WorkError("inputMismatch", "Stage mmCIF coordinates disagree with the identified State")

    state_cell = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom),
                            dtype=float)
    def cell_values(vectors: Any) -> Any:
        if vectors is None:
            return None
        return np.asarray([vector.value_in_unit(unit.angstrom) for vector in vectors], dtype=float)

    for name, vectors in (("mmCIF", pdb.coordinate_box_vectors),
                          ("bonded topology sidecar", pdb.topology.getPeriodicBoxVectors()),
                          ("System", system.getDefaultPeriodicBoxVectors())):
        actual = cell_values(vectors)
        if (actual is None or actual.shape != (3, 3) or not np.isfinite(actual).all() or
                state_cell.shape != (3, 3) or not np.isfinite(state_cell).all() or
                np.max(np.linalg.norm(actual - state_cell, axis=1)) > 0.001):
            raise WorkError("inputMismatch", f"Stage {name} cell disagrees with the identified State")


def observe_stage(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm import VerletIntegrator, unit
    from .local_state_observations import observe_local_state
    from ProteinInMembraneSystem.worker.geometry_observations import observe_stage_protein_geometry

    pdb, system, state = _load_stage(directory, payload, "stateXmlPath")
    _verify_stage_readback(pdb, system, state)
    sidecar = _read_topology_data(work_path(directory, payload.get("topologyJsonPath"), "topologyJsonPath"))
    atom_order_matched = _full_atom_sequence(pdb.topology) == _full_atom_sequence(sidecar)
    bonds_matched = (_bond_indices(pdb.topology) == _bond_indices(sidecar) and
                     _system_bonds_match(pdb.topology, system))
    if not atom_order_matched or not bonds_matched:
        raise WorkError("inputMismatch", "Stage atom order or bonded System identity disagrees")
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = _state_context(system, state, integrator)
    observed = _final_state(context)
    measurements = _stage_measurements(observed, system.getNumParticles())
    correspondence_path = work_path(directory, payload.get("correspondencePath"), "correspondencePath")
    verify_sha256(correspondence_path, require_text(payload.get("correspondenceSha256"),
                                                    "correspondenceSha256"), "correspondenceSha256")
    correspondence = require_mapping(json.loads(correspondence_path.read_text(encoding="utf-8")),
                                     "stage correspondence")
    pdb.topology.setPeriodicBoxVectors(observed.getPeriodicBoxVectors())
    local_state = observe_local_state(pdb.topology, observed.getPositions(), correspondence,
                                      payload.get("localObservationSpec"))
    protein_geometry = observe_stage_protein_geometry(
        pdb.topology, observed.getPositions(), correspondence,
        payload.get("stageProteinGeometrySpec"))
    del context, integrator
    progress("stageMeasurementsObserved", {"stageKind": payload.get("stageKind")})
    return {"artifacts": [],
            "observations": {"atomCount": system.getNumParticles(), "measurements": measurements,
                             "atomOrderMatched": atom_order_matched,
                             "bondsMatched": bonds_matched,
                             "contactWarnings": local_state["limitations"],
                             "structuralWarnings": local_state["limitations"] + protein_geometry["limitations"],
                             "numericalWarnings": [], "localState": local_state,
                             "proteinGeometry": protein_geometry},
            "provider": _provider("OpenMM stage observation")}


def _native_representation(directory: Path, raw: Any, species: str, category: str,
                           force_field_hashes: set[str]) -> tuple[dict[str, Any], Any]:
    representation = require_mapping(raw, f"{species} molecular representation")
    if (representation.get("speciesId") != species or
            representation.get("category") != category or
            not require_text(representation.get("chemistryId"), f"{species} chemistryId")):
        raise WorkError("unsupportedPolicy", f"Native construction needs the exact {species} {category} representation")
    template = work_path(directory, representation.get("templatePath"), f"{species} templatePath")
    template_sha = require_text(representation.get("templateSha256"), f"{species} templateSha256")
    verify_sha256(template, template_sha, f"{species} templateSha256")
    if template_sha.lower() not in force_field_hashes:
        raise WorkError("inputMismatch", f"{species} molecular parameter asset is absent from selected force fields")
    coordinates = work_path(directory, representation.get("coordinateTemplatePath"),
                            f"{species} coordinateTemplatePath")
    verify_sha256(coordinates, require_text(representation.get("coordinateTemplateSha256"),
                                            f"{species} coordinateTemplateSha256"),
                  f"{species} coordinateTemplateSha256")
    reference = _coordinate_file(coordinates)
    atoms = list(reference.topology.atoms())
    if (len(list(reference.topology.residues())) != 1 or
            len(atoms) != require_integer(representation.get("atomCount"), f"{species} atomCount", 1) or
            any(atom.element is None for atom in atoms)):
        raise WorkError("invalidTemplate", f"{species} coordinate reference lacks one exact molecular identity")
    require_number(representation.get("netChargeElementary"), f"{species} netChargeElementary")
    return representation, reference


def _native_atom_name(name: str, species: str) -> str:
    # OpenMM's PDB reader moves the final digit of long Lipid21 acyl names to
    # the front. The exact reference atom order and bonds are checked below.
    if species in {"DMPC", "POPC"}:
        match = re.fullmatch(r"([0-9])C(21|31)", name)
        if match:
            return f"C{match.group(2)}{match.group(1)}"
    if species in {"NA", "CL"}:
        return name.upper()
    return name


def _native_signed_volume(points: list[tuple[float, float, float]],
                          names: dict[str, int], raw: Any, species: str) -> float:
    check = require_mapping(raw, "molecular stereo check")
    if check.get("kind") != "tetrahedral":
        raise WorkError("unsupportedPolicy", f"Native {species} stereo policy needs a tetrahedral descriptor")
    requested = check.get("atomNames")
    if (not isinstance(requested, list) or len(requested) != 4 or
            len(set(requested)) != 4 or any(name not in names for name in requested)):
        raise WorkError("invalidTemplate", f"{species} stereo descriptor does not identify four distinct atoms")
    value = _signed_volume_at(points, tuple(names[name] for name in
                                            (requested[3], requested[0], requested[1], requested[2])))
    expected = check.get("expected")
    if expected not in {"negative", "positive"} or not math.isfinite(value):
        raise WorkError("invalidTemplate", f"{species} stereo descriptor or observed volume is invalid")
    if (expected == "negative" and value >= -1.0 or
            expected == "positive" and value <= 1.0):
        raise WorkError("providerMismatch", f"A {species} glycerol stereocenter differs from the declared chemistry")
    return value


def _native_alkene_cosine(points: list[tuple[float, float, float]],
                          names: dict[str, int], bonds: set[tuple[int, int]],
                          raw: Any) -> float:
    check = require_mapping(raw, "molecular alkene check")
    requested = check.get("atomNames")
    if (check.get("kind") != "alkene" or check.get("expected") != "cis" or
            not isinstance(requested, list) or len(requested) != 4 or
            len(set(requested)) != 4 or any(name not in names for name in requested)):
        raise WorkError("invalidTemplate", "POPC needs its identified cis-alkene descriptor")
    indices = [names[name] for name in requested]
    if any(tuple(sorted((first, second))) not in bonds
           for first, second in zip(indices, indices[1:])):
        raise WorkError("invalidTemplate", "POPC cis descriptor does not follow bonded atoms")
    first, second, third, fourth = (points[index] for index in indices)
    axis = tuple(third[i] - second[i] for i in range(3))
    axis_squared = sum(value * value for value in axis)
    if not math.isfinite(axis_squared) or axis_squared <= 1e-8:
        raise WorkError("providerMismatch", "A POPC alkene has unresolved central-bond geometry")
    left = tuple(first[i] - second[i] for i in range(3))
    right = tuple(fourth[i] - third[i] for i in range(3))
    dot_left = sum(left[i] * axis[i] for i in range(3)) / axis_squared
    dot_right = sum(right[i] * axis[i] for i in range(3)) / axis_squared
    left_projected = tuple(left[i] - dot_left * axis[i] for i in range(3))
    right_projected = tuple(right[i] - dot_right * axis[i] for i in range(3))
    numerator = sum(left_projected[i] * right_projected[i] for i in range(3))
    denominator = math.sqrt(sum(value * value for value in left_projected) *
                            sum(value * value for value in right_projected))
    cosine = numerator / denominator if denominator > 1e-8 else math.nan
    if not math.isfinite(cosine) or cosine <= 0.5:
        raise WorkError("providerMismatch", "A POPC cis alkene differs from the declared chemistry")
    return cosine


def _native_molecule_bonds(topology: Any, residue: Any) -> set[tuple[int, int]]:
    atoms = list(residue.atoms())
    local = {atom.index: index for index, atom in enumerate(atoms)}
    return {tuple(sorted((local[first.index], local[second.index])))
            for first, second in topology.bonds()
            if first.residue is residue and second.residue is residue}


def _native_molecule_matches(residue: Any, reference: Any, species: str, positions: Any,
                             stereo_checks: Any, actual_bonds: set[tuple[int, int]],
                             expected_bonds: set[tuple[int, int]]) -> list[Any]:
    from openmm import unit

    actual = list(residue.atoms())
    expected = list(reference.topology.atoms())
    if (len(actual) != len(expected) or
            any(_native_atom_name(a.name, species) != e.name or a.element != e.element
                for a, e in zip(actual, expected)) or
            actual_bonds != expected_bonds):
        raise WorkError("providerMismatch", f"Native {species} atom order, element, or molecular bonds differ")
    if species in {"DMPC", "POPC"}:
        expected_count = 1 if species == "DMPC" else 2
        if not isinstance(stereo_checks, list) or len(stereo_checks) != expected_count:
            raise WorkError("invalidTemplate", f"{species} lacks its qualified stereo descriptors")
        points = [tuple(float(value) for value in positions[atom.index].value_in_unit(unit.angstrom))
                  for atom in actual]
        names = {expected_atom.name: index for index, expected_atom in enumerate(expected)}
        _native_signed_volume(points, names, stereo_checks[0], species)
        if species == "POPC":
            _native_alkene_cosine(points, names, actual_bonds, stereo_checks[1])
    return actual


def _native_lipid_side(atoms: list[Any], points: list[tuple[float, float, float]],
                       center_angstrom: float, species: str) -> str:
    if species not in {"DMPC", "POPC"}:
        raise WorkError("unsupportedPolicy", "Native leaflet geometry has no identified lipid")
    phosphorus = [atom for atom in atoms if atom.name == "P"]
    if len(phosphorus) != 1:
        raise WorkError("providerMismatch", f"Native {species} lacks its one phosphate side witness")
    head_z = points[phosphorus[0].index][2]
    if head_z == center_angstrom:
        raise WorkError("providerMismatch", f"A {species} phosphate lies on the physical side divider")
    side = "upper" if head_z > center_angstrom else "lower"
    tail_names = {"C214", "C314"} if species == "DMPC" else {"C218", "C316"}
    tails = [atom for atom in atoms if _native_atom_name(atom.name, species) in tail_names]
    if len(tails) != 2 or (head_z - sum(points[atom.index][2] for atom in tails) / 2) * (
            1 if side == "upper" else -1) <= 0:
        raise WorkError("providerMismatch", f"A native {species} has the wrong leaflet head-to-tail orientation")
    return side


def _native_provider_identity(payload: dict[str, Any]) -> tuple[Path, str, str]:
    import importlib.metadata
    import openmm
    import openmm.app

    name = require_text(payload.get("providerName"), "providerName")
    version = require_text(payload.get("providerVersion"), "providerVersion")
    if (name != "OpenMM Modeller.addMembrane" or
            version != openmm.version.full_version or
            importlib.metadata.version("openmm") != "8.6.0"):
        raise WorkError("providerMismatch", "Installed OpenMM provider is not the identified native 8.6 build")
    species = payload.get("lipidTypeArgument")
    if species not in {"DMPC", "POPC"}:
        raise WorkError("unsupportedPolicy", "Native construction covers identified DMPC or POPC only")
    mode = payload.get("nativePatchMode", "installed")
    expected_path = (Path(openmm.app.__file__).resolve().parent / "data" / f"{species}.pdb").resolve()
    if mode == "installed":
        declared_path = Path(require_text(payload.get("nativePatchPath"), "nativePatchPath")).resolve()
        if (declared_path != expected_path or not expected_path.is_file() or
                payload.get("nativeSourcePatchPath") is not None or
                payload.get("nativeSourcePatchSha256") is not None or
                payload.get("removedNativeLipidResidueIds") not in (None, [])):
            raise WorkError("providerMismatch", f"Declared {species} patch is not the installed OpenMM resource")
        expected_sha = require_text(payload.get("nativePatchSha256"), "nativePatchSha256")
        verify_sha256(expected_path, expected_sha, "nativePatchSha256")
        return expected_path, expected_sha, version
    if (mode != "popc-62-109-deletion" or species != "POPC" or
            payload.get("removedNativeLipidResidueIds") != ["62", "109"]):
        raise WorkError("unsupportedPolicy", "No identified custom native patch derivation matches this request")
    declared_source = Path(require_text(payload.get("nativeSourcePatchPath"),
                                        "nativeSourcePatchPath")).resolve()
    if declared_source != expected_path or not expected_path.is_file():
        raise WorkError("providerMismatch", "Custom POPC source is not the installed OpenMM resource")
    source_sha = require_text(payload.get("nativeSourcePatchSha256"), "nativeSourcePatchSha256")
    verify_sha256(expected_path, source_sha, "nativeSourcePatchSha256")
    derived = Path(require_text(payload.get("nativePatchPath"), "nativePatchPath")).resolve()
    if derived == expected_path or not derived.is_file():
        raise WorkError("providerMismatch", "Custom POPC patch must be a separate identified resource")
    derived_sha = require_text(payload.get("nativePatchSha256"), "nativePatchSha256")
    verify_sha256(derived, derived_sha, "nativePatchSha256")
    return derived, derived_sha, version


def _native_verified_custom_popc_patch(source: Any, derived: Any, reference: Any,
                                       stereo_checks: Any,
                                       removed_ids: list[str]) -> None:
    """Prove the custom patch only deletes the named opposing lipids."""
    from collections import Counter, defaultdict
    from openmm import unit

    original_residues = list(source.topology.residues())
    derived_residues = list(derived.topology.residues())
    omitted = [residue for residue in original_residues
               if residue.name == "POP" and residue.id in removed_ids]
    survivors = [residue for residue in original_residues if residue not in omitted]
    if (len(omitted) != 2 or [residue.id for residue in omitted] != removed_ids or
            len(survivors) != len(derived_residues) or
            Counter(residue.name for residue in derived_residues) !=
                Counter({"POP": 126, "HOH": 5120})):
        raise WorkError("providerMismatch", "Custom POPC patch does not have the exact two-lipid deletion")
    source_cell = source.topology.getPeriodicBoxVectors()
    derived_cell = derived.topology.getPeriodicBoxVectors()
    if (source_cell is None or derived_cell is None or
            any(abs(source_cell[i][axis].value_in_unit(unit.angstrom) -
                    derived_cell[i][axis].value_in_unit(unit.angstrom)) > 1e-6
                for i in range(3) for axis in range(3))):
        raise WorkError("providerMismatch", "Custom POPC patch changed the source periodic cell")

    def bonds_by_residue(topology: Any) -> dict[int, set[tuple[int, int]]]:
        local = {atom.index: index for residue in topology.residues()
                 for index, atom in enumerate(residue.atoms())}
        bonds: dict[int, set[tuple[int, int]]] = defaultdict(set)
        for first, second in topology.bonds():
            if first.residue is not second.residue:
                raise WorkError("providerMismatch", "Custom POPC patch has a cross-residue bond")
            bonds[first.residue.index].add(tuple(sorted((local[first.index], local[second.index]))))
        return bonds

    source_bonds = bonds_by_residue(source.topology)
    derived_bonds = bonds_by_residue(derived.topology)
    for original, retained in zip(survivors, derived_residues):
        original_atoms = list(original.atoms())
        retained_atoms = list(retained.atoms())
        if ((original.name, original.id, original.chain.id, original.insertionCode) !=
                (retained.name, retained.id, retained.chain.id, retained.insertionCode) or
                len(original_atoms) != len(retained_atoms) or
                source_bonds[original.index] != derived_bonds[retained.index]):
            raise WorkError("providerMismatch", "Custom POPC patch changed a survivor identity or bond")
        for before, after in zip(original_atoms, retained_atoms):
            before_identity = before.element, before.name, before.formalCharge
            after_identity = after.element, after.name, after.formalCharge
            if (before_identity != after_identity or
                    max(abs((source.positions[before.index][axis] -
                             derived.positions[after.index][axis]).value_in_unit(unit.angstrom))
                        for axis in range(3)) > 1e-4):
                raise WorkError("providerMismatch", "Custom POPC patch changed survivor atom identity or coordinate")

    expected_bonds = _native_molecule_bonds(reference.topology, next(reference.topology.residues()))
    points = [tuple(float(value) for value in point.value_in_unit(unit.angstrom))
              for point in derived.positions]
    center = derived_cell[2][2].value_in_unit(unit.angstrom) / 2
    sides: Counter[str] = Counter()
    for residue in derived_residues:
        if residue.name != "POP":
            continue
        atoms = _native_molecule_matches(residue, reference, "POPC", derived.positions,
                                         stereo_checks, derived_bonds[residue.index], expected_bonds)
        sides[_native_lipid_side(atoms, points, center, "POPC")] += 1
    if sides != Counter({"upper": 63, "lower": 63}):
        raise WorkError("providerMismatch", "Custom POPC patch lacks exact opposing 63/63 leaflets")


def _native_system_settings(ff: Any, topology: Any, settings_raw: Any,
                            cell: list[float], water_atoms: set[int]) -> Any:
    from openmm import CMMotionRemover, NonbondedForce, unit
    from openmm.app import HBonds, PME

    settings = require_mapping(settings_raw, "systemSettings")
    if settings.get("nonbondedMethod") != "PME" or settings.get("constraints") != "HBonds":
        raise WorkError("unsupportedPolicy", "Native full-system route requires PME and HBonds")
    cutoff = require_number(settings.get("nonbondedCutoffNanometers"),
                            "nonbondedCutoffNanometers", 0.000001)
    if cutoff >= min(cell) / 20:
        raise WorkError("invalidGeometry", "PME cutoff is at least half the native cell's shortest dimension")
    rigid_water = settings.get("rigidWater")
    dispersion = settings.get("useDispersionCorrection")
    remove_cm = settings.get("removeCMMotion")
    if any(not isinstance(value, bool) for value in (rigid_water, dispersion, remove_cm)):
        raise WorkError("invalidRequest", "OpenMM system booleans must be explicitly declared")
    ewald = require_number(settings.get("ewaldErrorTolerance"), "ewaldErrorTolerance", 0.000000001)
    switch_raw = settings.get("switchDistanceNanometers")
    switch = None if switch_raw is None else require_number(switch_raw, "switchDistanceNanometers", 0.000001)
    if switch is not None and switch >= cutoff:
        raise WorkError("invalidRequest", "Nonbonded switching must begin before the cutoff")
    hydrogen_raw = settings.get("hydrogenMassDaltons")
    hydrogen_mass = None if hydrogen_raw is None else require_number(
        hydrogen_raw, "hydrogenMassDaltons", 0.000001) * unit.dalton
    system = ff.createSystem(topology, nonbondedMethod=PME,
                             nonbondedCutoff=cutoff * unit.nanometer,
                             constraints=HBonds, rigidWater=rigid_water,
                             ewaldErrorTolerance=ewald,
                             switchDistance=None if switch is None else switch * unit.nanometer,
                             removeCMMotion=remove_cm, hydrogenMass=hydrogen_mass)
    if system.getNumParticles() != topology.getNumAtoms() or not _system_bonds_match(topology, system):
        raise WorkError("providerMismatch", "Native topology lacks complete combined System parameterization")
    forces = [system.getForce(index) for index in range(system.getNumForces())]
    nonbonded = [force for force in forces if isinstance(force, NonbondedForce)]
    if len(nonbonded) != 1:
        raise WorkError("providerMismatch", "Native full System needs exactly one nonbonded force")
    nonbonded[0].setUseDispersionCorrection(dispersion)
    if (nonbonded[0].getNonbondedMethod() != NonbondedForce.PME or
            abs(nonbonded[0].getCutoffDistance().value_in_unit(unit.nanometer) - cutoff) > 1e-10 or
            abs(nonbonded[0].getEwaldErrorTolerance() - ewald) > 1e-12 or
            nonbonded[0].getUseSwitchingFunction() != (switch is not None) or
            (switch is not None and abs(nonbonded[0].getSwitchingDistance().value_in_unit(unit.nanometer) - switch) > 1e-10) or
            nonbonded[0].getUseDispersionCorrection() != dispersion or
            sum(isinstance(force, CMMotionRemover) for force in forces) != int(remove_cm)):
        raise WorkError("providerMismatch", "Native System settings differ from the selected policy")
    constrained = {tuple(sorted(system.getConstraintParameters(index)[:2]))
                   for index in range(system.getNumConstraints())}
    for first, second in topology.bonds():
        if ((first.element.symbol == "H" or second.element.symbol == "H") and
                tuple(sorted((first.index, second.index))) not in constrained):
            raise WorkError("providerMismatch", "HBonds setting left a bonded hydrogen unconstrained")
    if rigid_water:
        for residue in topology.residues():
            indices = [atom.index for atom in residue.atoms() if atom.index in water_atoms and atom.element.symbol == "H"]
            if len(indices) == 2 and tuple(sorted(indices)) not in constrained:
                raise WorkError("providerMismatch", "Rigid-water setting left a water hydrogen pair unconstrained")
    if hydrogen_mass is not None:
        target = hydrogen_mass.value_in_unit(unit.dalton)
        for first, second in topology.bonds():
            for atom in (first, second):
                if (atom.element.symbol == "H" and atom.index not in water_atoms and
                        abs(system.getParticleMass(atom.index).value_in_unit(unit.dalton) - target) > 1e-6):
                    raise WorkError("providerMismatch", "Hydrogen mass differs from selected repartitioning")
    return system


def construct_system(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    """Observe one exact native membrane build; the product judges its result."""
    from collections import Counter, defaultdict
    from openmm import Context, Platform, VerletIntegrator, XmlSerializer, unit
    from openmm.app import ForceField, Modeller, PDBFile, PDBxFile
    from .local_state_observations import observe_local_state

    require_text(payload.get("studyRevisionId"), "studyRevisionId")
    require_text(payload.get("attemptId"), "attemptId")
    oriented_path = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    verify_sha256(oriented_path, require_text(payload.get("orientedPdbSha256"),
                                              "orientedPdbSha256"), "orientedPdbSha256")
    prepared_path = work_path(directory, payload.get("preparedPdbPath"), "preparedPdbPath")
    prepared_sha = require_text(payload.get("preparedPdbSha256"), "preparedPdbSha256")
    verify_sha256(prepared_path, prepared_sha, "preparedPdbSha256")
    graph_path = work_path(directory, payload.get("preparedBondGraphPath"), "preparedBondGraphPath")
    verify_sha256(graph_path, require_text(payload.get("preparedBondGraphSha256"),
                                          "preparedBondGraphSha256"), "preparedBondGraphSha256")
    lipid_type = payload.get("lipidTypeArgument")
    patch_mode = payload.get("nativePatchMode", "installed")
    installed_patch, patch_sha, provider_version = _native_provider_identity(payload)
    if (lipid_type not in {"DMPC", "POPC"} or
            payload.get("positiveIonArgument") != "Na+" or
            payload.get("negativeIonArgument") != "Cl-"):
        raise WorkError("unsupportedPolicy", "This native route covers identified pure DMPC or POPC with NaCl only")
    center = require_number(payload.get("membraneCenterZNanometers"), "membraneCenterZNanometers")
    padding = require_number(payload.get("minimumPaddingNanometers"), "minimumPaddingNanometers", 0.000001)
    ionic_strength = require_number(payload.get("ionicStrengthMolar"), "ionicStrengthMolar", 0)
    maximum_atoms = require_integer(payload.get("maximumAtomCount"), "maximumAtomCount", 1)
    maximum_cell = require_number(payload.get("maximumCellDimensionAngstrom"),
                                  "maximumCellDimensionAngstrom", 0.000001)
    if (center != 0 or padding != 1 or ionic_strength != 0.15 or
            maximum_atoms > 2_000_000 or maximum_cell > 10_000):
        raise WorkError("unsupportedPolicy", "Native membrane invocation differs from the bounded method")

    force_field_raw = payload.get("forceFieldFiles")
    ff_files = _force_field_files(directory, force_field_raw)
    ff_hashes = {require_text(raw.get("sha256"), "force-field sha256").lower()
                 for raw in force_field_raw}
    representations = {}
    references = {}
    for key, species, category in (("lipid", lipid_type, "lipid"), ("water", "HOH", "water"),
                                   ("sodium", "NA", "ion"), ("chloride", "CL", "ion")):
        representation, reference = _native_representation(
            directory, payload.get(key), species, category, ff_hashes)
        representations[species], references[species] = representation, reference
    lipid = representations[lipid_type]
    heads = lipid.get("headAtomIndices")
    if heads != [20] or list(references[lipid_type].topology.atoms())[19].name != "P":
        raise WorkError("invalidTemplate", f"The native {lipid_type} route requires the declared phosphate head index")

    oriented = PDBFile(str(oriented_path))
    prepared = PDBFile(str(prepared_path))
    approved_graph = _read_topology_data(graph_path)
    oriented_atoms = list(oriented.topology.atoms())
    prepared_atoms = list(prepared.topology.atoms())
    if len(oriented_atoms) > maximum_atoms:
        raise WorkError("resourceRefused", "The prepared protein alone exceeds the native atom bound")
    if (not oriented_atoms or
            _full_atom_sequence(oriented.topology) != _full_atom_sequence(prepared.topology) or
            _full_atom_sequence(oriented.topology) != _full_atom_sequence(approved_graph) or
            _bond_indices(oriented.topology) != _bond_indices(approved_graph)):
        raise WorkError("inputMismatch", "Oriented protein differs from the approved prepared identity and bonds")
    mapping = require_mapping(payload.get("preparedCorrespondence"), "preparedCorrespondence")
    mapped_atoms = mapping.get("atoms")
    if (mapping.get("complete") is not True or mapping.get("resultId") != prepared_sha or
            not isinstance(mapped_atoms, list) or len(mapped_atoms) != len(oriented_atoms)):
        raise WorkError("inputMismatch", "Prepared protein correspondence is incomplete or addresses another source")
    mapped_atoms = sorted((require_mapping(item, "prepared atom mapping") for item in mapped_atoms),
                          key=lambda item: require_integer(item.get("resultAtomIndex"), "resultAtomIndex"))
    if [item["resultAtomIndex"] for item in mapped_atoms] != list(range(len(oriented_atoms))):
        raise WorkError("inputMismatch", "Prepared mapping does not identify every ordered protein atom")
    for index, (atom, mapped) in enumerate(zip(prepared_atoms, mapped_atoms)):
        residue = atom.residue
        atom_id = f"{index}:{residue.chain.id}:{residue.id}:{residue.insertionCode.strip()}:{atom.name}"
        if (mapped.get("moleculeRole") != "protein" or mapped.get("element") != atom.element.symbol or
                mapped.get("resultAtomId") != atom_id or mapped.get("role") not in {"source", "generated"} or
                (mapped.get("role") == "source") != (mapped.get("sourceAtomId") is not None)):
            raise WorkError("inputMismatch", f"Prepared protein atom {index} lacks exact source correspondence")
    oriented_points = [tuple(float(value) for value in position.value_in_unit(unit.angstrom))
                       for position in oriented.positions]
    if any(not all(math.isfinite(value) for value in point) for point in oriented_points):
        raise WorkError("invalidStructure", "Oriented protein has a nonfinite coordinate")

    ff = ForceField(*ff_files)
    # The approved protein graph, rather than the PDB's implicit CONECT and
    # residue heuristics, supplies the protein charge and bonded identity.
    approved_graph.setUnitCellDimensions(None)
    protein_charge = _nonbonded_charge(ff.createSystem(approved_graph))
    if abs(protein_charge - round(protein_charge)) > 1e-4:
        raise WorkError("unsupportedPolicy", "Native monovalent neutralization requires integral protein charge")
    patch = PDBFile(str(installed_patch))
    native_residue_name = {"DMPC": "DMP", "POPC": "POP"}[lipid_type]
    patch_lipid = next((residue for residue in patch.topology.residues()
                        if residue.name == native_residue_name), None)
    if patch_lipid is None:
        raise WorkError("providerMismatch", f"Pinned OpenMM patch contains no {native_residue_name} molecular reference")
    reference_bonds = {species: _native_molecule_bonds(reference.topology,
                        next(reference.topology.residues())) for species, reference in references.items()}
    patch_bonds = _native_molecule_bonds(patch.topology, patch_lipid)
    _native_molecule_matches(patch_lipid, references[lipid_type], lipid_type, patch.positions,
                             lipid.get("stereoChecks"), patch_bonds, reference_bonds[lipid_type])
    if patch_mode == "popc-62-109-deletion":
        source_patch = PDBFile(require_text(payload.get("nativeSourcePatchPath"),
                                            "nativeSourcePatchPath"))
        _native_verified_custom_popc_patch(source_patch, patch, references[lipid_type],
                                           lipid.get("stereoChecks"),
                                           payload["removedNativeLipidResidueIds"])
    # A source PDB may carry a placeholder 1 A CRYST1 record.  Modeller uses
    # the input cell's Z when it exists, so it must be cleared before this call.
    oriented.topology.setUnitCellDimensions(None)
    modeller = Modeller(oriented.topology, oriented.positions)
    if modeller.topology.getUnitCellDimensions() is not None:
        raise WorkError("providerMismatch", "The placeholder oriented-protein cell was not cleared")
    progress("nativeMembraneStarted", {"providerVersion": provider_version,
                                       "nativePatchSha256": patch_sha.lower(),
                                       "nativePatchMode": patch_mode})
    modeller.addMembrane(ff, lipidType=patch if patch_mode == "popc-62-109-deletion" else lipid_type,
                        membraneCenterZ=center * unit.nanometer,
                        minimumPadding=padding * unit.nanometer,
                        positiveIon="Na+", negativeIon="Cl-",
                        ionicStrength=ionic_strength * unit.molar,
                        platform=Platform.getPlatformByName("CPU"))
    verify_sha256(installed_patch, patch_sha, "nativePatchSha256")
    if patch_mode == "popc-62-109-deletion":
        verify_sha256(Path(payload["nativeSourcePatchPath"]),
                      payload["nativeSourcePatchSha256"], "nativeSourcePatchSha256")
    for raw in force_field_raw:
        asset = require_mapping(raw, "force-field asset")
        verify_sha256(work_path(directory, asset.get("path"), "force-field path"),
                      require_text(asset.get("sha256"), "force-field sha256"),
                      "force-field sha256")
    for species, representation in representations.items():
        for path_key, digest_key in (("templatePath", "templateSha256"),
                                     ("coordinateTemplatePath", "coordinateTemplateSha256")):
            verify_sha256(work_path(directory, representation.get(path_key), f"{species} {path_key}"),
                          require_text(representation.get(digest_key), f"{species} {digest_key}"),
                          f"{species} {digest_key}")
    progress("nativeMembraneReturned", {"atomCount": modeller.topology.getNumAtoms()})

    topology = modeller.topology
    positions = modeller.positions
    atoms = list(topology.atoms())
    protein_count = len(oriented_atoms)
    if len(atoms) > maximum_atoms or len(atoms) != len(positions):
        raise WorkError("resourceRefused", "Native membrane exceeds the declared atom count or lacks coordinates")
    if _full_atom_sequence(topology)[:protein_count] != _full_atom_sequence(oriented.topology):
        raise WorkError("providerMismatch", "Native builder changed prepared protein atom identity or order")
    if ({pair for pair in _bond_indices(topology) if pair[0] < protein_count or pair[1] < protein_count} !=
            _bond_indices(approved_graph)):
        raise WorkError("providerMismatch", "Native builder changed approved protein bonds or added a cross-bond")
    points = [tuple(float(value) for value in point.value_in_unit(unit.angstrom)) for point in positions]
    if any(not all(math.isfinite(value) for value in point) for point in points):
        raise WorkError("nonfiniteObservation", "Native membrane contains a nonfinite coordinate")
    protein_deviation = max(math.dist(original, actual) for original, actual in
                            zip(oriented_points, points[:protein_count]))
    if protein_deviation > 1e-6:
        raise WorkError("providerMismatch", "Native builder moved the approved oriented protein")
    box = topology.getPeriodicBoxVectors()
    if box is None:
        raise WorkError("invalidGeometry", "Native membrane returned no periodic cell")
    vectors = [[float(value) for value in vector.value_in_unit(unit.angstrom)] for vector in box]
    if any(not math.isfinite(value) for vector in vectors for value in vector) or any(
            abs(vectors[i][j]) > 1e-6 for i in range(3) for j in range(3) if i != j):
        raise WorkError("invalidGeometry", "Native route requires a finite orthorhombic cell")
    cell = [vectors[i][i] for i in range(3)]
    if any(length <= 0 or length > maximum_cell for length in cell):
        raise WorkError("invalidGeometry", "Native cell is outside the declared finite bounds")
    extent = _coordinate_extent(oriented_points)
    gaps = [cell[0] - (extent[1] - extent[0]), cell[1] - (extent[3] - extent[2]),
            cell[2] - (extent[5] - extent[4])]
    if any(gap < 20 * padding - 1e-4 for gap in gaps):
        raise WorkError("invalidGeometry", "Native cell lacks required periodic image separation for the protein")

    # One pass over bonds gives exact per-residue molecular graphs without a
    # quadratic walk over the tens of thousands of generated residues.
    bonds_by_residue: dict[int, set[tuple[int, int]]] = defaultdict(set)
    local_index = {atom.index: index for residue in topology.residues()
                   for index, atom in enumerate(residue.atoms())}
    for first, second in topology.bonds():
        if first.residue is second.residue:
            bonds_by_residue[first.residue.index].add(tuple(sorted(
                (local_index[first.index], local_index[second.index]))))
        elif first.index >= protein_count or second.index >= protein_count:
            raise WorkError("providerMismatch", "Native builder joined generated molecules across residues")

    species_counts: Counter[tuple[str, str, str]] = Counter()
    role_by_index: list[tuple[str, str] | None] = [None] * len(atoms)
    side_by_index: list[str | None] = [None] * len(atoms)
    generated_by_index: list[tuple[str, str, int] | None] = [None] * len(atoms)
    water_atoms: set[int] = set()
    for index, item in enumerate(mapped_atoms):
        role_by_index[index] = ("protein", require_text(item.get("atomRole"), "protein atomRole"))
    for residue in topology.residues():
        residue_atoms = list(residue.atoms())
        if not residue_atoms or residue_atoms[0].index < protein_count:
            if any(atom.index >= protein_count for atom in residue_atoms):
                raise WorkError("providerMismatch", "A generated atom shares a prepared-protein residue")
            continue
        species = {native_residue_name: lipid_type, "HOH": "HOH", "NA": "NA", "CL": "CL"}.get(residue.name)
        if species is None:
            raise WorkError("providerMismatch", f"Native builder generated unsupported residue {residue.name}")
        actual = _native_molecule_matches(residue, references[species], species, positions,
                                          representations[species].get("stereoChecks"),
                                          bonds_by_residue[residue.index], reference_bonds[species])
        if species == lipid_type:
            side = _native_lipid_side(actual, points, 10 * center, lipid_type)
            role, atom_role = "lipid", "body"
        elif species == "HOH":
            side = "upper" if points[actual[0].index][2] >= 10 * center else "lower"
            role, atom_role = "water", "body"
            water_atoms.update(atom.index for atom in actual)
        else:
            side = "upper" if points[actual[0].index][2] >= 10 * center else "lower"
            role = "positiveIon" if species == "NA" else "negativeIon"
            atom_role = "body"
        species_counts[(role, side, species)] += 1
        for atom in actual:
            role_by_index[atom.index] = ("ion" if role in {"positiveIon", "negativeIon"} else role,
                                        "head" if species == lipid_type and atom.name == "P" else atom_role)
            side_by_index[atom.index] = side
            generated_by_index[atom.index] = (species, role, residue.index)
    if (any(role is None for role in role_by_index) or
            species_counts[("lipid", "upper", lipid_type)] == 0 or
            species_counts[("lipid", "lower", lipid_type)] == 0 or
            species_counts[("lipid", "upper", lipid_type)] != species_counts[("lipid", "lower", lipid_type)]):
        raise WorkError("providerMismatch", f"Native result lacks complete pure symmetric {lipid_type} leaflets")

    system = _native_system_settings(ff, topology, payload.get("systemSettings"), cell, water_atoms)
    from openmm import NonbondedForce
    nonbonded = next(force for force in (system.getForce(index)
                       for index in range(system.getNumForces())) if isinstance(force, NonbondedForce))
    observed_charges = [float(nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge))
                        for index in range(len(atoms))]
    if any(not math.isfinite(value) for value in observed_charges):
        raise WorkError("nonfiniteObservation", "A parameterized atom charge is nonfinite")
    for residue in topology.residues():
        residue_atoms = list(residue.atoms())
        if residue_atoms[0].index < protein_count:
            continue
        species = generated_by_index[residue_atoms[0].index][0]
        charge = sum(observed_charges[atom.index] for atom in residue_atoms)
        declared = require_number(representations[species].get("netChargeElementary"),
                                  f"{species} netChargeElementary")
        if abs(charge - declared) > 1e-4:
            raise WorkError("providerMismatch", f"Generated {species} charge differs from selected chemistry")
    if abs(sum(observed_charges[:protein_count]) - protein_charge) > 1e-4:
        raise WorkError("providerMismatch", "Protein charge changed in the combined parameterized System")
    charge = sum(observed_charges)
    # OpenMM's PDB reader spells an absent insertion code as one blank; the
    # mmCIF reader returns the same absence as an empty string.  Canonicalize
    # that one field before serializing the bonded stage identity.
    for residue in topology.residues():
        if residue.insertionCode == " ":
            residue.insertionCode = ""
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator)
    context.setPositions(positions)
    initial = context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
    energy = initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        raise WorkError("nonfiniteObservation", "Native parameterized candidate has nonfinite initial energy")

    topology_path = directory / "constructed-topology.cif"
    with topology_path.open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(topology, positions, stream, keepIds=True)
    topology_json = directory / "constructed-topology.json"
    topology_json.write_text(json.dumps(_topology_data(topology), separators=(",", ":")), encoding="utf-8")
    readback = PDBxFile(str(topology_path))
    if (_full_atom_sequence(readback.topology) != _full_atom_sequence(topology) or
            len(readback.positions) != len(atoms)):
        raise WorkError("providerMismatch", "Native mmCIF readback lost ordered atom identities")
    readback_deviation = max(math.dist(actual, tuple(float(value) for value in exported.value_in_unit(unit.angstrom)))
                             for actual, exported in zip(points, readback.positions))
    if readback_deviation > 0.0002:
        raise WorkError("providerMismatch", "Native mmCIF coordinate readback exceeds its precision bound")
    system_path = directory / "constructed-system.xml"
    state_path = directory / "constructed-state.xml"
    system_path.write_text(XmlSerializer.serialize(system), encoding="utf-8")
    state_path.write_text(XmlSerializer.serialize(initial), encoding="utf-8")

    correspondence_atoms = []
    result_ids = set()
    for index, atom in enumerate(atoms):
        role, atom_role = role_by_index[index]
        if index < protein_count:
            source = mapped_atoms[index]
            result_id = "protein:" + json.dumps([index, atom.name, atom.residue.name], separators=(",", ":"))
            provenance = {"sourceAtomId": source.get("sourceAtomId"), "role": source["role"],
                          "sourceResidue": source.get("sourceResidue"),
                          "approvedChangeId": source.get("approvedChangeId"),
                          "generatedSpeciesId": None, "generatedComponentRole": None}
        else:
            species, component_role, residue_index = generated_by_index[index]
            result_id = "component:" + json.dumps(["OpenMM-native", residue_index, index, atom.name],
                                                   separators=(",", ":"))
            provenance = {"sourceAtomId": None, "role": "generated", "sourceResidue": None,
                          "approvedChangeId": None, "generatedSpeciesId": species,
                          "generatedComponentRole": component_role}
        if result_id in result_ids:
            raise WorkError("correspondenceFailed", "Native candidate has duplicate result atom identities")
        result_ids.add(result_id)
        correspondence_atoms.append({"resultAtomIndex": index, "resultAtomId": result_id,
                                     "moleculeRole": role, "atomRole": atom_role,
                                     "physicalSide": side_by_index[index],
                                     "element": atom.element.symbol, **provenance})
    correspondence_path = directory / "constructed-correspondence.json"
    correspondence = {"sourceId": sha256(oriented_path), "resultId": sha256(topology_path),
                      "atoms": correspondence_atoms, "complete": len(correspondence_atoms) == len(atoms)}
    correspondence_path.write_text(json.dumps(correspondence, separators=(",", ":")), encoding="utf-8")
    local_state = observe_local_state(topology, initial.getPositions(), correspondence,
                                      payload.get("localObservationSpec"))
    del context, integrator
    progress("constructedCoordinatesObserved", {"atomCount": len(atoms),
                                                "actualCellAngstrom": cell})
    counts = [{"role": role, "physicalSide": side, "speciesId": species, "count": count}
              for (role, side, species), count in sorted(species_counts.items())]
    return {"artifacts": [artifact(directory, topology_path, "topologyCif"),
                          artifact(directory, topology_json, "topologyJson"),
                          artifact(directory, system_path, "systemXml"),
                          artifact(directory, state_path, "stateXml"),
                          artifact(directory, correspondence_path, "correspondenceJson")],
            "observations": {"atomCount": len(atoms), "speciesCounts": counts,
                             "actualCellAngstrom": cell,
                             "waterCount": sum(count for (role, _, _), count in species_counts.items()
                                               if role == "water"),
                             "positiveIonCount": sum(count for (role, _, _), count in species_counts.items()
                                                     if role == "positiveIon"),
                             "negativeIonCount": sum(count for (role, _, _), count in species_counts.items()
                                                     if role == "negativeIon"),
                             "netChargeElementary": charge,
                             "proteinNetChargeElementary": protein_charge,
                             "maximumProteinCoordinateDeviationAngstrom": protein_deviation,
                             "proteinIdentityAndBondsPreserved": True,
                             "proteinPeriodicImageGapsAngstrom": gaps,
                             "nativePatchSha256": patch_sha.lower(),
                             "nativePatchMode": patch_mode,
                             "nativeSourcePatchSha256": payload.get("nativeSourcePatchSha256"),
                             "initialPotentialEnergyKjMol": energy,
                             "correspondedResultAtomCount": len(correspondence_atoms),
                             "contactWarnings": local_state["limitations"],
                             "geometryWarnings": local_state["limitations"],
                             "parameterWarnings": [], "localState": local_state},
            "provider": _provider("OpenMM Modeller.addMembrane", provider_version)}
