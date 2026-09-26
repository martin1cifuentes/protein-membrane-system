"""Owner-neutral local representation and stage mechanics used by established boundaries.

These helpers reconstruct and measure provider artifacts; they do not decide product
qualification or establish a separate scientific authority.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from ProteinInMembraneSystem.worker.exchange import (WorkError, require_integer, require_mapping,
                                        require_text, verify_sha256, work_path)

def _force_field_files(directory: Path, values: Any) -> list[str]:
    import xml.etree.ElementTree as ET

    if not isinstance(values, list) or not values:
        raise WorkError("missingPolicy", "An explicit nonempty force-field asset list is required")
    paths = []
    identities: dict[str, tuple[str, str, str, str]] = {}
    hash_identities: dict[str, tuple[str, str, str, str]] = {}
    for raw in values:
        asset = require_mapping(raw, "forceFieldFiles asset")
        identity = tuple(require_text(asset.get(field), f"force-field {field}")
                         for field in ("id", "version", "family", "sha256"))
        path = work_path(directory, asset.get("path"), "force-field path")
        verify_sha256(path, identity[3], "force-field sha256")
        existing = identities.get(str(path))
        if existing is not None:
            if existing != identity:
                raise WorkError("inputMismatch", "One force-field XML path has conflicting declared identities")
            continue
        same_bytes = hash_identities.get(identity[3].lower())
        if same_bytes is not None:
            if same_bytes != identity:
                raise WorkError("inputMismatch", "Identical force-field bytes have conflicting declared identities")
            identities[str(path)] = identity
            continue
        identities[str(path)] = identity
        hash_identities[identity[3].lower()] = identity
        paths.append(str(path))
    known = set(paths)
    for raw in paths:
        try:
            includes = ET.parse(raw).getroot().findall("Include")
        except ET.ParseError as exc:
            raise WorkError("invalidForceField", "A force-field asset is not parseable XML") from exc
        for included in includes:
            relative = included.attrib.get("file")
            if not relative:
                raise WorkError("invalidForceField", "A force-field Include has no file name")
            target = (Path(raw).parent / relative).resolve()
            if directory not in target.parents or str(target) not in known:
                raise WorkError("missingPolicy", "An included force-field XML lacks a staged, hashed asset",
                                {"includingAsset": raw, "includedFile": relative})
    return paths


def _provider(name: str, version: str | None = None) -> dict[str, str]:
    import importlib.metadata

    return {"name": name, "version": version or importlib.metadata.version("openmm")}




def _copy_molecule(target: Any, source: Any, *, chain_id: str | None = None) -> list[Any]:
    """Copy exact atom/bond order from an OpenMM topology into another."""
    from openmm.app import Topology

    mapping = {}
    for chain in source.chains():
        out_chain = target.addChain(chain_id)
        for residue in chain.residues():
            out_residue = target.addResidue(residue.name, out_chain, id=residue.id,
                                            insertionCode=residue.insertionCode)
            for atom in residue.atoms():
                mapping[atom] = target.addAtom(atom.name, atom.element, out_residue, id=atom.id)
    for bond in source.bonds():
        target.addBond(mapping[bond[0]], mapping[bond[1]],
                       type=getattr(bond, "type", None), order=getattr(bond, "order", None))
    return [mapping[a] for a in source.atoms()]


def _coordinate_file(path: Path) -> Any:
    from openmm.app import PDBFile, PDBxFile

    return PDBxFile(str(path)) if path.suffix.lower() in {".cif", ".mmcif"} else PDBFile(str(path))


def _topology_data(topology: Any) -> dict[str, Any]:
    from openmm import unit

    chains, residues, atoms = [], [], []
    for chain in topology.chains():
        chains.append({"id": chain.id})
        for residue in chain.residues():
            residues.append({"chainIndex": chain.index, "id": residue.id,
                             "name": residue.name, "insertionCode": residue.insertionCode})
            for atom in residue.atoms():
                atoms.append({"residueIndex": residue.index, "id": atom.id,
                              "name": atom.name,
                              "element": atom.element.symbol if atom.element else None})
    box = topology.getPeriodicBoxVectors()
    bonds = sorted(({"atomIndices": list(sorted((bond[0].index, bond[1].index))),
                     "order": getattr(bond, "order", None)} for bond in topology.bonds()),
                   key=lambda bond: bond["atomIndices"])
    return {"formatVersion": 1, "chains": chains, "residues": residues, "atoms": atoms,
            "bonds": bonds,
            "boxVectorsAngstrom": ([[float(value) for value in vector.value_in_unit(unit.angstrom)]
                                    for vector in box] if box is not None else None)}


def _read_topology_data(path: Path) -> Any:
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    data = require_mapping(json.loads(path.read_text(encoding="utf-8")), "topology sidecar")
    if data.get("formatVersion") != 1:
        raise WorkError("inputMismatch", "Unknown topology sidecar format")
    topology = Topology()
    chains = [topology.addChain(require_text(item.get("id"), "chain id"))
              for item in data.get("chains", [])]
    residues = []
    for item in data.get("residues", []):
        chain_index = require_integer(item.get("chainIndex"), "chainIndex")
        if chain_index >= len(chains):
            raise WorkError("inputMismatch", "Topology sidecar refers to an absent chain")
        residues.append(topology.addResidue(require_text(item.get("name"), "residue name"),
                                            chains[chain_index], id=str(item.get("id")),
                                            insertionCode=str(item.get("insertionCode") or "")))
    atoms = []
    for item in data.get("atoms", []):
        residue_index = require_integer(item.get("residueIndex"), "residueIndex")
        if residue_index >= len(residues):
            raise WorkError("inputMismatch", "Topology sidecar refers to an absent residue")
        symbol = require_text(item.get("element"), "atom element")
        atoms.append(topology.addAtom(require_text(item.get("name"), "atom name"),
                                      element.get_by_symbol(symbol), residues[residue_index],
                                      id=str(item.get("id"))))
    for bond in data.get("bonds", []):
        bond = require_mapping(bond, "topology bond")
        pair = bond.get("atomIndices")
        if not isinstance(pair, list) or len(pair) != 2 or any(
                not isinstance(index, int) or index < 0 or index >= len(atoms) for index in pair):
            raise WorkError("inputMismatch", "Topology sidecar contains an invalid bond")
        order = bond.get("order")
        if order is not None and (not isinstance(order, int) or order < 1):
            raise WorkError("inputMismatch", "Topology sidecar contains an invalid bond order")
        topology.addBond(atoms[pair[0]], atoms[pair[1]], order=order)
    box = data.get("boxVectorsAngstrom")
    if box is not None:
        if not isinstance(box, list) or len(box) != 3 or any(not isinstance(vector, list) or len(vector) != 3 for vector in box):
            raise WorkError("inputMismatch", "Topology sidecar has an invalid periodic cell")
        topology.setPeriodicBoxVectors(tuple(Vec3(*(float(value) / 10 for value in vector)) * unit.nanometer
                                             for vector in box))
    rebuilt = _topology_data(topology)
    if any(rebuilt.get(key) != data.get(key) for key in
           ("formatVersion", "chains", "residues", "atoms", "bonds")):
        raise WorkError("inputMismatch", "Topology sidecar did not reconstruct identities and bonds exactly")
    if box is not None and any(abs(rebuilt["boxVectorsAngstrom"][i][axis] - box[i][axis]) > 1e-8
                               for i in range(3) for axis in range(3)):
        raise WorkError("inputMismatch", "Topology sidecar did not reconstruct its periodic cell")
    return topology


def _nonbonded_charge(system: Any) -> float:
    from openmm import NonbondedForce, unit

    forces = [system.getForce(index) for index in range(system.getNumForces())]
    nonbonded = [force for force in forces if isinstance(force, NonbondedForce)]
    if len(nonbonded) != 1:
        raise WorkError("unobservedCharge", "Exactly one parameterized NonbondedForce is needed to observe charge")
    charge = sum(nonbonded[0].getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
                 for index in range(system.getNumParticles()))
    if not math.isfinite(charge):
        raise WorkError("nonfiniteObservation", "Parameterized net charge is not finite")
    return charge




def _full_atom_sequence(topology: Any) -> list[tuple[str, str, str, str, str, str]]:
    return [(atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode,
             atom.residue.name, atom.name, atom.element.symbol if atom.element else "")
            for atom in topology.atoms()]


def _bond_indices(topology: Any) -> set[tuple[int, int]]:
    return {tuple(sorted((a.index, b.index))) for a, b in topology.bonds()}


def _system_bonds_match(topology: Any, system: Any) -> bool:
    """Compare explicit topology bonds with realized System bond/constraint terms.

    Rigid water adds an H–H distance constraint even though the two hydrogen
    atoms are not directly bonded; that one same-residue constraint is allowed.
    This is a check of the selected Amber/OpenMM route, not a general parser for
    every possible force-field functional form.
    """
    from openmm import CustomBondForce, HarmonicBondForce

    atoms = list(topology.atoms())
    if len(atoms) != system.getNumParticles():
        return False
    declared = _bond_indices(topology)
    realized = set()
    for force_index in range(system.getNumForces()):
        force = system.getForce(force_index)
        if isinstance(force, (HarmonicBondForce, CustomBondForce)):
            for bond_index in range(force.getNumBonds()):
                first, second = force.getBondParameters(bond_index)[:2]
                realized.add(tuple(sorted((int(first), int(second)))))
    for constraint_index in range(system.getNumConstraints()):
        first, second, _ = system.getConstraintParameters(constraint_index)
        pair = tuple(sorted((int(first), int(second))))
        if pair not in declared:
            a, b = (atoms[index] for index in pair)
            if not (a.element.symbol == b.element.symbol == "H" and a.residue is b.residue):
                return False
        else:
            realized.add(pair)
    return realized == declared


def _load_stage(directory: Path, payload: dict[str, Any], state_key: str):
    from types import SimpleNamespace
    from openmm import XmlSerializer
    from openmm.app import PDBxFile

    topology_path = work_path(directory, payload.get("topologyCifPath"), "topologyCifPath")
    topology_json = work_path(directory, payload.get("topologyJsonPath"), "topologyJsonPath")
    system_path = work_path(directory, payload.get("systemXmlPath"), "systemXmlPath")
    state_path = work_path(directory, payload.get(state_key), state_key)
    state_hash_key = state_key.removesuffix("Path") + "Sha256"
    for path, hash_key in ((topology_path, "topologyCifSha256"),
                           (topology_json, "topologyJsonSha256"),
                           (system_path, "systemXmlSha256"),
                           (state_path, state_hash_key)):
        verify_sha256(path, require_text(payload.get(hash_key), hash_key), hash_key)
    cif = PDBxFile(str(topology_path))
    topology = _read_topology_data(topology_json)
    system = XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    state = XmlSerializer.deserialize(state_path.read_text(encoding="utf-8"))
    count = len(list(topology.atoms()))
    if count != system.getNumParticles() or len(state.getPositions()) != count:
        raise WorkError("inputMismatch", "Topology, System, and State disagree on ordered particle count")
    if topology.getPeriodicBoxVectors() is None:
        raise WorkError("inputMismatch", "Full-system stage topology has no periodic cell")
    if _full_atom_sequence(cif.topology) != _full_atom_sequence(topology):
        raise WorkError("inputMismatch", "mmCIF and bond-topology sidecar disagree on ordered atom identities")
    if len(cif.positions) != count:
        raise WorkError("inputMismatch", "mmCIF coordinate count differs from full bonded topology")
    return SimpleNamespace(topology=topology, positions=cif.positions,
                           coordinate_box_vectors=cif.topology.getPeriodicBoxVectors()), system, state


def _rms_force(state: Any, particle_count: int) -> float:
    from openmm import unit

    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    value = math.sqrt(float((forces * forces).sum()) / (3 * particle_count))
    if not math.isfinite(value):
        raise WorkError("nonfiniteObservation", "Stage force observation is nonfinite")
    return value


def _stage_measurements(state: Any, particle_count: int) -> list[dict[str, Any]]:
    from openmm import unit

    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not (math.isfinite(potential) and math.isfinite(kinetic)):
        raise WorkError("nonfiniteObservation", "Stage energy observation is nonfinite")
    values = [
        {"name": "potentialEnergy", "value": potential, "unit": "kJ/mol", "scope": "wholeSystem"},
        {"name": "kineticEnergy", "value": kinetic, "unit": "kJ/mol", "scope": "wholeSystem"},
        {"name": "rmsForce", "value": _rms_force(state, particle_count), "unit": "kJ/mol/nm", "scope": "wholeSystem"},
    ]
    vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    for name, vector in zip(("cellX", "cellY", "cellZ"), vectors):
        length = math.sqrt(float((vector * vector).sum()))
        values.append({"name": name, "value": length, "unit": "angstrom", "scope": "periodicCell"})
    return values


def _state_context(system: Any, state: Any, integrator: Any):
    from openmm import Context

    context = Context(system, integrator)
    context.setState(state)
    return context


def _final_state(context: Any):
    return context.getState(getPositions=True, getVelocities=True, getForces=True,
                            getEnergy=True, enforcePeriodicBox=False)


def _write_stage_coordinates(path: Path, topology: Any, positions: Any, box_vectors: Any) -> None:
    from openmm.app import PDBxFile

    # The barostat may have changed the periodic cell since the constructed
    # topology was written. The coordinate artifact must describe this stage's
    # actual State cell, not the original construction cell.
    topology.setPeriodicBoxVectors(box_vectors)
    with path.open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(topology, positions, stream, keepIds=True)
