"""Owner-local completed-stage export and independent read-back observations.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import WorkError, artifact, require_number
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _provider, _topology_data, _read_topology_data, _full_atom_sequence,
    _bond_indices, _system_bonds_match, _load_stage, _write_stage_coordinates)


def _gemmi_atoms(structure: Any) -> list[tuple[str, int, str, str, str, str, float, float, float]]:
    result = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    result.append((chain.name, int(residue.seqid.num), str(residue.seqid.icode).strip(),
                                   residue.name, atom.name, atom.element.name,
                                   atom.pos.x, atom.pos.y, atom.pos.z))
    return result


def verify_export(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    import gemmi
    from openmm import unit
    from openmm.app import PDBxFile

    stage, system, state = _load_stage(directory, payload, "stateXmlPath")
    tolerance = require_number(payload.get("coordinateReadBackToleranceAngstrom"),
                               "coordinateReadBackToleranceAngstrom", 0)
    cell_length_tolerance = require_number(payload.get("cellLengthReadBackToleranceAngstrom"),
                                           "cellLengthReadBackToleranceAngstrom", 0)
    cell_angle_tolerance = require_number(payload.get("cellAngleReadBackToleranceDegrees"),
                                          "cellAngleReadBackToleranceDegrees", 0)
    stage.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    mmcif = directory / "prepared-stage.cif"
    _write_stage_coordinates(mmcif, stage.topology, state.getPositions(), state.getPeriodicBoxVectors())
    openmm_readback = PDBxFile(str(mmcif))
    after = gemmi.read_structure(str(mmcif))
    if len(after) != 1:
        raise WorkError("invalidStructure", "Export stage must contain exactly one coordinate model")
    source_atoms = []
    for atom, position in zip(stage.topology.atoms(), state.getPositions()):
        x, y, z = position.value_in_unit(unit.angstrom)
        source_atoms.append((atom.residue.chain.id, int(atom.residue.id), atom.residue.insertionCode,
                             atom.residue.name, atom.name, atom.element.symbol if atom.element else "",
                             float(x), float(y), float(z)))
    exported_atoms = _gemmi_atoms(after)
    identity_before = [atom[:6] for atom in source_atoms]
    identity_after = [atom[:6] for atom in exported_atoms]
    count = len(source_atoms)
    corresponding = sum(a == b for a, b in zip(identity_before, identity_after))
    order_matched = (identity_before == identity_after and count == system.getNumParticles()
                     and _full_atom_sequence(openmm_readback.topology) == _full_atom_sequence(stage.topology))
    coordinate_deviation = (max(math.dist(a[6:], b[6:]) for a, b in zip(source_atoms, exported_atoms))
                            if source_atoms and exported_atoms else math.inf)
    vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    desired_lengths = [math.sqrt(float((vector * vector).sum())) for vector in vectors]
    def angle(first: Any, second: Any) -> float:
        cosine = float((first * second).sum()) / (
            math.sqrt(float((first * first).sum())) * math.sqrt(float((second * second).sum())))
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    desired_angles = [angle(vectors[1], vectors[2]), angle(vectors[0], vectors[2]),
                      angle(vectors[0], vectors[1])]
    cell_matched = (all(abs(actual - expected) <= cell_length_tolerance for actual, expected in
                        zip((after.cell.a, after.cell.b, after.cell.c), desired_lengths))
                    and all(abs(actual - expected) <= cell_angle_tolerance for actual, expected in
                            zip((after.cell.alpha, after.cell.beta, after.cell.gamma), desired_angles)))
    bond_file = directory / "stage-bond-graph.json"
    bond_file.write_text(json.dumps(_topology_data(stage.topology), separators=(",", ":")), encoding="utf-8")
    loaded_topology = _read_topology_data(bond_file)
    bonds_matched = (_bond_indices(loaded_topology) == _bond_indices(stage.topology)
                     and _full_atom_sequence(loaded_topology) == _full_atom_sequence(stage.topology)
                     and _system_bonds_match(loaded_topology, system))
    readback = (order_matched and corresponding == count and cell_matched and bonds_matched
                and math.isfinite(coordinate_deviation) and coordinate_deviation <= tolerance)
    warnings = []
    if not readback:
        warnings.append("Read-back did not preserve every required atom, coordinate, cell, and bond-sidecar fact.")
    progress("exportReadBackObserved", {"readBackMatched": readback})
    return {"artifacts": [artifact(directory, mmcif, "stageMmcif"),
                          artifact(directory, bond_file, "stageBondGraph")],
            "observations": {"sourceAtomCount": count, "exportedAtomCount": len(exported_atoms),
                             "correspondingElementAndResidueCount": corresponding,
                             "atomOrderMatched": order_matched, "bondsMatched": bonds_matched,
                             "cellMatched": cell_matched, "readBackMatched": readback,
                             "coordinateMaxDeviationAngstrom": coordinate_deviation if math.isfinite(coordinate_deviation) else 1e300,
                             "warnings": warnings},
            "provider": _provider("Gemmi/OpenMM export read-back",
                                   f"{gemmi.__version__}/{_provider('OpenMM')['version']}")}
