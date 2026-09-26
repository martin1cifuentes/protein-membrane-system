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


def _gemmi_atoms(atom_sites: Any) -> list[tuple[str, str, str, str, str, str, float, float, float]]:
    """Read mmCIF atom-site rows in particle order, including repeated residue IDs.

    Gemmi's Structure hierarchy groups atoms with the same chain/residue ID;
    those groups need not preserve the input row order for generated waters.
    """
    result = []
    for row in atom_sites:
        insertion = row["_atom_site.pdbx_PDB_ins_code"]
        result.append((row["_atom_site.auth_asym_id"], row["_atom_site.auth_seq_id"],
                       "" if insertion in {".", "?"} else insertion,
                       row["_atom_site.auth_comp_id"], row["_atom_site.auth_atom_id"],
                       row["_atom_site.type_symbol"],
                       float(row["_atom_site.Cartn_x"]), float(row["_atom_site.Cartn_y"]),
                       float(row["_atom_site.Cartn_z"])))
    return result


def verify_export(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    import gemmi
    import numpy as np
    from openmm import unit
    from openmm.app import PDBxFile

    stage, system, state = _load_stage(directory, payload, "stateXmlPath")
    tolerance = require_number(payload.get("coordinateReadBackToleranceAngstrom"),
                               "coordinateReadBackToleranceAngstrom", 0)
    cell_length_tolerance = require_number(payload.get("cellLengthReadBackToleranceAngstrom"),
                                           "cellLengthReadBackToleranceAngstrom", 0)
    cell_angle_tolerance = require_number(payload.get("cellAngleReadBackToleranceDegrees"),
                                          "cellAngleReadBackToleranceDegrees", 0)
    selected_xyz = np.asarray(stage.positions.value_in_unit(unit.angstrom), dtype=float)
    state_xyz = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom), dtype=float)
    selected_coordinate_deviation = (float(np.max(np.linalg.norm(selected_xyz - state_xyz, axis=1)))
                                     if selected_xyz.shape == state_xyz.shape and len(state_xyz) else math.inf)
    selected_coordinates_matched = (math.isfinite(selected_coordinate_deviation)
                                    and selected_coordinate_deviation <= tolerance)

    def cell_dimensions(vectors: Any) -> tuple[list[float], list[float]] | None:
        values = np.asarray(vectors, dtype=float)
        if values.shape != (3, 3) or not np.isfinite(values).all():
            return None
        determinant = float(np.linalg.det(values))
        if not math.isfinite(determinant) or determinant <= 0:
            return None
        lengths = np.linalg.norm(values, axis=1)
        if not np.isfinite(lengths).all() or np.any(lengths <= 0):
            return None
        pairs = ((1, 2), (0, 2), (0, 1))
        angles = [math.degrees(math.acos(max(-1.0, min(1.0,
                  float(np.dot(values[i], values[j]) / (lengths[i] * lengths[j]))))))
                  for i, j in pairs]
        return lengths.tolist(), angles

    state_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    desired_cell = cell_dimensions(state_vectors)
    selected_box = stage.coordinate_box_vectors
    selected_cell = (cell_dimensions(selected_box.value_in_unit(unit.angstrom))
                     if selected_box is not None else None)
    sidecar_box = stage.topology.getPeriodicBoxVectors()
    sidecar_vectors = ([vector.value_in_unit(unit.angstrom) for vector in sidecar_box]
                       if sidecar_box is not None else None)
    sidecar_cell = cell_dimensions(sidecar_vectors) if sidecar_vectors is not None else None
    system_vectors = [vector.value_in_unit(unit.angstrom)
                      for vector in system.getDefaultPeriodicBoxVectors()]
    system_cell = cell_dimensions(system_vectors)

    def cell_matched(actual: tuple[list[float], list[float]] | None,
                     reference: tuple[list[float], list[float]] | None) -> bool:
        return (actual is not None and reference is not None and
                all(abs(value - expected) <= cell_length_tolerance for value, expected in
                    zip(actual[0], reference[0])) and
                all(abs(value - expected) <= cell_angle_tolerance for value, expected in
                    zip(actual[1], reference[1])))

    def vectors_matched(actual: Any, reference: Any) -> bool:
        actual_values = np.asarray(actual, dtype=float)
        reference_values = np.asarray(reference, dtype=float)
        return (actual_values.shape == (3, 3) and reference_values.shape == (3, 3) and
                np.isfinite(actual_values).all() and np.isfinite(reference_values).all() and
                float(np.max(np.abs(actual_values - reference_values))) <= cell_length_tolerance)

    selected_cell_matched = cell_matched(selected_cell, desired_cell)
    sidecar_cell_matched = (cell_matched(sidecar_cell, desired_cell)
                            and cell_matched(sidecar_cell, selected_cell)
                            and vectors_matched(sidecar_vectors, state_vectors))
    system_cell_matched = (cell_matched(system_cell, desired_cell) and
                           vectors_matched(system_vectors, state_vectors))
    stage.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    mmcif = directory / "prepared-stage.cif"
    _write_stage_coordinates(mmcif, stage.topology, state.getPositions(), state.getPeriodicBoxVectors())
    openmm_readback = PDBxFile(str(mmcif))
    after = gemmi.read_structure(str(mmcif))
    if len(after) != 1:
        raise WorkError("invalidStructure", "Export stage must contain exactly one coordinate model")
    atom_sites = gemmi.cif.read_file(str(mmcif)).sole_block().find_mmcif_category("_atom_site.")
    if len({row["_atom_site.pdbx_PDB_model_num"] for row in atom_sites}) != 1:
        raise WorkError("invalidStructure", "Export atom sites must contain exactly one coordinate model")
    source_atoms = []
    for atom, position in zip(stage.topology.atoms(), state.getPositions()):
        x, y, z = position.value_in_unit(unit.angstrom)
        source_atoms.append((atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode,
                             atom.residue.name, atom.name, atom.element.symbol if atom.element else "",
                             float(x), float(y), float(z)))
    exported_atoms = _gemmi_atoms(atom_sites)
    identity_before = [atom[:6] for atom in source_atoms]
    identity_after = [atom[:6] for atom in exported_atoms]
    count = len(source_atoms)
    corresponding = sum(a == b for a, b in zip(identity_before, identity_after))
    order_matched = (identity_before == identity_after and count == system.getNumParticles()
                     and _full_atom_sequence(openmm_readback.topology) == _full_atom_sequence(stage.topology))
    output_coordinate_deviation = (max(math.dist(a[6:], b[6:]) for a, b in zip(source_atoms, exported_atoms))
                                   if source_atoms and exported_atoms else math.inf)
    coordinate_deviation = max(selected_coordinate_deviation, output_coordinate_deviation)
    exported_cell = ([after.cell.a, after.cell.b, after.cell.c],
                     [after.cell.alpha, after.cell.beta, after.cell.gamma])
    output_cell_matched = cell_matched(exported_cell, desired_cell)
    all_cells_matched = (selected_cell_matched and sidecar_cell_matched and
                         system_cell_matched and output_cell_matched)
    bond_file = directory / "stage-bond-graph.json"
    bond_file.write_text(json.dumps(_topology_data(stage.topology), separators=(",", ":")), encoding="utf-8")
    loaded_topology = _read_topology_data(bond_file)
    bonds_matched = (_bond_indices(loaded_topology) == _bond_indices(stage.topology)
                     and _full_atom_sequence(loaded_topology) == _full_atom_sequence(stage.topology)
                     and _system_bonds_match(loaded_topology, system))
    readback = (order_matched and corresponding == count and all_cells_matched and bonds_matched
                and selected_coordinates_matched and
                math.isfinite(output_coordinate_deviation) and output_coordinate_deviation <= tolerance)
    warnings = []
    if not selected_coordinates_matched:
        warnings.append("The selected stage mmCIF coordinates disagree with its identified State.")
    if not selected_cell_matched:
        warnings.append("The selected stage mmCIF cell disagrees with its identified State.")
    if not sidecar_cell_matched:
        warnings.append("The bonded topology sidecar cell disagrees with the selected stage mmCIF or State.")
    if not system_cell_matched:
        warnings.append("The identified System default periodic cell disagrees with the selected stage State.")
    if not (order_matched and corresponding == count and output_cell_matched and bonds_matched
            and math.isfinite(output_coordinate_deviation) and output_coordinate_deviation <= tolerance):
        warnings.append("Read-back did not preserve every required atom, coordinate, cell, and bond-sidecar fact.")
    progress("exportReadBackObserved", {"readBackMatched": readback})
    return {"artifacts": [artifact(directory, mmcif, "stageMmcif"),
                          artifact(directory, bond_file, "stageBondGraph")],
            "observations": {"sourceAtomCount": count, "exportedAtomCount": len(exported_atoms),
                             "correspondingElementAndResidueCount": corresponding,
                             "atomOrderMatched": order_matched, "bondsMatched": bonds_matched,
                             "cellMatched": all_cells_matched, "readBackMatched": readback,
                             "coordinateMaxDeviationAngstrom": coordinate_deviation if math.isfinite(coordinate_deviation) else 1e300,
                             "warnings": warnings},
            "provider": _provider("Gemmi/OpenMM export read-back",
                                   f"{gemmi.__version__}/{_provider('OpenMM')['version']}")}
