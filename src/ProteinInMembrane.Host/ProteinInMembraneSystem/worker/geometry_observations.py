"""Measured protein geometry from an exact selected or prepared atom graph.

These are bounded measurements, not universal geometric validity thresholds.
The host's identified structural-assessment policy interprets them.
"""

from __future__ import annotations

import math
import heapq
from collections import defaultdict
from typing import Any

from ProteinInMembraneSystem.worker.exchange import WorkError, require_integer, require_mapping, require_number
from ProteinInMembraneSystem.worker.protein_identity import _atom_address, _icode

_KINDS = ("covalentBond", "chainContinuity", "nonbondedDistance")


def unavailable_geometry(spec: Any, reason: str) -> dict[str, Any]:
    names = spec.get("requiredKinds", []) if isinstance(spec, dict) else []
    return {"standing": "Unavailable",
            "kinds": [{"kind": name, "standing": "Unavailable", "eligibleCount": 0,
                       "measuredCount": 0, "minimumDistanceAngstrom": None,
                       "maximumDistanceAngstrom": None, "unavailableReason": reason}
                      for name in names],
            "locatedDistances": [], "limitations": [reason]}


def observe_protein_geometry(topology: Any, positions: Any, spec: Any,
                             chain_map: list[dict[str, str]] | None, model_index: int,
                             atom_addresses: list[dict[str, Any]] | None = None,
                             periodic_lengths: tuple[float, float, float] | None = None) -> dict[str, Any]:
    from openmm import unit
    policy = require_mapping(spec, "geometrySpec")
    names = policy.get("requiredKinds")
    if not isinstance(names, list) or not names or len(set(names)) != len(names) or any(name not in _KINDS for name in names):
        raise WorkError("unsupportedObservation", "Protein geometry requests an unsupported or repeated kind")
    radius = require_number(policy.get("neighborSearchRadiusAngstrom"), "neighborSearchRadiusAngstrom", 0.000001)
    hops = require_integer(policy.get("excludedBondHops"), "excludedBondHops")
    maximum = require_integer(policy.get("maximumReportedPairs"), "maximumReportedPairs", 1)
    if hops > 8 or maximum > 10_000:
        raise WorkError("resourceRefused", "Protein geometry graph/report bounds exceed the supported observation size")
    radii = require_mapping(policy.get("atomRadiusByElementAngstrom"), "atomRadiusByElementAngstrom")
    atoms = list(topology.atoms())
    xyz = positions.value_in_unit(unit.angstrom)
    if len(xyz) != len(atoms) or any(not all(math.isfinite(float(value)) for value in point) for point in xyz):
        return unavailable_geometry(policy, "Atom positions are incomplete or nonfinite.")
    points = [tuple(float(value) for value in point) for point in xyz]
    if atom_addresses is not None and len(atom_addresses) != len(atoms):
        return unavailable_geometry(policy, "Protein atom-address mapping is incomplete for the stage graph.")
    if periodic_lengths is not None and any(not math.isfinite(value) or value <= 2 * radius
                                            for value in periodic_lengths):
        return unavailable_geometry(policy, "Protein geometry periodic cell is too small for unambiguous near-pair measurement.")

    def distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
        differences = [first[axis] - second[axis] for axis in range(3)]
        if periodic_lengths is not None:
            differences = [value - round(value / periodic_lengths[axis]) * periodic_lengths[axis]
                           for axis, value in enumerate(differences)]
        return math.sqrt(sum(value * value for value in differences))

    addresses = []
    element_radii = []
    for atom in atoms:
        residue = atom.residue
        key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode), atom.name)
        addresses.append(atom_addresses[atom.index] if atom_addresses is not None else
                         _atom_address(key, chain_map, model_index))
        symbol = atom.element.symbol if atom.element else ""
        if symbol not in radii:
            return unavailable_geometry(policy, f"No declared observation radius exists for element {symbol!r}.")
        element_radii.append(require_number(radii[symbol], f"radius {symbol}", 0.000001))

    bonds = []
    neighbors: list[set[int]] = [set() for _ in atoms]
    for first, second in topology.bonds():
        i, j = first.index, second.index
        bonds.append((i, j))
        neighbors[i].add(j)
        neighbors[j].add(i)
    excluded: list[set[int]] = []
    for i in range(len(atoms)):
        reached = {i}
        frontier = {i}
        for _ in range(hops):
            frontier = set().union(*(neighbors[index] for index in frontier)) - reached if frontier else set()
            reached.update(frontier)
        excluded.append(reached)

    stats: dict[str, dict[str, Any]] = {name: {"count": 0, "minimum": None,
                                               "maximum": None, "reported": []} for name in _KINDS}

    def add(name: str, distance: float, first: int, second: int) -> None:
        result = stats[name]
        result["count"] += 1
        result["minimum"] = distance if result["minimum"] is None else min(distance, result["minimum"])
        result["maximum"] = distance if result["maximum"] is None else max(distance, result["maximum"])
        entry = (-distance, first, second)
        if len(result["reported"]) < maximum:
            heapq.heappush(result["reported"], entry)
        elif entry > result["reported"][0]:
            heapq.heapreplace(result["reported"], entry)

    for i, j in bonds:
        add("covalentBond", distance(points[i], points[j]), i, j)
    for chain in topology.chains():
        residues = list(chain.residues())
        for previous, following in zip(residues, residues[1:]):
            first = next((atom for atom in previous.atoms() if atom.name == "C"), None)
            second = next((atom for atom in following.atoms() if atom.name == "N"), None)
            if first is None or second is None:
                return unavailable_geometry(policy, "A chain-continuity atom was not present in the selected graph.")
            add("chainContinuity", distance(points[first.index], points[second.index]),
                first.index, second.index)
    # Cell-linked candidate search covers every pair inside the declared
    # radius without making a quadratic all-pairs claim for large proteins.
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    bin_counts = (tuple(max(1, math.floor(length / radius)) for length in periodic_lengths)
                  if periodic_lengths is not None else None)
    for i, point in enumerate(points):
        if bin_counts is None:
            cell = tuple(math.floor(value / radius) for value in point)
        else:
            cell = tuple(min(bin_counts[axis] - 1,
                             math.floor((point[axis] % periodic_lengths[axis]) /
                                        periodic_lengths[axis] * bin_counts[axis])) for axis in range(3))
        neighbor_cells = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor_cell = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    if bin_counts is not None:
                        neighbor_cell = tuple(neighbor_cell[axis] % bin_counts[axis] for axis in range(3))
                    neighbor_cells.add(neighbor_cell)
        for neighbor_cell in neighbor_cells:
            for j in cells.get(neighbor_cell, []):
                if j in excluded[i]:
                    continue
                observed_distance = distance(point, points[j])
                if observed_distance <= radius:
                    add("nonbondedDistance", observed_distance, j, i)
        cells[cell].append(i)
    kinds = []
    located = []
    for name in names:
        observed = stats[name]
        count = observed["count"]
        kinds.append({"kind": name, "standing": "Observed" if count else "NotApplicable",
                      "eligibleCount": count, "measuredCount": count,
                      "minimumDistanceAngstrom": observed["minimum"],
                      "maximumDistanceAngstrom": observed["maximum"],
                      "unavailableReason": None if count else "No applicable pair was present in the declared observation scope."})
        for negative_distance, i, j in sorted(observed["reported"], reverse=True):
            distance = -negative_distance
            located.append({"kind": name, "first": addresses[i], "second": addresses[j],
                            "distanceAngstrom": distance,
                            "radiusSumAngstrom": element_radii[i] + element_radii[j] if name == "nonbondedDistance" else None})
    return {"standing": "Observed", "kinds": kinds,
            "locatedDistances": located, "limitations": []}


def observe_stage_protein_geometry(topology: Any, positions: Any,
                                   correspondence: Any, spec: Any) -> dict[str, Any]:
    """Measure only exact protein atoms of a completed stage, in source identity."""
    from openmm import Vec3, unit
    from openmm.app import Topology

    mapping = require_mapping(correspondence, "stage correspondence")
    raw = mapping.get("atoms")
    source_atoms = list(topology.atoms())
    if mapping.get("complete") is not True or not isinstance(raw, list) or len(raw) != len(source_atoms):
        return unavailable_geometry(spec, "Completed-stage correspondence is incomplete.")
    if sorted(item.get("resultAtomIndex") for item in raw if isinstance(item, dict)) != list(range(len(source_atoms))):
        return unavailable_geometry(spec, "Completed-stage atom order cannot be mapped exactly.")
    raw.sort(key=lambda item: item["resultAtomIndex"])
    selected = [index for index, item in enumerate(raw) if item.get("moleculeRole") == "protein"]
    if not selected:
        return unavailable_geometry(spec, "The completed stage has no identified protein atoms.")
    source_xyz = positions.value_in_unit(unit.angstrom)
    if len(source_xyz) != len(source_atoms):
        return unavailable_geometry(spec, "Completed-stage positions do not match source atom count.")
    selected_topology = Topology()
    selected_chains = {}
    selected_residues = {}
    selected_atoms = {}
    addresses = []
    selected_points = []
    for old_index in selected:
        atom = source_atoms[old_index]
        identity = raw[old_index]
        source_residue = identity.get("sourceResidue")
        if (not isinstance(source_residue, dict) or
                any(name not in source_residue for name in
                    ("model", "chain", "residue", "insertionCode", "copyId")) or
                identity.get("element") != (atom.element.symbol if atom.element else "")):
            return unavailable_geometry(spec, "A completed-stage protein atom lacks exact source identity.")
        chain_key = atom.residue.chain.index
        if chain_key not in selected_chains:
            selected_chains[chain_key] = selected_topology.addChain(atom.residue.chain.id)
        residue_key = atom.residue.index
        if residue_key not in selected_residues:
            selected_residues[residue_key] = selected_topology.addResidue(
                atom.residue.name, selected_chains[chain_key],
                id=atom.residue.id, insertionCode=atom.residue.insertionCode)
        selected_atoms[old_index] = selected_topology.addAtom(
            atom.name, atom.element, selected_residues[residue_key], id=atom.id)
        addresses.append({"residue": source_residue, "atomName": atom.name})
        selected_points.append(Vec3(*[float(value) for value in source_xyz[old_index]]))
    for first, second in topology.bonds():
        if first.index in selected_atoms and second.index in selected_atoms:
            selected_topology.addBond(selected_atoms[first.index], selected_atoms[second.index])
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return unavailable_geometry(spec, "Completed-stage protein geometry has no periodic cell.")
    box = [[float(value) for value in vector.value_in_unit(unit.angstrom)] for vector in vectors]
    if any(abs(box[i][j]) > 1e-6 for i in range(3) for j in range(3) if i != j):
        return unavailable_geometry(spec, "Stage protein-geometry measurement supports orthorhombic cells only.")
    lengths = tuple(box[i][i] for i in range(3))
    return observe_protein_geometry(selected_topology,
                                    unit.Quantity(selected_points, unit.angstrom), spec,
                                    None, 0, atom_addresses=addresses,
                                    periodic_lengths=lengths)
