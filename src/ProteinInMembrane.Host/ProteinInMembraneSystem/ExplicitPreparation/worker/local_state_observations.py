"""Actual local-state geometry for constructed and completed explicit systems.

The measurements are deliberately narrower than a support verdict. In
particular, an empty contact list means no pair inside the declared search
radius, not a scientifically accepted membrane-protein assembly.
"""

from __future__ import annotations

import math
import heapq
from collections import defaultdict
from typing import Any

from ProteinInMembraneSystem.worker.exchange import WorkError, require_integer, require_mapping, require_number, require_text

_METRICS = frozenset({"minimumIntermolecularDistanceAngstrom", "upperLipidHeadMeanZAngstrom",
                      "lowerLipidHeadMeanZAngstrom", "leafletHeadSeparationAngstrom",
                      "proteinBilayerMidplaneOffsetAngstrom"})
_ROLES = frozenset({"protein", "retainedPartner", "lipid", "water", "ion"})


def _unavailable(reason: str, covered: list[str] | None = None) -> dict[str, Any]:
    return {"standing": "Unavailable", "unavailableReason": reason,
            "measurements": [], "locatedContacts": [], "coveredRolePairs": covered or [],
            "rolePairMeasurements": [],
            "limitations": [reason]}


def observe_local_state(topology: Any, positions: Any, correspondence: Any,
                        specification: Any) -> dict[str, Any]:
    from openmm import unit

    spec = require_mapping(specification, "localObservationSpec")
    requested = spec.get("requiredMetricNames")
    if not isinstance(requested, list) or not requested or len(set(requested)) != len(requested) or any(
            name not in _METRICS for name in requested):
        raise WorkError("unsupportedObservation", "Local-state policy requests an unsupported or repeated metric")
    radius = require_number(spec.get("contactSearchRadiusAngstrom"), "contactSearchRadiusAngstrom", 0.000001)
    maximum = require_integer(spec.get("maximumReportedPairs"), "maximumReportedPairs", 1)
    if maximum > 10_000:
        raise WorkError("resourceRefused", "Local contact report bound exceeds supported observation size")
    use_periodic = spec.get("usePeriodicBoundary")
    if not isinstance(use_periodic, bool):
        raise WorkError("invalidRequest", "usePeriodicBoundary must be explicit")
    declared_pairs = spec.get("contactRolePairs")
    if not isinstance(declared_pairs, list) or not declared_pairs:
        raise WorkError("invalidRequest", "At least one local-state role pair is required")
    role_pairs = []
    for item in declared_pairs:
        pair = require_mapping(item, "contact role pair")
        first = require_text(pair.get("firstMoleculeRole"), "firstMoleculeRole")
        second = require_text(pair.get("secondMoleculeRole"), "secondMoleculeRole")
        if first not in _ROLES or second not in _ROLES or (first, second) in role_pairs:
            raise WorkError("unsupportedObservation", "Contact role pair is unknown or repeated")
        role_pairs.append((first, second))
    radii = require_mapping(spec.get("atomRadiusByElementAngstrom"), "atomRadiusByElementAngstrom")
    atoms = list(topology.atoms())
    mapped = require_mapping(correspondence, "stage correspondence")
    raw_mapping = mapped.get("atoms")
    if mapped.get("complete") is not True or not isinstance(raw_mapping, list) or len(raw_mapping) != len(atoms):
        return _unavailable("Stage correspondence is incomplete for ordered atoms.")
    if sorted(item.get("resultAtomIndex") for item in raw_mapping if isinstance(item, dict)) != list(range(len(atoms))):
        return _unavailable("Stage correspondence does not identify each ordered atom exactly once.")
    raw_mapping.sort(key=lambda item: item["resultAtomIndex"])
    xyz = positions.value_in_unit(unit.angstrom)
    if len(xyz) != len(atoms) or any(not all(math.isfinite(float(value)) for value in point) for point in xyz):
        return _unavailable("A stage atom coordinate is absent or nonfinite.")
    points = [tuple(float(value) for value in point) for point in xyz]
    observed_radii = []
    roles = []
    for atom, identity in zip(atoms, raw_mapping):
        role = identity.get("moleculeRole")
        symbol = atom.element.symbol if atom.element else ""
        if role not in _ROLES or symbol != identity.get("element"):
            return _unavailable("Stage molecular role or element differs from ordered topology identity.")
        if symbol not in radii:
            return _unavailable(f"Local-state policy has no declared radius for element {symbol!r}.")
        observed_radii.append(require_number(radii[symbol], f"radius {symbol}", 0.000001))
        roles.append(role)

    lengths = None
    if use_periodic:
        vectors = topology.getPeriodicBoxVectors()
        if vectors is None:
            return _unavailable("Periodic observation was requested but the stage lacks a cell.")
        box = [[float(value) for value in vector.value_in_unit(unit.angstrom)] for vector in vectors]
        if any(abs(box[i][j]) > 1e-6 for i in range(3) for j in range(3) if i != j):
            return _unavailable("The bounded local-state search supports orthorhombic cells only.")
        lengths = [box[i][i] for i in range(3)]
        if any(not math.isfinite(length) or length <= 2 * radius for length in lengths):
            return _unavailable("Periodic cell is too small for the declared unambiguous contact search.")

    # Connected components are exact molecular memberships from the bond graph.
    parent = list(range(len(atoms)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in topology.bonds():
        a, b = find(first.index), find(second.index)
        if a != b:
            parent[b] = a
    molecules = [find(index) for index in range(len(atoms))]
    if lengths is None:
        cell_of = lambda point: tuple(math.floor(value / radius) for value in point)
        counts = None
    else:
        counts = [max(1, math.floor(length / radius)) for length in lengths]
        cell_of = lambda point: tuple(min(counts[axis] - 1,
                                          math.floor((point[axis] % lengths[axis]) / lengths[axis] * counts[axis]))
                                      for axis in range(3))
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for index, point in enumerate(points):
        cells[cell_of(point)].append(index)
    present_roles = set(roles)
    if any(first not in present_roles or second not in present_roles for first, second in role_pairs):
        return _unavailable("A declared contact role has no identified atoms in this stage.")
    pair_set = set(role_pairs)
    observed_contacts: list[tuple[float, int, int, str, str]] = []
    # These aggregates cover every declared role-pair contact inside the
    # search radius. Only locatedContacts below is capped for presentation.
    pair_measurements = {pair: {"count": 0, "minimum": None} for pair in role_pairs}

    def retain(distance: float, i: int, j: int, first: str, second: str) -> None:
        result = pair_measurements[(first, second)]
        result["count"] += 1
        result["minimum"] = distance if result["minimum"] is None else min(result["minimum"], distance)
        entry = (-distance, i, j, first, second)
        if len(observed_contacts) < maximum:
            heapq.heappush(observed_contacts, entry)
        elif entry > observed_contacts[0]:
            heapq.heapreplace(observed_contacts, entry)

    minimum_all: float | None = None
    for i, point in enumerate(points):
        cell = cell_of(point)
        neighbor_cells = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    if counts is not None:
                        neighbor = tuple(neighbor[axis] % counts[axis] for axis in range(3))
                    neighbor_cells.add(neighbor)
        for neighbor in neighbor_cells:
            for j in cells.get(neighbor, []):
                if j <= i or molecules[i] == molecules[j]:
                    continue
                delta = [point[axis] - points[j][axis] for axis in range(3)]
                if lengths is not None:
                    delta = [value - round(value / lengths[axis]) * lengths[axis]
                             for axis, value in enumerate(delta)]
                distance = math.sqrt(sum(value * value for value in delta))
                if distance > radius:
                    continue
                if minimum_all is None or distance < minimum_all:
                    minimum_all = distance
                direct = (roles[i], roles[j])
                reverse = (roles[j], roles[i])
                if direct in pair_set:
                    retain(distance, i, j, *direct)
                if reverse in pair_set and reverse != direct:
                    retain(distance, j, i, *reverse)

    if "minimumIntermolecularDistanceAngstrom" in requested and minimum_all is None:
        return _unavailable("No intermolecular atom pair was found inside the declared search radius.",
                            [f"{a}|{b}" for a, b in role_pairs])
    upper = [points[i][2] for i, item in enumerate(raw_mapping)
             if item.get("moleculeRole") == "lipid" and item.get("physicalSide") == "upper"
             and item.get("atomRole") == "head"]
    lower = [points[i][2] for i, item in enumerate(raw_mapping)
             if item.get("moleculeRole") == "lipid" and item.get("physicalSide") == "lower"
             and item.get("atomRole") == "head"]
    backbone = [points[i][2] for i, item in enumerate(raw_mapping)
                if item.get("moleculeRole") == "protein" and item.get("atomRole") == "backbone"]
    if (("upperLipidHeadMeanZAngstrom" in requested or "leafletHeadSeparationAngstrom" in requested or
         "proteinBilayerMidplaneOffsetAngstrom" in requested) and not upper or
        ("lowerLipidHeadMeanZAngstrom" in requested or "leafletHeadSeparationAngstrom" in requested or
         "proteinBilayerMidplaneOffsetAngstrom" in requested) and not lower or
        "proteinBilayerMidplaneOffsetAngstrom" in requested and not backbone):
        return _unavailable("A required protein-backbone or leaflet-head atom group is absent.",
                            [f"{a}|{b}" for a, b in role_pairs])
    if lengths is not None and backbone:
        # A stage can drift through a periodic boundary. Use the backbone's
        # circular z center as the image anchor, then place both leaflets in
        # that same local image before evaluating their relative arrangement.
        length_z = lengths[2]
        sine = sum(math.sin(2 * math.pi * z / length_z) for z in backbone)
        cosine = sum(math.cos(2 * math.pi * z / length_z) for z in backbone)
        if math.hypot(sine, cosine) / len(backbone) < 1e-6:
            return _unavailable("Protein z image is ambiguous across the periodic boundary.",
                                [f"{a}|{b}" for a, b in role_pairs])
        anchor = (math.atan2(sine, cosine) / (2 * math.pi) * length_z) % length_z

        def local_image(values: list[float]) -> list[float]:
            return [z - round((z - anchor) / length_z) * length_z for z in values]

        backbone, upper, lower = local_image(backbone), local_image(upper), local_image(lower)
    measurements = []
    for name in requested:
        if name == "minimumIntermolecularDistanceAngstrom":
            value, scope = minimum_all, "wholeSystem"
        elif name == "upperLipidHeadMeanZAngstrom":
            value, scope = sum(upper) / len(upper), "upperLeaflet"
        elif name == "lowerLipidHeadMeanZAngstrom":
            value, scope = sum(lower) / len(lower), "lowerLeaflet"
        elif name == "leafletHeadSeparationAngstrom":
            value = sum(upper) / len(upper) - sum(lower) / len(lower)
            scope = "bilayer"
        else:
            value = sum(backbone) / len(backbone) - (sum(upper) / len(upper) + sum(lower) / len(lower)) / 2
            scope = "proteinVsBilayer"
        measurements.append({"name": name, "value": value, "unit": "angstrom", "scope": scope})
    contacts = [{"firstAtomIndex": i, "secondAtomIndex": j,
                 "firstMoleculeRole": first, "secondMoleculeRole": second,
                 "distanceAngstrom": distance,
                 "radiusSumAngstrom": observed_radii[i] + observed_radii[j]}
                for negative_distance, i, j, first, second in sorted(observed_contacts, reverse=True)
                for distance in (-negative_distance,)]
    return {"standing": "Observed", "unavailableReason": None, "measurements": measurements,
            "locatedContacts": contacts,
            "coveredRolePairs": [f"{first}|{second}" for first, second in role_pairs],
            "rolePairMeasurements": [
                {"firstMoleculeRole": first, "secondMoleculeRole": second,
                 "pairsWithinSearchRadius": pair_measurements[(first, second)]["count"],
                 "minimumDistanceAngstrom": pair_measurements[(first, second)]["minimum"]}
                for first, second in role_pairs],
            "limitations": []}
