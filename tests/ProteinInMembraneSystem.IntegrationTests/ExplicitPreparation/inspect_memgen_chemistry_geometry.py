"""Independent chemistry and periodic-geometry read-back of a saved Memgen trial.

Usage: out/python/bin/python inspect_memgen_chemistry_geometry.py CASE_DIRECTORY

This test-side oracle reads the pinned request, Amber topology/restart and
source-versioned case manifest. It never calls the product's mapping, stereo or
contact routines. Its 1.5 A floor applies only to the conditioned handoff.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree


class InspectionError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise InspectionError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def components(adjacency: list[set[int]]) -> list[list[int]]:
    seen = set()
    answer = []
    for start in range(len(adjacency)):
        if start in seen:
            continue
        group, queue = [], [start]
        seen.add(start)
        while queue:
            node = queue.pop()
            group.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        answer.append(sorted(group))
    return answer


def topology_graph(topology) -> tuple[list, list[set[int]]]:
    atoms = list(topology.atoms())
    adjacency = [set() for _ in atoms]
    for left, right in topology.bonds():
        adjacency[left.index].add(right.index)
        adjacency[right.index].add(left.index)
    return atoms, adjacency


def local_graph(atoms: list, adjacency: list[set[int]], indices: list[int]):
    local = {global_index: index for index, global_index in enumerate(indices)}
    elements = [atoms[index].element.symbol for index in indices]
    names = [atoms[index].name for index in indices]
    edges = [{local[neighbor] for neighbor in adjacency[index] if neighbor in local}
             for index in indices]
    return elements, names, edges


def graph_colors(first_elements, first_edges, second_elements, second_edges):
    def initial(elements, edges):
        return [(elements[i], len(edges[i]), tuple(sorted(elements[j] for j in edges[i])))
                for i in range(len(elements))]
    first, second = initial(first_elements, first_edges), initial(second_elements, second_edges)
    for _ in range(max(len(first), 1)):
        palette = {key: number for number, key in enumerate(sorted(set(first + second)))}
        first_colors = [palette[key] for key in first]
        second_colors = [palette[key] for key in second]
        require(Counter(first_colors) == Counter(second_colors),
                "The final lipid graph differs from its selected template")
        next_first = [(first_colors[i], tuple(sorted(first_colors[j] for j in first_edges[i])))
                      for i in range(len(first))]
        next_second = [(second_colors[i], tuple(sorted(second_colors[j] for j in second_edges[i])))
                       for i in range(len(second))]
        if next_first == first and next_second == second:
            return first_colors, second_colors
        first, second = next_first, next_second
    return first_colors, second_colors


def graph_map(reference, actual) -> tuple[dict[int, int], set[int]]:
    """Map by independent graph refinement; reject ambiguous stereo anchors."""
    ref_elements, ref_names, ref_edges = reference
    act_elements, act_names, act_edges = actual
    require(len(ref_elements) == len(act_elements), "Lipid atom count changed")
    ref_colors, act_colors = graph_colors(ref_elements, ref_edges, act_elements, act_edges)
    ref_by_name = Counter(ref_names)
    act_by_name = Counter(act_names)
    ref_groups = defaultdict(list)
    act_groups = defaultdict(list)
    for i, color in enumerate(ref_colors):
        ref_groups[color].append(i)
    for i, color in enumerate(act_colors):
        act_groups[color].append(i)
    require({key: len(value) for key, value in ref_groups.items()} ==
            {key: len(value) for key, value in act_groups.items()},
            "Lipid graph color inventory changed")
    # Exact unique shared names anchor Lipid21's retained head/sterol identity;
    # same-spelled names at different chemical sites are not accepted as anchors.
    anchors = {}
    for i, name in enumerate(ref_names):
        if ref_by_name[name] == act_by_name[name] == 1:
            j = act_names.index(name)
            if ref_colors[i] == act_colors[j] and ref_elements[i] == act_elements[j]:
                anchors[i] = j
    require(len(set(anchors.values())) == len(anchors), "Lipid graph anchors collide")
    mapping = dict(anchors)
    used = set(mapping.values())
    order = sorted((i for i in range(len(ref_elements)) if i not in mapping),
                   key=lambda i: (len(act_groups[ref_colors[i]]), -len(ref_edges[i]), i))

    def match(position):
        if position == len(order):
            return True
        i = order[position]
        for j in act_groups[ref_colors[i]]:
            if j in used or any((other in ref_edges[i]) != (mapped in act_edges[j])
                                for other, mapped in mapping.items()):
                continue
            mapping[i] = j
            used.add(j)
            if match(position + 1):
                return True
            used.remove(j)
            del mapping[i]
        return False

    require(match(0), "The final lipid bond graph is not the selected template graph")
    require(all({mapping[j] for j in ref_edges[i]} == act_edges[mapping[i]]
                for i in range(len(ref_elements))), "Lipid bond graph mapping is incomplete")
    # A topologically symmetric site without an exact name anchor cannot prove
    # which stereochemical descriptor atom it represents.
    unique = {i for i, color in enumerate(ref_colors) if len(ref_groups[color]) == 1}
    return mapping, unique | set(anchors)


def unwrap(points: np.ndarray, edges: list[set[int]], lengths: np.ndarray) -> np.ndarray:
    result = np.empty_like(points, dtype=float)
    seen = {0}
    result[0] = points[0]
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in edges[i]:
            step = points[j] - points[i]
            step -= np.rint(step / lengths) * lengths
            expected = result[i] + step
            if j not in seen:
                result[j] = expected
                seen.add(j)
                queue.append(j)
            else:
                require(np.linalg.norm(result[j] - expected) < 0.05,
                        "The lipid cannot be unwrapped consistently through its bonds")
    require(len(seen) == len(points), "The lipid bond graph is disconnected")
    return result


def stereo_value(kind: str, points: np.ndarray) -> float:
    require(points.shape == (4, 3), "A stereo descriptor needs four mapped atoms")
    if kind == "tetrahedral":
        return float(np.linalg.det(np.stack([points[i] - points[3] for i in range(3)])))
    if kind == "alkene":
        axis = points[2] - points[1]
        norm = float(axis @ axis)
        require(norm > 1e-8, "The alkene bond is degenerate")
        left = points[0] - points[1]
        right = points[3] - points[2]
        left -= axis * ((left @ axis) / norm)
        right -= axis * ((right @ axis) / norm)
        denominator = np.linalg.norm(left) * np.linalg.norm(right)
        require(denominator > 1e-8, "The alkene torsion is degenerate")
        return float((left @ right) / denominator)
    raise InspectionError(f"Unknown stereo descriptor {kind}")


def check_stereo(species: str, checks: list[dict], reference, actual,
                 actual_points: np.ndarray, lengths: np.ndarray) -> int:
    mapping, identified = graph_map(reference, actual)
    ref_elements, ref_names, _ = reference
    _, _, actual_edges = actual
    points = unwrap(actual_points, actual_edges, lengths)
    checked = 0
    for descriptor in checks:
        names = descriptor["atomNames"]
        require(len(names) == 4 and len(set(names)) == 4 and all(
            ref_names.count(name) == 1 for name in names),
            f"{species} has an ambiguous selected stereo descriptor")
        indices = [ref_names.index(name) for name in names]
        require(all(index in identified for index in indices),
                f"{species} stereo atom has an unresolved graph automorphism")
        if descriptor["kind"] == "alkene":
            require(all(indices[i + 1] in reference[2][indices[i]] for i in range(3)),
                    f"{species} selected alkene is not a bonded quartet")
        value = stereo_value(descriptor["kind"], points[[mapping[i] for i in indices]])
        expected = descriptor["expected"]
        valid = ((expected == "positive" and value > 1.0) or
                 (expected == "negative" and value < -1.0) or
                 (expected == "cis" and value > 0.5))
        require(valid, f"{species} {descriptor['kind']} stereochemistry changed: {value:.6g}")
        checked += 1
    return checked


def periodic_contacts(points: np.ndarray, lengths: np.ndarray, elements: list[str],
                      roles: list[str], molecule_ids: list[int], floor: float = 1.5) -> dict:
    """Independent 27-image heavy-contact screen on the conditioned handoff."""
    require(points.shape == (len(elements), 3) and len(roles) == len(elements) == len(molecule_ids),
            "Periodic contact arrays have different atom counts")
    require(np.isfinite(points).all() and np.isfinite(lengths).all() and
            np.all(lengths > 2 * floor), "Periodic cell or coordinates are invalid")
    heavy = np.asarray([i for i, element in enumerate(elements) if element.upper() not in {"H", "D"}],
                       dtype=int)
    wrapped = points[heavy] % lengths
    pairs = cKDTree(wrapped, boxsize=lengths).query_pairs(floor, output_type="ndarray")
    by_class = Counter()
    violations = []
    for first_local, second_local in pairs:
        first, second = int(heavy[first_local]), int(heavy[second_local])
        if molecule_ids[first] == molecule_ids[second]:
            continue  # Same physical component; self-images are checked after graph unwrapping.
        delta = wrapped[first_local] - wrapped[second_local]
        image = np.rint(delta / lengths).astype(int)
        distance = float(np.linalg.norm(delta - image * lengths))
        if distance >= floor:
            continue
        role_pair = tuple(sorted((roles[first], roles[second])))
        image_kind = ("interior", "face", "edge", "corner")[np.count_nonzero(image)]
        by_class[(role_pair, image_kind)] += 1
        violations.append({"atoms": [first, second], "roles": role_pair,
                           "imageKind": image_kind, "distanceAngstrom": distance})
    return {"violations": violations,
            "counts": {"|".join((*pair, image)): count for (pair, image), count in by_class.items()},
            "heavyAtomCount": len(heavy)}


def component_image_contacts(points: np.ndarray, lengths: np.ndarray,
                             adjacency: list[set[int]], elements: list[str],
                             roles: list[str], groups: list[list[int]], floor: float = 1.5) -> list[dict]:
    """Search copies of each unwrapped component, all 26 shifts."""
    violations = []
    for group in groups:
        local = {global_index: i for i, global_index in enumerate(group)}
        local_edges = [{local[j] for j in adjacency[index] if j in local} for index in group]
        local_points = unwrap(points[group], local_edges, lengths)
        heavy = np.asarray([i for i, index in enumerate(group)
                            if elements[index].upper() not in {"H", "D"}], dtype=int)
        if len(heavy) == 0:
            continue
        xyz = local_points[heavy]
        # A component narrower than L-floor in every axis cannot contact an image.
        if np.all(lengths - np.ptp(xyz, axis=0) >= floor):
            continue
        tree = cKDTree(xyz)
        for nx in (-1, 0, 1):
            for ny in (-1, 0, 1):
                for nz in (-1, 0, 1):
                    if (nx, ny, nz) == (0, 0, 0):
                        continue
                    shift = np.asarray([nx, ny, nz]) * lengths
                    hits = tree.query_ball_tree(cKDTree(xyz + shift), floor)
                    for i, neighbors in enumerate(hits):
                        for j in neighbors:
                            distance = float(np.linalg.norm(xyz[i] - xyz[j] - shift))
                            if distance < floor:
                                violations.append({"atoms": [group[heavy[i]], group[heavy[j]]],
                                    "role": roles[group[heavy[i]]],
                                    "imageKind": ("face", "edge", "corner")[
                                        sum(value != 0 for value in (nx, ny, nz)) - 1],
                                    "distanceAngstrom": distance})
    return violations


def inspect(case: Path) -> dict:
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, PDBxFile

    case = case.resolve()
    request_path = case / "request.json"
    request = json.loads(request_path.read_text())
    payload = request["payload"]
    manifest = json.loads((case.parent / "manifest.json").read_text())
    selected = [item for item in manifest["cases"]
                if Path(item["requestPath"]).resolve() == request_path]
    require(len(selected) == 1 and sha256(request_path) == selected[0]["requestSha256"],
            "The request is not bound by exactly one external fixture plan")
    terminals = [json.loads(line) for line in (case / "events.jsonl").read_text().splitlines()
                 if json.loads(line).get("kind") in {"result", "error"}]
    require(len(terminals) == 1 and terminals[0]["kind"] == "result",
            "The Memgen trial has no single observed result")
    result = terminals[0]["payload"]
    require(result.get("standing") == "observed" and
            result.get("requestId") == request["requestId"] and
            result.get("attemptId") == payload["attemptId"] and
            result.get("studyRevisionId") == payload["studyRevisionId"],
            "The observed Memgen trial is not the requested attempt and revision")
    artifacts = {item["role"]: item for item in result["artifacts"]}
    require(len(artifacts) == len(result["artifacts"]), "Provider artifact roles are duplicated")
    for role in ("amberTopology", "amberFinalRestart", "correspondenceJson", "providerPacked"):
        require(role in artifacts and sha256(Path(artifacts[role]["path"])) == artifacts[role]["sha256"],
                f"The {role} bytes changed")
    representations = payload["selectedSpeciesRepresentations"]
    references = {}
    for item in representations:
        source = Path(item["coordinateTemplatePath"])
        require(sha256(source) == item["coordinateTemplateSha256"],
                f"Selected {item['speciesId']} coordinate template changed")
        reference_atoms, reference_edges = topology_graph(PDBxFile(str(source)).topology)
        references[item["speciesId"]] = (item, local_graph(
            reference_atoms, reference_edges, list(range(len(reference_atoms)))))
    restart = AmberInpcrdFile(artifacts["amberFinalRestart"]["path"])
    topology = AmberPrmtopFile(artifacts["amberTopology"]["path"],
                               periodicBoxVectors=restart.boxVectors).topology
    atoms, adjacency = topology_graph(topology)
    points = np.asarray(restart.positions.value_in_unit(unit.angstrom), dtype=float)
    box = np.asarray([[float(value.value_in_unit(unit.angstrom)) for value in vector]
                      for vector in restart.boxVectors], dtype=float)
    require(points.shape == (len(atoms), 3) and box.shape == (3, 3) and
            np.max(np.abs(box - np.diag(np.diag(box)))) < 1e-6,
            "This independent periodic oracle requires the selected orthorhombic Memgen cell")
    lengths = np.diag(box)
    correspondence = json.loads(Path(artifacts["correspondenceJson"]["path"]).read_text())
    mapped = correspondence["atoms"]
    require(len(mapped) == len(atoms) and all(item["resultAtomIndex"] == i and
            item["element"] == atoms[i].element.symbol for i, item in enumerate(mapped)),
            "Final atom order or correspondence element changed")
    groups = components(adjacency)
    molecule_ids = [0] * len(atoms)
    for number, group in enumerate(groups):
        for index in group:
            molecule_ids[index] = number
    roles = [item["moleculeRole"] for item in mapped]
    elements = [atom.element.symbol for atom in atoms]
    observations = {"case": str(case),
                    "scope": "Independent selected-lipid graph/stereo and conditioned Amber periodic heavy-contact oracle",
                    "requestSha256": sha256(request_path),
                    "amberTopologySha256": artifacts["amberTopology"]["sha256"],
                    "amberFinalRestartSha256": artifacts["amberFinalRestart"]["sha256"],
                    "stereoChecksBySpecies": Counter(), "lipidMoleculesBySpecies": Counter(),
                    "periodicContactFloorAngstrom": 1.5,
                    "unavailableCoverage": [
                        "This focused oracle does not inspect retained source chemistry, provider controls or conditioning logs, salt-volume derivation, Amber-to-OpenMM parameter terms, final minimization, Host publication or export; separate inspectors cover those boundaries."]}
    counts_by_side = Counter()
    for group in groups:
        actual = local_graph(atoms, adjacency, group)
        formula = Counter(actual[0])
        group_roles = {roles[index] for index in group}
        candidates = [(species, item, reference) for species, (item, reference) in references.items()
                      if len(group) == len(reference[0]) and formula == Counter(reference[0])]
        if not candidates and "lipid" not in group_roles:
            continue
        require(len(candidates) == 1, "Selected lipid formulas do not uniquely identify a component")
        require(group_roles == {"lipid"}, "A lipid component has the wrong molecule role")
        declared_species = {mapped[index].get("generatedSpeciesId") for index in group}
        declared_sides = {mapped[index].get("physicalSide") for index in group}
        require(len(declared_species) == len(declared_sides) == 1 and
                next(iter(declared_sides)) in {"upper", "lower"},
                "A lipid has ambiguous generated species or physical side")
        species, item, reference = candidates[0]
        require(declared_species == {species}, "Lipid chemistry contradicts its generated species")
        observed = check_stereo(species, item["stereoChecks"], reference, actual,
                                points[group], lengths)
        observations["stereoChecksBySpecies"][species] += observed
        observations["lipidMoleculesBySpecies"][species] += 1
        counts_by_side[(species, next(iter(declared_sides)))] += 1
    require(all(observations["lipidMoleculesBySpecies"][species] > 0 for species in references),
            "A selected lipid was absent from the final Amber result")
    declared_counts = {(item["speciesId"], item["physicalSide"]): item["count"]
                       for item in result["observations"]["speciesCounts"]
                       if item["role"] == "lipid"}
    require(counts_by_side == declared_counts,
            "Actual lipid component counts do not match the provider's reported leaflet counts")
    ordinary = periodic_contacts(points, lengths, elements, roles, molecule_ids)
    images = component_image_contacts(points, lengths, adjacency, elements, roles, groups)
    observations["intermolecularHeavyContactScreen"] = ordinary
    observations["componentImageHeavyContactScreen"] = images
    observations["lipidMoleculesBySide"] = {"|".join(key): value for key, value in counts_by_side.items()}
    require(not ordinary["violations"] and not images,
            "The conditioned handoff has a severe periodic heavy-atom contact")
    observations["stereoChecksBySpecies"] = dict(observations["stereoChecksBySpecies"])
    observations["lipidMoleculesBySpecies"] = dict(observations["lipidMoleculesBySpecies"])
    observations["standing"] = "checksPassed"
    return observations


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: inspect_memgen_chemistry_geometry.py CASE_DIRECTORY")
    try:
        report = inspect(Path(sys.argv[1]))
    except Exception as error:
        report = {"standing": "failed", "case": sys.argv[1],
                  "failure": f"{type(error).__name__}: {error}"}
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["standing"] == "checksPassed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
