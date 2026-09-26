"""Read-only, independent inspection of saved Slice 4 construction and final State.

This is a test evidence tool. It neither calls the scientific worker nor decides
whether the product may qualify a completed stage. The expected file must bind
all seven saved artifacts by SHA-256 and supply the expected ordered atom count.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from openmm import XmlSerializer, unit
from openmm.app import PDBxFile


FILES = {
    "topologyCif": "constructed-topology.cif",
    "topologyJson": "constructed-topology.json",
    "systemXml": "constructed-system.xml",
    "stateXml": "constructed-state.xml",
    "correspondenceJson": "constructed-correspondence.json",
    "minimizedCif": "minimized-coordinates.cif",
    "minimizedStateXml": "minimized-state.xml",
}
CONTACT_CLASSES = (("water", "water"), ("water", "ion"),
                   ("ion", "ion"), ("protein", "lipid"))
VALID_ROLES = {"protein", "lipid", "water", "ion"}


class InspectionError(ValueError):
    pass


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise InspectionError(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _number(value, label: str) -> float:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool),
             f"{label} must be numeric")
    number = float(value)
    _require(math.isfinite(number), f"{label} is nonfinite")
    return number


def _load_paths(constructed_dir: Path, minimized_dir: Path, expected: dict):
    _require(isinstance(expected, dict), "Expected metadata must be a JSON object")
    count = expected.get("atomCount")
    _require(isinstance(count, int) and not isinstance(count, bool) and count > 0,
             "Expected positive atomCount is required")
    hashes = expected.get("sha256ByRole")
    _require(isinstance(hashes, dict) and set(hashes) == set(FILES),
             "Expected sha256ByRole must identify all seven artifact roles exactly")
    paths = {}
    actual_hashes = {}
    for role, filename in FILES.items():
        directory = minimized_dir if role.startswith("minimized") else constructed_dir
        path = directory / filename
        _require(path.is_file(), f"Missing {role} artifact: {path}")
        declared = hashes[role]
        _require(isinstance(declared, str) and len(declared) == 64 and
                 all(letter in "0123456789abcdefABCDEF" for letter in declared),
                 f"Invalid expected SHA-256 for {role}")
        actual = _sha256(path)
        _require(actual == declared.lower(), f"SHA-256 mismatch for {role}")
        paths[role] = path
        actual_hashes[role] = actual
    return count, paths, actual_hashes


def _identity(cif: PDBxFile, sidecar: dict, label: str) -> None:
    _require(sidecar.get("formatVersion") == 1, "Unknown topology sidecar format")
    chains, residues, atoms = (sidecar.get(name) for name in ("chains", "residues", "atoms"))
    _require(all(isinstance(group, list) for group in (chains, residues, atoms)),
             "Malformed ordered topology identity")
    observed = list(cif.topology.atoms())
    _require(len(observed) == len(atoms), f"{label} CIF atom count differs from topology sidecar")
    for index, (actual, declared) in enumerate(zip(observed, atoms)):
        _require(isinstance(declared, dict), f"Malformed sidecar atom {index}")
        residue_index = declared.get("residueIndex")
        _require(isinstance(residue_index, int) and 0 <= residue_index < len(residues),
                 f"Invalid sidecar residue index for atom {index}")
        residue = residues[residue_index]
        chain_index = residue.get("chainIndex") if isinstance(residue, dict) else None
        _require(isinstance(chain_index, int) and 0 <= chain_index < len(chains),
                 f"Invalid sidecar chain index for atom {index}")
        chain = chains[chain_index]
        # OpenMM's mmCIF writer/readback assigns sequential _atom_site.id
        # values to generated atoms. The sidecar keeps native atom IDs; those
        # technical serials are not the atom identity or correspondence key.
        _require(isinstance(declared.get("id"), str) and declared["id"] and
                 isinstance(actual.id, str) and actual.id,
                 f"{label} atom serial is missing at index {index}")
        expected_identity = (str(chain.get("id")), str(residue.get("id")),
                             str(residue.get("insertionCode") or ""),
                             str(residue.get("name")), str(declared.get("name")),
                             str(declared.get("element")))
        actual_identity = (actual.residue.chain.id, actual.residue.id,
                           actual.residue.insertionCode, actual.residue.name,
                           actual.name,
                           actual.element.symbol if actual.element else "")
        _require(actual_identity == expected_identity,
                 f"{label} CIF ordered identity differs at atom {index}")


def _bond_components(sidecar: dict, count: int) -> np.ndarray:
    bonds = sidecar.get("bonds")
    _require(isinstance(bonds, list), "Topology sidecar lacks bonds")
    parent = list(range(count))

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    seen = set()
    for ordinal, item in enumerate(bonds):
        pair = item.get("atomIndices") if isinstance(item, dict) else None
        _require(isinstance(pair, list) and len(pair) == 2 and
                 all(isinstance(i, int) and not isinstance(i, bool) and 0 <= i < count
                     for i in pair) and pair[0] != pair[1],
                 f"Invalid topology bond {ordinal}")
        ordered = tuple(sorted(pair))
        _require(ordered not in seen, f"Duplicate topology bond {ordinal}")
        seen.add(ordered)
        a, b = root(pair[0]), root(pair[1])
        parent[b] = a
    return np.asarray([root(i) for i in range(count)], dtype=np.int32)


def _correspondence(data: dict, sidecar: dict, count: int, expected: dict,
                    topology_hash: str):
    _require(isinstance(data, dict) and data.get("complete") is True,
             "Construction correspondence is incomplete")
    _require(data.get("resultId") == topology_hash,
             "Correspondence resultId differs from constructed CIF SHA-256")
    if "sourceId" in expected:
        _require(data.get("sourceId") == expected["sourceId"],
                 "Correspondence sourceId differs from expected input")
    mapping = data.get("atoms")
    _require(isinstance(mapping, list) and len(mapping) == count,
             "Correspondence does not cover every ordered atom")
    indices = [item.get("resultAtomIndex") if isinstance(item, dict) else None
               for item in mapping]
    _require(indices == list(range(count)), "Correspondence atom indices are not complete and ordered")
    ids = [item.get("resultAtomId") for item in mapping]
    _require(all(isinstance(value, str) and value for value in ids) and len(set(ids)) == count,
             "Correspondence result atom IDs are missing or repeated")
    roles = []
    for index, item in enumerate(mapping):
        role = item.get("moleculeRole")
        _require(role in VALID_ROLES, f"Unknown molecule role at atom {index}")
        _require(item.get("element") == sidecar["atoms"][index]["element"],
                 f"Correspondence element differs at atom {index}")
        roles.append(role)
    return mapping, np.asarray(roles, dtype="U8")


def _population(sidecar: dict, mapping: list[dict], molecules: np.ndarray) -> dict:
    """Count actual bonded molecules, using correspondence only for their roles."""
    groups = defaultdict(list)
    for index, root in enumerate(molecules):
        groups[int(root)].append(index)
    lipids = defaultdict(int)
    protein_atoms = water = sodium = chloride = 0
    expected_residue_name = {"DMPC": "DMP", "HOH": "HOH", "NA": "NA", "CL": "CL"}
    for indices in groups.values():
        identities = [mapping[index] for index in indices]
        roles = {item["moleculeRole"] for item in identities}
        _require(len(roles) == 1, "A bonded molecule mixes molecular roles")
        role = next(iter(roles))
        if role == "protein":
            _require(all(item.get("generatedSpeciesId") is None and
                         item.get("generatedComponentRole") is None for item in identities),
                     "A source protein component carries a generated identity")
            protein_atoms += len(indices)
            continue
        species = {item.get("generatedSpeciesId") for item in identities}
        components = {item.get("generatedComponentRole") for item in identities}
        residues = {sidecar["atoms"][index]["residueIndex"] for index in indices}
        _require(len(species) == len(components) == len(residues) == 1,
                 "A generated molecule has mixed species, component roles or residues")
        name = next(iter(species))
        component = next(iter(components))
        residue_name = sidecar["residues"][next(iter(residues))]["name"]
        _require(name in expected_residue_name and residue_name == expected_residue_name[name],
                 "Generated correspondence species differs from topology residue chemistry")
        if role == "lipid":
            sides = {item.get("physicalSide") for item in identities}
            _require(name == "DMPC" and component == "lipid" and len(sides) == 1 and
                     next(iter(sides)) in {"upper", "lower"},
                     "Lipid molecule lacks exact DMPC species or physical leaflet identity")
            heads = [index for index in indices if mapping[index].get("atomRole") == "head"]
            _require(len(heads) == 1 and sidecar["atoms"][heads[0]]["element"] == "P",
                     "DMPC molecule lacks one identified phosphate head")
            lipids[(name, next(iter(sides)))] += 1
        elif role == "water":
            _require(name == "HOH" and component == "water",
                     "Water molecule lacks exact HOH species identity")
            water += 1
        elif role == "ion":
            _require(len(indices) == 1 and
                     ((name == "NA" and component == "positiveIon" and
                       sidecar["atoms"][indices[0]]["element"] == "Na") or
                      (name == "CL" and component == "negativeIon" and
                       sidecar["atoms"][indices[0]]["element"] == "Cl")),
                     "Ion molecule lacks exact NA/CL identity or charge role")
            if name == "NA":
                sodium += 1
            else:
                chloride += 1
        else:
            raise InspectionError(f"Unsupported generated molecular role {role}")
    _require(protein_atoms > 0 and lipids[("DMPC", "upper")] > 0 and
             lipids[("DMPC", "lower")] > 0 and water > 0,
             "Required protein, both DMPC leaflets or water population is absent")
    return {"proteinAtomCount": protein_atoms,
            "lipidCounts": [{"speciesId": species, "physicalSide": side, "count": count}
                            for (species, side), count in sorted(lipids.items())],
            "waterCount": water, "sodiumCount": sodium, "chlorideCount": chloride,
            "bondedComponentCount": len(groups)}


def _state_xyz(state, count: int, label: str) -> np.ndarray:
    try:
        xyz = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom), dtype=float)
    except Exception as error:
        raise InspectionError(f"{label} State has no positions") from error
    _require(xyz.shape == (count, 3) and np.isfinite(xyz).all(),
             f"{label} State coordinates are absent or nonfinite")
    return xyz


def _box(state, sidecar: dict, label: str) -> np.ndarray:
    try:
        vectors = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom),
                             dtype=float)
    except Exception as error:
        raise InspectionError(f"{label} State has no periodic cell") from error
    _require(vectors.shape == (3, 3) and np.isfinite(vectors).all(),
             f"{label} periodic cell is invalid")
    _require(np.max(np.abs(vectors - np.diag(np.diag(vectors)))) <= 1e-6 and
             np.all(np.diag(vectors) > 0),
             f"{label} periodic cell is not positive orthorhombic")
    if label == "initial":
        declared = np.asarray(sidecar.get("boxVectorsAngstrom"), dtype=float)
        _require(declared.shape == (3, 3) and np.isfinite(declared).all() and
                 np.max(np.abs(declared - vectors)) <= 1e-6,
                 "Initial State cell differs from topology sidecar")
    return np.diag(vectors)


def _coordinates_match(cif: PDBxFile, xyz: np.ndarray, cell: np.ndarray,
                       label: str) -> float:
    exported = np.asarray(cif.positions.value_in_unit(unit.angstrom), dtype=float)
    _require(exported.shape == xyz.shape and np.isfinite(exported).all(),
             f"{label} CIF coordinates are absent or nonfinite")
    deviation = float(np.max(np.linalg.norm(exported - xyz, axis=1)))
    _require(deviation <= 0.0002, f"{label} CIF and State coordinates differ by {deviation:.6g} A")
    cif_cell = np.asarray(cif.topology.getPeriodicBoxVectors().value_in_unit(unit.angstrom), dtype=float)
    _require(cif_cell.shape == (3, 3) and
             np.max(np.abs(cif_cell - np.diag(cell))) <= 0.0011,
             f"{label} CIF and State cells differ")
    return deviation


def _force_projection(system, final_state, xyz_angstrom: np.ndarray) -> dict:
    count = len(xyz_angstrom)
    _require(system.getNumConstraints() > 0,
             "Final System lacks the expected constrained degrees of freedom")
    try:
        forces = np.asarray(final_state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer), dtype=float)
    except Exception as error:
        raise InspectionError("Final State has no saved physical forces") from error
    _require(forces.shape == (count, 3) and np.isfinite(forces).all(),
             "Final State forces are absent or nonfinite")
    xyz_nm = xyz_angstrom / 10.0
    parent = list(range(count))

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    edges = []
    seen = set()
    max_error = 0.0
    for ordinal in range(system.getNumConstraints()):
        i, j, target = system.getConstraintParameters(ordinal)
        i, j = int(i), int(j)
        wanted = float(target.value_in_unit(unit.nanometer))
        _require(0 <= i < count and 0 <= j < count and i != j and
                 math.isfinite(wanted) and wanted > 0, f"Invalid System constraint {ordinal}")
        pair = tuple(sorted((i, j)))
        _require(pair not in seen, f"Duplicate System constraint {ordinal}")
        seen.add(pair)
        delta = xyz_nm[i] - xyz_nm[j]
        observed = float(np.linalg.norm(delta))
        _require(observed > 0 and math.isfinite(observed),
                 f"Invalid final constrained distance {ordinal}")
        max_error = max(max_error, abs(observed - wanted) / wanted)
        edges.append((i, j, delta / observed))
        parent[root(j)] = root(i)
    groups = defaultdict(list)
    for edge in edges:
        groups[root(edge[0])].append(edge)
    tangent = forces.copy()
    max_orthogonality = 0.0
    for group_edges in groups.values():
        indices = sorted({atom for i, j, _ in group_edges for atom in (i, j)})
        _require(len(indices) <= 16 and len(group_edges) <= 16,
                 "Constraint component exceeds bounded HBonds inspection")
        lookup = {atom: ordinal for ordinal, atom in enumerate(indices)}
        jacobian = np.zeros((len(group_edges), 3 * len(indices)))
        for row, (i, j, direction) in enumerate(group_edges):
            jacobian[row, 3 * lookup[i]:3 * lookup[i] + 3] = direction
            jacobian[row, 3 * lookup[j]:3 * lookup[j] + 3] = -direction
        # An SVD gives an orthonormal basis of the constraint-normal space.
        # This independently checks the worker's least-squares projection.
        _, singular, right = np.linalg.svd(jacobian, full_matrices=False)
        rank = int(np.count_nonzero(singular > max(singular[0] * 1e-12, 1e-12)))
        _require(rank == len(group_edges), "Constraint Jacobian is rank deficient")
        normal_basis = right[:rank]
        local = forces[indices].reshape(-1)
        remainder = local - normal_basis.T @ (normal_basis @ local)
        max_orthogonality = max(max_orthogonality,
                                float(np.max(np.abs(jacobian @ remainder))))
        tangent[indices] = remainder.reshape(len(indices), 3)
    _require(max_orthogonality <= 1e-6 * max(1, float(np.max(np.abs(forces)))),
             "Constraint-tangent projection residual is unresolved")
    energy = float(final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    _require(math.isfinite(energy), "Final State energy is nonfinite")
    return {"method": "independent-svd-constraint-tangent-per-particle-v1",
            "constraintCount": len(edges), "constraintComponentCount": len(groups),
            "maximumRelativeConstraintError": max_error,
            "maximumTangentConstraintDotResidual": max_orthogonality,
            "rawComponentRmsKjMolNm": float(np.sqrt(np.mean(forces * forces))),
            "tangentPerParticleRmsKjMolNm": float(np.sqrt(np.sum(tangent * tangent) / count)),
            "finalPotentialEnergyKjMol": energy}


def _relative_heads(xyz: np.ndarray, box: np.ndarray, mapping: list[dict]) -> dict:
    upper = [index for index, item in enumerate(mapping)
             if item["moleculeRole"] == "lipid" and item.get("atomRole") == "head"
             and item.get("physicalSide") == "upper"]
    lower = [index for index, item in enumerate(mapping)
             if item["moleculeRole"] == "lipid" and item.get("atomRole") == "head"
             and item.get("physicalSide") == "lower"]
    backbone = [index for index, item in enumerate(mapping)
                if item["moleculeRole"] == "protein" and item.get("atomRole") == "backbone"]
    _require(bool(upper and lower and backbone),
             "Required protein backbone or physical leaflet head group is absent")
    length = box[2]
    phase = np.exp(2j * np.pi * xyz[backbone, 2] / length).mean()
    _require(abs(phase) >= 1e-6, "Protein periodic z image is ambiguous")
    anchor = np.angle(phase) / (2 * np.pi) * length

    def local_mean(indices):
        values = xyz[indices, 2]
        return float(np.mean(values - np.rint((values - anchor) / length) * length))

    upper_z, lower_z, protein_z = (local_mean(group) for group in (upper, lower, backbone))
    return {"upperHeadMeanZAngstrom": upper_z, "lowerHeadMeanZAngstrom": lower_z,
            "leafletHeadSeparationAngstrom": upper_z - lower_z,
            "proteinBilayerMidplaneOffsetAngstrom": protein_z - (upper_z + lower_z) / 2}


def _contacts(xyz: np.ndarray, box: np.ndarray, roles: np.ndarray,
              elements: np.ndarray, molecules: np.ndarray, radius: float) -> dict:
    _require(np.all(box > 2 * radius), "Cell is too small for complete periodic contact search")
    counts = np.floor(box / radius).astype(int)
    wrapped = xyz % box
    cell_indices = np.floor(wrapped / box * counts).astype(int)
    bins = defaultdict(list)
    for atom, cell in enumerate(cell_indices):
        bins[tuple(int(value) for value in cell)].append(atom)
    bins = {key: np.asarray(value, dtype=np.int32) for key, value in bins.items()}
    heavy_report = {"|".join(pair): {"heavyPairsWithinRadius": 0,
                                     "heavyPairsBelow1_5Angstrom": 0,
                                     "heavyPairsBelow2_2Angstrom": 0,
                                     "periodicImageHeavyPairsWithinRadius": 0,
                                     "periodicImageHeavyPairsBelow1_5Angstrom": 0,
                                     "periodicImageHeavyPairsBelow2_2Angstrom": 0,
                                     "nearestHeavyAngstrom": None}
                    for pair in CONTACT_CLASSES}
    all_report = {name: {"allAtomPairsWithinRadius": 0,
                         "allAtomPairsBelow1_5Angstrom": 0,
                         "allAtomPairsBelow2_2Angstrom": 0,
                         "periodicImageAllAtomPairsWithinRadius": 0,
                         "periodicImageAllAtomPairsBelow1_5Angstrom": 0,
                         "periodicImageAllAtomPairsBelow2_2Angstrom": 0,
                         "nearestAllAtomAngstrom": None}
                  for name in ("wholeSystem", *("|".join(pair) for pair in CONTACT_CLASSES))}

    def add_all(item: dict, distances: np.ndarray, across: np.ndarray) -> None:
        item["allAtomPairsWithinRadius"] += int(len(distances))
        item["allAtomPairsBelow1_5Angstrom"] += int(np.count_nonzero(distances < 1.5))
        item["allAtomPairsBelow2_2Angstrom"] += int(np.count_nonzero(distances < 2.2))
        item["periodicImageAllAtomPairsWithinRadius"] += int(np.count_nonzero(across))
        item["periodicImageAllAtomPairsBelow1_5Angstrom"] += int(np.count_nonzero(
            across & (distances < 1.5)))
        item["periodicImageAllAtomPairsBelow2_2Angstrom"] += int(np.count_nonzero(
            across & (distances < 2.2)))
        nearest = float(np.min(distances))
        if item["nearestAllAtomAngstrom"] is None or nearest < item["nearestAllAtomAngstrom"]:
            item["nearestAllAtomAngstrom"] = nearest

    def add_heavy(item: dict, distances: np.ndarray, across: np.ndarray) -> None:
        item["heavyPairsWithinRadius"] += int(len(distances))
        item["heavyPairsBelow1_5Angstrom"] += int(np.count_nonzero(distances < 1.5))
        item["heavyPairsBelow2_2Angstrom"] += int(np.count_nonzero(distances < 2.2))
        item["periodicImageHeavyPairsWithinRadius"] += int(np.count_nonzero(across))
        item["periodicImageHeavyPairsBelow1_5Angstrom"] += int(np.count_nonzero(
            across & (distances < 1.5)))
        item["periodicImageHeavyPairsBelow2_2Angstrom"] += int(np.count_nonzero(
            across & (distances < 2.2)))
        nearest = float(np.min(distances))
        if item["nearestHeavyAngstrom"] is None or nearest < item["nearestHeavyAngstrom"]:
            item["nearestHeavyAngstrom"] = nearest

    radius_sq = radius * radius
    tile = 128  # At most 128 x 128 pairs and three coordinates per distance block.
    for key, first_bin in bins.items():
        neighbors = {tuple((key[axis] + shift[axis]) % counts[axis] for axis in range(3))
                     for shift in product((-1, 0, 1), repeat=3)}
        for other in neighbors:
            if other < key or other not in bins:
                continue
            second_bin = bins[other]
            for first_start in range(0, len(first_bin), tile):
                first = first_bin[first_start:first_start + tile]
                for second_start in range(0, len(second_bin), tile):
                    second = second_bin[second_start:second_start + tile]
                    raw = wrapped[first, None, :] - wrapped[None, second, :]
                    images = np.rint(raw / box)
                    delta = raw - images * box
                    squared = np.einsum("ijk,ijk->ij", delta, delta)
                    inside = squared <= radius_sq
                    if other == key:
                        inside &= first[:, None] < second[None, :]
                    row, column = np.nonzero(inside)
                    if len(row) == 0:
                        continue
                    a, b = first[row], second[column]
                    distinct = molecules[a] != molecules[b]
                    a, b, row, column = a[distinct], b[distinct], row[distinct], column[distinct]
                    if len(a) == 0:
                        continue
                    distances = np.sqrt(squared[row, column])
                    # Distance is wrap invariant. The image flag describes
                    # the saved State coordinate representation, whose raw
                    # coordinates may lie on opposite periodic sides.
                    across = np.any(np.rint((xyz[a] - xyz[b]) / box) != 0, axis=1)
                    add_all(all_report["wholeSystem"], distances, across)
                    heavy = (elements[a] != "H") & (elements[b] != "H")
                    for left, right in CONTACT_CLASSES:
                        selected = (roles[a] == left) & (roles[b] == right)
                        if left != right:
                            selected |= (roles[a] == right) & (roles[b] == left)
                        if not np.any(selected):
                            continue
                        name = f"{left}|{right}"
                        add_all(all_report[name], distances[selected], across[selected])
                        selected_heavy = selected & heavy
                        if np.any(selected_heavy):
                            add_heavy(heavy_report[name], distances[selected_heavy],
                                      across[selected_heavy])
    return {"heavy": heavy_report, "allAtom": all_report}


def inspect(constructed_dir: Path, minimized_dir: Path, expected: dict,
            *, contact_radius_angstrom: float = 6.0) -> dict:
    radius = _number(contact_radius_angstrom, "contact radius")
    _require(radius > 0, "Contact radius must be positive")
    count, paths, hashes = _load_paths(constructed_dir, minimized_dir, expected)
    try:
        sidecar = json.loads(paths["topologyJson"].read_text(encoding="utf-8"))
        correspondence = json.loads(paths["correspondenceJson"].read_text(encoding="utf-8"))
        initial_cif = PDBxFile(str(paths["topologyCif"]))
        final_cif = PDBxFile(str(paths["minimizedCif"]))
        system = XmlSerializer.deserialize(paths["systemXml"].read_text(encoding="utf-8"))
        initial = XmlSerializer.deserialize(paths["stateXml"].read_text(encoding="utf-8"))
        final = XmlSerializer.deserialize(paths["minimizedStateXml"].read_text(encoding="utf-8"))
    except Exception as error:
        raise InspectionError(f"Saved artifact is malformed: {error}") from error
    _require(isinstance(sidecar, dict), "Topology sidecar must be an object")
    _identity(initial_cif, sidecar, "Initial")
    _identity(final_cif, sidecar, "Final")
    _require(len(sidecar["atoms"]) == count and system.getNumParticles() == count,
             "Expected, sidecar and System particle counts differ")
    mapping, roles = _correspondence(correspondence, sidecar, count, expected,
                                      hashes["topologyCif"])
    elements = np.asarray([item["element"] for item in sidecar["atoms"]], dtype="U2")
    molecules = _bond_components(sidecar, count)
    population = _population(sidecar, mapping, molecules)
    initial_xyz = _state_xyz(initial, count, "Initial")
    final_xyz = _state_xyz(final, count, "Final")
    initial_box = _box(initial, sidecar, "initial")
    final_box = _box(final, sidecar, "final")
    default_box = np.asarray([vector.value_in_unit(unit.angstrom)
                              for vector in system.getDefaultPeriodicBoxVectors()], dtype=float)
    _require(default_box.shape == (3, 3) and
             np.max(np.abs(default_box - np.diag(initial_box))) <= 1e-6 and
             np.max(np.abs(final_box - initial_box)) <= 1e-6,
             "Saved System, initial State and final State periodic cells differ")
    initial_energy = float(initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    _require(math.isfinite(initial_energy), "Initial State energy is nonfinite")
    initial_roundtrip = _coordinates_match(initial_cif, initial_xyz, initial_box, "Initial")
    final_roundtrip = _coordinates_match(final_cif, final_xyz, final_box, "Final")
    initial_contacts = _contacts(initial_xyz, initial_box, roles, elements, molecules, radius)
    final_contacts = _contacts(final_xyz, final_box, roles, elements, molecules, radius)
    force = _force_projection(system, final, final_xyz)
    stages = {
        "initial": {"cellAngstrom": initial_box.tolist(),
                    "contactRadiusAngstrom": radius,
                    "heavyPeriodicContacts": initial_contacts["heavy"],
                    "allAtomPeriodicContacts": initial_contacts["allAtom"],
                    "relativeOrganization": _relative_heads(initial_xyz, initial_box, mapping)},
        "final": {"cellAngstrom": final_box.tolist(),
                  "contactRadiusAngstrom": radius,
                  "heavyPeriodicContacts": final_contacts["heavy"],
                  "allAtomPeriodicContacts": final_contacts["allAtom"],
                  "relativeOrganization": _relative_heads(final_xyz, final_box, mapping)},
    }
    comparison = {name: {
        "changeInHeavyPairsBelow2_2Angstrom":
            final_contacts["heavy"][name]["heavyPairsBelow2_2Angstrom"] -
            initial_contacts["heavy"][name]["heavyPairsBelow2_2Angstrom"],
        "changeInPeriodicImageHeavyPairsBelow2_2Angstrom":
            final_contacts["heavy"][name]["periodicImageHeavyPairsBelow2_2Angstrom"] -
            initial_contacts["heavy"][name]["periodicImageHeavyPairsBelow2_2Angstrom"]}
        for name in initial_contacts["heavy"]}
    all_comparison = {name: {
        "changeInAllAtomPairsBelow1_5Angstrom":
            final_contacts["allAtom"][name]["allAtomPairsBelow1_5Angstrom"] -
            initial_contacts["allAtom"][name]["allAtomPairsBelow1_5Angstrom"],
        "changeInPeriodicImageAllAtomPairsBelow1_5Angstrom":
            final_contacts["allAtom"][name]["periodicImageAllAtomPairsBelow1_5Angstrom"] -
            initial_contacts["allAtom"][name]["periodicImageAllAtomPairsBelow1_5Angstrom"]}
        for name in initial_contacts["allAtom"]}
    return {"status": "observed", "inspectorVersion": 1,
            "artifactSha256": hashes, "atomCount": count,
            "population": population,
            "topologyBondCount": len(sidecar["bonds"]),
            "maximumCifStateDeviationAngstrom": {"initial": initial_roundtrip,
                                                  "final": final_roundtrip},
            "finalForce": force, "stages": stages, "contactChange": comparison,
            "allAtomContactChange": all_comparison}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constructed-dir", required=True, type=Path)
    parser.add_argument("--minimized-dir", type=Path,
                        help="Defaults to the constructed artifact directory")
    parser.add_argument("--expected", required=True, type=Path,
                        help="JSON with atomCount and sha256ByRole for all seven artifacts")
    parser.add_argument("--output", type=Path, help="Write report JSON here; otherwise print stdout")
    parser.add_argument("--contact-radius-angstrom", type=float, default=6.0)
    args = parser.parse_args(argv)
    try:
        expected = json.loads(args.expected.read_text(encoding="utf-8"))
        result = inspect(args.constructed_dir, args.minimized_dir or args.constructed_dir,
                         expected, contact_radius_angstrom=args.contact_radius_angstrom)
        rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.output:
            args.output.write_text(rendered, encoding="utf-8")
        else:
            sys.stdout.write(rendered)
        return 0
    except Exception as error:
        sys.stderr.write(json.dumps({"status": "failed", "reason": str(error)}) + "\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
