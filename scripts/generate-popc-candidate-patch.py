#!/usr/bin/env python3
"""Regenerate and verify the bounded OpenMM 8.6 POPC patch candidate.

This only removes a confirmed inverted lipid and its closest opposite-leaflet
partner.  It does not qualify the resulting patch for membrane construction.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys

import openmm
import openmm.app
from openmm import unit
from openmm.app import Modeller, PDBFile, PDBxFile


ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = ROOT / "config/policies/protein-membrane-slice4.json"
OUTPUT = ROOT / "config/policies/membrane-templates/POPC-OpenMM86-63x63.pdb"
SOURCE_SHA256 = "a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29"
DERIVED_SHA256 = "92a19c3470605d44ca779772f306ad693afdd11f1f78176655a3f1d2ee3a77b4"
REMOVED_IDS = frozenset({"62", "109"})


def require(condition: bool, explanation: str) -> None:
    if not condition:
        raise ValueError(explanation)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalized_name(name: str) -> str:
    # OpenMM's PDB reader rotates the final digit of five-character acyl names.
    match = re.fullmatch(r"([0-9])C(21|31)", name)
    return f"C{match.group(2)}{match.group(1)}" if match else name


def bonds(topology, residue) -> set[tuple[int, int]]:
    local = {atom.index: index for index, atom in enumerate(residue.atoms())}
    return {
        tuple(sorted((local[first.index], local[second.index])))
        for first, second in topology.bonds()
        if first.residue is residue and second.residue is residue
    }


def signed_volume(a, b, c, d) -> float:
    left = [a[i] - d[i] for i in range(3)]
    middle = [b[i] - d[i] for i in range(3)]
    right = [c[i] - d[i] for i in range(3)]
    return sum(left[i] * (middle[(i + 1) % 3] * right[(i + 2) % 3] -
                          middle[(i + 2) % 3] * right[(i + 1) % 3]) for i in range(3))


def lipid_observation(pdb, residue, reference, expected_bonds, stereo_checks, center):
    actual = list(residue.atoms())
    expected = list(reference.topology.atoms())
    require(len(actual) == len(expected) == 134, f"POP {residue.id}: atom count")
    require(all(normalized_name(a.name) == e.name and a.element == e.element
                for a, e in zip(actual, expected)), f"POP {residue.id}: atom identity/order")
    actual_bonds = bonds(pdb.topology, residue)
    require(actual_bonds == expected_bonds, f"POP {residue.id}: molecular bonds")
    require(len(stereo_checks) == 2 and stereo_checks[0]["kind"] == "tetrahedral" and
            stereo_checks[0]["atomNames"] == ["O21", "C1", "C3", "HS"] and
            stereo_checks[0]["expected"] == "negative" and
            stereo_checks[1]["kind"] == "alkene" and
            stereo_checks[1]["atomNames"] == ["C28", "C29", "C210", "C211"] and
            stereo_checks[1]["expected"] == "cis", "POPC catalogue stereo descriptors changed")

    positions = pdb.positions.value_in_unit(unit.angstrom)
    points = [positions[atom.index] for atom in actual]
    names = {normalized_name(atom.name): index for index, atom in enumerate(actual)}
    chiral = signed_volume(*(points[names[name]] for name in ("O21", "C1", "C3", "HS")))
    require(math.isfinite(chiral) and abs(chiral) > 1.0,
            f"POP {residue.id}: unresolved glycerol stereocenter")

    i1, i2, i3, i4 = (names[name] for name in ("C28", "C29", "C210", "C211"))
    require(all(tuple(sorted(pair)) in actual_bonds for pair in ((i1, i2), (i2, i3), (i3, i4))),
            f"POP {residue.id}: cis descriptor does not follow bonds")
    first, second, third, fourth = (points[index] for index in (i1, i2, i3, i4))
    axis = [third[i] - second[i] for i in range(3)]
    axis_sq = sum(value * value for value in axis)
    require(axis_sq > 1e-8, f"POP {residue.id}: degenerate double bond")
    left = [first[i] - second[i] for i in range(3)]
    right = [fourth[i] - third[i] for i in range(3)]
    ld = sum(left[i] * axis[i] for i in range(3)) / axis_sq
    rd = sum(right[i] * axis[i] for i in range(3)) / axis_sq
    lp = [left[i] - ld * axis[i] for i in range(3)]
    rp = [right[i] - rd * axis[i] for i in range(3)]
    denom = math.sqrt(sum(value * value for value in lp) * sum(value * value for value in rp))
    require(denom > 1e-8 and sum(lp[i] * rp[i] for i in range(3)) / denom > 0.5,
            f"POP {residue.id}: cis alkene")

    phosphorus = points[names["P"]]
    side = "upper" if phosphorus[2] > center else "lower"
    require(phosphorus[2] != center, f"POP {residue.id}: phosphate on side divider")
    tail_z = sum(points[names[name]][2] for name in ("C218", "C316")) / 2
    require((phosphorus[2] - tail_z) * (1 if side == "upper" else -1) > 0,
            f"POP {residue.id}: leaflet orientation")
    return chiral, side, (float(phosphorus[0]), float(phosphorus[1]))


def verify_patch(pdb, reference, stereo_checks, *, source: bool):
    box = pdb.topology.getUnitCellDimensions()
    require(box is not None, "Patch has no unit cell")
    dimensions = box.value_in_unit(unit.angstrom)
    require(all(math.isfinite(float(value)) and value > 0 for value in dimensions),
            "Patch cell is not finite and positive")
    lipids = [residue for residue in pdb.topology.residues() if residue.name == "POP"]
    waters = [residue for residue in pdb.topology.residues() if residue.name == "HOH"]
    require(len(waters) == 5120 and len(lipids) == (128 if source else 126) and
            pdb.topology.getNumResidues() == len(waters) + len(lipids), "Patch species/count")
    expected_bonds = bonds(reference.topology, next(reference.topology.residues()))
    sides = Counter()
    inverted = []
    positions = {}
    for residue in lipids:
        chiral, side, phosphorus_xy = lipid_observation(
            pdb, residue, reference, expected_bonds, stereo_checks, float(dimensions[2]) / 2)
        sides[side] += 1
        positions[residue.id] = phosphorus_xy
        if chiral > 0:
            inverted.append(residue.id)
    require(sides == Counter({"upper": 64 if source else 63,
                              "lower": 64 if source else 63}), "Patch leaflet counts")
    require(inverted == (["109"] if source else []), "Patch POPC stereochemistry")
    return lipids, positions, dimensions


def verify_exact_survivors(source, derived) -> None:
    old_atoms = [atom for atom in source.topology.atoms()
                 if atom.residue.name != "POP" or atom.residue.id not in REMOVED_IDS]
    new_atoms = list(derived.topology.atoms())
    require(len(old_atoms) == len(new_atoms) == 32244, "Readback atom count")
    identity = lambda a: (a.residue.name, a.residue.id, a.name, a.element.symbol)
    require(all(identity(a) == identity(b) for a, b in zip(old_atoms, new_atoms)),
            "Readback changed survivor atom order or identity")
    old_points = source.positions.value_in_unit(unit.angstrom)
    new_points = derived.positions.value_in_unit(unit.angstrom)
    require(all(all(float(old_points[a.index][i]) == float(new_points[b.index][i])
                    for i in range(3)) for a, b in zip(old_atoms, new_atoms)),
            "Readback changed survivor coordinates")
    remap = {atom.index: index for index, atom in enumerate(old_atoms)}
    expected_bonds = {
        tuple(sorted((remap[a.index], remap[b.index])))
        for a, b in source.topology.bonds()
        if a.index in remap and b.index in remap
    }
    actual_bonds = {tuple(sorted((a.index, b.index))) for a, b in derived.topology.bonds()}
    require(actual_bonds == expected_bonds and len(actual_bonds) == 26998,
            "Readback changed survivor molecular bonds")
    require(source.topology.getUnitCellDimensions() == derived.topology.getUnitCellDimensions(),
            "Readback changed unit cell")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true", help="write the pinned candidate asset")
    action.add_argument("--check", action="store_true", help="verify the existing asset")
    args = parser.parse_args()

    require(openmm.version.full_version == "8.6.0.dev-c6173db",
            "Installed OpenMM build differs from the pinned provider")
    patch_path = Path(openmm.app.__file__).resolve().parent / "data/POPC.pdb"
    require(patch_path.is_file() and digest(patch_path.read_bytes()) == SOURCE_SHA256,
            "Installed POPC source patch differs from the pinned bytes")
    catalogue = json.loads(CATALOGUE.read_text())
    lipid = next(item for item in catalogue["lipids"] if item["speciesId"] == "POPC")
    reference_path = CATALOGUE.parent / lipid["coordinateTemplatePath"]
    require(digest(reference_path.read_bytes()) == lipid["coordinateTemplateSha256"],
            "POPC catalogue reference bytes differ")
    reference = PDBxFile(str(reference_path))
    source = PDBFile(str(patch_path))
    lipids, xy, box = verify_patch(source, reference, lipid["stereoChecks"], source=True)

    defective = xy["109"]
    upper = [residue for residue in lipids if residue.id != "109" and
             source.positions[next(atom.index for atom in residue.atoms() if atom.name == "P")]
             .value_in_unit(unit.angstrom)[2] > float(box[2]) / 2]
    def xy_distance(residue):
        delta = [xy[residue.id][i] - defective[i] for i in range(2)]
        delta = [value - round(value / float(box[i])) * float(box[i])
                 for i, value in enumerate(delta)]
        return math.hypot(*delta)
    paired = min(upper, key=xy_distance)
    require(paired.id == "62" and abs(xy_distance(paired) - 2.225238751) < 0.001,
            "The paired upper POP is no longer the nearest opposite leaflet site")
    remove = [residue for residue in lipids if residue.id in REMOVED_IDS]
    require(len(remove) == 2, "The precise deletion set is absent")
    modeller = Modeller(source.topology, source.positions)
    modeller.delete(remove)
    buffer = io.StringIO()
    PDBFile.writeFile(modeller.topology, modeller.positions, buffer, keepIds=True)
    derived_bytes = buffer.getvalue().encode("utf-8")
    require(digest(derived_bytes) == DERIVED_SHA256,
            "Derived candidate bytes differ from the pinned recipe")

    if args.write:
        OUTPUT.write_bytes(derived_bytes)
    else:
        require(OUTPUT.is_file() and OUTPUT.read_bytes() == derived_bytes,
                "Candidate asset is missing or differs from reproducible bytes")
    readback = PDBFile(str(OUTPUT))
    verify_exact_survivors(source, readback)
    verify_patch(readback, reference, lipid["stereoChecks"], source=False)
    print(json.dumps({"sourceSha256": SOURCE_SHA256, "derivedSha256": DERIVED_SHA256,
                      "deletedPopResidues": sorted(REMOVED_IDS, key=int),
                      "pairedPhosphateXyDistanceAngstrom": round(xy_distance(paired), 6),
                      "retainedPopcPerLeaflet": 63, "retainedWaters": 5120,
                      "readbackAtoms": 32244, "readbackBonds": 26998,
                      "standing": "candidate; not scientifically qualified"}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, StopIteration) as error:
        sys.exit(f"POPC candidate patch refusal: {error}")
