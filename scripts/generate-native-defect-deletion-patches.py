#!/usr/bin/env python3
"""Reproduce bounded DOPC/DPPC native-patch deletion candidates.

The assets are candidate starting patches.  Reproduction and exact chemical
identity do not qualify their physical starting state or downstream route.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import sys

import openmm
import openmm.app
from openmm import unit
from openmm.app import Modeller, PDBFile, PDBxFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))
from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import (
    _native_lipid_side, _native_molecule_bonds, _native_molecule_matches,
)
from ProteinInMembraneSystem.worker.exchange import WorkError

SETTINGS = {
    "DOPC": {"residueName": "DOP", "removed": ("4", "90"),
             "sourceSha256": "05bc09ae7c7ac6cbdb221868fede000f310b86999105fda1cfbbb262135c7f10",
             "derivedSha256": "f1fe152d3c9bf22ac5eec70ff3d207d7506816ee03dcc910766b2f29d7fbddee"},
    "DPPC": {"residueName": "DPP", "removed": ("6", "70"),
             "sourceSha256": "6020f03c08d38b8469c3074da826fdf979a1ea627354459a092e27558af503d6",
             "derivedSha256": "c961627244da03a88aad66124fb8e877f84fa750e80b4f430f20137a975e1926"},
}
PINNED_HEADER = "REMARK   1 CREATED WITH OPENMM 8.6, 2026-09-27\n"


def require(yes: bool, message: str) -> None:
    if not yes:
        raise ValueError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify_species(species: str, pdb: PDBFile, reference: PDBxFile, representation: dict,
                   allowed_bad: set[str]) -> tuple[dict[str, str], int]:
    raw = SETTINGS[species]
    residues = [residue for residue in pdb.topology.residues()
                if residue.name == raw["residueName"]]
    require(len(residues) == (128 if allowed_bad else 126),
            f"{species} has an unexpected source/derived lipid count")
    expected_bonds = _native_molecule_bonds(reference.topology,
                                             next(reference.topology.residues()))
    positions = [tuple(float(value) for value in point.value_in_unit(unit.angstrom))
                 for point in pdb.positions]
    center = float(pdb.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)[2]) / 2
    actual_bad = set()
    sides = Counter()
    for residue in residues:
        try:
            atoms = _native_molecule_matches(
                residue, reference, species, pdb.positions,
                representation["stereoChecks"],
                _native_molecule_bonds(pdb.topology, residue), expected_bonds)
        except WorkError as error:
            require(error.code == "providerMismatch", f"{species} has unexpected chemistry failure")
            actual_bad.add(residue.id)
        else:
            sides[_native_lipid_side(atoms, positions, center, species)] += 1
    require(actual_bad == allowed_bad, f"{species} defects differ: {sorted(actual_bad)}")
    if not allowed_bad:
        require(sides == Counter({"upper": 63, "lower": 63}),
                f"{species} derived patch lost balanced molecular orientation")
    return {residue.id: "bad" if residue.id in actual_bad else "good" for residue in residues}, len(residues)


def verify_survivors(source: PDBFile, derived: PDBFile, removed: set[str], residue_name: str) -> None:
    remove_indices = {atom.index for residue in source.topology.residues()
                      if residue.name == residue_name and residue.id in removed
                      for atom in residue.atoms()}
    old_atoms = [atom for atom in source.topology.atoms() if atom.index not in remove_indices]
    new_atoms = list(derived.topology.atoms())
    require(len(old_atoms) == len(new_atoms), "Derivative changed surviving atom count")
    for old, new in zip(old_atoms, new_atoms):
        require((old.name, old.element.symbol, old.residue.name, old.residue.id,
                 old.residue.chain.id) ==
                (new.name, new.element.symbol, new.residue.name, new.residue.id,
                 new.residue.chain.id), "Derivative changed a surviving atom identity")
    old_points = source.positions.value_in_unit(unit.angstrom)
    new_points = derived.positions.value_in_unit(unit.angstrom)
    require(all(all(float(old_points[atom.index][axis]) == float(new_points[index][axis])
                    for axis in range(3)) for index, atom in enumerate(old_atoms)),
            "Derivative changed surviving coordinates")
    mapping = {atom.index: index for index, atom in enumerate(old_atoms)}
    old_bonds = {tuple(sorted((mapping[first.index], mapping[second.index])))
                 for first, second in source.topology.bonds()
                 if first.index in mapping and second.index in mapping}
    new_bonds = {tuple(sorted((first.index, second.index)))
                 for first, second in derived.topology.bonds()}
    require(old_bonds == new_bonds, "Derivative changed surviving bonds")
    require(source.topology.getPeriodicBoxVectors() == derived.topology.getPeriodicBoxVectors(),
            "Derivative changed the periodic cell")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("species", choices=sorted(SETTINGS))
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    require(openmm.version.full_version == "8.6.0.dev-c6173db", "Installed OpenMM version changed")
    setting = SETTINGS[args.species]
    path = Path(openmm.app.__file__).resolve().parent / "data" / f"{args.species}.pdb"
    require(path.is_file() and digest(path.read_bytes()) == setting["sourceSha256"],
            "Installed native source bytes changed")
    catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
    representation = next(item for item in catalogue["lipids"]
                          if item["speciesId"] == args.species)
    reference_path = ROOT / "config/policies" / representation["coordinateTemplatePath"]
    require(digest(reference_path.read_bytes()) == representation["coordinateTemplateSha256"],
            "Exact single-molecule reference bytes changed")
    reference = PDBxFile(str(reference_path))
    source = PDBFile(str(path))
    removed = set(setting["removed"])
    verify_species(args.species, source, reference, representation, removed)
    residues = [residue for residue in source.topology.residues()
                if residue.name == setting["residueName"] and residue.id in removed]
    require(len(residues) == 2, "The exact source defects are missing")
    modeller = Modeller(source.topology, source.positions)
    modeller.delete(residues)
    output = io.StringIO()
    PDBFile.writeFile(modeller.topology, modeller.positions, output, keepIds=True)
    generated = output.getvalue()
    header, newline, rest = generated.partition("\n")
    require(header.startswith("REMARK   1 CREATED WITH OPENMM 8.6, ") and newline,
            "OpenMM PDB header changed")
    data = (PINNED_HEADER + rest).encode("utf-8")
    output_path = ROOT / "config/policies/membrane-templates" / f"{args.species}-OpenMM86-63x63.pdb"
    if args.write:
        output_path.write_bytes(data)
    else:
        require(output_path.is_file() and output_path.read_bytes() == data,
                "Derived patch bytes differ from the reproducible candidate")
    actual_digest = digest(data)
    if setting["derivedSha256"] is not None:
        require(actual_digest == setting["derivedSha256"], "Pinned derivative digest changed")
    derivative = PDBFile(str(output_path))
    verify_survivors(source, derivative, removed, setting["residueName"])
    verify_species(args.species, derivative, reference, representation, set())
    print(json.dumps({"species": args.species, "sourceSha256": setting["sourceSha256"],
                      "derivedSha256": actual_digest, "removedResidueIds": setting["removed"],
                      "retainedPerLeaflet": 63, "standing": "candidate; not qualified"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
