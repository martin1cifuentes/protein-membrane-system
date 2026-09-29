#!/usr/bin/env python3
"""Reproduce the pinned one-residue POPC addMembrane patch from its GRO source.

Run with the repository's identified OpenMM environment. This is an asset
conversion, not a membrane equilibration or a replacement for chemistry checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from openmm import Vec3, unit
from openmm.app import PDBxFile, Topology, element


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "config/policies/membrane-templates"
MANIFEST = ASSETS / "POPC-Lipid21-Zenodo-64x64-map.json"
TEMPLATE = ASSETS / "POPC.cif"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output: Path) -> None:
    manifest = json.loads(MANIFEST.read_text())
    source = ASSETS / manifest["sourceFile"]
    if digest(source) != manifest["sourceSha256"] or digest(TEMPLATE) != manifest[
            "selectedCoordinateTemplateSha256"]:
        raise ValueError("The exact source or selected Lipid21 reference template changed")
    lines = source.read_text().splitlines()
    atom_count = int(lines[1])
    lipid_count = manifest["lipids"]
    water_count = manifest["waters"]
    per_lipid = manifest["mappedAtomCountPerLipid"]
    if atom_count != lipid_count * per_lipid + water_count * 3 or len(lines) != atom_count + 3:
        raise ValueError("The source molecule and atom populations changed")
    box = [10 * float(value) for value in lines[-1].split()]
    if len(box) != 3 or any(abs(a - b) > 0.0001 for a, b in zip(box,
                                                                  manifest["cellAngstrom"])):
        raise ValueError("The source periodic cell changed")
    atom_lines = lines[2:-1]
    names = [line[10:15].strip() for line in atom_lines]
    source_lipid_names = names[:per_lipid]
    if any(names[index * per_lipid:(index + 1) * per_lipid] != source_lipid_names
           for index in range(lipid_count)):
        raise ValueError("The source lipid atom ordering changed")
    if any(names[lipid_count * per_lipid + index * 3:
                 lipid_count * per_lipid + (index + 1) * 3] != ["O", "H1", "H2"]
           for index in range(water_count)):
        raise ValueError("The source water atom ordering changed")
    points = [tuple(10 * float(line[start:end]) for start, end in
                    ((20, 28), (28, 36), (36, 44))) for line in atom_lines]
    permutation = manifest["sourceAtomIndexForTemplateAtom"]
    if len(permutation) != per_lipid or sorted(permutation) != list(range(per_lipid)):
        raise ValueError("The recorded POPC atom permutation is incomplete")
    template = PDBxFile(str(TEMPLATE))
    target_atoms = list(template.topology.atoms())
    if len(target_atoms) != per_lipid:
        raise ValueError("The selected POPC template atom count changed")
    template_bonds = [(bond[0].index, bond[1].index, bond.type, bond.order)
                      for bond in template.topology.bonds()]
    topology = Topology()
    coordinates = []
    lipid_chain = topology.addChain("L")
    for index in range(lipid_count):
        residue = topology.addResidue("POPC", lipid_chain, str(index + 1))
        atoms = [topology.addAtom(atom.name, atom.element, residue)
                 for atom in target_atoms]
        for first, second, kind, order in template_bonds:
            topology.addBond(atoms[first], atoms[second], kind, order)
        coordinates.extend(points[index * per_lipid + source_index]
                           for source_index in permutation)
    water_chain = topology.addChain("W")
    for index in range(water_count):
        residue = topology.addResidue("HOH", water_chain, str(index + 1))
        atoms = [topology.addAtom(name, element.get_by_symbol(symbol), residue)
                 for name, symbol in (("O", "O"), ("H1", "H"), ("H2", "H"))]
        topology.addBond(atoms[0], atoms[1])
        topology.addBond(atoms[0], atoms[2])
        first = lipid_count * per_lipid + index * 3
        coordinates.extend(points[first:first + 3])
    topology.setUnitCellDimensions(Vec3(*box) * unit.angstrom)
    with output.open("w") as stream:
        PDBxFile.writeFile(topology, coordinates * unit.angstrom, stream, keepIds=True)
    # OpenMM writes the current date in the mmCIF header. Reproducibility of
    # this qualified molecular asset requires the recorded conversion date.
    written = output.read_text()
    written, replacements = re.subn(
        r"(?m)^# Created with OpenMM 8\.6, \d{4}-\d{2}-\d{2}$",
        f"# Created with OpenMM 8.6, {manifest['mappedHeaderDate']}", written)
    if replacements != 1:
        raise ValueError("The converted mmCIF header differs from its recorded OpenMM format")
    output.write_text(written)
    if digest(output) != manifest["mappedSha256"]:
        raise ValueError("The converted bytes differ from the qualified mapped patch")
    print(f"Verified {output}: {manifest['mappedSha256']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
