#!/usr/bin/env python3
"""Generate exact one-molecule Amber19 TIP3P/ion identity references.

The source water coordinates and parameter file come from the pinned OpenMM
installation. Ion reference positions are arbitrary single-atom origins;
OpenMM addMembrane supplies the actual water and ion coordinates. The generated
files are checked against the matching Amber19 TIP3P/Joung-Cheatham templates.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import openmm.app as app
from openmm import Vec3, unit
from openmm.app import ForceField, NoCutoff, PDBFile, PDBxFile, Topology
from openmm.app.element import Element


ROOT = Path(__file__).resolve().parents[1]
OPENMM_DATA = Path(app.__file__).resolve().parent / "data"
WATER_SOURCE = OPENMM_DATA / "tip3p.pdb"
ION_PARAMETERS = OPENMM_DATA / "amber19" / "tip3p.xml"
DESTINATION = ROOT / "config" / "policies" / "membrane-templates"
WATER_SOURCE_SHA256 = "14fd37900d627c0e258d6086a14c6084e4bec1422e9d20fccfc83c3f814fd7ee"
ION_PARAMETERS_SHA256 = "3f4b188dbcb6c02863230eaca231e927fb6bf3307ce947d8a50d0f46f6dd83d9"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit(name: str, residue_name: str, atoms: list[tuple[str, str]],
         coordinates: list[Vec3], bonds: list[tuple[int, int]]) -> None:
    topology = Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue(residue_name, chain)
    inserted = [topology.addAtom(atom_name, Element.getBySymbol(symbol), residue)
                for atom_name, symbol in atoms]
    for first, second in bonds:
        topology.addBond(inserted[first], inserted[second])
    forcefield = ForceField(str(ION_PARAMETERS))
    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
    if system.getNumParticles() != len(atoms):
        raise RuntimeError(f"Amber19 TIP3P did not parameterize {name} completely")
    DESTINATION.mkdir(parents=True, exist_ok=True)
    target = DESTINATION / name
    with target.open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(topology, coordinates * unit.angstrom, stream, keepIds=True)
    observed = PDBxFile(str(target))
    if [(atom.name, atom.element.symbol) for atom in observed.topology.atoms()] != atoms:
        raise RuntimeError(f"The {name} mmCIF round trip changed atom identity")
    print(f"{target.relative_to(ROOT)} sha256={digest(target)} atom_count={len(atoms)}")


def main() -> None:
    if digest(WATER_SOURCE) != WATER_SOURCE_SHA256 or digest(ION_PARAMETERS) != ION_PARAMETERS_SHA256:
        raise RuntimeError("Pinned OpenMM water source or Amber19 TIP3P parameter identity changed")
    source = PDBFile(str(WATER_SOURCE))
    water = next(source.topology.residues())
    water_atoms = list(water.atoms())
    if [(atom.name, atom.element.symbol) for atom in water_atoms] != [
        ("O", "O"), ("H1", "H"), ("H2", "H")
    ]:
        raise RuntimeError("Pinned OpenMM TIP3P source water atom identities changed")
    positions = [Vec3(*position.value_in_unit(unit.angstrom))
                 for position in source.positions[:3]]
    anchor = positions[0]
    centered = [position - anchor for position in positions]
    print(f"OpenMM TIP3P source sha256={digest(WATER_SOURCE)}")
    print(f"Amber19 TIP3P parameters sha256={digest(ION_PARAMETERS)}")
    emit("TIP3P-HOH.cif", "HOH", [("O", "O"), ("H1", "H"), ("H2", "H")],
         centered, [(0, 1), (0, 2)])
    emit("JC-NA.cif", "NA", [("NA", "Na")], [Vec3(0, 0, 0)], [])
    emit("JC-CL.cif", "CL", [("CL", "Cl")], [Vec3(0, 0, 0)], [])


if __name__ == "__main__":
    main()
