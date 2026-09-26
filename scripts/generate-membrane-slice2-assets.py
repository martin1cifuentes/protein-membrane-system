#!/usr/bin/env python3
"""Rebuild the pinned slice-2 coordinate catalogue from primary local assets.

Run with the pinned out/python/bin/python after scripts/build-local.sh.  The output
is deterministic for the identified OpenMM 8.6.0 PDB/force-field bytes and the
CC0 RCSB cholesterol CCD definition kept under config/policies/source-assets.
"""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import gemmi
import numpy as np
from openmm import Vec3, unit
from openmm.app import ForceField, PDBFile, PDBxFile, Topology, element


ROOT = Path(__file__).resolve().parent.parent
POLICIES = ROOT / "config" / "policies"
OPENMM_DATA = ROOT / "out" / "python" / "lib" / "python3.11" / "site-packages" / "openmm" / "app" / "data"
LIPID21_XML = OPENMM_DATA / "amber19" / "lipid21.xml"
FORCE_FIELD_SHA = "4ed8dc7485a3df14f55e1d9d12fb2a99813495293ba07c2d910fe0f66185e23a"
CLR_CIF = POLICIES / "source-assets" / "CLR.cif"
CLR_SHA = "0ee512527d26bd087cb19ffe657dada7aae0e35af4595249d126f4b928e6a811"

# The PDBs are OpenMM's seven bundled membrane patches.  Their hashes also fix
# the 2017 source geometries rather than assuming future wheels are identical.
PATCHES = {
    "DLPC": ("59d37cac8e13fdd6d812b1e195706e67f76eb548ed027f1d5639379a70aaa075", 933.23),
    "DLPE": ("ad4d1243e8bad41639dd648d2252407e3c945d72e0599f22bbc112ac8805b4c1", 907.0),
    "DMPC": ("19d8513a792ddf6f34edc16028da4410dc6d686945bede006ab1e5a27accb4ba", 1039.32),
    "DPPC": ("6020f03c08d38b8469c3074da826fdf979a1ea627354459a092e27558af503d6", 1163.50),
    "DOPC": ("05bc09ae7c7ac6cbdb221868fede000f310b86999105fda1cfbbb262135c7f10", 1232.84),
    "POPC": ("a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29", 1190.32),
    "POPE": ("71314d07ab9a28a293f724869c876ccf541c5c5eb35aa66f8357b4bb08884c36", 1129.94),
}
VOLUME_CONTEXT = {
    "DLPC": "Lipid21 pure DLPC simulation at 303 K",
    "DLPE": "fully hydrated fluid DLPE at 35 °C in Nagle and Tristram-Nagle Table 6",
    "DMPC": "Lipid21 pure DMPC simulation at 303 K",
    "DPPC": "Lipid21 pure DPPC simulation at 323 K",
    "DOPC": "Lipid21 pure DOPC simulation at 303 K",
    "POPC": "Lipid21 pure POPC simulation at 303 K",
    "POPE": "Lipid21 pure POPE simulation at 310 K",
}


def verify_hash(path: Path, expected: str) -> None:
    actual = sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"Pinned asset digest changed: {path}: {actual}")


def template_records(species: str) -> tuple[dict[str, tuple[str, float]], set[tuple[str, str]]]:
    root = ET.parse(LIPID21_XML).getroot()
    residue = next((item for item in root.find("Residues") if item.get("name") == species), None)
    if residue is None:
        raise ValueError(f"The pinned Lipid21 XML has no complete {species} template")
    atoms = {item.get("name"): (item.get("type"), float(item.get("charge")))
             for item in residue.findall("Atom")}
    bonds = {tuple(sorted((item.get("atomName1"), item.get("atomName2"))))
             for item in residue.findall("Bond")}
    return atoms, bonds


def lipid_atom_name(original: str) -> str:
    # Older PDB patches wrote five-character acyl atom names into a four-column
    # field: C210 was parsed as 0C21.  This is the uniquely identified inverse.
    match = re.fullmatch(r"([0-9])C([23])1", original)
    return f"C{match.group(2)}1{match.group(1)}" if match else original


def make_topology(species: str, atoms: list[tuple[str, str]],
                  bonds: list[tuple[str, str, int | None]]) -> Topology:
    topology = Topology()
    residue = topology.addResidue(species, topology.addChain("A"), id="1")
    output = {name: topology.addAtom(name, element.get_by_symbol(symbol), residue)
              for name, symbol in atoms}
    if len(output) != len(atoms):
        raise ValueError(f"Duplicate molecular atom names in {species}")
    for a, b, order in bonds:
        topology.addBond(output[a], output[b], order=order)
    return topology


def tetrahedral_sign(coordinates: dict[str, np.ndarray], names: list[str]) -> str:
    vectors = [coordinates[name] - coordinates[names[3]] for name in names[:3]]
    volume = float(np.linalg.det(np.array(vectors)))
    if abs(volume) < 1.0:
        raise ValueError(f"Degenerate tetrahedral geometry: {names}: {volume}")
    return "positive" if volume > 0 else "negative"


def cis_or_trans(coordinates: dict[str, np.ndarray], names: list[str]) -> str:
    before, first, second, after = [coordinates[name] for name in names]
    axis = second - first
    left = np.cross(axis, before - first)
    right = np.cross(axis, after - second)
    cosine = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
    if abs(cosine) < 0.7:
        raise ValueError(f"Indeterminate alkene geometry: {names}: {cosine}")
    return "cis" if cosine > 0 else "trans"


def stereo_checks(species: str, coordinates: dict[str, np.ndarray],
                  ccd_centers: dict[str, str] | None = None,
                  adjacency: dict[str, list[str]] | None = None) -> list[dict]:
    if species == "CHL1":
        assert ccd_centers is not None and adjacency is not None
        checks = []
        for center, configuration in ccd_centers.items():
            names = sorted(adjacency[center])
            if len(names) != 4:
                raise ValueError(f"CCD center {center} does not have four substituents")
            checks.append({"kind": "tetrahedral", "atomNames": names,
                           "expected": tetrahedral_sign(coordinates, names),
                           "centerAtomName": center, "sourceConfiguration": configuration})
        if len(checks) != 8:
            raise ValueError("The CCD cholesterol identity must have eight chiral centers")
        return checks
    checks = [{"kind": "tetrahedral", "atomNames": ["O21", "C1", "C3", "HS"],
               "expected": "negative", "centerAtomName": "C2", "sourceConfiguration": "R"}]
    if tetrahedral_sign(coordinates, checks[0]["atomNames"]) != "negative":
        raise ValueError(f"{species} does not have the selected sn-glycerol C2 configuration")
    if species in {"DOPC", "POPC", "POPE"}:
        names = ["C28", "C29", "C210", "C211"]
        if cis_or_trans(coordinates, names) != "cis":
            raise ValueError(f"{species} oleoyl chain is not 9Z")
        checks.append({"kind": "alkene", "atomNames": names, "expected": "cis"})
    if species == "DOPC":
        names = ["C38", "C39", "C310", "C311"]
        if cis_or_trans(coordinates, names) != "cis":
            raise ValueError("DOPC second oleoyl chain is not 9Z")
        checks.append({"kind": "alkene", "atomNames": names, "expected": "cis"})
    return checks


def check_and_write(species: str, topology: Topology, positions: list[Vec3],
                    expected_stereo: list[dict]) -> tuple[str, int, int]:
    expected_atoms, expected_bonds = template_records(species)
    actual_atoms = {atom.name for atom in topology.atoms()}
    actual_bonds = {tuple(sorted((a.name, b.name))) for a, b in topology.bonds()}
    if actual_atoms != set(expected_atoms) or actual_bonds != expected_bonds:
        raise ValueError(f"{species} atom names or bonds disagree with exact Lipid21 XML")
    force_field = ForceField(str(LIPID21_XML))
    matches = force_field.getMatchingTemplates(topology)
    if len(matches) != 1 or matches[0].name != species:
        raise ValueError(f"{species} did not match its exact parameter template")
    system = force_field.createSystem(topology)
    if system.getNumParticles() != len(actual_atoms):
        raise ValueError(f"{species} parameterized particle count differs")
    coordinate_map = {atom.name: np.array(positions[atom.index], dtype=float)
                      for atom in topology.atoms()}
    for check in expected_stereo:
        observed = (tetrahedral_sign(coordinate_map, check["atomNames"])
                    if check["kind"] == "tetrahedral" else
                    cis_or_trans(coordinate_map, check["atomNames"]))
        if observed != check["expected"]:
            raise ValueError(f"{species} stereochemistry disagrees with identified source")
    output = POLICIES / "membrane-templates" / f"{species}.cif"
    buffer = StringIO()
    PDBxFile.writeFile(topology, unit.Quantity(positions, unit.angstrom), buffer, keepIds=True)
    output.write_text(buffer.getvalue(), encoding="utf-8")
    reread = PDBxFile(str(output))
    if ({atom.name for atom in reread.topology.atoms()} != actual_atoms or
            {tuple(sorted((a.name, b.name))) for a, b in reread.topology.bonds()} != expected_bonds or
            force_field.getMatchingTemplates(reread.topology)[0].name != species):
        raise ValueError(f"{species} mmCIF did not round-trip exactly")
    return (sha256(output.read_bytes()).hexdigest(), len(actual_atoms),
            next((atom.index + 1 for atom in topology.atoms() if atom.name == "P"), 0))


def make_lipid(species: str) -> dict:
    source_sha, volume = PATCHES[species]
    source = OPENMM_DATA / f"{species}.pdb"
    verify_hash(source, source_sha)
    pdb = PDBFile(str(source))
    residues = [residue for residue in pdb.topology.residues() if residue.name != "HOH"]
    if len(residues) != 128:
        raise ValueError(f"{species} source patch has an unexpected lipid population")
    first = residues[0]
    source_atoms = list(first.atoms())
    names = {atom: lipid_atom_name(atom.name) for atom in source_atoms}
    atoms = [(names[atom], atom.element.symbol) for atom in source_atoms]
    bonds = [(names[a], names[b], getattr(bond, "order", None))
             for bond in pdb.topology.bonds() for a, b in [bond]
             if a in names and b in names]
    topology = make_topology(species, atoms, bonds)
    positions = [Vec3(*pdb.positions[atom.index].value_in_unit(unit.angstrom))
                 for atom in source_atoms]
    coordinates = {names[atom]: np.array(positions[index], dtype=float)
                   for index, atom in enumerate(source_atoms)}
    checks = stereo_checks(species, coordinates)
    digest, atom_count, head_index = check_and_write(species, topology, positions, checks)
    box = pdb.topology.getPeriodicBoxVectors()
    if box is None:
        raise ValueError(f"{species} source patch has no measured XY area")
    axes = [vector.value_in_unit(unit.angstrom) for vector in box]
    area = round(float(np.linalg.norm(np.cross(axes[0], axes[1]))) / 64.0, 3)
    return representation(species, f"Lipid21:{species}:sn-2R", "lipid", digest,
                          atom_count, area, volume, [head_index], checks,
                          ["Coordinates are one selected molecule from the pinned OpenMM 8.6.0 patch; exact Lipid21 template and sn-2R geometry were checked.",
                           "The XY footprint is the source patch's area divided by 64 lipids per leaflet, an initial packing estimate, not a validated area of this species in every mixture or at 303 K.",
                           f"The listed molecular volume is from {VOLUME_CONTEXT[species]}; it is not a measured mixture partial volume, a 303 K value for other temperatures, or a completed-system observation."])


def make_cholesterol() -> dict:
    verify_hash(CLR_CIF, CLR_SHA)
    block = gemmi.cif.read_file(str(CLR_CIF)).sole_block()
    rows = list(block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
                            "_chem_comp_atom.pdbx_stereo_config",
                            "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
                            "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
                            "_chem_comp_atom.pdbx_model_Cartn_z_ideal"]))
    raw = [(row[0], row[1], row[2], Vec3(*(float(row[i]) for i in (3, 4, 5))))
           for row in rows]
    source_bonds = [(row[0], row[1], 2 if row[2] == "DOUB" else 1)
                    for row in block.find(["_chem_comp_bond.atom_id_1",
                                           "_chem_comp_bond.atom_id_2",
                                           "_chem_comp_bond.value_order"])]
    atom_elements = {name: symbol for name, symbol, _, _ in raw}
    adjacency = defaultdict(list)
    for a, b, _ in source_bonds:
        adjacency[a].append(b)
        adjacency[b].append(a)
    target_atoms, target_bonds = template_records("CHL1")
    target_hydrogens = defaultdict(list)
    for a, b in target_bonds:
        for hydrogen, parent in ((a, b), (b, a)):
            if hydrogen.startswith("H") and not parent.startswith("H"):
                target_hydrogens[parent].append(hydrogen)
    mapping = {name: ("O3" if name == "O1" else name)
               for name, symbol, _, _ in raw if symbol != "H"}
    for parent, source_names in adjacency.items():
        if atom_elements[parent] == "H":
            continue
        source_hydrogens = sorted(name for name in source_names if atom_elements[name] == "H")
        target_names = sorted(target_hydrogens[mapping[parent]])
        if len(source_hydrogens) != len(target_names):
            raise ValueError(f"Cholesterol hydrogen mapping differs at {parent}")
        mapping.update(zip(source_hydrogens, target_names))
    if len(mapping) != 74 or len(set(mapping.values())) != 74:
        raise ValueError("Cholesterol atom-name mapping is incomplete or ambiguous")
    atoms = [(mapping[name], symbol) for name, symbol, _, _ in raw]
    topology = make_topology("CHL1", atoms,
                             [(mapping[a], mapping[b], order) for a, b, order in source_bonds])
    positions = [coords for _, _, _, coords in raw]
    coordinate_map = {mapping[name]: np.array(coords, dtype=float)
                      for name, _, _, coords in raw}
    mapped_adjacency = {mapping[center]: [mapping[name] for name in neighbors]
                        for center, neighbors in adjacency.items()}
    centers = {mapping[name]: configuration for name, _, configuration, _ in raw
               if configuration in {"R", "S"}}
    checks = stereo_checks("CHL1", coordinate_map, centers, mapped_adjacency)
    digest, atom_count, _ = check_and_write("CHL1", topology, positions, checks)
    return representation("CHL1", "Lipid21:CHL1:cholesterol:CCD-CLR", "sterol",
                          digest, atom_count, 0.0, 630.0, [], checks,
                          ["Coordinates derive from the CC0 RCSB CCD CLR ideal model; its eight declared stereocenters are preserved by the atom-name map and mmCIF round trip.",
                           "A transferable cholesterol footprint for mixed or asymmetric packing is not established; areaPerMoleculeAngstromSquared is zero and construction must not infer a viable packing area from it.",
                           "630 Å³ is an experimental approximate partial molecular volume for cholesterol in fluid phosphatidylcholine bilayers, not a universal local or mixture volume."])


def representation(species: str, chemistry: str, category: str, digest: str,
                   atoms: int, area: float, volume: float, head_indices: list[int],
                   checks: list[dict], limitations: list[str]) -> dict:
    return {"speciesId": species, "chemistryId": chemistry, "category": category,
            "templatePath": "../../out/python/lib/python3.11/site-packages/openmm/app/data/amber19/lipid21.xml",
            "templateSha256": FORCE_FIELD_SHA,
            "coordinateTemplatePath": f"membrane-templates/{species}.cif",
            "coordinateTemplateSha256": digest, "atomCount": atoms,
            "netChargeElementary": 0.0, "areaPerMoleculeAngstromSquared": area,
            "volumeAngstromCubed": volume, "headAtomIndices": head_indices,
            "forceFieldFamily": "Lipid21", "forceFieldVersion": "8.6.0",
            "stereoChecks": [{key: value for key, value in check.items()
                              if key in {"kind", "atomNames", "expected"}}
                             for check in checks],
            "limitations": limitations}


def main() -> None:
    verify_hash(LIPID21_XML, FORCE_FIELD_SHA)
    species = [make_lipid(name) for name in PATCHES]
    species.append(make_cholesterol())
    catalogue = json.loads((POLICIES / "protein-slice1.json").read_text(encoding="utf-8"))
    catalogue["version"] = "membrane-slice2-1.0.0"
    catalogue["evidenceReferences"].extend([
        "https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addMembrane",
        "https://www.rcsb.org/ligand/CLR",
        "https://www.rcsb.org/pages/usage-policy",
        "https://pubs.acs.org/doi/10.1021/acs.jctc.1c01217"])
    catalogue["lipids"] = species
    catalogue["membranePolicies"] = [{
        "id": "lipid21-neutral-planar-physical-leaflets", "version": "1.0.0",
        "coveredSpeciesIds": [item["speciesId"] for item in species],
        "allowsMixedLeaflets": True, "allowsAsymmetricLeaflets": True,
        "coversArbitraryCoherentFractions": True,
        "evidenceReferences": [
            "https://pubs.acs.org/doi/10.1021/acs.jctc.1c01217",
            "https://lipid.phys.cmu.edu/bba/bba.pdf",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2695672/",
            "https://docs.openmm.org/latest/userguide/application/02_running_sims.html#force-fields",
            "docs/Protein-in-membrane-system/requirements-elicitation/Research-Family-0002-Membrane-Model-and-Study-Conditions/2-Valid State Space.md",
            "config/policies/README.md"],
        "limitations": [
            "Membrane-local assessment establishes exact neutral lipid or sterol molecular representations and coherent independently specified physical leaflets only; it does not predict phase, miscibility, native biological sidedness, leaflet stress, or completed-system qualification.",
            "Arbitrary coherent fractions are accepted as stated intentions and exact representable combinations, not as blanket evidence that every mixture or asymmetry is scientifically suitable.",
            "CHL1 has no qualified construction footprint in this catalogue; a cholesterol-containing assessed model needs further attempt-specific packing qualification before construction.",
            "The fixed target is 0.15 M NaCl and optional equilibration target is 303 K; published single-component volume or area references at other temperatures cannot be treated as those achieved conditions."]
    }]
    output = POLICIES / "protein-membrane-slice2.json"
    output.write_text(json.dumps(catalogue, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for item in species:
        print(item["speciesId"], item["atomCount"], item["coordinateTemplateSha256"],
              len(item["stereoChecks"]))
    print(output)


if __name__ == "__main__":
    main()
