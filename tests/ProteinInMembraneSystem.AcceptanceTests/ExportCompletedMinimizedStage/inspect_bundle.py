"""Independent read-back of a downloaded completed-stage bundle.

Uses Gemmi and OpenMM installed in the identified scientific environment. It
does not call the product's export worker or trust its readBackMatched flag.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import gemmi
import numpy as np
from openmm import HarmonicBondForce, XmlSerializer, unit


REQUIRED = {"structure.cif", "topology.json", "system.xml", "state.xml", "manifest.json"}
SOURCE_6QWR_SHA256 = "f4c1503a60321c0cfe513e8211e43e20e2199aa5f71f22429e0e0ae8c97b0779"
DMPC_PATCH_SHA256 = "19d8513a792ddf6f34edc16028da4410dc6d686945bede006ab1e5a27accb4ba"


def cell_lengths_angles(vectors: np.ndarray):
    lengths = np.linalg.norm(vectors, axis=1)
    angles = []
    for left, right in ((1, 2), (0, 2), (0, 1)):
        cosine = np.dot(vectors[left], vectors[right]) / (lengths[left] * lengths[right])
        angles.append(math.degrees(math.acos(float(np.clip(cosine, -1, 1)))))
    return lengths, np.array(angles)


def inspect(path: Path, stage_id: str, attempt_id: str, revision_id: str,
            assessment_id: str) -> dict:
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        assert REQUIRED <= names, f"Missing bundle entries: {sorted(REQUIRED - names)}"
        assert archive.testzip() is None, "An archived entry failed CRC read-back"
        entries = {name: archive.read(name) for name in REQUIRED}

    manifest = json.loads(entries["manifest.json"])
    assert manifest["stage"]["id"] == stage_id
    assert manifest["attempt"]["id"] == attempt_id
    assert manifest["assessment"]["id"] == assessment_id
    assert manifest["assessment"]["stageId"] == stage_id
    assert manifest["assessment"]["qualification"] == "Indeterminate"
    # The portable account must identify both original study and molecular
    # inputs; the exact nested attribution is checked against the owner schema.
    assert manifest["stage"]["kind"] == "Minimization"
    assert manifest["stage"]["policyId"]
    assert manifest["study"]["id"] == revision_id
    assert manifest["preparedProtein"] and manifest["membrane"] and manifest["placement"]
    assert manifest["construction"] and manifest["forceFieldAssets"]
    assert manifest["lineage"] and manifest["attribution"]
    assert manifest["findings"] is not None and manifest["limitations"] is not None
    lineage = manifest["lineage"]
    assert lineage["sourceCoordinateSha256"] == SOURCE_6QWR_SHA256
    assert lineage["sourceModelIndex"] == 0  # First RCSB model, zero-indexed in the product.
    assert lineage["sourceAccession"] == "6QWR"
    assert manifest["preparedProtein"]["source"]["sha256"] == SOURCE_6QWR_SHA256
    assert manifest["attempt"]["nativePatchSha256"] == DMPC_PATCH_SHA256
    assert manifest["construction"]["provider"]["providerVersion"] == "8.6.0.dev-c6173db"
    asset_ids = {item["id"] for item in manifest["forceFieldAssets"]}
    assert asset_ids == {"openmm-amber19-protein-ff19sb",
                         "openmm-amber19-lipid21",
                         "openmm-amber19-tip3p-jc-ions"}
    attribution = manifest["attribution"]
    data = attribution["incorporatedData"]
    references = attribution["referenceInputs"]
    tools = attribution["toolsUsed"]
    assert len(data) >= 5 and len(references) >= 1 and len(tools) >= 2
    for item in data + references + tools:
        for field in ("role", "name", "version", "sourceUri", "rights",
                      "rightsUrl", "citation", "useLimitations", "relationship"):
            assert isinstance(item[field], str) and item[field].strip(), (item.get("name"), field)
    for item in data + references:
        assert isinstance(item["sha256"], str) and len(item["sha256"]) == 64
        assert all(character in "0123456789abcdef" for character in item["sha256"].lower())
    for item in tools:
        digest = item["sha256"]
        assert digest is None or (isinstance(digest, str) and len(digest) == 64 and
                                  all(character in "0123456789abcdef" for character in digest.lower()))
    source_attribution = next(item for item in data
                              if item["role"] == "protein source coordinates")
    assert source_attribution["name"] == "6QWR"
    assert source_attribution["sha256"] == SOURCE_6QWR_SHA256
    assert source_attribution["sourceUri"] == "https://www.rcsb.org/structure/6QWR"
    assert "CC0" in source_attribution["rights"]
    assert any(item["sha256"] == DMPC_PATCH_SHA256 and
               item["role"] == "native lipid coordinate patch" for item in data)
    assert any(item["role"] == "assessed molecular reference" and
               "not inserted" in item["useLimitations"].lower() for item in references)
    assert any(item["role"] == "orientation tool" and item["name"] == "PPM 2.0" and
               item["sha256"] is not None
               for item in tools)
    assert any("OpenMM" in item["name"] and item["role"] == "construction tool" and
               item["sha256"] is None and "not recorded" in item["useLimitations"]
               for item in tools)

    topology = json.loads(entries["topology.json"])
    atoms = topology["atoms"]
    residues = topology["residues"]
    bonds = topology["bonds"]
    atom_count = len(atoms)
    assert atom_count > 0 and len(bonds) > 0
    assert manifest["stage"]["atomCount"] == atom_count
    assert manifest["molecularIdentity"]["atomCount"] == atom_count
    mapped_atoms = manifest["correspondence"]["atoms"]
    assert len(mapped_atoms) == atom_count
    artifact_entries = {item["name"]: item for item in manifest["artifacts"]}
    assert set(artifact_entries) == REQUIRED - {"manifest.json"}
    for name, declared in artifact_entries.items():
        assert declared["sha256"] == hashlib.sha256(entries[name]).hexdigest()
        assert declared["byteLength"] == len(entries[name])
    assert manifest["molecularIdentity"]["coordinatesSha256"] == \
        artifact_entries["structure.cif"]["sha256"]
    bond_pairs = set()
    for bond in bonds:
        left, right = bond["atomIndices"]
        assert 0 <= left < atom_count and 0 <= right < atom_count and left != right
        pair = (min(left, right), max(left, right))
        assert pair not in bond_pairs, f"Duplicate topology bond: {pair}"
        bond_pairs.add(pair)

    block = gemmi.cif.read_string(entries["structure.cif"].decode("utf-8")).sole_block()
    table = block.find("_atom_site.", [
        "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num",
    ])
    assert len(table) == atom_count
    cif_positions = np.empty((atom_count, 3), dtype=float)
    models = set()
    labelled_atoms = defaultdict(list)
    molecule_roles = Counter()
    for index, row in enumerate(table):
        atom = atoms[index]
        residue = residues[atom["residueIndex"]]
        chain = topology["chains"][residue["chainIndex"]]
        assert (row[0], row[1], row[2]) == (
            atom["element"], atom["name"], residue["name"]), \
            f"Atom order or identity changed at row {index}"
        assert (row[3], row[4], "" if row[5] in (".", "?") else row[5]) == (
            chain["id"], residue["id"], residue["insertionCode"]), \
            f"Atom membership changed at row {index}"
        mapped = mapped_atoms[index]
        assert mapped["resultAtomIndex"] == index and mapped["element"] == atom["element"]
        assert mapped["resultAtomId"] and mapped["moleculeRole"]
        molecule_roles[mapped["moleculeRole"]] += 1
        labelled_atoms[(row[3], row[4], row[1])].append(index)
        cif_positions[index] = (float(row[6]), float(row[7]), float(row[8]))
        models.add(row[9])
    assert models == {"1"}, f"Unexpected molecular models: {models}"
    assert all(molecule_roles[role] > 0 for role in ("protein", "lipid", "water", "ion"))
    connections = block.find("_struct_conn.", [
        "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
        "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id",
    ])
    assert len(connections) > 0
    cif_bonds = set()
    for row in connections:
        left = labelled_atoms[(row[0], row[1], row[2])]
        right = labelled_atoms[(row[3], row[4], row[5])]
        assert len(left) == len(right) == 1, "An mmCIF bond endpoint is missing or ambiguous"
        pair = (min(left[0], right[0]), max(left[0], right[0]))
        assert pair not in cif_bonds
        cif_bonds.add(pair)
    assert cif_bonds <= bond_pairs, "An mmCIF bond contradicts the explicit topology"

    system = XmlSerializer.deserialize(entries["system.xml"].decode("utf-8"))
    state = XmlSerializer.deserialize(entries["state.xml"].decode("utf-8"))
    assert system.getNumParticles() == atom_count
    constrained_pairs = {tuple(sorted(system.getConstraintParameters(index)[:2]))
                         for index in range(system.getNumConstraints())}
    parameter_bonds = set()
    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            parameter_bonds.update(tuple(sorted(force.getBondParameters(index)[:2]))
                                   for index in range(force.getNumBonds()))
    realized_pairs = constrained_pairs | parameter_bonds
    assert bond_pairs <= realized_pairs, "A topology bond lacks a System bond or constraint"
    additional_constraints = realized_pairs - bond_pairs
    water_count = sum(residue["name"] == "HOH" for residue in residues)
    assert len(additional_constraints) == water_count
    for left, right in additional_constraints:
        left_residue = atoms[left]["residueIndex"]
        right_residue = atoms[right]["residueIndex"]
        assert ((left, right) in constrained_pairs and left_residue == right_residue and
                residues[left_residue]["name"] == "HOH" and
                atoms[left]["element"] == atoms[right]["element"] == "H"), \
            "An unexplained System constraint is absent from the topology"
    state_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    assert len(state_positions) == atom_count
    deviations = np.linalg.norm(cif_positions - state_positions, axis=1)
    max_deviation = float(np.max(deviations))
    assert math.isfinite(max_deviation) and max_deviation <= 0.01, max_deviation

    vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    topology_vectors = np.asarray(topology["boxVectorsAngstrom"], dtype=float)
    assert topology_vectors.shape == (3, 3)
    max_topology_vector_difference = float(np.max(np.abs(topology_vectors - vectors)))
    assert max_topology_vector_difference <= 0.01, max_topology_vector_difference
    lengths, angles = cell_lengths_angles(vectors)
    topology_lengths, topology_angles = cell_lengths_angles(topology_vectors)
    assert float(np.max(np.abs(topology_lengths - lengths))) <= 0.01
    assert float(np.max(np.abs(topology_angles - angles))) <= 0.1
    cif_lengths = np.array([float(block.find_value(f"_cell.length_{name}"))
                            for name in ("a", "b", "c")])
    cif_angles = np.array([float(block.find_value(f"_cell.angle_{name}"))
                           for name in ("alpha", "beta", "gamma")])
    max_length_difference = float(np.max(np.abs(cif_lengths - lengths)))
    max_angle_difference = float(np.max(np.abs(cif_angles - angles)))
    assert max_length_difference <= 0.01, max_length_difference
    assert max_angle_difference <= 0.1, max_angle_difference

    return {
        "stageId": stage_id,
        "attemptId": attempt_id,
        "studyRevisionId": revision_id,
        "assessmentId": assessment_id,
        "qualification": manifest["assessment"]["qualification"],
        "atomCount": atom_count,
        "bondCount": len(bonds),
        "cifBondCount": len(cif_bonds),
        "waterCount": water_count,
        "moleculeRoleCounts": dict(molecule_roles),
        "cellLengthsAngstrom": lengths.tolist(),
        "cellAnglesDegrees": angles.tolist(),
        "coordinateMaxDeviationAngstrom": max_deviation,
        "cellMaxLengthDifferenceAngstrom": max_length_difference,
        "cellMaxAngleDifferenceDegrees": max_angle_difference,
        "topologyMaxVectorDifferenceAngstrom": max_topology_vector_difference,
        "bundleSha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "entrySha256": {name: hashlib.sha256(data).hexdigest()
                        for name, data in entries.items()},
    }


if __name__ == "__main__":
    result = inspect(Path(sys.argv[1]), *sys.argv[2:6])
    print(json.dumps(result))
