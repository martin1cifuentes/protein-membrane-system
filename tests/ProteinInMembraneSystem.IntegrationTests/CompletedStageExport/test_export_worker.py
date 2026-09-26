"""One-request completed-stage export verification against actual molecular files."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

import gemmi
import numpy as np
from openmm import (Context, HarmonicBondForce, Platform, System, Vec3,
                    VerletIntegrator, XmlSerializer, unit)
from openmm.app import PDBxFile, Topology, element


ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src/ProteinInMembrane.Host"
PYTHON = Path(os.environ.get("PIM_TEST_PYTHON", ROOT / "out/python/bin/python"))
SAVED = ROOT / "out/browser-acceptance/slice4/final-state-artifacts"
SAVED_HASHES = ROOT / "out/browser-acceptance/slice4/final-state-expected.json"
TOLERANCES = {
    "coordinateReadBackToleranceAngstrom": 0.01,
    "cellLengthReadBackToleranceAngstrom": 0.01,
    "cellAngleReadBackToleranceDegrees": 0.1,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def invoke(work: Path, payload: dict, request_id: str, timeout: int = 120):
    request = {"requestId": request_id, "operation": "verify_export",
               "workingDirectory": str(work), "payload": payload}
    environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT),
                       PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT), OPENMM_CPU_THREADS="2")
    started = time.monotonic()
    completed = subprocess.run([str(PYTHON), "-m", "ProteinInMembraneSystem.worker"],
                               input=json.dumps(request) + "\n", text=True,
                               capture_output=True, cwd=SOURCE_ROOT, env=environment,
                               timeout=timeout, check=False)
    events = [json.loads(line) for line in completed.stdout.splitlines()]
    assert events, completed.stderr
    terminal = [event for event in events if event["kind"] in {"result", "error"}]
    assert len(terminal) == 1 and events[-1] is terminal[0], completed.stdout + completed.stderr
    assert terminal[0]["requestId"] == request_id
    return completed, terminal[0], time.monotonic() - started


def payload_for(work: Path) -> dict:
    names = {
        "topologyCif": "stage.cif", "topologyJson": "topology.json",
        "systemXml": "system.xml", "stateXml": "state.xml",
    }
    payload = {"studyRevisionId": "revision-1", "attemptId": "attempt-1",
               "stageId": "minimized-stage-1", **TOLERANCES}
    for role, name in names.items():
        payload[f"{role}Path"] = str(work / name)
        payload[f"{role}Sha256"] = sha256(work / name)
    return payload


def write_small_stage(work: Path) -> dict:
    work.mkdir(parents=True)
    topology = Topology()
    box = tuple(Vec3(*point) for point in
                ((4, 0, 0), (0, 4, 0), (0, 0, 4))) * unit.nanometer
    topology.setPeriodicBoxVectors(box)
    chain = topology.addChain("A")
    residue = topology.addResidue("GLY", chain, id="1")
    carbon = topology.addAtom("C", element.carbon, residue, id="1")
    hydrogen = topology.addAtom("H", element.hydrogen, residue, id="2")
    topology.addBond(carbon, hydrogen)
    system = System()
    system.addParticle(12 * unit.dalton)
    system.addParticle(1 * unit.dalton)
    system.setDefaultPeriodicBoxVectors(*box)
    bonds = HarmonicBondForce()
    bonds.addBond(0, 1, 0.1 * unit.nanometer,
                  1000 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(bonds)
    positions = [Vec3(0.1, 0.2, 0.3), Vec3(0.2, 0.2, 0.3)] * unit.nanometer
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    state = context.getState(getPositions=True, getEnergy=True)
    del context, integrator
    with (work / "stage.cif").open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(topology, state.getPositions(), stream, keepIds=True)
    sidecar = {
        "formatVersion": 1,
        "chains": [{"id": "A"}],
        "residues": [{"chainIndex": 0, "id": "1", "name": "GLY", "insertionCode": ""}],
        "atoms": [
            {"residueIndex": 0, "id": "1", "name": "C", "element": "C"},
            {"residueIndex": 0, "id": "2", "name": "H", "element": "H"},
        ],
        "bonds": [{"atomIndices": [0, 1], "order": None}],
        "boxVectorsAngstrom": [[40.0, 0.0, 0.0], [0.0, 40.0, 0.0], [0.0, 0.0, 40.0]],
    }
    (work / "topology.json").write_text(json.dumps(sidecar), encoding="utf-8")
    (work / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
    (work / "state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")
    return payload_for(work)


def rewrite_cif(work: Path, *, first_x_angstrom: float = 0.0,
                cell_a_angstrom: float = 0.0, cell_gamma_degrees: float = 0.0) -> None:
    cif = PDBxFile(str(work / "stage.cif"))
    points = list(cif.positions.value_in_unit(unit.nanometer))
    points[0] = points[0] + Vec3(first_x_angstrom / 10, 0, 0)
    original = cif.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
    gamma_shift = math.tan(math.radians(cell_gamma_degrees)) * 4.0
    box = (Vec3(original[0][0] + cell_a_angstrom / 10, 0, 0),
           Vec3(gamma_shift, 4, 0), Vec3(0, 0, 4))
    cif.topology.setPeriodicBoxVectors(box * unit.nanometer)
    with (work / "stage.cif").open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(cif.topology, points * unit.nanometer, stream, keepIds=True)


def rewrite_state(work: Path, *, first_x_angstrom: float = 0.0,
                  cell_a_angstrom: float = 0.0) -> None:
    system = XmlSerializer.deserialize((work / "system.xml").read_text(encoding="utf-8"))
    original = XmlSerializer.deserialize((work / "state.xml").read_text(encoding="utf-8"))
    positions = list(original.getPositions().value_in_unit(unit.nanometer))
    positions[0] = positions[0] + Vec3(first_x_angstrom / 10, 0, 0)
    box = np.asarray(original.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    box[0][0] += cell_a_angstrom / 10
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPeriodicBoxVectors(*(Vec3(*vector) * unit.nanometer for vector in box))
    context.setPositions(positions * unit.nanometer)
    state = context.getState(getPositions=True, getEnergy=True)
    del context, integrator
    (work / "state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")


def rewrite_whole_stage_cell(work: Path, lengths_angstrom: tuple[float, float, float]) -> None:
    """Represent a pressure-changed final cell in all four stage-specific artifacts."""
    box = tuple(Vec3(*(length if axis == vector else 0
                       for axis, length in enumerate(lengths_angstrom)))
                for vector in range(3)) * unit.angstrom
    system = XmlSerializer.deserialize((work / "system.xml").read_text(encoding="utf-8"))
    system.setDefaultPeriodicBoxVectors(*box)
    old_state = XmlSerializer.deserialize((work / "state.xml").read_text(encoding="utf-8"))
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPeriodicBoxVectors(*box)
    context.setPositions(old_state.getPositions())
    state = context.getState(getPositions=True, getEnergy=True)
    del context, integrator
    cif = PDBxFile(str(work / "stage.cif"))
    cif.topology.setPeriodicBoxVectors(box)
    with (work / "stage.cif").open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(cif.topology, state.getPositions(), stream, keepIds=True)
    sidecar_path = work / "topology.json"
    sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
    sidecar["boxVectorsAngstrom"] = [list(vector.value_in_unit(unit.angstrom)) for vector in box]
    sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
    (work / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
    (work / "state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")


class SmallStageExportWorkerCases(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="pim-export-small-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name) / "work"
        self.payload = write_small_stage(self.work)

    def verify(self, case: str):
        return invoke(self.work, self.payload, f"export-{case}")

    def assert_unmatched(self, case: str) -> dict:
        completed, terminal, _ = self.verify(case)
        self.assertEqual(0, completed.returncode, completed.stdout + completed.stderr)
        self.assertEqual("result", terminal["kind"])
        observed = terminal["payload"]["observations"]
        self.assertFalse(observed["readBackMatched"])
        self.assertTrue(observed["warnings"])
        return observed

    def assert_input_mismatch(self, case: str) -> None:
        completed, terminal, _ = self.verify(case)
        self.assertNotEqual(0, completed.returncode)
        self.assertEqual("error", terminal["kind"])
        self.assertEqual("inputMismatch", terminal["payload"]["failureCode"])

    def test_complete_selected_stage_is_read_back_with_one_correlated_result(self):
        completed, terminal, _ = self.verify("small-positive")
        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertEqual("result", terminal["kind"])
        result = terminal["payload"]
        self.assertEqual(("revision-1", "attempt-1", "minimized-stage-1"),
                         (result["studyRevisionId"], result["attemptId"], result["stageId"]))
        self.assertTrue(result["observations"]["readBackMatched"])
        self.assertEqual(2, result["observations"]["exportedAtomCount"])
        self.assertEqual(0, len(result["observations"]["warnings"]))
        paths = {item["role"]: Path(item["path"]) for item in result["artifacts"]}
        self.assertEqual({"stageMmcif", "stageBondGraph"}, set(paths))
        for item in result["artifacts"]:
            self.assertEqual(sha256(Path(item["path"])), item["sha256"])
        exported = PDBxFile(str(paths["stageMmcif"]))
        self.assertEqual(("C", "H"), tuple(atom.name for atom in exported.topology.atoms()))
        bond_graph = json.loads(paths["stageBondGraph"].read_text(encoding="utf-8"))
        self.assertEqual({(0, 1)}, {tuple(bond["atomIndices"]) for bond in bond_graph["bonds"]})
        self.assertEqual(2, len(exported.positions))

    def test_equilibrated_stage_with_changed_cell_in_all_four_artifacts_reads_back(self):
        previous = {name: sha256(self.work / name)
                    for name in ("stage.cif", "topology.json", "system.xml", "state.xml")}
        rewrite_whole_stage_cell(self.work, (42.0, 39.0, 43.0))
        self.payload = payload_for(self.work)
        self.payload["stageId"] = "equilibrated-stage-1"
        for name, old_digest in previous.items():
            self.assertNotEqual(old_digest, sha256(self.work / name), name)
        completed, terminal, _ = self.verify("equilibrated-new-cell")
        self.assertEqual(0, completed.returncode, completed.stdout + completed.stderr)
        self.assertEqual("result", terminal["kind"])
        result = terminal["payload"]
        self.assertEqual("equilibrated-stage-1", result["stageId"])
        self.assertTrue(result["observations"]["readBackMatched"])
        self.assertTrue(result["observations"]["cellMatched"])
        exported = PDBxFile(next(item["path"] for item in result["artifacts"]
                                 if item["role"] == "stageMmcif"))
        self.assertEqual((42.0, 39.0, 43.0),
                         tuple(round(length, 3) for length in
                               exported.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)))
        state = XmlSerializer.deserialize((self.work / "state.xml").read_text(encoding="utf-8"))
        system = XmlSerializer.deserialize((self.work / "system.xml").read_text(encoding="utf-8"))
        sidecar = json.loads((self.work / "topology.json").read_text(encoding="utf-8"))
        state_vectors = np.asarray(state.getPeriodicBoxVectors(asNumpy=True)
                                   .value_in_unit(unit.angstrom))
        np.testing.assert_allclose(sidecar["boxVectorsAngstrom"], state_vectors, atol=1e-8)
        np.testing.assert_allclose([vector.value_in_unit(unit.angstrom)
                                    for vector in system.getDefaultPeriodicBoxVectors()],
                                   state_vectors, atol=1e-8)
        np.testing.assert_allclose(exported.topology.getPeriodicBoxVectors()
                                   .value_in_unit(unit.angstrom), state_vectors, atol=0.01)

    def test_rehashed_equilibrated_sidecar_vector_shift_within_angle_tolerance_is_refused(self):
        rewrite_whole_stage_cell(self.work, (42.0, 39.0, 43.0))
        sidecar_path = self.work / "topology.json"
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        sidecar["boxVectorsAngstrom"][1][0] = 0.05
        sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("equilibrated-sidecar-vector-shift")
        self.assertFalse(observed["cellMatched"])
        self.assertTrue(any("bonded topology sidecar cell" in warning
                            for warning in observed["warnings"]))

    def test_rehashed_equilibrated_state_cell_shift_is_refused(self):
        rewrite_whole_stage_cell(self.work, (42.0, 39.0, 43.0))
        rewrite_state(self.work, cell_a_angstrom=0.02)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("equilibrated-state-cell-shift")
        self.assertFalse(observed["cellMatched"])

    def test_rehashed_equilibrated_system_default_cell_shift_is_refused(self):
        rewrite_whole_stage_cell(self.work, (42.0, 39.0, 43.0))
        system = XmlSerializer.deserialize((self.work / "system.xml").read_text(encoding="utf-8"))
        system.setDefaultPeriodicBoxVectors(Vec3(4.202, 0, 0) * unit.nanometer,
                                            Vec3(0, 3.9, 0) * unit.nanometer,
                                            Vec3(0, 0, 4.3) * unit.nanometer)
        (self.work / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("equilibrated-system-cell-shift")
        self.assertFalse(observed["cellMatched"])
        self.assertTrue(any("System default periodic cell" in warning
                            for warning in observed["warnings"]))

    def test_rehashed_selected_cif_coordinate_outside_tolerance_is_refused(self):
        rewrite_cif(self.work, first_x_angstrom=0.02)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("selected-coordinate-shift")
        self.assertGreater(observed["coordinateMaxDeviationAngstrom"], 0.01)
        self.assertIn("selected stage mmCIF coordinates", observed["warnings"][0])

    def test_rehashed_selected_cif_changes_within_declared_tolerances_are_permitted(self):
        rewrite_cif(self.work, first_x_angstrom=0.005, cell_a_angstrom=0.005,
                    cell_gamma_degrees=0.05)
        self.payload = payload_for(self.work)
        completed, terminal, _ = self.verify("selected-within-tolerance")
        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertEqual("result", terminal["kind"])
        self.assertTrue(terminal["payload"]["observations"]["readBackMatched"])
        self.assertLessEqual(terminal["payload"]["observations"]["coordinateMaxDeviationAngstrom"],
                             TOLERANCES["coordinateReadBackToleranceAngstrom"])

    def test_rehashed_selected_cif_cell_length_outside_tolerance_is_refused(self):
        rewrite_cif(self.work, cell_a_angstrom=0.02)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("selected-cell-length-shift")
        self.assertFalse(observed["cellMatched"])
        self.assertIn("selected stage mmCIF cell", observed["warnings"][0])

    def test_rehashed_selected_cif_cell_angle_outside_tolerance_is_refused(self):
        rewrite_cif(self.work, cell_gamma_degrees=0.2)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("selected-cell-angle-shift")
        self.assertFalse(observed["cellMatched"])
        self.assertIn("selected stage mmCIF cell", observed["warnings"][0])

    def test_rehashed_topology_sidecar_cell_outside_tolerance_is_refused(self):
        sidecar_path = self.work / "topology.json"
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        sidecar["boxVectorsAngstrom"][0][0] += 0.02
        sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("topology-sidecar-cell-shift")
        self.assertFalse(observed["cellMatched"])
        self.assertTrue(any("bonded topology sidecar cell" in warning for warning in observed["warnings"]))

    def test_rehashed_state_coordinate_outside_tolerance_is_refused(self):
        rewrite_state(self.work, first_x_angstrom=0.02)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("state-coordinate-shift")
        self.assertGreater(observed["coordinateMaxDeviationAngstrom"], 0.01)

    def test_rehashed_state_cell_outside_tolerance_is_refused(self):
        rewrite_state(self.work, cell_a_angstrom=0.02)
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("state-cell-shift")
        self.assertFalse(observed["cellMatched"])

    def test_changed_selected_stage_bytes_with_old_digest_are_refused(self):
        rewrite_cif(self.work, first_x_angstrom=0.001)
        completed, terminal, _ = self.verify("stale-selected-digest")
        self.assertNotEqual(0, completed.returncode)
        self.assertEqual("error", terminal["kind"])
        self.assertEqual("inputMismatch", terminal["payload"]["failureCode"])

    def test_rehashed_wrong_atom_identity_is_refused(self):
        text = (self.work / "stage.cif").read_text(encoding="utf-8")
        original = "ATOM      1 C   C    . GLY"
        self.assertIn(original, text)
        changed = text.replace(original, "ATOM      1 C   N    . GLY", 1)
        self.assertIn("GLY A    C     1", changed)
        changed = changed.replace("GLY A    C     1", "GLY A    N     1", 1)
        (self.work / "stage.cif").write_text(changed, encoding="utf-8")
        self.payload = payload_for(self.work)
        self.assert_input_mismatch("wrong-atom-identity")

    def test_rehashed_wrong_atom_order_is_refused(self):
        path = self.work / "stage.cif"
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        atom_lines = [index for index, line in enumerate(lines) if line.startswith("ATOM")]
        self.assertEqual(2, len(atom_lines))
        first, second = atom_lines
        lines[first], lines[second] = lines[second], lines[first]
        path.write_text("".join(lines), encoding="utf-8")
        self.payload = payload_for(self.work)
        self.assert_input_mismatch("wrong-atom-order")

    def test_rehashed_wrong_bond_is_refused(self):
        sidecar = json.loads((self.work / "topology.json").read_text(encoding="utf-8"))
        sidecar["bonds"] = []
        (self.work / "topology.json").write_text(json.dumps(sidecar), encoding="utf-8")
        self.payload = payload_for(self.work)
        observed = self.assert_unmatched("wrong-bond")
        self.assertFalse(observed["bondsMatched"])

    def test_rehashed_system_particle_count_mismatch_is_refused(self):
        system = XmlSerializer.deserialize((self.work / "system.xml").read_text(encoding="utf-8"))
        system.addParticle(12 * unit.dalton)
        (self.work / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
        self.payload = payload_for(self.work)
        self.assert_input_mismatch("system-particle-mismatch")


class SavedRealStageExportWorkerCase(unittest.TestCase):
    def test_full_85320_atom_stage_crosses_worker_and_independent_read_back(self):
        if not SAVED_HASHES.is_file() or not PYTHON.is_file():
            self.skipTest("Retained exact Slice 4 stage or scientific Python is unavailable")
        expected = json.loads(SAVED_HASHES.read_text(encoding="utf-8"))
        names = {"topologyCif": "minimized-coordinates.cif",
                 "topologyJson": "constructed-topology.json",
                 "systemXml": "constructed-system.xml", "stateXml": "minimized-state.xml"}
        expected_names = {"topologyCif": "minimizedCif", "topologyJson": "topologyJson",
                          "systemXml": "systemXml", "stateXml": "minimizedStateXml"}
        with tempfile.TemporaryDirectory(prefix="pim-export-real-") as temporary:
            work = Path(temporary)
            payload = {"studyRevisionId": "slice4-revision", "attemptId": "slice4-attempt",
                       "stageId": "slice4-minimized", **TOLERANCES}
            for role, name in names.items():
                source = SAVED / name
                self.assertEqual(expected["sha256ByRole"][expected_names[role]], sha256(source))
                target = work / name
                shutil.copyfile(source, target)
                payload[f"{role}Path"] = str(target)
                payload[f"{role}Sha256"] = sha256(target)
            completed, terminal, seconds = invoke(work, payload, "real-85320-export", timeout=600)
            self.assertEqual(0, completed.returncode, completed.stdout + completed.stderr)
            self.assertEqual("result", terminal["kind"])
            result = terminal["payload"]
            self.assertEqual(("slice4-revision", "slice4-attempt", "slice4-minimized"),
                             (result["studyRevisionId"], result["attemptId"], result["stageId"]))
            observed = result["observations"]
            self.assertTrue(observed["readBackMatched"], observed)
            self.assertEqual(85320, observed["sourceAtomCount"])
            self.assertEqual(85320, observed["exportedAtomCount"])
            self.assertLessEqual(observed["coordinateMaxDeviationAngstrom"], 0.01)
            artifacts = {item["role"]: item for item in result["artifacts"]}
            self.assertEqual({"stageMmcif", "stageBondGraph"}, set(artifacts))
            for item in artifacts.values():
                self.assertEqual(item["sha256"], sha256(Path(item["path"])))
            exported_path = Path(artifacts["stageMmcif"]["path"])
            gemmi_structure = gemmi.read_structure(str(exported_path))
            self.assertEqual(1, len(gemmi_structure))
            output = PDBxFile(str(exported_path))
            input_stage = PDBxFile(str(work / names["topologyCif"]))
            sidecar = json.loads(Path(artifacts["stageBondGraph"]["path"]).read_text(encoding="utf-8"))
            system = XmlSerializer.deserialize((work / names["systemXml"]).read_text(encoding="utf-8"))
            state = XmlSerializer.deserialize((work / names["stateXml"]).read_text(encoding="utf-8"))
            self.assertEqual(85320, system.getNumParticles())
            self.assertEqual(85320, len(sidecar["atoms"]))
            self.assertEqual(85320, len(list(output.topology.atoms())))
            identity = lambda topology: [(atom.residue.chain.id, atom.residue.id,
                                          atom.residue.name, atom.name, atom.element.symbol)
                                         for atom in topology.atoms()]
            self.assertEqual(identity(input_stage.topology), identity(output.topology))
            sidecar_identity = [
                (sidecar["chains"][sidecar["residues"][atom["residueIndex"]]["chainIndex"]]["id"],
                 sidecar["residues"][atom["residueIndex"]]["id"],
                 sidecar["residues"][atom["residueIndex"]]["name"], atom["name"], atom["element"])
                for atom in sidecar["atoms"]]
            self.assertEqual(sidecar_identity, identity(output.topology))
            atom_sites = gemmi.cif.read_file(str(exported_path)).sole_block().find_mmcif_category("_atom_site.")
            self.assertEqual(85320, len(atom_sites))
            self.assertEqual([atom[3] for atom in sidecar_identity],
                             [row["_atom_site.auth_atom_id"] for row in atom_sites])
            original_sidecar = json.loads((work / names["topologyJson"]).read_text(encoding="utf-8"))
            self.assertEqual({tuple(bond["atomIndices"]) for bond in original_sidecar["bonds"]},
                             {tuple(bond["atomIndices"]) for bond in sidecar["bonds"]})
            desired = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
            selected = np.asarray(input_stage.positions.value_in_unit(unit.angstrom))
            actual = np.asarray(output.positions.value_in_unit(unit.angstrom))
            self.assertLessEqual(float(np.max(np.linalg.norm(selected - desired, axis=1))), 0.01)
            self.assertLessEqual(float(np.max(np.linalg.norm(actual - desired, axis=1))), 0.01)
            state_cell = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
            selected_cell = np.asarray(input_stage.topology.getPeriodicBoxVectors().value_in_unit(unit.angstrom))
            output_cell = np.asarray(output.topology.getPeriodicBoxVectors().value_in_unit(unit.angstrom))
            self.assertLessEqual(float(np.max(np.abs(selected_cell - state_cell))), 0.01)
            self.assertLessEqual(float(np.max(np.abs(output_cell - state_cell))), 0.01)
            print(f"full saved-stage verify_export read-back: {seconds:.2f}s; "
                  f"{system.getNumParticles()} atoms; "
                  f"input CIF {payload['topologyCifSha256']}; State {payload['stateXmlSha256']}; "
                  f"output CIF {artifacts['stageMmcif']['sha256']}; "
                  f"bond graph {artifacts['stageBondGraph']['sha256']}; "
                  f"cell vectors A {state_cell.tolist()}; "
                  f"maximum deviation A {observed['coordinateMaxDeviationAngstrom']:.7g}; "
                  f"provider {result['provider']}")


if __name__ == "__main__":
    unittest.main()
