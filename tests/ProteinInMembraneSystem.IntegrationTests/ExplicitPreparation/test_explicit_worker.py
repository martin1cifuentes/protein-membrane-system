"""One-request native OpenMM construction and owner-local observation crossings.

The exact 6QWR source is prepared and oriented through the installed workers
before the native builder is invoked. No expected lipid count or cell is fed to
the builder; those are checked against the returned candidate from that call.
"""

from __future__ import annotations

import copy
from collections import defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src/ProteinInMembrane.Host"
CATALOGUE = ROOT / "config/policies/protein-membrane-slice4.json"
SOURCE = ROOT / "config/policies/source-assets/6QWR.pdb"
PYTHON = os.environ.get("PIM_TEST_PYTHON", str(ROOT / "out/python/bin/python"))
PPM = ROOT / "out/ppm2/immers"
PPM_LIBRARY = ROOT / "out/ppm2/res.lib"
PATCH = ROOT / "out/python/lib/python3.11/site-packages/openmm/app/data/DMPC.pdb"
PATCH_SHA256 = "19d8513a792ddf6f34edc16028da4410dc6d686945bede006ab1e5a27accb4ba"
PPM_SHA256 = "d58c7189f27b7e81150ba4f013cdc039aea9cd71c550e86b4796c5df9360f275"
LIBRARY_SHA256 = "26ef3b6ee3d237b0c29bcf549b721e33b0b157ccba198b728daa65dc2f817435"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def environment() -> dict[str, str]:
    return dict(os.environ, PYTHONPATH=str(SOURCE_ROOT), PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT),
                OPENMM_CPU_THREADS="2")


def invoke(work: Path, operation: str, payload: dict, request_id: str,
           timeout: int = 180) -> tuple[subprocess.CompletedProcess, list[dict], float]:
    request = {"requestId": request_id, "operation": operation,
               "workingDirectory": str(work), "payload": payload}
    started = time.monotonic()
    completed = subprocess.run([PYTHON, "-m", "ProteinInMembraneSystem.worker"],
                               input=json.dumps(request) + "\n", text=True,
                               capture_output=True, cwd=SOURCE_ROOT,
                               env=environment(), timeout=timeout, check=False)
    events = [json.loads(line) for line in completed.stdout.splitlines()]
    return completed, events, time.monotonic() - started


def terminal(completed: subprocess.CompletedProcess, events: list[dict], kind: str) -> dict:
    assert events, completed.stderr
    assert [event["kind"] for event in events].count("result") + [event["kind"] for event in events].count("error") == 1
    assert events[-1]["kind"] == kind, completed.stdout + completed.stderr
    return events[-1]["payload"]


def assets(result: dict) -> dict[str, Path]:
    return {item["role"]: Path(item["path"]) for item in result["artifacts"]}


def heavy_periodic_pair_minima(topology, positions, correspondence: dict,
                                cell: list[float], radius: float = 2.2) -> tuple[dict[str, float], dict[str, float]]:
    """Independent water/ion heavy scan; separately report image crossings."""
    from openmm import unit

    atoms = list(topology.atoms())
    coordinates = positions.value_in_unit(unit.angstrom)
    selected = [index for index, item in enumerate(correspondence["atoms"])
                if item["moleculeRole"] in {"water", "ion"} and atoms[index].element.symbol != "H"]
    roles = [correspondence["atoms"][index]["moleculeRole"] for index in selected]
    points = [tuple(float(coordinates[index][axis]) % cell[axis] for axis in range(3))
              for index in selected]
    divisions = [max(3, math.floor(length / radius)) for length in cell]
    bins: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    def bin_of(point):
        return tuple(min(divisions[axis] - 1,
                         math.floor(point[axis] / cell[axis] * divisions[axis]))
                     for axis in range(3))
    for index, point in enumerate(points):
        bins[bin_of(point)].append(index)
    periodic = {}
    seams = {}
    for index, point in enumerate(points):
        address = bin_of(point)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = tuple((address[axis] + delta) % divisions[axis]
                                     for axis, delta in enumerate((dx, dy, dz)))
                    for other in bins.get(neighbor, ()):
                        if other <= index:
                            continue
                        raw = [point[axis] - points[other][axis] for axis in range(3)]
                        delta = [raw[axis] - round(raw[axis] / cell[axis]) * cell[axis]
                                 for axis in range(3)]
                        distance = math.sqrt(sum(value * value for value in delta))
                        if distance > radius:
                            continue
                        key = "|".join(sorted((roles[index], roles[other])))
                        periodic[key] = min(periodic.get(key, math.inf), distance)
                        if any(abs(raw[axis] - delta[axis]) > 1e-6 for axis in range(3)):
                            seams[key] = min(seams.get(key, math.inf), distance)
    return periodic, seams


class NativeExplicitWorkerCrossing(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for path, expected in ((SOURCE, "f4c1503a60321c0cfe513e8211e43e20e2199aa5f71f22429e0e0ae8c97b0779"),
                               (PATCH, PATCH_SHA256), (PPM, PPM_SHA256),
                               (PPM_LIBRARY, LIBRARY_SHA256)):
            if not path.is_file() or not Path(PYTHON).is_file():
                raise unittest.SkipTest(f"Installed scientific asset is unavailable: {path}")
            if digest(path) != expected:
                raise AssertionError(f"Identified scientific bytes changed: {path}")
        cls.catalogue = json.loads(CATALOGUE.read_text())
        cls.policy = cls.catalogue["preparationPolicies"][0]
        cls.temporary = tempfile.TemporaryDirectory(prefix="pim-native-explicit-")
        cls.root = Path(cls.temporary.name)
        cls._prepare_and_orient()
        cls.real_result = None

    @classmethod
    def real_candidate(cls) -> None:
        if cls.real_result is not None:
            return
        cls.real_work, cls.real_payload = cls.new_work("real-native")
        cls.real_completed, cls.real_events, cls.real_seconds = invoke(
            cls.real_work, "construct_system", cls.real_payload, "native-6qwr-candidate", timeout=450)
        if cls.real_completed.returncode != 0:
            raise AssertionError(cls.real_completed.stdout + "\n" + cls.real_completed.stderr)
        cls.real_result = terminal(cls.real_completed, cls.real_events, "result")

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "temporary"):
            cls.temporary.cleanup()

    @classmethod
    def _stage(cls, work: Path, source: Path) -> Path:
        target = work / source.name
        shutil.copyfile(source, target)
        return target

    @classmethod
    def _stage_asset(cls, work: Path, raw: dict) -> dict:
        source = (CATALOGUE.parent / raw["path"]).resolve()
        if not source.is_file() or digest(source) != raw["sha256"]:
            raise AssertionError(f"Changed selected asset: {source}")
        return dict(raw, path=str(cls._stage(work, source)))

    @classmethod
    def _prepare_and_orient(cls) -> None:
        from gemmi import read_structure

        work = cls.root / "prepare"
        work.mkdir()
        source = cls._stage(work, SOURCE)
        protein_asset = cls._stage_asset(work, cls.catalogue["proteinChemicalStates"][0]["forceFieldFiles"][0])
        variants = []
        structure = read_structure(str(source))
        for residue in structure[0][0]:
            if residue.name in {"ASP", "GLU", "CYS", "HIS", "LYS"}:
                variants.append({"residue": {"model": 0, "chain": "A", "residue": residue.seqid.num,
                                              "insertionCode": "", "copyId": "A"},
                                 "variant": "HID" if residue.name == "HIS" else residue.name})
        oxt = {"model": 0, "chain": "A", "residue": 211, "insertionCode": "", "copyId": "A"}
        preparation = {"studyRevisionId": "native-6qwr-revision", "nominalPh": 7.0,
                       "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                       "sourcePath": str(source), "sourceSha256": digest(source),
                       "modelIndex": 0, "assemblyId": None,
                       "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                       "altlocChoices": [], "retainedPartners": [],
                       "approvedHeavyAtoms": [{"residue": oxt, "atomName": "OXT",
                                               "decisionId": "native-6qwr-oxt-approval"}],
                       "approvedDisulfides": [], "residueVariants": variants,
                       "forceFieldFiles": [protein_asset],
                       "geometrySpec": cls.catalogue["proteinStructuralPolicies"][0]["measurement"]}
        completed, events, _ = invoke(work, "prepare_protein", preparation, "native-6qwr-preparation")
        if completed.returncode != 0:
            raise AssertionError(completed.stdout + "\n" + completed.stderr)
        prepared = assets(terminal(completed, events, "result"))
        cls.prepared = prepared
        cls.prepared_mapping = json.loads(prepared["correspondenceJson"].read_text())

        oriented_work = cls.root / "orient"
        oriented_work.mkdir()
        staged_prepared = cls._stage(oriented_work, prepared["preparedPdb"])
        placement = {"studyRevisionId": "native-6qwr-revision", "preparedProteinId": "native-6qwr-protein",
                     "preparedPdbPath": str(staged_prepared), "preparedSha256": digest(staged_prepared),
                     "ppmExecutablePath": str(PPM),
                     "ppmVersion": "2.0 GitLab 0de704acd10fbf7f4fbafb9bc925e01b06f38ceb",
                     "ppmExecutableSha256": PPM_SHA256,
                     "ppmResidueLibraryPath": str(PPM_LIBRARY),
                     "ppmResidueLibrarySha256": LIBRARY_SHA256,
                     "topologyKind": "membrane-spanning", "ppmNterminalSide": "in"}
        completed, events, _ = invoke(oriented_work, "place_ppm", placement,
                                      "native-6qwr-orientation")
        if completed.returncode != 0:
            raise AssertionError(completed.stdout + "\n" + completed.stderr)
        cls.oriented = assets(terminal(completed, events, "result"))["orientedPdb"]

    @classmethod
    def new_work(cls, name: str) -> tuple[Path, dict]:
        work = cls.root / name
        work.mkdir()
        prepared = cls._stage(work, cls.prepared["preparedPdb"])
        oriented = cls._stage(work, cls.oriented)
        graph = cls._stage(work, cls.prepared["preparedBondGraph"])
        force_fields = [cls._stage_asset(work, item) for item in cls.policy["forceFieldFiles"]]
        def representation(raw: dict) -> dict:
            item = dict(raw)
            for key in ("templatePath", "coordinateTemplatePath"):
                item[key] = str(cls._stage(work, (CATALOGUE.parent / raw[key]).resolve()))
            return item
        lipid = next(item for item in cls.catalogue["lipids"] if item["speciesId"] == "DMPC")
        construction = cls.policy["construction"]
        return work, {
            "studyRevisionId": "native-6qwr-revision", "attemptId": "native-6qwr-attempt",
            "orientedPdbPath": str(oriented), "orientedPdbSha256": digest(oriented),
            "preparedPdbPath": str(prepared), "preparedPdbSha256": digest(prepared),
            "preparedCorrespondence": copy.deepcopy(cls.prepared_mapping),
            "preparedBondGraphPath": str(graph), "preparedBondGraphSha256": digest(graph),
            "lipid": representation(lipid), "water": representation(cls.policy["water"]),
            "sodium": representation(cls.policy["sodium"]),
            "chloride": representation(cls.policy["chloride"]),
            "nativePatchPath": str(PATCH), "nativePatchSha256": PATCH_SHA256,
            "providerName": construction["providerName"],
            "providerVersion": construction["providerVersion"],
            "lipidTypeArgument": "DMPC", "positiveIonArgument": "Na+",
            "negativeIonArgument": "Cl-", "membraneCenterZNanometers": 0.0,
            "minimumPaddingNanometers": construction["minimumPaddingNanometers"],
            "ionicStrengthMolar": 0.15,
            "forceFieldFiles": force_fields,
            "systemSettings": copy.deepcopy(cls.policy["systemSettings"]),
            "localObservationSpec": copy.deepcopy(cls.policy["localStateObservation"]),
            "maximumAtomCount": construction["maximumAtomCount"],
            "maximumCellDimensionAngstrom": construction["maximumCellDimensionAngstrom"]}

    def test_real_native_candidate_is_correlated_complete_and_same_invocation(self) -> None:
        self.real_candidate()
        from openmm import XmlSerializer, unit
        from openmm.app import PDBxFile
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.worker.parameterized_structure import (
            _full_atom_sequence, _read_topology_data, _system_bonds_match)

        self.assertEqual(0, self.real_completed.returncode)
        self.assertEqual({"native-6qwr-candidate"}, {event["requestId"] for event in self.real_events})
        self.assertEqual(["started", "nativeMembraneStarted", "nativeMembraneReturned",
                          "constructedCoordinatesObserved"],
                         [event["payload"]["stage"] for event in self.real_events
                          if event["kind"] == "progress"])
        result = self.real_result
        observed = result["observations"]
        self.assertEqual({"name": "OpenMM Modeller.addMembrane",
                          "version": "8.6.0.dev-c6173db"}, result["provider"])
        self.assertEqual(PATCH_SHA256, observed["nativePatchSha256"])
        self.assertEqual(3205, len(self.prepared_mapping["atoms"]))
        self.assertGreater(observed["atomCount"], 3205)
        self.assertLessEqual(observed["atomCount"], self.policy["construction"]["maximumAtomCount"])
        self.assertEqual(observed["atomCount"], observed["correspondedResultAtomCount"])
        self.assertTrue(observed["proteinIdentityAndBondsPreserved"])
        self.assertLessEqual(observed["maximumProteinCoordinateDeviationAngstrom"], 1e-6)
        self.assertEqual(3, len(observed["proteinPeriodicImageGapsAngstrom"]))
        self.assertTrue(all(value >= 20.0 - 1e-4 for value in observed["proteinPeriodicImageGapsAngstrom"]))
        self.assertAlmostEqual(-8, observed["proteinNetChargeElementary"], places=4)
        self.assertAlmostEqual(0, observed["netChargeElementary"], places=4)
        self.assertTrue(math.isfinite(observed["initialPotentialEnergyKjMol"]))
        counts = {(item["role"], item["physicalSide"], item["speciesId"]): item["count"]
                  for item in observed["speciesCounts"]}
        self.assertEqual(counts["lipid", "upper", "DMPC"], counts["lipid", "lower", "DMPC"])
        self.assertGreater(counts["lipid", "upper", "DMPC"], 0)
        self.assertEqual(observed["waterCount"], sum(value for (role, _, _), value in counts.items()
                                                      if role == "water"))
        self.assertEqual(observed["positiveIonCount"], sum(value for (role, _, _), value in counts.items()
                                                            if role == "positiveIon"))
        self.assertEqual(observed["negativeIonCount"], sum(value for (role, _, _), value in counts.items()
                                                            if role == "negativeIon"))
        neutralizers = 8
        salt_pairs = observed["negativeIonCount"]
        equivalent = observed["waterCount"] + observed["positiveIonCount"] + observed["negativeIonCount"]
        self.assertEqual(neutralizers + salt_pairs, observed["positiveIonCount"])
        self.assertEqual(math.floor(0.5 + (equivalent - neutralizers) * 0.15 / 55.4), salt_pairs)
        self.assertEqual(["Observed"], [observed["localState"]["standing"]])
        measurements = {item["name"]: item["value"] for item in observed["localState"]["measurements"]}
        self.assertLess(measurements["minimumIntermolecularDistanceAngstrom"], 1.5)
        self.assertGreaterEqual(measurements["minimumIntermolecularHeavyAtomDistanceAngstrom"], 1.5)
        self.assertGreater(measurements["leafletHeadSeparationAngstrom"], 24.0)
        self.assertLess(measurements["leafletHeadSeparationAngstrom"], 46.0)
        self.assertLess(abs(measurements["proteinBilayerMidplaneOffsetAngstrom"]), 20.0)
        role_pairs = {(item["firstMoleculeRole"], item["secondMoleculeRole"]): item
                      for item in observed["localState"]["rolePairMeasurements"]}
        self.assertEqual({("protein", "lipid"), ("water", "water"),
                          ("water", "ion"), ("ion", "ion")}, set(role_pairs))
        self.assertGreater(role_pairs["protein", "lipid"]["pairsWithinSearchRadius"], 0)
        self.assertLess(role_pairs["protein", "lipid"]["minimumDistanceAngstrom"], 6.0)
        self.assertGreater(role_pairs["water", "water"]["pairsWithinSearchRadius"], 0)
        self.assertGreater(role_pairs["water", "ion"]["pairsWithinSearchRadius"], 0)
        self.assertGreater(role_pairs["ion", "ion"]["pairsWithinSearchRadius"], 0)

        artifacts = {item["role"]: item for item in result["artifacts"]}
        self.assertEqual({"topologyCif", "topologyJson", "systemXml", "stateXml",
                          "correspondenceJson"}, set(artifacts))
        for item in artifacts.values():
            self.assertEqual(item["sha256"], digest(Path(item["path"])))
        cif = PDBxFile(artifacts["topologyCif"]["path"])
        topology = _read_topology_data(Path(artifacts["topologyJson"]["path"]))
        system = XmlSerializer.deserialize(Path(artifacts["systemXml"]["path"]).read_text())
        state = XmlSerializer.deserialize(Path(artifacts["stateXml"]["path"]).read_text())
        self.assertEqual(observed["atomCount"], system.getNumParticles())
        self.assertEqual(observed["atomCount"], len(state.getPositions()))
        self.assertEqual(_full_atom_sequence(cif.topology), _full_atom_sequence(topology))
        self.assertTrue(_system_bonds_match(topology, system))
        actual_cell = [vector[i].value_in_unit(unit.angstrom)
                       for i, vector in enumerate(topology.getPeriodicBoxVectors())]
        for expected, actual in zip(observed["actualCellAngstrom"], actual_cell):
            self.assertAlmostEqual(expected, actual, places=6)
        correspondence = json.loads(Path(artifacts["correspondenceJson"]["path"]).read_text())
        self.assertTrue(correspondence["complete"])
        self.assertEqual(observed["atomCount"], len(correspondence["atoms"]))
        self.assertEqual(digest(Path(artifacts["topologyCif"]["path"])), correspondence["resultId"])
        self.assertTrue(all(item["resultAtomId"].startswith("protein:")
                            for item in correspondence["atoms"][:3205]))
        self.assertTrue(all(item["resultAtomId"].startswith("component:")
                            for item in correspondence["atoms"][3205:]))
        periodic_heavy, seams = heavy_periodic_pair_minima(
            cif.topology, cif.positions, correspondence, observed["actualCellAngstrom"])
        self.assertIn("water|water", periodic_heavy)
        self.assertLess(periodic_heavy["water|water"], 2.2)
        self.assertGreaterEqual(periodic_heavy["water|water"] + 1e-3,
                                measurements["minimumIntermolecularHeavyAtomDistanceAngstrom"])
        print(f"native 6QWR construction: {self.real_seconds:.2f}s, {observed['atomCount']} atoms; "
              f"heavy periodic minima {periodic_heavy}; image crossings {seams}", flush=True)

    def test_changed_provider_source_and_chemistry_are_distinct_pre_call_refusals(self) -> None:
        cases = [
            ("patch-digest", lambda p: p.update(nativePatchSha256="0" * 64), "inputMismatch"),
            ("patch-path", lambda p: p.update(nativePatchPath=p["preparedPdbPath"]), "providerMismatch"),
            ("provider-version", lambda p: p.update(providerVersion="8.5.0"), "providerMismatch"),
            ("other-lipid", lambda p: p.update(lipidTypeArgument="DOPC"), "unsupportedPolicy"),
            ("oriented-digest", lambda p: p.update(orientedPdbSha256="0" * 64), "inputMismatch"),
            ("prepared-digest", lambda p: p.update(preparedPdbSha256="0" * 64), "inputMismatch"),
            ("graph-digest", lambda p: p.update(preparedBondGraphSha256="0" * 64), "inputMismatch"),
            ("lipid-template", lambda p: p["lipid"].update(coordinateTemplateSha256="0" * 64), "inputMismatch"),
            ("force-field", lambda p: p["forceFieldFiles"][1].update(sha256="0" * 64), "inputMismatch"),
            ("insufficient-atom-cap", lambda p: p.update(maximumAtomCount=3204), "resourceRefused"),
            ("wrong-sodium-role", lambda p: p["sodium"].update(speciesId="CL"), "unsupportedPolicy"),
        ]
        for name, change, code in cases:
            with self.subTest(name=name):
                work, payload = self.new_work("refuse-" + name)
                change(payload)
                completed, events, _ = invoke(work, "construct_system", payload, name)
                self.assertNotEqual(0, completed.returncode)
                self.assertEqual(code, terminal(completed, events, "error")["failureCode"])
                self.assertFalse((work / "constructed-topology.cif").exists())
                self.assertFalse((work / "constructed-system.xml").exists())

    def test_identified_patch_maps_to_canonical_dmpc_graph_and_stereo(self) -> None:
        from openmm.app import PDBFile, PDBxFile
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import (
            _native_lipid_side, _native_molecule_bonds, _native_molecule_matches)
        from ProteinInMembraneSystem.worker.exchange import WorkError

        patch = PDBFile(str(PATCH))
        reference = PDBxFile(str((CATALOGUE.parent / next(item for item in self.catalogue["lipids"]
                   if item["speciesId"] == "DMPC")["coordinateTemplatePath"]).resolve()))
        molecule = next(item for item in patch.topology.residues() if item.name == "DMP")
        canonical = next(reference.topology.residues())
        stereo = next(item for item in self.catalogue["lipids"] if item["speciesId"] == "DMPC")["stereoChecks"]
        atoms = _native_molecule_matches(molecule, reference, "DMPC", patch.positions, stereo,
                                         _native_molecule_bonds(patch.topology, molecule),
                                         _native_molecule_bonds(reference.topology, canonical))
        self.assertEqual(118, len(atoms))
        self.assertEqual("4C21", atoms[77].name)
        from openmm import unit
        points = [tuple(float(value) for value in position.value_in_unit(unit.angstrom))
                  for position in patch.positions]
        side = _native_lipid_side(atoms, points, 0.0, "DMPC")
        self.assertIn(side, {"upper", "lower"})
        inverted_points = list(points)
        head = next(atom for atom in atoms if atom.name == "P")
        wrong_tail_z = points[head.index][2] + (10 if side == "upper" else -10)
        for atom in atoms:
            if atom.name in {"4C21", "4C31"}:
                x, y, _ = inverted_points[atom.index]
                inverted_points[atom.index] = (x, y, wrong_tail_z)
        with self.assertRaises(WorkError) as leaflet_refusal:
            _native_lipid_side(atoms, inverted_points, 0.0, "DMPC")
        self.assertEqual("providerMismatch", leaflet_refusal.exception.code)
        inverted = copy.deepcopy(stereo)
        inverted[0]["expected"] = "positive"
        with self.assertRaises(WorkError) as refusal:
            _native_molecule_matches(molecule, reference, "DMPC", patch.positions, inverted,
                                     _native_molecule_bonds(patch.topology, molecule),
                                     _native_molecule_bonds(reference.topology, canonical))
        self.assertEqual("providerMismatch", refusal.exception.code)

    def test_changed_candidate_bytes_refuse_stage_load_before_minimization(self) -> None:
        self.real_candidate()
        real = {item["role"]: item for item in self.real_result["artifacts"]}
        for role, field in (("topologyCif", "topologyCif"),
                            ("topologyJson", "topologyJson"),
                            ("systemXml", "systemXml"),
                            ("stateXml", "stateXml")):
            with self.subTest(role=role):
                work = self.root / ("changed-stage-" + role)
                work.mkdir()
                payload = {"studyRevisionId": "native-6qwr-revision",
                           "attemptId": "native-6qwr-attempt", "stageId": "candidate-stage",
                           "maxIterations": 1, "rmsForceTargetKjMolNm": 10.0}
                for artifact_role in ("topologyCif", "topologyJson", "systemXml", "stateXml"):
                    source = Path(real[artifact_role]["path"])
                    copied = self._stage(work, source)
                    payload[artifact_role + "Path"] = str(copied)
                    payload[artifact_role + "Sha256"] = real[artifact_role]["sha256"]
                changed = Path(payload[field + "Path"])
                with changed.open("ab") as stream:
                    stream.write(b"\nchanged-after-review\n")
                completed, events, _ = invoke(work, "minimize", payload,
                                              "changed-candidate-" + role, timeout=35)
                self.assertNotEqual(0, completed.returncode)
                self.assertEqual("inputMismatch", terminal(completed, events, "error")["failureCode"])
                self.assertFalse((work / "minimized-state.xml").exists())

    def test_postcall_changed_parameter_and_patch_bytes_refuse_candidate(self) -> None:
        from openmm.app import Modeller
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.worker import construction
        from ProteinInMembraneSystem.worker.exchange import WorkError

        work, payload = self.new_work("postcall-force-field")
        changed_asset = Path(payload["forceFieldFiles"][1]["path"])
        def changed_parameters(_modeller, *args, **kwargs):
            with changed_asset.open("ab") as stream:
                stream.write(b"\nchanged-during-native-call\n")
        with mock.patch.object(Modeller, "addMembrane", autospec=True,
                               side_effect=changed_parameters):
            with self.assertRaises(WorkError) as refused:
                construction.construct_system(work, payload, lambda _stage, _detail: None)
        self.assertEqual("inputMismatch", refused.exception.code)
        self.assertFalse((work / "constructed-topology.cif").exists())

        work, payload = self.new_work("postcall-patch")
        original_check = construction.verify_sha256
        patch_checks = 0
        def changed_patch(path, expected, name):
            nonlocal patch_checks
            if Path(path).resolve() == PATCH.resolve():
                patch_checks += 1
                if patch_checks == 2:
                    raise WorkError("inputMismatch", "Installed native patch changed during provider call")
            return original_check(path, expected, name)
        with mock.patch.object(Modeller, "addMembrane", autospec=True,
                               side_effect=lambda _modeller, *args, **kwargs: None), \
             mock.patch.object(construction, "verify_sha256", side_effect=changed_patch):
            with self.assertRaises(WorkError) as refused:
                construction.construct_system(work, payload, lambda _stage, _detail: None)
        self.assertEqual(2, patch_checks)
        self.assertEqual("inputMismatch", refused.exception.code)
        self.assertFalse((work / "constructed-topology.cif").exists())

    def test_malformed_native_identity_and_cell_refuse_before_candidate(self) -> None:
        from openmm import Vec3, unit
        from openmm.app import Modeller, Topology, element
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.worker import construction
        from ProteinInMembraneSystem.worker.exchange import WorkError

        def changed_protein_coordinate(modeller, *_args, **_kwargs):
            modeller.positions[0] = modeller.positions[0] + Vec3(0.00001, 0, 0) * unit.nanometer

        def changed_protein_bond(modeller, *_args, **_kwargs):
            atoms = list(modeller.topology.atoms())
            modeller.topology.addBond(atoms[0], atoms[100])

        def missing_cell(modeller, *_args, **_kwargs):
            modeller.topology.setPeriodicBoxVectors(None)

        def insufficient_cell(modeller, *_args, **_kwargs):
            modeller.topology.setPeriodicBoxVectors((Vec3(12, 0, 0), Vec3(0, 10, 0),
                                                     Vec3(0, 0, 1)) * unit.nanometer)

        def unidentified_species(modeller, *_args, **_kwargs):
            extra = Topology()
            residue = extra.addResidue("UNQ", extra.addChain("U"), id="1")
            extra.addAtom("O", element.oxygen, residue)
            modeller.add(extra, [Vec3(0, 0, 0)] * unit.nanometer)

        def wrong_ion_identity(modeller, *_args, **_kwargs):
            extra = Topology()
            residue = extra.addResidue("NA", extra.addChain("I"), id="1")
            extra.addAtom("O", element.oxygen, residue)
            modeller.add(extra, [Vec3(0, 0, 0)] * unit.nanometer)

        cases = [
            ("protein-coordinate", changed_protein_coordinate, "providerMismatch", "moved"),
            ("protein-bond", changed_protein_bond, "providerMismatch", "bonds"),
            ("missing-cell", missing_cell, "invalidGeometry", "periodic cell"),
            ("insufficient-cell", insufficient_cell, "invalidGeometry", "image separation"),
            ("unidentified-species", unidentified_species, "providerMismatch", "unsupported residue"),
            ("wrong-ion-identity", wrong_ion_identity, "providerMismatch", "Native NA atom order"),
        ]
        for name, change, code, reason in cases:
            with self.subTest(name=name):
                work, payload = self.new_work("native-malformed-" + name)
                def changed_output(modeller, *args, **kwargs):
                    modeller.topology.setPeriodicBoxVectors((Vec3(12, 0, 0), Vec3(0, 10, 0),
                                                             Vec3(0, 0, 12)) * unit.nanometer)
                    change(modeller, *args, **kwargs)
                with mock.patch.object(Modeller, "addMembrane", autospec=True,
                                       side_effect=changed_output):
                    with self.assertRaises(WorkError) as refused:
                        construction.construct_system(work, payload, lambda _stage, _detail: None)
                self.assertEqual(code, refused.exception.code)
                self.assertIn(reason, refused.exception.message)
                self.assertFalse((work / "constructed-topology.cif").exists())
                self.assertFalse((work / "constructed-system.xml").exists())

    def test_native_wrong_leaflet_orientation_and_parameter_failure_refuse(self) -> None:
        from openmm import Vec3, unit
        from openmm.app import ForceField, Modeller, PDBFile, Topology
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.worker import construction
        from ProteinInMembraneSystem.worker.exchange import WorkError

        patch = PDBFile(str(PATCH))
        patch_residues = list(patch.topology.residues())
        patch_bonds = list(patch.topology.bonds())

        def append_lipid(modeller, source_index: int, inverted_tail: bool = False) -> None:
            source = patch_residues[source_index]
            extra = Topology()
            residue = extra.addResidue("DMP", extra.addChain(str(source_index)), id="1")
            copied = {atom: extra.addAtom(atom.name, atom.element, residue)
                      for atom in source.atoms()}
            for first, second in patch_bonds:
                if first in copied and second in copied:
                    extra.addBond(copied[first], copied[second])
            positions = []
            head_z = next(patch.positions[atom.index].value_in_unit(unit.nanometer)[2]
                          for atom in source.atoms() if atom.name == "P") - 3.0
            for atom in source.atoms():
                point = patch.positions[atom.index] - Vec3(0, 0, 3) * unit.nanometer
                if inverted_tail and atom.name in {"4C21", "4C31"}:
                    point = Vec3(point[0].value_in_unit(unit.nanometer),
                                 point[1].value_in_unit(unit.nanometer),
                                 head_z + 1.0) * unit.nanometer
                positions.append(point)
            modeller.add(extra, positions)

        def native_output(modeller, *_args, inverted_tail=False, **_kwargs):
            append_lipid(modeller, 0, inverted_tail)
            if not inverted_tail:
                append_lipid(modeller, 64)
            modeller.topology.setPeriodicBoxVectors((Vec3(12, 0, 0), Vec3(0, 10, 0),
                                                     Vec3(0, 0, 12)) * unit.nanometer)

        work, payload = self.new_work("native-wrong-leaflet")
        with mock.patch.object(Modeller, "addMembrane", autospec=True,
                               side_effect=lambda modeller, *args, **kwargs:
                               native_output(modeller, *args, inverted_tail=True, **kwargs)):
            with self.assertRaises(WorkError) as refused:
                construction.construct_system(work, payload, lambda _stage, _detail: None)
        self.assertEqual("providerMismatch", refused.exception.code)
        self.assertIn("wrong leaflet head-to-tail orientation", refused.exception.message)
        self.assertFalse((work / "constructed-topology.cif").exists())

        work, payload = self.new_work("native-unparameterized")
        original_create_system = ForceField.createSystem
        def fail_generated_parameters(force_field, topology, *args, **kwargs):
            if topology.getNumAtoms() > len(self.prepared_mapping["atoms"]):
                raise ValueError("controlled generated molecule has no combined parameters")
            return original_create_system(force_field, topology, *args, **kwargs)
        with mock.patch.object(Modeller, "addMembrane", autospec=True, side_effect=native_output), \
             mock.patch.object(ForceField, "createSystem", autospec=True,
                               side_effect=fail_generated_parameters):
            with self.assertRaisesRegex(ValueError, "no combined parameters"):
                construction.construct_system(work, payload, lambda _stage, _detail: None)
        self.assertFalse((work / "constructed-topology.cif").exists())
        self.assertFalse((work / "constructed-system.xml").exists())

    def test_interrupted_native_provider_has_no_constructed_result(self) -> None:
        work, payload = self.new_work("interrupted-native")
        request = {"requestId": "native-interrupted", "operation": "construct_system",
                   "workingDirectory": str(work), "payload": payload}
        process = subprocess.Popen([PYTHON, "-m", "ProteinInMembraneSystem.worker"],
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True,
                                   cwd=SOURCE_ROOT, env=environment(), start_new_session=True)
        events = []
        try:
            assert process.stdin is not None and process.stdout is not None
            process.stdin.write(json.dumps(request) + "\n")
            process.stdin.flush()
            deadline = time.monotonic() + 35
            while time.monotonic() < deadline:
                readable, _, _ = select.select([process.stdout], [], [], 1)
                if not readable:
                    continue
                line = process.stdout.readline()
                if not line:
                    break
                event = json.loads(line)
                events.append(event)
                if event["kind"] == "progress" and event["payload"]["stage"] == "nativeMembraneStarted":
                    break
            self.assertTrue(any(event["kind"] == "progress" and
                                event["payload"]["stage"] == "nativeMembraneStarted" for event in events), events)
        finally:
            os.killpg(process.pid, signal.SIGKILL)
            process.communicate(timeout=10)
        self.assertFalse(any(event["kind"] == "result" for event in events))
        self.assertFalse((work / "constructed-topology.cif").exists())
        self.assertFalse((work / "constructed-system.xml").exists())


class NativeObservationMechanics(unittest.TestCase):
    def test_heavy_contact_metric_excludes_short_hydrogen_pair(self) -> None:
        from openmm import Vec3, unit
        from openmm.app import Topology, element
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.worker.local_state_observations import observe_local_state

        topology = Topology()
        atoms = []
        for index, symbol in enumerate(("C", "H", "O", "H")):
            chain = topology.addChain(str(index // 2))
            residue = topology.addResidue("MOL", chain, id="1")
            atoms.append(topology.addAtom(symbol, element.get_by_symbol(symbol), residue))
        topology.addBond(atoms[0], atoms[1])
        topology.addBond(atoms[2], atoms[3])
        mapped = {"complete": True, "atoms": [
            {"resultAtomIndex": index, "moleculeRole": "protein" if index < 2 else "water",
             "atomRole": "body", "physicalSide": None, "element": symbol}
            for index, symbol in enumerate(("C", "H", "O", "H"))]}
        positions = [Vec3(x, 0, 0) for x in (0.0, 1.0, 3.0, 1.3)] * unit.angstrom
        spec = {"requiredMetricNames": ["minimumIntermolecularDistanceAngstrom",
                 "minimumIntermolecularHeavyAtomDistanceAngstrom"],
                "contactSearchRadiusAngstrom": 4.0, "maximumReportedPairs": 4,
                "usePeriodicBoundary": False,
                "contactRolePairs": [{"firstMoleculeRole": "protein", "secondMoleculeRole": "water"}],
                "atomRadiusByElementAngstrom": {"C": 1.7, "H": 1.2, "O": 1.52}}
        observed = observe_local_state(topology, positions, mapped, spec)
        measurements = {item["name"]: item["value"] for item in observed["measurements"]}
        self.assertEqual("Observed", observed["standing"])
        self.assertAlmostEqual(0.3, measurements["minimumIntermolecularDistanceAngstrom"])
        self.assertAlmostEqual(3.0, measurements["minimumIntermolecularHeavyAtomDistanceAngstrom"])
        unavailable = observe_local_state(topology, positions, mapped,
                                          dict(spec, contactSearchRadiusAngstrom=2.0))
        self.assertEqual("Unavailable", unavailable["standing"])
        self.assertIn("heavy-atom pair", unavailable["unavailableReason"])

    def test_constraint_tangent_force_excludes_bond_normal_reaction(self) -> None:
        from openmm import Context, CustomExternalForce, HarmonicBondForce, System, Vec3, VerletIntegrator, unit
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.ExplicitPreparation.Minimization.worker.minimization import _constrained_tangent_force

        system = System()
        system.addParticle(12.0)
        system.addParticle(1.0)
        system.addConstraint(0, 1, 1.0)
        bond = HarmonicBondForce()
        bond.addBond(0, 1, 0.8, 1000.0)
        system.addForce(bond)
        transverse = CustomExternalForce("10*y")
        transverse.addParticle(0, [])
        system.addForce(transverse)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator)
        context.setPositions([Vec3(0, 0, 0), Vec3(1, 0, 0)] * unit.nanometer)
        tangent, raw, error = _constrained_tangent_force(
            system, context.getState(getPositions=True, getForces=True))
        self.assertGreater(raw, 100.0)
        self.assertAlmostEqual(math.sqrt(50), tangent, places=6)
        self.assertLess(error, 1e-10)
        del context, integrator

    def test_full_system_mmcif_representation_crosses_pdb_serial_capacity(self) -> None:
        from openmm import Vec3, unit
        from openmm.app import PDBxFile, Topology, element
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.worker.parameterized_structure import (
            _full_atom_sequence, _read_topology_data, _topology_data)

        count = 100_001
        topology = Topology()
        chain = topology.addChain("A")
        for index in range(count):
            residue = topology.addResidue("HOH", chain, id=str(index + 1))
            topology.addAtom("O", element.oxygen, residue)
        topology.setPeriodicBoxVectors((Vec3(20, 0, 0), Vec3(0, 20, 0),
                                        Vec3(0, 0, 20)) * unit.nanometer)
        positions = [Vec3(1, 1, 1) for _ in range(count)] * unit.nanometer
        with tempfile.TemporaryDirectory(prefix="pim-native-mmcif-capacity-") as temporary:
            path = Path(temporary) / "full-system.cif"
            with path.open("w", encoding="utf-8") as stream:
                PDBxFile.writeFile(topology, positions, stream, keepIds=True)
            sidecar = Path(temporary) / "full-system-topology.json"
            sidecar.write_text(json.dumps(_topology_data(topology), separators=(",", ":")))
            reread = PDBxFile(str(path))
            reconstructed = _read_topology_data(sidecar)
            self.assertEqual(count, len(list(reread.topology.atoms())))
            self.assertEqual("100001", list(reread.topology.residues())[-1].id)
            self.assertEqual(_full_atom_sequence(reread.topology),
                             _full_atom_sequence(reconstructed))
            self.assertEqual(count, len(list(reconstructed.atoms())))


if __name__ == "__main__":
    unittest.main()
