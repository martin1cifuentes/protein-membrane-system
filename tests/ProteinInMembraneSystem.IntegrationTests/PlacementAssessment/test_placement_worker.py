"""Slice-3 one-request placement worker and real local PPM crossing."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src" / "ProteinInMembrane.Host"
FIXTURE = Path(__file__).parent / "fixtures" / "6qwr-prepared.pdb"
PYTHON = os.environ.get("PIM_TEST_PYTHON", sys.executable)
PPM = Path(os.environ.get("PIM_PPM_EXECUTABLE", ROOT / "out" / "ppm2" / "immers"))
RESIDUE_LIBRARY = Path(os.environ.get("PIM_PPM_RESIDUE_LIBRARY", ROOT / "out" / "ppm2" / "res.lib"))
FIXTURE_SHA256 = "80e2dcade32491cd49f86750ea101b050899d900ae72190f5b97f79b2052abc6"
PPM_SHA256 = "d58c7189f27b7e81150ba4f013cdc039aea9cd71c550e86b4796c5df9360f275"
LIBRARY_SHA256 = "26ef3b6ee3d237b0c29bcf549b721e33b0b157ccba198b728daa65dc2f817435"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() != "DUM"]


def signature(line: str) -> tuple[str, ...]:
    return (line[:30], line[54:])


def xyz(lines: list[str]) -> np.ndarray:
    return np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                     for line in lines])


def payload(work: Path, *, executable: Path = PPM, library: Path = RESIDUE_LIBRARY) -> dict:
    prepared = work / "prepared.pdb"
    shutil.copyfile(FIXTURE, prepared)
    return {
        "studyRevisionId": "6qwr-revision", "preparedProteinId": "6qwr-prepared",
        "preparedPdbPath": str(prepared), "preparedSha256": digest(prepared),
        "ppmExecutablePath": str(executable), "ppmVersion": "2.0 GitLab 0de704acd10fbf7f4fbafb9bc925e01b06f38ceb",
        "ppmExecutableSha256": digest(executable) if executable.is_file() else PPM_SHA256,
        "ppmResidueLibraryPath": str(library),
        "ppmResidueLibrarySha256": digest(library) if library.is_file() else LIBRARY_SHA256,
        "topologyKind": "membrane-spanning", "ppmNterminalSide": "in",
    }


def invoke(work: Path, operation: str, data: dict, *, request_id: str = "placement-worker",
           timeout: float = 75) -> tuple[subprocess.CompletedProcess, list[dict]]:
    request = {"requestId": request_id, "operation": operation,
               "workingDirectory": str(work), "payload": data}
    environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT),
                       PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT))
    completed = subprocess.run(
        [PYTHON, "-m", "ProteinInMembraneSystem.worker"],
        input=json.dumps(request) + "\n", text=True, capture_output=True,
        cwd=SOURCE_ROOT, env=environment, check=False, timeout=timeout)
    return completed, [json.loads(line) for line in completed.stdout.splitlines()]


def terminal(events: list[dict], kind: str = "result") -> dict:
    assert events and len([event for event in events if event["kind"] in {"result", "error"}]) == 1
    assert events[-1]["kind"] == kind
    return events[-1]["payload"]


def fake_ppm(root: Path, *, output: str | None = None, exit_code: int = 0,
             sleep_seconds: float = 0) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    source = root / "supplied-output.pdb"
    if output is not None:
        source.write_text(output, encoding="utf-8")
    script = root / "fake-immers"
    script.write_text(
        "#!/usr/bin/env python3\nfrom pathlib import Path\nimport sys, time\n"
        f"time.sleep({sleep_seconds!r})\n"
        + (f"Path('proteinout.pdb').write_bytes(Path({str(source)!r}).read_bytes())\n"
           if output is not None else "")
        + f"sys.exit({exit_code})\n", encoding="utf-8")
    script.chmod(0o755)
    return script


class PlacementProviderCrossing(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for asset, expected in ((FIXTURE, FIXTURE_SHA256), (PPM, PPM_SHA256),
                                (RESIDUE_LIBRARY, LIBRARY_SHA256)):
            if not asset.is_file():
                raise unittest.SkipTest(f"Identified local PPM crossing asset is missing: {asset}")
            if digest(asset) != expected:
                raise AssertionError(f"Identified PPM crossing asset changed: {asset}")
        cls.temporary = tempfile.TemporaryDirectory(prefix="pim-placement-provider-")
        cls.root = Path(cls.temporary.name)
        cls.real_work = cls.root / "real"
        cls.real_work.mkdir()
        cls.real_payload = payload(cls.real_work)
        cls.real_completed, cls.real_events = invoke(cls.real_work, "place_ppm", cls.real_payload,
                                                     request_id="6qwr-real-ppm")
        if cls.real_completed.returncode != 0:
            raise AssertionError(cls.real_completed.stdout + "\n" + cls.real_completed.stderr)
        cls.real_result = terminal(cls.real_events)
        cls.raw = cls.real_work / "ppm" / "proteinout.pdb"

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "temporary"):
            cls.temporary.cleanup()

    def new_work(self, name: str, *, executable: Path = PPM,
                 library: Path = RESIDUE_LIBRARY) -> tuple[Path, dict]:
        work = self.root / name
        work.mkdir()
        return work, payload(work, executable=executable, library=library)

    def test_real_ppm_correspondence_and_all_atom_transform(self) -> None:
        self.assertEqual(0, self.real_completed.returncode)
        self.assertEqual("", self.real_completed.stderr)
        self.assertTrue(all(event["requestId"] == "6qwr-real-ppm" for event in self.real_events))
        self.assertEqual(["started", "orientationProviderStarted", "orientationCandidateObserved"],
                         [event["payload"]["stage"] for event in self.real_events if event["kind"] == "progress"])
        self.assertEqual("result", self.real_events[-1]["kind"])
        self.assertEqual("in", self.real_events[1]["payload"]["detail"]["ppmNterminalSide"])
        observed = self.real_result["observations"]
        self.assertEqual(3205, observed["alignedSourceAtomCount"])
        self.assertAlmostEqual(0, observed["midplaneAngstrom"], places=3)
        self.assertAlmostEqual(19.8, observed["thicknessAngstrom"], places=3)
        self.assertAlmostEqual(14.6, observed["tiltDegrees"], places=3)
        self.assertIn("symmetric DOPC", observed["assumedMembrane"])
        self.assertIn(PPM_SHA256, self.real_result["provider"]["version"])
        self.assertIn(LIBRARY_SHA256, self.real_result["provider"]["version"])
        self.assertEqual(LIBRARY_SHA256, digest(self.real_work / "ppm" / "res.lib"))
        self.assertEqual(PPM_SHA256, digest(self.real_work / "ppm" / "immers"))
        self.assertEqual(FIXTURE_SHA256, digest(self.real_work / "ppm" / "protein.pdb"))

        source_lines = atoms(FIXTURE)
        output = self.real_work / "oriented-protein.pdb"
        final_lines = atoms(output)
        provider_lines = atoms(self.raw)
        heavy_source = [line for line in source_lines if line[76:78].strip() not in {"H", "D"}]
        self.assertEqual((3205, 1629, 1629), (len(source_lines), len(heavy_source), len(provider_lines)))
        self.assertEqual([signature(line) for line in source_lines],
                         [signature(line) for line in final_lines])
        original_text_lines = FIXTURE.read_text(encoding="utf-8").splitlines()
        final_text_lines = output.read_text(encoding="utf-8").splitlines()
        self.assertEqual([["1.000", "1.000", "1.000"]],
                         [line.split()[1:4] for line in original_text_lines if line.startswith("CRYST1")])
        self.assertFalse(any(line.startswith("CRYST1") for line in final_text_lines))
        original_text_lines = [line for line in original_text_lines if not line.startswith("CRYST1")]
        self.assertEqual(len(original_text_lines), len(final_text_lines))
        self.assertTrue(all((before[:30] + before[54:] == after[:30] + after[54:]
                             if before.startswith(("ATOM  ", "HETATM")) else before == after)
                            for before, after in zip(original_text_lines, final_text_lines)))
        self.assertEqual([signature(line) for line in heavy_source],
                         [signature(line) for line in provider_lines])
        self.assertEqual(1576, sum(line[76:78].strip() in {"H", "D"} for line in final_lines))
        self.assertFalse(any(line[17:20].strip() == "DUM" for line in final_lines))
        self.assertEqual(1466, len(observed["planeMarkerIds"]))
        marker_lines = [line for line in self.raw.read_text(encoding="utf-8").splitlines()
                        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "DUM"]
        self.assertEqual({-9.9, 9.9}, {float(line[46:54]) for line in marker_lines})
        self.assertEqual({"N", "O"}, {line[12:16].strip() for line in marker_lines})
        self.assertTrue((self.real_work / "ppm" / "datapar1").is_file())
        self.assertTrue((self.real_work / "ppm" / "datasub1").is_file())

        # Independent least-squares affine fit, unlike the worker's SVD
        # rigid fit: the observed matrix must itself be orthogonal and proper.
        design = np.column_stack((xyz(heavy_source), np.ones(len(heavy_source))))
        fitted, _, _, _ = np.linalg.lstsq(design, xyz(provider_lines), rcond=None)
        matrix = fitted[:3, :]
        self.assertLess(np.linalg.norm(matrix.T @ matrix - np.eye(3)), 0.0001)
        self.assertAlmostEqual(1.0, np.linalg.det(matrix), places=4)
        self.assertLess(np.linalg.norm(design @ fitted - xyz(provider_lines), axis=1).max(), 0.005)
        all_fitted = np.column_stack((xyz(source_lines), np.ones(len(source_lines)))) @ fitted
        self.assertLess(np.linalg.norm(all_fitted - xyz(final_lines), axis=1).max(), 0.005)
        roles = {item["role"]: item for item in self.real_result["artifacts"]}
        self.assertEqual({"orientedPdb", "ppmRawOutput", "ppmStdout", "ppmStderr"}, set(roles))
        for item in roles.values():
            self.assertEqual(digest(Path(item["path"])), item["sha256"])
        self.assertIn("thickn= 19.8", (self.real_work / "ppm" / "ppm-stdout.txt").read_text())

    def test_identified_input_refusals(self) -> None:
        cases = {
            "changed-executable": ("ppmExecutableSha256", "0" * 64, "inputMismatch"),
            "changed-library": ("ppmResidueLibrarySha256", "0" * 64, "inputMismatch"),
            "changed-prepared": ("preparedSha256", "0" * 64, "inputMismatch"),
            "missing-prepared-digest": ("preparedSha256", None, "invalidRequest"),
            "missing-library": ("ppmResidueLibraryPath", str(self.root / "absent-res.lib"), "dependencyUnavailable"),
            "missing-binary": ("ppmExecutablePath", str(self.root / "absent-immers"), "dependencyUnavailable"),
            "wrong-side": ("ppmNterminalSide", "neither", "invalidSelection"),
        }
        for name, (field, value, code) in cases.items():
            with self.subTest(name=name):
                work, request = self.new_work(name)
                request[field] = value
                completed, events = invoke(work, "place_ppm", request)
                self.assertNotEqual(0, completed.returncode)
                self.assertEqual(code, terminal(events, "error")["failureCode"])
                self.assertFalse((work / "oriented-protein.pdb").exists())

    def test_controlled_provider_refusals_and_first_model_boundary(self) -> None:
        raw = self.raw.read_text(encoding="utf-8")
        lines = raw.splitlines(keepends=True)
        heavy_index = next(i for i, line in enumerate(lines) if line.startswith("ATOM  "))
        marker_index = next(i for i, line in enumerate(lines) if line[17:20].strip() == "DUM")
        changed_identity = lines.copy()
        line = changed_identity[heavy_index]
        changed_identity[heavy_index] = line[:12] + " ZZ " + line[16:]
        nonrigid = lines.copy()
        line = nonrigid[heavy_index]
        nonrigid[heavy_index] = line[:30] + f"{float(line[30:38]) + 0.1:8.3f}" + line[38:]
        missing_plane = [line for line in lines if line[17:20].strip() != "DUM"]
        changed_plane = lines.copy()
        line = changed_plane[marker_index]
        changed_plane[marker_index] = line[:46] + f"{float(line[46:54]) + 0.1:8.3f}" + line[54:]
        source_h = next(line for line in FIXTURE.read_text(encoding="utf-8").splitlines(keepends=True)
                        if line.startswith("ATOM  ") and line[76:78].strip() == "H")
        partial_h = lines.copy()
        partial_h.insert(marker_index, source_h)
        cases = {
            "identity": ("".join(changed_identity), "providerMismatch"),
            "nonrigid": ("".join(nonrigid), "providerMismatch"),
            "missing-plane": ("".join(missing_plane), "invalidProviderOutput"),
            "inconsistent-plane": ("".join(changed_plane), "invalidProviderOutput"),
            "partial-hydrogen": ("".join(partial_h), "providerMismatch"),
            "no-output": (None, "unobservedOutput"),
        }
        for name, (provided, code) in cases.items():
            with self.subTest(name=name):
                fake = fake_ppm(self.root / f"tool-{name}", output=provided)
                work, request = self.new_work(f"work-{name}", executable=fake)
                completed, events = invoke(work, "place_ppm", request)
                self.assertNotEqual(0, completed.returncode)
                self.assertEqual(code, terminal(events, "error")["failureCode"])
                self.assertFalse((work / "oriented-protein.pdb").exists())

        fake = fake_ppm(self.root / "tool-nonzero", exit_code=7)
        work, request = self.new_work("work-nonzero", executable=fake)
        completed, events = invoke(work, "place_ppm", request)
        self.assertEqual("providerFailed", terminal(events, "error")["failureCode"])
        self.assertEqual(7, terminal(events, "error")["details"]["exitCode"])

        # PPM may append an original-coordinate model after ENDMDL. Its atoms
        # must not enter the positioned artifact or alter the first-model fit.
        fake = fake_ppm(self.root / "tool-appended", output=raw + "ENDMDL\nMODEL        2\n" +
                        FIXTURE.read_text(encoding="utf-8"))
        work, request = self.new_work("work-appended", executable=fake)
        completed, events = invoke(work, "place_ppm", request)
        self.assertEqual(0, completed.returncode, completed.stdout)
        self.assertEqual(3205, terminal(events)["observations"]["alignedSourceAtomCount"])
        self.assertEqual(3205, len(atoms(work / "oriented-protein.pdb")))

    def test_timeout_and_cancellation_do_not_publish_a_candidate(self) -> None:
        fake = fake_ppm(self.root / "tool-timeout", sleep_seconds=1.0)
        work, request = self.new_work("work-timeout", executable=fake)
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.PlacementAssessment.worker import placement_assessment
        from ProteinInMembraneSystem.worker.exchange import WorkError
        try:
            with mock.patch.object(placement_assessment, "_PPM_TIMEOUT_SECONDS", 0.05):
                with self.assertRaises(WorkError) as caught:
                    placement_assessment.place_ppm(work, request, lambda _stage, _detail: None)
            self.assertEqual("providerTimeout", caught.exception.code)
            self.assertFalse((work / "oriented-protein.pdb").exists())
        finally:
            sys.path.remove(str(SOURCE_ROOT))

        fake = fake_ppm(self.root / "tool-cancel", sleep_seconds=5.0)
        work, request = self.new_work("work-cancel", executable=fake)
        envelope = {"requestId": "cancelled-ppm", "operation": "place_ppm",
                    "workingDirectory": str(work), "payload": request}
        environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT),
                           PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT))
        process = subprocess.Popen([PYTHON, "-m", "ProteinInMembraneSystem.worker"],
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, cwd=SOURCE_ROOT, env=environment, start_new_session=True)
        assert process.stdin is not None and process.stdout is not None
        process.stdin.write(json.dumps(envelope) + "\n")
        process.stdin.flush()
        progress = [json.loads(process.stdout.readline()), json.loads(process.stdout.readline())]
        self.assertEqual(["started", "orientationProviderStarted"],
                         [item["payload"]["stage"] for item in progress])
        os.killpg(process.pid, signal.SIGTERM)
        remaining, _ = process.communicate(timeout=5)
        self.assertNotEqual(0, process.returncode)
        self.assertFalse(any(item.get("kind") == "result" for item in
                             [json.loads(line) for line in remaining.splitlines()]))
        self.assertFalse((work / "oriented-protein.pdb").exists())

    def test_lost_final_hydrogen_fails_correspondence_guard(self) -> None:
        fake = fake_ppm(self.root / "tool-drop-h", output=self.raw.read_text(encoding="utf-8"))
        work, request = self.new_work("work-drop-h", executable=fake)
        sys.path.insert(0, str(SOURCE_ROOT))
        from ProteinInMembraneSystem.PlacementAssessment.worker import placement_assessment
        from ProteinInMembraneSystem.worker.exchange import WorkError

        def broken_reconstruction(source_lines, target, _rotation, _translation):
            target.write_text("".join(line for line in source_lines
                                      if not (line.startswith(("ATOM  ", "HETATM")) and
                                              line[76:78].strip() in {"H", "D"})), encoding="utf-8")
            return len(atoms(FIXTURE))  # A false count alone must not pass.

        try:
            with mock.patch.object(placement_assessment, "_write_oriented", broken_reconstruction):
                with self.assertRaises(WorkError) as caught:
                    placement_assessment.place_ppm(work, request, lambda _stage, _detail: None)
            self.assertEqual("providerMismatch", caught.exception.code)
        finally:
            sys.path.remove(str(SOURCE_ROOT))


def atom_line(serial: int, atom: str, residue: str, chain: str, number: int,
              position: tuple[float, float, float], element: str) -> str:
    x, y, z = position
    return (f"ATOM  {serial:5d} {atom:4s} {residue:3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n")


def small_placement() -> str:
    return (atom_line(1, "N", "ALA", "A", 1, (0, 0, -2), "N") +
            atom_line(2, "CA", "ALA", "A", 1, (1, 0, -0.5), "C") +
            atom_line(3, "H", "ALA", "A", 1, (0, 1, -0.5), "H") +
            atom_line(4, "N", "ALA", "B", 2, (0, 0, 0.5), "N") +
            atom_line(5, "CA", "ALA", "B", 2, (1, 0, 2), "C") +
            "END\n")


def source_addresses() -> list[dict]:
    return [{"model": 1, "chain": "X", "copyId": "X1", "residue": 1, "insertionCode": ""},
            {"model": 1, "chain": "Y", "copyId": "Y1", "residue": 2, "insertionCode": ""}]


def result_atom_ids() -> list[str]:
    return ["0:A:1::N", "1:A:1::CA", "2:A:1::H", "3:B:2::N", "4:B:2::CA"]


class PlacementGeometryCrossing(unittest.TestCase):
    def test_exact_residue_geometry_and_bounded_correction(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pim-placement-geometry-") as temporary:
            root = Path(temporary)
            input_pdb = root / "small.pdb"
            input_pdb.write_text(small_placement(), encoding="utf-8")
            measured, events = invoke(root, "measure_placement", {
                "studyRevisionId": "geometry-revision", "proposalId": "original",
                "orientedPdbPath": str(input_pdb), "orientedPdbSha256": digest(input_pdb),
                "membraneMidplaneAngstrom": 0,
                "coreLowerZAngstrom": -1, "coreUpperZAngstrom": 1,
                "residueAddressesInOrder": source_addresses(),
                "outputChainIdsInOrder": ["A", "B"],
                "expectedResultAtomIdsInOrder": result_atom_ids(),
            })
            self.assertEqual(0, measured.returncode, measured.stdout)
            observation = terminal(events)["observations"]
            self.assertEqual((5, 3, 1, 1),
                             (observation["atomCount"], observation["atomsWithinCore"],
                              observation["atomsAboveCore"], observation["atomsBelowCore"]))
            self.assertEqual((-2, 2), (observation["proteinZMinAngstrom"],
                                        observation["proteinZMaxAngstrom"]))
            first, second = observation["residues"]
            self.assertEqual((3, -2, -0.5, 2, 0, 1),
                             (first["atomCount"], first["minZAngstrom"], first["maxZAngstrom"],
                              first["atomsWithinCore"], first["atomsAboveCore"], first["atomsBelowCore"]))
            self.assertEqual((2, 0.5, 2, 1, 1, 0),
                             (second["atomCount"], second["minZAngstrom"], second["maxZAngstrom"],
                              second["atomsWithinCore"], second["atomsAboveCore"], second["atomsBelowCore"]))
            self.assertEqual((source_addresses(), ["A", "B"]),
                             ([item["address"] for item in observation["residues"]],
                              [item["outputChainId"] for item in observation["residues"]]))

            adjusted, adjust_events = invoke(root, "adjust_placement", {
                "studyRevisionId": "geometry-revision", "sourceProposalId": "original",
                "orientedPdbPath": str(input_pdb), "orientedPdbSha256": digest(input_pdb),
                "depthShiftAngstrom": 1,
                "tiltAboutXDegrees": 0, "tiltAboutYDegrees": 0,
                "rotationAboutNormalDegrees": 0, "rationale": "Independently check a 1 Å shift",
            })
            self.assertEqual(0, adjusted.returncode, adjusted.stdout)
            adjusted_pdb = Path(terminal(adjust_events)["artifacts"][0]["path"])
            self.assertEqual([signature(line) for line in atoms(input_pdb)],
                             [signature(line) for line in atoms(adjusted_pdb)])
            np.testing.assert_allclose(xyz(atoms(adjusted_pdb)) - xyz(atoms(input_pdb)),
                                       np.tile([0, 0, 1], (5, 1)), atol=1e-9)
            shifted, shifted_events = invoke(root, "measure_placement", {
                "studyRevisionId": "geometry-revision", "proposalId": "shifted",
                "orientedPdbPath": str(adjusted_pdb), "orientedPdbSha256": digest(adjusted_pdb),
                "membraneMidplaneAngstrom": 0,
                "coreLowerZAngstrom": -1, "coreUpperZAngstrom": 1,
                "residueAddressesInOrder": source_addresses(),
                "outputChainIdsInOrder": ["A", "B"],
                "expectedResultAtomIdsInOrder": result_atom_ids(),
            })
            self.assertEqual(0, shifted.returncode, shifted.stdout)
            shifted_observation = terminal(shifted_events)["observations"]
            self.assertEqual((3, 2, 0),
                             (shifted_observation["atomsWithinCore"], shifted_observation["atomsAboveCore"],
                              shifted_observation["atomsBelowCore"]))

            rotated_dir = root / "rotated"
            rotated_dir.mkdir()
            rotated_source = rotated_dir / "small.pdb"
            shutil.copyfile(input_pdb, rotated_source)
            rotated, rotate_events = invoke(rotated_dir, "adjust_placement", {
                "studyRevisionId": "geometry-revision", "sourceProposalId": "original",
                "orientedPdbPath": str(rotated_source), "orientedPdbSha256": digest(rotated_source),
                "depthShiftAngstrom": 0.25,
                "tiltAboutXDegrees": 90, "tiltAboutYDegrees": 0,
                "rotationAboutNormalDegrees": 90, "rationale": "Check finite orthogonal rotations",
            })
            self.assertEqual(0, rotated.returncode, rotated.stdout)
            rotated_pdb = Path(terminal(rotate_events)["artifacts"][0]["path"])
            before, after = xyz(atoms(rotated_source)), xyz(atoms(rotated_pdb))
            center = before.mean(axis=0)
            # Independent elementary rotation: X(90°), then Z(90°).
            expected = np.column_stack((before[:, 2] - center[2] + center[0],
                                        before[:, 0] - center[0] + center[1],
                                        before[:, 1] - center[1] + center[2] + 0.25))
            np.testing.assert_allclose(after, expected, atol=0.001)
            np.testing.assert_allclose(np.linalg.norm(before[:, None] - before[None, :], axis=2),
                                       np.linalg.norm(after[:, None] - after[None, :], axis=2), atol=0.002)

    def test_correspondence_and_structure_refusals(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pim-placement-refusals-") as temporary:
            root = Path(temporary)
            original = small_placement()
            mutations = {
                "swapped-output-chain": (original, source_addresses(), ["B", "A"], result_atom_ids(), "correspondenceFailed"),
                "conflicting-copy": (original, [source_addresses()[0],
                                                {**source_addresses()[1], "chain": "X", "copyId": "X1"}],
                                     ["A", "B"], result_atom_ids(), "correspondenceFailed"),
                "swapped-residue": (original, list(reversed(source_addresses())), ["A", "B"],
                                    result_atom_ids(), "correspondenceFailed"),
                "missing-hydrogen": (original.replace(atom_line(3, "H", "ALA", "A", 1,
                                                               (0, 1, -0.5), "H"), ""),
                                     source_addresses(), ["A", "B"], result_atom_ids(), "correspondenceFailed"),
                "changed-hydrogen-name": (original.replace(
                    atom_line(3, "H", "ALA", "A", 1, (0, 1, -0.5), "H"),
                    atom_line(3, "Q", "ALA", "A", 1, (0, 1, -0.5), "H")),
                                          source_addresses(), ["A", "B"], result_atom_ids(), "correspondenceFailed"),
                "marker-contamination": (original.replace("END\n", atom_line(6, "N", "DUM", "B", 3,
                                                                                 (0, 0, -1), "N") + "END\n"),
                                         source_addresses(), ["A", "B"], result_atom_ids(), "invalidStructure"),
                "malformed-coordinate": (original.replace("  -2.000", "    oops"),
                                         source_addresses(), ["A", "B"], result_atom_ids(), "invalidStructure"),
            }
            for name, (pdb_text, addresses, chains, atom_ids, code) in mutations.items():
                with self.subTest(name=name):
                    work = root / name
                    work.mkdir()
                    path = work / "tampered.pdb"
                    path.write_text(pdb_text, encoding="utf-8")
                    completed, events = invoke(work, "measure_placement", {
                        "studyRevisionId": "geometry-revision", "proposalId": "tampered",
                        "orientedPdbPath": str(path), "orientedPdbSha256": digest(path),
                        "membraneMidplaneAngstrom": 0,
                        "coreLowerZAngstrom": -1, "coreUpperZAngstrom": 1,
                        "residueAddressesInOrder": addresses,
                        "outputChainIdsInOrder": chains,
                        "expectedResultAtomIdsInOrder": atom_ids,
                    })
                    self.assertNotEqual(0, completed.returncode)
                    self.assertEqual(code, terminal(events, "error")["failureCode"])

    def test_stale_artifact_and_invalid_adjustment_refusals(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pim-placement-adjust-refusal-") as temporary:
            root = Path(temporary)
            input_pdb = root / "source.pdb"
            input_pdb.write_text(small_placement(), encoding="utf-8")
            base_adjust = {
                "studyRevisionId": "geometry-revision", "sourceProposalId": "exact-proposal",
                "orientedPdbPath": str(input_pdb), "orientedPdbSha256": digest(input_pdb),
                "depthShiftAngstrom": 0, "tiltAboutXDegrees": 0,
                "tiltAboutYDegrees": 0, "rotationAboutNormalDegrees": 0,
                "rationale": "Check the exact oriented artifact",
            }
            base_measure = {
                "studyRevisionId": "geometry-revision", "proposalId": "exact-proposal",
                "orientedPdbPath": str(input_pdb), "orientedPdbSha256": digest(input_pdb),
                "membraneMidplaneAngstrom": 0,
                "coreLowerZAngstrom": -1, "coreUpperZAngstrom": 1,
                "residueAddressesInOrder": source_addresses(),
                "outputChainIdsInOrder": ["A", "B"],
                "expectedResultAtomIdsInOrder": result_atom_ids(),
            }
            for operation, original in (("adjust_placement", base_adjust),
                                        ("measure_placement", base_measure)):
                for name, changes, code in (
                    ("digest", {"orientedPdbSha256": "0" * 64}, "inputMismatch"),
                    ("missing-digest", {"orientedPdbSha256": None}, "invalidRequest"),
                ):
                    with self.subTest(operation=operation, name=name):
                        completed, events = invoke(root, operation, dict(original, **changes))
                        self.assertNotEqual(0, completed.returncode)
                        self.assertEqual(code, terminal(events, "error")["failureCode"])
            for name, changes, code in (
                ("no-rationale", {"rationale": ""}, "invalidRequest"),
                ("nonfinite-depth", {"depthShiftAngstrom": float("nan")}, "invalidRequest"),
                ("unrepresentable-depth", {"depthShiftAngstrom": -1001}, "unsupportedRepresentation"),
            ):
                with self.subTest(name=name):
                    completed, events = invoke(root, "adjust_placement", dict(base_adjust, **changes))
                    self.assertNotEqual(0, completed.returncode)
                    self.assertEqual(code, terminal(events, "error")["failureCode"])
            self.assertFalse((root / "adjusted-placement.pdb").exists())


class PredictionRegionCrossing(unittest.TestCase):
    def test_asymmetric_directions_coverage_and_exact_evidence_identity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pim-placement-pae-") as temporary:
            root = Path(temporary)
            pae = root / "pae.json"
            pae.write_text(json.dumps([{"predicted_aligned_error": [[0, 2, 4], [7, 0, 5], [9, 11, 0]],
                                        "max_predicted_aligned_error": 20}]), encoding="utf-8")
            addresses = [{"model": 0, "chain": "A", "copyId": "A", "residue": number,
                          "insertionCode": ""} for number in (1, 2, 3)]
            mapping = root / "axis.json"
            mapping.write_text(json.dumps({"recordId": "AF-exact-record", "coordinateSha256": "a" * 64,
                                           "paeSha256": digest(pae), "axisResidues": addresses}), encoding="utf-8")
            request = {"studyRevisionId": "prediction-revision", "recordId": "AF-exact-record",
                       "coordinateSha256": "a" * 64, "paePath": str(pae), "paeSha256": digest(pae),
                       "paeMappingPath": str(mapping), "paeMappingSha256": digest(mapping),
                       "firstRegion": addresses[:2], "secondRegion": addresses[2:],
                       "maximumReportedPairs": 2}
            completed, events = invoke(root, "summarize_prediction_evidence", request)
            self.assertEqual(0, completed.returncode, completed.stdout)
            observation = terminal(events)["observations"]
            self.assertEqual("Observed", observation["standing"])
            forward, reverse = observation["firstAlignedOnSecond"], observation["secondAlignedOnFirst"]
            self.assertEqual((2, 2, 9, 11, 10),
                             (forward["possiblePairCount"], forward["validPairCount"],
                              forward["minimumAngstrom"], forward["maximumAngstrom"], forward["meanAngstrom"]))
            self.assertEqual((2, 2, 4, 5, 4.5),
                             (reverse["possiblePairCount"], reverse["validPairCount"],
                              reverse["minimumAngstrom"], reverse["maximumAngstrom"], reverse["meanAngstrom"]))
            self.assertEqual([11, 9], [item["predictedAlignedErrorAngstrom"]
                                       for item in forward["highestErrorPairs"]])
            missing = dict(request, firstRegion=[{**addresses[0], "residue": 99}])
            completed, events = invoke(root, "summarize_prediction_evidence", missing)
            self.assertEqual(0, completed.returncode)
            self.assertEqual("Unavailable", terminal(events)["observations"]["standing"])
            for name, changed in {
                "pae-hash": {"paeSha256": "0" * 64},
                "mapping-hash": {"paeMappingSha256": "0" * 64},
                "coordinate-hash": {"coordinateSha256": "b" * 64},
                "record": {"recordId": "AF-other-record"},
            }.items():
                with self.subTest(name=name):
                    completed, events = invoke(root, "summarize_prediction_evidence", dict(request, **changed))
                    self.assertNotEqual(0, completed.returncode)
                    self.assertEqual("inputMismatch", terminal(events, "error")["failureCode"])
