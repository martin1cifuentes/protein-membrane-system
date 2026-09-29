"""Unknown-provenance software construct through the named native DMPC route.

The two-alanine PDB is deliberately small. Its declared unknown provenance is
an intake and technical-readiness case, not a biological or predicted-protein
claim. The published browser, Host and worker must complete one authorized
native construction, final minimization, export and independent read-back.
"""

from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile
import unittest
import zipfile

from playwright.sync_api import expect, sync_playwright

from test_manual_position import (ROOT, POLICY, chromium, prepare_with_current_plan,
                                  running_host, state, two_alanines, until)

ARTIFACTS = ROOT / "out" / "browser-acceptance" / "unknown-source-native"
NATIVE_POLICY_ID = "canonical-protein-pure-dmpc-native-construction"
READ_BACK = r"""
import hashlib, json, math, sys, zipfile
from collections import Counter

import gemmi
import numpy as np
from openmm import HarmonicBondForce, XmlSerializer, unit

path, stage_id, attempt_id, revision_id, placement_proposal_id, source_path, source_sha, assessment_id, standing, count = sys.argv[1:]
count = int(count)
assert hashlib.sha256(open(source_path, "rb").read()).hexdigest() == source_sha
required = {"structure.cif", "topology.json", "system.xml", "state.xml"}
with zipfile.ZipFile(path) as archive:
    assert archive.testzip() is None
    assert required | {"manifest.json"} <= set(archive.namelist())
    content = {name: archive.read(name) for name in required}
    manifest = json.loads(archive.read("manifest.json"))

assert manifest["stage"]["id"] == stage_id
assert manifest["stage"]["kind"] == "Minimization"
assert manifest["stage"]["atomCount"] == count
assert manifest["attempt"]["id"] == attempt_id
assert manifest["attempt"]["policyId"] == "canonical-protein-pure-dmpc-native-construction"
assert manifest["attempt"]["placementId"] == manifest["placement"]["id"]
assert manifest["placement"]["proposal"]["id"] == placement_proposal_id
assert manifest["placement"]["studyRevisionId"] == revision_id
assert manifest["study"]["id"] == revision_id
assert manifest["assessment"]["id"] == assessment_id
assert manifest["assessment"]["stageId"] == stage_id
assert manifest["assessment"]["checkStanding"] == standing
assert manifest["preparedProtein"]["source"]["uploadProvenance"] == "unknown"
assert manifest["lineage"]["sourceCoordinateSha256"] == source_sha
assert manifest["construction"]["provider"]["route"] == "nativeOpenMm"
entries = {item["name"]: item for item in manifest["artifacts"]}
assert set(entries) == required
for name, data in content.items():
    assert entries[name]["sha256"] == hashlib.sha256(data).hexdigest()
    assert entries[name]["byteLength"] == len(data)
for name, field in (("structure.cif", "coordinatesSha256"),
                    ("topology.json", "topologySha256"),
                    ("system.xml", "systemXmlSha256"),
                    ("state.xml", "stateXmlSha256")):
    assert manifest["molecularIdentity"][field] == entries[name]["sha256"]

topology = json.loads(content["topology.json"])
atoms, residues, chains = topology["atoms"], topology["residues"], topology["chains"]
bonds = {tuple(sorted(bond["atomIndices"])) for bond in topology["bonds"]}
assert len(atoms) == count and len(bonds) == len(topology["bonds"]) and bonds
assert all(0 <= left < right < count for left, right in bonds)
mapped = manifest["correspondence"]["atoms"]
assert len(mapped) == count
roles = Counter()
source_atoms = {}
for line in open(source_path, encoding="utf-8"):
    if line.startswith("ATOM  ") and line[76:78].strip().upper() != "H":
        source_id = f"0:A:A:{int(line[22:26])}:{line[26].strip()}:{line[12:16].strip()}"
        assert source_id not in source_atoms
        source_atoms[source_id] = (line[12:16].strip(), line[76:78].strip().upper())
assert source_atoms
mapped_source = {}
block = gemmi.cif.read_string(content["structure.cif"].decode()).sole_block()
table = block.find("_atom_site.", ["type_symbol", "label_atom_id", "label_comp_id",
    "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
    "Cartn_z", "pdbx_PDB_model_num"])
assert len(table) == count
positions = np.empty((count, 3), dtype=float)
models = set()
for index, row in enumerate(table):
    atom = atoms[index]
    residue = residues[atom["residueIndex"]]
    chain = chains[residue["chainIndex"]]
    assert (row[0], row[1], row[2]) == (atom["element"], atom["name"], residue["name"])
    assert (row[3], row[4], "" if row[5] in (".", "?") else row[5]) == (
        chain["id"], residue["id"], residue["insertionCode"])
    assert mapped[index]["resultAtomIndex"] == index
    assert mapped[index]["element"] == atom["element"]
    roles[mapped[index]["moleculeRole"]] += 1
    if mapped[index]["sourceAtomId"] in source_atoms:
        source_id = mapped[index]["sourceAtomId"]
        assert source_id not in mapped_source
        assert mapped[index]["role"] == "source" and mapped[index]["moleculeRole"] == "protein"
        assert (atom["name"], atom["element"]) == source_atoms[source_id]
        mapped_source[source_id] = index
    positions[index] = [float(row[6]), float(row[7]), float(row[8])]
    models.add(row[9])
assert models == {"1"}
assert set(mapped_source) == set(source_atoms)
assert all(roles[role] > 0 for role in ("protein", "lipid", "water", "ion"))

system = XmlSerializer.deserialize(content["system.xml"].decode())
state = XmlSerializer.deserialize(content["state.xml"].decode())
assert system.getNumParticles() == count
assert len(state.getPositions()) == count
realized = {tuple(sorted(system.getConstraintParameters(i)[:2]))
    for i in range(system.getNumConstraints())}
for force in system.getForces():
    if isinstance(force, HarmonicBondForce):
        realized.update(tuple(sorted(force.getBondParameters(i)[:2]))
            for i in range(force.getNumBonds()))
assert bonds <= realized
state_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
coordinate_error = float(np.max(np.linalg.norm(positions - state_positions, axis=1)))
assert math.isfinite(coordinate_error) and coordinate_error <= 0.01
vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
topology_vectors = np.asarray(topology["boxVectorsAngstrom"], dtype=float)
cell_error = float(np.max(np.abs(vectors - topology_vectors)))
assert math.isfinite(cell_error) and cell_error <= 0.01
assert np.all(np.linalg.norm(vectors, axis=1) > 0)
lengths = np.linalg.norm(vectors, axis=1)
cif_lengths = np.asarray([float(block.find_value(f"_cell.length_{axis}"))
    for axis in ("a", "b", "c")])
assert float(np.max(np.abs(cif_lengths - lengths))) <= 0.01
angles = []
for left, right in ((1, 2), (0, 2), (0, 1)):
    cosine = np.dot(vectors[left], vectors[right]) / (lengths[left] * lengths[right])
    angles.append(math.degrees(math.acos(float(np.clip(cosine, -1, 1)))))
cif_angles = np.asarray([float(block.find_value(f"_cell.angle_{axis}"))
    for axis in ("alpha", "beta", "gamma")])
assert float(np.max(np.abs(cif_angles - angles))) <= 0.1
print(json.dumps({"stageId": stage_id, "attemptId": attempt_id,
    "studyRevisionId": revision_id, "assessmentId": assessment_id,
    "atomCount": count, "bondCount": len(bonds), "moleculeRoleCounts": dict(roles),
    "coordinateMaxDeviationAngstrom": coordinate_error,
    "cellMaxVectorDifferenceAngstrom": cell_error,
    "cellLengthsAngstrom": lengths.tolist(), "cellAnglesDegrees": angles,
    "bundleSha256": hashlib.sha256(open(path, "rb").read()).hexdigest()}, indent=2))
"""


class NativeConstructionBrowserTests(unittest.TestCase):
    def test_unknown_software_construct_completes_named_native_route_and_export(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        run_artifacts = Path(tempfile.mkdtemp(prefix="unknown-ala2-native-", dir=ARTIFACTS))
        # Retain the exact Host workspace if a late actor assertion fails.
        with nullcontext(str(run_artifacts)) as temporary:
            directory = Path(temporary)
            source = run_artifacts / "two-alanines-software-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n", encoding="utf-8")
            with running_host(directory / "workspace", policy=POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                page = None
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940}, accept_downloads=True)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("unknown")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    uploaded = until(page, lambda current: (current.get("study") or {}).get(
                        "uploadProvenance") == "unknown", "unknown-provenance source intake")
                    self.assertEqual(uploaded["study"]["selectedSourceKind"], "upload")
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    prepare_with_current_plan(page)
                    expect(page.locator(".work-pending")).to_have_count(0, timeout=120000)
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Membrane").click()
                    for side in ("Upper", "Lower"):
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                    use = page.get_by_role("button", name="Use this membrane")
                    expect(use).to_be_enabled(timeout=120000)
                    use.click()
                    until(page, lambda current: (current.get("membrane") or {}).get("status") == "assessed",
                          "checked DMPC")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Placement").click()
                    page.locator("#starting-position").select_option("upper")
                    positioned = until(page, lambda current: (current.get("placement") or {}).get("status") ==
                                       "supported" and (current["placement"].get("transform") or {})
                                       .get("startingPosition") == "upper", "checked upper-side pose")
                    self.assertEqual(positioned["placement"]["preparedProteinId"],
                                     positioned["protein"]["subjectId"])
                    self.assertEqual(positioned["placement"]["membraneModelId"],
                                     positioned["membrane"]["modelId"])
                    page.get_by_role("button", name="Use this position").click()
                    adopted = until(page, lambda current: current["study"]["adoptedPlacementProposalId"] ==
                                    positioned["placement"]["proposalId"], "adopted pose")
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=180000)
                    for width, height in ((1672, 941), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        expect(page.locator(".viewer-mount[data-camera-ready='true']"))\
                            .to_be_visible(timeout=120000)
                        page.evaluate("() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))")
                        expect(page.get_by_label("Position proposal and scientific support")).to_contain_text(
                            "Selected in the current study revision")
                        page.screenshot(path=str(run_artifacts / f"unknown-upper-placement-{width}.png"),
                                        full_page=True, animations="disabled")
                    page.set_viewport_size({"width": 1024, "height": 768})
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Preparation").click()
                    route_account = next(item for item in adopted["constructionRoutes"]
                                         if item["policyId"] == NATIVE_POLICY_ID)
                    self.assertTrue(route_account["available"], route_account)
                    self.assertIn("OpenMM native: DMPC", route_account["label"])
                    route = page.locator(".construction-route").filter(has_text="OpenMM native: DMPC")
                    expect(route.get_by_role("button", name="Build and minimize")).to_be_enabled()
                    route.get_by_role("button", name="Build and minimize").click()
                    admitted = until(page, lambda current: (current.get("attempt") or {}).get(
                        "attemptId") and current["attempt"]["policyId"] == NATIVE_POLICY_ID,
                        "one named native Build and minimize attempt", seconds=300)
                    attempt_id = admitted["attempt"]["attemptId"]
                    revision_id = adopted["study"]["id"]
                    self.assertEqual(admitted["attempt"]["studyRevisionId"], revision_id)
                    self.assertEqual(admitted["study"]["adoptedPlacementProposalId"],
                                     positioned["placement"]["proposalId"])
                    self.assertEqual(admitted["stages"], [])
                    self.assertFalse(any(action["kind"] == "continueMinimization"
                                         for action in admitted["actions"]))
                    candidate = until(page, lambda current: (current.get("attempt") or {}).get("constructed")
                                      is not None or (current.get("attempt") or {}).get("status") in
                                      ("failed", "resourceRefused", "unobserved"),
                                      "checked native construction or attributable failure", seconds=1100)
                    self.assertIsNotNone(candidate["attempt"]["constructed"], candidate["attempt"])
                    self.assertEqual(candidate["attempt"]["attemptId"], attempt_id)
                    self.assertIsNotNone(candidate["attempt"]["constructed"])
                    self.assertGreater(candidate["attempt"]["constructed"]["atomCount"],
                                       candidate["protein"]["atomCount"])
                    page.get_by_role("button", name="Review attempt").click()
                    page.get_by_role("button", name="Inspect verified constructed system").click()
                    until(page, lambda current: (current.get("inspection") or {}).get("subjectId") ==
                          candidate["attempt"]["constructed"]["subjectId"], "candidate inspection")
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    page.screenshot(path=str(run_artifacts / "unknown-native-checked-intermediate-1024.png"),
                                    full_page=True, animations="disabled")
                    if candidate["attempt"]["status"] != "completed":
                        self.assertFalse(any(item["status"] == "completed" for item in candidate["stages"]))
                    completed = until(page, lambda current: (current.get("attempt") or {}).get("status") in
                                      ("completed", "failed", "stopped", "unobserved") and
                                      (current.get("attempt") or {}).get("stageKind") == "Minimization",
                                      "real minimization outcome", seconds=1800)
                    self.assertEqual(completed["attempt"]["status"], "completed", completed["attempt"])
                    self.assertEqual(completed["attempt"]["attemptId"], attempt_id)
                    self.assertFalse(any(action["kind"] == "continueMinimization"
                                         for action in completed["actions"]))
                    minimized = next(item for item in completed["stages"]
                                     if item["attemptId"] == attempt_id and
                                     item["kind"] == "Minimization")
                    self.assertEqual(minimized["status"], "completed")
                    self.assertEqual(minimized["constructed"]["subjectId"],
                                     candidate["attempt"]["constructed"]["subjectId"])
                    self.assertEqual(minimized["observation"]["termination"], "converged")
                    self.assertEqual(minimized["studyRevisionId"], revision_id)
                    self.assertIsNotNone(minimized["assessment"])
                    self.assertTrue(minimized["assessment"]["currentlyApplicable"])
                    self.assertEqual(minimized["assessment"]["stageId"], minimized["stageId"])
                    self.assertIn(minimized["assessment"]["checkStanding"],
                                  ("checksPassed", "issuesFound", "checksIncomplete"))
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Results").click()
                    page.locator(f'.stage-card[title^="Stage {minimized["stageId"]} · "]').click()
                    until(page, lambda current: (current.get("inspection") or {}).get("subjectId") ==
                          minimized["stageId"], "minimized stage inspection")
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    self.assertTrue(next(action for action in state(page)["actions"]
                                         if action["kind"] == "exportStage" and
                                         action.get("subjectId") == minimized["stageId"])["enabled"])
                    for width, height in ((1672, 941), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        expect(page.locator(".viewer-mount[data-camera-ready='true']"))\
                            .to_be_visible(timeout=120000)
                        page.evaluate("() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))")
                        expect(page.get_by_role("button", name="Export this completed stage")).to_be_enabled()
                        page.screenshot(path=str(run_artifacts / f"unknown-native-result-{width}.png"),
                                        full_page=True, animations="disabled")
                    page.set_viewport_size({"width": 1024, "height": 768})
                    with page.expect_download(timeout=300000) as delivered:
                        page.get_by_role("button", name="Export this completed stage").click()
                    bundle = run_artifacts / "unknown-native-minimized-export.zip"
                    delivered.value.save_as(str(bundle))
                    export = until(page, lambda current: (next((item for item in current["stages"]
                          if item["stageId"] == minimized["stageId"]), {}).get("export") or {})
                          .get("status") == "verified", "verified exact-stage export")
                    exported = next(item for item in export["stages"] if item["stageId"] == minimized["stageId"])
                    self.assertEqual(hashlib.sha256(bundle.read_bytes()).hexdigest(), exported["export"]["sha256"])
                    with zipfile.ZipFile(bundle) as archive:
                        self.assertIsNone(archive.testzip())
                        self.assertGreaterEqual(set(archive.namelist()),
                            {"manifest.json", "structure.cif", "topology.json", "system.xml", "state.xml"})
                        manifest = json.loads(archive.read("manifest.json"))
                        self.assertEqual(manifest["stage"]["id"], minimized["stageId"])
                        self.assertEqual(manifest["attempt"]["id"], attempt_id)
                        self.assertEqual(manifest["study"]["id"], revision_id)
                        self.assertEqual(manifest["attempt"]["policyId"], NATIVE_POLICY_ID)
                        self.assertEqual(manifest["assessment"]["id"], minimized["assessment"]["id"])
                        self.assertEqual(manifest["assessment"]["checkStanding"],
                                         minimized["assessment"]["checkStanding"])
                        self.assertEqual(manifest["preparedProtein"]["source"]["uploadProvenance"],
                                         "unknown")
                        self.assertEqual(manifest["stage"]["atomCount"],
                                         candidate["attempt"]["constructed"]["atomCount"])
                        self.assertEqual(manifest["lineage"]["sourceCoordinateSha256"],
                                         hashlib.sha256(source.read_bytes()).hexdigest())
                        plan = manifest["preparedProtein"]["recommendationPlan"]
                        self.assertEqual(len(plan["planSha256"]), 64)
                        self.assertEqual(plan["checkedCandidateSha256"],
                                         manifest["lineage"]["preparedCoordinateSha256"])
                        self.assertEqual(plan["methodVersion"], "8.6.0")
                        self.assertIsInstance(plan["seed"], int)
                    scientific_readback = subprocess.run([str(ROOT / "out" / "python" / "bin" /
                        "python"), "-c", READ_BACK, str(bundle), minimized["stageId"], attempt_id,
                        revision_id, positioned["placement"]["proposalId"], str(source),
                        hashlib.sha256(source.read_bytes()).hexdigest(),
                        minimized["assessment"]["id"], minimized["assessment"]["checkStanding"],
                        str(candidate["attempt"]["constructed"]["atomCount"])],
                        capture_output=True, text=True, timeout=180)
                    self.assertEqual(scientific_readback.returncode, 0, scientific_readback.stderr)
                    readback = json.loads(scientific_readback.stdout)
                    self.assertEqual(readback["atomCount"],
                                     candidate["attempt"]["constructed"]["atomCount"])
                    (run_artifacts / "independent-readback.json").write_text(
                        json.dumps(readback, indent=2) + "\n", encoding="utf-8")
                    page.screenshot(path=str(run_artifacts / "unknown-native-export-1024.png"),
                                    full_page=True, animations="disabled")
                    (run_artifacts / "final-account.json").write_text(
                        json.dumps(export, indent=2) + "\n", encoding="utf-8")
                finally:
                    if page is not None:
                        try:
                            (run_artifacts / "last-account.json").write_text(
                                json.dumps(state(page), indent=2) + "\n", encoding="utf-8")
                            page.screenshot(path=str(run_artifacts / "last-browser.png"),
                                            full_page=True, animations="disabled")
                        except Exception as failure:
                            (run_artifacts / "last-capture-error.txt").write_text(
                                str(failure) + "\n", encoding="utf-8")
                    browser.close()


if __name__ == "__main__":
    unittest.main()
