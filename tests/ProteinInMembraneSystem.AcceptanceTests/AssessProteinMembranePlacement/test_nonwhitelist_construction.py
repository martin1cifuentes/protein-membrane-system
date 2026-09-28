"""Real OpenMM native DMPC candidate for a non-whitelisted software construct.

The two-alanine PDB is deliberately small and is not biological suitability
evidence. This route checks actual provider, identity and candidate handoffs.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import subprocess
import tempfile
import unittest
import zipfile

from playwright.sync_api import expect, sync_playwright

from test_manual_position import (ARTIFACTS, POLICY, chromium, prepare_with_current_plan,
                                  running_host, state, two_alanines, until)


class NativeConstructionBrowserTests(unittest.TestCase):
    def test_nonwhitelisted_construct_reaches_exact_candidate(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-software-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace", policy=POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940}, accept_downloads=True)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    prepare_with_current_plan(page)
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Membrane").click()
                    for side in ("Upper", "Lower"):
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                    page.get_by_role("button", name="Propose membrane model").click()
                    until(page, lambda current: (current.get("membrane") or {}).get("status") == "proposed",
                          "proposed DMPC")
                    page.get_by_role("button", name="Adopt and assess displayed proposal").click()
                    until(page, lambda current: (current.get("membrane") or {}).get("status") == "assessed",
                          "checked DMPC")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Placement").click()
                    page.locator("#starting-position").select_option("upper")
                    positioned = until(page, lambda current: (current.get("placement") or {}).get("status") ==
                                       "supported" and (current["placement"].get("transform") or {})
                                       .get("startingPosition") == "upper", "checked upper-side pose")
                    page.get_by_role("button", name="Use this position").click()
                    until(page, lambda current: current["study"]["adoptedPlacementProposalId"] ==
                          positioned["placement"]["proposalId"], "adopted pose")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Preparation").click()
                    expect(page.get_by_role("button", name="Construct system")).to_be_enabled()
                    page.get_by_role("button", name="Construct system").click()
                    candidate = until(page, lambda current: (current.get("attempt") or {}).get("status") in
                                      ("readyForMinimization", "failed", "resourceRefused", "unobserved"),
                                      "real native candidate or attributable failure", seconds=1100)
                    self.assertEqual(candidate["attempt"]["status"], "readyForMinimization",
                                     candidate["attempt"])
                    self.assertIsNotNone(candidate["attempt"]["constructed"])
                    self.assertGreater(candidate["attempt"]["constructed"]["atomCount"],
                                       candidate["protein"]["atomCount"])
                    until(page, lambda current: (current.get("inspection") or {}).get("subjectId") ==
                          candidate["attempt"]["constructed"]["subjectId"], "candidate inspection")
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    page.screenshot(path=str(ARTIFACTS / "nonwhitelist-native-candidate-1672.png"),
                                    full_page=True, animations="disabled")
                    continue_button = page.get_by_role("button", name="Authorize minimization of this candidate")
                    expect(continue_button).to_be_enabled()
                    continue_button.click()
                    completed = until(page, lambda current: (current.get("attempt") or {}).get("status") in
                                      ("completed", "failed", "stopped", "unobserved") and
                                      (current.get("attempt") or {}).get("stageKind") == "Minimization",
                                      "real minimization outcome", seconds=1800)
                    self.assertEqual(completed["attempt"]["status"], "completed", completed["attempt"])
                    minimized = next(item for item in completed["stages"]
                                     if item["attemptId"] == candidate["attempt"]["attemptId"] and
                                     item["kind"] == "Minimization")
                    self.assertEqual(minimized["status"], "completed")
                    self.assertEqual(minimized["constructed"]["subjectId"],
                                     candidate["attempt"]["constructed"]["subjectId"])
                    self.assertEqual(minimized["observation"]["termination"], "converged")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Results").click()
                    until(page, lambda current: (current.get("inspection") or {}).get("subjectId") ==
                          minimized["stageId"], "minimized stage inspection")
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    page.screenshot(path=str(ARTIFACTS / "nonwhitelist-minimized-1672.png"),
                                    full_page=True, animations="disabled")
                    with page.expect_download(timeout=300000) as delivered:
                        page.get_by_role("button", name="Export this completed stage").click()
                    bundle = ARTIFACTS / "nonwhitelist-minimized-export.zip"
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
                        self.assertEqual(manifest["attempt"]["id"], candidate["attempt"]["attemptId"])
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
                    scientific_readback = subprocess.run([str(Path(__file__).resolve().parents[3] /
                        "out" / "python" / "bin" / "python"), "-c", """
import io, sys, zipfile
from openmm import XmlSerializer
from openmm.app import PDBxFile
with zipfile.ZipFile(sys.argv[1]) as bundle:
    structure = PDBxFile(io.StringIO(bundle.read('structure.cif').decode()))
    system = XmlSerializer.deserialize(bundle.read('system.xml').decode())
    state = XmlSerializer.deserialize(bundle.read('state.xml').decode())
    assert structure.topology.getNumAtoms() == system.getNumParticles()
    assert len(state.getPositions()) == system.getNumParticles()
    assert structure.topology.getPeriodicBoxVectors() is not None
print(system.getNumParticles())
""", str(bundle)], capture_output=True, text=True, timeout=90)
                    self.assertEqual(scientific_readback.returncode, 0, scientific_readback.stderr)
                    self.assertEqual(int(scientific_readback.stdout.strip()),
                                     candidate["attempt"]["constructed"]["atomCount"])
                    page.screenshot(path=str(ARTIFACTS / "nonwhitelist-export-1672.png"),
                                    full_page=True, animations="disabled")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
