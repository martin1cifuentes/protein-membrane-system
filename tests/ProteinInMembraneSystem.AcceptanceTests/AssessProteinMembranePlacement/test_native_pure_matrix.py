"""Each installed pure Lipid21 patch through a real isolated native builder.

The non-whitelisted two-alanine input is a technical software fixture.
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

from test_manual_position import POLICY, chromium, prepare_with_current_plan, running_host, state, two_alanines, until


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out/reconciled-design/native-pure-full"


class NativePurePatchTests(unittest.TestCase):
    def run_species(self, species: str, position: str):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-software-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace", policy=POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1024, "height": 768}, accept_downloads=True)
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
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option(species)
                    page.get_by_role("button", name="Use this membrane").click()
                    until(page, lambda account: (account.get("membrane") or {}).get("status") == "assessed",
                          f"checked {species}")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Placement").click()
                    if position != "center":
                        page.locator("#starting-position").select_option(position)
                    positioned = until(page, lambda account: (account.get("placement") or {}).get("status") ==
                                       "supported" and (account["placement"].get("transform") or {})
                                       .get("startingPosition") == position, f"{species} checked {position} pose")
                    page.get_by_role("button", name="Use this position").click()
                    until(page, lambda account: account["study"]["adoptedPlacementProposalId"] ==
                          positioned["placement"]["proposalId"], "adopted exact pose")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Preparation").click()
                    route = page.locator(".construction-route").filter(
                        has_text=f"OpenMM native: {species}")
                    expect(route.get_by_role("button", name="Build and minimize")).to_be_enabled()
                    route.get_by_role("button", name="Build and minimize").click()
                    candidate = until(page, lambda account: (account.get("attempt") or {}).get("constructed")
                                      is not None or (account.get("attempt") or {}).get("status") in
                                      ("failed", "resourceRefused", "unobserved"),
                                      f"real {species} native construction", seconds=1100)
                    self.assertEqual(candidate["attempt"]["policyId"],
                                     f"canonical-protein-pure-{species.lower()}-native-construction")
                    constructed = candidate["attempt"]["constructed"]
                    self.assertIsNotNone(constructed)
                    self.assertEqual({item["speciesId"] for item in constructed["achievedComposition"]},
                                     {species})
                    self.assertEqual({item["physicalSide"] for item in constructed["achievedComposition"]},
                                     {"upper", "lower"})
                    self.assertGreater(constructed["waterCount"], 0)
                    self.assertGreater(constructed["atomCount"], candidate["protein"]["atomCount"])
                    self.assertTrue(all(0 < length <= 180 for length in constructed["cellAngstrom"]))
                    page.get_by_role("button", name="Review attempt").click()
                    page.get_by_role("button", name="Inspect verified constructed system").click()
                    until(page, lambda account: (account.get("inspection") or {}).get("subjectId") ==
                          constructed["subjectId"], f"{species} exact candidate inspection")
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    ARTIFACTS.mkdir(parents=True, exist_ok=True)
                    page.screenshot(path=str(ARTIFACTS / f"{species.lower()}-candidate-1024.png"),
                                    full_page=True, animations="disabled")
                    if candidate["attempt"]["status"] != "completed":
                        self.assertFalse(any(item["status"] == "completed" for item in candidate["stages"]))
                    completed = until(page, lambda account: (account.get("attempt") or {}).get("status") in
                                      ("completed", "failed", "stopped", "unobserved") and
                                      (account.get("attempt") or {}).get("stageKind") == "Minimization",
                                      f"{species} exact-candidate minimization", seconds=1800)
                    self.assertEqual(completed["attempt"]["status"], "completed", completed["attempt"])
                    minimized = next(item for item in completed["stages"]
                                     if item["attemptId"] == candidate["attempt"]["attemptId"] and
                                     item["kind"] == "Minimization")
                    self.assertEqual(minimized["constructed"]["subjectId"], constructed["subjectId"])
                    self.assertEqual(minimized["observation"]["termination"], "converged")
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Results").click()
                    page.locator(".result-stage-list button").first.click()
                    until(page, lambda account: (account.get("inspection") or {}).get("subjectId") ==
                          minimized["stageId"], f"{species} minimized-stage inspection")
                    expect(page.locator(".scene-loading")).not_to_be_visible(timeout=240000)
                    with page.expect_download(timeout=300000) as delivered:
                        page.get_by_role("button", name="Export this completed stage").click()
                    bundle = ARTIFACTS / f"{species.lower()}-minimized-export.zip"
                    delivered.value.save_as(str(bundle))
                    exported = until(page, lambda account: (next((item for item in account["stages"]
                          if item["stageId"] == minimized["stageId"]), {}).get("export") or {})
                          .get("status") == "verified", f"{species} exact-stage export")
                    stage = next(item for item in exported["stages"]
                                 if item["stageId"] == minimized["stageId"])
                    self.assertEqual(hashlib.sha256(bundle.read_bytes()).hexdigest(),
                                     stage["export"]["sha256"])
                    with zipfile.ZipFile(bundle) as archive:
                        self.assertIsNone(archive.testzip())
                        self.assertGreaterEqual(set(archive.namelist()),
                                                {"manifest.json", "structure.cif", "topology.json",
                                                 "system.xml", "state.xml"})
                        manifest = json.loads(archive.read("manifest.json"))
                        self.assertEqual(manifest["stage"]["id"], minimized["stageId"])
                        self.assertEqual(manifest["attempt"]["id"], candidate["attempt"]["attemptId"])
                        self.assertEqual(manifest["stage"]["atomCount"], constructed["atomCount"])
                        self.assertEqual(manifest["lineage"]["sourceCoordinateSha256"],
                                         hashlib.sha256(source.read_bytes()).hexdigest())
                    readback = subprocess.run([str(ROOT / "out/python/bin/python"), "-c", """
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
                    self.assertEqual(readback.returncode, 0, readback.stderr)
                    self.assertEqual(int(readback.stdout.strip()), constructed["atomCount"])
                    page.screenshot(path=str(ARTIFACTS / f"{species.lower()}-export-1024.png"),
                                    full_page=True, animations="disabled")
                finally:
                    browser.close()

    def test_dlpc_upper(self):
        self.run_species("DLPC", "upper")

    def test_dlpe_center(self):
        self.run_species("DLPE", "center")

    @unittest.skip("The installed native DOPC patch has two trans alkenes and is not enabled")
    def test_dopc_lower(self):
        self.run_species("DOPC", "lower")

    @unittest.skip("The installed native DPPC patch has inverted glycerol centres and is not enabled")
    def test_dppc_center(self):
        self.run_species("DPPC", "center")

    @unittest.skip("POPC uses the separately identified mapped Zenodo route")
    def test_popc_upper(self):
        self.run_species("POPC", "upper")

    def test_pope_center(self):
        self.run_species("POPE", "center")


if __name__ == "__main__":
    unittest.main()
