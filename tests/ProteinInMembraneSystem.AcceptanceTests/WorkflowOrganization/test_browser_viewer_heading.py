"""Viewer-heading presentation checks with identified, non-scientific fixtures.

The source and stage accounts below only verify wording and subject continuity;
their small coordinate fixture does not establish a prepared or minimized result.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out" / "ui-screenshots" / "viewer-heading" / "fixtures"
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import account, completed, controlled_account_server, inspection  # noqa: E402
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/WorkflowOrganization"))
from test_browser_connected_source_inspection import (  # noqa: E402
    open_browser, small_alanine_chain, source_account,
)


class ViewerHeadingBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_source_heading_tracks_exact_displayed_source_and_origin(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.structure_bytes = small_alanine_chain()
            host.replace(source_account("rcsb:1AAA", 11))
            browser, page = open_browser(playwright, 1280, 800)
            try:
                page.goto(base, wait_until="domcontentloaded")
                heading = page.locator(".scene-heading h1")
                subtitle = page.locator(".scene-subtitle")
                expect(heading).to_have_text("Source structure · 1AAA", timeout=60000)
                expect(subtitle).to_contain_text("RCSB PDB · Study revision 2 · current study")
                self.assertNotIn("Connected structural inspection", page.locator(".scene-heading").inner_text())
                self.assertNotIn("Experimental", subtitle.inner_text(),
                                 "An exact RCSB reference alone does not establish experimental origin")
                page.get_by_role("button", name="Collapse inputs").click()
                page.set_viewport_size({"width": 900, "height": 740})
                expect(heading).to_be_visible()
                self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth"), 900)
                page.screenshot(path=str(ARTIFACTS / "source-small-collapsed.png"))

                experimental = source_account("rcsb:1AAA", 11)
                experimental["sourceCandidates"][0]["provenance"] = \
                    "RCSB PDB experimental entry search"
                host.replace(experimental)
                expect(subtitle).to_contain_text("RCSB PDB · experimental entry")

                predicted = source_account("alphafold:AF-A0A6M8B2E4-F1", 12)
                predicted["study"]["selectedSourceKind"] = "alphafold"
                predicted["inspection"]["structureUrl"] += "&source=predicted"
                host.replace(predicted)
                expect(heading).to_have_text("Source structure")
                expect(subtitle).to_contain_text("AlphaFold DB prediction")
                expect(page.locator(".selected-source-account")).to_contain_text(
                    "Prediction coordinates; not experimental validation")
                page.screenshot(path=str(ARTIFACTS / "predicted-source-small-collapsed.png"))

                uploaded = source_account("upload:fixture-one", 13)
                uploaded["study"].update(selectedSourceKind="upload", uploadProvenance="experimental")
                uploaded["inspection"]["structureUrl"] += "&source=upload"
                host.replace(uploaded)
                expect(heading).to_have_text("Source structure")
                expect(subtitle).to_contain_text("Researcher upload · researcher-declared experimental origin")
                expect(page.locator(".selected-source-account")).to_contain_text(
                    "this label is not independently verified")
                self.assertEqual(host.commands, [], "Heading changes cannot submit scientific work")
            finally:
                browser.close()

    def test_prepared_constructed_and_historical_stage_headings_follow_inspection(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.structure_bytes = small_alanine_chain()
            prepared = account()
            prepared["attempt"] = None
            prepared["inspection"] = inspection("prepared-one", "preparedProtein")
            prepared["inspection"].update(structureUrl="/api/structures/tiny?format=pdb&subject=prepared",
                                           omittedMolecules=[])
            prepared["protein"] = {
                "subjectId": "prepared-one", "status": "prepared", "summary": "Controlled protein fixture",
                "atomCount": 41, "changes": [], "findings": [], "candidateId": None,
                "prediction": None, "geometry": None, "sourceGeometry": None,
            }
            host.replace(prepared)
            browser, page = open_browser(playwright, 1280, 800)
            try:
                page.goto(base, wait_until="domcontentloaded")
                heading = page.locator(".scene-heading h1")
                subtitle = page.locator(".scene-subtitle")
                expect(heading).to_have_text("Prepared protein", timeout=60000)
                expect(subtitle).to_contain_text("Assessed prepared coordinates · Study revision 1 · current study")
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                page.wait_for_function("""() => !!Reflect.get(document.querySelector('.viewer-mount'),
                  Symbol.for('molstar.viewer'))?.plugin.managers.structure.hierarchy.current.structures[0]
                  ?.cell.obj?.data""", timeout=60000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                page.wait_for_timeout(200)
                page.screenshot(path=str(ARTIFACTS / "prepared-protein-expanded.png"))

                construction = completed(account())
                construction["study"].update(id="revision-two", number=2)
                construction["inspection"].update(
                    structureUrl="/api/structures/tiny?format=pdb&subject=constructed",
                    omittedMolecules=[])
                construction["stages"][0]["assessment"]["currentlyApplicable"] = False
                host.replace(construction)
                page.get_by_role("button", name="Preparation", exact=True).click()
                expect(heading).to_have_text("Protein–membrane system")
                expect(subtitle).to_contain_text("Checked constructed candidate; minimization not complete")
                expect(subtitle).to_contain_text("Study revision 1 · historical revision")
                expect(page.locator(".execution-scene-label")).to_contain_text(
                    "Protein–membrane system · historical revision")
                page.wait_for_timeout(350)
                page.screenshot(path=str(ARTIFACTS / "historical-constructed-expanded.png"))

                minimized = deepcopy(construction)
                minimized["revision"] += 1
                minimized["inspection"] = inspection("stage-one", "completedStage",
                                                      minimized["stages"][0]["assessment"])
                minimized["inspection"].update(
                    structureUrl="/api/structures/tiny?format=pdb&subject=minimized",
                    omittedMolecules=[])
                host.replace(minimized)
                page.get_by_role("button", name="Results", exact=True).click()
                expect(heading).to_have_text("Minimized system")
                expect(subtitle).to_contain_text("Completed molecular result · Study revision 1 · historical revision")
                expect(page.locator(".execution-scene-label")).to_contain_text(
                    "Minimized system · historical revision")
                page.get_by_role("button", name="Collapse inputs").click()
                page.set_viewport_size({"width": 900, "height": 740})
                expect(page.locator(".execution-scene-label")).to_be_visible()
                self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth"), 900)
                page.screenshot(path=str(ARTIFACTS / "historical-minimized-small-collapsed.png"))

                equilibrated = deepcopy(minimized)
                equilibrated["revision"] += 1
                equilibrated["stages"][0].update(stageId="stage-two", kind="Equilibration",
                                                  sourceStageId="stage-one",
                                                  summary="Controlled equilibrated-stage fixture")
                equilibrated["stages"][0]["assessment"]["stageId"] = "stage-two"
                equilibrated["inspection"] = inspection("stage-two", "completedStage",
                                                         equilibrated["stages"][0]["assessment"])
                equilibrated["inspection"].update(
                    structureUrl="/api/structures/tiny?format=pdb&subject=equilibrated",
                    omittedMolecules=[])
                host.replace(equilibrated)
                expect(heading).to_have_text("Equilibrated system")
                expect(page.locator(".execution-scene-label")).to_contain_text(
                    "Equilibrated system · historical revision")

                opaque = deepcopy(equilibrated)
                opaque_id = "stage-" + "a" * 40
                opaque["revision"] += 1
                opaque["stages"][0]["stageId"] = opaque_id
                opaque["stages"][0]["assessment"]["stageId"] = opaque_id
                opaque["inspection"] = inspection(opaque_id, "completedStage",
                                                    opaque["stages"][0]["assessment"])
                opaque["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb&subject=opaque"
                host.replace(opaque)
                expect(page.locator(".execution-scene-label")).to_contain_text("Equilibrated system")
                expect(page.locator(".evidence-panel")).to_contain_text(opaque_id)
                self.assertTrue(page.locator(".execution-scene-label").evaluate(
                    "element => element.scrollWidth <= element.clientWidth + 1"),
                    "The compact viewer label must fit without raw IDs")
                self.assertTrue(all(kind == "selectInspectionSubject" for kind in host.commands),
                                "Presentation fixtures cannot submit scientific work")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
