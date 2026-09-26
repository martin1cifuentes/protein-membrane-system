"""Controlled presentation of a distinct optional stage and its exact predecessor.

The fixture exercises actor presentation and navigation. A published scientific
route and class-specific qualification are separate Slice 7 acceptance claims.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "ExportCompletedMinimizedStage"))
from test_browser_export import ControlledExport, selected_stage_account  # noqa: E402
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "InspectAndReattach"))
from test_browser_presentation import CHROMIUM, DIST, controlled_account_server, inspection  # noqa: E402


CAPTURES = ROOT / "out" / "browser-acceptance" / "slice7"


def equilibrated_account() -> dict:
    value = selected_stage_account()
    minimized = value["stages"][0]
    minimized["sourceStageId"] = None
    value["stages"] = [minimized]
    equilibrated = deepcopy(minimized)
    equilibrated.update(
        stageId="stage-equilibrated", kind="Equilibration", sourceStageId="stage-one",
        summary="Declared optional procedure completed; assessment remains distinct",
        assessment={"id": "assessment-stage-equilibrated",
                    "stageId": "stage-equilibrated", "qualification": "indeterminate",
                    "reason": "Required independent observations remain unavailable",
                    "currentlyApplicable": True, "evidence": [], "findings": [],
                    "limitations": ["Controlled presentation only"]},
    )
    value["stages"].append(equilibrated)
    value["attempt"].update(currentStageId=equilibrated["stageId"], stageKind="Equilibration",
                            status="completed", message="Observed optional procedure completion")
    value["inspection"] = inspection(equilibrated["stageId"], "completedStage",
                                     equilibrated["assessment"])
    value["inspection"].update(structureUrl="/api/structures/tiny?format=pdb",
                               omittedMolecules=[])
    value["actions"] = [item for item in value["actions"] if item["kind"] != "exportStage"]
    value["actions"].extend([
        {"kind": "exportStage", "subjectId": minimized["stageId"], "enabled": False,
         "reason": "Select this completed stage before export."},
        {"kind": "exportStage", "subjectId": equilibrated["stageId"], "enabled": True,
         "reason": None},
    ])
    return value


class EquilibratedStageBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (DIST / "index.html").is_file():
            raise AssertionError("Build browser/ before Slice 7 controlled presentation")

    def test_exact_minimized_predecessor_and_equilibrated_status(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(equilibrated_account())
            controlled = ControlledExport(host)
            browser = playwright.chromium.launch(
                executable_path=CHROMIUM, headless=True,
                args=["--disable-dev-shm-usage", "--use-gl=angle", "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.route("**/api/commands", controlled.command)
                page.goto(base, wait_until="domcontentloaded")
                review = page.get_by_label("Equilibration stage review and distinct scientific assessment")
                expect(review).to_contain_text("Selected stage")
                expect(review).to_contain_text("Equilibrated")
                expect(review).to_contain_text("Earlier stage")
                expect(review).to_contain_text("Minimized — available")
                expect(review).to_contain_text("Optional procedure")
                expect(review).to_contain_text("Completed")
                expect(review).to_contain_text("indeterminate")
                expect(review).to_contain_text(
                    "Procedure completion and scientific qualification are separate")
                expect(page.locator(".stage-strip")).to_contain_text("Minimized · Equilibrated")
                tabs = page.get_by_role("navigation", name="Jump to researcher workflow")
                expect(tabs.get_by_role("button", name="Minimized stage review"))\
                    .to_be_visible()
                expect(tabs.get_by_role("button", name="Equilibrated stage review"))\
                    .to_have_attribute("aria-current", "page")
                expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
                select_minimized = page.get_by_role("button", name="Select minimized stage")
                expect(select_minimized).to_be_enabled()
                CAPTURES.mkdir(parents=True, exist_ok=True)
                for width, height in ((1672, 941), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1)
                    expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                    expect(select_minimized).to_be_visible()
                    page.screenshot(path=str(CAPTURES / f"controlled-equilibrated-review-{width}.png"),
                                    full_page=True, animations="disabled")
                select_minimized.click()
                expect(page.get_by_label("Minimized stage review and distinct scientific assessment"))\
                    .to_contain_text("Required evidence unavailable")
                self.assertEqual(controlled.commands[-1],
                                 ("selectInspectionSubject", "stage-one"))
                self.assertEqual(page.request.get(base + "/api/state").json()["inspection"]["subjectId"],
                                 "stage-one")
                expect(tabs.get_by_role("button", name="Minimized stage review"))\
                    .to_have_attribute("aria-current", "page")
                expect(page.locator(".stage-strip")).to_contain_text("Minimized · Equilibrated")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
