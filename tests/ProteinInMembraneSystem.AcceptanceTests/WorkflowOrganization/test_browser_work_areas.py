"""Focused work-area and feedback checks.

Controlled accounts test presentation and command timing. The separate upload
case crosses the published host and real one-request worker in a fresh workspace.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out" / "workflow-review" / "final"
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import CHROMIUM, account, controlled_account_server  # noqa: E402
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/ExportCompletedMinimizedStage"))
from test_browser_export import ControlledExport, selected_stage_account  # noqa: E402
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein"))
from test_browser_route import running_host, two_alanines  # noqa: E402


def open_browser(playwright, width=1672, height=941):
    browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"])
    return browser, browser.new_page(viewport={"width": width, "height": height},
                                     accept_downloads=True)


def area(page, label):
    page.get_by_role("navigation", name="Research work areas")\
        .get_by_role("button", name=label, exact=True).click()
    expect(page.get_by_role("navigation", name="Research work areas")
           .get_by_role("button", name=label, exact=True)).to_have_attribute("aria-current", "page")


class WorkAreaBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_navigation_preserves_drafts_subject_and_running_progress(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            initial = account()
            initial["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                  "summary": "Controlled assessed protein", "atomCount": 200,
                                  "changes": [], "findings": [], "candidateId": None,
                                  "prediction": None, "geometry": None, "sourceGeometry": None}
            initial["membrane"] = {"modelId": "membrane-one", "status": "assessed",
                                   "upper": [{"speciesId": "DMPC", "fraction": 1.0}],
                                   "lower": [{"speciesId": "DMPC", "fraction": 1.0}],
                                   "scientificPurpose": "Controlled planar bilayer",
                                   "limitations": [], "reason": None,
                                   "policyId": "controlled-policy", "policyVersion": "1",
                                   "evidence": [], "speciesSupport": []}
            initial["placement"] = {"proposalId": "placement-one", "status": "supported",
                                    "preparedProteinId": "protein-one", "membraneModelId": "membrane-one",
                                    "topologyKind": "membrane-spanning", "physicalSide": "both",
                                    "midplaneAngstrom": 0, "thicknessAngstrom": 30,
                                    "depthAngstrom": 0, "tiltDegrees": 0, "sidedness": None,
                                    "contactingRegions": [], "limitations": [],
                                    "policyId": "controlled-policy", "policyVersion": "1",
                                    "witnessId": "controlled-witness", "reason": "Controlled presentation fixture",
                                    "evidence": [], "prediction": None}
            initial["study"]["adoptedPlacementProposalId"] = "placement-one"
            initial["actions"].extend([
                {"kind": "startPreparation", "subjectId": None, "enabled": False,
                 "reason": "The identified attempt is already running."},
                {"kind": "stopAttempt", "subjectId": None, "enabled": True,
                 "reason": None},
            ])
            host.replace(initial)
            browser, page = open_browser(playwright, 820, 760)
            try:
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Preparation")
                expect(page.locator(".attempt-running")).to_contain_text("Minimization in progress")
                expect(page.locator(".attempt-account progress")).to_have_count(1)
                page.screenshot(path=str(ARTIFACTS / "running-preparation-820x760.png"),
                                animations="disabled")

                area(page, "Protein")
                page.locator("#source-query").fill("draft protein query")
                self.assertEqual(page.locator(".rail-section:visible h3").all_text_contents(),
                                 ["Structural source", "Protein model"])
                area(page, "Membrane")
                page.get_by_label("Upper leaflet percentage 1").fill("65")
                self.assertEqual(page.locator(".rail-section:visible h3").all_text_contents(),
                                 ["Membrane ready for placement", "Membrane composition draft"])
                area(page, "Placement")
                expect(page.locator(".rail-section:visible h3"))\
                    .to_have_text("Position the prepared protein")
                area(page, "Preparation")
                expect(page.locator("#preparation-workflow .prerequisite"))\
                    .to_contain_text("already in progress")
                area(page, "Protein")
                expect(page.locator("#source-query")).to_have_value("draft protein query")
                area(page, "Membrane")
                expect(page.get_by_label("Upper leaflet percentage 1")).to_have_value("65")
                self.assertTrue(all(command == "selectInspectionSubject" for command in host.commands),
                                "Navigation and draft editing must not submit scientific work")
                current = page.request.get(base + "/api/state").json()
                self.assertEqual(current["revision"], initial["revision"])
                self.assertEqual(current["inspection"]["subjectId"], "constructed-one")
                expect(page.locator(".header-context")).to_contain_text("Membrane composition draft")
                self.assertLessEqual(page.evaluate("document.body.scrollHeight"), 760)
            finally:
                browser.close()

    def test_search_pending_no_results_and_failure_keep_input(self):
        intercept_search = """
        (() => {
          const nativeFetch = window.fetch.bind(window);
          window.__searchCalls = 0;
          window.fetch = (input, init) => {
            if (String(input) === '/api/commands' && init?.body &&
                JSON.parse(init.body).kind === 'searchSource') {
              window.__searchCalls += 1;
              return new Promise(resolve => {
                window.__respondSearch = (body, status) => resolve(new Response(
                  JSON.stringify(body), {status, headers: {'Content-Type': 'application/json'}}));
              });
            }
            return nativeFetch(input, init);
          };
        })();
        """
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            initial = account()
            initial["attempt"] = None
            initial["inspection"] = None
            initial["actions"].append({"kind": "searchSource", "subjectId": None,
                                       "enabled": True, "reason": None})
            host.replace(initial)
            browser, page = open_browser(playwright)
            try:
                page.add_init_script(intercept_search)
                page.goto(base, wait_until="domcontentloaded")
                area(page, "Protein")
                page.locator("#source-query").fill("unmatched structure")
                page.get_by_role("button", name="Find sources").click()
                expect(page.locator(".action-feedback.pending")).to_contain_text(
                    "Searching structural databases")
                expect(page.get_by_role("button", name="Find sources")).to_be_disabled()
                self.assertEqual(page.evaluate("window.__searchCalls"), 1)
                page.screenshot(path=str(ARTIFACTS / "search-pending-1672x941.png"),
                                animations="disabled")
                area(page, "Membrane")
                page.get_by_label("Upper leaflet percentage 1").fill("65")
                area(page, "Protein")
                result = deepcopy(initial)
                result["revision"] += 1
                page.evaluate("([body, status]) => window.__respondSearch(body, status)",
                              [result, 200])
                expect(page.locator(".action-feedback.warning")).to_contain_text(
                    "No matching structures found")
                expect(page.get_by_role("button", name="Find sources")).to_be_enabled()
                expect(page.locator("#source-query")).to_have_value("unmatched structure")
                area(page, "Membrane")
                expect(page.get_by_label("Upper leaflet percentage 1")).to_have_value("65")
                area(page, "Protein")
                page.get_by_role("button", name="Find sources").click()
                expect(page.locator(".action-feedback.pending")).to_be_visible()
                page.evaluate("([body, status]) => window.__respondSearch(body, status)",
                              [{"reason": "Structural database temporarily unavailable."}, 422])
                expect(page.locator(".action-feedback.error")).to_contain_text(
                    "Structural database temporarily unavailable")
                expect(page.locator("#source-query")).to_have_value("unmatched structure")
                self.assertEqual(page.evaluate("window.__searchCalls"), 2)
                self.assertEqual(host.commands, [], "Controlled timing must not modify the host account")
            finally:
                browser.close()

    def test_completed_stage_guidance_and_export_retry(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = selected_stage_account()
            value["actions"].append({"kind": "requestEquilibration", "subjectId": "stage-one",
                                     "enabled": False,
                                     "reason": "No validated applicable optional procedure is available for this completed minimized stage."})
            host.replace(value)
            controlled = ControlledExport(host)
            controlled.command_failure_once = True
            browser, page = open_browser(playwright)
            try:
                page.route("**/api/commands", controlled.command)
                page.route("**/api/export/*", controlled.download)
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".workspace-context-tabs button.active")).to_have_text("Results")
                expect(page.locator(".header-context")).to_contain_text("Minimized system")
                expect(page.locator(".header-context .context-chip[title='Stage stage-one']"))\
                    .to_contain_text("Minimization")
                expect(page.get_by_role("button", name="Request optional equilibration from this stage"))\
                    .to_be_disabled()
                expect(page.locator(".stage-export .prerequisite")).to_contain_text(
                    "No validated optional equilibration procedure")
                expect(page.locator(".viewer-mount[data-camera-ready='true']"))\
                    .to_be_visible(timeout=120000)
                page.wait_for_timeout(300)
                page.screenshot(path=str(ARTIFACTS / "results-stage-1672x941.png"),
                                animations="disabled")
                page.get_by_role("button", name="Export this completed stage").click()
                expect(page.locator(".stage-export")).to_contain_text("Export not delivered")
                expect(page.locator(".stage-export")).to_contain_text(
                    "Bundle correspondence not verified")
                self.assertEqual(controlled.gets, [])
                page.screenshot(path=str(ARTIFACTS / "export-failure-1672x941.png"),
                                animations="disabled")
                with page.expect_download(timeout=30000) as pending:
                    page.get_by_role("button", name="Export this completed stage").click()
                self.assertEqual(pending.value.suggested_filename,
                                 "protein-membrane-stage-one.zip")
                expect(page.locator(".stage-export")).to_contain_text("Verified bundle received")
                current = page.request.get(base + "/api/state").json()
                self.assertEqual(current["stages"][0]["assessment"]["qualification"],
                                 "indeterminate")
                page.set_viewport_size({"width": 820, "height": 760})
                page.wait_for_timeout(400)
                page.screenshot(path=str(ARTIFACTS / "results-stage-820x760.png"),
                                animations="disabled")
            finally:
                browser.close()

    def test_unsupported_membrane_guidance_points_to_its_own_reason(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            value["attempt"] = None
            value["inspection"] = None
            value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                "summary": "Prepared protein assessed", "atomCount": 23,
                                "changes": [], "findings": [], "candidateId": None,
                                "prediction": None, "geometry": None, "sourceGeometry": None}
            value["membrane"] = {"modelId": "membrane-one", "status": "notEstablished",
                                 "upper": [{"speciesId": "X", "fraction": 1.0}],
                                 "lower": [{"speciesId": "X", "fraction": 1.0}],
                                 "scientificPurpose": "Inspect an exact unsupported choice",
                                 "limitations": [], "reason": "No qualified assets for X",
                                 "policyId": None, "policyVersion": None,
                                 "evidence": [], "speciesSupport": []}
            value["actions"].extend([
                {"kind": "proposePlacement", "subjectId": None, "enabled": False,
                 "reason": "A corresponding assessed membrane model is required."},
                {"kind": "startPreparation", "subjectId": None, "enabled": False,
                 "reason": "A supported adopted placement and exact assets are required."},
                {"kind": "adoptMembrane", "subjectId": None, "enabled": False,
                 "reason": "This model has no qualified lipid assets."},
            ])
            host.replace(value)
            browser, page = open_browser(playwright)
            try:
                page.goto(base, wait_until="domcontentloaded")
                area(page, "Placement")
                expect(page.locator("#placement-workflow .prerequisite"))\
                    .to_contain_text("This membrane composition is unsupported")
                page.locator("#placement-workflow .text-link").click()
                expect(page.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Membrane")
                expect(page.locator("#membrane-workflow")).to_contain_text(
                    "No qualified assets for X")
                expect(page.get_by_role("region", name="Membrane task outcome"))\
                    .to_contain_text("Membrane support not established")
                area(page, "Preparation")
                expect(page.locator("#preparation-workflow .prerequisite"))\
                    .to_contain_text("This membrane composition is unsupported")
                self.assertTrue(all(command == "selectInspectionSubject" for command in host.commands))
            finally:
                browser.close()

    def test_real_uploaded_source_and_protein_assessment(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser, page = open_browser(playwright)
                try:
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.locator(".workspace-context-tabs button.active"))\
                        .to_have_text("Protein")
                    area(page, "Protein")
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text(
                        "Structure displayed", timeout=120000)
                    expect(page.get_by_role("button", name="Inspect source structure")).to_have_count(0)
                    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    expect(page.get_by_role("region", name="Protein task outcome")).to_contain_text(
                        "Protein prepared", timeout=120000)
                    state = page.request.get(base + "/api/state").json()
                    self.assertEqual(state["protein"]["status"], "assessed")
                    self.assertEqual(state["stages"], [])
                    page.get_by_role("button", name="View prepared protein").click()
                    expect(page.get_by_role("heading", name="Prepared protein · two-alanines.pdb"))\
                        .to_be_visible()
                    expect(page.locator(".evidence-panel")).to_contain_text(
                        "Preparation complete")
                    area(page, "Membrane")
                    expect(page.locator(".header-context")).to_contain_text(
                        "Membrane composition draft")
                finally:
                    browser.close()

    def test_area_view_restores_source_camera_and_empty_membrane_is_not_a_protein(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser, page = open_browser(playwright, 1024, 800)
                try:
                    page.goto(base, wait_until="domcontentloaded")
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    expect(page.locator(".chain-color-swatch")).to_have_count(1)
                    expect(page.get_by_role("button", name="Focus in viewer")).to_be_visible()
                    page.get_by_role("button", name="Focus in viewer").click()
                    page.wait_for_timeout(700)
                    before = page.evaluate("""() => {
                      const camera = document.querySelector('.viewer-mount')
                        [Symbol.for('molstar.viewer')].plugin.canvas3d.camera;
                      const snapshot = camera.getSnapshot();
                      camera.setState({position: [snapshot.position[0] + 4,
                        snapshot.position[1] + 3, snapshot.position[2] + 2]}, 0);
                      return [...camera.getSnapshot().position];
                    }""")
                    page.wait_for_timeout(120)
                    revision = page.request.get(base + "/api/state").json()["revision"]
                    area(page, "Membrane")
                    expect(page.locator(".scene-heading h1")).to_have_text("Membrane composition draft")
                    expect(page.locator(".bilayer-unassigned")).to_have_count(2)
                    expect(page.locator(".bilayer-unassigned").first).to_have_text("Not specified")
                    expect(page.locator(".viewer-mount")).to_have_count(0)
                    self.assertEqual(page.request.get(base + "/api/state").json()["revision"], revision)
                    area(page, "Protein")
                    expect(page.locator(".scene-heading h1")).to_contain_text("Source structure")
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=120000)
                    page.wait_for_timeout(900)
                    after = page.evaluate("""() => [...document.querySelector('.viewer-mount')
                      [Symbol.for('molstar.viewer')].plugin.canvas3d.camera.getSnapshot().position]""")
                    self.assertLess(max(abs(a - b) for a, b in zip(before, after)), 1e-5)
                    self.assertFalse(page.evaluate("document.documentElement.scrollWidth > window.innerWidth"))
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
