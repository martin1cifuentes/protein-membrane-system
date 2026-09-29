"""Focused work-area and feedback checks.

Controlled accounts test presentation and command timing. The separate upload
case crosses the published host and real one-request worker in a fresh workspace.
"""

from __future__ import annotations

from copy import deepcopy
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import tempfile
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out" / "workflow-review" / "final"
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import CHROMIUM, account, controlled_account_server, inspection  # noqa: E402
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/ExportCompletedMinimizedStage"))
from test_browser_export import ControlledExport, selected_stage_account  # noqa: E402
_select_route = ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein/test_browser_route.py"
_select_spec = spec_from_file_location("select_prepare_browser_route", _select_route)
assert _select_spec is not None and _select_spec.loader is not None
_select_module = module_from_spec(_select_spec)
_select_spec.loader.exec_module(_select_module)
running_host, two_alanines = _select_module.running_host, _select_module.two_alanines


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

    def test_diagnostic_download_requires_opaque_token_and_same_origin(self):
        """Published Host endpoint boundary; root owner tests bind and tamper a real token."""
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            with running_host(directory / "workspace") as base:
                def refused(path: str, headers: dict[str, str], expected: int):
                    with self.assertRaises(HTTPError) as failure:
                        urlopen(Request(base + path, headers=headers), timeout=5)
                    self.assertEqual(failure.exception.code, expected)

                refused("/api/diagnostics/unknown-token", {}, 404)
                refused("/api/diagnostics/unknown-token", {"Origin": "http://example.test"}, 403)
                refused("/api/diagnostics/unknown-token",
                        {"Origin": base, "Sec-Fetch-Site": "cross-site"}, 403)
                refused("/api/diagnostics/bad%3Atoken", {"Origin": base}, 404)

    def test_terminal_pre_admission_refusal_has_no_phantom_attempt(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            value["inspection"] = None
            value["attempt"].update(attemptId="", status="failed", stageKind=None,
                                    progress=None, message="The exact inputs cannot be bound.",
                                    studyRevisionId=None, policyId=None, policyVersion=None,
                                    constructed=None, failureCode="preAdmissionRefused")
            value["actions"].append({"kind": "stopAttempt", "subjectId": None,
                                     "enabled": False, "reason": "No attempt is running."})
            host.replace(value)
            browser, page = open_browser(playwright, 820, 760)
            try:
                page.goto(base, wait_until="domcontentloaded")
                area(page, "Preparation")
                expect(page.locator(".attempt-account")).to_contain_text(
                    "Construction request not admitted")
                expect(page.locator(".attempt-account")).to_contain_text(
                    "The exact inputs cannot be bound.")
                expect(page.locator(".attempt-account")).to_contain_text(
                    "correct the inputs or choose an available recipe")
                expect(page.locator(".attempt-running")).to_have_count(0)
                expect(page.get_by_role("button", name="Review attempt")).to_have_count(0)
                expect(page.get_by_role("button", name="Stop unfinished work")).to_be_disabled()
                expect(page.locator(".stage-strip")).to_contain_text("No attempt admitted")
                expect(page.locator(".header-context")).to_contain_text("Build request")
                page.locator(".attempt-account details").first.locator("summary").click()
                expect(page.locator(".attempt-account details")).to_contain_text("No attempt admitted")
                self.assertNotIn("Attempt ", page.locator(".attempt-account details p").inner_text())
                page.screenshot(path=str(ARTIFACTS / "pre-admission-refused-820x760.png"),
                                full_page=True, animations="disabled")

                refused = deepcopy(value)
                refused["revision"] += 1
                refused["attempt"].update(status="resourceRefused",
                                           message="The bounded memory budget is unavailable.",
                                           failureCode="resourceLimit")
                host.replace(refused)
                page.reload(wait_until="domcontentloaded")
                area(page, "Preparation")
                expect(page.locator(".attempt-account")).to_contain_text(
                    "The bounded memory budget is unavailable.")
                expect(page.locator(".stage-strip")).to_contain_text("No attempt admitted")
                expect(page.get_by_role("button", name="Stop unfinished work")).to_be_disabled()
            finally:
                browser.close()

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
                {"kind": "buildAndMinimize", "subjectId": None, "enabled": False,
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
                expect(page.locator(".attempt-running")).to_contain_text("Final minimization in progress")
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
                self.assertEqual(current["stages"][0]["assessment"]["checkStanding"],
                                 "checksIncomplete")
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
                {"kind": "buildAndMinimize", "subjectId": None, "enabled": False,
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
                    authorize = page.get_by_role("button", name="Prepare selected protein", exact=True).or_(
                        page.get_by_role("button", name="Prepare with recommendations", exact=True))
                    expect(authorize).to_be_visible(timeout=120000)
                    authorize.click()
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
                    page.evaluate("""() => {
                      const camera = document.querySelector('.viewer-mount')
                        [Symbol.for('molstar.viewer')].plugin.canvas3d.camera;
                      const snapshot = camera.getSnapshot();
                      camera.setState({position: [snapshot.position[0] + 4,
                        snapshot.position[1] + 3, snapshot.position[2] + 2]}, 1);
                    }""")
                    page.wait_for_timeout(120)
                    before = page.evaluate("""() => [...document.querySelector('.viewer-mount')
                      [Symbol.for('molstar.viewer')].plugin.canvas3d.camera.getSnapshot().position]""")
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

    def test_running_attempt_keeps_placement_and_preparation_cameras_distinct(self):
        """Fresh reconnect: one positioned protein, two honest work-area views."""
        protein_bytes = (ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein" /
                         "fixtures/1CRN-RCSB.pdb").read_bytes()

        def settled_scene(page):
            mount = page.locator(".viewer-mount")
            expect(mount.locator("canvas")).to_have_count(1, timeout=120000)
            page.wait_for_function("""() => {
              const mount = document.querySelector('.viewer-mount');
              const viewer = mount && Reflect.get(mount, Symbol.for('molstar.viewer'));
              const viewport = viewer?.plugin.canvas3d?.camera.viewport;
              const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]
                ?.cell.obj?.data;
              const bounds = mount?.getBoundingClientRect();
              return mount?.dataset.cameraReady === 'true' && structure && viewport && bounds &&
                Math.abs(viewport.width - bounds.width) <= 4 &&
                Math.abs(viewport.height - bounds.height) <= 4;
            }""", timeout=120000)
            page.evaluate("""() => new Promise(resolve => requestAnimationFrame(() =>
              requestAnimationFrame(() => requestAnimationFrame(resolve))))""")
            page.evaluate("window.__pimCameraStability = null")
            page.wait_for_function("""() => {
              const mount = document.querySelector('.viewer-mount');
              const camera = mount && Reflect.get(mount, Symbol.for('molstar.viewer'))
                ?.plugin.canvas3d?.camera;
              if (!camera) return false;
              const snapshot = camera.getSnapshot();
              const values = [...snapshot.target, ...snapshot.position, snapshot.radius];
              const prior = window.__pimCameraStability;
              if (!prior || Math.max(...values.map((value, i) =>
                  Math.abs(value - prior.values[i]))) > 1e-4) {
                window.__pimCameraStability = {values, since: performance.now()};
                return false;
              }
              return performance.now() - prior.since >= 300;
            }""", polling=50, timeout=10000)
            return page.evaluate("""() => {
              const mount = document.querySelector('.viewer-mount');
              const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
              const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0]
                ?.cell.obj?.data;
              const camera = viewer.plugin.canvas3d.camera.getSnapshot();
              if (!structure) throw new Error('The identified protein did not load');
              const box = structure.boundary.box;
              return {subject: mount.getAttribute('aria-label'), atoms: structure.elementCount,
                representations: viewer.plugin.canvas3d.reprCount.value,
                position: [...camera.position], target: [...camera.target], radius: camera.radius,
                distance: Math.hypot(...camera.position.map((value, i) =>
                  value - camera.target[i])),
                targetDistance: viewer.plugin.canvas3d.camera.getTargetDistance(camera.radius),
                center: [0, 1, 2].map(i => (box.min[i] + box.max[i]) / 2),
                zSpan: box.max[2] - box.min[2]};
            }""")

        for width, height, pending_stop in ((1672, 941, False), (1024, 768, False),
                                            (820, 760, True)):
            with self.subTest(viewport=width, pending_stop=pending_stop), \
                    controlled_account_server() as (host, base), sync_playwright() as playwright:
                value = account()
                value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                    "summary": "Controlled positioned protein", "atomCount": 327,
                                    "changes": [], "findings": [], "candidateId": None,
                                    "prediction": None, "geometry": None, "sourceGeometry": None}
                value["membrane"] = {"modelId": "membrane-one", "status": "assessed",
                                     "upper": [{"speciesId": "DMPC", "fraction": 1.0}],
                                     "lower": [{"speciesId": "DMPC", "fraction": 1.0}],
                                     "scientificPurpose": "Controlled intended bilayer",
                                     "limitations": [], "reason": None,
                                     "policyId": "controlled-policy", "policyVersion": "1",
                                     "evidence": [], "speciesSupport": []}
                value["placement"] = {"proposalId": "placement-one", "status": "supported",
                                      "preparedProteinId": "protein-one", "membraneModelId": "membrane-one",
                                      "topologyKind": "membrane-spanning", "physicalSide": "both",
                                      "midplaneAngstrom": 0, "thicknessAngstrom": 30,
                                      "depthAngstrom": 0, "tiltDegrees": 0, "sidedness": None,
                                      "contactingRegions": [], "limitations": [],
                                      "policyId": "controlled-policy", "policyVersion": "1",
                                      "witnessId": "controlled-witness", "reason": "Controlled presentation fixture",
                                      "evidence": [], "prediction": None}
                value["study"]["adoptedPlacementProposalId"] = "placement-one"
                value["inspection"] = inspection(
                    "placement-one", "oriented-protein-with-proposed-membrane-bounds")
                value["inspection"].update(structureUrl="/api/structures/tiny?format=pdb",
                                           omittedMolecules=[])
                value["attempt"].update(stageKind="Construction", phase="providerPacking",
                                        constructed=None, progress=0.2,
                                        stopRequested=pending_stop,
                                        message="Provider packing remains in progress")
                value["actions"].extend([
                    {"kind": "buildAndMinimize", "subjectId": None, "enabled": False,
                     "reason": "The identified attempt is already running."},
                    {"kind": "stopAttempt", "subjectId": None,
                     "enabled": not pending_stop, "reason": None},
                ])
                host.structure_bytes = protein_bytes
                host.replace(value)
                browser, page = open_browser(playwright, width, height)
                try:
                    # Reconnect with an already admitted attempt and placement inspection.
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.locator(".workspace-context-tabs button.active"))\
                        .to_have_text("Preparation")
                    expect(page.locator(".scene-panel[aria-label='Current attempt inspection']"))\
                        .to_be_visible()
                    expect(page.locator(".stage-card")).to_have_count(0)
                    expect(page.locator(".placement-scene-legend")).to_have_count(0)
                    if pending_stop:
                        expect(page.locator(".attempt-running").last).to_contain_text(
                            "Stop requested · waiting for the worker")
                    expect(page.locator(".scene-view-disclosure")).to_contain_text(
                        "Protein shown", timeout=120000)
                    initial_preparation = settled_scene(page)
                    page.screenshot(path=str(ARTIFACTS / f"running-preparation-initial-{width}.png"),
                                    animations="disabled")
                    self.assertGreater(initial_preparation["representations"], 0)
                    self.assertGreater(initial_preparation["distance"],
                                       initial_preparation["targetDistance"] * 0.55)
                    self.assertLess(max(abs(a - b) for a, b in zip(
                        initial_preparation["target"], initial_preparation["center"])), 1e-4)
                    self.assertGreater(initial_preparation["zSpan"] /
                                       (2 * initial_preparation["radius"]), 0.45)

                    area(page, "Placement")
                    expect(page.locator(".scene-panel[aria-label='Placed protein and intended bilayer']"))\
                        .to_be_visible()
                    expect(page.locator(".placement-scene-legend")).to_contain_text(
                        "Measured target planes · no packed lipids", timeout=120000)
                    first = settled_scene(page)
                    self.assertGreater(first["distance"], first["targetDistance"] * 0.55)
                    page.screenshot(path=str(ARTIFACTS / f"running-placement-initial-{width}.png"),
                                    animations="disabled")
                    self.assertEqual(first["subject"], "3D structure for Protein placement proposal")
                    self.assertGreater(first["atoms"], 300)
                    self.assertLess(abs(first["target"][2]), 5,
                                    "The placement overview must stay near the intended bilayer core")
                    page.evaluate("""() => {
                      const camera = Reflect.get(document.querySelector('.viewer-mount'),
                        Symbol.for('molstar.viewer')).plugin.canvas3d.camera;
                      const snapshot = camera.getSnapshot();
                      camera.setState({position: [snapshot.position[0] + 4,
                        snapshot.position[1] + 3, snapshot.position[2] + 2]}, 0);
                    }""")
                    page.evaluate("""() => new Promise(resolve => requestAnimationFrame(() =>
                      requestAnimationFrame(resolve)))""")
                    saved_position = settled_scene(page)["position"]
                    page.screenshot(path=str(ARTIFACTS / f"running-placement-camera-{width}.png"),
                                    animations="disabled")

                    area(page, "Preparation")
                    expect(page.locator(".scene-panel[aria-label='Current attempt inspection']"))\
                        .to_be_visible()
                    expect(page.locator(".placement-scene-legend")).to_have_count(0)
                    expect(page.locator(".stage-card")).to_have_count(0)
                    preparation = settled_scene(page)
                    self.assertGreater(preparation["distance"],
                                       preparation["targetDistance"] * 0.55)
                    self.assertEqual(preparation["subject"], first["subject"])
                    self.assertEqual(preparation["atoms"], first["atoms"])
                    self.assertEqual(preparation["representations"],
                                     initial_preparation["representations"])
                    self.assertLess(preparation["representations"], first["representations"],
                                    "The intended-bilayer representation must leave the attempt canvas")
                    self.assertLess(max(abs(a - b) for a, b in zip(
                        preparation["target"], preparation["center"])), 1e-4)
                    self.assertGreater(preparation["zSpan"] / (2 * preparation["radius"]), 0.45,
                                       "The protein must occupy a useful fraction of the attempt view")
                    page.screenshot(path=str(ARTIFACTS / f"running-preparation-camera-{width}.png"),
                                    animations="disabled")

                    area(page, "Placement")
                    expect(page.locator(".placement-scene-legend")).to_contain_text(
                        "Measured target planes · no packed lipids", timeout=120000)
                    restored = settled_scene(page)
                    self.assertEqual(restored["representations"], first["representations"])
                    self.assertLess(max(abs(a - b) for a, b in zip(
                        saved_position, restored["position"])), 1e-5)
                    self.assertEqual(page.request.get(base + "/api/state").json()["attempt"]["attemptId"],
                                     "attempt-one")
                    self.assertEqual(page.request.get(base + "/api/state").json()["stages"], [])
                    self.assertTrue(all(kind == "selectInspectionSubject" for kind in host.commands),
                                    "Work-area navigation must not submit scientific work")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
