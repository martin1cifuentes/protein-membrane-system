"""Focused browser presentation checks with a controlled host account.

Build `browser/` first. This exercises initial state, SSE loss/reconnection,
revision labels and subject evidence without launching a scientific provider.
The real Slice 5 acceptance route must separately exercise the published host.
"""

from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
import os
from pathlib import Path
import threading
import time
import unittest
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
DIST = ROOT / "browser" / "dist"
CHROMIUM = os.environ.get("PIM_BROWSER_CHROMIUM", "/home/marti/.cache/ms-playwright/chromium-1228/chrome-linux64/chrome")
TINY_STRUCTURE = ("HETATM    1  C1  DMP A   1       0.000   0.000   0.000  1.00  0.00           C\n"
                  "HETATM    2  O1  DMP A   1       1.500   0.000   0.000  1.00  0.00           O\n"
                  "CONECT    1    2\nEND\n").encode()


def evidence(subject: str) -> dict:
    return {"id": f"evidence-{subject}", "subjectId": subject,
            "source": "controlled host account", "method": "identified observation",
            "observation": f"Observation for {subject}", "applicability": f"Only {subject}",
            "uncertainty": "Controlled presentation fixture", "bearing": "context"}


def finding(subject: str) -> dict:
    return {"id": f"finding-{subject}", "subjectId": subject,
            "evidenceId": f"evidence-{subject}", "meaning": f"Finding for {subject}",
            "consequence": "Review this subject", "disposition": "context", "material": False}


def inspection(subject: str, kind: str, assessment: dict | None = None) -> dict:
    return {"subjectId": subject, "studyRevisionId": "revision-one",
            "studyRevisionNumber": 1, "structureUrl": None,
            "representationKind": kind, "omittedMolecules": ["water rendering unavailable"],
            "evidence": [evidence(subject)], "findings": [finding(subject)],
            "assessment": assessment, "focusId": None, "focus": None,
            "annotations": [{"id": "unlocated", "subjectPartId": subject,
                             "label": "Unlocated context", "meaning": "No verified atom focus",
                             "evidenceId": f"evidence-{subject}", "geometryFocus": None}],
            "metrics": [{"name": "observed count", "value": "12", "unit": "pairs",
                         "subjectPartId": subject, "evidenceId": f"evidence-{subject}"}]}


def account() -> dict:
    constructed = {"subjectId": "constructed-one", "attemptId": "attempt-one",
                   "atomCount": 1000,
                   "achievedComposition": [{"physicalSide": "upper", "speciesId": "DMPC", "count": 2,
                                            "intendedFraction": 1.0},
                                           {"physicalSide": "lower", "speciesId": "DMPC", "count": 2,
                                            "intendedFraction": 1.0}],
                   "cellAngstrom": [50, 50, 80], "waterCount": 100,
                   "sodiumCount": 2, "chlorideCount": 2,
                   "conditionsTreatment": "Observed construction conditions", "localState": None}
    return {"revision": 10,
            "study": {"id": "revision-one", "number": 1, "summary": "Controlled study",
                      "selectedSourceId": None, "modelIndex": None,
                      "adoptedPlacementProposalId": None, "biologicalAssemblyId": None,
                      "chainIds": [], "partners": [], "alternateLocations": [],
                      "conditions": {"nominalPh": 7.0, "targetNaClMolar": 0.15,
                                     "optionalTemperatureKelvin": 300}},
            "sourceCandidates": [], "sourceModels": [], "sourcePrediction": None,
            "availableLipids": [], "protein": None, "membrane": None, "placement": None,
            "attempt": {"attemptId": "attempt-one", "status": "running",
                        "stageKind": "Minimization", "progress": 0.35,
                        "message": "Observed minimization progress", "studyRevisionId": "revision-one",
                        "policyId": "policy-one", "policyVersion": "1", "currentStageId": None,
                        "derivation": None, "constructed": constructed},
            "stages": [], "inspection": inspection("constructed-one", "constructedSystem"),
            "actions": [{"kind": "selectInspectionSubject", "subjectId": None,
                         "enabled": True, "reason": None},
                        {"kind": "setInspectionFocus", "subjectId": None,
                         "enabled": True, "reason": None}],
            "notices": []}


def completed(account_value: dict) -> dict:
    result = deepcopy(account_value)
    result["revision"] += 1
    result["attempt"].update(status="completed", stageKind="Minimization", progress=1.0,
                             message="Completed minimized stage established", currentStageId="stage-one")
    assessment = {"qualification": "indeterminate", "reason": "Required evidence unavailable",
                  "evidence": [evidence("stage-one")], "findings": [finding("stage-one")],
                  "limitations": ["One observation remains unavailable"],
                  "currentlyApplicable": True}
    result["stages"] = [{"stageId": "stage-one", "attemptId": "attempt-one",
                         "studyRevisionId": "revision-one", "kind": "Minimization",
                         "status": "completed", "assessment": assessment,
                         "summary": "Completed minimized stage", "observation": None,
                         "constructed": result["attempt"]["constructed"]}]
    return result


class AccountServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self):
        self.account = account()
        self.version = 0
        self.commands: list[str] = []
        self.atom_requests: list[tuple[str, str, int]] = []
        self.lock = threading.Lock()
        super().__init__(("127.0.0.1", 0), AccountHandler)

    def replace(self, value: dict) -> None:
        with self.lock:
            self.account = value
            self.version += 1


class AccountHandler(BaseHTTPRequestHandler):
    server: AccountServer

    def log_message(self, *_args):
        pass

    def json_response(self, value: dict, status: int = 200):
        body = json.dumps(value).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlsplit(self.path).path
        if path == "/api/state":
            with self.server.lock:
                value = deepcopy(self.server.account)
            self.json_response(value)
            return
        if path == "/api/events":
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            last = -1
            try:
                for _ in range(240):
                    with self.server.lock:
                        current = self.server.version
                    if current != last:
                        self.wfile.write(b"data: refresh\n\n")
                        self.wfile.flush()
                        last = current
                    time.sleep(0.1)
            except (BrokenPipeError, ConnectionResetError):
                pass
            return
        if path == "/api/structures/tiny":
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb")
            self.send_header("Content-Length", str(len(TINY_STRUCTURE)))
            self.end_headers()
            self.wfile.write(TINY_STRUCTURE)
            return
        if path.startswith("/api/inspection/atoms/"):
            parts = path.split("/")
            if len(parts) == 7 and parts[4] == "constructed-one" and parts[5] == "tiny" and parts[6] in {"0", "1"}:
                row = int(parts[6])
                with self.server.lock:
                    self.server.atom_requests.append((parts[4], parts[5], row))
                self.json_response({"subjectId": parts[4], "studyRevisionId": "revision-one",
                                    "structureToken": parts[5], "atomSiteIndex": row,
                                    "atom": {"resultAtomIndex": row, "resultAtomId": f"lipid:[{row},DMP]",
                                             "sourceAtomId": None, "role": "generated",
                                             "moleculeRole": "lipid", "atomRole": "headGroup",
                                             "element": "C" if row == 0 else "O",
                                             "sourceResidue": None, "approvedChangeId": None,
                                             "physicalSide": "upper", "generatedSpeciesId": "DMPC",
                                             "generatedComponentRole": "lipid"}})
            else:
                self.send_error(404)
            return
        target = DIST / (path.removeprefix("/") or "index.html")
        if not target.is_file() or not target.resolve().is_relative_to(DIST.resolve()):
            self.send_error(404)
            return
        body = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", mimetypes.guess_type(target.name)[0] or "application/octet-stream")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        if urlsplit(self.path).path != "/api/commands":
            self.send_error(404)
            return
        command = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        with self.server.lock:
            self.server.commands.append(command["kind"])
            current = self.server.account
            if command["kind"] == "setInspectionFocus" and \
                    command["data"].get("annotationId") == "located" and \
                    command["expectedRevision"] == current["revision"] and \
                    current["inspection"]["annotations"][0]["id"] == "located":
                updated = deepcopy(current)
                annotation = updated["inspection"]["annotations"][0]
                updated["revision"] += 1
                updated["inspection"]["focusId"] = annotation["subjectPartId"]
                updated["inspection"]["focus"] = annotation["geometryFocus"]
                self.server.account = updated
                self.server.version += 1
                self.json_response(updated)
                return
            if command["kind"] != "selectInspectionSubject" or \
                    command["data"].get("subjectId") != "stage-one" or \
                    not current["stages"] or command["expectedRevision"] != current["revision"]:
                self.json_response({"reason": "Only exact stage inspection is available."}, 422)
                return
            updated = deepcopy(current)
            updated["revision"] += 1
            updated["inspection"] = inspection("stage-one", "completedStage",
                                                updated["stages"][0]["assessment"])
            self.server.account = updated
            self.server.version += 1
        self.json_response(updated)


@contextmanager
def controlled_account_server():
    server = AccountServer()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server, f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


class BrowserPresentationTests(unittest.TestCase):
    def assert_stage_hierarchy(self, page, viewport_height: int):
        labels = page.locator(".evidence-content > .account-card").evaluate_all(
            "elements => elements.map(element => element.getAttribute('aria-label'))")
        self.assertEqual(labels[:2], ["Completed stage information",
                                      "Minimized stage review and distinct scientific assessment"])
        self.assertNotIn("Selected subject and study revision", labels)
        expect(page.get_by_label("Selected atom correspondence")).to_have_count(0)
        stage_information = page.get_by_label("Completed stage information")
        stage_review = page.get_by_label("Minimized stage review and distinct scientific assessment")
        expect(stage_information).to_contain_text("Historical stage from study revision revision-one")
        expect(stage_information).to_contain_text("current revision 2")
        self.assertLess(stage_information.bounding_box()["y"], stage_review.bounding_box()["y"])
        self.assertLess(stage_review.bounding_box()["y"], viewport_height,
                        "The minimized review must begin in the initial viewport")

    def test_missing_stage_structure_keeps_stage_account_and_recovers_on_reload(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = completed(account())
            value["study"].update(id="revision-two", number=2)
            value["inspection"] = inspection("stage-one", "completedStage",
                                             value["stages"][0]["assessment"])
            value["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
            value["inspection"]["omittedMolecules"] = []
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1280, "height": 800})
                page.route("**/api/structures/*", lambda route: route.fulfill(
                    status=404, content_type="text/plain", body="Structure unavailable"))
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".scene-error")).to_contain_text("Structure unavailable", timeout=60000)
                expect(page.get_by_label("Completed stage information")).to_contain_text(
                    "Historical stage from study revision revision-one")
                expect(page.get_by_label("Minimized stage review and distinct scientific assessment"))\
                    .to_contain_text("Required evidence unavailable")
                expect(page.get_by_label("Selected subject findings and evidence")).to_contain_text(
                    "Observation for stage-one")
                expect(page.locator(".stage-strip")).to_contain_text("stage-one")
                self.assertEqual(host.commands, [], "A missing representation must not alter the stage")

                page.unroute("**/api/structures/*")
                page.reload(wait_until="domcontentloaded")
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                expect(page.locator(".scene-error")).to_have_count(0)
                self.assertTrue(page.locator(".viewer-mount").evaluate("""mount => {
                  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                  return !!viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
                }"""), "The recovered view must contain the selected structure")
                expect(page.get_by_label("Completed stage information")).to_contain_text("stage-one")
                self.assertEqual(host.commands, [], "Reload must not change the scientific stage")
            finally:
                browser.close()

    def test_located_evidence_focus_and_long_findings_remain_reachable(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            subject = value["inspection"]["subjectId"]
            value["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
            value["inspection"]["omittedMolecules"] = []
            value["inspection"]["annotations"] = [{
                "id": "located", "subjectPartId": "A:1", "label": "Located observed region",
                "meaning": "A measured subject part", "evidenceId": f"evidence-{subject}",
                "geometryFocus": {"authAsymId": "A", "authSeqId": 1,
                                  "insertionCode": None, "authAtomId": None},
            }]
            value["inspection"]["metrics"] = [{
                "name": "nearest distance", "value": "3.4", "unit": "Å",
                "subjectPartId": "A:1", "evidenceId": f"evidence-{subject}",
            }]
            value["inspection"]["findings"] = [
                {**finding(subject), "id": f"finding-{number}",
                 "meaning": f"Located finding {number} for {subject}"}
                for number in range(80)
            ]
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 820, "height": 760})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                linked = page.get_by_label("Selected subject findings and evidence")
                expect(linked).to_contain_text("Located finding 0")
                expect(linked).to_contain_text("identified observation")
                expect(linked).to_contain_text("Controlled presentation fixture")
                measured = page.locator(".execution-metrics")
                expect(measured).to_contain_text("nearest distance")
                expect(measured).to_contain_text("3.4 Å")
                page.get_by_role("button", name="Located observed region").click()
                expect(page.get_by_label("Selected subject and study revision")).to_contain_text(
                    "Current study revision")
                expect(page.get_by_role("button", name="Located observed region")).to_have_class(
                    "annotation-item selected")
                self.assertEqual(host.account["inspection"]["focusId"], "A:1")
                self.assertEqual(host.account["study"]["id"], "revision-one")
                self.assertEqual(host.account["attempt"]["status"], "running")
                self.assertEqual(host.account["stages"], [])
                self.assertEqual(host.commands, ["setInspectionFocus"])
                panel = page.locator(".evidence-content")
                panel.evaluate("element => { element.scrollTop = element.scrollHeight; }")
                expect(linked).to_contain_text("Located finding 79")
                last = linked.locator(".finding-item").last
                last.scroll_into_view_if_needed()
                expect(last).to_be_visible()
                panel_box = panel.bounding_box()
                last_box = last.bounding_box()
                self.assertIsNotNone(panel_box)
                self.assertIsNotNone(last_box)
                self.assertGreater(min(panel_box["y"] + panel_box["height"],
                                       last_box["y"] + last_box["height"]) -
                                   max(panel_box["y"], last_box["y"]), 0,
                                   "The last finding must intersect the evidence viewport")
                self.assertGreater(panel.evaluate("element => element.scrollTop"), 0)
                expect(page.locator(".execution-decision")).to_be_visible()
                expect(page.locator(".stage-strip")).to_be_visible()
                page.get_by_role("button", name="Models").click()
                page.get_by_role("button", name="Notes").click()
                page.locator(".placement-view-tools").get_by_role("button", name="Focus").click()
                self.assertEqual(host.account["study"]["id"], "revision-one")
                self.assertEqual(host.account["attempt"]["status"], "running")
                self.assertEqual(host.account["stages"], [])
                self.assertEqual(host.commands, ["setInspectionFocus"],
                                 "View and evidence tabs must not issue scientific commands")
            finally:
                browser.close()

    def test_molstar_click_resolves_exact_coordinate_row_through_host(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            account_value = account()
            account_value["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
            account_value["inspection"]["omittedMolecules"] = []
            account_value["attempt"]["constructed"].update(
                atomCount=2, waterCount=0, sodiumCount=0, chlorideCount=0)
            host.replace(account_value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1280, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                expect(page.locator(".scene-error")).to_have_count(0)
                expect(page.get_by_label("Selected atom correspondence")).to_have_count(0)
                for row in (0, 1):
                    selected = page.evaluate("""row => {
                      const mount = document.querySelector('.viewer-mount');
                      const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                      const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
                      if (!structure) return false;
                      for (const unit of structure.units) {
                        for (let offset = 0; offset < unit.elements.length; offset++) {
                          const element = unit.elements[offset];
                          if (unit.model.atomicHierarchy.atomSourceIndex.value(element) !== row) continue;
                          const loci = { kind: 'element-loci', structure,
                            elements: [{ unit, indices: Int32Array.of(offset) }] };
                          viewer.plugin.behaviors.interaction.click.next({
                            current: { loci, repr: null }, buttons: 1, button: 1,
                            modifiers: { alt: false, control: false, meta: false, shift: false }
                          });
                          return true;
                        }
                      }
                      return false;
                    }""", row)
                    self.assertTrue(selected, f"Mol* did not contain coordinate row {row}")
                    expect(page.get_by_label("Selected atom correspondence")).to_contain_text(
                        f"lipid:[{row},DMP]")
                    expect(page.get_by_label("Selected atom correspondence")).to_contain_text(
                        "DMPC")
                self.assertEqual(host.atom_requests, [("constructed-one", "tiny", 0),
                                                      ("constructed-one", "tiny", 1)])
                self.assertEqual(host.commands, [], "Mol* selection must not make a scientific command")
            finally:
                browser.close()

    def test_same_host_account_survives_view_loss_and_sse_interruption(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(page.locator(".execution-decision")).to_have_count(1)
                expect(page.get_by_label("Current preparation attempt")).to_contain_text("attempt-one")
                expect(page.get_by_label("Selected subject and study revision")).to_contain_text(
                    "Current study revision")
                expect(page.locator(".stage-strip")).to_contain_text("No completed stage")
                expect(page.get_by_label("Selected subject findings and evidence")).to_contain_text(
                    "Observation for constructed-one")
                expect(page.get_by_role("button", name="Unlocated context")).to_be_disabled()
                page.close()

                revised = deepcopy(host.account)
                revised["revision"] += 1
                revised["study"].update(id="revision-two", number=2)
                revised["attempt"].update(progress=0.62, message="Later observed progress")
                host.replace(revised)
                interrupted = browser.new_page(viewport={"width": 820, "height": 760})
                interrupted.route("**/api/events", lambda route: route.abort("failed"))
                interrupted.goto(base, wait_until="domcontentloaded")
                expect(interrupted.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(interrupted.locator(".execution-decision")).to_have_count(1)
                expect(interrupted.get_by_role("alert").filter(
                    has_text="The live connection is interrupted")).to_be_visible()
                expect(interrupted.get_by_label("Selected subject and study revision")).to_contain_text(
                    "Historical subject; current study is revision 2")
                expect(interrupted.get_by_label("Observed preparation progress")).to_contain_text(
                    "Later observed progress")
                expect(interrupted.locator(".stage-strip")).to_contain_text("No completed stage")
                interrupted.close()

                reconnected = browser.new_page(viewport={"width": 1672, "height": 941})
                reconnected.goto(base, wait_until="domcontentloaded")
                expect(reconnected.get_by_label("Observed preparation progress")).to_contain_text(
                    "Later observed progress")
                host.replace(completed(revised))
                expect(reconnected.locator(".stage-card")).to_have_count(1)
                expect(reconnected.locator(".stage-card")).to_contain_text("stage-one")
                reconnected.locator(".stage-card").click()
                expect(reconnected.get_by_label("Completed stage information")).to_contain_text(
                    "Historical stage from study revision revision-one")
                self.assert_stage_hierarchy(reconnected, 941)
                expect(reconnected.get_by_label("Selected subject findings and evidence")).to_contain_text(
                    "Observation for stage-one")
                reconnected.close()

                reopened = browser.new_page(viewport={"width": 820, "height": 760})
                reopened.goto(base, wait_until="domcontentloaded")
                expect(reopened.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(reopened.locator(".execution-decision")).to_have_count(1)
                expect(reopened.get_by_label("Completed stage information")).to_contain_text(
                    "Historical stage from study revision revision-one")
                self.assert_stage_hierarchy(reopened, 760)
                expect(reopened.locator(".stage-strip")).to_contain_text("indeterminate")
                self.assertEqual(host.commands, ["selectInspectionSubject"],
                                 "Reload and SSE must not start or stop scientific work")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
