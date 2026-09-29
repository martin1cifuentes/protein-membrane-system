"""Focused browser presentation checks with a controlled host account.

Build `browser/` first. This exercises initial state, SSE loss/reconnection,
earlier-input meaning and subject evidence without launching a scientific provider.
The real Slice 5 acceptance route must separately exercise the published host.
"""

from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
import base64
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


def wait_for_tiny_fixture_pixels(page) -> dict:
    """Wait for actual central red/green model pixels, excluding Mol*'s corner axes."""
    mount = page.locator(".viewer-mount")
    expect(mount.locator("canvas")).to_have_count(1, timeout=60000)
    page.wait_for_function("""() => {
      const mount = document.querySelector('.viewer-mount');
      const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
      return !!viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
    }""", timeout=60000)
    deadline = time.monotonic() + 30
    observed = {"red": 0, "green": 0}
    while time.monotonic() < deadline:
        encoded = base64.b64encode(mount.screenshot()).decode("ascii")
        observed = page.evaluate("""async encoded => {
          const image = new Image();
          image.src = `data:image/png;base64,${encoded}`;
          await image.decode();
          const surface = document.createElement('canvas');
          surface.width = image.width;
          surface.height = image.height;
          const context = surface.getContext('2d', {willReadFrequently: true});
          context.drawImage(image, 0, 0);
          const rgba = context.getImageData(0, 0, image.width, image.height).data;
          let red = 0, green = 0;
          for (let y = Math.floor(image.height * .2); y < Math.floor(image.height * .8); y += 2)
            for (let x = Math.floor(image.width * .2); x < Math.floor(image.width * .8); x += 2) {
              const pixel = (y * image.width + x) * 4;
              const r = rgba[pixel], g = rgba[pixel + 1], b = rgba[pixel + 2];
              if (r > 110 && r > g + 45 && r > b + 45) red++;
              if (g > 65 && g > r + 25 && g > b + 15) green++;
            }
          return {red, green};
        }""", encoded)
        if observed["red"] >= 30 and observed["green"] >= 30:
            return observed
        page.wait_for_timeout(150)
    raise AssertionError(f"Tiny fixture model pixels did not render: {observed}")


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
                        "stopRequested": False,
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
    assessment = {"checkStanding": "checksIncomplete", "reason": "Required evidence unavailable",
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
        self.command_handler = None
        self.atom_requests: list[tuple[str, str, int]] = []
        self.structure_bytes = TINY_STRUCTURE
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
            structure_bytes = self.server.structure_bytes
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb")
            self.send_header("Content-Length", str(len(structure_bytes)))
            self.end_headers()
            self.wfile.write(structure_bytes)
            return
        if path.startswith("/api/inspection/components/"):
            parts = path.split("/")
            with self.server.lock:
                selected = deepcopy(self.server.account["inspection"])
            if (len(parts) == 6 and parts[4] == selected["subjectId"] and
                    parts[5] == "tiny" and selected["representationKind"] in
                    {"constructedSystem", "completedStage"}):
                self.json_response({"subjectId": parts[4],
                                    "studyRevisionId": selected["studyRevisionId"],
                                    "structureToken": "tiny", "atomCount": 2,
                                    "runs": [{"start": 0, "endExclusive": 2,
                                              "role": "lipid", "sourceChain": None,
                                              "copyId": None}]})
            else:
                self.send_error(404)
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
        if self.server.command_handler and self.server.command_handler(self, command):
            return
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
            if command["kind"] == "selectInspectionSubject" and \
                    command["data"].get("subjectId") == "diagnostic-one" and \
                    command["expectedRevision"] == current["revision"] and \
                    any(artifact.get("subjectId") == "diagnostic-one"
                        for attempt in [current.get("attempt"), *current.get("priorAttempts", [])]
                        if attempt is not None
                        for trial in attempt.get("trials", [])
                        for artifact in trial.get("diagnosticArtifacts", [])):
                updated = deepcopy(current)
                updated["revision"] += 1
                updated["inspection"] = inspection("diagnostic-one", "providerDiagnostic")
                original = next(attempt for attempt in
                                [current.get("attempt"), *current.get("priorAttempts", [])]
                                if attempt is not None and any(
                                    artifact.get("subjectId") == "diagnostic-one"
                                    for trial in attempt.get("trials", [])
                                    for artifact in trial.get("diagnosticArtifacts", [])))
                updated["inspection"].update(structureUrl="/api/structures/tiny?format=pdb",
                                             assessment=None, omittedMolecules=[],
                                             studyRevisionId=original["studyRevisionId"],
                                             studyRevisionNumber=1 if original["studyRevisionId"] == "revision-one" else 0)
                self.server.account = updated
                self.server.version += 1
                self.json_response(updated)
                return
            selected_stage = next((stage for stage in current["stages"]
                                   if stage["stageId"] == command["data"].get("subjectId")), None)
            if command["kind"] != "selectInspectionSubject" or \
                    selected_stage is None or command["expectedRevision"] != current["revision"]:
                self.json_response({"reason": "Only exact stage inspection is available."}, 422)
                return
            updated = deepcopy(current)
            updated["revision"] += 1
            updated["inspection"] = inspection(selected_stage["stageId"], "completedStage",
                                                selected_stage["assessment"])
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
    def test_current_warnings_follow_phase_and_resolution_without_policy_version_clutter(self):
        """Controlled account checks display scope; the Host owns issue applicability."""
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = completed(account())
            another_stage = deepcopy(value["stages"][0])
            another_stage["stageId"] = "stage-two"
            another_stage["attemptId"] = "attempt-two"
            another_stage["assessment"]["reason"] = "This result has its own check outcome."
            value["stages"].append(another_stage)
            value["study"].update(modelIndex=0, chainIds=["A", "B"])
            value["notices"] = [
                {"id": "protein-repair", "conditionKey": "protein-repair:source-one",
                 "severity": "error", "subjectId": "source-one",
                 "message": "Required protein repair is still unresolved.",
                 "affectedAreas": ["protein", "placement", "preparation"], "correctionArea": "protein"},
                {"id": "placement-check", "conditionKey": "placement-check:pose-one",
                 "severity": "error", "subjectId": "pose-one",
                 "message": "The current position check did not pass.",
                 "affectedAreas": ["placement", "preparation"], "correctionArea": "placement"},
                {"id": "older-result", "conditionKey": "stage-check:stage-one",
                 "severity": "warning", "subjectId": "stage-one",
                 "message": "An earlier result has a separate finding.",
                 "affectedAreas": ["results"], "correctionArea": "results"},
            ]
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 940})
                page.goto(base, wait_until="domcontentloaded")
                warnings = page.get_by_label("Current warnings")
                expect(warnings).to_be_visible()
                expect(warnings).to_contain_text("Required protein repair is still unresolved.")
                expect(warnings).to_contain_text("The current position check did not pass.")
                expect(warnings).not_to_contain_text("An earlier result has a separate finding.")
                expect(warnings.get_by_text("Blocking error")).to_have_count(2)
                expect(page.locator(".compact-model-context")).to_have_count(0)
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Membrane").click()
                expect(warnings).to_have_count(0)
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Placement").click()
                expect(warnings).to_contain_text("Required protein repair is still unresolved.")
                expect(warnings).to_contain_text("The current position check did not pass.")
                pending = deepcopy(value)
                pending["revision"] += 1
                pending["notices"] = [value["notices"][1], value["notices"][2]]
                host.replace(pending)
                expect(warnings).not_to_contain_text("Required protein repair is still unresolved.")
                expect(warnings).to_contain_text("The current position check did not pass.")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(warnings).to_be_visible()
                expect(warnings).to_contain_text("The current position check did not pass.")
                output = ROOT / "out" / "browser-acceptance" / "feedback-continuation"
                output.mkdir(parents=True, exist_ok=True)
                for width, height in ((1672, 940), (1024, 768), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    page.screenshot(path=str(output / f"collapsed-active-warnings-{width}.png"))
                page.get_by_role("button", name="Show inputs").click()
                resolved = deepcopy(pending)
                resolved["revision"] += 1
                resolved["notices"] = [value["notices"][2]]
                host.replace(resolved)
                expect(warnings).to_have_count(0)
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Results").click()
                expect(warnings).to_contain_text("An earlier result has a separate finding.")
                page.locator(".stage-card[data-stage-id='stage-two']").click()
                expect(warnings).to_have_count(0)
                page.locator(".stage-card[data-stage-id='stage-one']").click()
                expect(warnings).to_contain_text("An earlier result has a separate finding.")
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Preparation").click()
                expect(warnings).to_have_count(0)
                expect(page.get_by_text("Method version 1")).to_have_count(0)
            finally:
                browser.close()

    def test_protein_measurements_are_grouped_addressed_and_not_called_passes(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            value["attempt"] = None
            value["inspection"] = inspection("protein-one", "preparedProtein")
            values = [1.4999999999, 1.5000000001, 31.669000000000004] + [1.43] * 77
            geometry = {"standing": "Observed", "kinds": [{
                "kind": "covalentBond", "standing": "Observed", "eligibleCount": 80,
                "measuredCount": 80, "minimumDistanceAngstrom": min(values),
                "maximumDistanceAngstrom": max(values), "unavailableReason": None,
            }], "locatedDistances": [{
                "kind": "covalentBond", "distanceAngstrom": distance,
                "radiusSumAngstrom": None,
                "first": {"residue": {"model": 0, "chain": "A", "residue": index + 1,
                                       "insertionCode": "", "copyId": "A"}, "atomName": "C"},
                "second": {"residue": {"model": 0, "chain": "A", "residue": index + 2,
                                        "insertionCode": "", "copyId": "A"}, "atomName": "N"},
            } for index, distance in enumerate(values)],
                "limitations": ["Controlled coordinates are presentation evidence only."]}
            value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                "summary": "Controlled prepared protein", "atomCount": 162,
                                "changes": [], "findings": [], "prediction": None,
                                "geometry": geometry, "sourceGeometry": None}
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Protein").click()
                geometry_account = page.get_by_label("Protein geometry measurements")
                geometry_account.get_by_text("Measurement details").click()
                expect(geometry_account).to_contain_text("80 of 80 applicable distances measured")
                expect(geometry_account).not_to_contain_text("80 checks passed")
                addressed = geometry_account.locator(".geometry-address-list")
                expect(addressed.locator("summary")).to_contain_text("Measured bond lengths · 80 addressed measurements")
                addressed.locator("summary").click()
                expect(addressed.locator("li")).to_have_count(80)
                expect(addressed).to_contain_text("1.4999999999 Å")
                expect(addressed).to_contain_text("1.5000000001 Å")
                expect(addressed).to_contain_text("31.669 Å")
                expect(addressed).to_contain_text("Model 1 · chain A · residue 80 · atom C")
                expect(addressed).to_contain_text("Method: distance between the identified atom coordinates")
                expect(addressed).to_contain_text("Scope: prepared protein coordinates")
                self.assertEqual(host.commands, [], "Reading measurements must not change the scientific account")
            finally:
                browser.close()

    def test_running_current_operation_displays_each_provider_phase(self):
        """Controlled account: live operation is distinct from diagnostic history."""
        phases = ("providerPacking", "providerCleanup", "amberParameterization",
                  "providerConditioningRestrained", "providerConditioningUnrestrained")
        labels = ("Packing — arranging molecules around the protein",
                  "Checking and cleaning the packed system",
                  "Preparing molecular parameters", "Conditioning the starting system",
                  "Conditioning the starting system")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            value["inspection"] = None
            value["attempt"].update(stageKind="Construction", constructed=None,
                                    phase=phases[0], progress=0.1,
                                    message="Controlled active provider operation")
            value["attempt"]["trials"] = [{
                "trialId": "trial-one", "trialIndex": 0, "standing": "running",
                "lateralPaddingAngstrom": 15.0, "aqueousPaddingAngstrom": 17.5,
                "proposedLipidCounts": [], "achievedLipidCounts": [],
                "cleanupRemovedLipidCounts": [], "proposedCellAngstrom": [],
                "actualCellAngstrom": [], "conditions": None, "failureCode": None,
                "message": "The same trial is still running.",
                "diagnosticArtifacts": [{
                    "role": f"diagnostic{index}", "phase": phase,
                    "sha256": f"{index + 1:x}" * 64,
                    "fileName": f"diagnostic-{index}.log", "downloadUrl": None,
                    "structureUrl": None, "subjectId": None,
                } for index, phase in enumerate(phases)],
            }]
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Preparation", exact=True).click()
                progress = page.get_by_label("Observed preparation progress")
                current_line = progress.locator(".execution-progress-body.is-running strong")
                expect(page.get_by_label("System preparation view")).to_contain_text(
                    "No checked molecular system is available yet")
                expect(page.locator(".viewer-mount")).to_have_count(0)
                expect(page.get_by_role("button", name="Review attempt")).to_have_count(0)
                expect(page.get_by_role("button", name="Models")).to_have_count(0)
                expect(page.get_by_role("button", name="Notes")).to_have_count(0)
                for phase, label in zip(phases, labels):
                    expect(current_line).to_have_text(label)
                    expect(progress).to_contain_text("Method-internal diagnostics (5)")
                    value = deepcopy(value)
                    value["revision"] += 1
                    value["attempt"]["phase"] = (
                        phases[phases.index(phase) + 1]
                        if phase != phases[-1] else "finalMinimization")
                    if phase == phases[-1]:
                        value["attempt"]["stageKind"] = "Minimization"
                    host.replace(value)
                expect(current_line).to_have_text(
                    "Minimizing — adjusting atomic positions to reduce energy")
                expect(page.locator(".stage-strip")).to_contain_text(
                    "No minimized or equilibrated system stage has completed.")
                self.assertEqual(host.commands, [])
            finally:
                browser.close()

    def test_earlier_failed_attempt_keeps_its_diagnostics_after_new_attempt(self):
        """Controlled account: history is read-only and distinct from current controls."""
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            prior = deepcopy(value["attempt"])
            prior.update(attemptId="prior-one", studyRevisionId="revision-zero",
                         status="failed", stageKind="Construction", progress=None,
                         message="Prior provider conditioning failed.",
                         failureCode="providerFailure", constructed=None,
                         trials=[{
                             "trialId": "prior-trial", "trialIndex": 0, "standing": "failed",
                             "lateralPaddingAngstrom": 12.0, "aqueousPaddingAngstrom": 20.0,
                             "proposedLipidCounts": [], "achievedLipidCounts": [],
                             "cleanupRemovedLipidCounts": [], "proposedCellAngstrom": [],
                             "actualCellAngstrom": [], "conditions": None,
                             "failureCode": "providerFailure", "message": "Prior trial failed.",
                             "diagnosticArtifacts": [{
                                 "role": "providerPacked", "sha256": "a" * 64,
                                 "fileName": "prior-packed.pdb", "phase": "providerPacking",
                                 "downloadUrl": "/api/diagnostics/prior-packed",
                                 "structureUrl": "/api/structures/tiny?format=pdb",
                                 "subjectId": "diagnostic-one",
                             }],
                         }])
            value["priorAttempts"] = [prior]
            value["constructionRoutes"] = [{"policyId": "policy-one", "label": "Controlled native DMPC",
                                            "route": "nativeOpenMM", "saltConvention": "backgroundPairsPlusNeutralization",
                                            "available": True, "reason": None}]
            value["actions"].append({"kind": "buildAndMinimize", "subjectId": "policy-one",
                                     "enabled": False, "reason": "The current attempt is still running."})
            value["actions"].append({"kind": "stopAttempt", "subjectId": None,
                                     "enabled": True, "reason": None})
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Preparation", exact=True).click()
                expect(page.get_by_role("button", name="Build and minimize")).to_be_disabled()
                expect(page.get_by_role("button", name="Stop unfinished work")).to_be_enabled()
                history = page.get_by_label("Earlier preparation attempts")
                expect(history.locator("summary")).to_be_visible()
                history.locator("summary").click()
                expect(history.get_by_role("button", name="Earlier attempt 1 · failed")).to_be_visible()
                history.get_by_role("button", name="Earlier attempt 1 · failed").click()
                earlier = page.get_by_role("region", name="Earlier preparation attempt", exact=True)
                expect(earlier).to_contain_text(
                    "Read-only account from a prior attempt")
                expect(earlier).to_contain_text(
                    "This result uses earlier inputs")
                expect(page.locator(".attempt-account")).to_contain_text("Current system attempt")
                expect(page.get_by_role("button", name="Stop unfinished work")).to_be_enabled()
                expect(page.get_by_role("button", name="Show current run")).to_be_visible()
                self.assertEqual(host.commands, [])
                page.get_by_text("Construction trials and actual causes").click()
                page.get_by_text("Method-internal diagnostics (1)").click()
                expect(page.get_by_role("link", name="Download diagnostic")).to_have_attribute(
                    "href", "/api/diagnostics/prior-packed")
                page.get_by_role("button", name="Inspect diagnostic structure").click()
                expect(page.locator(".execution-scene-label")).to_contain_text("earlier attempt")
                expect(page.locator(".scene-caption")).to_contain_text(
                    "Provider method diagnostic; no completed stage")
                expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                    "Prior provider conditioning failed")
                expect(page.get_by_label("Subject applicability and technical checks")).to_contain_text(
                    "This result uses earlier inputs")
                self.assertEqual(host.account["inspection"]["studyRevisionId"], "revision-zero")
                page.reload(wait_until="domcontentloaded")
                page.get_by_role("button", name="Preparation", exact=True).click()
                expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                    "Prior provider conditioning failed")
                expect(page.locator(".scene-caption")).to_contain_text(
                    "Provider method diagnostic; no completed stage")
                screenshots = ROOT / "out" / "actor-reconciliation" / "current"
                screenshots.mkdir(parents=True, exist_ok=True)
                desktop_pixels = wait_for_tiny_fixture_pixels(page)
                page.screenshot(path=str(screenshots / "controlled-fixture-prior-diagnostic-1672.png"), full_page=True)
                page.set_viewport_size({"width": 820, "height": 760})
                self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth - innerWidth"), 0)
                compact_pixels = wait_for_tiny_fixture_pixels(page)
                page.screenshot(path=str(screenshots / "controlled-fixture-prior-diagnostic-820.png"), full_page=True)
                (screenshots / "controlled-fixture-prior-diagnostic-pixels.json").write_text(json.dumps({
                    "fixture": "controlled two-atom provider diagnostic",
                    "criterion": "At least 30 red and 30 green sampled model pixels in the central viewer area",
                    "1672": desktop_pixels, "820": compact_pixels,
                }, indent=2) + "\n")
            finally:
                browser.close()

    def test_failed_provider_trial_diagnostics_are_inspectable_without_completed_stage_claim(self):
        """Controlled account: diagnostic links retain failed trial identity."""
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = account()
            value["attempt"].update(status="failed", stageKind="Construction", progress=None,
                                    phase="providerConditioningUnrestrained",
                                    message="The provider trial failed before a checked handoff.",
                                    failureCode="providerFailure", constructed=None,
                                    trials=[{
                                        "trialId": "trial-one", "trialIndex": 0, "standing": "failed",
                                        "lateralPaddingAngstrom": 12.0,
                                        "aqueousPaddingAngstrom": 20.0,
                                        "proposedLipidCounts": [], "achievedLipidCounts": [],
                                        "cleanupRemovedLipidCounts": [],
                                        "proposedCellAngstrom": [], "actualCellAngstrom": [],
                                        "conditions": None, "failureCode": "providerFailure",
                                        "message": "Provider conditioning stopped without a checked construct.",
                                        "diagnosticArtifacts": [
                                            {"role": "providerPacked", "sha256": "a" * 64,
                                             "fileName": "packed.pdb", "phase": "providerPacking",
                                             "downloadUrl": "/api/diagnostics/packed-one",
                                             "structureUrl": "/api/structures/tiny?format=pdb",
                                             "subjectId": "diagnostic-one"},
                                            {"role": "providerUnrestrainedLog", "sha256": "b" * 64,
                                             "fileName": "unrestrained.log",
                                             "phase": "providerConditioningUnrestrained",
                                             "downloadUrl": "/api/diagnostics/log-one",
                                             "structureUrl": None, "subjectId": None},
                                        ],
                                    }])
            value["inspection"] = None
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Preparation", exact=True).click()
                expect(page.get_by_label("System preparation view")).to_contain_text(
                    "No checked constructed system")
                expect(page.locator(".viewer-mount")).to_have_count(0)
                expect(page.get_by_role("button", name="Review attempt")).to_have_count(0)
                page.get_by_text("Construction trials and actual causes").click()
                trial = page.get_by_text("Trial 1 · failed")
                expect(trial).to_be_visible()
                page.get_by_text("Method-internal diagnostics (2)").click()
                expect(page.get_by_text("Provider Packed")).to_be_visible()
                expect(page.get_by_text("Provider Unrestrained Log")).to_be_visible()
                downloads = page.get_by_role("link", name="Download diagnostic")
                self.assertEqual(downloads.count(), 2)
                self.assertEqual([downloads.nth(index).get_attribute("href") for index in range(2)],
                                 ["/api/diagnostics/packed-one", "/api/diagnostics/log-one"])
                self.assertEqual(page.get_by_role("button", name="Inspect diagnostic structure").count(), 1)
                expect(page.locator(".stage-strip")).to_contain_text("No minimized or equilibrated system stage has completed")
                page.get_by_role("button", name="Inspect diagnostic structure").click()
                expect(page.locator(".execution-scene-label")).to_contain_text(
                    "Provider method diagnostic")
                expect(page.locator(".scene-caption")).to_contain_text(
                    "Provider method diagnostic; no completed stage")
                self.assertEqual(host.account["inspection"]["subjectId"], "diagnostic-one")
                self.assertIsNone(host.account["inspection"]["assessment"])
                page.reload(wait_until="domcontentloaded")
                expect(page.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Preparation")
                expect(page.locator(".scene-caption")).to_contain_text(
                    "Provider method diagnostic; no completed stage")
                expect(page.locator(".stage-strip")).to_contain_text("No minimized or equilibrated system stage has completed")
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-provider-diagnostic-1672.png"), full_page=True)
                page.set_viewport_size({"width": 820, "height": 760})
                self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth - innerWidth"), 0)
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-provider-diagnostic-820.png"), full_page=True)
            finally:
                browser.close()

    def test_checked_candidate_remains_an_intermediate_under_one_authorization(self):
        """Presentation fixture; owner composition tests establish command authority."""
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            ready = account()
            ready["attempt"].update(status="running", stageKind="Construction",
                                    phase="handoffChecks", progress=0.8,
                                    message="Presentation fixture: checking this construction before final minimization.")
            ready["attempt"]["constructed"]["maximumProteinCoordinateDeviationAngstrom"] = 2.125
            ready["constructionRoutes"] = [{"policyId": "policy-one", "label": "Controlled native DMPC",
                                            "route": "nativeOpenMM", "saltConvention": "backgroundPairsPlusNeutralization",
                                            "available": True, "reason": None}]
            ready["actions"].extend([
                {"kind": "buildAndMinimize", "subjectId": "policy-one",
                 "enabled": False, "reason": "The identified attempt is already running."},
                {"kind": "stopAttempt", "subjectId": None, "enabled": True, "reason": None},
            ])
            host.replace(ready)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 950})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Preparation", exact=True).click()
                candidate = page.get_by_label("Checked constructed candidate")
                expect(candidate).to_contain_text("1,000 atoms")
                expect(candidate).to_contain_text("upper DMPC 2")
                expect(candidate).to_contain_text("50.0 Å × 50.0 Å × 80.0 Å")
                expect(candidate).to_contain_text("The checked construction is an intermediate")
                expect(page.locator(".rail").get_by_role(
                    "button", name="Authorize minimization of this candidate")).to_have_count(0)
                expect(page.get_by_role("button", name="Build and minimize")).to_be_disabled()
                expect(page.locator(".attempt-running")).to_contain_text(
                    "Checking the constructed system")
                expect(page.locator(".stage-strip")).to_contain_text(
                    "No minimized or equilibrated system stage has completed.")
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-candidate-intermediate-1672.png"), full_page=True)
                page.set_viewport_size({"width": 820, "height": 720})
                self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth - innerWidth"), 0)
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-candidate-intermediate-820.png"), full_page=True)
                page.get_by_role("button", name="Collapse inputs").click()
                actual = page.get_by_label("Validated actual constructed system")
                expect(actual).to_contain_text("Upper physical leafletDMPC 2")
                expect(actual).to_contain_text("Lower physical leafletDMPC 2")
                expect(page.locator(".stage-strip")).to_contain_text(
                    "No minimized or equilibrated system stage has completed.")
                page.get_by_role("button", name="Show inputs").click()
                candidate.locator("details", has_text="Construction conditions").click()
                expect(candidate).to_contain_text("2.125 Å")
                failed = deepcopy(ready)
                failed["revision"] += 1
                failed["attempt"].update(status="failed", stageKind="Minimization",
                                         phase=None, progress=None,
                                         message="Required final minimization failed before a completed stage.")
                host.replace(failed)
                expect(candidate).to_contain_text("The candidate is retained for inspection")
                expect(candidate).not_to_contain_text("authorization continues")
                expect(page.locator(".stage-strip")).to_contain_text(
                    "No minimized or equilibrated system stage has completed.")
                finished = completed(ready)
                host.replace(finished)
                expect(candidate).to_contain_text("Final minimization established a completed stage")
            finally:
                browser.close()

    def test_requested_stop_is_not_presented_as_observed_termination(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            pending = account()
            pending["attempt"]["stopRequested"] = True
            pending["notices"] = [{"id": "fixture-stop", "severity": "information", "subjectId": "attempt-one",
                                   "message": "Presentation fixture: stop requested; worker outcome not yet observed."}]
            host.replace(pending)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 950})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".attempt-account")).to_contain_text(
                    "Stop requested — awaiting observed stop")
                expect(page.locator(".attempt-account")).to_contain_text("running")
                expect(page.locator(".attempt-account")).not_to_contain_text("Work stopped")
                expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                    "The worker’s observed outcome is pending")
                expect(page.locator(".stage-strip")).to_contain_text("None yet")
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-stop-requested-1672.png"), full_page=True)
                stopped = deepcopy(pending)
                stopped["revision"] += 1
                stopped["attempt"].update(status="stopped", stopRequested=False,
                                          message="Worker cancellation observed; no stage completed.")
                host.replace(stopped)
                expect(page.locator(".attempt-account")).to_contain_text("Work stopped")
                expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                    "Worker cancellation observed")
                expect(page.locator(".attempt-account")).not_to_contain_text("Stop requested")
                expect(page.locator(".stage-strip")).to_contain_text("None yet")
                unknown = deepcopy(stopped)
                unknown["revision"] += 1
                unknown["attempt"].update(status="unobserved", message="Worker outcome unavailable")
                host.replace(unknown)
                expect(page.locator(".attempt-account")).to_contain_text("Worker outcome unavailable")
                expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                    "Worker outcome unavailable")
                expect(page.locator(".attempt-account")).not_to_contain_text("Work stopped")
                page.set_viewport_size({"width": 820, "height": 720})
                page.screenshot(path=str(ROOT / "out" / "actor-reconciliation" / "current" /
                                         "controlled-fixture-stop-unobserved-820.png"), full_page=True)
            finally:
                browser.close()

    def assert_stage_hierarchy(self, page, viewport_height: int):
        labels = page.locator(".evidence-content > .account-card").evaluate_all(
            "elements => elements.map(element => element.getAttribute('aria-label'))")
        self.assertEqual(labels[:2], ["Completed stage information",
                                      "Minimized stage review and technical checks"])
        self.assertNotIn("Selected subject and study revision", labels)
        expect(page.get_by_label("Inspection selection")).to_have_count(0)
        stage_information = page.get_by_label("Completed stage information")
        stage_review = page.get_by_label("Minimized stage review and technical checks")
        expect(stage_information).to_contain_text("This result uses earlier inputs")
        expect(stage_information).not_to_contain_text("revision-one")
        expect(stage_information).not_to_contain_text("revision-two")
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
            value["attempt"]["constructed"].update(atomCount=2, waterCount=0,
                                                      sodiumCount=0, chlorideCount=0)
            value["stages"][0]["constructed"] = value["attempt"]["constructed"]
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
                    "This result uses earlier inputs")
                expect(page.get_by_label("Minimized stage review and technical checks"))\
                    .to_contain_text("Required evidence unavailable")
                expect(page.get_by_label("Stage-specific findings and evidence")).to_contain_text(
                    "Observation for stage-one")
                expect(page.locator(".stage-card.selected")).to_be_visible()
                self.assertEqual(host.commands, [], "A missing representation must not alter the stage")

                page.unroute("**/api/structures/*")
                page.reload(wait_until="domcontentloaded")
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                expect(page.locator(".scene-error")).to_have_count(0)
                self.assertTrue(page.locator(".viewer-mount").evaluate("""mount => {
                  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                  return !!viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
                }"""), "The recovered view must contain the selected structure")
                expect(page.get_by_label("Completed stage information")).to_contain_text(
                    "Verified protein + bilayer + water and ions")
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
            value["attempt"]["constructed"].update(atomCount=2, waterCount=0,
                                                      sodiumCount=0, chlorideCount=0)
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
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                linked = page.get_by_label("Selected subject findings and evidence")
                expect(linked).to_contain_text("Located finding 0")
                expect(linked).to_contain_text("identified observation")
                expect(linked).to_contain_text("Controlled presentation fixture")
                measured = page.locator(".execution-metrics")
                expect(measured).to_contain_text("Nearest Distance")
                expect(measured).to_contain_text("3.4 Å")
                page.get_by_role("button", name="Located observed region").click()
                expect(page.get_by_label("Inspection selection")).to_contain_text(
                    "Located observed region")
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
                expect(page.get_by_role("button", name="Models")).to_have_count(0)
                expect(page.get_by_role("button", name="Notes")).to_have_count(0)
                expect(page.get_by_label("Details")).to_be_visible()
                # Choosing a located finding already opens its local neighborhood.
                expect(page.get_by_role("button", name="Inspect locally")).to_be_disabled()
                expect(page.locator(".scene-local-state")).to_contain_text("5 Å")
                self.assertEqual(host.account["study"]["id"], "revision-one")
                self.assertEqual(host.account["attempt"]["status"], "running")
                self.assertEqual(host.account["stages"], [])
                self.assertEqual(host.commands, ["setInspectionFocus"],
                                 "Viewing details and focus must not issue other scientific commands")
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
                expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
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
                    if row == 0:
                        page.get_by_text("Exact coordinate correspondence").click()
                    expect(page.get_by_label("Inspection selection")).to_contain_text(
                        f"lipid:[{row},DMP]")
                    expect(page.get_by_label("Inspection selection")).to_contain_text(
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
                expect(page.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Preparation")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(page.locator(".execution-decision")).to_have_count(1)
                expect(page.get_by_label("Current preparation attempt")).to_contain_text(
                    "Validated protein + bilayer + water and ions")
                expect(page.get_by_label("Current attempt inspection")).to_contain_text(
                    "Protein–membrane system · current inputs")
                expect(page.get_by_label("Subject applicability and technical checks")).to_have_count(0)
                self.assertEqual(host.account["attempt"]["attemptId"], "attempt-one")
                expect(page.locator(".stage-strip")).to_contain_text("None yet")
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
                expect(interrupted.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Preparation")
                interrupted.get_by_role("button", name="Collapse inputs").click()
                expect(interrupted.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(interrupted.locator(".execution-decision")).to_have_count(1)
                expect(interrupted.get_by_role("alert").filter(
                    has_text="The live connection is interrupted")).to_be_visible()
                expect(interrupted.get_by_label("Subject applicability and technical checks"))\
                    .to_contain_text("This result uses earlier inputs")
                expect(interrupted.get_by_label("Observed preparation progress")).to_contain_text(
                    "62% reported for the current operation")
                self.assertEqual(host.account["attempt"]["message"], "Later observed progress")
                expect(interrupted.locator(".stage-strip")).to_contain_text("None yet")
                interrupted.close()

                reconnected = browser.new_page(viewport={"width": 1672, "height": 941})
                reconnected.goto(base, wait_until="domcontentloaded")
                expect(reconnected.get_by_label("Observed preparation progress")).to_contain_text(
                    "62% reported for the current operation")
                host.replace(completed(revised))
                expect(reconnected.locator(".stage-card")).to_have_count(1)
                expect(reconnected.locator(".stage-card")).to_be_visible()
                reconnected.locator(".stage-card").click()
                reconnected.get_by_role("button", name="Collapse inputs").click()
                expect(reconnected.get_by_label("Completed stage information")).to_contain_text(
                    "This result uses earlier inputs")
                self.assert_stage_hierarchy(reconnected, 941)
                expect(reconnected.get_by_label("Stage-specific findings and evidence")).to_contain_text(
                    "Finding for stage-one")
                expect(reconnected.get_by_label("Selected subject findings and evidence")).to_have_count(0)
                reconnected.close()

                reopened = browser.new_page(viewport={"width": 820, "height": 760})
                reopened.goto(base, wait_until="domcontentloaded")
                expect(reopened.locator(".workspace-context-tabs button.active"))\
                    .to_have_text("Results")
                reopened.get_by_role("button", name="Collapse inputs").click()
                expect(reopened.locator(".workspace.workflow-closed")).to_have_count(1)
                expect(reopened.locator(".execution-decision")).to_have_count(1)
                expect(reopened.get_by_label("Completed stage information")).to_contain_text(
                    "This result uses earlier inputs")
                self.assert_stage_hierarchy(reopened, 760)
                expect(reopened.locator(".stage-strip")).to_contain_text("Checks incomplete")
                self.assertEqual(host.commands, ["selectInspectionSubject"],
                                 "Reload and SSE must not start or stop scientific work")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
