"""Controlled browser proof for exact-stage export delivery and retry.

Build browser/ first. The published host/worker route separately proves the
real bundle's molecular content; this fixture isolates actor presentation and
the command-to-HTTP-byte crossing without repeating minimization.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from zipfile import ZIP_DEFLATED, ZipFile

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
CAPTURES = ROOT / "out" / "browser-acceptance" / "slice6"
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "InspectAndReattach"))
from test_browser_presentation import (  # noqa: E402
    CHROMIUM, DIST, account, completed, controlled_account_server, inspection,
)
from test_export import capture as connected_capture, choose_stage  # noqa: E402 — same connected view checks


@lru_cache(maxsize=None)
def bundle(stage_id: str = "stage-one") -> bytes:
    buffer = BytesIO()
    with ZipFile(buffer, "w", ZIP_DEFLATED) as archive:
        archive.writestr("manifest.json", json.dumps({"stageId": stage_id}))
        archive.writestr("coordinates.cif", f"data_{stage_id.replace('-', '_')}\n#\n")
    return buffer.getvalue()


BUNDLE = bundle("stage-one")
DIGEST = sha256(BUNDLE).hexdigest()


def assessment(stage_id: str, check_standing: str, reason: str) -> dict:
    return {"id": f"assessment-{stage_id}", "stageId": stage_id,
            "checkStanding": check_standing, "reason": reason,
            "evidence": [], "findings": [], "limitations": ["Controlled status"],
            "currentlyApplicable": True}


def stage_account(stage_id: str, origin: str, check_standing: str, reason: str,
                  constructed: dict) -> dict:
    return {"stageId": stage_id, "attemptId": "attempt-one",
            "studyRevisionId": origin, "kind": "Minimization", "status": "completed",
            "assessment": assessment(stage_id, check_standing, reason),
            "summary": f"Completed minimized stage {stage_id}",
            "observation": None, "constructed": constructed, "export": None,
            "originProteinLabel": "1CRN · coordinate model 1 · deposited coordinates · chains A"
            if stage_id == "stage-one" else
            "1UBQ · coordinate model 1 · deposited coordinates · chains A",
            "originMethodLabel": "PACKMOL-Memgen 1.0",
            "runStartedAt": "2026-09-28T09:00:00Z" if stage_id == "stage-one"
            else "2026-09-27T11:00:00Z"}


def export_account(stage_id: str, status: str, reason: str | None = None) -> dict:
    payload = bundle(stage_id)
    return {"stageId": stage_id, "assessmentId": f"assessment-{stage_id}",
            "status": status, "reason": reason,
            "sha256": sha256(payload).hexdigest() if status == "verified" else None,
            "byteLength": len(payload) if status == "verified" else None}


def selected_stage_account() -> dict:
    value = completed(account())
    value["stages"][0]["assessment"].update(id="assessment-stage-one", stageId="stage-one")
    value["inspection"] = inspection("stage-one", "completedStage",
                                     value["stages"][0]["assessment"])
    value["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
    value["inspection"]["omittedMolecules"] = []
    value["stages"][0]["export"] = None
    value["stages"][0].update(originProteinLabel="1CRN · coordinate model 1 · deposited coordinates · chains A",
                               originMethodLabel="PACKMOL-Memgen 1.0",
                               runStartedAt="2026-09-28T09:00:00Z")
    old_constructed = deepcopy(value["attempt"]["constructed"])
    old_constructed.update(subjectId="constructed-zero", attemptId="attempt-zero")
    older = stage_account("stage-two", "revision-zero", "issuesFound",
                                         "A separate observed condition failed",
                                         old_constructed)
    older["attemptId"] = "attempt-zero"
    value["stages"].append(older)
    value["actions"].extend([
        {"kind": "exportStage", "subjectId": "stage-one", "enabled": True,
         "reason": None},
        {"kind": "exportStage", "subjectId": "stage-two", "enabled": False,
         "reason": "Select this completed stage before export."},
    ])
    return value


class ControlledExport:
    def __init__(self, host):
        self.host = host
        self.commands: list[tuple[str, str]] = []
        self.gets: list[tuple[str, str | None]] = []
        self.command_failure_once = False
        self.select_failure_once = False
        self.stale_revision_once = False
        self.revision_checks: list[tuple[int | None, int]] = []
        self.get_failure_once: str | None = None

    def _replace(self, value: dict):
        value["revision"] += 1
        self.host.replace(value)

    def command(self, route):
        request = json.loads(route.request.post_data or "{}")
        kind = request.get("kind")
        subject = request.get("data", {}).get("stageId") or request.get("data", {}).get("subjectId")
        if kind not in {"exportStage", "selectInspectionSubject"}:
            route.continue_()
            return
        self.commands.append((kind, subject))
        with self.host.lock:
            if self.stale_revision_once:
                self.host.account["revision"] += 1
                self.stale_revision_once = False
            current = deepcopy(self.host.account)
        self.revision_checks.append((request.get("expectedRevision"), current["revision"]))
        if request.get("expectedRevision") != current["revision"]:
            route.fulfill(status=422, content_type="application/json",
                          body=json.dumps({"reason": "The workspace changed."}))
            return
        stage = next((item for item in current["stages"] if item["stageId"] == subject), None)
        if stage is None:
            route.fulfill(status=422, content_type="application/json",
                          body=json.dumps({"reason": "The selected stage is unavailable."}))
            return
        if kind == "selectInspectionSubject":
            if self.select_failure_once:
                self.select_failure_once = False
                route.fulfill(status=422, content_type="application/json",
                              body=json.dumps({"reason": "This result could not be opened."}))
                return
            current["inspection"] = inspection(subject, "completedStage", stage["assessment"])
            current["inspection"]["studyRevisionId"] = stage["studyRevisionId"]
            current["inspection"]["studyRevisionNumber"] = 1 if subject == "stage-one" else 0
            current["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
            current["inspection"]["omittedMolecules"] = []
            for item in current["actions"]:
                if item["kind"] == "exportStage":
                    item["enabled"] = item["subjectId"] == subject
                    item["reason"] = None if item["enabled"] else "Select this completed stage before export."
            self._replace(current)
            route.fulfill(status=200, content_type="application/json", body=json.dumps(current))
            return
        if current["inspection"]["subjectId"] != subject:
            route.fulfill(status=422, content_type="application/json",
                          body=json.dumps({"reason": "Select this completed stage before export."}))
            return
        if self.command_failure_once:
            self.command_failure_once = False
            stage["export"] = export_account(subject, "failed", "Bundle correspondence not verified.")
            self._replace(current)
            route.fulfill(status=422, content_type="application/json",
                          body=json.dumps({"reason": "Bundle correspondence not verified."}))
            return
        stage["export"] = export_account(subject, "verified")
        self._replace(current)
        route.fulfill(status=200, content_type="application/json", body=json.dumps(current))

    def download(self, route):
        stage_id = route.request.url.split("/api/export/", 1)[-1]
        if_match = route.request.headers.get("if-match")
        self.gets.append((stage_id, if_match))
        with self.host.lock:
            current = deepcopy(self.host.account)
        if self.get_failure_once == "server":
            self.get_failure_once = None
            stage = next(item for item in current["stages"] if item["stageId"] == stage_id)
            stage["export"] = export_account(stage_id, "failed", "Bundle bytes changed after validation.")
            self._replace(current)
            route.fulfill(status=409, content_type="application/json",
                          body=json.dumps({"reason": "Bundle bytes changed after validation."}))
            return
        if self.get_failure_once == "missing":
            self.get_failure_once = None
            route.fulfill(status=404, content_type="application/json",
                          body=json.dumps({"reason": "The verified bundle is unavailable."}))
            return
        body = bundle(stage_id)
        digest = sha256(body).hexdigest()
        if if_match != f'"{digest}"':
            route.fulfill(status=412, content_type="application/json",
                          body=json.dumps({"reason": "The bundle identity changed."}))
            return
        if self.get_failure_once == "changed-response":
            self.get_failure_once = None
            body += b"changed after HTTP validation"
        route.fulfill(status=200, body=body, headers={
            "Content-Type": "application/zip",
            "Content-Disposition": f'attachment; filename="protein-membrane-{stage_id}.zip"',
            "Cache-Control": "no-store", "ETag": f'"{digest}"',
            "X-Content-SHA256": digest,
        })


class ExportBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (DIST / "index.html").is_file():
            raise AssertionError("Build browser/ before Slice 6 controlled browser tests")

    def open_page(self, playwright, base: str, width: int, height: int):
        browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                             args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                   "--use-angle=swiftshader"])
        page = browser.new_page(viewport={"width": width, "height": height}, accept_downloads=True)
        return browser, page

    def capture(self, page, name: str, width: int):
        expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
        CAPTURES.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(CAPTURES / f"controlled-{name}-{width}.png"),
                        animations="disabled")

    def assert_unchanged_stage(self, page, stage_id="stage-one"):
        expect(page.get_by_label("Completed stage information")).to_contain_text("Verified protein + bilayer + water and ions")
        expect(page.get_by_label("Export validation and unchanged stage review"))\
            .to_contain_text("Required evidence unavailable")
        expect(page.locator(".stage-strip")).to_contain_text("Minimized")
        expect(page.locator(".stage-strip")).to_contain_text("Checks incomplete")
        with self.subTest("scientific account untouched"):
            value = page.request.get(page.url + "api/state").json()
            stage = next(item for item in value["stages"] if item["stageId"] == stage_id)
            self.assertEqual(stage["status"], "completed")
            self.assertEqual(stage["assessment"]["checkStanding"], "checksIncomplete")
            self.assertTrue(stage["assessment"]["currentlyApplicable"])

    def test_collapsed_review_geometry_at_required_viewports(self):
        """Controlled stage checks the connected route's layout assertions cheaply."""
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(selected_stage_account())
            browser, page = self.open_page(playwright, base, 1672, 941)
            try:
                page.goto(base, wait_until="domcontentloaded")
                for width, height in ((1672, 941), (1024, 768), (820, 760)):
                    connected_capture(page, CAPTURES, "controlled-review-ready", width, height,
                                      "stage-one", "assessment-stage-one", failed=False)
            finally:
                browser.close()

    def test_unfinished_attempts_have_no_completed_stage_export_action(self):
        cases = (
            ("pending", False),
            ("running", False),
            ("running", True),
            ("failed", True),
            ("stopped", True),
            ("unobserved", True),
        )
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            browser, first_page = self.open_page(playwright, base, 820, 760)
            first_page.close()
            try:
                for status, has_constructed_candidate in cases:
                    with self.subTest(status=status):
                        value = account()
                        value["attempt"].update(status=status,
                                                progress=None if status != "running" else 0.35)
                        if not has_constructed_candidate:
                            value["attempt"]["constructed"] = None
                            value["inspection"] = None
                        host.replace(value)
                        controlled = ControlledExport(host)
                        page = browser.new_page(viewport={"width": 820, "height": 760},
                                                accept_downloads=True)
                        try:
                            page.route("**/api/commands", controlled.command)
                            page.route("**/api/export/*", controlled.download)
                            page.goto(base, wait_until="domcontentloaded")
                            expect(page.get_by_label("System stages and attempts"))\
                                .to_contain_text("No minimized or equilibrated system stage has completed.")
                            expect(page.get_by_role("button", name="Export this completed stage"))\
                                .to_have_count(0)
                            expect(page.get_by_role("button", name="Export with status"))\
                                .to_have_count(0)
                            self.assertEqual(page.request.get(base + "/api/state").json()["stages"], [])
                            self.assertEqual(controlled.commands, [])
                            self.assertEqual(controlled.gets, [])
                        finally:
                            page.close()
            finally:
                browser.close()

    def test_absent_or_noncurrent_assessment_disables_selected_stage_export(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            browser, first_page = self.open_page(playwright, base, 820, 760)
            first_page.close()
            try:
                for reason in ("absent", "noncurrent"):
                    with self.subTest(reason=reason):
                        value = selected_stage_account()
                        stage = value["stages"][0]
                        if reason == "absent":
                            stage["assessment"] = None
                            value["inspection"]["assessment"] = None
                        else:
                            stage["assessment"]["currentlyApplicable"] = False
                            value["inspection"]["assessment"]["currentlyApplicable"] = False
                        export = next(item for item in value["actions"] if item["kind"] == "exportStage"
                                      and item["subjectId"] == "stage-one")
                        export.update(enabled=False,
                                      reason="The selected completed stage has no currently applicable assessment.")
                        host.replace(value)
                        controlled = ControlledExport(host)
                        page = browser.new_page(viewport={"width": 820, "height": 760},
                                                accept_downloads=True)
                        try:
                            page.route("**/api/commands", controlled.command)
                            page.route("**/api/export/*", controlled.download)
                            page.goto(base, wait_until="domcontentloaded")
                            expect(page.get_by_label("Completed stage information"))\
                                .to_contain_text("Stage information")
                            expect(page.get_by_role("button", name="Export this completed stage"))\
                                .to_be_disabled()
                            expect(page.get_by_role("button", name="Export with status"))\
                                .to_have_count(0)
                            expect(page.locator(".stage-export"))\
                                .to_contain_text("no currently applicable assessment")
                            self.assertEqual(controlled.commands, [])
                            self.assertEqual(controlled.gets, [])
                        finally:
                            page.close()
            finally:
                browser.close()

    def test_stale_actor_revision_refuses_before_get_and_refreshes_for_same_stage_retry(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(selected_stage_account())
            controlled = ControlledExport(host)
            browser, page = self.open_page(playwright, base, 820, 760)
            try:
                page.route("**/api/commands", controlled.command)
                page.route("**/api/export/*", controlled.download)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
                controlled.stale_revision_once = True
                page.get_by_role("button", name="Export with status").click()
                failure = page.get_by_label("Export validation and unchanged stage review")
                expect(failure).to_contain_text("The workspace changed.")
                self.assert_unchanged_stage(page)
                self.assertEqual(controlled.commands, [("exportStage", "stage-one")])
                self.assertNotEqual(*controlled.revision_checks[0])
                self.assertEqual(controlled.gets, [], "A stale command cannot request bundle bytes")
                with page.expect_download(timeout=30000) as pending:
                    page.get_by_role("button", name="Retry export").click()
                with tempfile.TemporaryDirectory() as temporary:
                    saved = Path(temporary) / "export.zip"
                    pending.value.save_as(saved)
                    self.assertEqual(sha256(saved.read_bytes()).hexdigest(), DIGEST)
                expect(failure).to_have_count(0)
                self.assertEqual(controlled.commands,
                                 [("exportStage", "stage-one"), ("exportStage", "stage-one")])
                self.assertEqual(*controlled.revision_checks[1])
                self.assertEqual(controlled.gets, [("stage-one", f'"{DIGEST}"')])
            finally:
                browser.close()

    def test_command_refusal_then_same_stage_retry_and_historical_selection(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(selected_stage_account())
            controlled = ControlledExport(host)
            controlled.command_failure_once = True
            browser, page = self.open_page(playwright, base, 1672, 941)
            try:
                page.route("**/api/commands", controlled.command)
                page.route("**/api/export/*", controlled.download)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
                self.capture(page, "export-ready", 1672)
                page.get_by_role("button", name="Export with status").click()
                failure = page.get_by_label("Export validation and unchanged stage review")
                expect(failure).to_contain_text("Export not delivered")
                expect(failure).to_contain_text("Bundle correspondence not verified")
                expect(page.get_by_label("Details")).to_be_visible()
                self.capture(page, "export-failure", 1672)
                self.assert_unchanged_stage(page)
                self.assertEqual(controlled.gets, [], "A refused command cannot begin ZIP delivery")
                expect(page.get_by_role("button", name="Retry export")).to_be_enabled()
                choose_stage(page, "stage-two")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.get_by_label("Export validation and unchanged stage review"))\
                    .to_have_count(0)
                expect(page.get_by_label("Minimized stage review and technical checks"))\
                    .to_contain_text("A separate observed condition failed")
                expect(page.get_by_label("Completed stage information"))\
                    .to_contain_text("This result uses earlier inputs")
                choose_stage(page, "stage-one")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(failure).to_contain_text("Bundle correspondence not verified")
                with page.expect_download(timeout=30000) as pending:
                    page.get_by_role("button", name="Retry export").click()
                download = pending.value
                self.assertEqual(download.suggested_filename, "protein-membrane-stage-one.zip")
                with tempfile.TemporaryDirectory() as temporary:
                    saved = Path(temporary) / "export.zip"
                    download.save_as(saved)
                    self.assertEqual(sha256(saved.read_bytes()).hexdigest(), DIGEST)
                expect(page.get_by_label("Export validation and unchanged stage review"))\
                    .to_have_count(0)
                expect(page.get_by_label("Minimized stage review and technical checks"))\
                    .to_contain_text("Required evidence unavailable")
                self.assertEqual(controlled.commands, [
                    ("exportStage", "stage-one"),
                    ("selectInspectionSubject", "stage-two"),
                    ("selectInspectionSubject", "stage-one"),
                    ("exportStage", "stage-one"),
                ])
                self.assertEqual(controlled.gets, [("stage-one", f'"{DIGEST}"')])
            finally:
                browser.close()

    def test_refused_stage_switch_keeps_each_export_fault_with_its_inspected_stage(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = selected_stage_account()
            value["stages"][0]["export"] = export_account("stage-one", "failed", "First result fault.")
            value["stages"][1]["export"] = export_account("stage-two", "failed", "Second result fault.")
            host.replace(value)
            controlled = ControlledExport(host)
            browser, page = self.open_page(playwright, base, 1024, 800)
            try:
                page.route("**/api/commands", controlled.command)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Results").click()
                choices = page.get_by_label("Completed stage choices")
                controlled.select_failure_once = True
                choices.locator('[data-stage-id="stage-two"]').click()
                expect(page.get_by_label("Completed stage information")).to_contain_text("1CRN")
                detail = page.get_by_label("Export validation and unchanged stage review")
                expect(detail).to_contain_text("First result fault.")
                expect(detail).not_to_contain_text("Second result fault.")
                expect(page.locator(".stage-export")).to_contain_text("Second result fault.")
                self.assertEqual(host.account["inspection"]["subjectId"], "stage-one")

                choices.locator('[data-stage-id="stage-two"]').click()
                expect(page.get_by_label("Completed stage information")).to_contain_text("1UBQ")
                controlled.select_failure_once = True
                choices.locator('[data-stage-id="stage-one"]').click()
                expect(detail).to_contain_text("Second result fault.")
                expect(detail).not_to_contain_text("First result fault.")
                expect(page.locator(".stage-export")).to_contain_text("First result fault.")
                self.assertEqual(host.account["inspection"]["subjectId"], "stage-two")
            finally:
                browser.close()

    def test_two_similar_earlier_results_keep_origin_and_exact_export_selection(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = selected_stage_account()
            earlier = value["stages"][1]
            earlier["assessment"] = assessment("stage-two", "checksPassed", "Controlled checks passed.")
            third = deepcopy(earlier)
            third.update(stageId="stage-three", attemptId="attempt-three",
                         studyRevisionId="revision-three", summary="Another earlier minimized result",
                         originProteinLabel="2L6W · coordinate model 1 · deposited coordinates · chains A",
                         runStartedAt="2026-09-26T11:00:00Z")
            third["assessment"] = assessment("stage-three", "checksPassed", "Controlled checks passed.")
            third["constructed"] = deepcopy(earlier["constructed"])
            third["constructed"].update(subjectId="constructed-three", attemptId="attempt-three")
            value["stages"].append(third)
            value["actions"].append({"kind": "exportStage", "subjectId": "stage-three",
                                      "enabled": False, "reason": "Select this completed stage before export."})
            host.replace(value)
            controlled = ControlledExport(host)
            browser, page = self.open_page(playwright, base, 1024, 800)
            try:
                page.route("**/api/commands", controlled.command)
                page.route("**/api/export/*", controlled.download)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Results").click()
                choices = page.get_by_label("Completed stage choices")
                expect(choices.locator('[data-stage-id="stage-two"]')).to_contain_text("1UBQ")
                expect(choices.locator('[data-stage-id="stage-three"]')).to_contain_text("2L6W")
                for stage_id, source in (("stage-two", "1UBQ"), ("stage-three", "2L6W")):
                    choices.locator(f'[data-stage-id="{stage_id}"]').click()
                    expect(page.get_by_label("Completed stage information")).to_contain_text(source)
                    expect(page.get_by_label("Completed stage information"))\
                        .to_contain_text("PACKMOL-Memgen 1.0")
                    with page.expect_download(timeout=30000) as pending:
                        page.get_by_role("button", name="Export this completed stage").click()
                    with tempfile.TemporaryDirectory() as temporary:
                        saved = Path(temporary) / "result.zip"
                        pending.value.save_as(saved)
                        with ZipFile(saved) as archive:
                            self.assertEqual(json.loads(archive.read("manifest.json"))["stageId"], stage_id)
                    self.assertEqual(host.account["inspection"]["subjectId"], stage_id)
                self.assertIn(("exportStage", "stage-two"), controlled.commands)
                self.assertIn(("exportStage", "stage-three"), controlled.commands)
                self.assertEqual([stage_id for stage_id, _ in controlled.gets],
                                 ["stage-two", "stage-three"])
            finally:
                browser.close()

    def test_http_and_byte_faults_keep_retry_reachable_at_reduced_size(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(selected_stage_account())
            controlled = ControlledExport(host)
            controlled.get_failure_once = "server"
            browser, page = self.open_page(playwright, base, 820, 760)
            try:
                page.route("**/api/commands", controlled.command)
                page.route("**/api/export/*", controlled.download)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
                self.capture(page, "export-ready", 820)
                page.get_by_role("button", name="Export with status").click()
                failure = page.get_by_label("Export validation and unchanged stage review")
                expect(failure).to_contain_text("Bundle bytes changed after validation")
                expect(page.get_by_label("Details")).to_be_visible()
                self.capture(page, "export-failure", 820)
                self.assert_unchanged_stage(page)
                retry = page.get_by_role("button", name="Retry export")
                expect(retry).to_be_visible()
                self.assertIsNotNone(retry.bounding_box())
                self.assertLessEqual(retry.bounding_box()["y"] + retry.bounding_box()["height"],
                                     760, "Retry must stay in the reduced viewport")
                self.assertEqual(page.locator(".stage-strip").bounding_box()["y"] +
                                 page.locator(".stage-strip").bounding_box()["height"] <= 762,
                                 True)
                controlled.get_failure_once = "changed-response"
                retry.click()
                expect(failure).to_contain_text("downloaded bytes did not match", timeout=30000)
                self.assert_unchanged_stage(page)
                with page.expect_download(timeout=30000) as pending:
                    retry.click()
                with tempfile.TemporaryDirectory() as temporary:
                    saved = Path(temporary) / "export.zip"
                    pending.value.save_as(saved)
                    self.assertEqual(sha256(saved.read_bytes()).hexdigest(), DIGEST)
                expect(failure).to_have_count(0)
                self.assertEqual([subject for subject, _ in controlled.gets],
                                 ["stage-one", "stage-one", "stage-one"])
                self.assertEqual([subject for kind, subject in controlled.commands if kind == "exportStage"],
                                 ["stage-one", "stage-one", "stage-one"])
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
