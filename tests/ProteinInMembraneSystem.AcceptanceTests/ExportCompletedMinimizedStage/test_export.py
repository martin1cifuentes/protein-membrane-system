"""Slice 6 completed-stage export through the published browser, host and worker.

Run after scripts/build-local.sh with the installed identified OpenMM/PPM assets.
This one long route constructs and minimizes the real 6QWR/DMPC system once.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "ConstructAndMinimizeExplicitSystem"))
from test_construct_and_minimize import (  # noqa: E402
    HOST, POLICY, PPM, SOURCE, await_reviewable_candidate, await_state,
    chromium, current_state, prepare_adopted_alkl_dmpc, running_host,
)

CAPTURES = ROOT / "out" / "browser-acceptance" / "slice6"
INSPECTOR = Path(__file__).with_name("inspect_bundle.py")


def milestone(label: str):
    print(f"Slice 6 route: {label}", flush=True)


def selected_stage(state: dict, stage_id: str) -> dict:
    return next(item for item in state["stages"] if item["stageId"] == stage_id)


def capture(page, name: str, width: int, height: int, stage_id: str,
            assessment_id: str, *, failed: bool):
    page.set_viewport_size({"width": width, "height": height})
    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
    expect(page.locator(".scene-error")).to_have_count(0)
    page.wait_for_function("""() => {
        const mount = document.querySelector('.viewer-mount');
        if (!mount || mount.dataset.cameraReady !== 'true') return false;
        const bounds = mount.getBoundingClientRect();
        return Math.abs(Number(mount.dataset.cameraViewportWidth) - bounds.width) <= 4 &&
            Math.abs(Number(mount.dataset.cameraViewportHeight) - bounds.height) <= 4;
    }""", timeout=30000)
    expect(page.locator(".execution-scene-label")).to_be_visible()
    expect(page.locator(".stage-strip")).to_contain_text("Minimized")
    expect(page.locator(".stage-strip")).to_contain_text("Indeterminate")
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    assert scene is not None and evidence is not None and stages is not None
    assert scene["width"] >= (350 if width <= 820 else width * .64)
    assert evidence["width"] >= (300 if width <= 820 else width * .25)
    assert evidence["x"] > scene["x"] and abs(evidence["y"] - scene["y"]) <= 3
    assert stages["y"] + stages["height"] <= height + 2
    state = current_state(page)
    stage = selected_stage(state, stage_id)
    assert stage["status"] == "completed"
    assert stage["assessment"]["id"] == assessment_id
    assert stage["assessment"]["qualification"] == "Indeterminate"
    assert state["inspection"]["subjectId"] == stage_id
    if failed:
        expect(page.locator(".export-validation")).to_contain_text("Export not delivered")
        expect(page.get_by_role("button", name="Inspect export issue")).to_be_visible()
        expect(page.get_by_role("button", name="Retry export")).to_be_visible()
        assert stage["export"]["status"] == "failed"
    else:
        expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
    CAPTURES.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(CAPTURES / f"{name}-{width}.png"), full_page=True,
                    animations="disabled")


class ExportCompletedMinimizedStageRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, POLICY, PPM, SOURCE, INSPECTOR,
                     ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Slice 6 browser prerequisite unavailable: {path}")
        CAPTURES.mkdir(parents=True, exist_ok=True)

    def test_real_export_fault_and_retry_same_completed_stage(self):
        with tempfile.TemporaryDirectory(prefix="real-route-", dir=CAPTURES) as temporary:
            workspace = Path(temporary) / "workspace"
            downloaded = Path(temporary) / "downloaded-stage.zip"
            milestone(f"retained diagnostic workspace {temporary}")
            with running_host(workspace, POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941},
                                            accept_downloads=True)
                    downloads = []
                    page.on("download", lambda download: downloads.append(download))
                    page.goto(base, wait_until="domcontentloaded")
                    adopted = prepare_adopted_alkl_dmpc(page)
                    revision_id = adopted["study"]["id"]
                    milestone(f"adopted exact study revision {revision_id}")
                    if page.locator(".workspace.workflow-closed").count():
                        page.get_by_role("button", name="Show workflow").click()
                    page.get_by_role("button", name="Construct system").click()
                    candidate = await_reviewable_candidate(page)
                    attempt_id = candidate["attempt"]["attemptId"]
                    milestone(f"native construction candidate {attempt_id} ready for review")
                    self.assertFalse(candidate["stages"])
                    page.get_by_role("button", name="Continue minimization").first.click()
                    milestone(f"required minimization requested for {attempt_id}")
                    completed = await_state(page, lambda state: any(
                        item["attemptId"] == attempt_id and item["kind"] == "Minimization" and
                        item["status"] == "completed" for item in state["stages"]),
                        "factually completed minimized stage", timeout=7200)
                    stage = next(item for item in completed["stages"]
                                 if item["attemptId"] == attempt_id and
                                 item["kind"] == "Minimization")
                    stage_id = stage["stageId"]
                    assessment = stage["assessment"]
                    self.assertIsNotNone(assessment)
                    self.assertEqual(assessment["qualification"], "Indeterminate")
                    self.assertEqual(stage["studyRevisionId"], revision_id)
                    self.assertEqual(completed["attempt"]["currentStageId"], stage_id)
                    if completed["inspection"] is None or \
                            completed["inspection"]["subjectId"] != stage_id:
                        page.locator(".stage-card").filter(has_text=stage_id).click()
                    await_state(page, lambda state: state["inspection"] is not None and
                                state["inspection"]["subjectId"] == stage_id,
                                "selected exact minimized stage")
                    capture(page, "export-ready", 1672, 941, stage_id,
                            assessment["id"], failed=False)
                    capture(page, "export-ready", 820, 760, stage_id,
                            assessment["id"], failed=False)
                    milestone("actual completed stage and export-ready view captured")

                    # Corrupt exactly the ZIP that the real owner just published,
                    # before the browser's first GET reaches the real host.
                    tampered = {"path": None, "requests": 0}

                    def corrupt_first_delivery(route):
                        tampered["requests"] += 1
                        if tampered["path"] is None:
                            zip_paths = list((workspace / "exports").rglob("*.zip"))
                            self.assertEqual(len(zip_paths), 1,
                                             "The command must publish one identified ZIP")
                            tampered["path"] = zip_paths[0]
                            with zip_paths[0].open("rb") as bundle:
                                expected_digest = hashlib.file_digest(bundle, "sha256").hexdigest()
                            self.assertEqual(route.request.headers.get("if-match"),
                                             f'"{expected_digest}"',
                                             "The GET must bind the command's verified bundle")
                            with zip_paths[0].open("ab") as bundle:
                                bundle.write(b"\x00")
                        route.continue_()

                    page.route("**/api/export/*", corrupt_first_delivery)
                    page.get_by_role("button", name="Export with status").click()
                    failed = await_state(page, lambda state: (
                        (selected_stage(state, stage_id).get("export") or {}).get("status") == "failed"),
                        "stage-bound export delivery fault", timeout=180)
                    self.assertEqual(tampered["requests"], 1)
                    self.assertIsNotNone(tampered["path"])
                    self.assertEqual(downloads, [], "A corrupt ZIP must not start a browser download")
                    self.assertEqual(selected_stage(failed, stage_id)["assessment"]["id"],
                                     assessment["id"])
                    self.assertEqual(failed["attempt"]["currentStageId"], stage_id)
                    expect(page.locator(".export-validation")).to_contain_text(
                        "Export not delivered", timeout=30000)
                    capture(page, "export-failed", 1672, 941, stage_id,
                            assessment["id"], failed=True)
                    capture(page, "export-failed", 820, 760, stage_id,
                            assessment["id"], failed=True)
                    milestone("corrupt host ZIP refused; scientific stage unchanged")

                    page.unroute("**/api/export/*", corrupt_first_delivery)
                    with page.expect_download(timeout=600000) as transfer:
                        page.get_by_role("button", name="Retry export").click()
                    download = transfer.value
                    self.assertEqual(download.suggested_filename,
                                     f"protein-membrane-{stage_id}.zip")
                    download.save_as(downloaded)
                    self.assertTrue(downloaded.is_file())
                    delivered = await_state(page, lambda state: (
                        (selected_stage(state, stage_id).get("export") or {}).get("status") == "verified"),
                        "verified retry for the same completed stage", timeout=180)
                    delivered_stage = selected_stage(delivered, stage_id)
                    self.assertEqual(delivered_stage["assessment"]["id"], assessment["id"])
                    self.assertEqual(delivered_stage["export"]["assessmentId"], assessment["id"])
                    self.assertEqual(delivered_stage["export"]["stageId"], stage_id)
                    self.assertEqual(delivered_stage["export"]["sha256"],
                                     hashlib.sha256(downloaded.read_bytes()).hexdigest())
                    self.assertEqual(delivered_stage["export"]["byteLength"],
                                     downloaded.stat().st_size)
                    self.assertEqual(delivered["inspection"]["subjectId"], stage_id)
                    self.assertEqual(delivered["attempt"]["attemptId"], attempt_id)
                    wrong_digest = page.request.get(
                        base + f"/api/export/{stage_id}",
                        headers={"If-Match": '"wrong-bundle-digest"'})
                    self.assertEqual(wrong_digest.status, 412)
                    foreign_origin = page.request.get(
                        base + f"/api/export/{stage_id}",
                        headers={
                            "If-Match": f'"{delivered_stage["export"]["sha256"]}"',
                            "Origin": "http://unrelated.invalid",
                        })
                    self.assertEqual(foreign_origin.status, 403)
                    self.assertEqual(page.request.get(base + "/api/export/unknown-stage").status,
                                     404)
                    self.assertEqual(selected_stage(current_state(page), stage_id)["export"]["status"],
                                     "verified", "A wrong client digest must not invalidate the bundle")
                    milestone("same-stage retry downloaded verified bundle")
                finally:
                    browser.close()

            # The saved file remains usable after the browser and host are gone.
            checked = subprocess.run([
                str(ROOT / "out" / "python" / "bin" / "python"), str(INSPECTOR),
                str(downloaded), stage_id, attempt_id, revision_id, assessment["id"]],
                cwd=ROOT, text=True, capture_output=True, timeout=600, check=False)
            self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
            evidence = json.loads(checked.stdout)
            self.assertEqual(evidence["stageId"], stage_id)
            self.assertEqual(evidence["attemptId"], attempt_id)
            self.assertEqual(evidence["studyRevisionId"], revision_id)
            self.assertEqual(evidence["atomCount"], stage["constructed"]["atomCount"])
            (CAPTURES / "independent-bundle-readback.json").write_text(
                json.dumps(evidence, indent=2) + "\n")
            milestone("downloaded bundle independently read back after host stopped")


if __name__ == "__main__":
    unittest.main()
