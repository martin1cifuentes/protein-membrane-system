"""One connected published-browser route: Build and minimize, reattach, export/read-back.

This is the only long 6QWR/DMPC actor route. It selects the qualified named
native recipe and uses the actual host, worker, source, browser and saved ZIP.
Controlled presentation cases in sibling slices cover faster state variants.
"""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "ConstructAndMinimizeExplicitSystem"))
from test_construct_and_minimize import (  # noqa: E402
    HOST, POLICY, PPM, SOURCE, await_completed_stage, await_state, chromium,
    current_state, prepare_adopted_alkl_dmpc, running_host, start_named_route,
)
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "InspectAndReattachToCurrentWork"))
from test_reattach import atom_site_elements, click_molstar_atom, selected_atom  # noqa: E402

CAPTURES = ROOT / "out" / "browser-acceptance" / "connected-native"
INSPECTOR = Path(__file__).with_name("inspect_bundle.py")
NATIVE_POLICY_ID = "canonical-protein-pure-dmpc-native-construction"
NATIVE_ROUTE_LABEL = "OpenMM native: DMPC"
STANDINGS = {"checksPassed", "issuesFound", "checksIncomplete"}


def stage_in(state: dict, stage_id: str) -> dict:
    return next(item for item in state["stages"] if item["stageId"] == stage_id)


def choose_stage(page, stage_id: str):
    page.locator(f'.stage-card[data-stage-id="{stage_id}"]').click()


@contextmanager
def retained_run_directory():
    # Keep exact host state and artifacts on a late browser assertion failure.
    # The route is expensive, and a retained completed stage can be reattached.
    yield tempfile.mkdtemp(prefix="real-native-", dir=CAPTURES)


def capture(page, folder: Path, name: str, width: int, height: int, stage_id: str,
            assessment_id: str, *, failed: bool):
    page.set_viewport_size({"width": width, "height": height})
    if page.locator(".workspace.workflow-closed").count() == 0:
        page.get_by_role("button", name="Collapse inputs").click()
    expect(page.locator(".workspace.workflow-closed.execution-review")).to_have_count(1)
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
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    toolbar = page.locator(".scene-inspection-toolbar").bounding_box()
    scene_label = page.locator(".execution-scene-label").bounding_box()
    assert all(box is not None for box in (scene, evidence, stages, toolbar, scene_label))
    folder.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(folder / f"{name}-{width}.png"), full_page=True,
                    animations="disabled")
    (folder / f"{name}-{width}-layout.json").write_text(json.dumps({
        "fixtureKind": "connected real route" if not name.startswith("controlled-") else "controlled fixture",
        "stageId": stage_id, "assessmentId": assessment_id,
        "viewport": {"width": width, "height": height},
        "scene": scene, "evidence": evidence, "stages": stages,
        "toolbar": toolbar, "sceneLabel": scene_label,
    }, indent=2) + "\n")
    assert scene["width"] >= (350 if width <= 820 else 600 if width <= 1150 else width * .64)
    assert evidence["width"] >= (330 if width <= 820 else 390 if width <= 1150 else width * .25)
    assert evidence["x"] >= scene["x"] + scene["width"] - 2
    assert abs(evidence["y"] - scene["y"]) <= 3
    assert scene_label["y"] >= toolbar["y"] + toolbar["height"] - 2
    assert stages["y"] + stages["height"] <= height + 2
    if width > 1150:
        scene_share = scene["width"] / (scene["width"] + evidence["width"])
        assert .64 <= scene_share <= .77
    state = current_state(page)
    stage = stage_in(state, stage_id)
    assert stage["status"] == "completed"
    assert stage["assessment"]["id"] == assessment_id
    assert stage["assessment"]["checkStanding"] in STANDINGS
    assert state["inspection"]["subjectId"] == stage_id
    if failed:
        expect(page.locator(".export-validation")).to_contain_text("Export not delivered")
        expect(page.get_by_role("button", name="Inspect export issue").first).to_be_visible()
        expect(page.get_by_role("button", name="Retry export")).to_be_visible()
        # A changed published ZIP is rehashed by the Host at delivery and
        # invalidates only this stage's export, not its scientific assessment.
        assert stage["export"]["status"] == "failed"
    else:
        expect(page.get_by_role("button", name="Export with status")).to_be_enabled()


class ConnectedBuildReattachExportRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, POLICY, PPM, SOURCE, INSPECTOR,
                     ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Connected browser prerequisite unavailable: {path}")
        CAPTURES.mkdir(parents=True, exist_ok=True)

    def test_one_real_native_build_reattach_export_fault_retry_and_readback(self):
        with retained_run_directory() as temporary:
            workspace = Path(temporary) / "workspace"
            downloaded = Path(temporary) / "downloaded-stage.zip"
            print(f"Connected native route workspace: {temporary}", flush=True)
            with running_host(workspace, POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941},
                                            accept_downloads=True)
                    downloads = []
                    page.on("download", lambda transfer: downloads.append(transfer))
                    page.goto(base, wait_until="domcontentloaded")
                    adopted = prepare_adopted_alkl_dmpc(page)
                    revision_id = adopted["study"]["id"]
                    admitted = start_named_route(page, NATIVE_ROUTE_LABEL, NATIVE_POLICY_ID)
                    attempt_id = admitted["attempt"]["attemptId"]
                    self.assertEqual(admitted["attempt"]["studyRevisionId"], revision_id)
                    self.assertEqual(admitted["attempt"]["policyId"], NATIVE_POLICY_ID)
                    self.assertEqual(admitted["stages"], [])
                    duplicate = page.request.post(base + "/api/commands", data=json.dumps({
                        "kind": "buildAndMinimize", "data": {"policyId": NATIVE_POLICY_ID},
                        "expectedRevision": current_state(page)["revision"]}), headers={
                            "Content-Type": "application/json", "Origin": base})
                    self.assertEqual(duplicate.status, 422, duplicate.text())
                    self.assertEqual(current_state(page)["attempt"]["attemptId"], attempt_id)
                    print(f"One authorized native attempt: {attempt_id}", flush=True)

                    # The host-owned attempt survives browser closure and an interrupted
                    # live connection. Reopening must not start a second attempt.
                    page.close()
                    page = browser.new_page(viewport={"width": 1672, "height": 941},
                                            accept_downloads=True)
                    page.on("download", lambda transfer: downloads.append(transfer))
                    blocked_events = []
                    def block_events(route):
                        blocked_events.append(route.request.url)
                        route.abort()
                    page.route("**/api/events", block_events)
                    page.goto(base, wait_until="domcontentloaded")
                    reopened = await_state(page, lambda state: state["attempt"] is not None and
                                           state["attempt"]["attemptId"] == attempt_id,
                                           "same attempt after browser reattachment", timeout=120)
                    self.assertEqual(reopened["study"]["id"], revision_id)
                    self.assertIn(reopened["attempt"]["status"], ("running", "completed"))
                    if reopened["attempt"]["status"] == "running":
                        self.assertEqual(reopened["stages"], [])
                    else:
                        self.assertTrue(any(stage["attemptId"] == attempt_id and
                                            stage["status"] == "completed"
                                            for stage in reopened["stages"]))
                    self.assertEqual(len(list((workspace / "attempts").iterdir())), 1)
                    expect(page.get_by_role("alert").filter(
                        has_text="The live connection is interrupted")).to_be_visible(timeout=10000)
                    self.assertTrue(blocked_events)
                    page.unroute("**/api/events")
                    page.reload(wait_until="domcontentloaded")

                    # A checked constructed account remains an intermediate under
                    # the one upfront authorization, with no midpoint action/stage.
                    handoff = await_state(page, lambda state: state["attempt"] is not None and
                                          state["attempt"]["attemptId"] == attempt_id and
                                          state["attempt"].get("constructed") is not None,
                                          "checked construction handoff", timeout=1200)
                    constructed = handoff["attempt"]["constructed"]
                    self.assertEqual(constructed["attemptId"], attempt_id)
                    self.assertFalse(any(item["kind"] == "continueMinimization"
                                         for item in handoff["actions"]))
                    if handoff["attempt"]["status"] == "running":
                        self.assertEqual(handoff["stages"], [])
                    self.assertTrue((workspace / "attempts" / attempt_id /
                                     "constructed-correspondence.json").is_file())
                    print("Checked constructed intermediate observed under same attempt", flush=True)

                    completed, stage = await_completed_stage(page, attempt_id)
                    stage_id = stage["stageId"]
                    assessment = stage["assessment"]
                    self.assertEqual(stage["studyRevisionId"], revision_id)
                    self.assertEqual(stage["constructed"]["subjectId"], constructed["subjectId"])
                    self.assertIsNotNone(stage["observation"])
                    self.assertEqual(stage["observation"]["termination"], "converged")
                    final_force = next(value for value in stage["observation"]["measurements"]
                                       if value["name"] == "finalRmsForce" and
                                       value["scope"] ==
                                       "final unrestrained constraint tangent; per particle")
                    self.assertLessEqual(final_force["value"], 10.0)
                    self.assertIn(assessment["checkStanding"], STANDINGS)
                    self.assertTrue(assessment["reason"])
                    self.assertTrue(assessment["currentlyApplicable"])
                    print(f"Completed {stage_id}: {assessment['checkStanding']}", flush=True)

                    # Select and reopen the actual completed stage. Coordinate-row
                    # identity is read from the saved independent correspondence file.
                    choose_stage(page, stage_id)
                    selected = await_state(page, lambda state: state["inspection"] is not None and
                                           state["inspection"]["subjectId"] == stage_id,
                                           "exact completed-stage inspection", timeout=120)
                    self.assertEqual(selected["inspection"]["studyRevisionId"], revision_id)
                    self.assertEqual(selected["inspection"]["assessment"]["id"], assessment["id"])
                    page.close()
                    page = browser.new_page(viewport={"width": 1672, "height": 941},
                                            accept_downloads=True)
                    page.on("download", lambda transfer: downloads.append(transfer))
                    export_commands = []
                    page.on("request", lambda request: export_commands.append(
                        request.post_data_json) if request.url.endswith("/api/commands") and
                        request.post_data else None)
                    page.goto(base, wait_until="domcontentloaded")
                    attached = await_state(page, lambda state: state["attempt"] is not None and
                                           state["attempt"]["attemptId"] == attempt_id and
                                           any(item["stageId"] == stage_id for item in state["stages"]),
                                           "same completed stage after fresh browser open", timeout=120)
                    self.assertEqual(attached["study"]["id"], revision_id)
                    self.assertEqual(len(list((workspace / "attempts").iterdir())), 1)
                    if attached["inspection"] is None or attached["inspection"]["subjectId"] != stage_id:
                        choose_stage(page, stage_id)
                    attached = await_state(page, lambda state: state["inspection"] is not None and
                                           state["inspection"]["subjectId"] == stage_id,
                                           "reattached exact completed stage", timeout=120)
                    exact_stage = page.get_by_label("Completed stage information").locator("details")
                    exact_stage.locator("summary").click()
                    expect(exact_stage).to_contain_text(revision_id)
                    mapping = json.loads((workspace / "attempts" / attempt_id /
                                          "constructed-correspondence.json").read_text())
                    structure = page.request.get(base + attached["inspection"]["structureUrl"])
                    self.assertTrue(structure.ok)
                    elements = atom_site_elements(structure.body())
                    self.assertEqual(len(elements), len(mapping["atoms"]))
                    chosen = {}
                    for atom in mapping["atoms"]:
                        role = atom["moleculeRole"]
                        key = (role, atom.get("physicalSide") if role == "lipid" else None)
                        chosen.setdefault(key, atom)
                    for key in (("protein", None), ("lipid", "upper"),
                                ("lipid", "lower"), ("water", None), ("ion", None)):
                        expected = chosen[key]
                        self.assertEqual(selected_atom(base, page, attached,
                                                       expected["resultAtomIndex"]), expected)
                        self.assertEqual(elements[expected["resultAtomIndex"]], expected["element"])
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
                    click_molstar_atom(page, chosen[("protein", None)])
                    click_molstar_atom(page, chosen[("lipid", "upper")])
                    capture(page, Path(temporary), "export-ready", 1672, 941, stage_id,
                            assessment["id"], failed=False)
                    capture(page, Path(temporary), "export-ready", 820, 760, stage_id,
                            assessment["id"], failed=False)
                    print("Reattached stage, atom identity and browser view checked", flush=True)

                    # Corrupt the exact ZIP after the owner publishes it and before
                    # the browser's first GET. The stage/assessment must survive.
                    tampered = {"path": None, "requests": 0}
                    def corrupt_first_delivery(route):
                        tampered["requests"] += 1
                        if tampered["path"] is None:
                            zip_paths = list((workspace / "exports").rglob("*.zip"))
                            self.assertEqual(len(zip_paths), 1)
                            tampered["path"] = zip_paths[0]
                            with zip_paths[0].open("rb") as bundle:
                                expected_digest = hashlib.file_digest(bundle, "sha256").hexdigest()
                            self.assertEqual(route.request.headers.get("if-match"),
                                             f'"{expected_digest}"')
                            with zip_paths[0].open("ab") as bundle:
                                bundle.write(b"\x00")
                        route.continue_()
                    page.route("**/api/export/*", corrupt_first_delivery)
                    page.get_by_role("button", name="Export with status").click()
                    failed = await_state(page, lambda state: (
                        (stage_in(state, stage_id).get("export") or {}).get("status") == "failed"),
                        "exact-stage export delivery integrity failure", timeout=180)
                    expect(page.locator(".export-validation")).to_contain_text(
                        "Export not delivered", timeout=30000)
                    expect(page.get_by_role("button", name="Retry export")).to_be_visible()
                    self.assertEqual(tampered["requests"], 1)
                    self.assertEqual(downloads, [])
                    self.assertEqual(stage_in(failed, stage_id)["export"]["status"], "failed")
                    self.assertEqual(stage_in(failed, stage_id)["assessment"]["id"], assessment["id"])
                    self.assertEqual(failed["attempt"]["currentStageId"], stage_id)
                    capture(page, Path(temporary), "export-failed", 1672, 941, stage_id,
                            assessment["id"], failed=True)
                    capture(page, Path(temporary), "export-failed", 820, 760, stage_id,
                            assessment["id"], failed=True)
                    page.unroute("**/api/export/*", corrupt_first_delivery)
                    with page.expect_download(timeout=600000) as transfer:
                        page.get_by_role("button", name="Retry export").click()
                    download = transfer.value
                    self.assertEqual(download.suggested_filename,
                                     f"protein-membrane-{stage_id}.zip")
                    download.save_as(downloaded)
                    delivered = await_state(page, lambda state: (
                        (stage_in(state, stage_id).get("export") or {}).get("status") == "verified"),
                        "verified same-stage export retry", timeout=180)
                    delivered_stage = stage_in(delivered, stage_id)
                    self.assertEqual(delivered_stage["assessment"]["id"], assessment["id"])
                    self.assertEqual(delivered_stage["export"]["stageId"], stage_id)
                    self.assertEqual(delivered_stage["export"]["assessmentId"], assessment["id"])
                    self.assertEqual(delivered_stage["export"]["sha256"],
                                     hashlib.sha256(downloaded.read_bytes()).hexdigest())
                    self.assertEqual(delivered_stage["export"]["byteLength"], downloaded.stat().st_size)
                    self.assertEqual([item["data"]["stageId"] for item in export_commands
                                      if item["kind"] == "exportStage"], [stage_id, stage_id])
                    self.assertEqual(page.request.get(base + f"/api/export/{stage_id}",
                                      headers={"If-Match": '"wrong-bundle-digest"'}).status, 412)
                    self.assertEqual(page.request.get(base + f"/api/export/{stage_id}", headers={
                        "If-Match": f'"{delivered_stage["export"]["sha256"]}"',
                        "Origin": "http://unrelated.invalid"}).status, 403)
                    self.assertEqual(page.request.get(base + "/api/export/unknown-stage").status, 404)
                    print("Same-stage export retry and HTTP identity checked", flush=True)

                    # Controlled browser cases cover historical selection and
                    # structure-view recovery. Keep this costly route focused on
                    # the real delivered bundle and independent read-back.
                finally:
                    try:
                        with urlopen(base + "/api/state", timeout=10) as response:
                            account_snapshot = json.load(response)
                        (Path(temporary) / "final-account.json").write_text(
                            json.dumps(account_snapshot, indent=2) + "\n")
                    except Exception as failure:
                        (Path(temporary) / "account-read-error.txt").write_text(str(failure) + "\n")
                    browser.close()

            # The downloaded result is inspected outside the host workspace,
            # after its browser and server have stopped.
            checked = subprocess.run([
                str(ROOT / "out" / "python" / "bin" / "python"), str(INSPECTOR),
                str(downloaded), stage_id, attempt_id, revision_id, assessment["id"],
                assessment["checkStanding"]],
                cwd=ROOT, text=True, capture_output=True, timeout=600, check=False)
            self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
            evidence = json.loads(checked.stdout)
            self.assertEqual(evidence["stageId"], stage_id)
            self.assertEqual(evidence["attemptId"], attempt_id)
            self.assertEqual(evidence["studyRevisionId"], revision_id)
            self.assertEqual(evidence["checkStanding"], assessment["checkStanding"])
            self.assertEqual(evidence["atomCount"], stage["constructed"]["atomCount"])
            (Path(temporary) / "independent-bundle-readback.json").write_text(
                json.dumps(evidence, indent=2) + "\n")
            print("Downloaded bundle independently read back after host stop", flush=True)


if __name__ == "__main__":
    unittest.main()
