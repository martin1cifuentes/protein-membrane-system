"""Slice 5 browser reattachment through the published host and actual worker.

Run after scripts/build-local.sh with the identified Slice 4 catalogue and
PIM_BROWSER_CHROMIUM set to the local Chromium executable. The long positive
performs one real 6QWR/DMPC construction and minimization.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from urllib.parse import urlsplit
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "ConstructAndMinimizeExplicitSystem"))
from test_construct_and_minimize import (  # noqa: E402
    HOST, POLICY, PPM, SOURCE, await_reviewable_candidate, await_state,
    chromium, current_state, prepare_adopted_alkl_dmpc, running_host,
)

CAPTURES = ROOT / "out" / "browser-acceptance" / "slice5"


def milestone(label: str):
    print(f"Slice 5 route: {label}", flush=True)


def capture(page, name: str, width: int, height: int, expected_status: str):
    page.set_viewport_size({"width": width, "height": height})
    page.wait_for_timeout(150)
    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
    expect(page.locator(".scene-error")).to_have_count(0)
    expect(page.locator(".execution-scene-label")).to_be_visible()
    expect(page.locator(".execution-decision")).to_be_visible()
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    assert scene is not None and evidence is not None and stages is not None
    assert scene["width"] >= (350 if width <= 820 else width * .64)
    assert evidence["width"] >= (300 if width <= 820 else width * .25)
    assert evidence["x"] > scene["x"] and abs(evidence["y"] - scene["y"]) <= 3
    assert stages["y"] + stages["height"] <= height + 2
    state = current_state(page)
    assert state["attempt"]["status"] == expected_status
    if expected_status == "running":
        assert state["attempt"]["stageKind"] == "Minimization" and not state["stages"]
    else:
        assert any(stage["kind"] == "Minimization" and stage["status"] == "completed"
                   for stage in state["stages"])
    CAPTURES.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(CAPTURES / f"{name}-{width}.png"), full_page=True,
                    animations="disabled")
    assert current_state(page)["attempt"]["status"] == expected_status


def selected_atom(base: str, page, account: dict, index: int):
    inspection = account["inspection"]
    token = urlsplit(inspection["structureUrl"]).path.rsplit("/", 1)[1]
    response = page.request.get(base + f"/api/inspection/atoms/"
                                f"{inspection['subjectId']}/{token}/{index}")
    assert response.ok, f"Unverified atom row {index}: {response.status} {response.text()}"
    value = response.json()
    assert value["subjectId"] == inspection["subjectId"]
    assert value["studyRevisionId"] == inspection["studyRevisionId"]
    assert value["structureToken"] == token and value["atomSiteIndex"] == index
    assert value["atom"]["resultAtomIndex"] == index
    return value["atom"]


def atom_site_elements(mmcif: bytes):
    lines = mmcif.decode("utf-8").splitlines()
    headers = []
    for line_index, line in enumerate(lines):
        if line.startswith("_atom_site."):
            headers.append(line.strip())
            continue
        if headers and line and not line.startswith("_atom_site."):
            if "_atom_site.type_symbol" not in headers:
                raise AssertionError("Inspected mmCIF has no atom-site element column")
            element_column = headers.index("_atom_site.type_symbol")
            elements = []
            for atom_line in lines[line_index:]:
                if atom_line.startswith(("#", "loop_")):
                    break
                if atom_line.startswith(("ATOM", "HETATM")):
                    elements.append(atom_line.split()[element_column])
            if elements:
                return elements
    raise AssertionError("Inspected mmCIF has no atom-site rows")


def click_molstar_atom(page, expected: dict):
    index = expected["resultAtomIndex"]
    emitted = page.locator(".viewer-mount").evaluate("""(mount, atomSiteIndex) => {
      const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
      const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!structure) return false;
      for (const unit of structure.units) {
        if (!unit.model?.atomicHierarchy?.atomSourceIndex) continue;
        for (let offset = 0; offset < unit.elements.length; offset++) {
          const modelElement = unit.elements[offset];
          if (unit.model.atomicHierarchy.atomSourceIndex.value(modelElement) !== atomSiteIndex) continue;
          const loci = {kind: 'element-loci', structure,
            elements: [{unit, indices: Int32Array.of(offset)}]};
          viewer.plugin.behaviors.interaction.click.next({
            current: {loci, repr: null}, buttons: 1, button: 1,
            modifiers: {alt: false, control: false, meta: false, shift: false},
          });
          return true;
        }
      }
      return false;
    }""", index)
    assert emitted, f"Mol* did not contain expected atom-site row {index}"
    card = page.get_by_role("region", name="Selected atom correspondence")
    expect(card).to_contain_text(expected["resultAtomId"], timeout=30000)
    expect(card).to_contain_text(expected["moleculeRole"])
    if expected["generatedSpeciesId"]:
        expect(card).to_contain_text(expected["generatedSpeciesId"])
    if expected["physicalSide"]:
        expect(card).to_contain_text(expected["physicalSide"])
    if expected["sourceAtomId"]:
        expect(card).to_contain_text(expected["sourceAtomId"])


class InspectAndReattachRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, POLICY, PPM, SOURCE,
                     ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Slice 5 browser prerequisite unavailable: {path}")
        CAPTURES.mkdir(parents=True, exist_ok=True)

    def test_real_running_and_completed_reattachment(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "workspace"
            with running_host(workspace, POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    adopted = prepare_adopted_alkl_dmpc(page)
                    origin_id = adopted["study"]["id"]
                    if page.locator(".workspace.workflow-closed").count():
                        page.get_by_role("button", name="Show workflow").click()
                    page.get_by_role("button", name="Construct system").click()
                    candidate = await_reviewable_candidate(page)
                    attempt_id = candidate["attempt"]["attemptId"]
                    constructed_id = candidate["attempt"]["constructed"]["subjectId"]
                    page.get_by_role("button", name="Continue minimization").first.click()
                    await_state(page, lambda state: state["attempt"] is not None and
                                state["attempt"]["attemptId"] == attempt_id and
                                state["attempt"]["stageKind"] == "Minimization" and
                                state["attempt"]["status"] == "running" and
                                not state["stages"], "identified live minimization", timeout=120)
                    current = current_state(page)
                    if current["inspection"] is None or current["inspection"]["subjectId"] != constructed_id:
                        page.get_by_role("button", name="Inspect verified constructed system").click()
                    inspected = await_state(page, lambda state: state["inspection"] is not None and
                                            state["inspection"]["subjectId"] == constructed_id,
                                            "exact constructed inspection during minimization")
                    self.assertEqual(inspected["inspection"]["studyRevisionId"], origin_id)
                    self.assertEqual(inspected["inspection"]["representationKind"], "constructedSystem")
                    self.assertIsNotNone(inspected["inspection"]["structureUrl"])
                    capture(page, "running-before-reattach", 1672, 941, "running")
                    capture(page, "running-before-reattach", 820, 760, "running")

                    # Browser presentation ends; the identified host and worker survive.
                    page.close()
                    interrupted = browser.new_page(viewport={"width": 1672, "height": 941})
                    interrupted.add_init_script("""(() => {
                      window.__slice5Sse = {opens: 0, messages: 0};
                      const NativeEventSource = window.EventSource;
                      window.EventSource = class extends NativeEventSource {
                        constructor(...args) {
                          super(...args);
                          window.__slice5Sse.opens += 1;
                          this.addEventListener('message', () => { window.__slice5Sse.messages += 1; });
                        }
                      };
                    })()""")
                    blocked_events = []
                    def block_events(route):
                        blocked_events.append(route.request.url)
                        route.abort()
                    interrupted.route("**/api/events", block_events)
                    interrupted.goto(base, wait_until="domcontentloaded")
                    reacquired = await_state(interrupted, lambda state: state["attempt"] is not None and
                                             state["attempt"]["attemptId"] == attempt_id and
                                             state["attempt"]["status"] == "running" and
                                             state["attempt"]["stageKind"] == "Minimization" and
                                             state["inspection"] is not None and
                                             state["inspection"]["subjectId"] == constructed_id and
                                             not state["stages"],
                                             "same running attempt from initial state without SSE", timeout=120)
                    self.assertEqual(reacquired["inspection"]["studyRevisionId"], origin_id)
                    interrupted.wait_for_function(
                        "() => window.__slice5Sse?.opens > 0", timeout=10000)
                    self.assertTrue(blocked_events, "The actor's SSE request must actually be interrupted")
                    expect(interrupted.locator(".attempt-account")).to_contain_text(attempt_id)
                    expect(interrupted.locator(".inspection-origin-account")).to_contain_text(origin_id)
                    expect(interrupted.locator(".execution-scene-label")).to_contain_text(
                        "Current constructed system")
                    expect(interrupted.locator(".stage-strip")).to_contain_text("None")
                    self.assertEqual(len(list((workspace / "attempts").iterdir())), 1,
                                     "Reopening a view must not start a second attempt")
                    capture(interrupted, "running-sse-interrupted", 1672, 941, "running")
                    capture(interrupted, "running-sse-interrupted", 820, 760, "running")
                    interrupted.unroute("**/api/events")
                    interrupted.reload(wait_until="domcontentloaded")
                    reconnected = await_state(interrupted, lambda state: state["attempt"] is not None and
                                              state["attempt"]["attemptId"] == attempt_id and
                                              state["inspection"] is not None and
                                              state["inspection"]["subjectId"] == constructed_id,
                                              "same attempt after SSE reconnect", timeout=120)
                    self.assertEqual(reconnected["study"]["id"], origin_id)
                    self.assertEqual(reconnected["inspection"]["studyRevisionId"], origin_id)
                    expect(interrupted.locator(".attempt-account")).to_contain_text(attempt_id)
                    expect(interrupted.locator(".inspection-origin-account")).to_contain_text(origin_id)
                    mapping_path = workspace / "attempts" / attempt_id / "constructed-correspondence.json"
                    mapping = json.loads(mapping_path.read_text())
                    structure = interrupted.request.get(
                        base + reconnected["inspection"]["structureUrl"])
                    self.assertTrue(structure.ok)
                    elements = atom_site_elements(structure.body())
                    self.assertEqual(len(elements), len(mapping["atoms"]))
                    chosen = {}
                    for expected in mapping["atoms"]:
                        role = expected["moleculeRole"]
                        key = (role, expected.get("physicalSide") if role == "lipid" else None)
                        if key not in chosen:
                            chosen[key] = expected
                    for key in (("protein", None), ("lipid", "upper"),
                                ("lipid", "lower"), ("water", None), ("ion", None)):
                        expected = chosen[key]
                        index = expected["resultAtomIndex"]
                        actual = selected_atom(base, interrupted, reconnected, index)
                        self.assertEqual(actual, expected,
                                         f"Mol* coordinate row {index} must resolve exact {key} identity")
                        self.assertEqual(actual["element"], elements[index])
                    click_molstar_atom(interrupted, chosen[("protein", None)])
                    click_molstar_atom(interrupted, chosen[("lipid", "upper")])

                    completed = await_state(interrupted, lambda state: any(
                        stage["attemptId"] == attempt_id and stage["kind"] == "Minimization" and
                        stage["status"] == "completed" for stage in state["stages"]),
                        "observed completed minimized stage", timeout=2400)
                    stage = next(stage for stage in completed["stages"]
                                 if stage["attemptId"] == attempt_id and stage["kind"] == "Minimization")
                    self.assertEqual(completed["attempt"]["currentStageId"], stage["stageId"])
                    self.assertEqual(completed["attempt"]["status"], "completed")
                    self.assertEqual(stage["studyRevisionId"], origin_id)
                    self.assertIsNotNone(stage["assessment"])
                    interrupted.wait_for_function(
                        "() => window.__slice5Sse?.messages > 0", timeout=30000)
                    interrupted.locator(".stage-card").filter(has_text=stage["stageId"]).click()
                    before_close = await_state(interrupted, lambda state: state["inspection"] is not None and
                                               state["inspection"]["subjectId"] == stage["stageId"],
                                               "completed stage selected before browser close", timeout=120)
                    self.assertEqual(before_close["inspection"]["studyRevisionId"], origin_id)
                    self.assertEqual(before_close["inspection"]["assessment"]["id"],
                                     stage["assessment"]["id"])
                    expect(interrupted.locator(".execution-identity-account")).to_contain_text(
                        stage["stageId"])
                    capture(interrupted, "minimized-before-reattach", 1672, 941, "completed")
                    capture(interrupted, "minimized-before-reattach", 820, 760, "completed")
                    milestone("completed stage captured before reopening")
                    interrupted.close()

                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    reopened = await_state(page, lambda state: state["attempt"] is not None and
                                           state["attempt"]["attemptId"] == attempt_id and
                                           any(item["stageId"] == stage["stageId"] for item in state["stages"]),
                                           "completed stage after fresh browser open", timeout=120)
                    self.assertEqual(reopened["attempt"]["status"], "completed")
                    self.assertEqual(reopened["study"]["id"], origin_id)
                    expect(page.locator(".stage-strip")).to_contain_text(stage["stageId"])
                    if reopened["inspection"] is None or \
                            reopened["inspection"]["subjectId"] != stage["stageId"]:
                        page.locator(".stage-card").filter(has_text=stage["stageId"]).click()
                    selected = await_state(page, lambda state: state["inspection"] is not None and
                                           state["inspection"]["subjectId"] == stage["stageId"],
                                           "selected exact completed stage", timeout=120)
                    self.assertEqual(selected["inspection"]["studyRevisionId"], origin_id)
                    self.assertEqual(selected["inspection"]["assessment"]["id"],
                                     stage["assessment"]["id"])
                    stage_elements = atom_site_elements(page.request.get(
                        base + selected["inspection"]["structureUrl"]).body())
                    self.assertEqual(len(stage_elements), len(mapping["atoms"]))
                    revision_before_display = selected["revision"]
                    page.locator(".execution-display-tool").get_by_role("button", name="Display").click()
                    options = page.locator("#execution-display-options")
                    expect(options).to_contain_text("Full model:")
                    options.get_by_label("Hidden").check()
                    expect(page.locator(".execution-display-disclosure")).to_contain_text("water hidden")
                    options.get_by_label("Na⁺ / Cl⁻ ions").uncheck()
                    expect(page.locator(".execution-display-disclosure")).to_contain_text("ions hidden")
                    self.assertEqual(current_state(page)["revision"], revision_before_display,
                                     "Molecule display choices must not alter scientific state")
                    options.get_by_label("All").check()
                    options.get_by_label("Na⁺ / Cl⁻ ions").check()
                    for expected in chosen.values():
                        index = expected["resultAtomIndex"]
                        actual = selected_atom(base, page, selected, index)
                        self.assertEqual(actual, expected,
                                         "Completed-stage atom identity must retain every source/generated field")
                        self.assertEqual(actual["element"], stage_elements[index])
                        click_molstar_atom(page, expected)
                    options.get_by_label("Representative sample").check()
                    page.locator(".execution-display-tool").get_by_role("button", name="Display").click()
                    expect(page.locator(".execution-identity-account")).to_contain_text(origin_id)
                    capture(page, "minimized-after-reattach", 1672, 941, "completed")
                    capture(page, "minimized-after-reattach", 820, 760, "completed")
                    milestone("completed stage captured after reopening")

                    # A new source choice is a new study revision. The completed stage
                    # remains inspectable and labelled with its original revision.
                    page.set_viewport_size({"width": 1672, "height": 941})
                    if page.locator(".workspace.workflow-closed").count():
                        page.get_by_role("button", name="Show workflow").click()
                    page.locator("#source-upload").set_input_files(str(SOURCE))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Inspect uploaded source").click()
                    revised = await_state(page, lambda state: state["study"] is not None and
                                          state["study"]["id"] != origin_id,
                                          "new study revision from new source choice", timeout=180)
                    milestone("new source revision established")
                    self.assertTrue(any(item["stageId"] == stage["stageId"] for item in revised["stages"]))
                    page.get_by_role("button", name="Review attempt").click()
                    page.get_by_role("button", name="Inspect verified constructed system").click()
                    earlier_candidate = await_state(page, lambda state: state["inspection"] is not None and
                                                    state["inspection"]["subjectId"] == constructed_id and
                                                    state["inspection"]["studyRevisionId"] == origin_id,
                                                    "earlier constructed subject retained after study revision",
                                                    timeout=120)
                    self.assertNotEqual(earlier_candidate["study"]["id"], origin_id)
                    self.assertEqual(earlier_candidate["inspection"]["representationKind"],
                                     "constructedSystem")
                    self.assertTrue(page.request.get(
                        base + earlier_candidate["inspection"]["structureUrl"]).ok)
                    expect(page.locator(".inspection-origin-account")).to_contain_text("Historical subject")
                    milestone("historical constructed subject reselected")
                    page.locator(".stage-card").filter(has_text=stage["stageId"]).click()
                    historical = await_state(page, lambda state: state["inspection"] is not None and
                                             state["inspection"]["subjectId"] == stage["stageId"] and
                                             state["inspection"]["studyRevisionId"] == origin_id,
                                             "historical stage with its originating revision", timeout=120)
                    self.assertNotEqual(historical["study"]["id"], historical["inspection"]["studyRevisionId"])
                    expect(page.locator(".execution-identity-account")).to_contain_text(
                        "Historical stage")
                    expect(page.locator(".execution-identity-account")).to_contain_text(
                        f"current revision {historical['study']['number']}")
                    self.assertEqual(historical["inspection"]["assessment"]["id"],
                                     stage["assessment"]["id"])
                    self.assertEqual(
                        [item["id"] for item in historical["inspection"]["evidence"]],
                        [item["id"] for item in selected["inspection"]["evidence"]],
                        "Reselecting a historical stage must preserve evidence identities")
                    self.assertEqual(
                        [item["id"] for item in historical["inspection"]["findings"]],
                        [item["id"] for item in selected["inspection"]["findings"]],
                        "Reselecting a historical stage must preserve finding identities")
                    self.assertEqual(selected_atom(base, page, historical, 0)["moleculeRole"], "protein")
                    milestone("historical stage and evidence continuity verified")
                    page.route("**/api/structures/*",
                               lambda route: route.fulfill(status=404, body="Structure unavailable"))
                    page.reload(wait_until="domcontentloaded")
                    expect(page.locator(".scene-error")).to_contain_text("Structure unavailable",
                                                                       timeout=180000)
                    milestone("unavailable view reported")
                    expect(page.locator(".execution-identity-account")).to_contain_text(origin_id)
                    expect(page.locator(".execution-assessment-account")).to_contain_text(
                        stage["assessment"]["reason"])
                    self.assertEqual(current_state(page)["stages"][0]["stageId"], stage["stageId"])
                    page.unroute("**/api/structures/*")
                    page.reload(wait_until="domcontentloaded")
                    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
                    expect(page.locator(".scene-error")).to_have_count(0)
                    milestone("structure restored after view failure")
                    (CAPTURES / "final-account.json").write_text(
                        json.dumps(historical, indent=2) + "\n")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
