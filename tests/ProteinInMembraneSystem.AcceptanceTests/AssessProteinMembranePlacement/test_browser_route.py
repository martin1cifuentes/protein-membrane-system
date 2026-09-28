"""Slice 3 through the published launcher, real worker, local PPM and browser.

Build with scripts/build-local.sh. Then run with the Playwright test environment:

    PIM_BROWSER_CHROMIUM=/path/to/chrome out/browser-test-python/bin/python \
      -m unittest discover \
      -s tests/ProteinInMembraneSystem.AcceptanceTests/AssessProteinMembranePlacement \
      -p 'test_*.py' -v

The qualified combined catalogue may be overridden with
PIM_PLACEMENT_POLICY_CATALOGUE. Captures are kept in out/browser-acceptance/slice3.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import tempfile
import time
import unittest
from urllib.error import URLError
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
HOST = Path(os.environ.get("PIM_BROWSER_HOST", str(ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"))).resolve()
WORKER = ROOT / "out" / "python" / "bin" / "python"
LAUNCHER = ROOT / "scripts" / "start-local.sh"
SOURCE = ROOT / "config" / "policies" / "source-assets" / "6QWR.pdb"
POLICY = Path(os.environ.get(
    "PIM_PLACEMENT_POLICY_CATALOGUE",
    str(ROOT / "config" / "policies" / "protein-membrane-slice3.json"),
)).resolve()
PPM = Path(os.environ.get("PIM_PPM_EXECUTABLE", str(ROOT / "out" / "ppm2" / "immers"))).resolve()
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "slice3"


@contextmanager
def running_host(workspace: Path, catalogue: Path, ppm: Path | None = PPM):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    environment = dict(os.environ, PIM_PORT=str(port), ASPNETCORE_URLS=base,
                       PIM_WORKSPACE_ROOT=str(workspace),
                       PIM_POLICY_CATALOGUE=str(catalogue),
                       PIM_WORKER_PYTHON=str(WORKER),
                       PIM_OWNER_SOURCE_ROOT=str(ROOT / "src" / "ProteinInMembrane.Host"))
    if ppm is None:
        environment.pop("PIM_PPM_EXECUTABLE", None)
    else:
        environment["PIM_PPM_EXECUTABLE"] = str(ppm)
    log_path = workspace.parent / "placement-host.log"
    with log_path.open("wb") as log:
        command = ["dotnet", str(HOST)] if os.environ.get("PIM_BROWSER_HOST") else [str(LAUNCHER)]
        process = subprocess.Popen(command, cwd=ROOT, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT)
        try:
            for _ in range(150):
                if process.poll() is not None:
                    raise AssertionError("Local launcher exited: " + log_path.read_text())
                try:
                    with urlopen(base + "/api/state", timeout=2) as response:
                        if response.status == 200:
                            break
                except URLError:
                    time.sleep(0.1)
            else:
                raise AssertionError("Local host did not become ready: " + log_path.read_text())
            yield base
        except BaseException:
            print("\nHost output tail:\n" + log_path.read_text()[-2500:])
            raise
        finally:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)


def chromium(playwright):
    executable = os.environ.get("PIM_BROWSER_CHROMIUM")
    return playwright.chromium.launch(
        executable_path=executable or None, headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"],
    )


def current_state(page):
    return page.evaluate("async () => (await fetch('/api/state', {cache: 'no-store'})).json()")


def await_state(page, predicate, label: str, timeout: float = 180):
    deadline = time.monotonic() + timeout
    state = None
    while time.monotonic() < deadline:
        state = current_state(page)
        if predicate(state):
            return state
        time.sleep(0.15)
    summary = None if state is None else {
        "revision": state["revision"],
        "protein": None if state["protein"] is None else {
            "status": state["protein"]["status"],
            "changes": [(item["id"], item["kind"], item["residue"]["residue"])
                        for item in state["protein"]["changes"]],
        },
        "membrane": None if state["membrane"] is None else {
            "status": state["membrane"]["status"], "reason": state["membrane"].get("reason")},
        "placement": None if state["placement"] is None else {
            "status": state["placement"]["status"], "reason": state["placement"].get("reason")},
        "notices": state["notices"][-3:],
    }
    raise AssertionError(f"Timed out waiting for {label}: {summary}")


def prepare_exact_protein(page):
    page.get_by_role("navigation", name="Research work areas")\
        .get_by_role("button", name="Protein").click()
    page.locator("#source-upload").set_input_files(str(SOURCE))
    page.locator("#upload-provenance").select_option("experimental")
    page.get_by_role("button", name="Upload source").click()
    expect(page.locator(".source-context")).to_contain_text(
        "Structure displayed", timeout=180000)
    if page.locator("#model-index").count():
        page.locator("#model-index").select_option("0")
    else:
        expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
    page.locator("#assembly-choice").select_option("deposited")
    page.get_by_label("Chain A").check()
    page.get_by_role("button", name="Assess selected protein").click()
    assessed_or_changes = await_state(page,
        lambda state: state["protein"] is not None and (
            state["protein"]["status"] == "assessed" or state["protein"]["changes"]),
        "assessed protein or located change proposals")
    changes = assessed_or_changes["protein"]["changes"]
    for residue_number, expected_change in ((108, "HID"), (211, "OXT")):
        if not any(item["residue"]["residue"] == residue_number and
                   expected_change.lower() in (item["proposedChange"] + " " + item["rationale"]).lower()
                   for item in changes):
            raise AssertionError(f"Expected {expected_change} at residue {residue_number}: " +
                                 str([(item["kind"], item["residue"]["residue"], item["proposedChange"])
                                      for item in changes]))
    for decision in assessed_or_changes["preparationReview"]["decisions"]:
        if decision["standing"] != "pending":
            continue
        kind = decision["kind"]
        address = decision["residue"]
        site = f"{address['residue']} · Chain {address['chain']}"
        if kind in ("residueState", "alternateLocation"):
            page.locator(".review-site-list .review-site").filter(has_text=site).click()
            option = next((item for item in decision["options"]
                           if item["proposedChange"] == "HID"), decision["options"][0])
            page.locator(".review-option").filter(has_text=option["proposedChange"])\
                .locator("input").check()
            button = "Confirm state" if kind == "residueState" else "Confirm conformer"
        elif kind in ("heavyAtom", "disulfide"):
            if page.locator(".review-other").get_attribute("open") is None:
                page.locator(".review-other summary").click()
            page.locator(".review-other .review-site").filter(has_text=site).click()
            option = decision["options"][0]
            page.locator(".review-option").filter(has_text=option["proposedChange"])\
                .locator("input").check()
            button = "Approve repair" if kind == "heavyAtom" else "Confirm possible bond"
        else:
            raise AssertionError(f"Unexpected preparation decision kind: {kind}")
        expect(page.get_by_role("button", name=button)).to_be_enabled()
        page.get_by_role("button", name=button).click()
        await_state(page, lambda state: next((item["standing"] for item in
                    state["preparationReview"]["decisions"] if item["id"] == decision["id"]), None)
                    == "confirmed", f"recorded {kind} at {site}")
    assessed = await_state(page, lambda state: state["protein"] is not None and
                           state["protein"]["status"] == "assessed", "assessed prepared protein")
    if page.locator(".workspace.workflow-closed").count():
        page.get_by_role("button", name="Show inputs").click()
    return assessed


def choose_exact_membrane(page, species_id="DMPC"):
    page.get_by_role("navigation", name="Research work areas")\
        .get_by_role("button", name="Membrane").click()
    for side in ("Upper", "Lower"):
        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option(species_id)
        page.get_by_label(f"{side} leaflet percentage 1", exact=True).fill("100")
    page.get_by_role("button", name="Propose membrane model").click()
    proposed = await_state(page, lambda state: state["membrane"] is not None and
                           state["membrane"]["status"] == "proposed" and
                           all(fraction["speciesId"] == species_id for fraction in
                               state["membrane"]["upper"] + state["membrane"]["lower"]),
                           f"proposed {species_id} membrane")
    page.get_by_role("button", name="Adopt and assess displayed proposal").click()
    assessed = await_state(page, lambda state: state["membrane"] is not None and
                           state["membrane"]["status"] == "assessed" and
                           state["membrane"]["policyId"] and
                           state["membrane"]["modelId"] == proposed["membrane"]["modelId"],
                           f"assessed {species_id} membrane")
    assert assessed["membrane"]["modelId"] == proposed["membrane"]["modelId"]
    return assessed


def request_placement(page):
    page.get_by_role("navigation", name="Research work areas")\
        .get_by_role("button", name="Placement").click()
    page.locator("#topology-kind").select_option("membrane-spanning")
    page.locator("#ppm-nterminal-side").select_option("in")
    page.locator("#biological-sidedness").fill("extracellular-upper; periplasmic-lower")
    expect(page.get_by_role("button", name="Obtain and assess position")).to_be_enabled()
    before = current_state(page)
    page.get_by_role("button", name="Obtain and assess position").click()
    state = await_state(page, lambda account:
                        account["placement"] is not None and
                        account["placement"]["status"] not in ("proposed", "Proposed")
                        and account["placement"]["reason"] or
                        account["revision"] > before["revision"] and
                        any(item["severity"] == "warning" for item in
                            account["notices"][len(before["notices"]):]),
                        "terminal placement account or route failure", timeout=300)
    if state["placement"] is None:
        raise AssertionError("Placement route produced no proposal: " +
                             str([item["message"] for item in
                                  state["notices"][len(before["notices"]):]]))
    return state


def inspect_placement(page, proposal_id: str):
    page_errors = []
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    selected = current_state(page)["inspection"]
    if selected is None or selected["subjectId"] != proposal_id:
        if not page.get_by_role("button", name="View positioned protein").is_visible():
            page.get_by_role("button", name="Show inputs").click()
        page.get_by_role("button", name="View positioned protein").click()
    try:
        expect(page.locator(".workspace")).to_have_class(re.compile(r"\bworkflow-open\b"))
    except AssertionError as failure:
        raise AssertionError(f"Placement review failed to render: {page_errors}; "
                             f"page text: {page.locator('body').inner_text()[:600]}") from failure
    expect(page.locator(".scene-heading h1")).to_contain_text("Protein placement proposal")
    expect(page.locator(".evidence-panel")).to_contain_text(proposal_id)
    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=120000)
    expect(page.locator(".scene-error")).to_have_count(0)
    expect(page.locator(".scene-loading")).to_have_count(0, timeout=120000)


def capture_review(test: unittest.TestCase, page, name: str, width: int, height: int):
    page.set_viewport_size({"width": width, "height": height})
    page.evaluate("() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))")
    page.locator(".evidence-content").evaluate("element => { element.scrollTop = 0; }")
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    test.assertIsNotNone(scene)
    test.assertIsNotNone(evidence)
    test.assertIsNotNone(stages)
    test.assertGreaterEqual(scene["width"], 350 if width <= 820 else width * 0.48)
    test.assertGreaterEqual(evidence["width"], 210 if width <= 820 else width * 0.19)
    test.assertGreaterEqual(scene["height"], height * 0.55)
    test.assertGreater(evidence["x"], scene["x"])
    test.assertAlmostEqual(evidence["y"], scene["y"], delta=3)
    test.assertLessEqual(stages["y"] + stages["height"], height + 2)
    overflowing = page.evaluate("""() => [...document.querySelectorAll(
        '.workspace, .workspace-header, .rail, .scene-panel, .evidence-panel, .evidence-content')]
      .filter(element => element.scrollWidth > element.clientWidth + 1)
      .map(element => ({region: element.className, width: element.clientWidth,
                        contentWidth: element.scrollWidth}))""")
    test.assertEqual(overflowing, [])
    expect(page.get_by_role("img", name=re.compile(r"intended.*bilayer.*positioned protein", re.I))).to_be_visible()
    expect(page.locator(".placement-scene-legend")).to_contain_text("Measured target planes")
    expect(page.locator(".stage-strip")).to_contain_text("None")
    expect(page.locator(".stage-strip")).to_contain_text("No minimized or equilibrated system stage has completed.")
    page.screenshot(path=str(ARTIFACTS / f"{name}-{width}.png"),
                    full_page=True, animations="disabled")
    if width <= 820:
        lower = page.locator(".placement-orientation-list button").filter(has_text="Lower physical side")
        lower.scroll_into_view_if_needed()
        expect(lower).to_be_visible()
        page.get_by_label("Placement standing and uncertainty").scroll_into_view_if_needed()
        expect(page.get_by_label("Placement standing and uncertainty")).to_be_visible()
        page.locator(".evidence-content").evaluate("element => { element.scrollTop = 0; }")


def catalogue_without_placement_policy(source: Path, destination: Path):
    catalogue = json.loads(source.read_text())
    source_base = source.parent

    def absolute_paths(value):
        if isinstance(value, list):
            return [absolute_paths(item) for item in value]
        if isinstance(value, dict):
            return {key: str((source_base / item).resolve())
                    if key in {"path", "templatePath", "coordinateTemplatePath", "coordinatePath",
                               "paePath", "mappingPath", "ppmResidueLibraryPath"} and isinstance(item, str) and item and
                    not Path(item).is_absolute() else absolute_paths(item)
                    for key, item in value.items()}
        return value

    catalogue = absolute_paths(catalogue)
    catalogue["placementPolicies"] = []
    catalogue["placementWitnesses"] = []
    catalogue["version"] += "-test-no-placement-policy"
    destination.write_text(json.dumps(catalogue, indent=2) + "\n")


class BrowserPlacementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, WORKER, LAUNCHER, SOURCE, POLICY, PPM,
                     HOST.parent / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Placement route prerequisite unavailable: {path}")
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_exact_real_candidate_can_be_inspected_and_only_supported_one_adopted(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    protein = prepare_exact_protein(page)["protein"]["subjectId"]
                    membrane = choose_exact_membrane(page)["membrane"]["modelId"]
                    state = request_placement(page)
                    placement = state["placement"]
                    self.assertEqual(placement["status"], "supported", placement["reason"])
                    self.assertEqual(placement["preparedProteinId"], protein)
                    self.assertEqual(placement["membraneModelId"], membrane)
                    self.assertEqual(placement["physicalSide"], "both")
                    self.assertEqual(placement["policyId"], "alkl-6qwr-dmpc-spanning-placement")
                    self.assertEqual(placement["policyVersion"], "1.0.0")
                    self.assertEqual(placement["witnessId"],
                                     "alkl-6qwr-model1-dmpc-regional-topology")
                    self.assertIsNotNone(placement["midplaneAngstrom"])
                    self.assertGreater(placement["thicknessAngstrom"], 0)
                    self.assertTrue(placement["evidence"])
                    self.assertTrue(any("PPM" in item["source"] for item in placement["evidence"]))
                    self.assertTrue(any(item["bearing"].lower() == "supports" and
                                        item["method"] == "Independently witnessed placement relationship"
                                        for item in placement["evidence"]),
                                    "PPM context alone cannot establish support")
                    self.assertTrue(any("implicit" in item["uncertainty"].lower()
                                        for item in placement["evidence"]))
                    self.assertEqual(state["stages"], [])
                    inspect_placement(page, placement["proposalId"])
                    inspected = current_state(page)["inspection"]
                    upper_anchor = next(item for item in inspected["annotations"]
                                        if item["label"] == "upper side")
                    self.assertIsNotNone(upper_anchor["geometryFocus"])
                    self.assertTrue(any(item["evidenceId"] == upper_anchor["evidenceId"]
                                        for item in inspected["metrics"]))
                    page.locator(".placement-orientation-list button").filter(
                        has_text="Upper physical side").click()
                    focused = await_state(page, lambda value: value["inspection"] is not None and
                                          value["inspection"]["focusId"] == upper_anchor["subjectPartId"],
                                          "linked upper-side spatial and numerical focus")
                    self.assertEqual(focused["inspection"]["focus"]["authSeqId"],
                                     upper_anchor["geometryFocus"]["authSeqId"])
                    expect(page.locator(".scene-local-state")).to_have_count(1)
                    local_radius = page.evaluate("""() => document.querySelector('.viewer-mount')
                      ?.[Symbol.for('molstar.viewer')]?.plugin.canvas3d?.camera.getSnapshot().radius""")
                    page.get_by_role("button", name="Clear selected part").click()
                    await_state(page, lambda value: value["inspection"] is not None and
                                value["inspection"]["focusId"] is None, "cleared spatial focus")
                    expect(page.locator(".scene-local-state")).to_have_count(0)
                    page.wait_for_function("""localRadius => document.querySelector('.viewer-mount')
                      ?.[Symbol.for('molstar.viewer')]?.plugin.canvas3d?.camera.getSnapshot().radius > localRadius * 2""",
                      arg=local_radius)
                    panel = page.locator(".evidence-panel")
                    expect(panel).to_contain_text(protein)
                    expect(panel).to_contain_text(membrane)
                    expect(panel).to_contain_text("PPM")
                    expect(panel).to_contain_text("DOPC")
                    expect(panel).to_contain_text("physical", ignore_case=True)
                    expect(panel).to_contain_text("supported", ignore_case=True)
                    expect(page.get_by_role("button", name="Adopt supported position")).to_be_enabled()
                    capture_review(self, page, "supported", 1672, 941)
                    capture_review(self, page, "supported", 820, 760)
                    page.get_by_role("button", name="Adopt supported position").click()
                    adopted = await_state(page, lambda value: value["study"]["number"] ==
                                          state["study"]["number"] + 1 and value["placement"] is not None and
                                          value["study"]["adoptedPlacementProposalId"] == placement["proposalId"],
                                          "adopted placement")
                    self.assertEqual(adopted["placement"]["status"], "supported")
                    self.assertEqual(adopted["placement"]["proposalId"], placement["proposalId"])
                    self.assertEqual(adopted["study"]["adoptedPlacementProposalId"], placement["proposalId"])
                    self.assertFalse(next(item["enabled"] for item in adopted["actions"]
                                          if item["kind"] == "adoptPlacement"),
                                     "The same proposal cannot be adopted twice")
                    self.assertIsNone(adopted["attempt"],
                                      "Adopting placement must not itself start construction")
                    self.assertEqual(adopted["stages"], [])
                    self.assertEqual(adopted["inspection"]["subjectId"], placement["proposalId"])
                    expect(page.get_by_label("Position proposal and scientific support"))\
                        .to_contain_text("Adopted in the current study revision")
                    capture_review(self, page, "adopted", 1672, 941)
                    capture_review(self, page, "adopted", 820, 760)
                    page.set_viewport_size({"width": 1672, "height": 941})
                    changed = choose_exact_membrane(page, "POPC")
                    self.assertEqual(changed["study"]["number"], adopted["study"]["number"] + 1)
                    self.assertIsNone(changed["placement"],
                                      "A consequential membrane edit must withdraw current placement reliance")
                    self.assertFalse(next(item["enabled"] for item in changed["actions"]
                                          if item["kind"] == "adoptPlacement"))
                    self.assertFalse(next(item["enabled"] for item in changed["actions"]
                                          if item["kind"] == "startPreparation"))
                finally:
                    browser.close()

    def test_real_reoriented_candidate_contradicts_exact_witness_and_withdraws_old_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    prepare_exact_protein(page)
                    choose_exact_membrane(page)
                    before = request_placement(page)
                    original = before["placement"]
                    self.assertEqual(original["status"], "supported", original["reason"])
                    inspect_placement(page, original["proposalId"])
                    original_url = current_state(page)["inspection"]["structureUrl"]
                    page.get_by_label("Tilt about X (°)").fill("180")
                    page.get_by_role("button", name="Reassess corrected placement").click()
                    revised = await_state(page, lambda state: state["placement"] is not None and
                                          state["placement"]["proposalId"] != original["proposalId"] and
                                          state["placement"]["status"] not in ("proposed", "Proposed"),
                                          "freshly assessed reversed placement", timeout=300)
                    placement = revised["placement"]
                    self.assertEqual(placement["status"], "unsupported", placement["reason"])
                    self.assertEqual(placement["preparedProteinId"], original["preparedProteinId"])
                    self.assertEqual(placement["membraneModelId"], original["membraneModelId"])
                    self.assertEqual(revised["study"]["number"], before["study"]["number"])
                    self.assertTrue(any(item["bearing"].lower() == "contradicts"
                                        for item in placement["evidence"]))
                    self.assertTrue(set(item["id"] for item in placement["evidence"]) -
                                    set(item["id"] for item in original["evidence"]),
                                    "The revised relationship needs fresh observed evidence")
                    self.assertFalse(next(item["enabled"] for item in revised["actions"]
                                          if item["kind"] == "adoptPlacement"))
                    self.assertFalse(next(item["enabled"] for item in revised["actions"]
                                          if item["kind"] == "startPreparation"))
                    inspect_placement(page, placement["proposalId"])
                    self.assertNotEqual(current_state(page)["inspection"]["structureUrl"], original_url)
                    panel = page.locator(".evidence-panel")
                    expect(panel).to_contain_text("Unsupported", ignore_case=True)
                    expect(panel).to_contain_text(placement["reason"])
                    expect(panel).to_contain_text("physical", ignore_case=True)
                    expect(page.get_by_role("button", name="Reassess corrected placement"))\
                        .to_be_visible()
                    capture_review(self, page, "unsupported", 1672, 941)
                    capture_review(self, page, "unsupported", 820, 760)
                finally:
                    browser.close()

    def test_real_candidate_without_placement_policy_is_inspectable_not_established(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            no_policy = directory / "without-placement-policy.json"
            catalogue_without_placement_policy(POLICY, no_policy)
            with running_host(directory / "workspace", no_policy) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    protein = prepare_exact_protein(page)["protein"]["subjectId"]
                    membrane = choose_exact_membrane(page)["membrane"]["modelId"]
                    state = request_placement(page)
                    placement = state["placement"]
                    self.assertEqual(placement["status"], "notEstablished", placement["reason"])
                    self.assertEqual(placement["preparedProteinId"], protein)
                    self.assertEqual(placement["membraneModelId"], membrane)
                    self.assertRegex(placement["reason"], r"(?i)policy|witness")
                    self.assertFalse(next(item["enabled"] for item in state["actions"]
                                          if item["kind"] == "adoptPlacement"))
                    self.assertFalse(next(item["enabled"] for item in state["actions"]
                                          if item["kind"] == "startPreparation"))
                    self.assertEqual(state["stages"], [])
                    inspect_placement(page, placement["proposalId"])
                    panel = page.locator(".evidence-panel")
                    expect(panel).to_contain_text("not established", ignore_case=True)
                    expect(panel).to_contain_text(placement["reason"])
                    expect(panel).to_contain_text("PPM")
                    expect(panel).to_contain_text("DOPC")
                    expect(page.get_by_role("button", name="Reassess corrected placement"))\
                        .to_be_visible()
                    capture_review(self, page, "not-established", 1672, 941)
                    capture_review(self, page, "not-established", 820, 760)
                finally:
                    browser.close()

    def test_missing_local_ppm_does_not_claim_a_position_or_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            missing_binary = directory / "unavailable-immers"
            with running_host(directory / "workspace", POLICY, missing_binary) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    prepare_exact_protein(page)
                    choose_exact_membrane(page)
                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Placement").click()
                    page.locator("#topology-kind").select_option("membrane-spanning")
                    page.locator("#ppm-nterminal-side").select_option("in")
                    page.locator("#biological-sidedness").fill(
                        "extracellular-upper; periplasmic-lower")
                    state = current_state(page)
                    self.assertIsNone(state["placement"])
                    self.assertFalse(next(item["enabled"] for item in state["actions"]
                                          if item["kind"] == "proposePlacement"))
                    self.assertFalse(next(item["enabled"] for item in state["actions"]
                                          if item["kind"] == "adoptPlacement"))
                    self.assertFalse(next(item["enabled"] for item in state["actions"]
                                          if item["kind"] == "startPreparation"))
                    expect(page.get_by_role("button", name="Obtain and assess position"))\
                        .to_be_disabled()
                    expect(page.locator("#placement-workflow")).to_contain_text("hash-verified PPM")
                    expect(page.locator(".stage-strip")).to_contain_text("None")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
