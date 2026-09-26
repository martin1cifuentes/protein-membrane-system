"""Slice 4 official browser route through the published local host and real worker.

Build with scripts/build-local.sh, then run:

    PIM_BROWSER_CHROMIUM=/path/to/chrome out/browser-test-python/bin/python \
      -m unittest discover \
      -s tests/ProteinInMembraneSystem.AcceptanceTests/ConstructAndMinimizeExplicitSystem \
      -p 'test_*.py' -v

The real positive needs config/policies/protein-membrane-slice4.json and its
identified OpenMM assets. Captures are retained under
out/browser-acceptance/slice4.
"""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest
from urllib.error import URLError
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "AssessProteinMembranePlacement"))
from test_browser_route import (  # noqa: E402 — reuse the exact accepted 6QWR actor setup
    await_state, choose_exact_membrane, chromium, current_state,
    inspect_placement, request_placement,
)

HOST = ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"
WORKER = ROOT / "out" / "python" / "bin" / "python"
LAUNCHER = ROOT / "scripts" / "start-local.sh"
POLICY = Path(os.environ.get("PIM_PREPARATION_POLICY_CATALOGUE", str(
    ROOT / "config" / "policies" / "protein-membrane-slice4.json"))).resolve()
PPM = Path(os.environ.get("PIM_PPM_EXECUTABLE", str(ROOT / "out" / "ppm2" / "immers"))).resolve()
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "slice4"
SOURCE = ROOT / "config" / "policies" / "source-assets" / "6QWR.pdb"
NATIVE_PATCH = ROOT / "out" / "python" / "lib" / "python3.11" / "site-packages" / "openmm" / "app" / "data" / "DMPC.pdb"
TERMINAL_ATTEMPT_STATUSES = {"failed", "resourceRefused", "unobserved", "stopped"}


@contextmanager
def running_host(workspace: Path, catalogue: Path):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    environment = dict(os.environ, PIM_PORT=str(port), PIM_WORKSPACE_ROOT=str(workspace),
                       PIM_POLICY_CATALOGUE=str(catalogue), PIM_PPM_EXECUTABLE=str(PPM))
    log_path = workspace.parent / "preparation-host.log"
    with log_path.open("wb") as log:
        process = subprocess.Popen([str(LAUNCHER)], cwd=ROOT, env=environment,
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
            ARTIFACTS.mkdir(parents=True, exist_ok=True)
            diagnostic = ARTIFACTS / "last-failed-route"
            diagnostic.mkdir(parents=True, exist_ok=True)
            shutil.copy2(log_path, diagnostic / "preparation-host.log")
            try:
                with urlopen(base + "/api/state", timeout=5) as response:
                    account = json.load(response)
                (diagnostic / "account.json").write_text(json.dumps(account, indent=2))
                attempt = account.get("attempt")
                print("\nFinal attempt account: " + json.dumps(None if attempt is None else {
                    key: attempt.get(key) for key in
                    ("attemptId", "status", "stageKind", "message", "currentStageId")
                }))
            except Exception as failure:
                print("\nFinal attempt account unavailable: " + str(failure))
            retained_names = {
                "request.json", "outcome.json", "worker.stderr.log", "worker.stdout.jsonl",
                "constructed-topology.cif", "constructed-topology.json",
                "constructed-system.xml", "constructed-state.xml",
                "constructed-correspondence.json", "minimized-coordinates.cif",
                "minimized-state.xml",
            }
            captured = []
            for source in workspace.rglob("*"):
                if source.is_file() and source.name in retained_names:
                    relative = source.relative_to(workspace)
                    target = diagnostic / relative
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
                    captured.append({"path": str(relative),
                                     "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
            (diagnostic / "captured-files.json").write_text(json.dumps(
                {"workspace": str(workspace), "captured": sorted(captured, key=lambda item: item["path"])},
                indent=2) + "\n")
            print("\nHost output tail:\n" + log_path.read_text()[-4000:])
            raise
        finally:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)


def prepare_exact_protein(page):
    page.locator("#source-upload").set_input_files(str(SOURCE))
    page.locator("#upload-provenance").select_option("experimental")
    page.get_by_role("button", name="Inspect uploaded source").click()
    expect(page.locator(".source-context")).to_contain_text("Selected structural source", timeout=180000)
    page.locator("#model-index").select_option("0")
    page.get_by_label("A · copy A").check()
    page.get_by_role("button", name="Assess selected protein").click()
    proposed = await_state(page, lambda item: item["protein"] is not None and
                           len(item["protein"]["changes"]) > 0,
                           "source-backed protein change proposals", timeout=180)
    for residue_number, variant in ((108, "HID"), (211, "OXT")):
        change = next((item for item in proposed["protein"]["changes"]
                       if item["residue"]["residue"] == residue_number and
                       item["proposedChange"].strip() == variant), None)
        if change is None:
            raise AssertionError(f"Missing source-backed {variant} proposal at residue {residue_number}")
        if page.locator(".workspace.workflow-closed").count():
            page.get_by_role("button", name="Show workflow").click()
        card = page.locator(".change-card").filter(has_text=change["id"])
        inspect = card.get_by_role("button", name="Inspect change and evidence")
        expect(inspect).to_be_enabled(timeout=120000)
        inspect.click()
        decision = page.locator(".proposal-decision")
        expect(page.locator(".proposal-account")).to_contain_text(change["id"])
        rationale = decision.locator(f"#decision-{change['id']}")
        if rationale.count() > 0:
            rationale.fill(f"Accept identified {variant} for first-model 6QWR AlkL at pH 7.")
        approve = decision.get_by_role("button", name="Approve change")
        expect(approve).to_be_enabled(timeout=120000)
        approve.click()
        await_state(page, lambda item: any(action["kind"] == "approvePreparationChange" and
                    action["subjectId"] == change["id"] and not action["enabled"] for action in item["actions"]),
                    f"approved {variant} at residue {residue_number}", timeout=180)
    return await_state(page, lambda item: item["protein"] is not None and
                       item["protein"]["status"] == "assessed",
                       "assessed exact prepared protein", timeout=180)


def prepare_adopted_alkl_dmpc(page):
    prepared_id = prepare_exact_protein(page)["protein"]["subjectId"]
    if page.locator(".workspace.workflow-closed").count():
        page.get_by_role("button", name="Show workflow").click()
    membrane_id = choose_exact_membrane(page, "DMPC")["membrane"]["modelId"]
    placement_state = request_placement(page)
    placement = placement_state["placement"]
    assert placement["status"] == "supported", placement["reason"]
    assert placement["preparedProteinId"] == prepared_id
    assert placement["membraneModelId"] == membrane_id
    assert placement["policyId"] == "alkl-6qwr-dmpc-spanning-placement"
    inspect_placement(page, placement["proposalId"])
    page.get_by_role("button", name="Adopt supported placement").click()
    state = await_state(page, lambda item: item["study"] is not None and
                        item["study"]["adoptedPlacementProposalId"] == placement["proposalId"],
                        "exact supported placement adopted")
    return state


def await_reviewable_candidate(page, attempt_id: str | None = None):
    def observed(item):
        attempt = item["attempt"]
        if attempt is None or (attempt_id is not None and attempt["attemptId"] != attempt_id):
            return False
        return attempt["status"] in TERMINAL_ATTEMPT_STATUSES or (
            attempt["status"] == "readyForMinimization" and
            attempt["derivation"] is not None and attempt["constructed"] is not None and
            not item["stages"])

    # The provider's declared construction bound is 900 s; allow host validation
    # to report the terminal outcome instead of ending a valid in-flight request.
    account = await_state(page, observed, "reviewable or terminal native construction", timeout=1000)
    attempt = account["attempt"]
    if attempt["status"] != "readyForMinimization":
        raise AssertionError("Native construction ended without a reviewable candidate: " +
                             json.dumps(attempt))
    return account


def observe_live_sse_accounts(page):
    page.evaluate("""() => {
      window.__slice4SseAccounts = [];
      window.__slice4Sse = new EventSource('/api/events');
      window.__slice4Sse.onmessage = async () => {
        try {
          const state = await (await fetch('/api/state', {cache: 'no-store'})).json();
          window.__slice4SseAccounts.push({revision: state.revision,
            attemptId: state.attempt?.attemptId, status: state.attempt?.status,
            stageKind: state.attempt?.stageKind, currentStageId: state.attempt?.currentStageId});
        } catch { /* A later live notice may still be observed. */ }
      };
    }""")


def require_live_sse_account(page, attempt_id: str, status: str, stage_kind: str | None = None):
    deadline = time.monotonic() + 30
    observed = None
    while time.monotonic() < deadline:
        observed = page.evaluate("""() => ({readyState: window.__slice4Sse?.readyState,
          accounts: window.__slice4SseAccounts})""")
        if any(item["attemptId"] == attempt_id and item["status"] == status and
               (stage_kind is None or item["stageKind"] == stage_kind)
               for item in observed["accounts"] or []):
            return
        time.sleep(0.15)
    raise AssertionError(f"No correlated live SSE account for {status}/{stage_kind}: " +
                         json.dumps(observed))


def inspect_saved_final_state(test: unittest.TestCase, workspace: Path,
                              attempt_id: str, constructed_account: dict,
                              completed_account: dict, planned: dict,
                              prepared_protein_atom_count: int, final_force: float,
                              final_coordinate_bytes: bytes, observed_local_state: dict):
    atom_count = constructed_account["atomCount"]
    constructed_subject_id = constructed_account["subjectId"]
    attempt_dir = workspace / "attempts" / attempt_id
    filenames = {
        "topologyCif": "constructed-topology.cif",
        "topologyJson": "constructed-topology.json",
        "systemXml": "constructed-system.xml",
        "stateXml": "constructed-state.xml",
        "correspondenceJson": "constructed-correspondence.json",
        "minimizedCif": "minimized-coordinates.cif",
        "minimizedStateXml": "minimized-state.xml",
    }
    paths = {role: attempt_dir / name for role, name in filenames.items()}
    for role, path in paths.items():
        test.assertTrue(path.is_file(), f"Missing saved {role} artifact")
    digests = {role: hashlib.sha256(path.read_bytes()).hexdigest()
               for role, path in paths.items()}
    test.assertEqual(digests["topologyCif"], constructed_subject_id,
                     "The host candidate identity must bind the independently inspected CIF bytes")
    test.assertEqual(paths["minimizedCif"].read_bytes(), final_coordinate_bytes,
                     "The browser's completed-stage coordinates must be the saved final artifact")
    expected = ARTIFACTS / "final-state-expected.json"
    report_path = ARTIFACTS / "final-state-report.json"
    expected.write_text(json.dumps({"atomCount": atom_count, "sha256ByRole": digests}, indent=2))
    inspector = (ROOT / "tests" / "ProteinInMembraneSystem.IntegrationTests" /
                 "ExplicitPreparation" / "inspect_final_state.py")
    completed = subprocess.run([str(WORKER), str(inspector), "--constructed-dir", str(attempt_dir),
                                "--expected", str(expected), "--output", str(report_path)],
                               capture_output=True, text=True, timeout=300)
    test.assertEqual(completed.returncode, 0, completed.stderr[-4000:])
    report = json.loads(report_path.read_text())
    test.assertEqual(report["status"], "observed")
    test.assertEqual(report["atomCount"], atom_count)
    population = report["population"]
    test.assertEqual(population["proteinAtomCount"], prepared_protein_atom_count,
                     "Saved bonded topology must retain the assessed protein atom population")
    observed_lipids = {(item["physicalSide"], item["speciesId"], item["count"])
                       for item in population["lipidCounts"]}
    for label, account in (("reviewed candidate", constructed_account),
                           ("completed stage", completed_account)):
        test.assertEqual(account["subjectId"], constructed_subject_id)
        test.assertEqual(observed_lipids,
                         {(item["physicalSide"], item["speciesId"], item["count"])
                          for item in account["achievedComposition"]},
                         f"Saved bonded lipid populations differ from the {label} account")
        for saved_name, account_name in (("waterCount", "waterCount"),
                                         ("sodiumCount", "sodiumCount"),
                                         ("chlorideCount", "chlorideCount")):
            test.assertEqual(population[saved_name], account[account_name],
                             f"Saved bonded {saved_name} differs from the {label} account")
    test.assertEqual(observed_lipids,
                     {(item["physicalSide"], item["speciesId"], item["count"])
                      for item in planned["lipidCounts"]})
    for name in ("waterCount", "sodiumCount", "chlorideCount"):
        test.assertEqual(population[name], planned[name])
    test.assertAlmostEqual(report["finalForce"]["tangentPerParticleRmsKjMolNm"],
                           final_force, delta=0.001,
                           msg="Independent final-State force must agree with the host's stage observation")
    test.assertLessEqual(report["finalForce"]["maximumRelativeConstraintError"], 1e-5)
    test.assertEqual(set(report["stages"]["final"]["heavyPeriodicContacts"]),
                     {"water|water", "water|ion", "ion|ion", "protein|lipid"})
    initial_contacts = report["stages"]["initial"]["heavyPeriodicContacts"]
    final_contacts = report["stages"]["final"]["heavyPeriodicContacts"]
    test.assertGreater(initial_contacts["water|water"]["heavyPairsBelow2_2Angstrom"], 0,
                       "The representative starting State must expose the correctable solvent seam")
    for role_pair in ("water|water", "water|ion", "ion|ion"):
        test.assertEqual(final_contacts[role_pair]["heavyPairsBelow2_2Angstrom"], 0,
                         f"The exact first-class minimized State retains a {role_pair} severe seam")
    initial_all = report["stages"]["initial"]["allAtomPeriodicContacts"]
    final_all = report["stages"]["final"]["allAtomPeriodicContacts"]
    test.assertGreater(initial_all["wholeSystem"]["allAtomPairsBelow1_5Angstrom"], 0)
    for role_pair in ("water|water", "water|ion", "ion|ion"):
        test.assertEqual(final_all[role_pair]["allAtomPairsBelow1_5Angstrom"], 0,
                         f"The minimized State retains an all-atom {role_pair} seam below 1.5 Å")
    test.assertGreater(final_all["wholeSystem"]["nearestAllAtomAngstrom"],
                       initial_all["wholeSystem"]["nearestAllAtomAngstrom"])
    measured_nearest = next(item["value"] for item in observed_local_state["measurements"]
                            if item["name"] == "minimumIntermolecularDistanceAngstrom" and
                            item["scope"] == "wholeSystem")
    test.assertAlmostEqual(final_all["wholeSystem"]["nearestAllAtomAngstrom"],
                           measured_nearest, delta=0.01,
                           msg="Independent final all-atom minimum must agree with the worker observation")
    test.assertGreater(final_contacts["protein|lipid"]["heavyPairsWithinRadius"], 0)
    test.assertLessEqual(final_contacts["protein|lipid"]["nearestHeavyAngstrom"], 6.0)
    initial_relative = report["stages"]["initial"]["relativeOrganization"]
    final_relative = report["stages"]["final"]["relativeOrganization"]
    test.assertGreater(final_relative["leafletHeadSeparationAngstrom"], 0)
    test.assertTrue(all(math.isfinite(value) for value in final_relative.values()))
    test.assertTrue(any(abs(final_relative[key] - initial_relative[key]) > 1e-5
                        for key in final_relative),
                    "Final protein and leaflet organization must be measured from changed coordinates")
    retained = ARTIFACTS / "final-state-artifacts"
    retained.mkdir(parents=True, exist_ok=True)
    for role, path in paths.items():
        copy = retained / path.name
        shutil.copy2(path, copy)
        test.assertEqual(hashlib.sha256(copy.read_bytes()).hexdigest(), digests[role],
                         f"Retained {role} evidence differs from the inspected stage")


def capture_review(test: unittest.TestCase, page, name: str, width: int, height: int):
    before = current_state(page)
    if name == "candidate":
        test.assertEqual(before["attempt"]["status"], "readyForMinimization")
        test.assertIsNotNone(before["attempt"]["constructed"])
        test.assertEqual(before["stages"], [])
    elif name == "running":
        test.assertEqual(before["attempt"]["status"], "running")
        test.assertEqual(before["attempt"]["stageKind"], "Minimization")
        test.assertIsNotNone(before["attempt"]["constructed"])
        test.assertEqual(before["stages"], [])
    else:
        test.assertTrue(any(stage["status"] == "completed" and stage["kind"] == "Minimization"
                            for stage in before["stages"]))
    page.set_viewport_size({"width": width, "height": height})
    page.evaluate("() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))")
    # Let the canvas resize observation settle before the actor's Focus action.
    page.wait_for_timeout(150)
    page.locator(".placement-view-tools").get_by_role("button", name="Focus").click()
    page.locator(".evidence-content").evaluate("element => { element.scrollTop = 0; }")
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    test.assertIsNotNone(scene)
    test.assertIsNotNone(evidence)
    test.assertIsNotNone(stages)
    test.assertGreaterEqual(scene["width"], 350 if width <= 820 else width * 0.64)
    test.assertGreaterEqual(evidence["width"], 300 if width <= 820 else width * 0.25)
    test.assertGreaterEqual(scene["height"], height * 0.55)
    test.assertGreater(evidence["x"], scene["x"])
    test.assertAlmostEqual(evidence["y"], scene["y"], delta=3)
    test.assertLessEqual(stages["y"] + stages["height"], height + 2)
    overflowing = page.evaluate("""() => [...document.querySelectorAll(
        '.workspace, .workspace-header, .scene-panel, .evidence-panel, .evidence-content')]
      .filter(element => element.scrollWidth > element.clientWidth + 1)
      .map(element => ({region: element.className, width: element.clientWidth,
                        contentWidth: element.scrollWidth}))""")
    test.assertEqual(overflowing, [])
    expect(page.locator(".execution-scene-label")).to_be_visible()
    expect(page.locator(".execution-decision")).to_be_visible()
    page.screenshot(path=str(ARTIFACTS / f"{name}-{width}.png"),
                    full_page=True, animations="disabled")
    after = current_state(page)
    if name == "candidate":
        test.assertEqual(after["attempt"]["status"], "readyForMinimization")
        test.assertEqual(after["stages"], [])
    elif name == "running":
        test.assertEqual(after["attempt"]["status"], "running",
                         "A running-state capture must coincide with real ongoing work")
        test.assertEqual(after["stages"], [])
    if width <= 820:
        if name in ("running", "candidate"):
            details = page.locator(".execution-actual-values")
            if details.get_attribute("open") is None:
                details.locator("summary").first.click()
            row = details.locator("dt").filter(has_text="Actual cell")
        else:
            row = page.locator(".execution-findings-account details p").filter(
                has_text="Fresh protein geometry")
        row.scroll_into_view_if_needed()
        expect(row).to_be_visible()
        decision = page.locator(".execution-decision")
        expect(decision).to_be_visible()
        row_box = row.bounding_box()
        decision_box = decision.bounding_box()
        test.assertIsNotNone(row_box)
        test.assertIsNotNone(decision_box)
        test.assertLessEqual(row_box["y"] + row_box["height"], decision_box["y"] + 1,
                             "Exact evidence must be scrollable above the standing/action footer")
        page.locator(".evidence-content").evaluate("element => { element.scrollTop = 0; }")


def catalogue_without_preparation_policy(destination: Path):
    catalogue = json.loads(POLICY.read_text())
    source_base = POLICY.parent

    def absolute_paths(value):
        if isinstance(value, list):
            return [absolute_paths(item) for item in value]
        if isinstance(value, dict):
            return {key: str((source_base / item).resolve())
                    if key in {"path", "templatePath", "coordinateTemplatePath", "coordinatePath",
                               "paePath", "mappingPath", "ppmResidueLibraryPath", "nativePatchPath"}
                    and isinstance(item, str) and item and not Path(item).is_absolute()
                    else absolute_paths(item)
                    for key, item in value.items()}
        return value

    catalogue = absolute_paths(catalogue)
    catalogue["preparationPolicies"] = []
    catalogue["version"] += "-test-no-preparation-policy"
    destination.write_text(json.dumps(catalogue, indent=2) + "\n")


def atom_coordinates(mmcif: bytes):
    """Read the atom-site Cartesian loop independently of the product account."""
    lines = mmcif.decode("utf-8").splitlines()
    headers = []
    for index, line in enumerate(lines):
        if line.startswith("_atom_site."):
            headers.append(line.strip())
            continue
        if headers and line and not line.startswith("_atom_site."):
            required = ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
            if any(name not in headers for name in required):
                raise AssertionError("The inspected mmCIF has no complete atom-site coordinates")
            columns = [headers.index(name) for name in required]
            positions = []
            for atom_line in lines[index:]:
                if atom_line.startswith("#") or atom_line.startswith("loop_"):
                    break
                if not atom_line.startswith(("ATOM", "HETATM")):
                    continue
                fields = atom_line.split()
                positions.append(tuple(float(fields[column]) for column in columns))
            if not positions:
                raise AssertionError("The inspected mmCIF has no atom-site rows")
            return positions
    raise AssertionError("The inspected mmCIF has no atom-site loop")


class BrowserConstructionAndMinimizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, WORKER, LAUNCHER, PPM, NATIVE_PATCH, POLICY,
                     ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Slice 4 official route prerequisite unavailable: {path}")
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_real_constructed_system_running_and_completed_minimized_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    adopted = prepare_adopted_alkl_dmpc(page)
                    self.assertTrue(any(line.startswith("CRYST1") and
                                        line[6:15].strip() == "1.000" for line in SOURCE.read_text().splitlines()),
                                    "This representative source must exercise the placeholder-cell correction")
                    oriented_files = list((directory / "workspace" / "placement").rglob(
                        "oriented-protein.pdb"))
                    self.assertEqual(len(oriented_files), 1)
                    self.assertNotIn("CRYST1", oriented_files[0].read_text(),
                                     "The isolated oriented protein must not carry 6QWR's dummy 1 Å cell")
                    self.assertEqual(adopted["stages"], [])
                    page.get_by_role("button", name="Show workflow").click()
                    observe_live_sse_accounts(page)
                    start = page.get_by_role("button", name="Construct system")
                    expect(start).to_be_enabled()
                    start.click()
                    admitted = await_state(page, lambda item: item["attempt"] is not None and
                                           item["attempt"]["attemptId"] and item["attempt"]["studyRevisionId"] ==
                                           adopted["study"]["id"], "identified accepted attempt", timeout=300)
                    attempt_id = admitted["attempt"]["attemptId"]
                    self.assertEqual(admitted["attempt"]["policyId"],
                                     "alkl-6qwr-dmpc-explicit-minimization")
                    self.assertEqual(admitted["attempt"]["policyVersion"], "1.2.0")
                    self.assertEqual(admitted["stages"], [])
                    candidate = await_reviewable_candidate(page, attempt_id)
                    planned = candidate["attempt"]["derivation"]
                    constructed = candidate["attempt"]["constructed"]
                    require_live_sse_account(page, attempt_id, "running")
                    require_live_sse_account(page, attempt_id, "readyForMinimization")
                    self.assertEqual(planned["attemptId"], attempt_id)
                    self.assertEqual(constructed["attemptId"], attempt_id)
                    self.assertTrue(planned["approximations"])
                    self.assertGreater(constructed["atomCount"], 0)
                    self.assertEqual({item["speciesId"] for item in constructed["achievedComposition"]}, {"DMPC"})
                    self.assertEqual({item["physicalSide"] for item in constructed["achievedComposition"]},
                                     {"upper", "lower"})
                    self.assertEqual(
                        {(item["physicalSide"], item["speciesId"], item["count"])
                         for item in constructed["achievedComposition"]},
                        {(item["physicalSide"], item["speciesId"], item["count"])
                         for item in planned["lipidCounts"]},
                        "Validated achieved leaflets must correspond to the proposed finite counts")
                    self.assertGreater(constructed["waterCount"], 0)
                    self.assertEqual(constructed["waterCount"], planned["waterCount"])
                    self.assertEqual(constructed["sodiumCount"], planned["sodiumCount"])
                    self.assertEqual(constructed["chlorideCount"], planned["chlorideCount"])
                    self.assertEqual(candidate["stages"], [], "A constructed candidate is not a completed stage")
                    self.assertTrue(any(action["kind"] == "continueMinimization" and action["enabled"]
                                        for action in candidate["actions"]),
                                    "The exact current candidate needs an explicit continuation action")
                    def refused_continue(data, expected_reason, stale=False):
                        for _ in range(3):
                            revision = page.request.get(base + "/api/state").json()["revision"]
                            reply = page.request.post(base + "/api/commands",
                                                      data=json.dumps({"kind": "continueMinimization",
                                                                       "data": data,
                                                                       "expectedRevision": revision - int(stale)}),
                                                      headers={"Content-Type": "application/json",
                                                               "Origin": base})
                            self.assertEqual(reply.status, 422, reply.text())
                            reason = reply.json()["reason"]
                            if expected_reason in reason:
                                break
                            if stale or "workspace changed" not in reason:
                                self.fail(f"Unexpected continuation refusal: {reason}")
                        else:
                            self.fail("The workspace revision changed during every exact-candidate refusal")
                        unchanged = page.request.get(base + "/api/state").json()
                        self.assertEqual(unchanged["attempt"]["status"], "readyForMinimization")
                        self.assertEqual(unchanged["stages"], [])

                    refused_continue({"attemptId": attempt_id,
                                      "constructedSubjectId": "wrong-candidate"},
                                     "exact current constructed candidate")
                    refused_continue({"attemptId": "wrong-attempt",
                                      "constructedSubjectId": constructed["subjectId"]},
                                     "exact current constructed candidate")
                    refused_continue({"attemptId": attempt_id,
                                      "constructedSubjectId": constructed["subjectId"]},
                                     "workspace changed", stale=True)
                    expect(page.locator(".execution-identity-account")).to_contain_text(
                        "Constructed system review", timeout=120000)
                    capture_review(self, page, "candidate", 1672, 941)
                    capture_review(self, page, "candidate", 820, 760)
                    page.set_viewport_size({"width": 1672, "height": 941})
                    if page.locator(".workspace.workflow-closed").count():
                        page.get_by_role("button", name="Show workflow").click()
                    continue_button = page.get_by_role("button", name="Continue minimization").first
                    expect(continue_button).to_be_enabled(timeout=120000)
                    continue_button.click()
                    running = await_state(page, lambda item: item["attempt"] is not None and
                                          item["attempt"]["attemptId"] == attempt_id and
                                          item["attempt"]["stageKind"] == "Minimization" and
                                          item["attempt"]["status"] == "running" and
                                          item["attempt"]["constructed"] is not None and not item["stages"],
                                          "explicitly continued minimization of reviewed candidate", timeout=120)
                    self.assertEqual(running["attempt"]["constructed"]["subjectId"], constructed["subjectId"])
                    require_live_sse_account(page, attempt_id, "running", "Minimization")
                    expect(page.locator(".execution-identity-account")).to_contain_text(
                        "Current minimization", timeout=120000)
                    selected_inspection = current_state(page)["inspection"]
                    if selected_inspection is None or selected_inspection["subjectId"] != constructed["subjectId"]:
                        inspect_constructed = page.get_by_role("button", name="Inspect verified constructed system")
                        expect(inspect_constructed).to_be_enabled(timeout=120000)
                        inspect_constructed.click()
                    inspected = await_state(page, lambda item: item["inspection"] is not None and
                                            item["inspection"]["subjectId"] == constructed["subjectId"] and
                                            item["inspection"]["representationKind"] == "constructedSystem",
                                            "verified full-system inspection", timeout=120)
                    self.assertTrue(inspected["inspection"]["structureUrl"])
                    constructed_response = page.request.get(base + inspected["inspection"]["structureUrl"])
                    self.assertTrue(constructed_response.ok)
                    constructed_coordinates = atom_coordinates(constructed_response.body())
                    self.assertEqual(len(constructed_coordinates), constructed["atomCount"])
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
                    expect(page.locator(".scene-error")).to_have_count(0)
                    expect(page.locator(".execution-scene-label")).to_contain_text("Current constructed system")
                    expect(page.locator(".stage-strip")).to_contain_text("None")
                    expect(page.locator(".stage-strip")).to_contain_text("Not applicable")
                    capture_review(self, page, "running", 1672, 941)
                    capture_review(self, page, "running", 820, 760)
                    completed = await_state(page, lambda item: any(stage["attemptId"] == attempt_id and
                                            stage["kind"] == "Minimization" and stage["status"] == "completed"
                                            for stage in item["stages"]), "observed completed minimized stage", timeout=1200)
                    stage = next(item for item in completed["stages"] if item["attemptId"] == attempt_id and
                                 item["kind"] == "Minimization")
                    require_live_sse_account(page, attempt_id, "completed", "Minimization")
                    self.assertEqual(stage["studyRevisionId"], adopted["study"]["id"])
                    self.assertIsNotNone(stage["constructed"])
                    self.assertEqual(stage["constructed"]["subjectId"], constructed["subjectId"])
                    self.assertEqual(stage["constructed"]["cellAngstrom"], constructed["cellAngstrom"])
                    self.assertIsNotNone(stage["observation"])
                    self.assertEqual(stage["observation"]["termination"], "converged")
                    final_force = next(value for value in stage["observation"]["measurements"]
                                       if value["name"] == "finalRmsForce" and
                                       value["scope"] == "final unrestrained constraint tangent; per particle")
                    self.assertEqual(final_force["unit"], "kJ mol^-1 nm^-1")
                    self.assertLessEqual(final_force["value"], 10.0)
                    self.assertIsNotNone(stage["assessment"])
                    self.assertTrue(stage["assessment"]["currentlyApplicable"])
                    actual_qualification = stage["assessment"]["qualification"].casefold()
                    self.assertIn(actual_qualification, ("indeterminate", "notqualified"),
                                  "Without positive criteria, only an attributable defect may establish NotQualified")
                    self.assertTrue(stage["assessment"]["reason"])
                    if actual_qualification == "notqualified":
                        reason = stage["assessment"]["reason"].lower()
                        if "material finding" in reason:
                            self.assertTrue(any(item["material"] and
                                                item["disposition"].lower() == "disqualifies"
                                                for item in stage["assessment"]["findings"]))
                        elif "protein intramolecular geometry" in reason:
                            self.assertEqual(stage["observation"]["proteinGeometry"]["standing"], "Observed")
                        elif "molecular contact" in reason:
                            self.assertEqual(stage["observation"]["localState"]["standing"], "Observed")
                        else:
                            self.fail("NotQualified lacks an inspectable measured-defect reason: " + reason)
                    self.assertIn("localState", stage["observation"])
                    self.assertIn("proteinGeometry", stage["observation"])
                    self.assertEqual(stage["observation"]["localState"]["standing"], "Observed")
                    self.assertEqual(stage["observation"]["proteinGeometry"]["standing"], "Observed")
                    page.locator(".stage-card").filter(has_text=stage["stageId"]).click()
                    selected = await_state(page, lambda item: item["inspection"] is not None and
                                           item["inspection"]["subjectId"] == stage["stageId"] and
                                           item["inspection"]["representationKind"] == "completedStage",
                                           "completed stage selected for exact inspection")
                    self.assertTrue(selected["inspection"]["structureUrl"])
                    final_response = page.request.get(base + selected["inspection"]["structureUrl"])
                    self.assertTrue(final_response.ok)
                    final_coordinates = atom_coordinates(final_response.body())
                    self.assertEqual(len(final_coordinates), len(constructed_coordinates))
                    self.assertTrue(any(any(abs(before - after) > 1e-4 for before, after in zip(start, end))
                                        for start, end in zip(constructed_coordinates, final_coordinates)),
                                    "The selected final coordinates must reflect actual minimization")
                    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
                    expect(page.locator(".scene-error")).to_have_count(0)
                    expect(page.locator(".execution-assessment-account")).to_contain_text("Scientific assessment")
                    expect(page.locator(".execution-assessment-account")).to_contain_text("Unrestrained convergence")
                    displayed_force = page.locator(".execution-assessment-account dt").filter(
                        has_text="Final RMS force").locator("xpath=following-sibling::dd[1]").inner_text()
                    self.assertAlmostEqual(float(displayed_force.split()[0].replace(",", "")),
                                           final_force["value"], delta=0.001,
                                           msg="The stage review must display the physical tangent force, not raw constraint reaction")
                    expect(page.locator(".execution-assessment-account")).to_contain_text(
                        stage["assessment"]["reason"])
                    details = page.locator(".execution-findings-account details")
                    details.locator("summary").click()
                    expect(details).to_contain_text("Fresh local-state observation")
                    expect(details).to_contain_text("Fresh protein geometry")
                    self.assertTrue(any(action["kind"] == "exportStage" and
                                        action["subjectId"] == stage["stageId"] and
                                        not action["enabled"] for action in completed["actions"]),
                                    "Unverified Slice 6 export must not be offered for this completed stage")
                    expect(page.get_by_role("button", name="Export with status")).to_have_count(0)
                    capture_review(self, page, "minimized", 1672, 941)
                    capture_review(self, page, "minimized", 820, 760)
                    inspect_saved_final_state(self, directory / "workspace", attempt_id,
                                              constructed, stage["constructed"], planned,
                                              adopted["protein"]["atomCount"],
                                              final_force["value"], final_response.body(),
                                              stage["observation"]["localState"])
                finally:
                    browser.close()

    def test_missing_preparation_policy_refuses_before_an_attempt(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            no_policy = directory / "without-preparation-policy.json"
            catalogue_without_preparation_policy(no_policy)
            with running_host(directory / "workspace", no_policy) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    prepare_adopted_alkl_dmpc(page)
                    page.get_by_role("button", name="Show workflow").click()
                    expect(page.get_by_role("button", name="Construct system")).to_be_disabled()
                    account = current_state(page)
                    self.assertIsNone(account["attempt"])
                    self.assertEqual(account["stages"], [])
                    self.assertTrue(next(item["reason"] for item in account["actions"]
                                         if item["kind"] == "startPreparation"))
                    expect(page.locator("#researcher-workflow")).to_contain_text("policy", ignore_case=True)
                finally:
                    browser.close()

    def test_reviewed_candidate_can_be_declined_without_a_minimized_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    prepare_adopted_alkl_dmpc(page)
                    page.get_by_role("button", name="Show workflow").click()
                    page.get_by_role("button", name="Construct system").click()
                    ready = await_reviewable_candidate(page)
                    attempt_id = ready["attempt"]["attemptId"]
                    candidate_id = ready["attempt"]["constructed"]["subjectId"]
                    self.assertEqual(ready["stages"], [])
                    page.close()
                    with urlopen(base + "/api/state", timeout=5) as response:
                        surviving = json.load(response)
                    self.assertEqual(surviving["attempt"]["status"], "readyForMinimization")
                    self.assertEqual(surviving["attempt"]["constructed"]["subjectId"], candidate_id)
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    interrupted_events = []
                    def lose_sse(route):
                        interrupted_events.append(route.request.url)
                        route.abort("failed")
                    page.route("**/api/events", lose_sse)
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.get_by_role("alert").filter(
                        has_text="The live connection is interrupted")).to_be_visible(timeout=10000)
                    self.assertTrue(interrupted_events)
                    reopened = current_state(page)
                    self.assertEqual(reopened["attempt"]["status"], "readyForMinimization")
                    self.assertEqual(reopened["attempt"]["constructed"]["subjectId"], candidate_id)
                    expect(page.get_by_role("button", name="Decline candidate").first).to_be_enabled()
                    page.get_by_role("button", name="Decline candidate").first.click()
                    stopped = await_state(page, lambda item: item["attempt"] is not None and
                                          item["attempt"]["attemptId"] == attempt_id and
                                          item["attempt"]["status"] == "stopped",
                                          "declined exact candidate", timeout=120)
                    self.assertEqual(stopped["stages"], [])
                    self.assertEqual(stopped["attempt"]["constructed"]["subjectId"], candidate_id)
                    self.assertFalse(any(action["kind"] == "continueMinimization" and action["enabled"]
                                         for action in stopped["actions"]))
                    page.reload(wait_until="domcontentloaded")
                    persisted = current_state(page)
                    self.assertEqual(persisted["attempt"]["attemptId"], attempt_id)
                    self.assertEqual(persisted["attempt"]["status"], "stopped")
                    self.assertEqual(persisted["stages"], [])
                finally:
                    browser.close()

    def test_exact_running_attempt_can_be_stopped_without_a_completed_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    prepare_adopted_alkl_dmpc(page)
                    page.get_by_role("button", name="Show workflow").click()
                    observe_live_sse_accounts(page)
                    page.get_by_role("button", name="Construct system").click()
                    running = await_state(page, lambda item: item["attempt"] is not None and
                                          item["attempt"]["attemptId"] and
                                          item["attempt"]["status"] == "running" and
                                          next((action["enabled"] for action in item["actions"]
                                                if action["kind"] == "stopAttempt"), False),
                                          "stoppable accepted attempt", timeout=300)
                    attempt_id = running["attempt"]["attemptId"]
                    require_live_sse_account(page, attempt_id, "running")
                    self.assertEqual(running["stages"], [])
                    # A view and its SSE subscription may disappear while the
                    # application-owned attempt continues. Reconnect without
                    # SSE and inspect the same host account before requesting
                    # the explicit stop below.
                    page.close()
                    with urlopen(base + "/api/state", timeout=5) as response:
                        surviving = json.load(response)
                    self.assertEqual(surviving["attempt"]["attemptId"], attempt_id)
                    self.assertEqual(surviving["attempt"]["status"], "running")
                    self.assertEqual(surviving["stages"], [])
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    interrupted_events = []
                    def lose_sse(route):
                        interrupted_events.append(route.request.url)
                        route.abort("failed")
                    page.route("**/api/events", lose_sse)
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.get_by_role("alert").filter(
                        has_text="The live connection is interrupted")).to_be_visible(timeout=10000)
                    self.assertTrue(interrupted_events, "The SSE route was not actually interrupted")
                    reopened = await_state(page, lambda item: item["attempt"] is not None and
                                           item["attempt"]["attemptId"] == attempt_id,
                                           "same admitted attempt after SSE loss", timeout=30)
                    self.assertEqual(reopened["attempt"]["status"], "running")
                    self.assertEqual(reopened["stages"], [])
                    expect(page.get_by_role("button", name="Stop unfinished work").last).to_be_enabled()
                    page.get_by_role("button", name="Stop unfinished work").last.click()
                    stopped = await_state(page, lambda item: item["attempt"] is not None and
                                          item["attempt"]["attemptId"] == attempt_id and
                                          item["attempt"]["status"] == "stopped",
                                          "stopped exact attempt", timeout=180)
                    self.assertEqual(stopped["stages"], [])
                    page.reload(wait_until="domcontentloaded")
                    expect(page.locator(".attempt-account")).to_contain_text("stopped")
                    expect(page.locator(".stage-strip")).to_contain_text("None")
                    self.assertFalse(any(stage["attemptId"] == attempt_id for stage in stopped["stages"]))
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
