"""Shared exact-actor setup for the single connected 6QWR/DMPC route.

Build with scripts/build-local.sh, then run:

    PIM_BROWSER_CHROMIUM=/path/to/chrome out/browser-test-python/bin/python \
      -m unittest discover \
      -s tests/ProteinInMembraneSystem.AcceptanceTests/ExportCompletedMinimizedStage \
      -p 'test_*.py' -v

The connected test in ExportCompletedMinimizedStage/test_export.py uses these
helpers. This module intentionally has no separate long-running test.
"""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
from urllib.error import URLError
from urllib.request import urlopen

from playwright.sync_api import expect


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "AssessProteinMembranePlacement"))
from test_browser_route import (  # noqa: E402 — reuse the exact accepted 6QWR actor setup
    await_state, choose_exact_membrane, chromium, current_state,
    request_placement,
)

HOST = ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"
LAUNCHER = ROOT / "scripts" / "start-local.sh"
POLICY = Path(os.environ.get("PIM_PREPARATION_POLICY_CATALOGUE", str(
    ROOT / "config" / "policies" / "protein-membrane-current.json"))).resolve()
PPM = Path(os.environ.get("PIM_PPM_EXECUTABLE", str(ROOT / "out" / "ppm2" / "immers"))).resolve()
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "slice4"
SOURCE = ROOT / "config" / "policies" / "source-assets" / "6QWR.pdb"


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
    page.get_by_role("button", name="Upload source").click()
    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=180000)
    if page.get_by_role("button", name="Show inputs").count():
        page.get_by_role("button", name="Show inputs").click()
    page.locator("#model-index").select_option("0")
    page.locator("#assembly-choice").select_option("deposited")
    page.get_by_label("Chain A").check()
    page.get_by_role("button", name="Assess selected protein").click()
    ready = await_state(page, lambda item: (item.get("preparationPlan") or {}).get("standing") ==
                        "ready", "jointly checked current protein plan", timeout=180)
    assert ready["preparationPlan"]["authorizationAvailable"]
    assert ready["preparationPlan"]["stateChoiceCount"] > 0
    assert ready["preparationPlan"]["heavyAtomCount"] == 1
    page.get_by_role("button", name="Prepare with recommendations").click()
    prepared = await_state(page, lambda item: item["protein"] is not None and
                           item["protein"]["status"] == "assessed" and
                           (item.get("proteinTask") or {}).get("standing") == "assessed",
                           "assessed exact prepared protein", timeout=240)
    assert prepared["preparationPlan"]["planSha256"] == ready["preparationPlan"]["planSha256"]
    return prepared


def prepare_adopted_alkl_dmpc(page):
    prepared_id = prepare_exact_protein(page)["protein"]["subjectId"]
    if page.get_by_role("button", name="Show inputs").count():
        page.get_by_role("button", name="Show inputs").click()
    membrane_id = choose_exact_membrane(page, "DMPC")["membrane"]["modelId"]
    placement_state = request_placement(page)
    placement = placement_state["placement"]
    assert placement["status"] == "supported", placement["reason"]
    assert placement["preparedProteinId"] == prepared_id
    assert placement["membraneModelId"] == membrane_id
    assert placement["policyId"]
    page.get_by_role("button", name="Use this position").click()
    state = await_state(page, lambda item: item["study"] is not None and
                        item["study"]["adoptedPlacementProposalId"] == placement["proposalId"],
                        "exact supported placement adopted")
    return state



def start_named_route(page, route_label: str, policy_id: str):
    """Authorize one explicitly named complete route from the actor's current account."""
    before = current_state(page)
    route = next(item for item in before["constructionRoutes"]
                 if item["policyId"] == policy_id)
    assert route["available"], route["reason"]
    assert route["label"] == route_label
    continue_to_build = page.get_by_role("button", name="Continue to system preparation")
    if continue_to_build.is_visible():
        continue_to_build.click()
    if page.locator(".workspace.workflow-closed").count():
        page.get_by_role("button", name="Show workflow").click()
    button = page.locator(".construction-route").filter(has_text=route_label).get_by_role(
        "button", name="Build and minimize")
    expect(button).to_be_enabled()
    button.click()
    admitted = await_state(page, lambda state: state["attempt"] is not None and
                           state["attempt"]["attemptId"] and
                           state["attempt"]["policyId"] == policy_id and
                           state["attempt"]["studyRevisionId"] == before["study"]["id"],
                           "one identified Build and minimize authorization", timeout=300)
    assert admitted["stages"] == []
    return admitted


def await_completed_stage(page, attempt_id: str, timeout: int = 2400):
    account = await_state(page, lambda state: any(
        stage["attemptId"] == attempt_id and stage["kind"] == "Minimization" and
        stage["status"] == "completed" for stage in state["stages"]),
        "completed minimized stage for one accepted attempt", timeout=timeout)
    assert account["attempt"]["attemptId"] == attempt_id
    assert account["attempt"]["status"] == "completed"
    stage = next(stage for stage in account["stages"] if stage["attemptId"] == attempt_id and
                 stage["kind"] == "Minimization")
    assert account["attempt"]["currentStageId"] == stage["stageId"]
    return account, stage
