"""Slice 1 browser acceptance through the published loopback Host and real worker.

Run after the local build with the acceptance-test environment:
    out/browser-test-python/bin/python -m unittest discover -s tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein -p 'test_*.py'

Set PIM_BROWSER_CHROMIUM to override Playwright's installed Chromium binary.
Screenshots are retained under out/browser-acceptance for visual review.
Set PIM_BROWSER_REMOTE=1 to include a live exact RCSB source review.
"""

from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
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
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.IntegrationTests" / "ProteinPreparation"))
from test_worker_exchange import two_alanines  # noqa: E402

HOST = ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"
WORKER = ROOT / "out" / "python" / "bin" / "python"
POLICY = ROOT / "config" / "policies" / "protein-slice1.json"
ARTIFACTS = ROOT / "out" / "browser-acceptance"


@contextmanager
def running_host(workspace: Path):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    environment = dict(os.environ,
                       ASPNETCORE_URLS=base,
                       PIM_WORKER_PYTHON=str(WORKER),
                       PIM_OWNER_SOURCE_ROOT=str(ROOT / "src" / "ProteinInMembrane.Host"),
                       PIM_WORKSPACE_ROOT=str(workspace),
                       PIM_POLICY_CATALOGUE=str(POLICY))
    environment.pop("PIM_PPM_EXECUTABLE", None)
    environment.pop("PIM_PACKMOL_EXECUTABLE", None)
    log_path = workspace.parent / "host.log"
    with log_path.open("wb") as log:
        process = subprocess.Popen(["dotnet", str(HOST)], cwd=ROOT,
                                   env=environment, stdout=log, stderr=subprocess.STDOUT)
        try:
            for _ in range(100):
                if process.poll() is not None:
                    raise AssertionError("Host exited during startup: " + log_path.read_text())
                try:
                    with urlopen(base + "/api/state", timeout=2) as response:
                        if response.status == 200:
                            break
                except URLError:
                    time.sleep(0.1)
            else:
                raise AssertionError("Host did not become ready: " + log_path.read_text())
            yield base
        except BaseException:
            print("\nHost output:\n" + log_path.read_text()[-12000:], file=sys.stderr)
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
        executable_path=executable or None,
        headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"],
    )


def upload_and_choose(page, source: Path):
    page.locator("#source-upload").set_input_files(str(source))
    page.locator("#upload-provenance").select_option("experimental")
    page.get_by_role("button", name="Inspect uploaded source").click()
    expect(page.locator(".source-context")).to_contain_text("Selected structural source", timeout=120000)
    expect(page.locator("#model-index")).to_be_visible()
    page.locator("#model-index").select_option("0")
    page.get_by_label("A · copy A").check()
    page.get_by_role("button", name="Assess selected protein").click()


def incomplete_alanine_source() -> str:
    # The second ALA retains its observed backbone and is missing only CB.
    lines = two_alanines(1, "A").splitlines(keepends=True)
    return "".join(line for line in lines if not (
        line.startswith("ATOM") and line[12:16].strip() == "CB" and
        line[22:26].strip() == "2")) + "END\n"


def horizontally_overflowing_regions(page):
    return page.evaluate("""() => [...document.querySelectorAll(
        '.workspace, .workspace-header, .rail, .scene-panel, .evidence-panel, .evidence-content')]
      .filter(element => element.scrollWidth > element.clientWidth + 1)
      .map(element => ({ region: element.className, width: element.clientWidth,
                         contentWidth: element.scrollWidth,
                         overflowingChildren: [...element.querySelectorAll('*')]
                           .filter(child => child.getBoundingClientRect().right >
                             element.getBoundingClientRect().right + 1)
                           .slice(0, 5).map(child => child.tagName + '.' + child.className) }))""")


def wait_for_selected_structure(page):
    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=120000)
    expect(page.locator(".scene-loading")).to_have_count(0, timeout=120000)
    expect(page.locator(".scene-error")).to_have_count(0)
    # Mol* completes the first WebGL frame after its structure-loading promise.
    page.wait_for_timeout(700)


class BrowserRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, WORKER, POLICY, ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise unittest.SkipTest(f"Local build prerequisite is missing: {path}")
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_uploaded_exact_model_reaches_assessed_protein_and_connected_inspection(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1440, "height": 900})
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.locator("#exact-identifier")).to_be_visible()
                    expect(page.get_by_role("button", name="Inspect exact source")).to_be_disabled()
                    upload_and_choose(page, source)
                    expect(page.locator(".protein-standing")).to_contain_text(
                        "Assessed prepared protein", timeout=120000)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["protein"]["status"], "assessed")
                    self.assertEqual(state["study"]["modelIndex"], 0)
                    self.assertEqual(state["study"]["chainIds"], ["A"])
                    self.assertIsNone(state["placement"])
                    self.assertEqual(state["stages"], [])
                    self.assertIn(state["study"]["selectedSourceId"],
                                  page.locator(".source-context").inner_text())
                    page.get_by_role("button", name="Inspect protein", exact=True).click()
                    expect(page.locator(".workspace")).to_have_class("workspace workflow-closed")
                    expect(page.locator(".scene-subtitle")).to_contain_text(
                        state["protein"]["subjectId"])
                    expect(page.locator(".evidence-panel")).to_contain_text(
                        state["protein"]["subjectId"])
                    expect(page.locator(".evidence-panel")).to_contain_text("Protein preparation")
                    expect(page.locator(".stage-strip")).to_contain_text(
                        "No completed stage is currently established")
                    expect(page.locator(".stage-strip")).to_contain_text("Stage assessment")
                    expect(page.locator(".source-account")).to_contain_text(
                        state["study"]["selectedSourceId"])
                    wait_for_selected_structure(page)
                    self.assertEqual(horizontally_overflowing_regions(page), [])
                    page.screenshot(path=str(ARTIFACTS / "assessed-1440.png"), full_page=True,
                                    animations="disabled")
                    page.get_by_role("button", name="Show workflow").click()
                    expect(page.locator("#researcher-workflow")).to_be_visible()
                    expect(page.locator("#researcher-workflow")).to_contain_text("3 · Planar membrane")
                    print("PASS: browser upload, exact model/chain, assessed protein and linked inspection")
                finally:
                    browser.close()

    def test_required_change_review_decline_and_reduced_desktop(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "missing-sidechain-cb.pdb"
            source.write_text(incomplete_alanine_source())
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1440, "height": 900})
                    page.goto(base, wait_until="domcontentloaded")
                    upload_and_choose(page, source)
                    proposal = page.locator(".change-card")
                    expect(proposal).to_have_count(1, timeout=120000)
                    expect(proposal).to_contain_text("Heavy atom · A2 · copy A")
                    before_review = page.evaluate("async () => (await fetch('/api/state')).json()")
                    pending_id = before_review["protein"]["changes"][0]["id"]
                    self.assertFalse(next(item["enabled"] for item in before_review["actions"]
                                          if item["kind"] == "approvePreparationChange"
                                          and item["subjectId"] == pending_id))
                    proposal.get_by_role("button", name="Inspect change and evidence").click()
                    expect(proposal).to_have_class("change-card selected")
                    expect(page.locator(".workspace")).to_have_class("workspace workflow-closed")
                    decision = page.locator(".proposal-decision")
                    approve = decision.get_by_role("button", name="Approve change")
                    expect(page.locator(".evidence-panel")).to_contain_text("Proposed change")
                    expect(page.locator(".evidence-panel")).to_contain_text("Affected region")
                    expect(page.locator(".evidence-panel")).to_contain_text("Source and assembly")
                    expect(page.locator(".evidence-panel")).to_contain_text("Located evidence and findings")
                    expect(approve).to_be_enabled()
                    wait_for_selected_structure(page)
                    proposal_id = page.evaluate("async () => (await fetch('/api/state')).json()")\
                        ["protein"]["changes"][0]["id"]
                    page.locator(".evidence-content").evaluate("element => element.scrollTop = 0")
                    page.screenshot(path=str(ARTIFACTS / "proposal-1440.png"), full_page=True,
                                    animations="disabled")

                    page.set_viewport_size({"width": 820, "height": 760})
                    page.wait_for_function("() => document.querySelector('.evidence-content')?.scrollTop > 100")
                    scene_box = page.locator(".scene-panel").bounding_box()
                    evidence_box = page.locator(".evidence-panel").bounding_box()
                    self.assertIsNotNone(scene_box)
                    self.assertIsNotNone(evidence_box)
                    self.assertGreaterEqual(scene_box["width"], 350)
                    self.assertGreaterEqual(scene_box["height"], 450)
                    self.assertGreater(evidence_box["x"], scene_box["x"])
                    self.assertAlmostEqual(evidence_box["y"], scene_box["y"], delta=2)
                    self.assertEqual(horizontally_overflowing_regions(page), [])
                    expect(decision.get_by_role("button", name="Decline")).to_be_visible()
                    expect(page.locator(".scene-surface")).to_be_visible()
                    expect(page.locator(".compact-model-context")).to_be_visible()
                    self.assertTrue(page.evaluate("""() => {
                        const item = document.querySelector('.annotation-item').getBoundingClientRect();
                        const pane = document.querySelector('.evidence-content').getBoundingClientRect();
                        return item.top < pane.bottom && item.bottom > pane.top;
                    }"""))
                    page.screenshot(path=str(ARTIFACTS / "proposal-820.png"), full_page=True,
                                    animations="disabled")

                    page.get_by_role("button", name="Focus affected region in structure").click()
                    expect(page.locator(".annotation-item.selected")).to_have_count(1)
                    page.get_by_role("button", name="Clear selected part").click()
                    expect(page.locator(".annotation-item.selected")).to_have_count(0)

                    decision.get_by_role("button", name="Decline").click()
                    expect(decision).to_contain_text(
                        "Prepared protein not established", timeout=120000)
                    expect(decision).to_contain_text(proposal_id)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["protein"]["status"], "declined")
                    self.assertIsNone(state["placement"])
                    self.assertEqual(state["stages"], [])
                    expect(approve).to_be_disabled()
                    expect(page.locator(".stage-strip")).to_contain_text(
                        "No completed stage is currently established")
                    page.screenshot(path=str(ARTIFACTS / "declined-820.png"), full_page=True,
                                    animations="disabled")
                    page.get_by_role("button", name="Show workflow").click()
                    expect(page.locator("#researcher-workflow")).to_be_visible()
                    expect(page.locator(".protein-standing")).to_contain_text(
                        "Prepared protein not established")
                    print("PASS: browser proposal inspection, required decline and reduced desktop standing")
                finally:
                    browser.close()

    @unittest.skipUnless(os.environ.get("PIM_BROWSER_REMOTE") == "1",
                         "Live source review is run with PIM_BROWSER_REMOTE=1")
    def test_live_exact_rcsb_reference_has_visible_source_identity_and_structure(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1440, "height": 900})
                    page.goto(base, wait_until="domcontentloaded")
                    page.locator("#exact-source-kind").select_option("rcsb")
                    page.locator("#exact-identifier").fill("1CRN")
                    page.get_by_role("button", name="Inspect exact source").click()
                    expect(page.locator(".source-context")).to_contain_text(
                        "rcsb:1CRN", timeout=120000)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["study"]["selectedSourceId"], "rcsb:1CRN")
                    self.assertGreater(state["sourceModels"][0]["atomCount"], 300)
                    self.assertIsNone(state["protein"])
                    page.get_by_role("button", name="Inspect source structure").click()
                    expect(page.locator(".scene-subtitle")).to_contain_text("rcsb:1CRN")
                    expect(page.locator(".evidence-panel")).to_contain_text("rcsb:1CRN")
                    wait_for_selected_structure(page)
                    self.assertEqual(horizontally_overflowing_regions(page), [])
                    page.screenshot(path=str(ARTIFACTS / "source-rcsb-1440.png"),
                                    full_page=True, animations="disabled")
                    print("PASS: browser exact RCSB reference, source identity and real structure")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
