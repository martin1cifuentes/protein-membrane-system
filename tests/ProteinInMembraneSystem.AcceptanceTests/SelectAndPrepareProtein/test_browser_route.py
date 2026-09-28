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

HOST = Path(os.environ.get("PIM_BROWSER_HOST", str(ROOT / "out" / "host" / "ProteinInMembrane.Host.dll")))
WORKER = ROOT / "out" / "python" / "bin" / "python"
POLICY = ROOT / "config" / "policies" / "protein-slice1.json"
ARTIFACTS = ROOT / "out" / "browser-acceptance"


@contextmanager
def running_host(workspace: Path, policy: Path = POLICY, ppm: Path | None = None,
                 worker: Path | None = None):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    environment = dict(os.environ,
                       ASPNETCORE_URLS=base,
                       PIM_WORKER_PYTHON=str(worker or WORKER),
                       PIM_OWNER_SOURCE_ROOT=str(ROOT / "src" / "ProteinInMembrane.Host"),
                       PIM_WORKSPACE_ROOT=str(workspace),
                       PIM_POLICY_CATALOGUE=str(policy))
    if ppm is None:
        environment.pop("PIM_PPM_EXECUTABLE", None)
    else:
        environment["PIM_PPM_EXECUTABLE"] = str(ppm)
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
    page.get_by_role("button", name="Upload source").click()
    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
    page.get_by_label("Chain A").check()
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
        for path in (HOST, WORKER, POLICY, HOST.parent / "wwwroot" / "index.html"):
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
                    expect(page.get_by_role("button", name="Load source", exact=True)).to_be_disabled()
                    upload_and_choose(page, source)
                    expect(page.get_by_role("region", name="Protein task outcome"))\
                        .to_contain_text("Protein prepared", timeout=120000)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["protein"]["status"], "assessed")
                    self.assertEqual(state["study"]["modelIndex"], 0)
                    self.assertEqual(state["study"]["chainIds"], ["A"])
                    self.assertIsNone(state["placement"])
                    self.assertEqual(state["stages"], [])
                    self.assertIn(source.name, page.locator(".source-context").inner_text())
                    expect(page.get_by_role("heading", name=f"Protein before preparation · {source.name}"))\
                        .to_be_visible()
                    page.get_by_role("button", name="View prepared protein").click()
                    expect(page.locator(".workspace")).to_have_class("workspace workflow-open")
                    expect(page.get_by_role("heading", name=f"Prepared protein · {source.name}"))\
                        .to_be_visible()
                    expect(page.locator(".evidence-panel")).to_contain_text("Preparation complete")
                    expect(page.locator(".stage-strip")).to_contain_text("No minimized or equilibrated")
                    self.assertEqual(page.evaluate("async () => (await fetch('/api/state')).json()")
                                     ["inspection"]["subjectId"], state["protein"]["subjectId"])
                    wait_for_selected_structure(page)
                    self.assertEqual(horizontally_overflowing_regions(page), [])
                    page.screenshot(path=str(ARTIFACTS / "assessed-1440.png"), full_page=True,
                                    animations="disabled")
                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Membrane").click()
                    expect(page.locator("#researcher-workflow")).to_be_visible()
                    expect(page.locator("#membrane-workflow")).to_be_visible()
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
                    expect(page.locator(".review-progress")).to_contain_text("0 site choices", timeout=120000)
                    expect(page.locator(".review-other .review-site")).to_have_count(1)
                    page.locator(".review-other .review-site").click()
                    expect(page.locator(".review-site-scope")).to_contain_text("Missing-atom repair")
                    before_review = page.evaluate("async () => (await fetch('/api/state')).json()")
                    pending_id = before_review["protein"]["changes"][0]["id"]
                    self.assertTrue(next(item["enabled"] for item in before_review["actions"]
                                          if item["kind"] == "approvePreparationChange"
                                          and item["subjectId"] == pending_id))
                    page.locator(".review-option input").check()
                    approve = page.get_by_role("button", name="Approve repair")
                    expect(approve).to_be_enabled()
                    expect(page.locator(".review-evidence-current")).to_be_visible()
                    page.locator(".review-option.picked").get_by_role("button", name="Focus on this residue").click()
                    wait_for_selected_structure(page)
                    proposal_id = page.evaluate("async () => (await fetch('/api/state')).json()")\
                        ["protein"]["changes"][0]["id"]
                    page.screenshot(path=str(ARTIFACTS / "proposal-1440.png"), full_page=True,
                                    animations="disabled")

                    page.set_viewport_size({"width": 820, "height": 760})
                    scene_box = page.locator(".scene-panel").bounding_box()
                    evidence_box = page.locator(".evidence-panel").bounding_box()
                    self.assertIsNotNone(scene_box)
                    self.assertIsNotNone(evidence_box)
                    self.assertGreaterEqual(scene_box["width"], 350)
                    self.assertGreaterEqual(scene_box["height"], 450)
                    self.assertGreater(evidence_box["x"], scene_box["x"])
                    self.assertAlmostEqual(evidence_box["y"], scene_box["y"], delta=2)
                    self.assertEqual(horizontally_overflowing_regions(page), [])
                    expect(page.get_by_role("button", name="Reject required repair")).to_be_visible()
                    expect(page.locator(".scene-surface")).to_be_visible()
                    expect(page.locator(".compact-model-context")).to_be_visible()
                    approve.scroll_into_view_if_needed()
                    expect(approve).to_be_in_viewport()
                    page.screenshot(path=str(ARTIFACTS / "proposal-actions-820.png"), full_page=True,
                                    animations="disabled")
                    page.locator(".rail").evaluate("element => element.scrollTop = 0")
                    page.locator(".evidence-content").evaluate("element => element.scrollTop = 0")
                    page.screenshot(path=str(ARTIFACTS / "proposal-820.png"), full_page=True,
                                    animations="disabled")

                    page.get_by_role("button", name="Reject required repair").click()
                    expect(page.locator(".review-blocker-inline")).to_contain_text(
                        "required missing-atom repair was declined", timeout=120000)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["protein"]["status"], "declined")
                    self.assertEqual(state["preparationReview"]["decisions"][0]["options"][0]["disposition"], "declined")
                    self.assertEqual(state["preparationReview"]["decisions"][0]["options"][0]["proposalId"], proposal_id)
                    self.assertIsNone(state["placement"])
                    self.assertEqual(state["stages"], [])
                    expect(approve).to_have_count(0)
                    expect(page.locator(".stage-strip")).to_contain_text(
                        "No minimized or equilibrated")
                    page.screenshot(path=str(ARTIFACTS / "declined-820.png"), full_page=True,
                                    animations="disabled")
                    expect(page.locator("#researcher-workflow")).to_be_visible()
                    expect(page.get_by_role("region", name="Protein task outcome"))\
                        .to_contain_text("Protein preparation blocked")
                    print("PASS: browser exact repair evidence, required decline and reduced desktop standing")
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
                    page.get_by_role("button", name="Load source", exact=True).click()
                    expect(page.locator(".source-context")).to_contain_text(
                        "rcsb:1CRN", timeout=120000)
                    state = page.evaluate("async () => (await fetch('/api/state')).json()")
                    self.assertEqual(state["study"]["selectedSourceId"], "rcsb:1CRN")
                    self.assertGreater(state["sourceModels"][0]["atomCount"], 300)
                    self.assertIsNone(state["protein"])
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    expect(page.get_by_role("button", name="Inspect source structure")).to_have_count(0)
                    expect(page.locator(".scene-heading h1")).to_have_text("Source structure · 1CRN")
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
