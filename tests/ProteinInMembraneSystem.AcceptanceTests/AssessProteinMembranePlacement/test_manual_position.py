"""Real small-construct manual placement through the isolated Host and browser.

The two-alanine PDB is a software fixture for coordinate and interaction checks;
this test makes no claim about biological membrane suitability.
"""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import time
import unittest

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" / "SelectAndPrepareProtein"))
from test_browser_route import chromium, running_host, two_alanines  # noqa: E402

POLICY = ROOT / "config" / "policies" / "protein-membrane-current.json"
ARTIFACTS = ROOT / "out" / "reconciled-design" / "manual-placement"


def state(page):
    return page.evaluate("async () => (await (await fetch('/api/state', {cache: 'no-store'})).json())")


def until(page, predicate, description, seconds=120):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        current = state(page)
        if predicate(current):
            return current
        time.sleep(0.15)
    raise AssertionError(f"Timed out waiting for {description}: {current.get('placementTask')}; "
                         f"notices={current.get('notices', [])[-3:]}")


def prepare_with_current_plan(page):
    planned = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing")
                    in ("ready", "failed", "partial"), "protein starting-state plan")
    if planned["preparationPlan"]["standing"] != "ready":
        raise AssertionError(f"No complete checked starting-state plan: {planned['preparationPlan']}")
    page.get_by_role("button", name="Prepare with recommendations").click()
    outcome = until(page, lambda current: (current.get("proteinTask") or {}).get("standing")
                    in ("assessed", "failed", "unavailable"), "authorized prepared protein")
    if outcome["proteinTask"]["standing"] != "assessed":
        raise AssertionError(f"Protein preparation did not establish the result: {outcome['proteinTask']}")
    return outcome


class ManualPlacementBrowserTests(unittest.TestCase):
    def test_position_check_is_automatic_and_adoption_uses_current_transform(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-software-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            with running_host(directory / "workspace", policy=POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    prepare_with_current_plan(page)
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Membrane").click()
                    for side in ("Upper", "Lower"):
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                    page.get_by_role("button", name="Use this membrane").click()
                    until(page, lambda current: (current.get("membrane") or {}).get("status") == "assessed",
                          "checked membrane")
                    page.get_by_role("navigation", name="Research work areas").get_by_role("button", name="Placement").click()
                    initial = until(page, lambda current: current.get("placement") is not None and
                                    current["placement"]["status"] in ("supported", "unsupported", "notEstablished"),
                                    "automatic position check")
                    self.assertEqual(initial["placement"]["status"], "supported", initial["placement"])
                    self.assertEqual(initial["placement"]["transform"]["startingPosition"], "center")
                    expect(page.get_by_role("button", name="Check position now")).to_have_count(0)
                    page.get_by_label("Move X (Å)").fill("")
                    expect(page.get_by_label("Current placement status")).to_contain_text("Incomplete input")
                    expect(page.get_by_label("Current placement status")).to_contain_text("Move X")
                    expect(page.get_by_role("button", name="Use this position")).to_have_count(0)
                    expect(page.get_by_role("button", name="Use earlier checked position")).to_be_enabled()
                    page.wait_for_timeout(650)
                    self.assertEqual(state(page)["placement"]["proposalId"], initial["placement"]["proposalId"])
                    page.locator("#starting-position").select_option("upper")
                    page.get_by_label("Move X (Å)").fill("3.5")
                    expect(page.get_by_role("button", name="Use this position")).to_have_count(0)
                    changed = until(page, lambda current: current.get("placement") is not None and
                                    current["placement"]["proposalId"] != initial["placement"]["proposalId"] and
                                    current["placement"]["status"] in ("supported", "unsupported", "notEstablished"),
                                    "changed exact position")
                    self.assertEqual(changed["placement"]["status"], "supported", changed["placement"])
                    self.assertEqual(changed["placement"]["transform"]["startingPosition"], "upper")
                    self.assertAlmostEqual(changed["placement"]["transform"]["offsetXAngstrom"], 3.5)
                    expect(page.get_by_role("button", name="Use this position")).to_be_enabled()
                    expect(page.get_by_role("button", name="View positioned protein")).to_have_count(0)
                    page.get_by_role("button", name="Use this position").click()
                    adopted = until(page, lambda current: current["study"]["adoptedPlacementProposalId"] ==
                                    changed["placement"]["proposalId"], "adopted exact position")
                    self.assertEqual(adopted["placement"]["proposalId"], changed["placement"]["proposalId"])
                    expect(page.locator("#starting-position")).to_have_value("upper")
                    expect(page.get_by_label("Move X (Å)")).to_have_value("3.5")
                    self.assertTrue(any(item["enabled"] for item in adopted["actions"]
                                        if item["kind"] == "buildAndMinimize"),
                                    "At least one identified construction method should be available")
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Preparation").click()
                    expect(page.get_by_text("Method: PACKMOL-Memgen")).to_be_visible()
                    page.locator(".construction-route-list").get_by_text("Method details").first.click()
                    expect(page.locator(".construction-route-list")).to_contain_text("900 seconds")
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Placement").click()
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                            .to_have_count(1, timeout=60000)
                        page.wait_for_timeout(250)
                        guide_and_legend = page.evaluate("""() => {
                          const mount = document.querySelector('.placement-scene-surface .viewer-mount');
                          const legend = document.querySelector('.placement-scene-legend');
                          const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                          const camera = viewer.plugin.canvas3d.camera;
                          const lowerCorners = [[-18, -18, -23], [-18, 18, -23],
                            [18, -18, -23], [18, 18, -23]];
                          const mountTop = mount.getBoundingClientRect().top;
                          return {planeBottom: mountTop + Math.max(...lowerCorners.map(point =>
                              camera.project(new Float32Array(4), new Float32Array(point))[1])),
                            legendTop: legend.getBoundingClientRect().top};
                        }""")
                        self.assertGreater(guide_and_legend["legendTop"],
                                           guide_and_legend["planeBottom"] + 2,
                                           f"The preview legend obscured a boundary at {width}px: {guide_and_legend}")
                        page.screenshot(path=str(ARTIFACTS / f"manual-adopted-{width}.png"),
                                        animations="disabled")
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
