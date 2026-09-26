"""A completed-system camera refits after the actor changes desktop width."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import CHROMIUM, TINY_STRUCTURE, controlled_account_server  # noqa: E402
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/ExportCompletedMinimizedStage"))
from test_browser_export import selected_stage_account  # noqa: E402


# A wide, identified bonded pair makes the camera-radius error observable
# without requiring another 85,000-atom render or scientific provider run.
WIDE_STRUCTURE = TINY_STRUCTURE.replace(b"   1.500", b" 120.000")


def fitted_camera(page) -> dict:
    page.wait_for_function("""() => {
        const mount = document.querySelector('.viewer-mount');
        if (!mount || mount.dataset.cameraReady !== 'true') return false;
        const bounds = mount.getBoundingClientRect();
        return Math.abs(Number(mount.dataset.cameraViewportWidth) - bounds.width) <= 4 &&
            Math.abs(Number(mount.dataset.cameraViewportHeight) - bounds.height) <= 4;
    }""", timeout=30000)
    return page.evaluate("""() => {
        const mount = document.querySelector('.viewer-mount');
        const camera = mount[Symbol.for('molstar.viewer')].plugin.canvas3d.camera;
        return {radius: camera.state.radius, viewport: camera.viewport};
    }""")


class CompletedStageCameraResizeTests(unittest.TestCase):
    def test_widening_after_reduced_desktop_refits_selected_system(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(selected_stage_account())
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                args=["--disable-dev-shm-usage", "--use-gl=angle", "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.route("**/api/structures/tiny?format=pdb", lambda route: route.fulfill(
                    status=200, content_type="chemical/x-pdb", body=WIDE_STRUCTURE))
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                expect(page.locator(".viewer-mount canvas")).to_have_count(1)
                before = fitted_camera(page)
                page.set_viewport_size({"width": 820, "height": 760})
                narrow = fitted_camera(page)
                page.set_viewport_size({"width": 1672, "height": 941})
                after = fitted_camera(page)
                self.assertTrue(all(math.isfinite(item["radius"]) and item["radius"] > 0
                                    for item in (before, narrow, after)))
                # A narrow viewport may clip the distant patch edges; it must
                # retain a prominent central system instead of shrinking it.
                self.assertLessEqual(narrow["radius"], before["radius"] * 1.1)
                self.assertAlmostEqual(after["radius"], before["radius"], delta=0.1)
                self.assertGreater(after["viewport"]["width"], narrow["viewport"]["width"] * 2)
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
