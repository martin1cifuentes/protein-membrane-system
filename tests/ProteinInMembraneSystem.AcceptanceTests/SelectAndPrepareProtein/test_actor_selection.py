"""Actor selection identity through the real source worker and isolated browser host.

The mmCIF is an authored software fixture with two coordinate models, an
assembly, and two hemes. It establishes interface/command identity only; it is
not evidence that the construct is biologically or chemically suitable.
"""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from urllib.parse import urlencode

from playwright.sync_api import expect, sync_playwright

from test_browser_route import ROOT, chromium, running_host, horizontally_overflowing_regions


FIXTURE = Path(__file__).parent / "fixtures" / "actor-selection-identities.cif"
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "actor-selection"


class ActorSelectionBrowserTests(unittest.TestCase):
    def test_model_assembly_membership_partner_scope_and_submitted_identity(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            with running_host(Path(temporary) / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    submitted = []
                    page.on("request", lambda request: submitted.append(json.loads(request.post_data))
                            if request.url.endswith("/api/commands") and request.post_data else None)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(FIXTURE))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    expect(page.get_by_role("heading", name=f"Source structure · {FIXTURE.name}"))\
                        .to_be_visible()
                    models = page.locator("#model-index")
                    expect(models).to_have_value("")
                    expect(models.locator("option[value='0']")).to_contain_text("Source model 7")
                    expect(models.locator("option[value='1']")).to_contain_text("Source model 9")
                    expect(page.locator("#assembly-choice")).to_have_count(0)

                    models.select_option("0")
                    expect(page.locator("#assembly-choice")).to_have_value("")
                    expect(page.locator("#assembly-choice option[value='']")).to_be_disabled()
                    expect(page.get_by_role("button", name="Assess selected protein")).to_have_count(0)
                    page.locator("#assembly-choice").select_option("1")
                    expect(page.get_by_label("Chain A · copy A1")).to_be_visible()
                    expect(page.get_by_label("Chain A · copy A2")).to_be_visible()
                    expect(page.locator(".scene-subtitle")).to_contain_text(
                        "Source model 7 · biological assembly 1", timeout=30000)
                    expect(page.locator(".chain-color-swatch")).to_have_count(2, timeout=30000)
                    selected_source = page.request.get(base + "/api/state").json()["study"]["selectedSourceId"]
                    preview = page.request.get(base + "/api/source-preview?" + urlencode({
                        "sourceId": selected_source, "modelIndex": 0, "assemblyId": "1"})).json()
                    self.assertEqual(preview["chainIds"], ["A1", "B1", "A2", "B2"])
                    self.assertEqual(page.request.get(base + preview["structureUrl"]).status, 200)
                    camera_target = lambda: page.locator(".viewer-mount").evaluate(
                        "el => Array.from(el[Symbol.for('molstar.viewer')].plugin.canvas3d.camera.getSnapshot().target)")
                    page.locator(".chain-choice").filter(has_text="copy A1").get_by_role(
                        "button", name="Focus in viewer").click()
                    page.wait_for_timeout(600)
                    first_copy_target = camera_target()
                    page.locator(".chain-choice").filter(has_text="copy A2").get_by_role(
                        "button", name="Focus in viewer").click()
                    page.wait_for_timeout(600)
                    second_copy_target = camera_target()
                    self.assertGreater(max(abs(a - b) for a, b in zip(first_copy_target,
                                                                          second_copy_target)), 0.1)
                    expect(page.get_by_label("Chain B · copy B1")).to_have_count(0)
                    expect(page.get_by_label("Chain B · copy B2")).to_have_count(0)
                    expect(page.locator(".partner-choice")).to_have_count(1)
                    expect(page.locator(".partner-choice")).to_contain_text("Heme (HEM) · Chain B · Residue 147")
                    page.get_by_label("Chain A · copy A1").check()
                    expect(page.locator(".partner-choice")).to_have_count(2)
                    heme_a = page.locator(".partner-choice").filter(has_text="Chain A · Residue 142")
                    heme_b = page.locator(".partner-choice").filter(has_text="Chain B · Residue 147")
                    expect(heme_a).to_contain_text("Applies to copy A1")
                    expect(heme_b).to_contain_text("every assembly copy (B1, B2)")
                    page.get_by_label("Chain A · copy A2").check()
                    expect(heme_a).to_contain_text("every assembly copy (A1, A2)")
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        self.assertEqual(horizontally_overflowing_regions(page), [])
                        page.screenshot(path=str(ARTIFACTS / f"assembly-partners-{width}.png"),
                                        full_page=True, animations="disabled")
                        heme_b.scroll_into_view_if_needed()
                        page.screenshot(path=str(ARTIFACTS / f"assembly-partner-decisions-{width}.png"),
                                        full_page=True, animations="disabled")
                        page.locator(".rail").evaluate("element => element.scrollTop = 0")

                    # Changing models discards a noncorresponding draft. Deposited coordinates
                    # remain a deliberate alternative to the authored assembly.
                    models.select_option("1")
                    expect(page.locator("#assembly-choice")).to_have_value("")
                    expect(page.get_by_label("Chain A · copy A1")).to_have_count(0)
                    page.locator("#assembly-choice").select_option("deposited")
                    page.get_by_label("Chain A", exact=True).check()
                    expect(page.locator(".partner-choice")).to_have_count(2)
                    for card in page.locator(".partner-choice").all():
                        card.get_by_label("Exclude").check()
                    expect(page.get_by_label("Protein membership to assess"))\
                        .to_contain_text("Source model 9 · deposited coordinates")
                    expect(page.get_by_label("Protein membership to assess"))\
                        .to_contain_text("Heme (HEM) · Chain A · Residue 142")
                    expect(page.get_by_label("Protein membership to assess"))\
                        .to_contain_text("Heme (HEM) · Chain B · Residue 147")
                    page.screenshot(path=str(ARTIFACTS / "deposited-summary-820.png"),
                                    full_page=True, animations="disabled")
                    with page.expect_request(lambda request: request.url.endswith("/api/commands") and
                                             request.post_data_json.get("kind") == "selectProteinModel"):
                        page.get_by_role("button", name="Assess selected protein").click()
                    command = next(request for request in submitted if request["kind"] == "selectProteinModel")
                    self.assertEqual(command["data"]["modelIndex"], 1)
                    self.assertIsNone(command["data"]["biologicalAssemblyId"])
                    self.assertEqual(command["data"]["chains"], [{"sourceChain": "A", "copyId": "A"}])
                    self.assertEqual({item["sourceId"] for item in command["data"]["partners"]},
                                     {"1:A:142::HEM", "1:B:147::HEM"})
                    self.assertTrue(all(not item["retain"] for item in command["data"]["partners"]))
                    page.wait_for_function("async () => (await (await fetch('/api/state')).json()).study.modelIndex === 1",
                                           timeout=120000)
                    state = page.request.get(base + "/api/state").json()
                    self.assertEqual(state["study"]["modelIndex"], 1)
                    self.assertEqual(state["study"]["biologicalAssemblyId"], None)
                    self.assertEqual(state["study"]["chainIds"], ["A"])
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
