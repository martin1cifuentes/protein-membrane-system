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

from test_browser_route import ROOT, chromium, running_host, horizontally_overflowing_regions, two_alanines


FIXTURE = Path(__file__).parent / "fixtures" / "actor-selection-identities.cif"
WATER_FIXTURE = Path(__file__).parent / "fixtures" / "actor-selection-ordinary-water.cif"
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "actor-selection"


def viewer_chain_colors(page) -> dict[str, str]:
    """Read the loaded Mol* cartoon theme for each actual source chain/copy."""
    return page.locator(".viewer-mount").evaluate("""element => {
      const viewer = element[Symbol.for('molstar.viewer')];
      const components = viewer.plugin.managers.structure.hierarchy.current.structures[0].components;
      const polymer = components.find(item => item.cell.transform.tags?.includes('structure-component-static-polymer'));
      const structure = polymer.cell.obj.data;
      const cartoon = polymer.representations.find(item => item.cell.params?.values?.type?.name === 'cartoon');
      const theme = cartoon.cell.obj.data.repr.theme.color;
      const colors = {};
      for (const unit of structure.units) for (const elementIndex of unit.elements) {
        const hierarchy = unit.model.atomicHierarchy;
        const chainIndex = hierarchy.chainAtomSegments.index[elementIndex];
        const chain = hierarchy.chains.auth_asym_id.value(chainIndex);
        if (colors[chain]) continue;
        const value = theme.color({kind: 'element-location', structure, unit, element: elementIndex}, false);
        colors[chain] = `rgb(${value >> 16 & 255}, ${value >> 8 & 255}, ${value & 255})`;
      }
      return colors;
    }""")


def mixed_partner_source() -> str:
    """Small source-only fixture with distinct water, ion and heme addresses."""
    partners = (
        "HETATM   13  O   HOH A  10       9.000   4.000   0.000  1.00 20.00           O\n"
        "HETATM   14  O   HOH B  20      10.000   4.000   0.000  1.00 20.00           O\n"
        "HETATM   15  NA   NA C  30      11.000   4.000   0.000  1.00 20.00          NA\n"
        "HETATM   16  FE  HEM D  40      12.000   4.000   0.000  1.00 20.00          FE\n"
        "HETATM   17  FE  HEM E  41      13.000   4.000   0.000  1.00 20.00          FE\n"
    )
    return two_alanines(1, "A").replace("ENDMDL\n", "") + partners + "ENDMDL\nEND\n"


def expand_partner_group(group) -> None:
    details = group.locator("details.partner-members")
    if not details.evaluate("element => element.open"):
        details.locator("summary").click()


class ActorSelectionBrowserTests(unittest.TestCase):
    def test_bulk_receipt_counts_54_current_waters_and_withdraws_on_override(self):
        """Controlled source-membership UI case; no preparation is claimed."""
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "54-waters-software-fixture.pdb"
            backbone = two_alanines(1, "A").replace("ENDMDL\n", "")
            waters = "".join(
                f"HETATM{serial:5d}  O   HOH A{serial:4d}    {float(serial):8.3f}{4.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}           O\n"
                for serial in range(20, 74)
            )
            source.write_text(backbone + waters + "ENDMDL\nEND\n")
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A", exact=True).check()
                    group = page.get_by_role("region", name="Source waters (54)")
                    expect(group).to_contain_text("0 kept · 0 excluded · 54 undecided")
                    group.get_by_role("button", name="Exclude all source waters").click()
                    expect(group.get_by_role("status")).to_have_text(
                        "54 source waters excluded from preparation.")
                    expect(group).to_contain_text("0 kept · 54 excluded · 0 undecided")
                    ARTIFACTS.mkdir(parents=True, exist_ok=True)
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        group.scroll_into_view_if_needed()
                        self.assertEqual(horizontally_overflowing_regions(page), [])
                        page.screenshot(path=str(ARTIFACTS / f"54-waters-receipt-{width}.png"))
                    group.get_by_role("button", name="Exclude all source waters").focus()
                    page.keyboard.press("Enter")
                    expect(group.get_by_role("status")).to_contain_text("No additional source waters excluded")
                    expand_partner_group(group)
                    group.locator(".partner-choice").first.get_by_label("Keep").check()
                    expect(group).to_contain_text("1 kept · 53 excluded · 0 undecided")
                    expect(group.get_by_role("status")).to_have_count(0)
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        self.assertEqual(horizontally_overflowing_regions(page), [])
                        page.screenshot(path=str(ARTIFACTS / f"54-waters-feedback-{width}.png"))
                finally:
                    browser.close()

    def test_category_actions_change_only_current_partner_draft_and_allow_one_override(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "mixed-partners-software-fixture.pdb"
            source.write_text(mixed_partner_source())
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    submitted = []
                    page.on("request", lambda request: submitted.append(json.loads(request.post_data))
                            if request.url.endswith("/api/commands") and request.post_data else None)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text(
                        "Structure displayed", timeout=120000)
                    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
                    page.get_by_label("Chain A", exact=True).check()
                    waters = page.get_by_role("region", name="Source waters (2)")
                    ions = page.get_by_role("region", name="Ions (1)")
                    ligands = page.get_by_role("region", name="Ligands and cofactors (2)")
                    waters.get_by_role("button", name="Exclude all source waters").click()
                    expect(waters.get_by_role("status")).to_have_text(
                        "2 source waters excluded from preparation.")
                    summary = page.get_by_label("Protein membership to assess")
                    expect(summary).to_contain_text("Source waters: 0 kept · 2 excluded · 0 to decide")
                    expect(waters).to_contain_text("0 kept · 2 excluded · 0 undecided")
                    waters.get_by_role("button", name="Exclude all source waters").click()
                    expect(waters.get_by_role("status")).to_contain_text("No additional source waters excluded")
                    expect(summary).to_contain_text("Ions: 0 kept · 0 excluded · 1 to decide")
                    expect(summary).to_contain_text("Ligands and cofactors: 0 kept · 0 excluded · 2 to decide")
                    ions.get_by_role("button", name="Keep all ions").click()
                    ligands.get_by_role("button", name="Exclude all ligands and cofactors").click()
                    expand_partner_group(waters)
                    expand_partner_group(ions)
                    expand_partner_group(ligands)
                    water_b = waters.locator(".partner-choice").filter(
                        has_text="Water (HOH) · Chain B · Residue 20")
                    water_b.get_by_label("Keep").check()
                    expect(waters.get_by_role("status")).to_have_count(0)
                    expect(ions).to_contain_text("Sodium ion (NA) · Chain C · Residue 30")
                    expect(ligands).to_contain_text("Heme (HEM) · Chain D · Residue 40")
                    expect(ligands).to_contain_text("Heme (HEM) · Chain E · Residue 41")
                    expect(summary).to_contain_text("Source waters: 1 kept · 1 excluded · 0 to decide")
                    expect(summary).to_contain_text("Ions: 1 kept · 0 excluded · 0 to decide")
                    expect(summary).to_contain_text("Ligands and cofactors: 0 kept · 2 excluded · 0 to decide")
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        self.assertEqual(horizontally_overflowing_regions(page), [])
                    with page.expect_response(lambda response: response.url.endswith("/api/commands") and
                                              response.request.post_data_json.get("kind") == "selectProteinModel") as selected:
                        page.get_by_role("button", name="Assess selected protein").click()
                    self.assertEqual(selected.value.status, 200, selected.value.text())
                    command = next(item for item in submitted if item["kind"] == "selectProteinModel")
                    self.assertEqual(command["data"]["chains"], [{"sourceChain": "A", "copyId": "A"}])
                    self.assertEqual({item["sourceId"]: item["retain"] for item in command["data"]["partners"]}, {
                        "0:A:10::HOH": False, "0:B:20::HOH": True, "0:C:30::NA": True,
                        "0:D:40::HEM": False, "0:E:41::HEM": False,
                    })
                finally:
                    browser.close()

    def test_partner_only_water_follows_selected_assembly_frames(self):
        """A controlled water fixture checks actor choice and submitted copy scope."""
        with tempfile.TemporaryDirectory() as temporary:
            with running_host(Path(temporary) / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1280, "height": 800})
                    submitted = []
                    page.on("request", lambda request: submitted.append(json.loads(request.post_data))
                            if request.url.endswith("/api/commands") and request.post_data else None)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(WATER_FIXTURE))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text(
                        "Structure displayed", timeout=120000)
                    page.locator("#model-index").select_option("0")
                    page.locator("#assembly-choice").select_option("1")
                    page.get_by_label("Chain A · copy A1").check()
                    page.get_by_label("Chain A · copy A2").check()
                    waters = page.get_by_role("region", name="Source waters (2)")
                    expect(waters).to_contain_text("Source waters (2)")
                    waters.get_by_role("button", name="Exclude all source waters").click()
                    expand_partner_group(waters)
                    partner_a = page.locator(".partner-choice").filter(
                        has_text="Water (HOH) · Chain A · Residue 142")
                    partner_b = page.locator(".partner-choice").filter(
                        has_text="Water (HOH) · Chain B · Residue 147")
                    expect(partner_b).to_contain_text("selected protein assembly frame")
                    expect(partner_b).to_contain_text("A1, A2")
                    expect(partner_a).to_contain_text("every assembly copy (A1, A2)")
                    expect(partner_a.get_by_label("Exclude")).to_be_checked()
                    expect(partner_b.get_by_label("Exclude")).to_be_checked()
                    partner_b.get_by_label("Keep").check()
                    expect(page.get_by_label("Protein membership to assess")).to_contain_text(
                        "Source waters: 1 kept · 1 excluded · 0 to decide")
                    with page.expect_response(lambda response: response.url.endswith("/api/commands") and
                                              response.request.post_data_json.get("kind") == "selectProteinModel") as selection_response:
                        page.get_by_role("button", name="Assess selected protein").click()
                    response = selection_response.value
                    self.assertEqual(response.status, 200, response.text())
                    command = next(request for request in submitted if request["kind"] == "selectProteinModel")
                    self.assertEqual(command["data"]["chains"], [
                        {"sourceChain": "A", "copyId": "A1"},
                        {"sourceChain": "A", "copyId": "A2"},
                    ])
                    self.assertEqual(command["data"]["biologicalAssemblyId"], "1")
                    self.assertEqual({item["sourceId"]: item["retain"]
                                      for item in command["data"]["partners"]},
                                     {"0:A:142::HOH": False, "0:B:147::HOH": True})
                finally:
                    browser.close()

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
                    mapped_colors = viewer_chain_colors(page)
                    for copy_id in ("A1", "A2"):
                        row = page.locator(".chain-choice").filter(has_text=f"copy {copy_id}")
                        expect(row.locator("input[type=checkbox]")).to_have_count(1)
                        expect(row.locator("label .chain-color-swatch")).to_have_count(0)
                        expect(row.locator(".chain-color-swatch")).to_have_count(1)
                        self.assertIn(copy_id, mapped_colors)
                        expect(row.locator(".chain-color-swatch")).to_have_css(
                            "background-color", mapped_colors[copy_id])
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
                    expect(page.get_by_role("button", name="Show whole system")).to_have_count(0)
                    page.get_by_role("button", name="Fit structure").click()
                    page.wait_for_function("""previous => {
                      const viewer = document.querySelector('.viewer-mount')[Symbol.for('molstar.viewer')];
                      const target = viewer.plugin.canvas3d.camera.getSnapshot().target;
                      return Math.max(...target.map((value, index) => Math.abs(value - previous[index]))) > 0.1;
                    }""", arg=second_copy_target)
                    whole_target = camera_target()
                    self.assertGreater(max(abs(a - b) for a, b in zip(whole_target,
                                                                          second_copy_target)), 0.1)
                    expect(page.get_by_label("Chain B · copy B1")).to_have_count(0)
                    expect(page.get_by_label("Chain B · copy B2")).to_have_count(0)
                    expect(page.get_by_role("region", name="Ligands and cofactors (1)"))\
                        .to_contain_text("Ligands and cofactors (1)")
                    expand_partner_group(page.get_by_role("region", name="Ligands and cofactors (1)"))
                    expect(page.locator(".partner-choice")).to_have_count(1)
                    expect(page.locator(".partner-choice")).to_contain_text("Heme (HEM) · Chain B · Residue 147")
                    page.get_by_label("Chain A · copy A1").check()
                    expect(page.locator(".partner-choice")).to_have_count(2)
                    heme_a = page.locator(".partner-choice").filter(has_text="Chain A · Residue 142")
                    heme_b = page.locator(".partner-choice").filter(has_text="Chain B · Residue 147")
                    expect(heme_a).to_contain_text("Applies to copy A1")
                    expect(heme_b).to_contain_text("selected protein assembly frame (A1)")
                    page.get_by_label("Chain A · copy A2").check()
                    expect(heme_a).to_contain_text("every assembly copy (A1, A2)")
                    expect(heme_b).to_contain_text("selected protein assembly frame (A1, A2)")
                    for width in (1024, 820):
                        page.set_viewport_size({"width": width, "height": 768})
                        expect(page.locator(".chain-choice")).to_have_count(2)
                        for copy_id in ("A1", "A2"):
                            row = page.locator(".chain-choice").filter(has_text=f"copy {copy_id}")
                            expect(row.locator("input[type=checkbox]")).to_be_checked()
                            expect(row.locator(".chain-color-swatch")).to_have_css(
                                "background-color", mapped_colors[copy_id])
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Membrane").click()
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Protein").click()
                    expect(page.locator(".chain-choice")).to_have_count(2)
                    expect(page.get_by_label("Chain A · copy A1")).to_be_checked()
                    expect(page.get_by_label("Chain A · copy A2")).to_be_checked()
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
                    ligands = page.get_by_role("region", name="Ligands and cofactors (2)")
                    expect(ligands).to_contain_text("Ligands and cofactors (2)")
                    expect(page.get_by_label("Protein membership to assess")).to_contain_text(
                        "Ligands and cofactors: 0 kept · 0 excluded · 2 to decide")
                    ligands.get_by_role("button", name="Exclude all ligands and cofactors").click()
                    expand_partner_group(ligands)
                    expect(page.locator(".partner-choice")).to_have_count(2)
                    expect(page.get_by_label("Protein membership to assess"))\
                        .to_contain_text("Source model 9 · deposited coordinates")
                    expect(page.get_by_label("Protein membership to assess")).to_contain_text(
                        "Ligands and cofactors: 0 kept · 2 excluded · 0 to decide")
                    expect(ligands.locator(".partner-choice").filter(has_text="Chain A · Residue 142"))\
                        .to_contain_text("Heme (HEM)")
                    expect(ligands.locator(".partner-choice").filter(has_text="Chain B · Residue 147"))\
                        .to_contain_text("Heme (HEM)")
                    page.screenshot(path=str(ARTIFACTS / "deposited-summary-820.png"),
                                    full_page=True, animations="disabled")
                    with page.expect_response(lambda response: response.url.endswith("/api/commands") and
                                              response.request.post_data_json.get("kind") == "selectProteinModel") as selection_response:
                        page.get_by_role("button", name="Assess selected protein").click()
                    response = selection_response.value
                    self.assertEqual(response.status, 200, response.text())
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
