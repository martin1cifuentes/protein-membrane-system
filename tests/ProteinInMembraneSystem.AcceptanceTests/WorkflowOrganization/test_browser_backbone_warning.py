"""Rendered guidance for observed protein-backbone completeness.

These controlled accounts exercise browser applicability only. The source
worker and preparation owner remain responsible for atom observations and
scientific refusal; this fixture does not establish a prepared protein.
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import (  # noqa: E402
    CHROMIUM, DIST, account, controlled_account_server, inspection,
)


def residue(model: int, chain: str, number: int, kind: str, complete: bool, name: str) -> dict:
    return {"address": {"model": model, "chain": chain, "residue": number,
                        "insertionCode": "", "copyId": ""},
            "name": name, "residueKind": kind,
            "backboneHeavyAtomsComplete": complete, "alternateLocations": []}


def source_model(index: int, chain_names: list[str], residues: list[dict]) -> dict:
    return {"index": index,
            "chains": [{"name": chain, "residueCount": sum(r["address"]["chain"] == chain
                                                         for r in residues), "atomCount": 4}
                       for chain in chain_names],
            "assemblies": [], "partners": [], "residues": residues,
            "atomCount": len(residues) * 4}


def warning_account() -> dict:
    value = account()
    value["attempt"] = None
    value["stages"] = []
    value["study"]["selectedSourceId"] = "fixture:backbone-guidance"
    value["study"]["selectedSourceKind"] = "upload"
    value["sourceCandidates"] = [{"id": "fixture:backbone-guidance",
                                  "label": "Backbone guidance presentation fixture",
                                  "kind": "upload", "provenance": "Controlled account fixture",
                                  "limitations": []}]
    value["inspection"] = inspection("fixture:backbone-guidance", "structuralSource")
    value["inspection"]["omittedMolecules"] = []
    value["sourceModels"] = [
        source_model(0, ["A"], [residue(0, "A", 1, "protein", True, "ALA")]),
        source_model(1, ["A"], [
            residue(1, "A", 1, "protein", True, "ALA"),
            residue(1, "A", 2, "solvent", False, "HOH"),
            residue(1, "A", 3, "heterogen", False, "LIG"),
        ]),
        source_model(2, ["A", "B"], [
            residue(2, "A", 1, "protein", True, "ALA"),
            # Chain B's observed ALA lacks its required backbone O atom.
            residue(2, "B", 1, "protein", False, "ALA"),
        ]),
    ]
    return value


def open_browser(playwright):
    browser = playwright.chromium.launch(
        executable_path=CHROMIUM, headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"],
    )
    return browser, browser.new_page(viewport={"width": 1280, "height": 800})


class BackboneWarningBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert (DIST / "index.html").is_file(), "Build browser/ before these focused checks"

    def test_complete_protein_with_and_without_nonprotein_components_has_no_warning(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(warning_account())
            browser, page = open_browser(playwright)
            try:
                page.goto(base, wait_until="domcontentloaded")
                warning = page.locator(".notice.warning").filter(has_text="backbone heavy-atom")
                page.locator("#model-index").select_option("0")
                expect(warning).to_have_count(0)
                expect(page.get_by_text("Choose at least one observed chain copy.")).to_be_visible()
                page.get_by_label("Chain A", exact=True).check()
                expect(warning).to_have_count(0)

                page.locator("#model-index").select_option("1")
                expect(warning).to_have_count(0)
                expect(page.get_by_label("Chain A", exact=True)).not_to_be_checked()
                page.get_by_label("Chain A", exact=True).check()
                expect(warning).to_have_count(0)
                self.assertEqual(host.commands, [], "Draft selection cannot submit assessment")
            finally:
                browser.close()

    def test_missing_backbone_warning_tracks_selected_chains_and_model_changes(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(warning_account())
            browser, page = open_browser(playwright)
            try:
                page.goto(base, wait_until="domcontentloaded")
                warning = page.locator(".notice.warning").filter(has_text="backbone heavy-atom")
                page.locator("#model-index").select_option("2")
                expect(warning).to_have_count(0)
                expect(page.get_by_text("Choose at least one observed chain copy.")).to_be_visible()

                page.get_by_label("Chain A", exact=True).check()
                expect(warning).to_have_count(0)
                page.get_by_label("Chain B", exact=True).check()
                expect(warning).to_contain_text("selected protein chains")
                expect(warning).to_contain_text("no missing backbone is silently reconstructed")
                page.get_by_label("Chain B", exact=True).uncheck()
                expect(warning).to_have_count(0)
                page.get_by_label("Chain A", exact=True).uncheck()
                page.get_by_label("Chain B", exact=True).check()
                expect(warning).to_have_count(1)

                page.locator("#model-index").select_option("0")
                expect(warning).to_have_count(0)
                page.get_by_label("Chain A", exact=True).check()
                expect(warning).to_have_count(0)
                page.locator("#model-index").select_option("2")
                expect(warning).to_have_count(0)
                expect(page.get_by_label("Chain B", exact=True)).not_to_be_checked()
                page.get_by_label("Chain B", exact=True).check()
                expect(warning).to_have_count(1)
                self.assertEqual(host.commands, [], "Warning presentation cannot submit assessment")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
