"""Connected 1L2Y Memgen routes and independent ZIP read-back.

The uploaded fixture is the exact RCSB 1L2Y solution-NMR PDB download. The
researcher explicitly selects deposited model 1 (account index 0), chain A,
and an upper-side or centered manual four-component or pure-DOPC pose. The
shared runner and inspector own the Build and minimize attempt, completed
stage, export and molecular read-back; this test supplies only source-specific
browser preparation. A passing preflight alone is not a provider result.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright

from test_memgen_connected import (
    ARTIFACTS, HOST, POLICY, ROOT, ROUTE_ID, ROUTE_LABEL,
    ConnectedMemgenCrambinRoute, await_state, choose_membrane_and_manual_pose,
    chromium, current_state, running_host_without_ppm,
)
from test_memgen_trpcage_preflight import (
    MODEL_ATOM_COUNT, MODEL_COUNT, MODEL_HYDROGEN_COUNT, PROVENANCE_NOTE,
    SOURCE, SOURCE_SHA256, selected_source_heavy_atom_ids,
)


def prepare_exact_trpcage(page):
    """Prepare exactly one deposited 1L2Y model in this connected Host run."""
    page.locator("#source-upload").set_input_files(str(SOURCE))
    page.locator("#upload-provenance").select_option("experimental")
    page.locator("#upload-provenance-note").fill(PROVENANCE_NOTE)
    page.get_by_role("button", name="Upload source").click()
    expect(page.locator(".source-context")).to_contain_text(
        "Structure displayed", timeout=180000)
    source = current_state(page)
    assert source["study"]["uploadProvenance"] == "experimental"
    assert source["study"]["uploadProvenanceNote"] == PROVENANCE_NOTE
    assert len(source["sourceModels"]) == MODEL_COUNT
    assert [model["index"] for model in source["sourceModels"]] == list(range(MODEL_COUNT))
    first = source["sourceModels"][0]
    assert first["sourceModelId"] == "1" and first["atomCount"] == MODEL_ATOM_COUNT
    assert first["chains"] == [{"name": "A", "residueCount": 20,
                                 "atomCount": MODEL_ATOM_COUNT}]
    assert first["partners"] == []
    expect(page.locator("#model-index")).to_have_value("")
    page.locator("#model-index").select_option("0")
    page.locator("#assembly-choice").select_option("deposited")
    page.get_by_label("Chain A", exact=True).check()
    expect(page.locator(".partner-choice")).to_have_count(0)
    with page.expect_request(lambda request: request.url.endswith("/api/commands") and
                             request.post_data_json.get("kind") == "selectProteinModel") as submitted:
        page.get_by_role("button", name="Assess selected protein").click()
    selection = submitted.value.post_data_json["data"]
    assert selection["modelIndex"] == 0 and selection["biologicalAssemblyId"] is None
    assert selection["chains"] == [{"sourceChain": "A", "copyId": "A"}]
    assert selection["partners"] == []

    planned = await_state(page, lambda state: (state.get("preparationPlan") or {}).get(
        "standing") in ("ready", "partial", "failed"),
        "1L2Y model-1 recommendations in the connected route", timeout=300)
    plan = planned["preparationPlan"]
    assert plan["standing"] == "ready" and plan["authorizationAvailable"], plan
    assert plan["method"] and plan["nominalPh"] == 7.0
    assert plan["heavyAtomCount"] == 0
    assert plan["removedSourceHydrogenCount"] == MODEL_HYDROGEN_COUNT
    assert not any(item["kind"] == "disulfide" for item in
                   planned["preparationReview"]["decisions"])
    expect(page.get_by_role("button", name="Prepare with recommendations"))\
        .to_be_enabled()
    page.get_by_role("button", name="Prepare with recommendations").click()
    prepared = await_state(page, lambda state: (state.get("proteinTask") or {}).get(
        "standing") in ("assessed", "failed", "unavailable"),
        "prepared 1L2Y model-1 outcome in the connected route", timeout=300)
    assert prepared["proteinTask"]["standing"] == "assessed", prepared["proteinTask"]
    assert prepared["preparationPlan"]["standing"] == "applied"
    assert prepared["protein"]["status"] == "assessed"
    assert prepared["study"]["modelIndex"] == 0
    assert prepared["study"]["biologicalAssemblyId"] is None
    assert prepared["study"]["chainIds"] == ["A"] and prepared["study"]["partners"] == []
    return prepared


class ConnectedMemgenTrpCageRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert SOURCE.stat().st_size == 959202
        assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
        assert len(selected_source_heavy_atom_ids()) == 154
        for path in (HOST, POLICY, ROOT / "out/host/wwwroot/index.html"):
            assert path.is_file(), path
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_00_centered_four_component_preflight(self):
        """Establish a centered four-component pose and route without Build."""
        self._preflight_centered_route("unequal-center", {
            ("upper", "POPC"): 0.8, ("upper", "CHL1"): 0.2,
            ("lower", "DLPE"): 0.7, ("lower", "DLPC"): 0.3,
        })

    def test_01_pure_dopc_centered_preflight(self):
        """Establish a centered pure-DOPC pose and route without Build."""
        self._preflight_centered_route("pure-dopc", {
            ("upper", "DOPC"): 1.0, ("lower", "DOPC"): 1.0,
        })

    def _preflight_centered_route(self, route_key: str,
                                  expected: dict[tuple[str, str], float]):
        folder = Path(tempfile.mkdtemp(
            prefix=f"preflight-1l2y-{route_key}-", dir=ARTIFACTS))
        with patch.dict(os.environ, {"PIM_BROWSER_HOST": str(HOST),
                                     "PIM_AMBERTOOLS_HOME": str(ROOT / "out/ambertools26")}), \
                running_host_without_ppm(folder / "workspace", POLICY, ppm=None) as base, \
                sync_playwright() as playwright:
            browser = chromium(playwright)
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                command_kinds = []
                page.on("request", lambda request: command_kinds.append(
                    request.post_data_json["kind"]) if request.url.endswith("/api/commands") else None)
                page.goto(base, wait_until="domcontentloaded")
                prepared = prepare_exact_trpcage(page)
                adopted = choose_membrane_and_manual_pose(
                    page, prepared["protein"]["subjectId"], route_key)
                assert {(side, item["speciesId"]): item["fraction"]
                        for side in ("upper", "lower") for item in
                        adopted["membrane"][side]} == expected
                placement = adopted["placement"]
                assert placement["status"] == "supported", placement
                assert placement["preparedProteinId"] == prepared["protein"]["subjectId"]
                assert placement["membraneModelId"] == adopted["membrane"]["modelId"]
                assert placement["transform"]["startingPosition"] == "center"
                assert placement["transform"]["offsetXAngstrom"] == 0.0
                assert adopted["study"]["adoptedPlacementProposalId"] == placement["proposalId"]
                assert adopted["attempt"] is None and adopted["stages"] == []
                ppm = next(item for item in adopted["placementMethods"]
                           if item["method"] == "PPM")
                assert ppm["standing"] == "unavailable", ppm
                route = next(item for item in adopted["constructionRoutes"]
                             if item["policyId"] == ROUTE_ID)
                assert route["label"] == ROUTE_LABEL and route["available"], route
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Preparation", exact=True).click()
                recipe = page.locator(".construction-route").filter(has_text=ROUTE_LABEL)
                expect(recipe).to_have_count(1)
                expect(recipe.get_by_role("button", name="Build and minimize"))\
                    .to_be_enabled()
                page.screenshot(path=str(folder / f"centered-{route_key}-route-1024.png"),
                                full_page=True, animations="disabled")
                assert "buildAndMinimize" not in command_kinds, command_kinds
                final = current_state(page)
                assert final["attempt"] is None and final["stages"] == []
            finally:
                try:
                    with urlopen(base + "/api/state", timeout=10) as response:
                        (folder / "preflight-account.json").write_text(
                            json.dumps(json.load(response), indent=2) + "\n")
                except Exception as failure:
                    (folder / "account-read-error.txt").write_text(str(failure) + "\n")
                browser.close()

    def test_one_unequal_leaflet_build_to_export_readback(self):
        """One real 1L2Y attempt; the shared inspector reads the downloaded ZIP."""
        ConnectedMemgenCrambinRoute._run_connected_route(
            self, route_key="unequal", source_case="trpcage",
            prepare_source=prepare_exact_trpcage)

    def test_two_centered_unequal_leaflet_build_to_export_readback(self):
        """A second real 1L2Y attempt uses the centered four-component pose."""
        ConnectedMemgenCrambinRoute._run_connected_route(
            self, route_key="unequal-center", source_case="trpcage",
            prepare_source=prepare_exact_trpcage)

    def test_three_pure_dopc_build_to_export_readback(self):
        """A real 1L2Y pure-DOPC attempt uses the centered manual pose."""
        ConnectedMemgenCrambinRoute._run_connected_route(
            self, route_key="pure-dopc", source_case="trpcage",
            prepare_source=prepare_exact_trpcage)


if __name__ == "__main__":
    unittest.main()
