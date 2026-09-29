"""Published-browser preflight for the independently sourced RCSB 1L2Y ensemble.

The fixture is the exact RCSB PDB download (CC0),
https://files.rcsb.org/download/1L2Y.pdb, SHA256
5d1bbb545a312dfff1ae1e64b6d8addecb2f561ddc4011aeb5bee9d1dfcd4438.
It contains 38 deposited solution-NMR coordinate models. This test explicitly
chooses model 1 (account index 0), deposited chain A, and an upper-side manual
pose. It checks only prerequisites and route availability. It never authorizes
construction or claims a Memgen result.
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
    ARTIFACTS, HOST, POLICY, ROOT, ROUTE_ID, ROUTE_LABEL, await_state,
    choose_membrane_and_manual_pose, chromium, current_state,
    running_host_without_ppm,
)


SOURCE = Path(__file__).parent / "fixtures" / "1L2Y-RCSB.pdb"
SOURCE_SHA256 = "5d1bbb545a312dfff1ae1e64b6d8addecb2f561ddc4011aeb5bee9d1dfcd4438"
MODEL_COUNT = 38
MODEL_ATOM_COUNT = 304
MODEL_HEAVY_ATOM_COUNT = 154
MODEL_HYDROGEN_COUNT = 150
PROVENANCE_NOTE = "RCSB PDB 1L2Y; deposited solution NMR ensemble, 38 models"


def selected_source_heavy_atom_ids() -> set[str]:
    """Derive expected identities from deposited model 1, not worker output."""
    model_lines = []
    in_first_model = False
    for line in SOURCE.read_text().splitlines():
        if line.startswith("MODEL "):
            if in_first_model:
                break
            in_first_model = True
        elif in_first_model and line.startswith("ENDMDL"):
            break
        elif in_first_model and line.startswith("ATOM  "):
            model_lines.append(line)
    assert len(model_lines) == MODEL_ATOM_COUNT
    assert all(line[21] == "A" and line[16] == " " for line in model_lines)
    assert len({(line[22:27], line[12:16]) for line in model_lines}) == MODEL_ATOM_COUNT
    heavy = [line for line in model_lines if line[76:78].strip() not in ("H", "D")]
    assert len(heavy) == MODEL_HEAVY_ATOM_COUNT
    assert len(model_lines) - len(heavy) == MODEL_HYDROGEN_COUNT
    return {f"0:A:A:{int(line[22:26])}:{line[26].strip()}:{line[12:16].strip()}"
            for line in heavy}


class ConnectedMemgenTrpCagePreflight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
        lines = SOURCE.read_text().splitlines()
        assert sum(line.startswith("MODEL ") for line in lines) == MODEL_COUNT
        assert sum(line.startswith("ENDMDL") for line in lines) == MODEL_COUNT
        assert any(line.startswith("EXPDTA") and "SOLUTION NMR" in line for line in lines)
        assert not any(line.startswith(("HETATM", "SSBOND")) for line in lines)
        selected_source_heavy_atom_ids()
        for path in (HOST, POLICY, ROOT / "out/host/wwwroot/index.html"):
            assert path.is_file(), path
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_00_exact_nmr_preparation_pose_and_named_route_without_build(self):
        folder = Path(tempfile.mkdtemp(prefix="preflight-1l2y-", dir=ARTIFACTS))
        workspace = folder / "workspace"
        with patch.dict(os.environ, {"PIM_BROWSER_HOST": str(HOST),
                                     "PIM_AMBERTOOLS_HOME": str(ROOT / "out/ambertools26")}), \
                running_host_without_ppm(workspace, POLICY, ppm=None) as base, \
                sync_playwright() as playwright:
            browser = chromium(playwright)
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                command_kinds = []
                page.on("request", lambda request: command_kinds.append(
                    request.post_data_json["kind"]) if request.url.endswith("/api/commands") else None)
                page.goto(base, wait_until="domcontentloaded")
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
                assert [model["sourceModelId"] for model in source["sourceModels"]] == [
                    str(number) for number in range(1, MODEL_COUNT + 1)]
                assert all(model["atomCount"] == MODEL_ATOM_COUNT and model["partners"] == []
                           for model in source["sourceModels"])
                first = source["sourceModels"][0]
                assert first["chains"] == [{"name": "A", "residueCount": 20,
                                             "atomCount": MODEL_ATOM_COUNT}]
                assert len(first["residues"]) == 20
                expect(page.locator("#model-index")).to_have_value("")
                page.locator("#model-index").select_option("0")
                expect(page.locator("#model-index")).to_have_value("0")
                page.locator("#assembly-choice").select_option("deposited")
                page.get_by_label("Chain A", exact=True).check()
                expect(page.locator(".partner-choice")).to_have_count(0)
                with page.expect_request(lambda request: request.url.endswith("/api/commands") and
                                         request.post_data_json.get("kind") == "selectProteinModel") as submitted:
                    page.get_by_role("button", name="Assess selected protein").click()
                selection = submitted.value.post_data_json["data"]
                assert selection["modelIndex"] == 0
                assert selection["biologicalAssemblyId"] is None
                assert selection["chains"] == [{"sourceChain": "A", "copyId": "A"}]
                assert selection["partners"] == []

                planned = await_state(page, lambda state: (state.get("preparationPlan") or {}).get(
                    "standing") in ("ready", "partial", "failed"),
                    "1L2Y model-1 preparation recommendations", timeout=300)
                plan = planned["preparationPlan"]
                assert plan["standing"] == "ready", plan
                assert plan["authorizationAvailable"] and plan["method"]
                assert plan["nominalPh"] == 7.0
                assert plan["removedSourceHydrogenCount"] == MODEL_HYDROGEN_COUNT
                assert not any(item["kind"] == "disulfide" for item in
                               planned["preparationReview"]["decisions"])
                expect(page.get_by_role("button", name="Prepare with recommendations"))\
                    .to_be_enabled()
                page.get_by_role("button", name="Prepare with recommendations").click()
                prepared = await_state(page, lambda state: (state.get("proteinTask") or {}).get(
                    "standing") in ("assessed", "failed", "unavailable"),
                    "prepared 1L2Y model-1 outcome", timeout=300)
                assert prepared["proteinTask"]["standing"] == "assessed", prepared["proteinTask"]
                assert prepared["preparationPlan"]["standing"] == "applied"
                assert prepared["protein"]["status"] == "assessed"
                assert prepared["study"]["modelIndex"] == 0
                assert prepared["study"]["chainIds"] == ["A"]
                assert prepared["study"]["partners"] == []

                correspondence_paths = list(workspace.rglob("protein-correspondence.json"))
                assert len(correspondence_paths) == 1, correspondence_paths
                correspondence = json.loads(correspondence_paths[0].read_text())
                assert correspondence["complete"]
                assert correspondence["sourceId"] == SOURCE_SHA256
                assert {atom["sourceAtomId"] for atom in correspondence["atoms"]
                        if atom["sourceAtomId"]} == selected_source_heavy_atom_ids()

                adopted = choose_membrane_and_manual_pose(
                    page, prepared["protein"]["subjectId"], "unequal")
                assert adopted["placement"]["transform"]["startingPosition"] == "upper"
                assert adopted["study"]["adoptedPlacementProposalId"] == \
                    adopted["placement"]["proposalId"]
                assert adopted["attempt"] is None and adopted["stages"] == []
                ppm = next(item for item in adopted["placementMethods"]
                           if item["method"] == "PPM")
                assert ppm["standing"] == "unavailable", ppm
                route = next(item for item in adopted["constructionRoutes"]
                             if item["policyId"] == ROUTE_ID)
                assert route["label"] == ROUTE_LABEL and route["available"], route
                page.screenshot(path=str(folder / "adopted-manual-1024.png"),
                                full_page=True, animations="disabled")
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Preparation", exact=True).click()
                recipe = page.locator(".construction-route").filter(has_text=ROUTE_LABEL)
                expect(recipe).to_have_count(1)
                expect(recipe.get_by_role("button", name="Build and minimize"))\
                    .to_be_enabled()
                page.screenshot(path=str(folder / "named-route-1024.png"),
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


if __name__ == "__main__":
    unittest.main()
