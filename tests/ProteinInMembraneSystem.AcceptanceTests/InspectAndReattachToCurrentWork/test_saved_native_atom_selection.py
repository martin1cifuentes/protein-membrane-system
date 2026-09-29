"""Controlled Mol* selection replay of one immutable, real minimized native stage.

The original Host closed after exporting this stage. This test serves its saved
account, CIF and correspondence through a controlled HTTP fixture; it does not
recover the Host, run construction, or establish a new scientific result. The
connected native route separately checked these atom rows against its live Host.
"""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re
import sys
import threading
import unittest
from urllib.parse import urlsplit
import zipfile

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" /
                       "InspectAndReattach"))
from test_browser_presentation import AccountHandler, AccountServer, CHROMIUM  # noqa: E402
from test_reattach import atom_site_elements, click_molstar_atom  # noqa: E402


SAVED = ROOT / "out" / "browser-acceptance" / "connected-native" / "real-native-65fk_lkd"
ATTEMPT_ID = "beb5f13ee4db4a78848caea311174ed9"
STAGE_ID = "3484b1b930cb4299bc5938735fb80d9f"
REVISION_ID = "c934dcf2efb2449783c3567f31dc4639"
ASSESSMENT_ID = "da312a3da0ae40e88ab6b9677cc4b96a"
STRUCTURE_TOKEN = "4a846749d7f84e319b96ecf177667629"
PINNED_SHA256 = {
    "final-account.json": "bdec54a0bbbdf90297a4a3365a81ba7f7d74e2f75b1309e68e3bb8bad77cb4b6",
    "minimized-coordinates.cif": "57c7edccbf833b7b387799fe1176957ffbcf87a71df417a7b5cae71cffe788b1",
    "constructed-correspondence.json": "1e1f1bc87b1a5a9006e3d98deb462f8997e2f581b77e0304043e0c5a96c852c1",
    "downloaded-stage.zip": "48e075b5839bca9b36de410da6eb450aa98e39ea4b7e0e735d8e26b07ca25103",
}
# Distinguishable real rows, selected before the controlled browser starts.
SELECTED_ROWS = (
    ("protein", 0, "protein", None, None),
    ("upper DMPC", 3205, "lipid", "upper", "DMPC"),
    ("lower DMPC", 9223, "lipid", "lower", "DMPC"),
    ("lower HOH", 29637, "water", "lower", "HOH"),
    ("lower Na", 85557, "ion", "lower", "NA"),
    ("upper Cl", 85616, "ion", "upper", "CL"),
)


def verified_saved_stage() -> tuple[dict, dict, bytes]:
    paths = {
        "final-account.json": SAVED / "final-account.json",
        "minimized-coordinates.cif": SAVED / "workspace" / "attempts" / ATTEMPT_ID /
            "minimized-coordinates.cif",
        "constructed-correspondence.json": SAVED / "workspace" / "attempts" / ATTEMPT_ID /
            "constructed-correspondence.json",
        "downloaded-stage.zip": SAVED / "downloaded-stage.zip",
    }
    for name, path in paths.items():
        assert path.is_file(), f"Pinned saved-stage artifact unavailable: {path}"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == PINNED_SHA256[name], name
    account = json.loads(paths["final-account.json"].read_text())
    mapping = json.loads(paths["constructed-correspondence.json"].read_text())
    cif = paths["minimized-coordinates.cif"].read_bytes()
    stage = next(item for item in account["stages"] if item["stageId"] == STAGE_ID)
    assert account["attempt"]["attemptId"] == ATTEMPT_ID
    assert account["study"]["id"] == REVISION_ID
    assert account["inspection"]["subjectId"] == STAGE_ID
    assert account["inspection"]["studyRevisionId"] == REVISION_ID
    assert account["inspection"]["structureUrl"] == \
        f"/api/structures/{STRUCTURE_TOKEN}?format=mmcif"
    assert account["inspection"]["representationKind"] == "completedStage"
    assert stage["attemptId"] == ATTEMPT_ID and stage["studyRevisionId"] == REVISION_ID
    assert stage["kind"] == "Minimization" and stage["status"] == "completed"
    assert stage["assessment"]["id"] == ASSESSMENT_ID
    assert stage["assessment"]["checkStanding"] == "checksPassed"
    assert stage["constructed"]["atomCount"] == 85667
    assert mapping["complete"] and len(mapping["atoms"]) == 85667
    assert mapping["sourceId"] == "eb2064ca8b53b56ab5c4bda5eec3e398434b72c610b4171990212ed051e3e628"
    assert mapping["resultId"] == "a64a105b7b30d0f90607d7ffd3065703b78cc432241b075786db9162a5eecdd0"
    elements = atom_site_elements(cif)
    assert len(elements) == len(mapping["atoms"])
    with zipfile.ZipFile(paths["downloaded-stage.zip"]) as bundle:
        manifest = json.loads(bundle.read("manifest.json"))
        assert hashlib.sha256(bundle.read("structure.cif")).hexdigest() == PINNED_SHA256[
            "minimized-coordinates.cif"]
        assert manifest["stage"]["id"] == STAGE_ID
        assert manifest["attempt"]["id"] == ATTEMPT_ID
        assert manifest["study"]["id"] == REVISION_ID
        assert manifest["assessment"]["id"] == ASSESSMENT_ID
        assert manifest["correspondence"]["atoms"] == mapping["atoms"]
    for label, index, role, side, species in SELECTED_ROWS:
        atom = mapping["atoms"][index]
        assert atom["resultAtomIndex"] == index, label
        assert atom["moleculeRole"] == role, label
        assert atom.get("physicalSide") == side, label
        assert atom.get("generatedSpeciesId") == species, label
        assert elements[index] == atom["element"], label
    return account, mapping, cif


class SavedStageHandler(AccountHandler):
    def do_GET(self):
        path = urlsplit(self.path).path
        if path.startswith("/api/structures/"):
            if path != f"/api/structures/{STRUCTURE_TOKEN}":
                self.send_error(404)
                return
            body = self.server.structure_bytes
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-mmcif")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        if path.startswith("/api/inspection/atoms/"):
            parts = path.split("/")
            if (len(parts) != 7 or parts[4:6] != [STAGE_ID, STRUCTURE_TOKEN] or
                    not parts[6].isdecimal() or int(parts[6]) >= len(self.server.saved_atoms)):
                self.send_error(404)
                return
            index = int(parts[6])
            self.server.atom_requests.append((STAGE_ID, STRUCTURE_TOKEN, index))
            self.json_response({"subjectId": STAGE_ID, "studyRevisionId": REVISION_ID,
                                "structureToken": STRUCTURE_TOKEN, "atomSiteIndex": index,
                                "atom": self.server.saved_atoms[index]})
            return
        super().do_GET()


@contextmanager
def saved_stage_server(account: dict, mapping: dict, cif: bytes):
    server = AccountServer()
    server.RequestHandlerClass = SavedStageHandler
    server.account = account
    server.structure_bytes = cif
    server.saved_atoms = mapping["atoms"]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server, f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


class SavedNativeAtomSelectionTests(unittest.TestCase):
    def test_identifies_lower_lipid_water_and_both_ions_in_saved_real_stage(self):
        account, mapping, cif = verified_saved_stage()
        with saved_stage_server(account, mapping, cif) as (host, base), sync_playwright() as playwright:
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".workspace.execution-review")).to_have_count(1, timeout=30000)
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
                expect(page.locator(".scene-error")).to_have_count(0)
                page.wait_for_function("""() => {
                  const mount = document.querySelector('.viewer-mount');
                  const viewer = mount && Reflect.get(mount, Symbol.for('molstar.viewer'));
                  const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
                  return structure && structure.elementCount === 85667;
                }""", timeout=240000)
                expect(page.locator(".stage-strip")).to_contain_text("Minimized")
                page.get_by_role("button", name=re.compile("Components")).click()
                components = page.get_by_label("Visible molecular components")
                components.get_by_role("radio", name="All").check()
                expect(page.locator(".execution-display-disclosure")).to_contain_text("water all shown")
                expect(page.locator(".execution-display-disclosure")).to_contain_text("ions shown")
                expect(components.get_by_role("checkbox", name=re.compile("Lipids"))).to_be_checked()
                expect(components.get_by_role("checkbox", name=re.compile("Ions"))).to_be_checked()
                page.wait_for_function("""() => {
                  const mount = document.querySelector('.viewer-mount');
                  const viewer = mount && Reflect.get(mount, Symbol.for('molstar.viewer'));
                  const components = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.components;
                  return components?.some(component =>
                    component.cell.transform.tags?.includes('structure-component-static-water') &&
                    !component.cell.state.isHidden);
                }""", timeout=30000)
                page.get_by_label("Protein representation").select_option("sticks")
                for label, index, _, _, _ in SELECTED_ROWS:
                    with self.subTest(label=label):
                        click_molstar_atom(page, mapping["atoms"][index])
                        selected = page.get_by_role("region", name="Inspection selection")
                        expect(selected).to_contain_text(mapping["atoms"][index]["resultAtomId"])
                        expect(selected).to_contain_text(mapping["atoms"][index]["moleculeRole"])
                        if label == "lower HOH":
                            # All water was visible for its Mol* selection; restore
                            # the lighter overview before checking the two ions.
                            components.get_by_role("radio", name="Hidden").check()
                self.assertEqual([row for _, _, row in host.atom_requests],
                                 [index for _, index, _, _, _ in SELECTED_ROWS])
                self.assertEqual(host.commands, [])
                self.assertEqual(host.account["study"]["id"], REVISION_ID)
                self.assertEqual(host.account["stages"][0]["assessment"]["id"], ASSESSMENT_ID)
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
