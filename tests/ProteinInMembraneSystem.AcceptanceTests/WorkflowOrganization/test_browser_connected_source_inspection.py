"""Focused browser checks for source continuation and local display.

The controlled account and coordinates below are presentation fixtures. The
separate real-host upload test covers the worker and published intake route.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out" / "ui-screenshots" / "revised-inspection" / "fixtures"
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import (  # noqa: E402
    CHROMIUM, DIST, account, controlled_account_server, inspection,
)
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein"))
from test_browser_route import two_alanines  # noqa: E402


def source_account(subject: str | None = None, revision: int = 10) -> dict:
    value = account()
    value["attempt"] = None
    value["revision"] = revision
    value["study"].update(
        id=f"source-revision-{revision}", number=revision - 9,
        selectedSourceId=subject, selectedSourceKind="rcsb" if subject else None,
    )
    value["sourceCandidates"] = [
        {"id": item, "label": label, "kind": "rcsb",
         "provenance": "Controlled source fixture", "limitations": []}
        for item, label in (("rcsb:1AAA", "PDB 1AAA"), ("rcsb:2BBB", "PDB 2BBB"))
    ]
    value["actions"].append(
        {"kind": "selectSource", "subjectId": None, "enabled": True, "reason": None})
    if subject:
        selected = inspection(subject, "structuralSource")
        selected.update(
            studyRevisionId=value["study"]["id"],
            studyRevisionNumber=value["study"]["number"],
            structureUrl="/api/structures/tiny?format=pdb",
            evidence=[], findings=[], annotations=[], metrics=[], omittedMolecules=[],
        )
        value["inspection"] = selected
    else:
        value["inspection"] = None
    return value


def small_alanine_chain() -> bytes:
    """Eight joined fixture residues, solely to exercise polymer display."""
    original = [line for line in two_alanines(1, "A").splitlines()
                if line.startswith("ATOM")]
    lines = ["SEQRES   1 A    8  ALA ALA ALA ALA ALA ALA ALA ALA", "MODEL        1"]
    serial = 0
    for pair in range(4):
        for line in original:
            if line[12:16].strip() == "OXT" and pair < 3:
                continue
            serial += 1
            residue = int(line[22:26]) + pair * 2
            x = float(line[30:38]) + pair * 6.6
            lines.append(f"{line[:6]}{serial:5d}{line[11:22]}{residue:4d}{line[26:30]}"
                         f"{x:8.3f}{line[38:]}")
    return ("\n".join(lines + ["TER", "ENDMDL", "END", ""])).encode()


def small_chain_with_ligand() -> bytes:
    """Display fixture with an independently hideable non-protein component."""
    return small_alanine_chain().replace(
        b"TER\n", b"HETATM   42  C1  LIG B   1      10.000  10.000   3.000  1.00 20.00           C\nTER\n")


INTERCEPT_SOURCE = """
(() => {
  const nativeFetch = window.fetch.bind(window);
  window.__sourceCalls = [];
  window.fetch = (input, init) => {
    if (String(input) === '/api/commands' && init?.body &&
        JSON.parse(init.body).kind === 'selectSource') {
      return new Promise(resolve => {
        window.__sourceCalls.push({ data: JSON.parse(init.body), resolve });
      });
    }
    return nativeFetch(input, init);
  };
  window.__answerSource = (index, body, status) => {
    window.__sourceCalls[index].resolve(new Response(JSON.stringify(body), {
      status, headers: { 'Content-Type': 'application/json' }
    }));
  };
})();
"""


PICK_FIRST_ATOM = """() => {
  const mount = document.querySelector('.viewer-mount');
  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
  const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!structure) return false;
  const unit = structure.units.find(item => item.model?.atomicHierarchy?.atomSourceIndex && item.elements.length);
  if (!unit) return false;
  const loci = { kind: 'element-loci', structure,
    elements: [{ unit, indices: Int32Array.of(0) }] };
  viewer.plugin.behaviors.interaction.click.next({
    current: { loci, repr: null }, buttons: 1, button: 1,
    modifiers: { alt: false, control: false, meta: false, shift: false }
  });
  return true;
}"""


CAMERA_STATE = """() => {
  const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
  const camera = viewer.plugin.canvas3d.camera;
  const state = camera.getSnapshot();
  const cells = [...viewer.plugin.state.data.cells.values()];
  const focusSelections = cells.filter(cell => cell.transform.tags?.some(tag =>
    tag === 'structure-focus-target-sel' || tag === 'structure-focus-surr-sel'));
  return { target: [...state.target], position: [...state.position], radius: state.radius,
    near: camera.near, far: camera.far, viewport: camera.viewport,
    focusActive: !!viewer.plugin.managers.structure.focus.current,
    focusAtoms: focusSelections.map(cell => cell.obj?.data?.elementCount ?? 0),
    transparencyLayers: cells.filter(cell => cell.transform.tags?.includes('transparency-controls')).length };
}"""


def open_browser(playwright, width=1280, height=800):
    browser = playwright.chromium.launch(
        executable_path=CHROMIUM, headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"],
    )
    return browser, browser.new_page(viewport={"width": width, "height": height})


class ConnectedSourceInspectionBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert (DIST / "index.html").is_file(), "Build browser/ before these focused checks"
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_switching_proteins_replaces_the_old_scene_with_a_loading_screen(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(source_account("rcsb:1AAA", 11))
            browser, page = open_browser(playwright)
            try:
                page.add_init_script(INTERCEPT_SOURCE)
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".source-context")).to_contain_text(
                    "Structure displayed", timeout=60000)
                expect(page.locator(".viewer-mount")).to_be_visible()
                held_structure_requests = []
                page.route("**/api/structures/tiny*",
                           lambda route: held_structure_requests.append(route))

                page.locator(".source-explorer > summary").click()
                page.get_by_role("button", name="PDB 2BBB").click()
                expect(page.get_by_role("status", name="Loading source PDB 2BBB"))\
                    .to_contain_text("Retrieving source")
                expect(page.locator(".viewer-mount")).to_have_count(0)
                expect(page.locator(".scene-subtitle")).to_contain_text("PDB 2BBB")
                expect(page.locator(".scene-subtitle")).not_to_contain_text("rcsb:1AAA")
                expect(page.get_by_role("button", name="PDB 1AAA")).to_contain_text(
                    "Current study source")
                expect(page.locator(".source-context")).not_to_contain_text(
                    "Source coordinates are displayed")
                expect(page.locator(".source-evidence-pending")).to_contain_text(
                    "evidence below still belongs to Source structure · 1AAA")
                page.screenshot(path=str(ARTIFACTS / "protein-switch-loading.png"),
                                animations="disabled")

                new_source = source_account("rcsb:2BBB", 12)
                new_source["inspection"]["structureUrl"] += "&source=2BBB"
                page.evaluate("body => window.__answerSource(0, body, 200)", new_source)
                expect(page.get_by_role("status", name="Loading source PDB 2BBB"))\
                    .to_contain_text("Loading structure")
                expect(page.locator(".source-evidence-pending")).to_contain_text(
                    "Loading the coordinates for Source structure · 2BBB")
                for _ in range(100):
                    if held_structure_requests:
                        break
                    page.wait_for_timeout(100)
                self.assertTrue(held_structure_requests,
                                "The new coordinate response should be held during rendering")
                page.screenshot(path=str(ARTIFACTS / "protein-switch-rendering.png"),
                                animations="disabled")
                for request in held_structure_requests:
                    request.continue_()
                expect(page.locator(".source-context")).to_contain_text(
                    "Structure displayed", timeout=60000)
                expect(page.locator(".scene-source-loading")).to_have_count(0)
                expect(page.locator(".viewer-mount")).to_be_visible()
                expect(page.locator(".scene-heading h1")).to_have_text("Source structure · 2BBB")
                expect(page.locator(".selected-source-account")).to_contain_text("rcsb:2BBB")
                self.assertEqual(host.commands, [], "The response is fixture data")
            finally:
                browser.close()

    def test_newer_source_selection_wins_after_earlier_request_completes(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(source_account())
            browser, page = open_browser(playwright)
            try:
                page.add_init_script(INTERCEPT_SOURCE)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="PDB 1AAA").click()
                expect(page.locator(".source-context")).to_contain_text("Retrieving source")
                page.get_by_role("button", name="PDB 1AAA").click()
                self.assertEqual(page.evaluate("window.__sourceCalls.length"), 1,
                                 "Repeated selection must not submit the same source twice")
                page.get_by_role("button", name="PDB 2BBB").click()
                expect(page.locator(".source-context")).to_contain_text("PDB 2BBB")
                page.evaluate("body => window.__answerSource(0, body, 200)",
                              source_account("rcsb:1AAA", 11))
                page.wait_for_function("window.__sourceCalls.length === 2")
                expect(page.locator(".scene-heading h1")).to_have_text("Loading source structure")
                expect(page.locator(".viewer-mount")).to_have_count(0)
                expect(page.locator(".source-context")).to_contain_text("PDB 2BBB")
                page.evaluate("body => window.__answerSource(1, body, 200)",
                              source_account("rcsb:2BBB", 12))
                expect(page.locator(".source-context")).to_contain_text(
                    "Structure displayed", timeout=60000)
                expect(page.locator(".selected-source-account")).to_contain_text("rcsb:2BBB")
                expect(page.locator(".scene-heading h1")).to_have_text("Source structure · 2BBB")
                expect(page.get_by_role("button", name="Inspect source structure")).to_have_count(0)
                self.assertEqual(host.commands, [], "The delayed responses are fixture data")
            finally:
                browser.close()

    def test_retrieval_and_render_failures_keep_source_and_offer_the_right_retry(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(source_account())
            browser, page = open_browser(playwright)
            try:
                page.add_init_script(INTERCEPT_SOURCE)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="PDB 1AAA").click()
                page.screenshot(path=str(ARTIFACTS / "source-retrieving.png"), animations="disabled")
                page.evaluate("() => window.__answerSource(0, {reason: 'Source archive unavailable'}, 503)")
                expect(page.locator(".source-context")).to_contain_text("Retrieval failed")
                expect(page.locator(".source-context")).to_contain_text("Source archive unavailable")
                page.screenshot(path=str(ARTIFACTS / "source-retrieval-failed.png"), animations="disabled")
                page.route("**/api/structures/tiny*", lambda route: route.fulfill(
                    status=503, content_type="text/plain", body="Representation unavailable"))
                page.get_by_role("button", name="Retry retrieval").click()
                page.wait_for_function("window.__sourceCalls.length === 2")
                page.evaluate("body => window.__answerSource(1, body, 200)",
                              source_account("rcsb:1AAA", 11))
                expect(page.locator(".source-context")).to_contain_text(
                    "Visualization unavailable", timeout=60000)
                expect(page.locator(".scene-error")).to_contain_text("Structure unavailable")
                page.screenshot(path=str(ARTIFACTS / "source-render-failed.png"), animations="disabled")
                page.unroute("**/api/structures/tiny*")
                page.get_by_role("button", name="Retry visualization").click()
                expect(page.locator(".source-context")).to_contain_text(
                    "Structure displayed", timeout=60000)
                expect(page.locator(".selected-source-account")).to_contain_text("rcsb:1AAA")
                self.assertEqual(host.commands, [], "These responses are fixture data")
            finally:
                browser.close()

    def test_failed_newer_selection_identifies_the_earlier_established_view(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.replace(source_account())
            browser, page = open_browser(playwright)
            try:
                page.add_init_script(INTERCEPT_SOURCE)
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="PDB 1AAA").click()
                page.get_by_role("button", name="PDB 2BBB").click()
                page.evaluate("body => window.__answerSource(0, body, 200)",
                              source_account("rcsb:1AAA", 11))
                page.wait_for_function("window.__sourceCalls.length === 2")
                page.evaluate("() => window.__answerSource(1, {reason: 'Source archive unavailable'}, 503)")
                expect(page.locator(".source-context")).to_contain_text("PDB 2BBB")
                expect(page.locator(".source-context")).to_contain_text(
                    "View and evidence remain on rcsb:1AAA")
                expect(page.locator(".selected-source-account")).to_contain_text("rcsb:1AAA")
                expect(page.locator(".scene-heading h1")).to_have_text("Source structure · 1AAA")
                self.assertEqual(host.commands, [], "Delayed responses are fixture data")
            finally:
                browser.close()

    def test_local_selection_and_collapse_preserve_draft_subject_and_viewer(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.structure_bytes = small_alanine_chain()
            host.replace(source_account("rcsb:1AAA", 11))
            browser, page = open_browser(playwright, 1024, 720)
            try:
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".source-context")).to_contain_text(
                    "Structure displayed", timeout=60000)
                page.locator(".source-explorer > summary").click()
                page.locator("#source-query").fill("unfinished search")
                self.assertTrue(page.evaluate("() => { window.__viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')); return !!window.__viewer; }"))
                self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                expect(page.get_by_label("Inspection selection")).to_contain_text("chain A")
                self.assertEqual(host.atom_requests, [], "Raw source selection uses its displayed coordinates")
                expect(page.locator(".scene-local-state")).to_have_count(0)
                page.get_by_role("button", name="Inspect locally").click()
                expect(page.locator(".scene-local-state")).to_contain_text("5 Å")
                page.get_by_role("button", name="Components").click()
                expect(page.get_by_label("Visible molecular components")).to_contain_text(
                    "Visibility and rendering never change coordinates")
                page.get_by_role("button", name="Components").click()
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.locator(".rail")).to_be_hidden()
                expect(page.get_by_label("Evidence and assessment")).to_be_visible()
                expect(page.get_by_label("Inspection selection")).to_contain_text("chain A")
                self.assertTrue(page.evaluate("Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')) === window.__viewer"))
                page.get_by_role("button", name="Show inputs").click()
                expect(page.locator("#source-query")).to_have_value("unfinished search")
                expect(page.locator(".scene-local-state")).to_contain_text("5 Å")
                page.get_by_label("Protein representation").select_option("sticks")
                expect(page.get_by_label("Protein representation")).to_have_value("sticks")
                page.wait_for_function("""() => {
                  const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                  const polymer = viewer.plugin.managers.structure.hierarchy.current.structures[0]
                    ?.components.find(item => item.cell.transform.tags?.includes('structure-component-static-polymer'));
                  return polymer?.representations.some(item => item.cell.params?.values?.type?.name === 'ball-and-stick' &&
                    !item.cell.state.isHidden);
                }""")
                page.get_by_role("button", name="Components").click()
                protein_visibility = page.locator(".scene-component-options label")\
                    .filter(has_text="Protein chains").locator("input")
                protein_visibility.uncheck()
                self.assertTrue(page.evaluate("""() => {
                  const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                  return viewer.plugin.managers.structure.hierarchy.current.structures[0]
                    ?.components.find(item => item.cell.transform.tags?.includes('structure-component-static-polymer'))
                    ?.cell.state.isHidden === true;
                }"""))
                protein_visibility.check()
                page.get_by_role("button", name="Components").click()
                page.get_by_role("button", name="Show whole system").click()
                expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
                expect(page.locator(".scene-local-state")).to_have_count(0)
                self.assertEqual(host.commands, [], "Inspection controls cannot submit scientific work")
                self.assertLessEqual(page.evaluate("document.body.scrollHeight"), 720)
                page.screenshot(path=str(ARTIFACTS / "small-restored-source.png"), animations="disabled")
            finally:
                browser.close()

    def test_clear_selection_restores_camera_only_when_leaving_local_inspection(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            host.structure_bytes = small_chain_with_ligand()
            host.replace(source_account("rcsb:1AAA", 11))
            browser, page = open_browser(playwright, 1024, 720)
            try:
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=60000)
                page.wait_for_function("() => document.querySelector('.viewer-mount')?.dataset.cameraReady === 'true'")
                page.get_by_label("Protein representation").select_option("sticks")
                page.wait_for_function("""() => {
                  const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                  const polymer = viewer.plugin.managers.structure.hierarchy.current.structures[0]
                    ?.components.find(item => item.cell.transform.tags?.includes('structure-component-static-polymer'));
                  return polymer?.representations.some(item => item.cell.params?.values?.type?.name === 'ball-and-stick' &&
                    !item.cell.state.isHidden);
                }""")
                page.get_by_role("button", name="Components").click()
                ligand_visibility = page.locator(".scene-component-options label")\
                    .filter(has_text="Ligands").locator("input")
                expect(ligand_visibility).to_be_enabled()
                ligand_visibility.uncheck()
                page.get_by_role("button", name="Components").click()
                page.wait_for_timeout(250)

                # An ordinary pick and clear preserve a researcher-positioned camera.
                page.evaluate("""() => {
                  const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                  const camera = viewer.plugin.canvas3d.camera;
                  const state = camera.getSnapshot();
                  camera.setState({ target: state.target.map((value, axis) => value + (axis === 0 ? 3 : 0)),
                    position: state.position.map((value, axis) => value + (axis === 0 ? 3 : 0)) }, 0);
                }""")
                ordinary_camera = page.evaluate(CAMERA_STATE)
                self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                expect(page.get_by_label("Inspection selection")).to_contain_text("chain A")
                page.get_by_role("region", name="Molecular structure").get_by_role("button", name="Clear selection").click()
                expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
                after_ordinary = page.evaluate(CAMERA_STATE)
                for key in ("target", "position"):
                    for actual, expected in zip(after_ordinary[key], ordinary_camera[key]):
                        self.assertAlmostEqual(actual, expected, places=3)
                self.assertAlmostEqual(after_ordinary["radius"], ordinary_camera["radius"], places=3)

                page.get_by_role("button", name="Show whole system").click()
                page.wait_for_function("""shifted => {
                  const camera = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                    .plugin.canvas3d.camera;
                  return Math.abs(camera.state.target[0] - shifted) > 1;
                }""", arg=ordinary_camera["target"][0])
                overview = page.evaluate(CAMERA_STATE)

                for collapsed in (False, True):
                    if collapsed:
                        page.get_by_role("button", name="Collapse inputs").click()
                        expect(page.locator(".rail")).to_be_hidden()
                        page.wait_for_function("() => document.querySelector('.viewer-mount')?.dataset.cameraReady === 'true'")
                        overview = page.evaluate(CAMERA_STATE)
                    self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                    page.get_by_role("button", name="Inspect locally").click()
                    page.wait_for_function("""limits => {
                      const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                      const camera = viewer.plugin.canvas3d.camera;
                      return !!document.querySelector('.scene-local-state') &&
                        camera.state.radius < limits.radius && camera.far - camera.near < limits.depth;
                    }""", arg={"radius": overview["radius"] * 0.8,
                                "depth": (overview["far"] - overview["near"]) * 0.8})
                    local_camera = page.evaluate(CAMERA_STATE)
                    self.assertLess(local_camera["far"] - local_camera["near"],
                                    overview["far"] - overview["near"])
                    page.get_by_role("region", name="Molecular structure").get_by_role("button", name="Clear selection").click()
                    page.wait_for_function("""limits => {
                      const viewer = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'));
                      const cells = [...viewer.plugin.state.data.cells.values()];
                      return !document.querySelector('.scene-local-state') &&
                        viewer.plugin.canvas3d.camera.state.radius > limits.radius &&
                        viewer.plugin.canvas3d.camera.far - viewer.plugin.canvas3d.camera.near > limits.depth &&
                        !viewer.plugin.managers.structure.focus.current &&
                        cells.filter(cell => cell.transform.tags?.some(tag =>
                          tag === 'structure-focus-target-sel' || tag === 'structure-focus-surr-sel'))
                          .every(cell => !cell.obj?.data?.elementCount) &&
                        !cells.some(cell => cell.transform.tags?.includes('transparency-controls'));
                    }""", arg={"radius": local_camera["radius"] * 1.5,
                                "depth": (local_camera["far"] - local_camera["near"]) * 1.5})
                    cleared = page.evaluate(CAMERA_STATE)
                    self.assertAlmostEqual(cleared["radius"], overview["radius"], delta=1)
                    self.assertGreater(cleared["far"] - cleared["near"],
                                       local_camera["far"] - local_camera["near"])
                    self.assertEqual(cleared["focusAtoms"], [0, 0])
                    self.assertEqual(cleared["transparencyLayers"], 0)
                    self.assertAlmostEqual(cleared["viewport"]["width"],
                                           page.locator(".viewer-mount").bounding_box()["width"], delta=4)
                    expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
                    expect(page.get_by_label("Protein representation")).to_have_value("sticks")
                    page.get_by_role("button", name="Components").click()
                    expect(ligand_visibility).not_to_be_checked()
                    page.get_by_role("button", name="Components").click()
                    page.screenshot(path=str(ARTIFACTS / ("clear-collapsed.png" if collapsed else "clear-expanded.png")))

                page.set_viewport_size({"width": 820, "height": 760})
                page.wait_for_function("() => document.querySelector('.viewer-mount')?.dataset.cameraReady === 'true'")
                page.get_by_role("button", name="Show whole system").click()
                page.wait_for_timeout(300)
                small_overview = page.evaluate(CAMERA_STATE)
                self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                page.get_by_role("button", name="Inspect locally").click()
                page.wait_for_function("limit => Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')).plugin.canvas3d.camera.state.radius < limit", arg=small_overview["radius"] * 0.8)
                page.get_by_role("region", name="Molecular structure").get_by_role("button", name="Clear selection").click()
                page.wait_for_function("limit => Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')).plugin.canvas3d.camera.state.radius > limit", arg=small_overview["radius"] * 0.8)
                expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
                page.screenshot(path=str(ARTIFACTS / "clear-small-collapsed.png"))
                page.get_by_role("button", name="Show inputs").click()
                expect(page.locator(".rail")).to_be_visible()
                page.wait_for_function("() => document.querySelector('.viewer-mount')?.dataset.cameraReady === 'true'")
                page.get_by_role("button", name="Show whole system").click()
                page.wait_for_timeout(300)
                small_expanded = page.evaluate(CAMERA_STATE)
                self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                page.get_by_role("button", name="Inspect locally").click()
                page.wait_for_function("limit => Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')).plugin.canvas3d.camera.state.radius < limit", arg=small_expanded["radius"] * 0.8)
                page.get_by_role("region", name="Molecular structure").get_by_role("button", name="Clear selection").click()
                page.wait_for_function("limit => Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer')).plugin.canvas3d.camera.state.radius > limit", arg=small_expanded["radius"] * 0.8)
                expect(page.get_by_label("Inspection selection")).to_contain_text("No residue selected")
                page.screenshot(path=str(ARTIFACTS / "clear-small-expanded.png"))

                # A pending Mol* focus update must not revive a cleared overlay.
                self.assertTrue(page.evaluate(PICK_FIRST_ATOM))
                page.get_by_role("button", name="Inspect locally").click()
                page.get_by_role("region", name="Molecular structure").get_by_role("button", name="Clear selection").click()
                page.wait_for_timeout(500)
                rapid_clear = page.evaluate(CAMERA_STATE)
                self.assertAlmostEqual(rapid_clear["radius"], small_expanded["radius"], delta=1)
                self.assertFalse(rapid_clear["focusActive"])
                self.assertEqual(rapid_clear["focusAtoms"], [0, 0])
                self.assertEqual(rapid_clear["transparencyLayers"], 0)
                self.assertEqual(host.commands, [], "Display changes cannot submit scientific work")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
