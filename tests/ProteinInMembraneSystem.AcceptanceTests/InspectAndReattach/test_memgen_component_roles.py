"""Local evidence replay: browser roles from saved real Memgen coordinates.

Build browser/ first. This controlled local account verifies rendering, not provider
construction; the coordinate and correspondence artifacts came from one real trial.
"""
from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
import json
from pathlib import Path
import threading
import unittest
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

from test_browser_presentation import AccountHandler, AccountServer, CHROMIUM, DIST, ROOT, account, completed, inspection


ARTIFACT = ROOT / "out/browser-acceptance/connected-memgen/real-memgen-trpcage-pure-dopc-wq3rtckx/workspace"
TRIAL = ARTIFACT / "attempts/639ffcd329d747c7b68cdbcfa90b850a/trials/feb8a4275f9a4e7a90127f8dfae0ceda"
CONSTRUCTED_CIF = TRIAL / "constructed-feb8a4275f9a4e7a90127f8dfae0ceda-topology.cif"
STAGE_CIF = ARTIFACT / "attempts/639ffcd329d747c7b68cdbcfa90b850a/minimized-coordinates.cif"
CORRESPONDENCE = TRIAL / "constructed-feb8a4275f9a4e7a90127f8dfae0ceda-correspondence.json"
PORTABLE_CIF = Path(__file__).with_name("portable-one-entity.cif")


def role_runs() -> tuple[int, list[dict], dict[str, int]]:
    atoms = json.loads(CORRESPONDENCE.read_text())["atoms"]
    runs: list[dict] = []
    counts: dict[str, int] = {}
    for index, atom in enumerate(atoms):
        assert atom["resultAtomIndex"] == index
        role = atom["moleculeRole"]
        if role not in {"protein", "lipid", "water", "ion"}:
            role = "retainedPartner"
        source = atom.get("sourceResidue")
        source_chain = source.get("chain") if source else None
        copy_id = source.get("copyId") if source else None
        counts[role] = counts.get(role, 0) + 1
        if runs and all((runs[-1][key] == value for key, value in
                         (("role", role), ("sourceChain", source_chain), ("copyId", copy_id)))):
            runs[-1]["endExclusive"] = index + 1
        else:
            runs.append({"start": index, "endExclusive": index + 1,
                         "role": role, "sourceChain": source_chain, "copyId": copy_id})
    return len(atoms), runs, counts


def portable_atom(index: int, run: dict) -> dict:
    role = run["role"]
    source = run["sourceChain"] is not None
    source_residue = ({"model": 0, "chain": run["sourceChain"],
                       "residue": 1 if index < 4 else 2, "insertionCode": "",
                       "copyId": run["copyId"]} if source else None)
    return {"resultAtomIndex": index, "resultAtomId": f"portable:{index}",
            "sourceAtomId": f"source:{index}" if source else None,
            "role": "source" if source else "generated", "moleculeRole": role,
            "atomRole": "mapped", "element": "C", "sourceResidue": source_residue,
            "approvedChangeId": None, "physicalSide": "upper" if role == "lipid" else None,
            "generatedSpeciesId": {"lipid": "LIP", "water": "HOH", "ion": "NA"}.get(role),
            "generatedComponentRole": role if not source else None}


class MemgenHandler(AccountHandler):
    def do_GET(self):
        path = urlsplit(self.path).path
        if path in {"/api/structures/constructed", "/api/structures/stage"}:
            body = self.server.structures_by_token[path.rsplit("/", 1)[1]]
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-mmcif")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        if path.startswith("/api/inspection/components/"):
            parts = path.split("/")
            with self.server.lock:
                selected = deepcopy(self.server.account["inspection"])
            if (len(parts) != 6 or parts[4] != selected["subjectId"] or
                    parts[5] != ("constructed" if selected["representationKind"] == "constructedSystem" else "stage")):
                self.send_error(404)
                return
            self.json_response({"subjectId": selected["subjectId"],
                                "studyRevisionId": selected["studyRevisionId"],
                                "structureToken": parts[5], "atomCount": self.server.atom_count,
                                "runs": self.server.runs})
            return
        if path.startswith("/api/inspection/atoms/"):
            parts = path.split("/")
            with self.server.lock:
                selected = deepcopy(self.server.account["inspection"])
            if len(parts) != 7 or parts[4] != selected["subjectId"] or not parts[6].isdigit():
                self.send_error(404)
                return
            row = int(parts[6])
            if row >= self.server.atom_count:
                self.send_error(404)
                return
            atom = self.server.atoms[row]
            with self.server.lock:
                self.server.atom_requests.append((parts[4], parts[5], row))
            self.json_response({"subjectId": parts[4], "studyRevisionId": selected["studyRevisionId"],
                                "structureToken": parts[5], "atomSiteIndex": row, "atom": atom})
            return
        super().do_GET()

    def do_POST(self):
        if not self.server.portable:
            super().do_POST()
            return
        if urlsplit(self.path).path != "/api/commands":
            self.send_error(404)
            return
        command = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        target = command.get("data", {}).get("subjectId")
        with self.server.lock:
            current = self.server.account
            self.server.commands.append(command.get("kind"))
            if (command.get("kind") != "selectInspectionSubject" or
                    command.get("expectedRevision") != current["revision"] or
                    target not in {"stage-one", "constructed-one"}):
                self.json_response({"reason": "Only the exact portable subjects are available."}, 422)
                return
            updated = deepcopy(current)
            updated["revision"] += 1
            if target == "stage-one":
                updated["inspection"] = inspection(target, "completedStage",
                    updated["stages"][0]["assessment"])
                token = "stage"
            else:
                updated["inspection"] = inspection(target, "constructedSystem")
                token = "constructed"
            updated["inspection"].update(structureUrl=f"/api/structures/{token}?format=mmcif",
                                         omittedMolecules=[])
            self.server.account = updated
            self.server.version += 1
        self.json_response(updated)


class MemgenServer(AccountServer):
    def __init__(self, portable: bool = False):
        self.portable = portable
        if portable:
            self.atom_count = 15
            self.runs = [
                {"start": 0, "endExclusive": 4, "role": "protein", "sourceChain": "A", "copyId": "A"},
                {"start": 4, "endExclusive": 8, "role": "protein", "sourceChain": "A", "copyId": "A-2"},
                {"start": 8, "endExclusive": 9, "role": "retainedPartner", "sourceChain": "A", "copyId": "A"},
                {"start": 9, "endExclusive": 11, "role": "lipid", "sourceChain": None, "copyId": None},
                {"start": 11, "endExclusive": 12, "role": "ion", "sourceChain": None, "copyId": None},
                {"start": 12, "endExclusive": 15, "role": "water", "sourceChain": None, "copyId": None},
            ]
            self.counts = {"protein": 8, "retainedPartner": 1, "lipid": 2, "ion": 1, "water": 3}
            self.atoms = [portable_atom(index, run) for run in self.runs
                          for index in range(run["start"], run["endExclusive"])]
            self.structures_by_token = {"constructed": PORTABLE_CIF.read_bytes(), "stage": PORTABLE_CIF.read_bytes()}
        else:
            self.atom_count, self.runs, self.counts = role_runs()
            self.atoms = json.loads(CORRESPONDENCE.read_text())["atoms"]
            self.structures_by_token = {"constructed": CONSTRUCTED_CIF.read_bytes(), "stage": STAGE_CIF.read_bytes()}
        super().__init__()
        self.RequestHandlerClass = MemgenHandler


@contextmanager
def memgen_server(portable: bool = False):
    server = MemgenServer(portable=portable)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server, f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def pick_coordinate_row(page, row: int) -> bool:
    return page.evaluate("""row => {
      const mount = document.querySelector('.viewer-mount');
      const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
      const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!structure) return false;
      for (const unit of structure.units) {
        for (let offset = 0; offset < unit.elements.length; offset++) {
          if (unit.model.atomicHierarchy.atomSourceIndex.value(unit.elements[offset]) !== row) continue;
          const loci = { kind: 'element-loci', structure,
            elements: [{ unit, indices: Int32Array.of(offset) }] };
          viewer.plugin.behaviors.interaction.click.next({
            current: { loci, repr: null }, buttons: 1, button: 1,
            modifiers: { alt: false, control: false, meta: false, shift: false }
          });
          return true;
        }
      }
      return false;
    }""", row)


class MemgenComponentRoleTests(unittest.TestCase):
    def test_portable_one_entity_role_components_and_disclosure(self):
        """Tracked 15-row mmCIF repeats the single label-entity classification hazard."""
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        coordinate_rows = [line.split() for line in PORTABLE_CIF.read_text().splitlines()
                           if line.startswith(("ATOM ", "HETATM "))]
        self.assertEqual(len(coordinate_rows), 15)
        self.assertEqual({row[7] for row in coordinate_rows}, {"1"},
                         "Every role must share one label entity in this controlled fixture")
        with memgen_server(portable=True) as (host, base), sync_playwright() as playwright:
            state = completed(account())
            state["attempt"]["constructed"].update(atomCount=15, waterCount=1,
                                                     sodiumCount=1, chlorideCount=0)
            state["stages"][0]["constructed"] = state["attempt"]["constructed"]
            state["inspection"].update(structureUrl="/api/structures/constructed?format=mmcif",
                                       omittedMolecules=[])
            host.replace(state)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1280, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                expect(page.locator(".scene-loading")).to_have_count(0, timeout=60000)
                expect(page.locator(".scene-error")).to_have_count(0)
                disclosure = page.locator(".execution-display-disclosure")
                expect(disclosure).to_be_visible(timeout=60000)
                observed = page.locator(".viewer-mount").evaluate("""mount => {
                  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                  const root = viewer?.plugin.managers.structure.hierarchy.current.structures[0];
                  return root?.components.map(component => ({
                    tags: component.cell.transform.tags ?? [],
                    count: component.cell.obj?.data.elementCount ?? 0,
                    hidden: component.cell.state.isHidden,
                  })) ?? [];
                }""")
                for role, count in host.counts.items():
                    selected = [item for item in observed
                                if f"structure-component-inspection-{role}" in item["tags"]]
                    self.assertEqual(len(selected), 1, (role, observed))
                    self.assertEqual(selected[0]["count"], count)
                presets = [item for item in observed
                           if not any(tag.startswith("structure-component-inspection-")
                                      for tag in item["tags"])]
                self.assertTrue(presets, "The fixture must exercise Mol* preset classification")
                self.assertTrue(any(item["count"] > max(host.counts.values()) for item in presets),
                                "At least one preset component must span multiple mapped roles")
                self.assertTrue(all(item["hidden"] for item in presets),
                                "Overlapping preset representations must be suppressed")
                expect(disclosure).to_contain_text("Protein shown")
                expect(disclosure).to_contain_text("retained partners shown")
                expect(disclosure).to_contain_text("lipids shown")
                expect(disclosure).to_contain_text("water sampled")
                expect(disclosure).to_contain_text("ions shown")
                expect(disclosure).to_contain_text("15 atoms")
                page.get_by_role("button", name="Components").click()
                options = page.get_by_label("Visible molecular components")
                expect(options).to_contain_text("Protein chains (A, A · copy A-2) · 8 atoms")
                expect(options).to_contain_text("Lipids · 2 atoms")
                expect(options).to_contain_text("Ions · 1 atoms")
                expect(options).to_contain_text("3 water atoms")
                for label, role, caption in (("Protein chains", "protein", "Protein hidden"),
                                             ("Ligands / retained partners", "retainedPartner", "retained partners hidden"),
                                             ("Lipids", "lipid", "lipids hidden"),
                                             ("Ions", "ion", "ions hidden")):
                    control = options.get_by_label(label)
                    control.uncheck()
                    expect(disclosure).to_contain_text(caption)
                    page.wait_for_function("""role => {
                      const root = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                        ?.plugin.managers.structure.hierarchy.current.structures[0];
                      return root?.components.find(c => c.cell.transform.tags?.includes(`structure-component-inspection-${role}`))
                        ?.cell.state.isHidden === true;
                    }""", arg=role)
                    control.check()
                    page.wait_for_function("""role => {
                      const root = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                        ?.plugin.managers.structure.hierarchy.current.structures[0];
                      return root?.components.find(c => c.cell.transform.tags?.includes(`structure-component-inspection-${role}`))
                        ?.cell.state.isHidden === false;
                    }""", arg=role)
                options.get_by_label("All", exact=True).check()
                expect(disclosure).to_contain_text("water all shown")
                options.get_by_label("Hidden", exact=True).check()
                expect(disclosure).to_contain_text("water hidden")
                options.get_by_label("Representative sample").check()
                expect(disclosure).to_contain_text("water sampled")
                self.assertEqual(host.commands, [], "Display controls must not issue scientific commands")
                areas = page.get_by_label("Research work areas")
                areas.get_by_role("button", name="Results").click()
                page.get_by_label("Completed stage choices").get_by_role("button").first.click()
                expect(page.locator(".viewer-mount")).to_have_attribute(
                    "aria-label", "3D structure for Minimized system", timeout=60000)
                areas.get_by_role("button", name="Preparation").click()
                page.get_by_role("button", name="Review attempt").click()
                expect(page.locator(".viewer-mount")).to_have_attribute(
                    "aria-label", "3D structure for Protein–membrane system", timeout=60000)
                expect(page.locator(".execution-display-disclosure")).to_contain_text("lipids shown")
                expect(page.get_by_label("Checked constructed candidate")).to_contain_text("15 atoms")
                self.assertEqual(host.account["inspection"]["subjectId"], "constructed-one")
            finally:
                browser.close()

    def test_real_memgen_constructed_and_stage_roles_are_exclusive_and_controllable(self):
        if not all(path.is_file() for path in (CONSTRUCTED_CIF, STAGE_CIF, CORRESPONDENCE)):
            self.skipTest("Local evidence replay requires ignored saved Memgen qualification artifacts")
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with memgen_server() as (host, base), sync_playwright() as playwright:
            self.assertEqual(host.atom_count, 23922)
            self.assertEqual(host.counts, {"protein": 304, "lipid": 12420, "ion": 17, "water": 11181})
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                for subject, token in (("constructed-one", "constructed"), ("stage-one", "stage")):
                    state = account() if token == "constructed" else completed(account())
                    state["attempt"]["constructed"].update(atomCount=host.atom_count,
                        waterCount=host.counts["water"] // 3, sodiumCount=8, chlorideCount=9)
                    if token == "stage":
                        state["stages"][0]["constructed"] = state["attempt"]["constructed"]
                        state["inspection"] = inspection("stage-one", "completedStage",
                                                          state["stages"][0]["assessment"])
                    state["inspection"].update(structureUrl=f"/api/structures/{token}?format=mmcif",
                                               omittedMolecules=[])
                    host.replace(state)
                    page = browser.new_page(viewport={"width": 1280, "height": 800})
                    try:
                        page.goto(base, wait_until="domcontentloaded")
                        expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=120000)
                        expect(page.locator(".scene-loading")).to_have_count(0, timeout=120000)
                        expect(page.locator(".scene-error")).to_have_count(0)
                        expect(page.locator(".execution-display-disclosure")).to_be_visible(timeout=120000)
                        observed = page.locator(".viewer-mount").evaluate("""mount => {
                          const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                          const root = viewer?.plugin.managers.structure.hierarchy.current.structures[0];
                          return root?.components.map(component => ({
                            tags: component.cell.transform.tags ?? [],
                            count: component.cell.obj?.data.elementCount ?? 0,
                            hidden: component.cell.state.isHidden,
                          })) ?? [];
                        }""")
                        for role, expected in host.counts.items():
                            rows = [item for item in observed if f"structure-component-inspection-{role}" in item["tags"]]
                            self.assertEqual(len(rows), 1, (subject, role, observed))
                            self.assertEqual(rows[0]["count"], expected, (subject, role))
                        self.assertTrue(all(item["hidden"] for item in observed
                                            if not any(tag.startswith("structure-component-inspection-")
                                                       for tag in item["tags"])), observed)
                        disclosure = page.locator(".execution-display-disclosure")
                        expect(disclosure).to_contain_text("Protein shown")
                        expect(disclosure).to_contain_text("lipids shown")
                        expect(disclosure).to_contain_text("water sampled")
                        expect(disclosure).to_contain_text("ions shown")
                        expect(disclosure).to_contain_text("23,922 atoms")
                        selection = page.get_by_label("Inspection selection")
                        pick_cases = [(0, "chain A", "copy A")]
                        pick_cases.extend((next(index for index, atom in enumerate(host.atoms)
                                                if atom["moleculeRole"] == "lipid" and atom["physicalSide"] == side),
                                           "DOPC", side)
                                          for side in ("lower", "upper"))
                        pick_cases.append((next(index for index, atom in enumerate(host.atoms)
                                                if atom["moleculeRole"] == "ion"), "ion", "Physical leaflet"))
                        for row, identity, detail in pick_cases:
                            self.assertTrue(pick_coordinate_row(page, row), f"Mol* omitted coordinate row {row}")
                            expect(selection).to_contain_text(identity)
                            if not selection.locator(".inspection-correspondence").evaluate("element => element.open"):
                                selection.get_by_text("Exact coordinate correspondence").click()
                            expect(selection.locator(".inspection-correspondence")).to_contain_text(
                                host.atoms[row]["resultAtomId"])
                            expect(selection.locator(".inspection-correspondence")).to_contain_text(detail)
                            self.assertIn((subject, token, row), host.atom_requests)
                        page.get_by_role("button", name="Components").click()
                        options = page.get_by_label("Visible molecular components")
                        expect(options).to_contain_text("Protein chains (A) · 304 atoms")
                        expect(options).to_contain_text("Lipids · 12,420 atoms")
                        expect(options).to_contain_text("Ions · 17 atoms")
                        expect(options).to_contain_text("11,181 water atoms")
                        options.get_by_label("Lipids").uncheck()
                        expect(disclosure).to_contain_text("lipids hidden")
                        page.wait_for_function("""() => {
                          const root = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                            ?.plugin.managers.structure.hierarchy.current.structures[0];
                          return root?.components.find(c => c.cell.transform.tags?.includes('structure-component-inspection-lipid'))
                            ?.cell.state.isHidden === true;
                        }""")
                        options.get_by_label("Ions").uncheck()
                        expect(disclosure).to_contain_text("ions hidden")
                        page.wait_for_function("""() => {
                          const root = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                            ?.plugin.managers.structure.hierarchy.current.structures[0];
                          return root?.components.find(c => c.cell.transform.tags?.includes('structure-component-inspection-ion'))
                            ?.cell.state.isHidden === true;
                        }""")
                        options.get_by_label("All", exact=True).check()
                        expect(disclosure).to_contain_text("water all shown")
                        page.wait_for_function("""() => {
                          const root = Reflect.get(document.querySelector('.viewer-mount'), Symbol.for('molstar.viewer'))
                            ?.plugin.managers.structure.hierarchy.current.structures[0];
                          const water = root?.components.find(c => c.cell.transform.tags?.includes('structure-component-inspection-water'));
                          const sample = root?.components.find(c => c.cell.transform.tags?.includes('structure-component-inspection-water-sample'));
                          return water?.cell.state.isHidden === false && sample?.cell.state.isHidden === true;
                        }""")
                        water_row = next(index for index, atom in enumerate(host.atoms)
                                         if atom["moleculeRole"] == "water")
                        self.assertTrue(pick_coordinate_row(page, water_row))
                        expect(selection).to_contain_text("water")
                        expect(selection.locator(".inspection-correspondence")).to_contain_text(
                            host.atoms[water_row]["resultAtomId"])
                        self.assertIn((subject, token, water_row), host.atom_requests)
                        options.get_by_label("Protein chains").uncheck()
                        expect(disclosure).to_contain_text("Protein hidden")
                        options.get_by_label("Protein chains").check()
                        options.get_by_label("Lipids").check()
                        options.get_by_label("Ions").check()
                        options.get_by_label("Representative sample").check()
                        expect(disclosure).to_contain_text("water sampled")
                        self.assertEqual(host.commands, [], "Display controls must not issue scientific commands")
                    finally:
                        page.close()
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
