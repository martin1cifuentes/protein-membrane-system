"""Exact atom-identity helpers for the single connected browser acceptance route.

The real 6QWR/DMPC construction runs once in
ExportCompletedMinimizedStage/test_export.py. The controlled reattachment cases
remain in InspectAndReattach/test_browser_presentation.py.
"""

from __future__ import annotations

from urllib.parse import urlsplit

from playwright.sync_api import expect


def selected_atom(base: str, page, account: dict, index: int):
    inspection = account["inspection"]
    token = urlsplit(inspection["structureUrl"]).path.rsplit("/", 1)[1]
    response = page.request.get(base + f"/api/inspection/atoms/"
                                f"{inspection['subjectId']}/{token}/{index}")
    assert response.ok, f"Unverified atom row {index}: {response.status} {response.text()}"
    value = response.json()
    assert value["subjectId"] == inspection["subjectId"]
    assert value["studyRevisionId"] == inspection["studyRevisionId"]
    assert value["structureToken"] == token and value["atomSiteIndex"] == index
    assert value["atom"]["resultAtomIndex"] == index
    return value["atom"]


def atom_site_elements(mmcif: bytes):
    lines = mmcif.decode("utf-8").splitlines()
    headers = []
    for line_index, line in enumerate(lines):
        if line.startswith("_atom_site."):
            headers.append(line.strip())
            continue
        if headers and line and not line.startswith("_atom_site."):
            if "_atom_site.type_symbol" not in headers:
                raise AssertionError("Inspected mmCIF has no atom-site element column")
            element_column = headers.index("_atom_site.type_symbol")
            elements = []
            for atom_line in lines[line_index:]:
                if atom_line.startswith(("#", "loop_")):
                    break
                if atom_line.startswith(("ATOM", "HETATM")):
                    elements.append(atom_line.split()[element_column])
            if elements:
                return elements
    raise AssertionError("Inspected mmCIF has no atom-site rows")


def click_molstar_atom(page, expected: dict):
    index = expected["resultAtomIndex"]
    emitted = page.locator(".viewer-mount").evaluate("""(mount, atomSiteIndex) => {
      const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
      const structure = viewer?.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!structure) return false;
      for (const unit of structure.units) {
        if (!unit.model?.atomicHierarchy?.atomSourceIndex) continue;
        for (let offset = 0; offset < unit.elements.length; offset++) {
          const modelElement = unit.elements[offset];
          if (unit.model.atomicHierarchy.atomSourceIndex.value(modelElement) !== atomSiteIndex) continue;
          const loci = {kind: 'element-loci', structure,
            elements: [{unit, indices: Int32Array.of(offset)}]};
          viewer.plugin.behaviors.interaction.click.next({
            current: {loci, repr: null}, buttons: 1, button: 1,
            modifiers: {alt: false, control: false, meta: false, shift: false},
          });
          return true;
        }
      }
      return false;
    }""", index)
    assert emitted, f"Mol* did not contain expected atom-site row {index}"
    card = page.get_by_role("region", name="Inspection selection")
    expect(card).to_contain_text(expected["resultAtomId"], timeout=30000)
    expect(card).to_contain_text(expected["moleculeRole"])
    if expected["generatedSpeciesId"]:
        expect(card).to_contain_text(expected["generatedSpeciesId"])
    if expected["physicalSide"]:
        expect(card).to_contain_text(expected["physicalSide"])
    if expected["sourceAtomId"]:
        expect(card).to_contain_text(expected["sourceAtomId"])
