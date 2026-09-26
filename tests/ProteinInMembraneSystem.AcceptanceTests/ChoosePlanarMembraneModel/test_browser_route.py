"""Slice 2 membrane choice through the published host and official browser route.

Build with scripts/build-local.sh, then run with the separate Playwright test venv:

    PIM_BROWSER_CHROMIUM=/path/to/chrome out/browser-test-python/bin/python \
      -m unittest discover \
      -s tests/ProteinInMembraneSystem.AcceptanceTests/ChoosePlanarMembraneModel \
      -p 'test_*.py' -v

PIM_MEMBRANE_POLICY_CATALOGUE may identify another qualified combined catalogue.
The default is config/policies/protein-membrane-slice2.json. Captures are retained in
out/browser-acceptance/slice2 for Movement D review.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import tempfile
import time
import unittest
from urllib.error import URLError
from urllib.request import urlopen

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
HOST = ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"
WORKER = ROOT / "out" / "python" / "bin" / "python"
LAUNCHER = ROOT / "scripts" / "start-local.sh"
POLICY = Path(os.environ.get(
    "PIM_MEMBRANE_POLICY_CATALOGUE",
    str(ROOT / "config" / "policies" / "protein-membrane-slice2.json"),
)).resolve()
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "slice2"
PURPOSE = "Compare a defined mixed asymmetric PC and PE bilayer for a placement study."


@contextmanager
def running_host(workspace: Path, catalogue: Path):
    """Use the same launcher, worker, catalogue and loopback origin as a researcher."""
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    environment = dict(
        os.environ,
        PIM_PORT=str(port),
        PIM_WORKSPACE_ROOT=str(workspace),
        PIM_POLICY_CATALOGUE=str(catalogue),
    )
    environment.pop("PIM_PPM_EXECUTABLE", None)
    log_path = workspace.parent / "membrane-host.log"
    with log_path.open("wb") as log:
        process = subprocess.Popen(
            [str(LAUNCHER)], cwd=ROOT, env=environment,
            stdout=log, stderr=subprocess.STDOUT,
        )
        try:
            for _ in range(150):
                if process.poll() is not None:
                    raise AssertionError("Local launcher exited during startup: " + log_path.read_text())
                try:
                    with urlopen(base + "/api/state", timeout=2) as response:
                        if response.status == 200:
                            break
                except URLError:
                    time.sleep(0.1)
            else:
                raise AssertionError("Local host did not become ready: " + log_path.read_text())
            yield base
        except BaseException:
            print("\nHost output tail:\n" + log_path.read_text()[-1500:])
            raise
        finally:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)


def chromium(playwright):
    executable = os.environ.get("PIM_BROWSER_CHROMIUM")
    return playwright.chromium.launch(
        executable_path=executable or None,
        headless=True,
        args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
              "--use-gl=angle", "--use-angle=swiftshader"],
    )


def current_state(page):
    return page.evaluate("async () => (await fetch('/api/state', {cache: 'no-store'})).json()")


def await_membrane_status(page, status: str):
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        state = current_state(page)
        membrane = state.get("membrane") or {}
        ready_status = membrane.get("status") == status
        ready_evidence = status != "assessed" or (membrane.get("policyId") and membrane.get("evidence"))
        ready_reason = status != "notEstablished" or membrane.get("reason")
        if ready_status and ready_evidence and ready_reason:
            break
        time.sleep(0.1)
    else:
        raise AssertionError(f"Timed out waiting for membrane {status}: {state.get('membrane')}")
    return state


def assert_no_horizontal_overflow(test: unittest.TestCase, page):
    overflowing = page.evaluate("""() => [...document.querySelectorAll(
        '.workspace, .workspace-header, .rail, .scene-panel, .evidence-panel, .evidence-content')]
      .filter(element => element.scrollWidth > element.clientWidth + 1)
      .map(element => ({region: element.className, width: element.clientWidth,
                        contentWidth: element.scrollWidth}))""")
    test.assertEqual(overflowing, [])


def assert_workspace_geometry(test: unittest.TestCase, page, width: int, height: int):
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    test.assertIsNotNone(scene)
    test.assertIsNotNone(evidence)
    test.assertIsNotNone(stages)
    test.assertGreaterEqual(scene["width"], 350 if width <= 820 else width * 0.55)
    test.assertGreaterEqual(evidence["width"], 300 if width <= 820 else width * 0.25)
    test.assertGreaterEqual(scene["height"], height * 0.55)
    test.assertGreater(evidence["x"], scene["x"])
    test.assertAlmostEqual(evidence["y"], scene["y"], delta=3)
    test.assertLessEqual(stages["y"] + stages["height"], height + 2)
    assert_no_horizontal_overflow(test, page)


def choose_fraction(page, side: str, index: int, species: str, percentage: str):
    page.get_by_label(f"{side} leaflet lipid {index}", exact=True).select_option(species)
    page.get_by_label(f"{side} leaflet percentage {index}", exact=True).fill(percentage)


def inspect_membrane(page):
    page.get_by_role("button", name="Inspect membrane").click()
    expect(page.locator(".workspace")).to_have_class(re.compile(r"\bworkflow-closed\b"))
    expect(page.locator(".evidence-panel")).to_be_visible()


def require_intended_preview(page):
    # A labeled spatial preview may be SVG, canvas or image. The accessible role
    # ensures the researcher can identify it without interpreting the pixels.
    expect(page.get_by_role(
        "img", name=re.compile(r"intended.*bilayer|bilayer.*intended", re.I),
    )).to_be_visible()


def assert_exact_leaflets(test: unittest.TestCase, panel_text: str,
                          upper: dict[str, int], lower: dict[str, int]):
    upper_match = re.search(r"Upper\s+(?:physical\s+)?leaflet", panel_text, re.I)
    lower_match = re.search(r"Lower\s+(?:physical\s+)?leaflet", panel_text, re.I)
    test.assertIsNotNone(upper_match, "Upper physical leaflet is absent from the evidence panel")
    test.assertIsNotNone(lower_match, "Lower physical leaflet is absent from the evidence panel")
    test.assertLess(upper_match.start(), lower_match.start())
    upper_text = panel_text[upper_match.end():lower_match.start()]
    lower_text = panel_text[lower_match.end():]
    for species, percent in upper.items():
        test.assertIn(species, upper_text)
        test.assertRegex(upper_text, rf"\b{percent}(?:\.0)?\s*%")
    for species, percent in lower.items():
        test.assertIn(species, lower_text)
        test.assertRegex(lower_text, rf"\b{percent}(?:\.0)?\s*%")


def capture_review(test: unittest.TestCase, page, name: str, width: int, height: int):
    page.set_viewport_size({"width": width, "height": height})
    page.locator(".evidence-content").evaluate("element => { element.scrollTop = 0; }")
    expect(page.locator(".scene-panel")).to_be_visible()
    expect(page.locator(".evidence-panel")).to_be_visible()
    expect(page.locator(".stage-strip")).to_be_visible()
    require_intended_preview(page)
    assert_workspace_geometry(test, page, width, height)
    page.screenshot(
        path=str(ARTIFACTS / f"{name}-{width}.png"),
        full_page=True, animations="disabled",
    )


def catalogue_with_absolute_asset_paths(source: Path):
    """Keep a test-local catalogue pointed at the original, identified assets."""
    catalogue = json.loads(source.read_text())
    source_base = source.parent

    def absolute_asset_paths(value):
        if isinstance(value, list):
            return [absolute_asset_paths(item) for item in value]
        if isinstance(value, dict):
            return {
                key: str((source_base / item).resolve())
                if key in {"path", "templatePath", "coordinateTemplatePath"}
                and isinstance(item, str) and item and not Path(item).is_absolute()
                else absolute_asset_paths(item)
                for key, item in value.items()
            }
        return value

    return absolute_asset_paths(catalogue)


def unsupported_policy_copy(source: Path, destination: Path):
    """Remove only membrane-local permission, keeping real identified assets."""
    catalogue = catalogue_with_absolute_asset_paths(source)
    catalogue["membranePolicies"] = []
    catalogue["version"] += "-test-without-membrane-policy"
    destination.write_text(json.dumps(catalogue, indent=2) + "\n")


def invalid_asset_catalogue_copy(source: Path, destination: Path, fault: str):
    """Invalidate one selected species' asset while preserving other asset paths."""
    catalogue = catalogue_with_absolute_asset_paths(source)
    popc = next(item for item in catalogue["lipids"] if item["speciesId"] == "POPC")
    if fault == "missing-coordinate":
        missing = destination.parent / "intentionally-absent-POPC.cif"
        assert not missing.exists()
        popc["coordinateTemplatePath"] = str(missing)
    else:
        digest_key = {
            "coordinate-digest": "coordinateTemplateSha256",
            "force-field-digest": "templateSha256",
        }[fault]
        popc[digest_key] = "0" * 64
    catalogue["version"] += f"-test-{fault}"
    destination.write_text(json.dumps(catalogue, indent=2) + "\n")


class BrowserMembraneChoiceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (HOST, WORKER, LAUNCHER, ROOT / "out" / "host" / "wwwroot" / "index.html"):
            if not path.is_file():
                raise AssertionError(f"Build prerequisite is missing: {path}; run scripts/build-local.sh")
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_mixed_asymmetric_choice_requires_adoption_and_shows_assessed_evidence(self):
        catalogue = json.loads(POLICY.read_text())
        membrane_policy = catalogue["membranePolicies"][0]
        representations = {item["speciesId"]: item for item in catalogue["lipids"]}
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    before = current_state(page)
                    self.assertIsNone(before["membrane"])
                    self.assertIsNone(before["placement"])
                    self.assertEqual(before["stages"], [])
                    candidates = {item["speciesId"] for item in before["availableLipids"]}
                    self.assertTrue({"POPC", "POPE", "CHL1"} <= candidates,
                                    f"Qualified POPC/POPE/CHL1 catalogue unavailable: {sorted(candidates)}")

                    choose_fraction(page, "Upper", 1, "POPC", "60")
                    page.get_by_role("button", name="Add lipid").nth(0).click()
                    choose_fraction(page, "Upper", 2, "POPE", "40")
                    choose_fraction(page, "Lower", 1, "POPC", "25")
                    page.get_by_role("button", name="Add lipid").nth(1).click()
                    choose_fraction(page, "Lower", 2, "POPE", "75")
                    page.locator("#scientific-purpose").fill(PURPOSE)
                    page.get_by_role("button", name="Propose membrane model").click()
                    proposed = await_membrane_status(page, "proposed")
                    model_id = proposed["membrane"]["modelId"]
                    self.assertEqual(proposed["study"]["number"], before["study"]["number"])
                    self.assertEqual(proposed["membrane"]["upper"], [
                        {"speciesId": "POPC", "fraction": 0.6},
                        {"speciesId": "POPE", "fraction": 0.4},
                    ])
                    self.assertEqual(proposed["membrane"]["lower"], [
                        {"speciesId": "POPC", "fraction": 0.25},
                        {"speciesId": "POPE", "fraction": 0.75},
                    ])
                    self.assertIsNone(proposed["placement"])
                    inspect_membrane(page)
                    panel = page.locator(".evidence-panel")
                    expect(panel).to_contain_text(re.compile(r"proposed", re.I))
                    assert_exact_leaflets(self, panel.inner_text(),
                                          {"POPC": 60, "POPE": 40},
                                          {"POPC": 25, "POPE": 75})
                    expect(panel).to_contain_text(PURPOSE)
                    expect(panel).to_contain_text("0.15 M")
                    expect(panel).to_contain_text("charge-neutralizing compatible monovalent counterions")
                    expect(panel).to_contain_text("303 K")
                    expect(page.locator(".stage-strip")).to_contain_text("None")
                    capture_review(self, page, "proposed", 1672, 941)
                    capture_review(self, page, "proposed", 820, 760)

                    page.get_by_role("button", name="Show workflow").click()
                    page.get_by_role("button", name="Adopt and assess model").click()
                    assessed = await_membrane_status(page, "assessed")
                    self.assertEqual(assessed["membrane"]["modelId"], model_id)
                    self.assertEqual(assessed["study"]["number"], before["study"]["number"] + 1)
                    self.assertEqual(assessed["membrane"]["policyId"], membrane_policy["id"],
                                     assessed["membrane"])
                    self.assertEqual(assessed["membrane"]["policyVersion"], membrane_policy["version"])
                    self.assertIn(membrane_policy["limitations"][0], assessed["membrane"]["limitations"])
                    support = {item["speciesId"]: item for item in assessed["membrane"]["speciesSupport"]}
                    self.assertEqual(set(support), {"POPC", "POPE"})
                    for species_id in support:
                        identified = representations[species_id]
                        self.assertEqual(support[species_id]["chemistryId"], identified["chemistryId"])
                        self.assertEqual(support[species_id]["coordinateSha256"], identified["coordinateTemplateSha256"])
                        self.assertEqual(support[species_id]["parameterSha256"], identified["templateSha256"])
                        self.assertEqual(support[species_id]["forceFieldFamily"], identified["forceFieldFamily"])
                        self.assertEqual(support[species_id]["forceFieldVersion"], identified["forceFieldVersion"])
                    self.assertTrue(assessed["membrane"]["evidence"])
                    self.assertIsNone(assessed["placement"])
                    self.assertEqual(assessed["stages"], [])
                    inspect_membrane(page)
                    expect(panel).to_contain_text(re.compile(r"assessed", re.I))
                    assert_exact_leaflets(self, panel.inner_text(),
                                          {"POPC": 60, "POPE": 40},
                                          {"POPC": 25, "POPE": 75})
                    expect(panel).to_contain_text(PURPOSE)
                    expect(panel).to_contain_text("Lipid21")
                    expect(panel).to_contain_text(membrane_policy["id"])
                    expect(panel).to_contain_text(membrane_policy["version"])
                    expect(panel).to_contain_text(membrane_policy["limitations"][0])
                    expect(panel).to_contain_text(representations["POPC"]["chemistryId"])
                    expect(panel).to_contain_text(representations["POPE"]["chemistryId"])
                    for summary in panel.get_by_text("Coordinate and parameter asset identities").all():
                        summary.click()
                    for species_id in ("POPC", "POPE"):
                        expect(panel).to_contain_text(representations[species_id]["coordinateTemplateSha256"])
                        expect(panel).to_contain_text(representations[species_id]["templateSha256"])
                    expect(panel).to_contain_text("Protein placement and an assembled membrane remain separate")
                    expect(panel).to_contain_text("These disclosed targets are not achieved conditions")
                    self.assertNotRegex(panel.inner_text(),
                                        r"(?i)\b(?:cytoplasmic|extracellular|periplasmic|luminal)\s+leaflet\b")
                    self.assertNotRegex(panel.inner_text(), r"(?i)\bachieved\s+(?:molecule count|composition)\b")
                    for summary in panel.get_by_text("Coordinate and parameter asset identities").all():
                        summary.click()
                    expect(page.locator(".stage-strip")).to_contain_text("Stage assessment")
                    capture_review(self, page, "assessed", 1672, 941)
                    capture_review(self, page, "assessed", 820, 760)

                    page.get_by_role("button", name="Show workflow").click()
                    page.get_by_label("Lower leaflet percentage 1", exact=True).fill("30")
                    page.get_by_label("Lower leaflet percentage 2", exact=True).fill("70")
                    page.get_by_role("button", name="Propose membrane model").click()
                    revised_proposal = await_membrane_status(page, "proposed")
                    revised_id = revised_proposal["membrane"]["modelId"]
                    self.assertNotEqual(revised_id, model_id)
                    self.assertEqual(revised_proposal["study"]["number"], assessed["study"]["number"])
                    self.assertEqual(revised_proposal["membrane"]["lower"], [
                        {"speciesId": "POPC", "fraction": 0.3},
                        {"speciesId": "POPE", "fraction": 0.7},
                    ])
                    self.assertIsNone(revised_proposal["membrane"]["policyId"])
                    self.assertEqual(revised_proposal["membrane"]["speciesSupport"], [])
                    self.assertEqual(revised_proposal["membrane"]["evidence"], [])
                    self.assertIsNone(revised_proposal["placement"])
                    page.get_by_role("button", name="Adopt and assess model").click()
                    revised_assessed = await_membrane_status(page, "assessed")
                    self.assertEqual(revised_assessed["membrane"]["modelId"], revised_id)
                    self.assertEqual(revised_assessed["study"]["number"], assessed["study"]["number"] + 1)
                    self.assertEqual(revised_assessed["membrane"]["lower"], revised_proposal["membrane"]["lower"])
                    self.assertEqual(revised_assessed["membrane"]["policyId"], membrane_policy["id"])
                    self.assertIsNone(revised_assessed["placement"])
                    self.assertEqual(revised_assessed["stages"], [])
                finally:
                    browser.close()

    def test_incoherent_and_unsupported_choices_do_not_establish_membrane(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            if not POLICY.is_file():
                self.fail(f"Qualified slice 2 catalogue is unavailable: {POLICY}")
            unsupported = directory / "without-membrane-policy.json"
            unsupported_policy_copy(POLICY, unsupported)
            with running_host(directory / "workspace", unsupported) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    before = current_state(page)
                    self.assertIn("POPC", {item["speciesId"] for item in before["availableLipids"]})
                    choose_fraction(page, "Upper", 1, "POPC", "60")
                    choose_fraction(page, "Lower", 1, "POPC", "100")
                    page.locator("#scientific-purpose").fill("A deliberately simple physical bilayer.")
                    expect(page.get_by_role("button", name="Propose membrane model")).to_be_disabled()
                    self.assertIsNone(current_state(page)["membrane"])
                    expect(page.locator("#researcher-workflow")).to_contain_text("totaling 100%")

                    page.get_by_label("Upper leaflet percentage 1", exact=True).fill("100")
                    page.get_by_role("button", name="Propose membrane model").click()
                    proposed = await_membrane_status(page, "proposed")
                    model_id = proposed["membrane"]["modelId"]
                    self.assertEqual(proposed["study"]["number"], before["study"]["number"])
                    page.get_by_role("button", name="Adopt and assess model").click()
                    rejected = await_membrane_status(page, "notEstablished")
                    self.assertEqual(rejected["membrane"]["modelId"], model_id)
                    self.assertEqual(rejected["study"]["number"], before["study"]["number"] + 1)
                    self.assertIsNone(rejected["placement"])
                    self.assertEqual(rejected["stages"], [])
                    self.assertFalse(next(item["enabled"] for item in rejected["actions"]
                                          if item["kind"] == "proposePlacement"))
                    inspect_membrane(page)
                    panel = page.locator(".evidence-panel")
                    expect(panel).to_contain_text("not established", ignore_case=True)
                    assert_exact_leaflets(self, panel.inner_text(), {"POPC": 100}, {"POPC": 100})
                    self.assertRegex(panel.inner_text(), r"(?i)no.*policy|policy.*unavailable|unsupported")
                    expect(page.locator(".stage-strip")).to_contain_text("None")
                    capture_review(self, page, "not-established", 1672, 941)
                    capture_review(self, page, "not-established", 820, 760)

                    page.get_by_role("button", name="Show workflow").click()
                    expect(page.get_by_role("button", name="Propose membrane model")).to_be_enabled()
                    page.locator("#scientific-purpose").fill(
                        "Reconsider the unsupported membrane premise before placement.")
                    page.get_by_role("button", name="Propose membrane model").click()
                    revised = await_membrane_status(page, "proposed")
                    self.assertNotEqual(revised["membrane"]["modelId"], model_id)
                    self.assertEqual(revised["study"]["number"], rejected["study"]["number"])
                    self.assertIsNone(revised["placement"])
                finally:
                    browser.close()

    def test_exact_unlisted_species_is_retained_and_refused_without_substitution(self):
        unlisted = "UNREPRESENTED-LIPID-0001"
        protein_only_catalogue = ROOT / "config" / "policies" / "protein-slice1.json"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with running_host(directory / "workspace", protein_only_catalogue) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 941})
                    page.goto(base, wait_until="domcontentloaded")
                    self.assertEqual(current_state(page)["availableLipids"], [])
                    for side in ("Upper", "Lower"):
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("__other__")
                        page.get_by_label(f"{side} leaflet exact species ID 1", exact=True).fill(unlisted)
                        page.get_by_label(f"{side} leaflet percentage 1", exact=True).fill("100")
                    page.get_by_role("button", name="Add lipid").nth(0).click()
                    page.locator("#scientific-purpose").fill(
                        "Inspect the exact requested chemical identity before choosing a representation.")
                    expect(page.get_by_role("button", name="Propose membrane model")).to_be_enabled()
                    page.get_by_role("button", name="Propose membrane model").click()
                    proposed = await_membrane_status(page, "proposed")
                    self.assertEqual(proposed["membrane"]["upper"], [
                        {"speciesId": unlisted, "fraction": 1.0},
                    ])
                    self.assertEqual(proposed["membrane"]["lower"], [
                        {"speciesId": unlisted, "fraction": 1.0},
                    ])
                    inspect_membrane(page)
                    expect(page.locator(".evidence-panel")).to_contain_text(unlisted)
                    page.get_by_role("button", name="Adopt and assess model").click()
                    rejected = await_membrane_status(page, "notEstablished")
                    self.assertEqual(rejected["membrane"]["upper"][0]["speciesId"], unlisted)
                    self.assertEqual(rejected["membrane"]["lower"][0]["speciesId"], unlisted)
                    self.assertIsNone(rejected["placement"])
                    self.assertFalse(next(item["enabled"] for item in rejected["actions"]
                                          if item["kind"] == "proposePlacement"))
                    page.get_by_role("button", name="Show workflow").click()
                    inspect_membrane(page)
                    expect(page.locator(".evidence-panel")).to_contain_text(unlisted)
                    expect(page.locator(".evidence-panel")).to_contain_text(
                        re.compile(r"not established|no evidence-backed policy|unsupported", re.I))
                finally:
                    browser.close()

    def test_unavailable_or_mismatched_assets_refuse_exact_species_without_substitution(self):
        for fault in ("coordinate-digest", "force-field-digest", "missing-coordinate"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                invalid = directory / f"{fault}.json"
                invalid_asset_catalogue_copy(POLICY, invalid, fault)
                with running_host(directory / "workspace", invalid) as base, sync_playwright() as playwright:
                    browser = chromium(playwright)
                    try:
                        page = browser.new_page(viewport={"width": 1672, "height": 941})
                        page.goto(base, wait_until="domcontentloaded")
                        before = current_state(page)
                        candidates = {item["speciesId"] for item in before["availableLipids"]}
                        self.assertNotIn("POPC", candidates,
                                         f"POPC with an invalid {fault} asset was offered as qualified")
                        self.assertIn("POPE", candidates,
                                      "An unaffected candidate should remain available")
                        for side in ("Upper", "Lower"):
                            page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("__other__")
                            page.get_by_label(f"{side} leaflet exact species ID 1", exact=True).fill("POPC")
                            page.get_by_label(f"{side} leaflet percentage 1", exact=True).fill("100")
                        page.locator("#scientific-purpose").fill(
                            f"Assess exactly POPC despite an identified {fault} asset failure.")
                        page.get_by_role("button", name="Propose membrane model").click()
                        proposed = await_membrane_status(page, "proposed")
                        self.assertEqual(proposed["membrane"]["upper"],
                                         [{"speciesId": "POPC", "fraction": 1.0}])
                        self.assertEqual(proposed["membrane"]["lower"],
                                         [{"speciesId": "POPC", "fraction": 1.0}])
                        page.get_by_role("button", name="Adopt and assess model").click()
                        refused = await_membrane_status(page, "notEstablished")
                        self.assertEqual(refused["membrane"]["modelId"], proposed["membrane"]["modelId"])
                        self.assertEqual(refused["membrane"]["upper"], proposed["membrane"]["upper"])
                        self.assertEqual(refused["membrane"]["lower"], proposed["membrane"]["lower"])
                        self.assertIn("POPC", refused["membrane"]["reason"])
                        self.assertRegex(refused["membrane"]["reason"], r"(?i)representation")
                        self.assertEqual(refused["membrane"]["speciesSupport"], [])
                        self.assertIsNone(refused["placement"])
                        self.assertEqual(refused["stages"], [])
                        self.assertFalse(next(item["enabled"] for item in refused["actions"]
                                              if item["kind"] == "proposePlacement"))
                        inspect_membrane(page)
                        panel = page.locator(".evidence-panel")
                        expect(panel).to_contain_text("POPC")
                        expect(panel).to_contain_text("not established", ignore_case=True)
                        expect(panel).to_contain_text(refused["membrane"]["reason"])
                    finally:
                        browser.close()


if __name__ == "__main__":
    unittest.main()
