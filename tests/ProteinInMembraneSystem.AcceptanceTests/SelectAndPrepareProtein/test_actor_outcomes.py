"""Connected task outcomes with real small protein work and a controlled PPM failure.

The two-alanine structure is a software fixture. The executable below deliberately
fails after a pause; its placement screens demonstrate UI behavior, not a scientific
placement calculation or support for that construct.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import time
import unittest

from playwright.sync_api import expect, sync_playwright

from test_browser_route import (ROOT, chromium, horizontally_overflowing_regions,
                                running_host, two_alanines)
from test_browser_decision_review import two_site_fixture


POLICY = ROOT / "config" / "policies" / "protein-membrane-slice3.json"
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "actor-outcomes"


def fake_provider_catalogue(destination: Path, executable: Path) -> Path:
    """Keep exact assets but identify this test-local failing executable honestly."""
    catalogue = json.loads(POLICY.read_text())
    def absolute_assets(value):
        if isinstance(value, list):
            return [absolute_assets(item) for item in value]
        if isinstance(value, dict):
            return {key: str((POLICY.parent / item).resolve())
                    if key in {"path", "coordinateTemplatePath", "templatePath", "ppmResidueLibraryPath"}
                    and isinstance(item, str) and item and not Path(item).is_absolute()
                    else absolute_assets(item) for key, item in value.items()}
        return value
    catalogue = absolute_assets(catalogue)
    catalogue["ppmExecutableSha256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
    catalogue["version"] += "-controlled-provider-failure"
    destination.write_text(json.dumps(catalogue, indent=2) + "\n")
    return destination


def account(page):
    return page.evaluate("async () => (await (await fetch('/api/state', {cache: 'no-store'})).json())")


def await_account(page, predicate, label: str, seconds: float = 120):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        state = account(page)
        if predicate(state):
            return state
        time.sleep(0.1)
    raise AssertionError(f"Timed out waiting for {label}: {state.get('proteinTask')}, {state.get('membrane')}, {state.get('placementTask')}")


class ActorOutcomeBrowserTests(unittest.TestCase):
    def test_final_confirmation_reports_preparation_without_stealing_the_view(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "one-chain-histidine-software-fixture.pdb"
            two_site_fixture(source)
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1024, "height": 768})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    await_account(page, lambda state: state["preparationReview"] is not None and
                                  len(state["preparationReview"]["decisions"]) == 2,
                                  "site and repair review")
                    page.locator(".review-option").filter(has_text="HIE").locator("input").check()
                    expect(page.get_by_role("button", name="Confirm state")).to_be_enabled()
                    page.get_by_role("button", name="Confirm state").click()
                    await_account(page, lambda state: state["preparationReview"]["confirmedCount"] == 1,
                                  "recorded histidine state")
                    next_repair = page.get_by_role("button", name="Next unresolved repair")
                    if next_repair.is_visible():
                        next_repair.click()
                    else:
                        details = page.locator(".review-other")
                        if details.get_attribute("open") is None:
                            details.locator("summary").click()
                        details.locator(".review-site").click()
                    expect(page.locator(".review-site-scope")).to_contain_text("Missing-atom repair")
                    page.locator(".review-option input").check()
                    final = page.get_by_role("button", name="Approve repair and prepare protein")
                    expect(final).to_be_enabled()
                    expect(page.locator(".review-consequence").first).to_contain_text(
                        "resulting protein checked")
                    before = account(page)
                    view_subject = before["inspection"]["subjectId"]
                    page.screenshot(path=str(ARTIFACTS / "final-confirmation-ready-1024.png"),
                                    full_page=True, animations="disabled")
                    final.click()
                    terminal = await_account(page, lambda state: state["proteinTask"] is not None and
                                             state["proteinTask"]["standing"] in ("assessed", "failed", "unavailable"),
                                             "preparation outcome")
                    self.assertEqual(terminal["preparationReview"]["confirmedCount"], 2)
                    self.assertEqual(terminal["inspection"]["subjectId"], view_subject)
                    if terminal["proteinTask"]["standing"] == "assessed":
                        expect(page.get_by_role("region", name="Protein task outcome"))\
                            .to_contain_text("Protein prepared")
                        expect(page.get_by_role("button", name="View prepared protein")).to_be_visible()
                    else:
                        expect(page.get_by_role("region", name="Protein task outcome"))\
                            .to_contain_text(terminal["proteinTask"]["message"])
                        expect(page.get_by_role("button", name="Review confirmed choices")).to_be_visible()
                    page.screenshot(path=str(ARTIFACTS / "final-confirmation-outcome-1024.png"),
                                    full_page=True, animations="disabled")
                finally:
                    browser.close()

    def test_no_choice_completion_elsewhere_membrane_and_failed_position(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-software-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            executable = directory / "controlled-failing-ppm"
            executable.write_text("#!/bin/sh\nsleep 2\nexit 2\n")
            executable.chmod(0o755)
            catalogue = fake_provider_catalogue(directory / "controlled-policy.json", executable)
            with running_host(directory / "workspace", catalogue, executable) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    commands = []
                    page.on("request", lambda request: commands.append(request.post_data_json)
                            if request.url.endswith("/api/commands") and request.post_data else None)
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Membrane").click()
                    expect(page.get_by_role("button", name="Membrane", exact=True)).to_have_attribute("aria-current", "page")
                    await_account(page, lambda state: state["proteinTask"] is not None and
                                  state["proteinTask"]["standing"] == "assessed", "assessed protein")
                    expect(page.locator("#membrane-workflow")).to_be_visible()
                    self.assertEqual(account(page)["proteinTask"]["standing"], "assessed")
                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Protein").click()
                    expect(page.get_by_role("region", name="Protein task outcome"))\
                        .to_contain_text("Protein prepared")
                    expect(page.get_by_role("heading", name=f"Protein before preparation · {source.name}"))\
                        .to_be_visible()
                    page.screenshot(path=str(ARTIFACTS / "no-choice-completed-elsewhere-1672.png"),
                                    full_page=True, animations="disabled")

                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Membrane").click()
                    for side in ("Upper", "Lower"):
                        page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                        page.get_by_label(f"{side} leaflet percentage 1", exact=True).fill("100")
                    page.get_by_role("button", name="Propose membrane model").click()
                    await_account(page, lambda state: state["membrane"] is not None and
                                  state["membrane"]["status"] == "proposed", "membrane proposal")
                    expect(page.locator("#membrane-workflow")).to_contain_text("Proposal awaiting adoption")
                    page.screenshot(path=str(ARTIFACTS / "membrane-proposed-1672.png"),
                                    full_page=True, animations="disabled")
                    page.get_by_role("button", name="Adopt and assess displayed proposal").click()
                    expect(page.locator("#membrane-workflow")).to_contain_text(
                        "assessing membrane support", timeout=15000)
                    await_account(page, lambda state: state["membrane"] is not None and
                                  state["membrane"]["status"] == "assessed", "assessed membrane")
                    expect(page.locator("#membrane-workflow")).to_contain_text(
                        "membrane support established")
                    page.get_by_role("button", name="View proposal account").click()
                    expect(page.get_by_role("heading", name="Intended membrane model")).to_be_visible()
                    expect(page.locator(".scene-subtitle")).not_to_contain_text(
                        account(page)["membrane"]["modelId"])
                    page.screenshot(path=str(ARTIFACTS / "membrane-assessed-1672.png"),
                                    full_page=True, animations="disabled")

                    page.get_by_role("navigation", name="Research work areas")\
                        .get_by_role("button", name="Placement").click()
                    page.locator("#topology-kind").select_option("membrane-spanning")
                    page.locator("#ppm-nterminal-side").select_option("in")
                    expect(page.get_by_role("button", name="Obtain and assess position")).to_be_enabled()
                    page.get_by_role("button", name="Obtain and assess position").click()
                    expect(page.get_by_role("region", name="Placement task outcome"))\
                        .to_contain_text("Obtaining a position", timeout=15000)
                    self.assertEqual(account(page)["placementTask"]["standing"], "obtaining")
                    page.screenshot(path=str(ARTIFACTS / "placement-obtaining-1672.png"),
                                    full_page=True, animations="disabled")
                    await_account(page, lambda state: state["placementTask"] is not None and
                                  state["placementTask"]["standing"] == "noProposal", "no-position outcome")
                    self.assertIsNone(account(page)["placement"])
                    expect(page.get_by_role("region", name="Placement task outcome"))\
                        .to_contain_text("No position established")
                    expect(page.get_by_label("Current placement outcome"))\
                        .to_contain_text("No reviewable position established")
                    expect(page.locator("#topology-kind")).to_have_value("membrane-spanning")
                    expect(page.locator("#ppm-nterminal-side")).to_have_value("in")
                    for width, height in ((1672, 940), (1024, 768), (820, 760)):
                        page.set_viewport_size({"width": width, "height": height})
                        self.assertEqual(horizontally_overflowing_regions(page), [])
                        page.screenshot(path=str(ARTIFACTS / f"placement-no-proposal-{width}.png"),
                                        full_page=True, animations="disabled")
                    expect(page.get_by_role("button", name="Obtain and assess position")).to_be_enabled()
                    page.get_by_role("button", name="Obtain and assess position").click()
                    await_account(page, lambda state: state["placementTask"] is not None and
                                  state["placementTask"]["standing"] == "noProposal", "retry outcome")
                    self.assertEqual(sum(command["kind"] == "proposePlacement" for command in commands), 2)
                    self.assertEqual(account(page)["stages"], [])
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
