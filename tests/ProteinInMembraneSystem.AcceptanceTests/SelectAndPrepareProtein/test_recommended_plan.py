"""Real host/worker/browser crossing for explicit starting-state plan authorization.

The two-histidine PDB is a small software fixture.  Its provider result checks
the implementation route, not a biological protonation-state claim.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

from playwright.sync_api import expect, sync_playwright

from test_browser_route import ROOT, chromium, running_host, two_alanines

sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" / "InspectAndReattach"))
from test_browser_presentation import account, controlled_account_server  # noqa: E402


POLICY = ROOT / "config/policies/protein-membrane-current.json"
CAPTURES = ROOT / "out/reconciled-design/recommendation-browser"


def atom_line(serial: int, name: str, residue: str, chain: str, number: int,
              xyz: tuple[float, float, float], element: str) -> str:
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:4s} {residue:3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n")


def two_histidine_fixture(path: Path) -> None:
    source = path.with_name("two-alanines-before-mutation.pdb")
    source.write_text(two_alanines(1, "A") + "END\n", encoding="utf-8")
    result = subprocess.run([str(ROOT / "out/python/bin/python"), "-c", """
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(filename=sys.argv[1])
fixer.applyMutations(['ALA-1-HIS', 'ALA-2-HIS'], 'A')
fixer.missingResidues = {}
fixer.findMissingAtoms()
fixer.addMissingAtoms()
with open(sys.argv[2], 'w', encoding='utf-8') as stream:
    PDBFile.writeFile(fixer.topology, fixer.positions, stream, keepIds=True)
""", str(source), str(path)], capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        raise AssertionError(result.stderr)


def state_of(page):
    return page.evaluate("async () => (await (await fetch('/api/state', {cache: 'no-store'})).json())")


def until(page, predicate, label: str, seconds: float = 120):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        state = state_of(page)
        if predicate(state):
            return state
        time.sleep(0.1)
    raise AssertionError(f"Timed out waiting for {label}: "
                         f"{state.get('preparationPlan')}, {state.get('proteinTask')}")


def controlled_plan(standing: str, *, authorization_available: bool = False,
                    has_overrides: bool = False) -> dict:
    value = account()
    value["attempt"] = None
    value["inspection"] = None
    value["preparationPlan"] = {
        "standing": standing, "planId": "controlled-plan", "planSha256": "a" * 64,
        "method": "OpenMM Modeller.addHydrogens starting-state rule · 8.1",
        "nominalPh": 7.0, "stateChoiceCount": 2, "heavyAtomCount": 1,
        "removedSourceHydrogenCount": 3, "choices": [],
        "message": "One exact conformer choice remains unresolved." if standing == "partial" else
                   "The suggestion operation could not finish." if standing == "failed" else
                   "Controlled plan standing",
        "authorizationAvailable": authorization_available, "hasOverrides": has_overrides,
    }
    value["actions"] = [{"kind": kind, "subjectId": subject, "enabled": True, "reason": None}
                        for kind, subject in (("authorizePreparationPlan", "controlled-plan"),
                                              ("retryPreparationPlan", None))]
    return value


class RecommendationBrowserTests(unittest.TestCase):
    def test_controlled_plan_states_keep_actions_and_method_details_distinct(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            browser = chromium(playwright)
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                scenarios = (
                    ("calculating", False, False, "Finding preparation suggestions…"),
                    ("ready", True, False, "Preparation suggestions ready"),
                    ("ready", True, True, "Updated choices ready"),
                    ("ready", False, False, "Suggestions need rechecking"),
                    ("partial", False, False, "Choices need your input"),
                    ("failed", False, False, "Recommendations unavailable"),
                )
                for standing, available, overridden, heading in scenarios:
                    with self.subTest(standing=standing, available=available, overridden=overridden):
                        host.replace(controlled_plan(standing, authorization_available=available,
                                                     has_overrides=overridden))
                        page.goto(base)
                        card = page.get_by_role("region", name="Preparation plan")
                        expect(card.get_by_role("heading", name=heading)).to_be_visible()
                        self.assertNotIn("a" * 64, card.inner_text())
                        if standing == "ready":
                            expect(card).to_contain_text("2 state choices")
                            expect(card).to_contain_text("1 atom repairs")
                            expect(card).to_contain_text("pH 7")
                            expect(card).to_contain_text("Changes not yet applied")
                            expect(card.get_by_role("button", name="Review or change choices"))\
                                .to_be_visible()
                            if available:
                                name = "Prepare with these choices" if overridden else "Prepare with recommendations"
                                expect(card.get_by_role("button", name=name)).to_be_enabled()
                            else:
                                expect(card).to_contain_text("no longer matches the current recorded choices")
                                expect(card.get_by_role("button", name="Prepare with recommendations"))\
                                    .to_have_count(0)
                            card.get_by_text("Method details").click()
                            expect(card).to_contain_text("OpenMM")
                            expect(card).to_contain_text("8.1")
                            self.assertNotIn("Modeller.addHydrogens", card.inner_text())
                            self.assertNotIn("a" * 64, card.inner_text())
                        elif standing == "calculating":
                            expect(card).to_contain_text("Finding suggested starting states")
                            expect(card.get_by_role("button", name="Prepare with recommendations"))\
                                .to_have_count(0)
                        elif standing == "partial":
                            expect(card).to_contain_text("One exact conformer choice remains unresolved")
                            expect(card.get_by_role("button", name="Prepare with recommendations"))\
                                .to_have_count(0)
                        else:
                            expect(card).to_contain_text("The suggestion operation could not finish")
                            expect(card.get_by_role("button", name="Retry suggestions")).to_be_enabled()
                            expect(card.get_by_role("button", name="Review choices manually")).to_be_visible()
            finally:
                browser.close()

    def test_failed_authorized_candidate_keeps_choices_and_retries_exact_plan(self):
        CAPTURES.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-histidines-failure-software-fixture.pdb"
            two_histidine_fixture(source)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
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
                    ready = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                  "ready", "checked plan before controlled artifact change")
                    digest = ready["preparationPlan"]["planSha256"]
                    correspondence = next((directory / "workspace").rglob("protein-correspondence.json"))
                    original = correspondence.read_bytes()
                    correspondence.write_bytes(original + b"\n")
                    page.get_by_role("button", name="Prepare with recommendations").click()
                    failed = until(page, lambda current: (current.get("proteinTask") or {}).get("standing") ==
                                   "failed", "controlled artifact-integrity failure")
                    self.assertIn("checked plan artifact changed", failed["proteinTask"]["message"])
                    self.assertEqual(failed["preparationPlan"]["standing"], "failedAfterAuthorization")
                    self.assertEqual(failed["preparationPlan"]["planSha256"], digest)
                    self.assertEqual(failed["preparationReview"]["confirmedCount"], 2)
                    self.assertIsNone(failed["proteinTask"]["preparedProteinId"])
                    expect(page.get_by_role("region", name="Current protein result"))\
                        .to_contain_text(failed["proteinTask"]["message"])
                    expect(page.get_by_role("button", name="Review confirmed choices")).to_be_visible()
                    page.screenshot(path=str(CAPTURES / "authorized-plan-failed-1024.png"),
                                    full_page=True, animations="disabled")
                    correspondence.write_bytes(original)
                    page.get_by_role("button", name="Retry protein preparation").click()
                    recovered = until(page, lambda current: (current.get("proteinTask") or {}).get("standing") ==
                                      "assessed", "same-plan retry after restoring checked artifact")
                    self.assertEqual(recovered["preparationPlan"]["planSha256"], digest)
                    self.assertEqual(recovered["preparationPlan"]["standing"], "applied")
                finally:
                    browser.close()

    def test_three_observed_bonds_need_exact_review_before_plan(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein" / \
                "fixtures/1CRN-RCSB.pdb"
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1024, "height": 768})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.locator("#assembly-choice").select_option("deposited")
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    partial = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                    "partial", "three bond exceptions")
                    self.assertEqual(len([item for item in partial["preparationReview"]["decisions"]
                                          if item["kind"] == "disulfide"]), 3)
                    plan_card = page.get_by_role("region", name="Preparation plan")
                    expect(plan_card.get_by_role("heading", name="Choices need your input"))\
                        .to_be_visible()
                    expect(plan_card.get_by_role("button", name="Prepare with recommendations"))\
                        .to_have_count(0)
                    expect(page.get_by_role("region", name="Protein preparation decisions"))\
                        .to_contain_text("3 structural exceptions")
                    CAPTURES.mkdir(parents=True, exist_ok=True)
                    for width in (1672, 1024, 820):
                        page.set_viewport_size({"width": width, "height": 768})
                        page.screenshot(path=str(CAPTURES / f"three-bond-exceptions-{width}.png"),
                                        full_page=True, animations="disabled")
                    page.set_viewport_size({"width": 1024, "height": 768})
                    for remaining in (2, 1, 0):
                        page.get_by_role("button", name="Confirm possible bond").click()
                        if remaining:
                            until(page, lambda current: sum(item["kind"] == "disulfide" and
                                item["standing"] == "pending" for item in
                                current["preparationReview"]["decisions"]) == remaining,
                                "next unresolved bond")
                            page.get_by_role("button", name="Next unresolved bond decision").click()
                    ready = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing")
                                  in ("ready", "failed"), "joint bond plan")
                    self.assertEqual(ready["preparationPlan"]["standing"], "ready",
                                     ready["preparationPlan"])
                    self.assertEqual(ready["preparationReview"]["confirmedCount"] >= 3, True)
                    self.assertIsNone(ready["proteinTask"]["preparedProteinId"])
                    page.screenshot(path=str(CAPTURES / "three-bond-plan-ready-1024.png"),
                                    full_page=True, animations="disabled")
                finally:
                    browser.close()

    def test_conformer_exception_recalculates_then_requires_plan_authorization(self):
        CAPTURES.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-alternate-software-fixture.pdb"
            original = atom_line(5, "CB", "ALA", "A", 1, (2, -0.77, -1.2), "C")
            shifted = atom_line(12, "CB", "ALA", "A", 1, (2.1, -0.77, -1.2), "C")
            source.write_text(two_alanines(1, "A").replace(
                original, original[:16] + "A" + original[17:] +
                shifted[:16] + "B" + shifted[17:]) + "END\n", encoding="utf-8")
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1024, "height": 768})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A").check()
                    page.locator(".altloc-row select").select_option(label="A")
                    page.get_by_role("button", name="Assess selected protein").click()
                    partial = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                    "partial", "conformer exception")
                    self.assertEqual(len([item for item in partial["preparationReview"]["decisions"]
                                          if item["kind"] == "alternateLocation"]), 1)
                    expect(page.get_by_role("region", name="Preparation plan")).to_contain_text(
                        "Choices need your input")
                    page.screenshot(path=str(CAPTURES / "conformer-exception-1024.png"),
                                    full_page=True, animations="disabled")
                    page.get_by_role("button", name="Confirm conformer").click()
                    ready = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                  "ready", "rechecked conformer plan")
                    self.assertIsNone(ready["proteinTask"]["preparedProteinId"])
                    self.assertEqual(ready["preparationReview"]["confirmedCount"], 1)
                    page.screenshot(path=str(CAPTURES / "conformer-plan-ready-1024.png"),
                                    full_page=True, animations="disabled")
                    page.get_by_role("button", name="Prepare with recommendations").click()
                    assessed = until(page, lambda current: (current.get("proteinTask") or {}).get("standing")
                                     in ("assessed", "failed"), "conformer plan outcome")
                    self.assertEqual(assessed["proteinTask"]["standing"], "assessed", assessed["proteinTask"])
                    self.assertEqual(assessed["preparationPlan"]["standing"], "applied")
                finally:
                    browser.close()

    def test_no_choice_plan_waits_for_explicit_authorization(self):
        CAPTURES.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-no-choice-fixture.pdb"
            source.write_text(two_alanines(1, "A") + "END\n", encoding="utf-8")
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
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
                    ready = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                  "ready", "no-choice plan")
                    self.assertEqual(ready["preparationPlan"]["stateChoiceCount"], 0)
                    self.assertEqual(ready["preparationPlan"]["heavyAtomCount"], 0)
                    self.assertIsNone(ready["proteinTask"]["preparedProteinId"])
                    expect(page.get_by_role("button", name="Prepare with recommendations")).to_be_enabled()
                    page.screenshot(path=str(CAPTURES / "no-choice-ready-1024.png"),
                                    full_page=True, animations="disabled")
                    view_subject = ready["inspection"]["subjectId"]
                    page.get_by_role("button", name="Prepare with recommendations").click()
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Membrane").click()
                    assessed = until(page, lambda current: (current.get("proteinTask") or {}).get("standing")
                                     in ("assessed", "failed"), "no-choice preparation outcome")
                    self.assertEqual(assessed["proteinTask"]["standing"], "assessed")
                    self.assertEqual(assessed["preparationPlan"]["standing"], "applied")
                    self.assertEqual(assessed["inspection"]["subjectId"], view_subject)
                    expect(page.get_by_role("button", name="Membrane", exact=True)).to_have_attribute(
                        "aria-current", "page")
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Protein").click()
                    expect(page.get_by_role("region", name="Protein task outcome"))\
                        .to_contain_text("Protein prepared")
                    page.screenshot(path=str(CAPTURES / "no-choice-applied-1024.png"),
                                    full_page=True, animations="disabled")
                finally:
                    browser.close()

    def test_recommend_override_authorize_and_preserve_work_area(self):
        CAPTURES.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-histidines-software-fixture.pdb"
            two_histidine_fixture(source)
            with running_host(directory / "workspace", POLICY) as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 940})
                    page.goto(base)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    ready = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing")
                                  in ("ready", "failed"), "recommendation outcome")
                    self.assertEqual(ready["preparationPlan"]["standing"], "ready", ready["preparationPlan"])
                    self.assertEqual(ready["preparationPlan"]["stateChoiceCount"], 2)
                    self.assertIsNone(ready["proteinTask"]["preparedProteinId"])
                    first_digest = ready["preparationPlan"]["planSha256"]
                    for width in (1672, 1024, 820):
                        page.set_viewport_size({"width": width, "height": 940 if width == 1672 else 768})
                        plan_card = page.get_by_role("region", name="Preparation plan")
                        expect(plan_card.get_by_role("heading", name="Preparation suggestions ready"))\
                            .to_be_visible()
                        expect(plan_card).to_contain_text("2 state choices")
                        expect(plan_card).to_contain_text("Changes not yet applied")
                        self.assertNotIn(first_digest, plan_card.inner_text())
                        page.screenshot(path=str(CAPTURES / f"ready-{width}.png"),
                                        full_page=True, animations="disabled")
                    page.set_viewport_size({"width": 1024, "height": 768})
                    plan_card.get_by_text("Method details").click()
                    self.assertNotIn("Modeller.addHydrogens", plan_card.inner_text())
                    self.assertNotIn(first_digest, plan_card.inner_text())
                    page.get_by_role("button", name="Review or change choices").click()
                    expect(page.locator(".review-site-list .review-site")).to_have_count(2)
                    expect(page.locator(".review-site-list")).to_contain_text("Suggested")
                    page.locator(".review-option").filter(has_text="HIP").locator("input").check()
                    expect(page.get_by_role("button", name="Use this choice")).to_be_enabled()
                    page.get_by_role("button", name="Use this choice").click()
                    updated = until(page, lambda current: (current.get("preparationPlan") or {}).get("standing") ==
                                    "ready" and current["preparationPlan"]["planSha256"] != first_digest,
                                    "rechecked override")
                    self.assertTrue(updated["preparationPlan"]["hasOverrides"])
                    self.assertEqual(updated["preparationPlan"]["choices"][0]["variant"], "HIP")
                    self.assertTrue(updated["preparationPlan"]["choices"][0]["overridden"])
                    page.screenshot(path=str(CAPTURES / "overridden-1024.png"),
                                    full_page=True, animations="disabled")
                    page.get_by_role("button", name="Hide choice review").click()
                    prepare = page.get_by_role("button", name="Prepare with these choices")
                    expect(prepare).to_be_enabled()
                    view_subject = updated["inspection"]["subjectId"]
                    prepare.click()
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Membrane").click()
                    completed = until(page, lambda current: (current.get("proteinTask") or {}).get("standing")
                                      in ("assessed", "failed", "unavailable"), "authorized protein outcome")
                    self.assertEqual(completed["proteinTask"]["standing"], "assessed", completed["proteinTask"])
                    self.assertEqual(completed["preparationPlan"]["standing"], "applied")
                    self.assertEqual(completed["preparationPlan"]["planSha256"],
                                     updated["preparationPlan"]["planSha256"])
                    self.assertEqual(completed["inspection"]["subjectId"], view_subject)
                    expect(page.get_by_role("button", name="Membrane", exact=True)).to_have_attribute(
                        "aria-current", "page")
                    page.get_by_role("navigation", name="Research work areas").get_by_role(
                        "button", name="Protein").click()
                    expect(page.get_by_role("region", name="Protein task outcome")).to_contain_text(
                        "Protein prepared")
                    page.screenshot(path=str(CAPTURES / "applied-1024.png"),
                                    full_page=True, animations="disabled")
                    correspondence = next((directory / "workspace").rglob("protein-correspondence.json"))
                    observed = json.loads(correspondence.read_text())
                    self.assertTrue(observed["complete"])
                    self.assertTrue(all(item["sourceResidue"]["copyId"] == "A"
                                        for item in observed["atoms"]))
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
