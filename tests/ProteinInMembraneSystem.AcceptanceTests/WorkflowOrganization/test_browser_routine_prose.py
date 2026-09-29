"""Small controlled accounts verify routine choices without mandatory prose.

The account and its placement are presentation fixtures, not scientific results.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import threading
import unittest

from playwright.sync_api import expect, sync_playwright


ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "routine-prose"
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/InspectAndReattach"))
from test_browser_presentation import CHROMIUM, account, controlled_account_server, inspection  # noqa: E402


class RoutineProseBrowserTests(unittest.TestCase):
    def test_membrane_edit_during_pending_use_keeps_submitted_snapshot(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            fixture = account()
            fixture["attempt"] = None
            fixture["inspection"] = None
            fixture["availableLipids"] = [
                {"speciesId": species, "displayName": species,
                 "chemistryId": species, "limitations": []}
                for species in ("DMPC", "POPC")
            ]
            fixture["actions"].append({"kind": "adoptMembrane", "subjectId": None,
                                       "enabled": True, "reason": None})
            host.replace(fixture)
            entered, release = threading.Event(), threading.Event()
            submitted = []

            def hold_use(handler, command):
                if command["kind"] != "adoptMembrane":
                    return False
                submitted.append(command)
                entered.set()
                release.wait(timeout=10)
                handler.json_response({"reason": "Controlled pending submission"}, 422)
                return True

            host.command_handler = hold_use
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                args=["--no-sandbox", "--disable-dev-shm-usage"])
            try:
                page = browser.new_page(viewport={"width": 1280, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Membrane", exact=True).click()
                for side in ("Upper", "Lower"):
                    page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                expect(page.locator(".membrane-draft-summary")).to_contain_text("DMPC 100.0%")
                self.assertIsNone(host.account["membrane"])
                page.get_by_role("button", name="Use this membrane").click()
                self.assertTrue(entered.wait(5), "The exact use command did not reach the controlled Host")
                page.get_by_label("Upper leaflet lipid 1", exact=True).select_option("POPC")
                expect(page.locator(".membrane-draft-summary")).to_contain_text("POPC 100.0%")
                self.assertEqual(submitted[0]["data"], {
                    "upper": [{"speciesId": "DMPC", "fraction": 1}],
                    "lower": [{"speciesId": "DMPC", "fraction": 1}],
                })
                self.assertEqual(len(submitted), 1)
                self.assertIsNone(host.account["membrane"])
                release.set()
                expect(page.locator(".action-feedback.error")).to_be_visible()
            finally:
                release.set()
                browser.close()

    def test_site_specific_approval_uses_inspected_choice_without_a_prose_field(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            fixture = account()
            fixture["attempt"] = None
            proposal_id = "fixture-histidine-choice"
            fixture["protein"] = {
                "subjectId": "fixture-protein", "status": "review",
                "summary": "One histidine state awaits review", "atomCount": None,
                "candidateId": None, "prediction": None, "geometry": None,
                "sourceGeometry": None, "findings": [],
                "changes": [{"id": proposal_id, "studyRevisionId": "revision-one",
                             "intendedProteinId": "fixture-intended", "kind": "residueState",
                             "proposedChange": "HID",
                             "rationale": "Inspect this exact site and its evidence.",
                             "limitations": [], "approvalRequired": True,
                             "residue": {"model": 0, "chain": "A", "residue": 7,
                                         "insertionCode": "", "copyId": "A"}}]
            }
            fixture["inspection"] = inspection(proposal_id, "preparationChange")
            fixture["preparationReview"] = {
                "studyRevisionId": "revision-one", "intendedProteinId": "fixture-intended",
                "confirmedCount": 0, "remainingCount": 1, "blockers": [],
                "preparationStanding": "awaitingDecisions", "preparationMessage": None,
                "decisions": [{"id": "fixture-histidine-site", "kind": "residueState",
                               "residue": {"model": 0, "chain": "A", "residue": 7,
                                           "insertionCode": "", "copyId": "A"},
                               "partnerResidue": None, "atomName": None,
                               "standing": "pending", "chosenProposalId": None, "blocker": None,
                               "options": [{"proposalId": proposal_id, "proposedChange": "HID",
                                            "disposition": "available", "decisionId": None,
                                            "evidenceCurrent": True, "confirmationBlocker": None,
                                            "startsPreparationOnConfirmation": True,
                                            "startsPreparationOnDecline": False,
                                            "evidence": []}]}]}
            fixture["actions"].extend([
                {"kind": "approvePreparationChange", "subjectId": proposal_id,
                 "enabled": True, "reason": None},
                {"kind": "declinePreparationChange", "subjectId": proposal_id,
                 "enabled": True, "reason": None},
            ])
            host.replace(fixture)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                args=["--no-sandbox", "--disable-dev-shm-usage"])
            page = browser.new_page(viewport={"width": 1672, "height": 941})
            sent = []
            page.on("request", lambda request: sent.append(json.loads(request.post_data))
                    if request.method == "POST" and request.url.endswith("/api/commands") else None)
            try:
                page.goto(base, wait_until="domcontentloaded")
                expect(page.get_by_text("Scientific rationale for this")).to_have_count(0)
                expect(page.get_by_placeholder("Explain why this is the intended model assumption"))\
                    .to_have_count(0)
                approve = page.get_by_role("button", name="Confirm state and prepare protein")
                expect(approve).to_be_enabled()
                page.screenshot(path=str(ARTIFACTS / "fixture-site-specific-review.png"),
                                animations="disabled")
                with page.expect_request("**/api/commands"):
                    approve.click()
                decision = next(item for item in sent if item["kind"] == "approvePreparationChange")
                self.assertEqual(decision["data"],
                                 {"proposalId": proposal_id, "approve": True})
            finally:
                browser.close()

    def test_structured_choices_and_correction_send_no_routine_prose(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            fixture = account()
            fixture["attempt"] = None
            fixture["study"]["selectedSourceId"] = "fixture:source"
            fixture["inspection"] = inspection("fixture:source", "structuralSource")
            fixture["sourceModels"] = [{
                "index": 0, "atomCount": 8,
                "chains": [{"name": "A", "residueCount": 1, "atomCount": 4}],
                "assemblies": [],
                "partners": [{"sourceId": "cofactor-1", "label": "Bound cofactor",
                              "kind": "heterogen", "atomCount": 4, "chain": "A", "residue": 1}],
                "residues": [{"address": {"model": 0, "chain": "A", "residue": 1,
                                           "insertionCode": "", "copyId": ""},
                              "name": "ALA", "residueKind": "protein",
                              "backboneHeavyAtomsComplete": True, "alternateLocations": []}]
            }]
            fixture["availableLipids"] = [{"speciesId": "DMPC", "displayName": "DMPC",
                                           "chemistryId": "DMPC", "limitations": []}]
            fixture["placement"] = {
                "proposalId": "fixture-placement", "status": "supported",
                "preparedProteinId": "fixture-protein", "membraneModelId": "fixture-membrane",
                "topologyKind": "membrane-spanning", "physicalSide": "both",
                "midplaneAngstrom": 0, "thicknessAngstrom": 30,
                "depthAngstrom": 0, "tiltDegrees": 0, "sidedness": None,
                "contactingRegions": [], "limitations": [], "policyId": "fixture-policy",
                "policyVersion": "1", "witnessId": "fixture-witness",
                "reason": "Controlled presentation fixture; no scientific support is asserted here.",
                "evidence": [], "prediction": None
            }
            fixture["actions"].extend({"kind": kind, "subjectId": None,
                                       "enabled": True, "reason": None} for kind in
                                      ("selectProteinModel", "adoptMembrane", "revisePlacement"))
            host.replace(fixture)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-webgl",
                      "--use-gl=angle", "--use-angle=swiftshader"])
            page = browser.new_page(viewport={"width": 1672, "height": 941})
            sent = []
            page.on("request", lambda request: sent.append(json.loads(request.post_data))
                    if request.method == "POST" and request.url.endswith("/api/commands") else None)
            try:
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".sole-model")).to_contain_text("Only coordinate model")
                page.get_by_label("Chain A", exact=True).check()
                page.locator(".partner-members summary").click()
                page.get_by_label("Exclude").check()
                expect(page.get_by_role("button", name="Assess selected protein")).to_be_enabled()
                expect(page.get_by_placeholder("Reason for this decision")).to_have_count(0)
                page.get_by_role("button", name="Assess selected protein").scroll_into_view_if_needed()
                page.screenshot(path=str(ARTIFACTS / "fixture-protein-partner.png"), animations="disabled")
                page.get_by_role("button", name="Assess selected protein").click()
                expect(page.locator(".action-feedback.error")).to_be_visible()
                partner = next(item for item in sent if item["kind"] == "selectProteinModel")
                self.assertEqual(partner["data"]["partners"],
                                 [{"sourceId": "cofactor-1", "retain": False}])

                page.get_by_role("button", name="Membrane", exact=True).click()
                for side in ("Upper", "Lower"):
                    page.get_by_label(f"{side} leaflet lipid 1", exact=True).select_option("DMPC")
                expect(page.locator("#scientific-purpose")).to_have_count(0)
                expect(page.get_by_role("button", name="Use this membrane")).to_be_enabled()
                page.screenshot(path=str(ARTIFACTS / "fixture-membrane.png"), animations="disabled")
                page.get_by_role("button", name="Use this membrane").click()
                expect(page.locator(".action-feedback.error")).to_be_visible()
                membrane = next(item for item in sent if item["kind"] == "adoptMembrane")
                self.assertNotIn("scientificPurpose", membrane["data"])
                self.assertEqual(membrane["data"]["upper"],
                                 [{"speciesId": "DMPC", "fraction": 1}])

                page.get_by_role("button", name="Placement", exact=True).click()
                page.get_by_label("Depth shift (Å)").fill("1")
                expect(page.locator("#placement-rationale")).to_have_count(0)
                expect(page.get_by_role("button", name="Reassess corrected placement")).to_be_enabled()
                page.get_by_role("button", name="Reassess corrected placement").scroll_into_view_if_needed()
                page.screenshot(path=str(ARTIFACTS / "fixture-placement-correction.png"), animations="disabled")
                page.get_by_role("button", name="Reassess corrected placement").click()
                expect(page.locator(".action-feedback.error")).to_be_visible()
                corrected = next(item for item in sent if item["kind"] == "revisePlacement")
                self.assertNotIn("rationale", corrected["data"])
                self.assertEqual(corrected["data"]["depthShiftAngstrom"], 1)
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
