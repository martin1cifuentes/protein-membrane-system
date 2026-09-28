"""Exact protein-review route and labelled density presentation through Chromium.

The 38-site account below is a presentation fixture. Its synthetic addresses and
decisions do not establish a scientifically suitable protein or provider output.
"""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import time
import unittest

from playwright.sync_api import expect, sync_playwright

from test_browser_route import ROOT, chromium, running_host


ARTIFACTS = ROOT / "out" / "browser-acceptance" / "decision-review"


def two_site_fixture(path: Path) -> None:
    """Two displaced copies of 6QWR model 1 for UI mechanics, not biology."""
    lines = (ROOT / "config" / "policies" / "source-assets" / "6QWR.pdb").read_text().splitlines()
    inside = False
    atoms = []
    for line in lines:
        if line.startswith("MODEL"):
            inside = True
            continue
        if inside and line.startswith("ENDMDL"):
            break
        if inside and line.startswith("ATOM"):
            atoms.append(line)
    assert sum(line[17:20] == "HIS" for line in atoms) > 0
    with path.open("w") as output:
        output.write("REMARK 900 TWO-SITE SOFTWARE FIXTURE; NOT SCIENTIFIC SUITABILITY EVIDENCE\n")
        serial = 0
        for chain, shift in (("A", 0), ("B", 80)):
            for line in atoms:
                serial += 1
                x = float(line[30:38]) + shift
                output.write(f"{line[:6]}{serial:5d}{line[11:21]}{chain}{line[22:30]}"
                             f"{x:8.3f}{line[38:]}\n")
            output.write("TER\n")
        output.write("END\n")


def density_state(base: dict, standing: str) -> dict:
    state = deepcopy(base)
    protein = state["protein"]
    review = state["preparationReview"]
    assert protein and review and state["inspection"]
    intended_id = review["intendedProteinId"]
    state["inspection"].update({"subjectId": intended_id, "representationKind": "intendedProtein",
                                "evidence": [], "annotations": [], "metrics": [],
                                "focusId": None, "focus": None})
    state["notices"] = [{"id": "fixture-density", "severity": "warning", "subjectId": None,
                         "message": "Presentation fixture: 38 synthetic site choices. Viewer coordinates do not establish these synthetic sites or scientific suitability."}]
    changes = []
    decisions = []
    for number in range(1, 39):
        address = {"model": 0, "chain": "A", "residue": 500 + number,
                   "insertionCode": "", "copyId": "A"}
        state["sourceModels"][0]["residues"].append({
            "address": dict(address, copyId=""), "name": "HIS", "residueKind": "protein",
            "backboneHeavyAtomsComplete": True, "alternateLocations": []})
        chosen = number == 1 and standing == "confirmed" or standing in ("preparing", "failed")
        blocked = number == 1 and standing == "blocked"
        options = []
        for variant in ("HID", "HIE", "HIP"):
            proposal_id = f"fixture-site-{number:02d}-{variant}"
            changes.append({"id": proposal_id, "studyRevisionId": review["studyRevisionId"],
                            "intendedProteinId": intended_id, "residue": address,
                            "kind": "residueState", "proposedChange": variant,
                            "rationale": "Fixture alternative", "limitations": [],
                            "approvalRequired": True, "partnerResidue": None})
            disposition = ("confirmed" if variant == "HID" else "notChosen") if chosen else "available"
            options.append({"proposalId": proposal_id, "proposedChange": variant,
                            "disposition": disposition,
                            "decisionId": f"fixture-decision-{number}" if variant == "HID" and chosen else None,
                            "evidenceCurrent": False, "confirmationBlocker": "Fixture evidence is not current.",
                            "startsPreparationOnConfirmation": False,
                            "startsPreparationOnDecline": False, "evidence": []})
        decisions.append({"id": f"fixture-site-{number:02d}", "kind": "residueState",
                          "residue": address, "partnerResidue": None, "atomName": None,
                          "standing": "blocked" if blocked else "confirmed" if chosen else "pending",
                          "chosenProposalId": options[0]["proposalId"] if chosen else None,
                          "options": options,
                          "blocker": "Fixture observation unavailable at this site." if blocked else None})
    protein["changes"] = changes
    protein["status"] = "review"
    protein["summary"] = "Presentation fixture; no scientific result is established."
    confirmed = sum(item["standing"] == "confirmed" for item in decisions)
    review.update({"decisions": decisions, "confirmedCount": confirmed,
                   "remainingCount": 38 - confirmed,
                   "blockers": ["Fixture observation unavailable at this site."] if standing == "blocked" else [],
                   "preparationStanding": standing if standing in ("preparing", "failed", "blocked") else "awaitingDecisions",
                   "preparationMessage": "Fixture: applying choices and checking the candidate…" if standing == "preparing"
                   else "Fixture provider did not return an assessed candidate; choices remain recorded." if standing == "failed" else None})
    state["proteinTask"].update({
        "standing": standing if standing in ("preparing", "failed", "blocked") else "awaitingDecisions",
        "message": "Fixture: applying choices and checking the candidate…" if standing == "preparing"
                   else "Fixture provider did not return an assessed candidate; choices remain recorded." if standing == "failed"
                   else "Fixture observation unavailable at this site." if standing == "blocked"
                   else f"{38 - confirmed} presentation choices remain.",
        "preparedProteinId": None, "retryAvailable": standing == "failed"})
    if standing == "failed":
        retry = next((item for item in state["actions"] if item["kind"] == "retryProteinPreparation"), None)
        if retry is None:
            state["actions"].append({"kind": "retryProteinPreparation", "subjectId": None,
                                     "enabled": True, "reason": None})
        else:
            retry.update({"enabled": True, "reason": None})
    return state


def distinct_meanings_state(base: dict, declined_bond: bool = False) -> dict:
    """Presentation-only options with different decision meanings and counts."""
    state = deepcopy(base)
    review = state["preparationReview"]
    protein = state["protein"]
    assert review and protein
    choices = [
        ("residueState", 901, "ASP", ["ASP", "ASH"], None),
        ("alternateLocation", 902, "SER", ["Conformer A", "Conformer B"], None),
        ("heavyAtom", 903, "ALA", ["Complete sidechain CB"], None),
        ("disulfide", 904, "CYS", ["Form CYS 904–905 bond"], 905),
    ]
    decisions = []
    changes = []
    for kind, number, residue_name, variants, partner_number in choices:
        address = {"model": 0, "chain": "A", "residue": number,
                   "insertionCode": "", "copyId": "A"}
        partner = dict(address, residue=partner_number) if partner_number else None
        for site_number, name in ((number, residue_name), (partner_number, "CYS")):
            if site_number is not None:
                state["sourceModels"][0]["residues"].append({
                    "address": dict(address, residue=site_number, copyId=""),
                    "name": name, "residueKind": "protein",
                    "backboneHeavyAtomsComplete": True,
                    "alternateLocations": ["A", "B"] if number == 902 else []})
        options = []
        for index, variant in enumerate(variants):
            proposal_id = f"fixture-{kind}-{number}-{index}"
            declined = kind == "disulfide" and declined_bond
            changes.append({"id": proposal_id, "studyRevisionId": review["studyRevisionId"],
                            "intendedProteinId": review["intendedProteinId"],
                            "residue": address, "partnerResidue": partner,
                            "kind": kind, "proposedChange": variant,
                            "rationale": "Presentation fixture; no molecular state is established.",
                            "limitations": [], "approvalRequired": True})
            options.append({"proposalId": proposal_id, "proposedChange": variant,
                            "disposition": "declined" if declined else "available",
                            "decisionId": f"fixture-declined-{number}" if declined else None,
                            "evidenceCurrent": False,
                            "confirmationBlocker": "Presentation fixture has no owner-verified confirmation evidence.",
                            "startsPreparationOnConfirmation": False,
                            "startsPreparationOnDecline": False,
                            "evidence": [{"id": f"fixture-evidence-{number}-{index}",
                                          "subjectId": proposal_id, "source": "Presentation fixture",
                                          "method": "Synthetic option wording",
                                          "observation": f"Option {variant} belongs only to this fixture.",
                                          "applicability": "UI mechanics only",
                                          "uncertainty": "No scientific suitability is established.",
                                          "bearing": "Unknown"}]})
        decisions.append({"id": f"fixture-{kind}-{number}", "kind": kind,
                          "residue": address, "partnerResidue": partner,
                          "atomName": "CB" if kind == "heavyAtom" else None,
                          "standing": "confirmed" if kind == "disulfide" and declined_bond else "pending",
                          "chosenProposalId": options[0]["proposalId"] if kind == "disulfide" and declined_bond else None,
                          "options": options, "blocker": None})
    protein.update(status="review", changes=changes,
                   findings=[{"id": "fixture-finding", "subjectId": review["intendedProteinId"],
                              "evidenceId": "fixture-finding-evidence",
                              "meaning": "Presentation fixture: informational geometry note",
                              "consequence": "Read this finding; it is not an approval obligation.",
                              "disposition": "Informational", "material": False}],
                   summary="Presentation fixture: synthetic decisions; displayed coordinates are not evidence for these sites.")
    review.update(decisions=decisions, confirmedCount=1 if declined_bond else 0,
                  remainingCount=3 if declined_bond else 4, blockers=[],
                  preparationStanding="awaitingDecisions", preparationMessage=None)
    state["notices"] = [{"id": "fixture-meanings", "severity": "warning", "subjectId": None,
                         "message": "Presentation fixture only; no scientific result is established."}]
    return state


class DecisionReviewBrowserTests(unittest.TestCase):
    def test_optional_focus_failure_and_late_response_do_not_change_the_choice(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-histidine-software-fixture.pdb"
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
                    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
                    page.get_by_label("Chain A").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    expect(page.locator(".review-progress")).to_contain_text("1 site choice", timeout=120000)
                    review = page.request.get(base + "/api/state").json()["preparationReview"]
                    site = next(item for item in review["decisions"] if item["kind"] == "residueState")
                    hid = next(item["proposalId"] for item in site["options"] if item["proposedChange"] == "HID")
                    hie = next(item["proposalId"] for item in site["options"] if item["proposedChange"] == "HIE")
                    page.locator(".review-option").filter(has_text="HID").locator("input").check()
                    confirm = page.get_by_role("button", name="Confirm state")
                    expect(confirm).to_be_enabled()
                    expect(page.locator(".viewer-mount[data-camera-ready='true']")).to_be_visible(timeout=120000)
                    before_camera = page.evaluate("""() => {
                      const viewer = document.querySelector('.viewer-mount')?.[Symbol.for('molstar.viewer')];
                      return viewer?.plugin.canvas3d?.camera.getSnapshot().target;
                    }""")

                    def refuse_evidence(route):
                        payload = route.request.post_data_json
                        if payload.get("kind") == "selectInspectionSubject" and payload.get("data", {}).get("subjectId") == hid:
                            route.fulfill(status=503, content_type="application/json",
                                          body=json.dumps({"reason": "Fixture evidence service unavailable."}))
                        else:
                            route.continue_()

                    page.route("**/api/commands", refuse_evidence)
                    page.locator(".review-option.picked").get_by_role("button", name="Focus on this residue").click()
                    expect(page.locator(".protein-review-panel .review-blocker-inline")).to_contain_text(
                        "Fixture evidence service unavailable.")
                    expect(confirm).to_be_enabled()
                    page.screenshot(path=str(ARTIFACTS / "real-evidence-failed-1024.png"), full_page=True)
                    page.unroute("**/api/commands", refuse_evidence)

                    page.evaluate("""id => {
                        window.__realReviewFetch = window.fetch.bind(window);
                        window.fetch = (...arguments_) => {
                            const options = arguments_[1];
                            if (String(arguments_[0]).endsWith('/api/commands') && options?.body) {
                                const request = JSON.parse(options.body);
                                if (request.kind === 'selectInspectionSubject' && request.data.subjectId === id)
                                    return new Promise(resolve => setTimeout(() =>
                                        resolve(window.__realReviewFetch(...arguments_)), 1200));
                            }
                            return window.__realReviewFetch(...arguments_);
                        };
                    }""", hid)
                    page.locator(".review-option.picked").get_by_role("button", name="Focus on this residue").click()
                    page.locator(".review-option").filter(has_text="HIE").locator("input").check()
                    expect(page.locator(".review-option.picked")).to_contain_text("HIE")
                    page.wait_for_timeout(1600)
                    expect(page.locator(".review-option.picked")).to_contain_text(
                        "Neutral histidine with its ring proton on NE2")
                    expect(page.get_by_label("Selected option explanation and evidence")).to_contain_text(
                        "has not measured which state is best")
                    expect(confirm).to_be_enabled()
                    self.assertEqual(page.request.get(base + "/api/state").json()["inspection"]["subjectId"], hid)
                    self.assertIsNone(page.request.get(base + "/api/state").json()["inspection"]["focusId"])
                    self.assertEqual(before_camera, page.evaluate("""() => {
                      const viewer = document.querySelector('.viewer-mount')?.[Symbol.for('molstar.viewer')];
                      return viewer?.plugin.canvas3d?.camera.getSnapshot().target;
                    }"""))
                    page.screenshot(path=str(ARTIFACTS / "real-stale-evidence-1024.png"), full_page=True)
                    page.evaluate("window.fetch = window.__realReviewFetch")
                    self.assertEqual(page.request.get(base + "/api/state").json()["preparationReview"]["confirmedCount"], 0)
                    confirm.click()
                    expect(page.locator(".review-receipt")).to_contain_text("HIE confirmed", timeout=120000)
                    self.assertEqual(page.request.get(base + "/api/state").json()["preparationReview"]["confirmedCount"], 1)
                finally:
                    browser.close()

    def test_two_site_continuation_and_labelled_density_states(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-histidine-software-fixture.pdb"
            two_site_fixture(source)
            with running_host(directory / "workspace") as base, sync_playwright() as playwright:
                browser = chromium(playwright)
                try:
                    page = browser.new_page(viewport={"width": 1672, "height": 950}, device_scale_factor=1)
                    requests = []
                    page.on("request", lambda request: requests.append(json.loads(request.post_data))
                            if request.url.endswith("/api/commands") and request.post_data else None)
                    page.goto(base)
                    expect(page.locator(".workspace-body")).to_be_visible(timeout=30000)
                    page.locator("#source-upload").set_input_files(str(source))
                    page.locator("#upload-provenance").select_option("experimental")
                    page.get_by_role("button", name="Upload source").click()
                    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=120000)
                    expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
                    page.get_by_label("Chain A").check()
                    page.get_by_label("Chain B").check()
                    page.get_by_role("button", name="Assess selected protein").click()
                    expect(page.locator(".review-progress")).to_contain_text("2 site choices", timeout=120000)
                    account = page.request.get(base + "/api/state").json()
                    review = account["preparationReview"]
                    self.assertEqual((len(account["protein"]["changes"]), len(review["decisions"])), (8, 4))
                    site_a = next(item for item in review["decisions"] if item["kind"] == "residueState" and item["residue"]["chain"] == "A")
                    site_b = next(item for item in review["decisions"] if item["kind"] == "residueState" and item["residue"]["chain"] == "B")
                    hid = next(item["proposalId"] for item in site_a["options"] if item["proposedChange"] == "HID")
                    hie = next(item["proposalId"] for item in site_b["options"] if item["proposedChange"] == "HIE")
                    page.screenshot(path=str(ARTIFACTS / "real-two-site-pending-1672.png"), full_page=True)

                    page.locator(".review-site-list .review-site").first.click()
                    page.locator(".review-option").filter(has_text="HID").locator("input").check()
                    confirm = page.get_by_role("button", name="Confirm state")
                    expect(confirm).to_be_enabled()
                    page.evaluate("window.__reviewViewer = document.querySelector('.viewer-mount')?.[Symbol.for('molstar.viewer')]")
                    confirm.click()
                    expect(page.locator(".review-receipt")).to_contain_text("HID confirmed", timeout=30000)
                    page.screenshot(path=str(ARTIFACTS / "real-two-site-confirmed-1672.png"), full_page=True)
                    review = page.request.get(base + "/api/state").json()["preparationReview"]
                    self.assertEqual((review["confirmedCount"], review["remainingCount"]), (1, 3))
                    self.assertEqual(sum(option["disposition"] == "notChosen" for option in
                                         next(item for item in review["decisions"] if item["id"] == site_a["id"])["options"]), 2)

                    page.get_by_role("button", name="Next unresolved site").click()
                    expect(page.locator(".review-site-title")).to_contain_text("Histidine 108 · Chain B")
                    radio = page.locator(".review-option").filter(has_text="HIE").locator("input")
                    radio.focus()
                    radio.press("Space")
                    expect(confirm).to_be_enabled()
                    self.assertTrue(page.evaluate("window.__reviewViewer === document.querySelector('.viewer-mount')?.[Symbol.for('molstar.viewer')]"))
                    confirm.press("Enter")
                    expect(page.locator(".review-receipt")).to_contain_text("HIE confirmed", timeout=30000)
                    self.assertEqual([(item["kind"], item["data"].get("proposalId")) for item in requests
                                      if item["kind"] == "approvePreparationChange"],
                                     [("approvePreparationChange", hid), ("approvePreparationChange", hie)])
                    self.assertEqual([item["data"]["subjectId"] for item in requests
                                      if item["kind"] == "selectInspectionSubject"], [])

                    page.get_by_role("button", name="Next unresolved repair").click()
                    expect(page.locator(".review-site-scope")).to_contain_text("Missing-atom repair")
                    expect(page.locator(".review-other")).to_have_attribute("open", "")
                    page.locator(".review-site-list .review-site").first.click()

                    slow_page = browser.new_page(viewport={"width": 1024, "height": 768})
                    slow_page.route("**/api/structures/**", lambda route: (time.sleep(1.2), route.continue_()))
                    slow_page.goto(base, wait_until="domcontentloaded")
                    expect(slow_page.locator(".scene-loading")).to_be_visible(timeout=5000)
                    slow_page.screenshot(path=str(ARTIFACTS / "real-two-site-loading-1024.png"), full_page=True)
                    expect(slow_page.locator(".viewer-mount[data-camera-ready='true']")).to_be_visible(timeout=120000)
                    slow_page.close()

                    failed_page = browser.new_page(viewport={"width": 1024, "height": 768})
                    failed_page.route("**/api/structures/**", lambda route: route.abort("failed"))
                    failed_page.goto(base)
                    expect(failed_page.locator(".scene-error")).to_contain_text("Structure unavailable", timeout=30000)
                    failed_page.screenshot(path=str(ARTIFACTS / "real-two-site-visualization-failed-1024.png"), full_page=True)
                    failed_page.locator(".review-site-list .review-site").first.click()
                    expect(failed_page.locator(".review-site-title")).to_contain_text("Histidine 108 · Chain A")
                    failed_page.close()

                    page.get_by_role("button", name="Reviewed sites (2)").click()
                    self.assertEqual(page.locator(".review-site-list .review-site").count(), 2)
                    page.get_by_role("button", name="Open sites (0)").click()
                    self.assertEqual(page.locator(".review-site-list .review-site").count(), 0)
                    page.get_by_role("button", name="All sites (2)").click()
                    page.locator(".review-site-list .review-site").first.press("Enter")
                    expect(page.locator(".review-receipt")).to_contain_text("HID confirmed")
                    page.set_viewport_size({"width": 1024, "height": 768})
                    page.screenshot(path=str(ARTIFACTS / "real-two-site-confirmed-1024.png"), full_page=True)
                    page.set_viewport_size({"width": 820, "height": 720})
                    page.screenshot(path=str(ARTIFACTS / "real-two-site-confirmed-820.png"), full_page=True)
                    self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth - innerWidth"), 0)
                    page.get_by_role("button", name="Collapse inputs").click()
                    page.screenshot(path=str(ARTIFACTS / "real-two-site-collapsed-820.png"), full_page=True)
                    page.get_by_role("button", name="Show inputs").click()
                    expect(page.locator(".review-receipt")).to_contain_text("HID confirmed")

                    base_account = page.request.get(base + "/api/state").json()
                    fixture = {"value": density_state(base_account, "pending")}
                    page.route("**/api/state", lambda route: route.fulfill(
                        status=200, content_type="application/json", body=json.dumps(fixture["value"])))
                    for standing in ("pending", "confirmed", "blocked", "preparing", "failed"):
                        fixture["value"] = density_state(base_account, standing)
                        for width, height in ((1672, 950), (1024, 768), (820, 720)):
                            page.set_viewport_size({"width": width, "height": height})
                            page.reload()
                            if standing in ("preparing", "failed"):
                                expect(page.get_by_role("region", name="Protein task outcome"))\
                                    .to_contain_text("Fixture", timeout=30000)
                                expect(page.locator(".review-progress")).to_have_count(0)
                            else:
                                expect(page.locator(".review-progress"))\
                                    .to_contain_text("38 site choices", timeout=30000)
                            if standing == "confirmed":
                                page.locator(".review-site-list .review-site").first.click()
                                expect(page.locator(".review-receipt")).to_contain_text("HID confirmed")
                            if standing == "blocked":
                                page.locator(".review-site-list .review-site").first.click()
                                expect(page.locator(".review-blocker")).to_be_visible()
                            if standing == "failed":
                                expect(page.get_by_role("button", name="Retry protein preparation")).to_be_enabled()
                            self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth - innerWidth"), 0)
                            page.screenshot(path=str(ARTIFACTS / f"fixture-38-{standing}-{width}.png"), full_page=True)
                            if standing == "pending" and width == 1672:
                                page.locator(".review-site-list .review-site").last.press("Enter")
                                expect(page.locator(".review-site-title")).to_contain_text("Histidine 538 · Chain A")
                                page.screenshot(path=str(ARTIFACTS / "fixture-38-last-site-1672.png"), full_page=True)
                    page.get_by_role("button", name="Review confirmed choices (38)").click()
                    self.assertEqual(page.locator(".review-site-list .review-site").count(), 38)

                    fixture["value"] = distinct_meanings_state(base_account)
                    page.set_viewport_size({"width": 1024, "height": 768})
                    page.reload()
                    expect(page.locator(".review-progress")).to_contain_text("2 site choices")
                    page.locator(".review-site-list .review-site").filter(has_text="Aspartate 901 · Chain A").click()
                    self.assertEqual(page.locator(".review-option").count(), 2)
                    page.locator(".review-option").filter(has_text="ASH").locator("input").check()
                    expect(page.locator(".review-option.picked"))\
                        .to_contain_text("Neutral, protonated aspartate")
                    self.assertEqual(page.get_by_role("button", name="Reject possible bond").count(), 0)
                    self.assertEqual(page.get_by_role("button", name="Reject required repair").count(), 0)
                    page.locator(".review-site-list .review-site").filter(has_text="Serine 902 · Chain A").click()
                    self.assertEqual(page.locator(".review-option").count(), 2)
                    page.locator(".review-option").filter(has_text="Conformer B").locator("input").check()
                    expect(page.get_by_label("Selected option explanation and evidence"))\
                        .to_contain_text("Option Conformer B belongs only to this fixture")
                    repairs = page.locator(".review-other").filter(has_text="Repairs and other findings")
                    repairs.locator("summary").click()
                    repairs.locator(".review-site").filter(has_text="ALA 903 · Chain A").click()
                    expect(page.get_by_role("button", name="Approve repair")).to_be_visible()
                    expect(page.get_by_role("button", name="Reject required repair")).to_be_visible()
                    repairs.locator(".review-site").filter(has_text="Cysteine 904 · Chain A").click()
                    expect(page.get_by_role("button", name="Confirm possible bond")).to_be_visible()
                    expect(page.get_by_role("button", name="Reject possible bond")).to_be_visible()
                    expect(page.locator(".review-site-scope")).to_contain_text("Bond partner: Cysteine 905, chain A")
                    self.assertEqual(page.locator(".review-other").filter(
                        has_text="informational geometry note").locator("button").count(), 0)
                    page.screenshot(path=str(ARTIFACTS / "fixture-distinct-decision-types-1024.png"), full_page=True)

                    fixture["value"] = distinct_meanings_state(base_account, declined_bond=True)
                    page.set_viewport_size({"width": 820, "height": 720})
                    page.reload()
                    repairs = page.locator(".review-other").filter(has_text="Repairs and other findings")
                    repairs.locator("summary").click()
                    repairs.locator(".review-site").filter(has_text="Cysteine 904 · Chain A").click()
                    expect(page.locator(".review-receipt")).to_contain_text("No bond confirmed for this pair")
                    self.assertEqual(page.get_by_role("button", name="Confirm possible bond").count(), 0)
                    self.assertEqual(page.get_by_role("button", name="Reject possible bond").count(), 0)
                    page.screenshot(path=str(ARTIFACTS / "fixture-declined-bond-820.png"), full_page=True)

                    no_choice = deepcopy(base_account)
                    no_choice["protein"]["changes"] = []
                    no_choice["protein"]["status"] = "review"
                    no_choice["inspection"].update({
                        "subjectId": no_choice["preparationReview"]["intendedProteinId"],
                        "representationKind": "intendedProtein", "evidence": [],
                        "annotations": [], "metrics": [], "focusId": None, "focus": None})
                    no_choice["preparationReview"].update({
                        "decisions": [], "confirmedCount": 0, "remainingCount": 0,
                        "blockers": [], "preparationStanding": "failed",
                        "preparationMessage": "Fixture: a no-choice preparation returned no assessed candidate."})
                    no_choice["proteinTask"].update({
                        "standing": "failed", "message": "Fixture: a no-choice preparation returned no assessed candidate.",
                        "preparedProteinId": None, "retryAvailable": True})
                    retry = next((item for item in no_choice["actions"]
                                  if item["kind"] == "retryProteinPreparation"), None)
                    if retry is None:
                        no_choice["actions"].append({"kind": "retryProteinPreparation", "subjectId": None,
                                                      "enabled": True, "reason": None})
                    else:
                        retry.update({"enabled": True, "reason": None})
                    no_choice["notices"] = [{"id": "fixture-no-choice", "severity": "warning",
                                              "subjectId": None,
                                              "message": "Presentation fixture: no-choice provider failure."}]
                    fixture["value"] = no_choice
                    page.set_viewport_size({"width": 820, "height": 720})
                    page.reload()
                    expect(page.get_by_role("button", name="Retry protein preparation")).to_be_enabled()
                    expect(page.locator(".review-progress")).to_have_count(0)
                    expect(page.locator(".viewer-mount[data-camera-ready='true']")).to_be_visible(timeout=120000)
                    page.screenshot(path=str(ARTIFACTS / "fixture-no-choice-failed-820.png"), full_page=True)
                finally:
                    browser.close()


if __name__ == "__main__":
    unittest.main()
