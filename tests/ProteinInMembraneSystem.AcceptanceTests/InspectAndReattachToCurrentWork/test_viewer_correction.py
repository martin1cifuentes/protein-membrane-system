"""Focused browser checks for the active Preparation subject and placement guide.

Controlled accounts exercise presentation without starting a scientific worker.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import threading
import unittest

from playwright.sync_api import expect, sync_playwright


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "InspectAndReattach"))
from test_browser_presentation import (  # noqa: E402
    CHROMIUM, DIST, account, completed, controlled_account_server, inspection,
)

VISUAL_EVIDENCE = Path(__file__).resolve().parents[3] / "out" / "ui-correction-visual"
TRPCAGE_SOURCE = (Path(__file__).resolve().parents[1] /
                  "ConstructAndMinimizeExplicitSystem" / "fixtures" / "1L2Y-RCSB.pdb")
SAVED_AB_87 = (Path(__file__).resolve().parents[3] / "out" / "workspace" / "placement" /
               "77cc43dcf04c429b87f0cf2bb1f66559" / "oriented-protein.pdb")


def positioned_account(position: str = "upper", z_translation: float = 29.0,
                       x_translation: float = 32.0) -> dict:
    value = account()
    value["placement"] = {
        "proposalId": "placement-one", "status": "supported",
        "preparedProteinId": "protein-one", "membraneModelId": "membrane-one",
        "topologyKind": "membrane-spanning", "physicalSide": position,
        "midplaneAngstrom": 0.0, "thicknessAngstrom": 46.0,
        "depthAngstrom": 30.0, "tiltDegrees": 0.0, "sidedness": None,
        "contactingRegions": [], "limitations": [], "policyId": None,
        "policyVersion": None, "witnessId": None,
        "reason": "Exact construct and placement frame are technically ready.",
        "evidence": [], "prediction": None,
        "transform": {"startingPosition": position, "offsetXAngstrom": x_translation,
                      "offsetYAngstrom": 0.0, "offsetZAngstrom": 0.0,
                      "rotationXDegrees": 0.0, "rotationYDegrees": 0.0,
                      "rotationZDegrees": 0.0, "appliedTranslationXAngstrom": x_translation,
                      "appliedTranslationYAngstrom": 0.0,
                      "appliedTranslationZAngstrom": z_translation,
                      "headgroupBoundaryAngstrom": 23.0},
    }
    value["inspection"] = inspection("placement-one", "oriented-protein-with-proposed-membrane-bounds")
    value["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
    value["inspection"]["omittedMolecules"] = []
    return value


def positioned_protein_bytes(z_translation: float = 29.0) -> bytes:
    atoms = (("N", 31.0, 0.0, z_translation), ("CA", 32.0, 0.0, z_translation + 1),
             ("C", 33.0, 0.0, z_translation + 2), ("O", 34.0, 0.0, z_translation + 3),
             ("N", 32.0, 1.0, z_translation + 4), ("CA", 33.0, 1.0, z_translation + 5),
             ("C", 34.0, 1.0, z_translation + 6), ("O", 35.0, 1.0, z_translation + 7))
    lines = [f"ATOM  {index:5d} {name:<4s} ALA A{1 if index <= 4 else 2:4d}    "
             f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}\n"
             for index, (name, x, y, z) in enumerate(atoms, 1)]
    return ("".join(lines) + "END\n").encode()


def trpcage_first_model_bytes() -> bytes:
    """Existing 1L2Y source fixture, first deposited NMR coordinate model only."""
    atoms: list[str] = []
    in_first_model = False
    for line in TRPCAGE_SOURCE.read_text().splitlines():
        if line.startswith("MODEL "):
            if in_first_model:
                break
            in_first_model = True
        elif in_first_model and line.startswith("ENDMDL"):
            break
        elif in_first_model and line.startswith(("ATOM  ", "HETATM")):
            atoms.append(line)
    if len(atoms) < 200:
        raise AssertionError("The first deposited 1L2Y model was unavailable")
    return ("\n".join(atoms) + "\nEND\n").encode()


class ViewerCorrectionTests(unittest.TestCase):
    def test_opm_no_match_keeps_earlier_pose_and_reports_independent_ppm_availability(self):
        """Controlled owner outcome; no external OPM search is claimed."""
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account("center", 0.0, 0.0)
            value["attempt"] = None
            value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                "summary": "Controlled prepared protein", "atomCount": 8,
                                "changes": [], "findings": [], "prediction": None,
                                "geometry": None, "sourceGeometry": None}
            value["membrane"] = {"modelId": "membrane-one", "status": "assessed",
                                 "upper": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "lower": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "scientificPurpose": None, "limitations": [], "reason": None,
                                 "policyId": None, "policyVersion": "2.0.0",
                                 "evidence": [], "speciesSupport": []}
            value["placement"]["policyVersion"] = "2.0.0"
            value["placementTask"] = {"studyRevisionId": value["study"]["id"],
                                      "preparedProteinId": "protein-one", "membraneModelId": "membrane-one",
                                      "proposalId": "placement-one", "standing": "supported",
                                      "message": "Earlier checked position", "adopted": False,
                                      "latestAttemptIssue": None, "routeOutcome": None}
            value["placementMethods"] = [
                {"method": "OPM", "standing": "lookupEligible", "reason": None},
                {"method": "PPM", "standing": "configured", "reason": None},
            ]
            value["actions"].extend([
                {"kind": "proposePlacement", "subjectId": None, "enabled": True, "reason": None},
                {"kind": "adoptPlacement", "subjectId": None, "enabled": True, "reason": None},
            ])
            host.structure_bytes = positioned_protein_bytes(0.0)
            host.replace(value)
            received = threading.Event()
            release = threading.Event()
            request_number = [0]

            def handle_command(handler, command):
                if command["kind"] == "adoptPlacement":
                    self.assertEqual(command["data"]["proposalId"], "placement-one")
                    selected = deepcopy(host.account)
                    selected["revision"] += 1
                    selected["study"]["id"] = "revision-adopted"
                    selected["study"]["number"] += 1
                    selected["study"]["adoptedPlacementProposalId"] = "placement-one"
                    selected["placementTask"]["studyRevisionId"] = "revision-adopted"
                    selected["placementTask"]["adopted"] = True
                    selected["inspection"]["studyRevisionId"] = "revision-adopted"
                    selected["inspection"]["studyRevisionNumber"] += 1
                    host.replace(selected)
                    handler.json_response(selected)
                    return True
                if command["kind"] != "proposePlacement":
                    return False
                self.assertEqual(command["data"]["orientationRoute"], "opm")
                self.assertEqual(command["data"]["physicalSide"], "both")
                request_number[0] += 1
                searching = deepcopy(host.account)
                searching["revision"] += 1
                searching["placementTask"]["standing"] = "obtaining"
                searching["placementTask"]["routeOutcome"] = {
                    "requestId": f"controlled-opm-request-{request_number[0]}",
                    "studyRevisionId": searching["study"]["id"],
                    "preparedProteinId": "protein-one", "membraneModelId": "membrane-one",
                    "route": "opm", "topologyKind": "membrane-spanning", "physicalSide": "both",
                    "ppmNterminalSide": None, "standing": "searching",
                    "message": "Searching OPM for the exact selected source…", "proposalId": None}
                host.replace(searching)
                received.set()
                if not release.wait(timeout=15):
                    handler.json_response({"reason": "Controlled OPM reply was not released"}, 500)
                    return True
                no_match = deepcopy(searching)
                no_match["revision"] += 1
                no_match["placementTask"]["standing"] = "noProposal"
                no_match["placementTask"]["routeOutcome"].update(
                    standing="noMatch", message="No matching OPM reference was found for this source.")
                host.replace(no_match)
                handler.json_response(no_match)
                return True

            host.command_handler = handle_command
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 940})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Placement").click()
                page.locator("#orientation-route").select_option("opm")
                topology = page.locator("#topology-kind")
                self.assertFalse(topology.locator("option[value='']").is_enabled())
                topology.select_option("membrane-spanning")
                physical = page.locator("#physical-side")
                self.assertEqual(physical.input_value(), "", "A physical side must remain unselected")
                self.assertFalse(physical.locator("option[value='']").is_enabled())
                physical.select_option("both")
                page.get_by_role("button", name="Calculate orientation estimate").click()
                self.assertTrue(received.wait(timeout=10))
                expect(page.locator(".method-route-outcome")).to_contain_text("Searching OPM")
                release.set()
                expect(page.locator(".method-route-outcome")).to_contain_text(
                    "No matching OPM reference found")
                expect(page.locator(".method-route-outcome")).to_contain_text(
                    "separately available local PPM method")
                expect(page.get_by_label("Current placement status")).to_contain_text("No position from OPM")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_enabled()
                expect(page.get_by_text("Placement assessed.")).to_have_count(0)
                expect(page.get_by_label("Placement proposal")).to_contain_text("Earlier checked position")
                expect(page.get_by_text("Earlier position checks")).to_have_count(1)
                self.assertEqual(host.account["placement"]["proposalId"], "placement-one")
                expect(page.get_by_text("Position check version 2.0.0")).to_have_count(0)
                output = VISUAL_EVIDENCE / "feedback-continuation"
                output.mkdir(parents=True, exist_ok=True)
                for width, height in ((1672, 940), (1024, 768), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                        .to_have_count(1, timeout=60000)
                    self.assert_framed(page, 0.0)
                    page.wait_for_timeout(180)
                    page.screenshot(path=str(output / f"opm-no-match-prior-pose-{width}.png"))
                unavailable = deepcopy(host.account)
                unavailable["revision"] += 1
                unavailable["placementMethods"][1] = {"method": "PPM", "standing": "unavailable",
                                                      "reason": "Local PPM is not verified."}
                host.replace(unavailable)
                expect(page.locator(".method-route-outcome")).not_to_contain_text(
                    "separately available local PPM method")
                expect(page.locator(".method-route-outcome")).to_contain_text("direct positioning")
                page.get_by_role("button", name="Use earlier checked position").click()
                expect(page.get_by_label("Current placement status")).to_contain_text("Position selected")
                expect(page.get_by_role("button", name="Review preparation requirements")).to_be_visible()
                self.assertEqual(host.account["study"]["adoptedPlacementProposalId"], "placement-one")
                page.get_by_role("button", name="Calculate orientation estimate").click()
                expect(page.locator(".method-route-outcome")).to_contain_text(
                    "No matching OPM reference found")
                expect(page.get_by_label("Current placement status")).to_contain_text("No position from OPM")
                expect(page.get_by_label("Current placement status")).to_contain_text("Position selected")
                expect(page.get_by_role("button", name="Review preparation requirements")).to_be_visible()
                for width, height in ((1672, 940), (1024, 768), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    page.get_by_label("Current placement status").scroll_into_view_if_needed()
                    self.assert_framed(page, 0.0)
                    page.wait_for_timeout(180)
                    page.screenshot(path=str(output / f"opm-no-match-selected-pose-{width}.png"))
                received.clear()
                release.clear()
                page.get_by_role("button", name="Calculate orientation estimate").click()
                self.assertTrue(received.wait(timeout=10))
                page.locator("#topology-kind").select_option("one-surface-associated")
                release.set()
                expect(page.get_by_text("The earlier orientation request used different method inputs.")).to_be_visible()
                expect(page.locator(".method-route-outcome")).to_have_count(0)
                expect(page.get_by_label("Current placement status")).to_contain_text("Position selected")
                expect(page.get_by_role("button", name="Review preparation requirements")).to_be_visible()
            finally:
                release.set()
                browser.close()

    @staticmethod
    def projected_guide_and_protein(page, z_translation: float) -> dict:
        return page.evaluate("""z => {
          const mount = document.querySelector('.placement-scene-surface .viewer-mount');
          const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
          const camera = viewer.plugin.canvas3d.camera;
          const points = [
            [-18, -18, -23], [-18, 18, -23], [18, -18, -23], [18, 18, -23],
            [-18, -18, 23], [-18, 18, 23], [18, -18, 23], [18, 18, 23],
            [31, 0, z], [35, 1, z + 7]
          ];
          return {width: camera.viewport.width, height: camera.viewport.height,
            pixels: points.map(point => Array.from(camera.project(
              new Float32Array(4), new Float32Array(point))).slice(0, 2))};
        }""", z_translation)

    def assert_framed(self, page, z_translation: float):
        projected = self.projected_guide_and_protein(page, z_translation)
        self.assertGreater(projected["width"], 0)
        self.assertGreater(projected["height"], 0)
        for x, y in projected["pixels"]:
            self.assertGreater(x, projected["width"] * 0.02)
            self.assertLess(x, projected["width"] * 0.98)
            self.assertGreater(y, projected["height"] * 0.02)
            self.assertLess(y, projected["height"] * 0.98)

    def test_preparation_has_no_protein_fallback_before_checked_construction(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            pending = positioned_account()
            pending["attempt"]["constructed"] = None
            pending["attempt"].update(stageKind="Construction", phase="providerPacking",
                                      message="Packing is in progress")
            host.replace(pending)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.get_by_label("System preparation view")).to_contain_text("Building system")
                expect(page.locator(".viewer-mount")).to_have_count(0)
                expect(page.get_by_role("button", name="Show whole system")).to_have_count(0)
                self.assertNotIn("selectInspectionSubject", host.commands)

                checked = deepcopy(pending)
                checked["revision"] += 1
                checked["attempt"]["constructed"] = account()["attempt"]["constructed"]
                checked["attempt"].update(stageKind="Minimization", phase="finalMinimization",
                                           message="Final minimization is running")
                checked["inspection"] = inspection("constructed-one", "constructedSystem")
                checked["inspection"]["structureUrl"] = "/api/structures/tiny?format=pdb"
                checked["inspection"]["omittedMolecules"] = []
                host.replace(checked)
                expect(page.get_by_label("Current attempt inspection").locator(".viewer-mount canvas"))\
                    .to_have_count(1, timeout=60000)
                expect(page.get_by_role("button", name="Show whole system")).to_be_enabled()
                expect(page.get_by_label("Visible molecular components")).to_have_count(0)
                page.get_by_role("button", name="Components").click()
                expect(page.get_by_label("Visible molecular components")).to_contain_text(
                    "Upper physical leaflet: DMPC 2")
                expect(page.get_by_label("Visible molecular components")).to_contain_text(
                    "Lower physical leaflet: DMPC 2")
                expect(page.locator(".stage-strip")).to_contain_text(
                    "No minimized or equilibrated system stage has completed.")
                self.assertNotIn("selectInspectionSubject", host.commands)
            finally:
                browser.close()

    def test_completed_stage_does_not_offer_measure_without_a_visible_target(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = completed(account())
            value["inspection"] = inspection("stage-one", "completedStage", value["stages"][0]["assessment"])
            value["inspection"].update(structureUrl="/api/structures/tiny?format=pdb", omittedMolecules=[])
            value["stages"][0]["observation"] = {
                "termination": "converged", "providerVersion": "OpenMM 8.6",
                "measurements": [{"name": "finalRmsForce", "value": 1.23,
                                  "unit": "kJ/mol/nm", "scope": "final stage"}],
                "localState": None, "proteinGeometry": None,
            }
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_role("button", name="Results").click()
                expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=60000)
                expect(page.get_by_role("button", name="Measure", exact=True)).to_have_count(0)
                observations = page.get_by_label("Stage-specific findings and evidence")
                observations.get_by_text("Observed final operation and measurements").click()
                expect(observations).to_contain_text("Final RMS force")
                expect(observations).to_contain_text("1.23")
            finally:
                browser.close()

    def test_complete_placement_draft_shows_actual_check_blocker(self):
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account("center", 0.0, 0.0)
            value["attempt"] = None
            value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                "summary": "Controlled prepared protein", "atomCount": 8,
                                "changes": [], "findings": [], "prediction": None,
                                "geometry": None, "sourceGeometry": None}
            value["membrane"] = {"modelId": "membrane-one", "status": "assessed",
                                 "upper": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "lower": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "scientificPurpose": None, "limitations": [], "reason": None,
                                 "policyId": None, "policyVersion": None,
                                 "evidence": [], "speciesSupport": []}
            value["actions"].append({"kind": "proposePlacement", "subjectId": None,
                                      "enabled": False, "reason": "A prepared protein is required before checking."})
            host.structure_bytes = positioned_protein_bytes(0.0)
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                page.get_by_label("Move X (Å)").fill("1")
                status = page.get_by_label("Current placement status")
                expect(status).to_contain_text("Check unavailable")
                expect(status).to_contain_text("A prepared protein is required before checking.")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_disabled()
                page.wait_for_timeout(650)
                self.assertNotIn("proposePlacement", host.commands)
            finally:
                browser.close()

    def test_late_position_check_cannot_ready_an_edited_draft(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account("center", 0.0, 0.0)
            value["attempt"] = None
            value["protein"] = {"subjectId": "protein-one", "status": "assessed",
                                "summary": "Controlled prepared protein", "atomCount": 8,
                                "changes": [], "findings": [], "prediction": None,
                                "geometry": None, "sourceGeometry": None}
            value["membrane"] = {"modelId": "membrane-one", "status": "assessed",
                                 "upper": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "lower": [{"speciesId": "DMPC", "fraction": 1.0}],
                                 "scientificPurpose": None, "limitations": [], "reason": None,
                                 "policyId": None, "policyVersion": None,
                                 "evidence": [], "speciesSupport": []}
            value["actions"].extend([
                {"kind": "proposePlacement", "subjectId": None, "enabled": True, "reason": None},
                {"kind": "adoptPlacement", "subjectId": None, "enabled": True, "reason": None},
            ])
            host.structure_bytes = positioned_protein_bytes(0.0)
            host.replace(value)
            received = [threading.Event(), threading.Event()]
            release = [threading.Event(), threading.Event()]
            requests: list[dict] = []

            def handle_command(handler, command):
                if command["kind"] != "proposePlacement":
                    return False
                index = len(requests)
                requests.append(command)
                received[index].set()
                if not release[index].wait(timeout=15):
                    handler.json_response({"reason": "Controlled response was not released"}, 500)
                    return True
                with host.lock:
                    updated = deepcopy(host.account)
                self.assertEqual(command["expectedRevision"], updated["revision"])
                updated["revision"] += 1
                updated["placement"]["proposalId"] = f"checked-{index + 1}"
                updated["placement"]["transform"]["offsetXAngstrom"] = command["data"]["offsetXAngstrom"]
                updated["placement"]["transform"]["appliedTranslationXAngstrom"] = command["data"]["offsetXAngstrom"]
                host.replace(updated)
                handler.json_response(updated)
                return True

            host.command_handler = handle_command
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                move_x = page.get_by_label("Move X (Å)")
                expect(page.get_by_role("button", name="Use this position")).to_be_enabled()
                move_x.fill("1")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_disabled()
                self.assertTrue(received[0].wait(timeout=10), "First automatic check did not start")
                move_x.fill("2")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_disabled()
                release[0].set()
                self.assertTrue(received[1].wait(timeout=10), "Edited draft did not get its own check")
                expect(page.get_by_label("Current placement status")).to_contain_text("Checking position")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_disabled()
                self.assertEqual([request["data"]["offsetXAngstrom"] for request in requests], [1, 2])
                release[1].set()
                expect(page.get_by_role("button", name="Use this position")).to_be_enabled(timeout=10000)
                move_x.click()
                move_x.press("ControlOrMeta+A")
                move_x.press("-")
                expect(page.get_by_label("Current placement status")).to_contain_text("Incomplete input")
                expect(page.get_by_role("button", name="Use earlier checked position")).to_be_enabled()
                expect(page.get_by_role("button", name="Use this position")).to_have_count(0)
                self.assertEqual(len(requests), 2)
            finally:
                release[0].set()
                release[1].set()
                browser.close()

    def test_placement_guide_and_translated_protein_fit_at_three_viewport_sizes(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        VISUAL_EVIDENCE.mkdir(parents=True, exist_ok=True)
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account()
            value["attempt"] = None
            host.structure_bytes = positioned_protein_bytes()
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 940})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.get_by_role("button", name="Fit view"))\
                    .to_be_enabled(timeout=60000)
                expect(page.get_by_role("button", name="Show whole system")).to_have_count(0)
                for width, height in ((1672, 940), (1024, 800), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                        .to_have_count(1, timeout=60000)
                    expect(page.locator(".placement-scene-legend")).to_contain_text(
                        "Membrane preview — lipids not yet built")
                    self.assert_framed(page, 29.0)
                    page.screenshot(path=str(VISUAL_EVIDENCE / f"controlled-fixture-placement-upper-{width}.png"),
                                    full_page=True)
                page.get_by_role("button", name="Collapse inputs").click()
                self.assert_framed(page, 29.0)
                page.screenshot(path=str(VISUAL_EVIDENCE / "controlled-fixture-placement-upper-820-collapsed.png"),
                                full_page=True)
                self.assertNotIn("selectInspectionSubject", host.commands)
            finally:
                browser.close()

    def test_center_lower_oblique_edge_on_and_zoom_are_visual_only(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        VISUAL_EVIDENCE.mkdir(parents=True, exist_ok=True)
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                for position, z_translation in (("center", 0.0), ("lower", -29.0)):
                    value = positioned_account(position, z_translation)
                    value["attempt"] = None
                    host.structure_bytes = positioned_protein_bytes(z_translation)
                    host.replace(value)
                    page = browser.new_page(viewport={"width": 1024, "height": 800})
                    page.goto(base, wait_until="domcontentloaded")
                    expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                        .to_have_count(1, timeout=60000)
                    self.assert_framed(page, z_translation)
                    page.locator(".placement-scene-surface").screenshot(
                        path=str(VISUAL_EVIDENCE / f"controlled-fixture-placement-{position}-frontal.png"))
                    for view, vector in (("oblique", [1, 1, 1]), ("edge-on", [1, 0, 0])):
                        page.evaluate("""direction => {
                          const mount = document.querySelector('.placement-scene-surface .viewer-mount');
                          const camera = Reflect.get(mount, Symbol.for('molstar.viewer')).plugin.canvas3d.camera;
                          const state = camera.getSnapshot();
                          const target = Array.from(state.target);
                          const distance = Math.hypot(...Array.from(state.position).map((v, i) => v - target[i]));
                          const length = Math.hypot(...direction);
                          state.position = new Float32Array(direction.map((v, i) => target[i] + v / length * distance));
                          state.up = new Float32Array([0, 0, 1]);
                          camera.setState(state, 0);
                        }""", vector)
                        page.wait_for_timeout(250)
                        self.assert_framed(page, z_translation)
                        page.locator(".placement-scene-surface").screenshot(
                            path=str(VISUAL_EVIDENCE / f"controlled-fixture-placement-{position}-{view}.png"))
                    surface = page.locator(".placement-scene-surface")
                    box = surface.bounding_box()
                    assert box is not None
                    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
                    page.mouse.wheel(0, -180)
                    page.wait_for_timeout(250)
                    surface.screenshot(path=str(VISUAL_EVIDENCE / f"controlled-fixture-placement-{position}-zoom.png"))
                    page.get_by_role("button", name="Rotate", exact=True).click()
                    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
                    page.mouse.down()
                    page.mouse.move(box["x"] + box["width"] / 2 + 45,
                                    box["y"] + box["height"] / 2 + 35, steps=5)
                    page.mouse.up()
                    page.wait_for_timeout(250)
                    surface.screenshot(path=str(VISUAL_EVIDENCE / f"controlled-fixture-placement-{position}-pointer-rotated.png"))
                    fit = page.get_by_role("button", name="Fit view")
                    fit.focus()
                    fit.press("Enter")
                    page.wait_for_timeout(250)
                    self.assert_framed(page, z_translation)
                    self.assertEqual(host.structure_bytes, positioned_protein_bytes(z_translation))
                    self.assertNotIn("selectInspectionSubject", host.commands)
                    page.close()
            finally:
                browser.close()

    def test_existing_trpcage_model_remains_legible_inside_centered_guide(self):
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        VISUAL_EVIDENCE.mkdir(parents=True, exist_ok=True)
        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account("center", 0.0, 0.0)
            value["attempt"] = None
            value["study"].update(selectedSourceId="rcsb:1L2Y", selectedSourceLabel="1L2Y first NMR model",
                                  modelIndex=0, chainIds=["A"])
            host.structure_bytes = trpcage_first_model_bytes()
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 800})
                page.goto(base, wait_until="domcontentloaded")
                expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                    .to_have_count(1, timeout=60000)
                for width, height in ((1672, 940), (1024, 800), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                        .to_have_count(1, timeout=60000)
                    page.locator(".placement-scene-surface").screenshot(
                        path=str(VISUAL_EVIDENCE / f"controlled-account-1L2Y-first-model-centered-{width}.png"))
                page.get_by_role("button", name="Collapse inputs").click()
                expect(page.locator(".placement-scene-surface .viewer-mount[data-camera-ready='true']"))\
                    .to_have_count(1, timeout=60000)
                page.locator(".placement-scene-surface").screenshot(
                    path=str(VISUAL_EVIDENCE / "controlled-account-1L2Y-first-model-centered-820-collapsed.png"))
                self.assertEqual(host.structure_bytes, trpcage_first_model_bytes())
                self.assertNotIn("selectInspectionSubject", host.commands)
            finally:
                browser.close()

    def test_saved_multichain_assisted_pose_uses_only_two_planes_and_fit_is_camera_only(self):
        if not SAVED_AB_87.is_file():
            self.skipTest("The saved two-chain PPM pose is unavailable in this workspace")
        self.assertTrue((DIST / "index.html").is_file(), "Build browser/ before this focused test")
        VISUAL_EVIDENCE.mkdir(parents=True, exist_ok=True)
        coordinate_bytes = SAVED_AB_87.read_bytes()
        atom_rows = [line for line in coordinate_bytes.splitlines()
                     if line.startswith((b"ATOM  ", b"HETATM"))]
        self.assertEqual(len(atom_rows), 4393)
        self.assertEqual({chr(line[21]) for line in atom_rows}, {"A", "B"})
        coordinates = [(float(line[30:38]), float(line[38:46]), float(line[46:54]))
                       for line in atom_rows]
        minima = [min(point[axis] for point in coordinates) for axis in range(3)]
        maxima = [max(point[axis] for point in coordinates) for axis in range(3)]
        half_x = max(18, (maxima[0] - minima[0]) / 2 + 8)
        half_y = max(18, (maxima[1] - minima[1]) / 2 + 8)
        centre_x = (minima[0] + maxima[0]) / 2
        centre_y = (minima[1] + maxima[1]) / 2
        guide_points = [[centre_x + sign_x * half_x, centre_y + sign_y * half_y, z]
                        for sign_x in (-1, 1) for sign_y in (-1, 1) for z in (-23, 23)]
        molecule_box = [[x, y, z] for x in (minima[0], maxima[0])
                        for y in (minima[1], maxima[1]) for z in (minima[2], maxima[2])]

        def assert_combined_frame(page):
            projected = page.evaluate("""points => {
              const mount = document.querySelector('.placement-scene-surface .viewer-mount');
              const camera = Reflect.get(mount, Symbol.for('molstar.viewer')).plugin.canvas3d.camera;
              return {width: camera.viewport.width, height: camera.viewport.height,
                pixels: points.map(point => Array.from(camera.project(
                  new Float32Array(4), new Float32Array(point))).slice(0, 2))};
            }""", guide_points + molecule_box)
            self.assertGreater(projected["width"], 0)
            self.assertGreater(projected["height"], 0)
            for x, y in projected["pixels"]:
                self.assertGreater(x, projected["width"] * 0.02)
                self.assertLess(x, projected["width"] * 0.98)
                self.assertGreater(y, projected["height"] * 0.02)
                self.assertLess(y, projected["height"] * 0.98)

        with controlled_account_server() as (host, base), sync_playwright() as playwright:
            value = positioned_account("lower", -29.0, 0.0)
            value["attempt"] = None
            # The coordinates and 87-degree PPM observation are saved facts.
            # This controlled 46 Å frame is a rendering fixture, not the
            # unavailable original checked placement account.
            value["placement"].update(status="notEstablished", transform=None,
                                      topologyKind="one-surface-associated",
                                      tiltDegrees=87.0,
                                      reason="Controlled frame for saved AB PPM pose; support is not claimed.")
            value["study"].update(selectedSourceLabel="Saved two-chain PPM pose",
                                  chainIds=["A", "B"])
            value["inspection"]["annotations"] = [{
                "id": "located", "subjectPartId": "A:1",
                "label": "Located AB region", "meaning": "A real saved protein residue",
                "evidenceId": "evidence-placement-one",
                "geometryFocus": {"authAsymId": "A", "authSeqId": 1,
                                  "insertionCode": None, "authAtomId": None},
            }]
            host.structure_bytes = coordinate_bytes
            host.replace(value)
            browser = playwright.chromium.launch(executable_path=CHROMIUM, headless=True,
                                                 args=["--disable-dev-shm-usage", "--use-gl=angle",
                                                       "--use-angle=swiftshader"])
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 940})
                page.goto(base, wait_until="domcontentloaded")
                mount = page.locator(".placement-scene-surface .viewer-mount")
                expect(mount.locator("canvas")).to_have_count(1, timeout=60000)
                expect(mount).to_have_attribute("data-camera-ready", "true", timeout=60000)
                expect(page.get_by_role("button", name="Fit view")).to_have_attribute(
                    "title", "Fit the protein and membrane preview in the window")
                expect(page.locator(".placement-scene-legend")).to_contain_text(
                    "Membrane preview — lipids not yet built")
                expect(page.locator(".scene-caption")).to_contain_text(
                    "Membrane preview — lipids not yet built")
                expect(page.locator(".placement-proposal-account")).to_contain_text(
                    "Placement frame thickness")
                loaded = page.evaluate("""() => {
                  const mount = document.querySelector('.placement-scene-surface .viewer-mount');
                  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                  const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
                  return {atomCount: structure.elementCount,
                    min: Array.from(structure.boundary.box.min), max: Array.from(structure.boundary.box.max)};
                }""")
                self.assertEqual(loaded["atomCount"], 4393)
                for width, height in ((1672, 940), (1024, 800), (820, 760)):
                    page.set_viewport_size({"width": width, "height": height})
                    expect(mount).to_have_attribute("data-camera-ready", "true", timeout=60000)
                    assert_combined_frame(page)
                    page.locator(".placement-scene-surface").screenshot(
                        path=str(VISUAL_EVIDENCE / f"saved-ab87-two-planes-{width}.png"))
                page.set_viewport_size({"width": 1024, "height": 800})
                for view, direction in (("reverse", [0, 1, 0]), ("oblique", [1, 1, 1]),
                                        ("edge-on", [1, 0, 0])):
                    page.evaluate("""direction => {
                      const mount = document.querySelector('.placement-scene-surface .viewer-mount');
                      const camera = Reflect.get(mount, Symbol.for('molstar.viewer')).plugin.canvas3d.camera;
                      const state = camera.getSnapshot();
                      const target = Array.from(state.target);
                      const distance = Math.hypot(...Array.from(state.position).map((v, i) => v - target[i]));
                      const length = Math.hypot(...direction);
                      state.position = new Float32Array(direction.map((v, i) => target[i] + v / length * distance));
                      state.up = new Float32Array([0, 0, 1]);
                      camera.setState(state, 0);
                    }""", direction)
                    page.wait_for_timeout(200)
                    assert_combined_frame(page)
                    page.locator(".placement-scene-surface").screenshot(
                        path=str(VISUAL_EVIDENCE / f"saved-ab87-two-planes-{view}.png"))
                page.get_by_role("button", name="Located AB region").click()
                expect(page.get_by_role("button", name="Located AB region")).to_have_class(
                    "annotation-item selected")
                self.assertEqual(host.account["inspection"]["focusId"], "A:1")
                commands_before_fit = list(host.commands)
                page.get_by_role("button", name="Fit view").click()
                page.locator(".placement-scene-surface").screenshot(
                    path=str(VISUAL_EVIDENCE / "saved-ab87-fit-with-selection.png"))
                assert_combined_frame(page)
                self.assertEqual(host.commands, commands_before_fit)
                self.assertEqual(host.account["inspection"]["focusId"], "A:1")
                self.assertEqual(host.account["study"]["adoptedPlacementProposalId"], None)
                self.assertEqual(host.structure_bytes, coordinate_bytes)
                after = page.evaluate("""() => {
                  const mount = document.querySelector('.placement-scene-surface .viewer-mount');
                  const viewer = Reflect.get(mount, Symbol.for('molstar.viewer'));
                  const structure = viewer.plugin.managers.structure.hierarchy.current.structures[0].cell.obj.data;
                  return {atomCount: structure.elementCount,
                    min: Array.from(structure.boundary.box.min), max: Array.from(structure.boundary.box.max)};
                }""")
                self.assertEqual(after, loaded)
                page.get_by_role("button", name="Collapse inputs").click()
                expect(mount).to_have_attribute("data-camera-ready", "true", timeout=60000)
                assert_combined_frame(page)
                page.locator(".placement-scene-surface").screenshot(
                    path=str(VISUAL_EVIDENCE / "saved-ab87-two-planes-collapsed.png"))
                self.assertEqual(host.account["inspection"]["focusId"], "A:1")
            finally:
                browser.close()


if __name__ == "__main__":
    unittest.main()
