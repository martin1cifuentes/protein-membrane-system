"""Connected, unrelated-protein Memgen routes through independent ZIP read-back.

Pinned source: RCSB 1CRN PDB, https://files.rcsb.org/download/1CRN.pdb,
wwPDB archive data distributed under CC0. The independently sourced tracked
fixture is the exact 49,491-byte downloaded file, SHA256
42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811.
It has one deposited protein chain, 327 source atoms, no HETATM partners and
three deposited disulfides. The browser explicitly reviews all three bonds.
The upload account truthfully remains an upload of the pinned RCSB file.
The expanded case chooses upper POPC/CHL1 80:20 and lower DLPE/DLPC 70:30
with an upper-side pose; the pure DOPC case uses a centered pose.

Run the long routes only after their isolated 1CRN preparation/placement
preflights and the current Memgen provider matrix have passed. Screenshots
are connected real-route evidence, not controlled presentation fixtures.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.request import urlopen
from zipfile import ZipFile

from playwright.sync_api import expect, sync_playwright

from test_construct_and_minimize import (  # same published actor and stage helpers
    HOST, POLICY, await_completed_stage, await_state, chromium, current_state,
    start_named_route,
)
from test_browser_route import running_host as running_host_without_ppm


ROOT = Path(__file__).resolve().parents[3]
SOURCE = (ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein" /
          "fixtures/1CRN-RCSB.pdb")
SOURCE_SHA256 = "42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811"
SOURCE_ATOM_COUNT = 327
SOURCE_DISULFIDES = {(3, 40), (4, 32), (16, 26)}
TRPCAGE_SOURCE = Path(__file__).parent / "fixtures" / "1L2Y-RCSB.pdb"
TRPCAGE_SOURCE_SHA256 = "5d1bbb545a312dfff1ae1e64b6d8addecb2f561ddc4011aeb5bee9d1dfcd4438"
SOURCE_CASES = {
    "crambin": {"path": SOURCE, "sha256": SOURCE_SHA256,
                "retainedAtomCount": SOURCE_ATOM_COUNT, "disulfides": SOURCE_DISULFIDES,
                "depositedHydrogenCount": 0},
    "trpcage": {"path": TRPCAGE_SOURCE, "sha256": TRPCAGE_SOURCE_SHA256,
                "retainedAtomCount": 154, "disulfides": set(),
                "depositedHydrogenCount": 150},
}
ROUTE_ID = "canonical-protein-general-memgen-construction"
ROUTE_LABEL = "PACKMOL-Memgen: general construction"
ARTIFACTS = ROOT / "out" / "browser-acceptance" / "connected-memgen"
STANDINGS = {"checksPassed", "issuesFound", "checksIncomplete"}
INTENDED = {("upper", "POPC"): 0.8, ("upper", "CHL1"): 0.2,
            ("lower", "DLPE"): 0.7,
            ("lower", "DLPC"): 0.3}
PURE_DOPC_INTENDED = {("upper", "DOPC"): 1.0, ("lower", "DOPC"): 1.0}
INTENDED_BY_ROUTE = {"unequal": INTENDED, "unequal-center": INTENDED,
                     "pure-dopc": PURE_DOPC_INTENDED}


def stage_in(account: dict, stage_id: str) -> dict:
    return next(item for item in account["stages"] if item["stageId"] == stage_id)


def selected_source_atom_lines(source_case: str, *, hydrogen: bool = False) -> list[str]:
    spec = SOURCE_CASES[source_case]
    selected = []
    in_model = source_case == "crambin"
    for line in spec["path"].read_text().splitlines():
        if source_case == "trpcage" and line.startswith("MODEL "):
            if in_model:
                break
            assert line[10:14].strip() == "1"
            in_model = True
            continue
        if source_case == "trpcage" and in_model and line.startswith("ENDMDL"):
            break
        if not in_model or not line.startswith("ATOM  "):
            continue
        is_hydrogen = line[76:78].strip().upper() in ("H", "D")
        if is_hydrogen == hydrogen:
            selected.append(line)
    count = spec["depositedHydrogenCount"] if hydrogen else spec["retainedAtomCount"]
    assert len(selected) == count, (source_case, hydrogen, len(selected), count)
    return selected


def source_atom_ids(source_case: str = "crambin") -> set[str]:
    expected = set()
    for line in selected_source_atom_lines(source_case):
        assert line[21] == "A" and line[16] == " "
        residue = int(line[22:26])
        insertion = line[26].strip()
        atom = line[12:16].strip()
        expected.add(f"0:A:A:{residue}:{insertion}:{atom}")
    assert len(expected) == SOURCE_CASES[source_case]["retainedAtomCount"]
    return expected


def observe_rendered_execution_phases(page):
    # Capture visible current-operation lines for the current running attempt.
    # Diagnostic artifact labels elsewhere in the review are not live progress.
    page.add_init_script("""(() => {
        window.__renderedMemgenPhaseEvents = [];
        window.__renderedFinalMinimizationAttemptIds = [];
        const capture = () => {
            const account = document.querySelector('[aria-label="Current preparation attempt"]');
            const review = account?.nextElementSibling;
            if (!review?.matches('[aria-label="Observed preparation progress"]')) return;
            const running = review.querySelector('.execution-progress-body.is-running');
            const operation = [...(running?.querySelectorAll('p.help-text') || [])]
                .find(line => line.textContent?.trim().startsWith('Current operation:'));
            if (!operation || operation.getClientRects().length === 0) return;
            const attemptId = account.querySelector('details.execution-exact-identity dd.tabular')
                ?.textContent?.trim();
            if (!attemptId || attemptId === 'Not yet accepted') return;
            const text = operation.textContent || '';
            for (const phase of ['provider Packing', 'provider Cleanup',
                                 'amber Parameterization',
                                 'provider Conditioning Restrained',
                                 'provider Conditioning Unrestrained']) {
                if (text.toLowerCase().includes(phase.toLowerCase()) &&
                    !window.__renderedMemgenPhaseEvents.some(event =>
                        event.attemptId === attemptId && event.phase === phase))
                    window.__renderedMemgenPhaseEvents.push({attemptId, phase});
            }
            if (text.toLowerCase().includes('final minimization') &&
                !window.__renderedFinalMinimizationAttemptIds.includes(attemptId))
                window.__renderedFinalMinimizationAttemptIds.push(attemptId);
        };
        const start = () => {
            new MutationObserver(capture).observe(document.documentElement,
                {subtree: true, childList: true, characterData: true});
            capture();
        };
        if (document.readyState === 'loading')
            document.addEventListener('DOMContentLoaded', start, {once: true});
        else start();
    })()""")


def prepare_exact_crambin(page):
    page.locator("#source-upload").set_input_files(str(SOURCE))
    page.locator("#upload-provenance").select_option("experimental")
    page.get_by_role("button", name="Upload source").click()
    expect(page.locator(".source-context")).to_contain_text("Structure displayed", timeout=180000)
    source = current_state(page)
    assert len(source["sourceModels"]) == 1
    model = source["sourceModels"][0]
    assert model["index"] == 0
    assert model["partners"] == []
    if page.locator("#model-index").count():
        page.locator("#model-index").select_option("0")
    else:
        expect(page.locator(".sole-model")).to_contain_text("selected for this draft")
    page.locator("#assembly-choice").select_option("deposited")
    page.get_by_label("Chain A", exact=True).check()
    expect(page.locator(".partner-choice")).to_have_count(0)
    with page.expect_request(lambda request: request.url.endswith("/api/commands") and
                             request.post_data_json.get("kind") == "selectProteinModel") as submitted:
        page.get_by_role("button", name="Assess selected protein").click()
    selection = submitted.value.post_data_json["data"]
    assert selection["modelIndex"] == 0 and selection["biologicalAssemblyId"] is None
    assert selection["chains"] == [{"sourceChain": "A", "copyId": "A"}]
    assert selection["partners"] == []
    plan = await_state(page, lambda state: (state.get("preparationPlan") or {}).get("standing")
                       in ("ready", "partial", "failed"), "1CRN preparation plan", timeout=240)
    assert plan["preparationPlan"]["standing"] == "partial", plan["preparationPlan"]
    bond_decisions = [item for item in plan["preparationReview"]["decisions"]
                      if item["kind"] == "disulfide"]
    assert len(bond_decisions) == len(SOURCE_DISULFIDES)
    observed_pairs = {tuple(sorted((item["residue"]["residue"],
                                    item["partnerResidue"]["residue"])))
                      for item in bond_decisions}
    assert observed_pairs == SOURCE_DISULFIDES
    for remaining in (2, 1, 0):
        page.get_by_role("button", name="Confirm possible bond").click()
        if remaining:
            await_state(page, lambda state: sum(item["kind"] == "disulfide" and
                        item["standing"] == "pending" for item in
                        state["preparationReview"]["decisions"]) == remaining,
                        "next unresolved deposited bond")
            page.get_by_role("button", name="Next unresolved bond decision").click()
    plan = await_state(page, lambda state: (state.get("preparationPlan") or {}).get("standing")
                       in ("ready", "failed"), "joint 1CRN bond plan", timeout=240)
    assert plan["preparationPlan"]["standing"] == "ready", plan["preparationPlan"]
    assert plan["preparationPlan"]["authorizationAvailable"]
    page.get_by_role("button", name="Prepare with recommendations").click()
    prepared = await_state(page, lambda state: (state.get("proteinTask") or {}).get("standing")
                           in ("assessed", "failed", "unavailable"),
                           "exact 1CRN protein preparation", timeout=300)
    assert prepared["proteinTask"]["standing"] == "assessed", prepared["proteinTask"]
    assert prepared["protein"]["status"] == "assessed"
    return prepared


def choose_membrane_and_manual_pose(page, prepared_id: str, route_key: str):
    expected = INTENDED_BY_ROUTE[route_key]
    page.get_by_role("navigation", name="Research work areas").get_by_role(
        "button", name="Membrane").click()
    page.get_by_label("Upper leaflet lipid 1", exact=True).select_option(
        "DOPC" if route_key == "pure-dopc" else "POPC")
    page.get_by_label("Upper leaflet percentage 1", exact=True).fill(
        "100" if route_key == "pure-dopc" else "80")
    page.get_by_label("Lower leaflet lipid 1", exact=True).select_option(
        "DOPC" if route_key == "pure-dopc" else "DLPE")
    page.get_by_label("Lower leaflet percentage 1", exact=True).fill(
        "100" if route_key == "pure-dopc" else "70")
    if route_key in ("unequal", "unequal-center"):
        page.get_by_role("button", name="Add lipid").nth(0).click()
        page.get_by_label("Upper leaflet lipid 2", exact=True).select_option("CHL1")
        page.get_by_label("Upper leaflet percentage 2", exact=True).fill("20")
        page.get_by_role("button", name="Add lipid").nth(1).click()
        page.get_by_label("Lower leaflet lipid 2", exact=True).select_option("DLPC")
        page.get_by_label("Lower leaflet percentage 2", exact=True).fill("30")
    page.get_by_role("button", name="Use this membrane").click()
    assessed = await_state(page, lambda state: (state.get("membrane") or {}).get("status")
                           == "assessed", "assessed unequal membrane")
    assert {(side, fraction["speciesId"]): fraction["fraction"]
            for side in ("upper", "lower") for fraction in assessed["membrane"][side]} == expected
    membrane_id = assessed["membrane"]["modelId"]
    page.get_by_role("navigation", name="Research work areas").get_by_role(
        "button", name="Placement").click()
    # No PPM executable is supplied by this Host. This is a researcher-defined
    # technical pose, with no claim of a native or biological orientation.
    starting_position = "center" if route_key in ("pure-dopc", "unequal-center") else "upper"
    offset_x = 0.0 if starting_position == "center" else 2.0
    page.locator("#starting-position").select_option(starting_position)
    page.get_by_label("Move X (Å)").fill(str(offset_x))
    positioned = await_state(page, lambda state: (state.get("placement") or {}).get("status")
                             in ("supported", "unsupported", "notEstablished") and
                             (state["placement"].get("transform") or {}).get("startingPosition")
                             == starting_position and abs(state["placement"]["transform"].get(
                                 "offsetXAngstrom", 0) - offset_x) < 1e-9,
                             f"checked researcher-defined {starting_position} pose")
    assert positioned["placement"]["status"] == "supported", positioned["placement"]
    assert positioned["placement"]["preparedProteinId"] == prepared_id
    assert positioned["placement"]["membraneModelId"] == membrane_id
    page.get_by_role("button", name="Use this position").click()
    adopted = await_state(page, lambda state: state["study"] is not None and
                          state["study"]["adoptedPlacementProposalId"] ==
                          positioned["placement"]["proposalId"], "adopted exact manual pose")
    return adopted


def capture_review(page, folder: Path, name: str, width: int, height: int,
                   stage_id: str, assessment_id: str):
    page.set_viewport_size({"width": width, "height": height})
    if not page.locator(".workspace.workflow-closed").count():
        page.get_by_role("button", name="Collapse inputs").click()
    expect(page.locator(".workspace.workflow-closed.execution-review")).to_have_count(1)
    expect(page.locator(".viewer-mount canvas")).to_have_count(1, timeout=180000)
    expect(page.locator(".scene-loading")).to_have_count(0, timeout=240000)
    expect(page.locator(".scene-error")).to_have_count(0)
    page.wait_for_function("""() => {
        const mount = document.querySelector('.viewer-mount');
        if (!mount || mount.dataset.cameraReady !== 'true') return false;
        const bounds = mount.getBoundingClientRect();
        return Math.abs(Number(mount.dataset.cameraViewportWidth) - bounds.width) <= 4 &&
            Math.abs(Number(mount.dataset.cameraViewportHeight) - bounds.height) <= 4;
    }""", timeout=30000)
    expect(page.get_by_role("button", name="Export with status")).to_be_enabled()
    expect(page.locator(".stage-strip")).to_contain_text("Minimized")
    scene = page.locator(".scene-panel").bounding_box()
    evidence = page.locator(".evidence-panel").bounding_box()
    stages = page.locator(".stage-strip").bounding_box()
    assert scene and evidence and stages
    folder.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(folder / f"{name}-{width}.png"), full_page=True,
                    animations="disabled")
    (folder / f"{name}-{width}-layout.json").write_text(json.dumps({
        "stageId": stage_id, "assessmentId": assessment_id,
        "viewport": {"width": width, "height": height},
        "scene": scene, "evidence": evidence, "stages": stages,
    }, indent=2) + "\n")
    assert scene["width"] >= 350 and evidence["width"] >= 330
    assert evidence["x"] >= scene["x"] + scene["width"] - 2
    assert abs(evidence["y"] - scene["y"]) <= 3
    assert stages["y"] + stages["height"] <= height + 2
    assert stage_in(current_state(page), stage_id)["assessment"]["id"] == assessment_id


def inspect_bundle(path: Path, stage_id: str, attempt_id: str, revision_id: str,
                   assessment_id: str, check_standing: str, route_key: str,
                   trial_id: str, case_folder: Path,
                   source_case: str = "crambin") -> dict:
    # Executed in the scientific Python environment, after Host/browser closure.
    case_folder = Path(case_folder)
    import gemmi
    import numpy as np
    from openmm import (CMAPTorsionForce, Context, Platform, VerletIntegrator,
                        XmlSerializer, unit)
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PDBxFile, PME
    from openmm.app.internal.amber_file_parser import PrmtopLoader
    from scipy.spatial import cKDTree
    sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.IntegrationTests" /
                           "ExplicitPreparation"))
    from inspect_memgen_chemistry_geometry import (
        check_stereo, component_image_contacts, local_graph, periodic_contacts,
        topology_graph,
    )
    source_spec = SOURCE_CASES[source_case]
    source_sha256 = source_spec["sha256"]
    disulfides = source_spec["disulfides"]
    selected_atom_ids = source_atom_ids(source_case)

    required = {"structure.cif", "topology.json", "system.xml", "state.xml", "manifest.json"}
    with ZipFile(path) as bundle:
        assert required <= set(bundle.namelist())
        assert bundle.testzip() is None
        data = {name: bundle.read(name) for name in required}
    manifest = json.loads(data["manifest.json"])
    assert manifest["stage"]["id"] == stage_id
    assert manifest["stage"]["kind"] == "Minimization"
    assert manifest["stage"]["policyId"] == ROUTE_ID
    assert manifest["attempt"]["id"] == attempt_id
    assert manifest["attempt"]["policyId"] == ROUTE_ID
    assert manifest["study"]["id"] == revision_id
    assert manifest["construction"]["derivation"]["selectedTrialId"] == trial_id
    assert manifest["assessment"]["id"] == assessment_id
    assert manifest["assessment"]["stageId"] == stage_id
    assert manifest["assessment"]["checkStanding"] == check_standing
    assert manifest["lineage"]["sourceCoordinateSha256"] == source_sha256
    assert manifest["lineage"]["sourceModelIndex"] == 0
    source = manifest["preparedProtein"]["source"]
    assert source["kind"] == "upload" and source["sha256"] == source_sha256
    assert source["uploadProvenance"] == "experimental"
    assert manifest["preparedProtein"]["partners"] == []
    assert manifest["preparedProtein"]["recommendationPlan"] is not None
    disulfide_changes = [item for item in manifest["preparedProtein"]["changes"]
                         if item["kind"] == "disulfide"]
    approved_bonds = {tuple(sorted((item["residue"]["residue"],
                                   item["partnerResidue"]["residue"])))
                      for item in disulfide_changes}
    assert approved_bonds == disulfides
    approved_change_ids = {item["subjectId"] for item in
                           manifest["approvedPreparationDecisions"]
                           if item["kind"] == "approvePreparationChange" and
                           item["chosenValue"] == "approved"}
    assert approved_change_ids >= {item["id"] for item in disulfide_changes}
    if source_case == "trpcage":
        assert manifest["preparedProtein"]["modelIndex"] == 0
        assert manifest["preparedProtein"]["sourceModelNumber"] == 1
        deposited_hydrogen_addresses = {
            (0, line[21], int(line[22:26]), line[26].strip(), "A", line[12:16].strip())
            for line in selected_source_atom_lines(source_case, hydrogen=True)
        }
        removed_source_hydrogens = manifest["preparedProtein"]["recommendationPlan"][
            "removedSourceHydrogens"]
        observed_hydrogen_addresses = {
            (item["residue"]["model"], item["residue"]["chain"],
             item["residue"]["residue"], item["residue"]["insertionCode"],
             item["residue"]["copyId"], item["atomName"])
            for item in removed_source_hydrogens
        }
        assert len(removed_source_hydrogens) == len(deposited_hydrogen_addresses) == 150
        assert observed_hydrogen_addresses == deposited_hydrogen_addresses
    provider = manifest["construction"]["provider"]
    assert provider["route"] == "packmolMemgen"
    assert provider["providerName"] == "PACKMOL-Memgen"
    assert provider["providerVersion"] == "2026.3.25"
    assert provider["saltConvention"] == "memgenChargeCompensated"
    catalogue = json.loads(POLICY.read_text())
    selected_policy = next(item for item in catalogue["preparationPolicies"]
                           if item["id"] == ROUTE_ID)
    selected_construction = selected_policy["construction"]
    assert (provider["route"], provider["providerName"], provider["providerVersion"],
            provider["saltConvention"]) == (
                selected_construction["route"], selected_construction["providerName"],
                selected_construction["providerVersion"],
                selected_construction["saltConvention"])

    def pinned_assets(entries: list[dict], keys: tuple[str, ...]) -> dict[str, tuple]:
        assets = {}
        for entry in entries:
            assert entry["id"] not in assets
            path = (POLICY.parent / entry["path"]).resolve()
            assert path.is_file(), path
            assert sha256(path.read_bytes()).hexdigest() == entry["sha256"].lower(), path
            assets[entry["id"]] = tuple(entry[key] for key in keys)
        return assets

    provider_assets = pinned_assets(selected_construction["providerAssets"],
                                    ("version", "sha256"))
    force_field_assets = pinned_assets(selected_policy["forceFieldFiles"],
                                       ("version", "family", "sha256"))
    assert len(provider["assets"]) == len(provider_assets)
    assert {item["id"]: (item["version"], item["sha256"])
            for item in provider["assets"]} == provider_assets
    assert len(manifest["forceFieldAssets"]) == len(force_field_assets)
    assert {item["id"]: (item["version"], item["family"], item["sha256"])
            for item in manifest["forceFieldAssets"]} == force_field_assets
    intended = manifest["membrane"]["intended"]
    requested = {(leaflet["physicalSide"], fraction["speciesId"]): fraction["fraction"]
                 for leaflet in (intended["upper"], intended["lower"])
                 for fraction in leaflet["fractions"]}
    expected = INTENDED_BY_ROUTE[route_key]
    assert requested == expected
    achieved = manifest["construction"]["achievedComposition"]
    observed = {(item["physicalSide"], item["speciesId"]): item["count"] for item in achieved}
    assert set(observed) == set(expected) and all(count > 0 for count in observed.values())
    artifacts = {entry["name"]: entry for entry in manifest["artifacts"]}
    assert set(artifacts) == required - {"manifest.json"}
    for name, entry in artifacts.items():
        assert entry["sha256"] == sha256(data[name]).hexdigest()
        assert entry["byteLength"] == len(data[name])

    topology = json.loads(data["topology.json"])
    atoms, residues = topology["atoms"], topology["residues"]
    correspondence = manifest["correspondence"]["atoms"]
    assert manifest["lineage"]["sourceToResult"] == manifest["correspondence"]
    assert manifest["correspondence"]["resultId"] == manifest["stage"]["moleculeId"] == stage_id
    assert len(atoms) == len(correspondence) == manifest["stage"]["atomCount"]
    mapped_source_ids = [atom["sourceAtomId"] for atom in correspondence
                         if atom["sourceAtomId"]]
    assert len(mapped_source_ids) == len(set(mapped_source_ids)) == len(selected_atom_ids)
    assert set(mapped_source_ids) == selected_atom_ids
    assert all(atom["resultAtomIndex"] == index for index, atom in enumerate(correspondence))
    roles = Counter(atom["moleculeRole"] for atom in correspondence)
    assert all(roles[role] > 0 for role in ("protein", "lipid", "water", "ion"))
    block = gemmi.cif.read_string(data["structure.cif"].decode()).sole_block()
    table = block.find("_atom_site.", ["type_symbol", "label_atom_id", "label_comp_id",
                                      "Cartn_x", "Cartn_y", "Cartn_z"])
    assert len(table) == len(atoms)
    positions = np.empty((len(atoms), 3), dtype=float)
    for index, row in enumerate(table):
        residue = residues[atoms[index]["residueIndex"]]
        assert (row[0], row[1], row[2]) == (
            atoms[index]["element"], atoms[index]["name"], residue["name"])
        positions[index] = [float(row[3]), float(row[4]), float(row[5])]
    system = XmlSerializer.deserialize(data["system.xml"].decode())
    state = XmlSerializer.deserialize(data["state.xml"].decode())
    assert system.getNumParticles() == len(atoms)
    assert any(isinstance(force, CMAPTorsionForce) and force.getNumTorsions() > 0
               for force in system.getForces())
    state_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    assert len(state_positions) == len(atoms)
    coordinate_error = float(np.max(np.linalg.norm(positions - state_positions, axis=1)))
    tolerance = manifest["preparationPolicy"]["exportCoordinateReadBackToleranceAngstrom"]
    assert math.isfinite(coordinate_error) and coordinate_error <= tolerance
    vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    topology_vectors = np.asarray(topology["boxVectorsAngstrom"], dtype=float)
    cell_error = float(np.max(np.abs(vectors - topology_vectors)))
    assert math.isfinite(cell_error) and cell_error <= \
        manifest["preparationPolicy"]["exportCellLengthReadBackToleranceAngstrom"]

    # The actor's checked trial is recorded outside the delivered ZIP. Read its
    # exact, hash-bound Amber artifacts from the retained Host attempt, then
    # compare them with the molecular result independently of the export worker.
    recorded = json.loads((case_folder / "observed-provider-phases.json").read_text())
    assert recorded["attemptId"] == attempt_id
    selected = [item for item in recorded["trials"] if item["trialId"] == trial_id]
    assert len(selected) == 1 and selected[0]["standing"] == "checked"
    manifest_trial = [item for item in manifest["construction"]["derivation"]["trials"]
                      if item["trialId"] == trial_id]
    assert len(manifest_trial) == 1 and manifest_trial[0]["standing"] == "checked"
    recorded_diagnostics = {item["role"]: item for item in selected[0]["diagnosticArtifacts"]}
    manifest_diagnostics = {item["role"]: item for item in
                            manifest_trial[0]["diagnosticArtifacts"]}
    assert len(recorded_diagnostics) == len(selected[0]["diagnosticArtifacts"])
    assert len(manifest_diagnostics) == len(manifest_trial[0]["diagnosticArtifacts"])
    trial_directory = (case_folder / "workspace" / "attempts" / attempt_id /
                       "trials" / trial_id).resolve()
    assert trial_directory.is_dir()

    def diagnostic(role: str) -> Path:
        account = recorded_diagnostics[role]
        assert manifest_diagnostics[role]["sha256"] == account["sha256"]
        assert manifest_diagnostics[role]["fileName"] == account["fileName"]
        candidates = [candidate for candidate in trial_directory.rglob(account["fileName"])
                      if candidate.is_file() and candidate.resolve().is_relative_to(trial_directory)
                      and sha256(candidate.read_bytes()).hexdigest() == account["sha256"]]
        assert len(candidates) == 1, (role, candidates)
        return candidates[0]

    def pdb_rows(file: Path, *, atom_only: bool = False,
                 lines: list[str] | None = None) -> list[dict]:
        return [{"chain": line[21],
                 "residue": int(line[22:26]) if line[22:26].strip().lstrip("-").isdigit() else None,
                 "insertion": line[26].strip(), "name": line[12:16].strip(),
                 "element": line[76:78].strip().upper(),
                 "xyz": tuple(float(line[i:i + 8]) for i in (30, 38, 46))}
                for line in (lines if lines is not None else file.read_text(encoding="ascii").splitlines())
                if line.startswith("ATOM  ") or not atom_only and line.startswith("HETATM")]

    def staged_pdb(expected_sha: str) -> Path:
        candidates = [item for item in (trial_directory / "inputs").glob("*.pdb")
                      if sha256(item.read_bytes()).hexdigest() == expected_sha.lower()]
        assert len(candidates) == 1, expected_sha
        return candidates[0]

    prepared = pdb_rows(staged_pdb(manifest["lineage"]["preparedCoordinateSha256"]))
    oriented = pdb_rows(staged_pdb(manifest["lineage"]["orientedCoordinateSha256"]))
    provider_input = pdb_rows(diagnostic("providerInput"))
    provider_protein = pdb_rows(diagnostic("providerProtein"))
    provider_packed = pdb_rows(diagnostic("providerPacked"))
    parameterized_path = diagnostic("providerParameterized")
    parameterized = pdb_rows(parameterized_path)
    assert len(prepared) == len(oriented) == len(provider_input) == len(provider_protein)
    assert all((left["chain"], left["residue"], left["insertion"], left["name"],
                left["element"]) == (right["chain"], right["residue"], right["insertion"],
                                     right["name"], right["element"])
               for left, right in zip(prepared, oriented))
    assert all(left["element"] == right["element"] and left["xyz"] == right["xyz"]
               for left, right in zip(oriented, provider_input))
    source_by_id = {}
    for atom in pdb_rows(source_spec["path"], atom_only=True,
                         lines=selected_source_atom_lines(source_case)):
        if atom["chain"] != "A":
            continue
        source_id = f"0:A:A:{atom['residue']}:{atom['insertion']}:{atom['name']}"
        assert source_id not in source_by_id
        source_by_id[source_id] = atom
    assert set(source_by_id) == selected_atom_ids
    prepared_index = {}
    for index, atom in enumerate(prepared):
        key = (atom["chain"], atom["residue"], atom["insertion"], atom["name"], atom["element"])
        assert key not in prepared_index
        prepared_index[key] = index

    def map_retained(previous: list[dict], later: list[dict], indices: dict[str, int]) -> dict[str, int]:
        tree = cKDTree(np.asarray([atom["xyz"] for atom in later], dtype=float))
        mapped = {}
        for source_id, before_index in indices.items():
            original = previous[before_index]
            matches = [index for index in tree.query_ball_point(original["xyz"], r=0.025)
                       if later[index]["element"] == original["element"]]
            assert len(matches) == 1, (source_id, len(matches))
            mapped[source_id] = matches[0]
        assert len(set(mapped.values())) == len(mapped)
        return mapped

    source_indices = {}
    for source_id, atom in source_by_id.items():
        key = (atom["chain"], atom["residue"], atom["insertion"], atom["name"], atom["element"])
        assert key in prepared_index, source_id
        source_indices[source_id] = prepared_index[key]
        if source_case == "trpcage":
            # All selected heavy atoms are retained in this complete NMR model.
            # Models 1 and 2 have the same atom addresses, so identity alone
            # cannot prove which deposited coordinates were prepared.
            assert prepared[source_indices[source_id]]["xyz"] == atom["xyz"], source_id
    source_indices = map_retained(provider_input, provider_protein, source_indices)
    source_indices = map_retained(provider_protein, provider_packed, source_indices)
    source_indices = map_retained(provider_packed, parameterized, source_indices)

    prmtop = diagnostic("amberTopology")
    restart = diagnostic("amberFinalRestart")
    diagnostic("providerRestrainedLog")
    diagnostic("providerRestrainedRestart")
    diagnostic("providerUnrestrainedLog")
    initial_restart = AmberInpcrdFile(str(diagnostic("amberCoordinates")))
    raw_amber = PrmtopLoader(str(prmtop))
    conditioned = AmberInpcrdFile(str(restart))
    amber = AmberPrmtopFile(str(prmtop), periodicBoxVectors=conditioned.boxVectors)
    amber_atoms = list(amber.topology.atoms())
    assert len(amber_atoms) == len(atoms) == len(conditioned.positions)
    assert len(parameterized) == len(amber_atoms)
    for index, (pdb_atom, amber_atom) in enumerate(zip(parameterized, amber_atoms)):
        first_amino_hydrogen_alias = (index == 1 and pdb_atom["name"] == "H1" and
                                      amber_atom.name == "H" and
                                      [item["name"] for item in parameterized[:4]] ==
                                      ["N", "H1", "H2", "H3"])
        assert (pdb_atom["name"] == amber_atom.name or first_amino_hydrogen_alias)
        assert pdb_atom["element"] == amber_atom.element.symbol.upper()
    for source_id, index in source_indices.items():
        source_atom = source_by_id[source_id]
        mapped = correspondence[index]
        assert mapped["resultAtomIndex"] == index and mapped["sourceAtomId"] == source_id
        assert mapped["role"] == "source" and mapped["moleculeRole"] == "protein"
        assert amber_atoms[index].name == source_atom["name"]
        assert amber_atoms[index].element.symbol.upper() == source_atom["element"]
        assert mapped["sourceResidue"] == {"model": 0, "chain": "A",
                                           "residue": source_atom["residue"],
                                           "insertionCode": source_atom["insertion"],
                                           "copyId": "A"}
    assert all((atom.name, atom.element.symbol, atom.residue.name) ==
               (atoms[index]["name"], atoms[index]["element"],
                residues[atoms[index]["residueIndex"]]["name"])
               for index, atom in enumerate(amber_atoms))
    amber_bonds = {tuple(sorted((left.index, right.index)))
                   for left, right in amber.topology.bonds()}
    exported_bonds = {tuple(sorted(item["atomIndices"])) for item in topology["bonds"]}
    assert amber_bonds == exported_bonds
    for first, second in disulfides:
        first_sg = source_indices[f"0:A:A:{first}::SG"]
        second_sg = source_indices[f"0:A:A:{second}::SG"]
        assert tuple(sorted((first_sg, second_sg))) in amber_bonds
    selected_species = {species for _, species in expected}
    membrane_species = {item["speciesId"]: item for item in manifest["membrane"]["species"]}
    assert set(membrane_species) == selected_species
    references = {}
    reference_entries = {}
    for species in selected_species:
        entry = next(item for item in catalogue["lipids"] if item["speciesId"] == species)
        template = (POLICY.parent / entry["coordinateTemplatePath"]).resolve()
        assert template.is_file() and sha256(template.read_bytes()).hexdigest() == \
            entry["coordinateTemplateSha256"] == membrane_species[species]["coordinateTemplateSha256"]
        assert entry["templateSha256"] == membrane_species[species]["templateSha256"]
        reference_atoms, reference_edges = topology_graph(PDBxFile(str(template)).topology)
        references[species] = local_graph(reference_atoms, reference_edges,
                                          list(range(len(reference_atoms))))
        reference_entries[species] = entry
    initial_positions = np.asarray(initial_restart.positions.value_in_unit(unit.angstrom))
    pdb_positions = np.asarray([item["xyz"] for item in parameterized], dtype=float)
    initial_amber_box = np.asarray([
        [float(value.value_in_unit(unit.angstrom)) for value in vector]
        for vector in initial_restart.boxVectors])
    conditioned_amber_box = np.asarray([
        [float(value.value_in_unit(unit.angstrom)) for value in vector]
        for vector in conditioned.boxVectors])
    initial_cell_lengths = np.diag(initial_amber_box)
    conditioned_cell_lengths = np.diag(conditioned_amber_box)
    assert np.all(initial_cell_lengths > 0) and len(initial_positions) == len(pdb_positions)
    assert np.array_equal(conditioned_amber_box, np.diag(conditioned_cell_lengths))
    assert np.array_equal(initial_amber_box, conditioned_amber_box)
    amber_shift = initial_positions[0] - pdb_positions[0]
    indexed_shift_error = (initial_positions - pdb_positions - amber_shift +
                           initial_cell_lengths / 2) % initial_cell_lengths - initial_cell_lengths / 2
    assert float(np.max(np.abs(indexed_shift_error))) <= 0.025
    output_midplane = manifest["placement"]["proposal"]["midplaneAngstrom"] + amber_shift[2]
    conditioned_positions = np.asarray(conditioned.positions.value_in_unit(unit.angstrom))
    final_cell_lengths = np.asarray([vectors[i][i] for i in range(3)], dtype=float)
    assert np.all(final_cell_lengths > 0)
    assert np.array_equal(conditioned_amber_box, vectors)
    _, amber_adjacency = topology_graph(amber.topology)
    parent = list(range(len(atoms)))

    def component(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left, right in amber_bonds:
        parent[component(right)] = component(left)
    molecular_groups = defaultdict(list)
    for index in range(len(atoms)):
        molecular_groups[component(index)].append(index)
    conditioned_lipid_counts = Counter()
    final_lipid_counts = Counter()
    conditioned_stereo_checks = Counter()
    final_stereo_checks = Counter()
    for indices in molecular_groups.values():
        members = [correspondence[index] for index in indices]
        if any(member["moleculeRole"] == "lipid" for member in members):
            actual = local_graph(amber_atoms, amber_adjacency, indices)
            candidates = [species for species, reference in references.items()
                          if len(actual[0]) == len(reference[0]) and
                          Counter(actual[0]) == Counter(reference[0])]
            assert len(candidates) == 1, candidates
            species = candidates[0]
            checked = check_stereo(species, reference_entries[species]["stereoChecks"],
                                   references[species], actual,
                                   conditioned_positions[indices], conditioned_cell_lengths)
            assert checked == len(reference_entries[species]["stereoChecks"])
            conditioned_stereo_checks[species] += checked
            final_checked = check_stereo(species, reference_entries[species]["stereoChecks"],
                                         references[species], actual,
                                         state_positions[indices], final_cell_lengths)
            assert final_checked == checked
            final_stereo_checks[species] += final_checked
            head_name = "O1" if species == "CHL1" else "P31"
            heads = [index for index in indices if amber_atoms[index].name == head_name]
            assert len(heads) == 1, (species, heads)
            conditioned_z = (conditioned_positions[heads[0], 2] - output_midplane +
                             conditioned_cell_lengths[2] / 2) % conditioned_cell_lengths[2] - \
                conditioned_cell_lengths[2] / 2
            final_z = (state_positions[heads[0], 2] - output_midplane +
                       final_cell_lengths[2] / 2) % final_cell_lengths[2] - final_cell_lengths[2] / 2
            assert abs(conditioned_z) > 1e-6 and abs(final_z) > 1e-6
            conditioned_side = "upper" if conditioned_z > 0 else "lower"
            final_side = "upper" if final_z > 0 else "lower"
            assert conditioned_side == final_side
            identities = {(member["moleculeRole"], member["physicalSide"],
                           member["generatedSpeciesId"]) for member in members}
            assert identities == {("lipid", final_side, species)}
            conditioned_lipid_counts[(conditioned_side, species)] += 1
            final_lipid_counts[(final_side, species)] += 1
    assert conditioned_lipid_counts == final_lipid_counts
    assert final_lipid_counts == Counter(observed)
    assert final_lipid_counts == Counter({(item["physicalSide"], item["speciesId"]): item["count"]
                                    for item in selected[0]["achievedLipidCounts"]})
    assert final_lipid_counts == Counter({(item["physicalSide"], item["speciesId"]): item["count"]
                                    for item in manifest["construction"]["achievedComposition"]})
    groups = list(molecular_groups.values())
    molecule_ids = [0] * len(atoms)
    for number, indices in enumerate(groups):
        for index in indices:
            molecule_ids[index] = number
    molecule_roles = [item["moleculeRole"] for item in correspondence]
    elements = [atom.element.symbol for atom in amber_atoms]
    conditioned_contacts = periodic_contacts(conditioned_positions, conditioned_cell_lengths,
                                             elements, molecule_roles, molecule_ids)
    conditioned_image_contacts = component_image_contacts(
        conditioned_positions, conditioned_cell_lengths, amber_adjacency, elements,
        molecule_roles, groups)
    assert not conditioned_contacts["violations"], conditioned_contacts
    assert not conditioned_image_contacts, conditioned_image_contacts
    final_contacts = periodic_contacts(state_positions, final_cell_lengths,
                                       elements, molecule_roles, molecule_ids)
    final_image_contacts = component_image_contacts(
        state_positions, final_cell_lengths, amber_adjacency, elements,
        molecule_roles, groups)
    final_local_state = manifest["observation"]["localState"]
    assert final_local_state["standing"] == "Observed"
    final_measurements = {(item["name"], item["scope"]): item["value"]
                          for item in manifest["observation"]["measurements"]}
    reported_heavy = final_measurements[("minimumIntermolecularHeavyAtomDistanceAngstrom",
                                         "wholeSystem")]
    assert math.isfinite(reported_heavy) and reported_heavy > 0
    heavy_indices = np.asarray([index for index, element in enumerate(elements)
                                if element.upper() not in {"H", "D"}], dtype=int)
    wrapped_heavy = state_positions[heavy_indices] % final_cell_lengths
    nearby = cKDTree(wrapped_heavy, boxsize=final_cell_lengths).query_pairs(
        reported_heavy + 1e-4, output_type="ndarray")
    independent_heavy = min((float(np.linalg.norm(
        (wrapped_heavy[first] - wrapped_heavy[second] + final_cell_lengths / 2) %
        final_cell_lengths - final_cell_lengths / 2))
        for first, second in nearby
        if molecule_ids[heavy_indices[first]] != molecule_ids[heavy_indices[second]]),
        default=None)
    assert independent_heavy is not None and abs(independent_heavy - reported_heavy) <= 1e-5
    protein_indices = sorted(index for index, item in enumerate(correspondence)
                             if item["moleculeRole"] == "protein")
    protein_set = set(protein_indices)
    protein_bonds = [(first, second) for first, second in amber_bonds
                     if first in protein_set and second in protein_set]
    protein_neighbors = {index: set() for index in protein_indices}
    for first, second in protein_bonds:
        protein_neighbors[first].add(second)
        protein_neighbors[second].add(first)

    def final_distance(first: int, second: int) -> float:
        delta = (state_positions[first] - state_positions[second] + final_cell_lengths / 2) % \
            final_cell_lengths - final_cell_lengths / 2
        return float(np.linalg.norm(delta))

    geometry_spec = selected_policy["stageProteinGeometryMeasurement"]
    final_geometry_distances = {
        "covalentBond": [final_distance(first, second) for first, second in protein_bonds],
        "chainContinuity": [], "nonbondedDistance": [],
    }
    protein_residue_atoms = defaultdict(dict)
    for index in protein_indices:
        protein_residue_atoms[atoms[index]["residueIndex"]][atoms[index]["name"]] = index
    protein_chains = defaultdict(list)
    for residue_index in sorted(protein_residue_atoms):
        protein_chains[residues[residue_index]["chainIndex"]].append(residue_index)
    for residue_indices in protein_chains.values():
        for earlier, later in zip(residue_indices, residue_indices[1:]):
            assert "C" in protein_residue_atoms[earlier]
            assert "N" in protein_residue_atoms[later]
            final_geometry_distances["chainContinuity"].append(final_distance(
                protein_residue_atoms[earlier]["C"], protein_residue_atoms[later]["N"]))
    excluded = {}
    for index in protein_indices:
        reached, frontier = {index}, {index}
        for _ in range(geometry_spec["excludedBondHops"]):
            frontier = set().union(*(protein_neighbors[item] for item in frontier)) - reached
            reached.update(frontier)
        excluded[index] = reached
    wrapped_protein = state_positions[protein_indices] % final_cell_lengths
    near_protein = cKDTree(wrapped_protein, boxsize=final_cell_lengths).query_pairs(
        geometry_spec["neighborSearchRadiusAngstrom"], output_type="ndarray")
    for first_local, second_local in near_protein:
        first, second = protein_indices[first_local], protein_indices[second_local]
        if second not in excluded[first]:
            distance = final_distance(first, second)
            if distance <= geometry_spec["neighborSearchRadiusAngstrom"]:
                final_geometry_distances["nonbondedDistance"].append(distance)
    geometry = manifest["observation"]["proteinGeometry"]
    assert geometry["standing"] == "Observed"
    reported_geometry = {item["kind"]: item for item in geometry["kinds"]}
    assert set(reported_geometry) == set(geometry_spec["requiredKinds"])
    for kind, distances in final_geometry_distances.items():
        assert distances, kind
        reported = reported_geometry[kind]
        assert reported["standing"] == "Observed"
        assert reported["measuredCount"] == reported["eligibleCount"] == len(distances)
        assert abs(reported["minimumDistanceAngstrom"] - min(distances)) <= 1e-5
        assert abs(reported["maximumDistanceAngstrom"] - max(distances)) <= 1e-5
    if check_standing == "checksPassed":
        for item in selected_policy["stageProteinGeometryCriteria"]:
            if item["stageKind"] != "Minimization":
                continue
            criterion = item["criterion"]
            distances = final_geometry_distances[criterion["kind"]]
            minimum = criterion["minimumObservedAngstrom"]
            maximum = criterion["maximumObservedAngstrom"]
            assert minimum is None or min(distances) >= minimum - 1e-6
            assert maximum is None or max(distances) <= maximum + 1e-6
    residue_counts = Counter(residue.name for residue in amber.topology.residues())
    conditions = manifest["construction"]["derivation"]["conditions"]
    assert conditions == selected[0]["conditions"]
    assert all(conditions[f"retained{species}Count"] == 0
               for species in ("Water", "Sodium", "Chloride"))
    assert (residue_counts["HOH"], residue_counts["Na+"], residue_counts["Cl-"]) == (
        conditions["finalWaterCount"], conditions["finalSodiumCount"],
        conditions["finalChlorideCount"])
    fresh = amber.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                               constraints=HBonds, rigidWater=True,
                               ewaldErrorTolerance=0.0005, removeCMMotion=True,
                               hydrogenMass=None)
    assert XmlSerializer.serialize(fresh) == XmlSerializer.serialize(system)
    amber_cmap_terms = int(raw_amber._raw_data.get("CMAP_COUNT", ["0"])[0])
    imported_cmap_terms = sum(force.getNumTorsions() for force in system.getForces()
                              if isinstance(force, CMAPTorsionForce))
    assert amber_cmap_terms == imported_cmap_terms and amber_cmap_terms > 0

    attempt_directory = trial_directory.parent.parent
    exact_stage_files = {
        "topology.json": trial_directory / f"constructed-{trial_id}-topology.json",
        "system.xml": trial_directory / f"constructed-{trial_id}-system.xml",
        "state.xml": attempt_directory / "minimized-state.xml",
        "structure.cif": attempt_directory / "minimized-coordinates.cif",
    }
    for exported_name, saved in exact_stage_files.items():
        assert saved.is_file() and sha256(saved.read_bytes()).digest() == \
            sha256(data[exported_name]).digest(), saved
    conditioned_state = XmlSerializer.deserialize(
        (trial_directory / f"constructed-{trial_id}-state.xml").read_text())
    conditioned_positions = conditioned.positions.value_in_unit(unit.angstrom)
    imported_positions = conditioned_state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    assert float(np.max(np.linalg.norm(conditioned_positions - imported_positions, axis=1))) < 1e-7
    imported_box = conditioned_state.getPeriodicBoxVectors(
        asNumpy=True).value_in_unit(unit.angstrom)
    assert np.array_equal(conditioned_amber_box, imported_box)

    # Re-evaluate the delivered final System at its delivered final positions.
    # Project fresh forces onto the HBonds constraint tangent rather than
    # accepting a serialized State force or an owner-reported convergence label.
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "1"})
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    evaluated = context.getState(getForces=True, getEnergy=True)
    force = np.asarray(evaluated.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    saved_force = np.asarray(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    assert float(np.max(np.abs(force - saved_force))) <= 0.05
    assert abs((evaluated.getPotentialEnergy() - state.getPotentialEnergy()).value_in_unit(
        unit.kilojoule_per_mole)) <= 0.1
    del context, integrator
    positions_nm = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    parent = list(range(len(atoms)))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    rows = []
    maximum_constraint_error = 0.0
    for index in range(system.getNumConstraints()):
        left, right, distance = system.getConstraintParameters(index)
        delta = positions_nm[left] - positions_nm[right]
        length = float(np.linalg.norm(delta))
        target = distance.value_in_unit(unit.nanometer)
        assert length > 0 and target > 0
        maximum_constraint_error = max(maximum_constraint_error, abs(length - target) / target)
        rows.append((left, right, delta / length))
        parent[root(right)] = root(left)
    groups = defaultdict(list)
    for row in rows:
        groups[root(row[0])].append(row)
    tangent = force.copy()
    for members in groups.values():
        indices = sorted({index for left, right, _ in members for index in (left, right)})
        local = {index: position for position, index in enumerate(indices)}
        jacobian = np.zeros((len(members), 3 * len(indices)))
        for row_index, (left, right, direction) in enumerate(members):
            jacobian[row_index, 3 * local[left]:3 * local[left] + 3] = direction
            jacobian[row_index, 3 * local[right]:3 * local[right] + 3] = -direction
        flat = force[indices].reshape(-1)
        multipliers = np.linalg.solve(jacobian @ jacobian.T, jacobian @ flat)
        tangent[indices] = (flat - jacobian.T @ multipliers).reshape((len(indices), 3))
    tangent_rms = math.sqrt(float(np.sum(tangent * tangent)) / len(atoms))
    measured = {(item["name"], item["scope"]): item["value"]
                for item in manifest["observation"]["measurements"]}
    reported_tangent = measured[("finalRmsForce", "final unrestrained constraint tangent; per particle")]
    reported_constraint = measured[("maximumRelativeConstraintError", "final HBonds constraints")]
    assert math.isfinite(tangent_rms) and tangent_rms <= 10.0
    assert abs(tangent_rms - reported_tangent) < 1e-3
    assert maximum_constraint_error <= 1e-5
    assert abs(maximum_constraint_error - reported_constraint) < 1e-7
    return {"stageId": stage_id, "attemptId": attempt_id, "studyRevisionId": revision_id,
            "assessmentId": assessment_id, "checkStanding": check_standing,
            "sourceSha256": source_sha256, "sourceCase": source_case,
            "atomCount": len(atoms),
            "selectedTrialId": trial_id, "amberTopologySha256": sha256(prmtop.read_bytes()).hexdigest(),
            "amberFinalRestartSha256": sha256(restart.read_bytes()).hexdigest(),
            "amberCmapTerms": amber_cmap_terms,
            "amberFrameMidplaneZAngstrom": float(output_midplane),
            "maximumIndexedAmberFrameDeviationAngstrom": float(np.max(np.abs(indexed_shift_error))),
            "sourceAtomsIndividuallyMapped": len(source_indices),
            "conditionedStereoChecksBySpecies": dict(conditioned_stereo_checks),
            "finalStereoChecksBySpecies": dict(final_stereo_checks),
            "conditionedPeriodicHeavyContactViolations": len(conditioned_contacts["violations"]),
            "conditionedComponentImageHeavyContactViolations": len(conditioned_image_contacts),
            "finalPeriodicHeavyContactViolations": len(final_contacts["violations"]),
            "finalComponentImageHeavyContactViolations": len(final_image_contacts),
            "finalMinimumIntermolecularHeavyDistanceAngstrom": independent_heavy,
            "finalProteinGeometryCounts": {kind: len(distances)
                                           for kind, distances in final_geometry_distances.items()},
            "finalStateLipidComponents": [
                {"physicalSide": side, "speciesId": species, "count": count}
                for (side, species), count in sorted(final_lipid_counts.items())],
            "pinnedProviderAssetCount": len(provider_assets),
            "pinnedForceFieldAssetCount": len(force_field_assets),
            "amberWaterIonCounts": {"water": residue_counts["HOH"],
                                    "sodium": residue_counts["Na+"],
                                    "chloride": residue_counts["Cl-"]},
            "constraintTangentRmsKjMolNm": tangent_rms,
            "maximumRelativeConstraintError": maximum_constraint_error,
            "achievedComposition": [
                {"physicalSide": side, "speciesId": species, "count": count}
                for (side, species), count in sorted(observed.items())],
            "moleculeRoles": dict(roles),
            "coordinateMaxDeviationAngstrom": coordinate_error,
            "cellVectorMaxDeviationAngstrom": cell_error,
            "bundleSha256": sha256(path.read_bytes()).hexdigest()}


class ConnectedMemgenCrambinRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert SOURCE.stat().st_size == 49491
        assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
        assert sum(line.startswith("ATOM  ") for line in SOURCE.read_text().splitlines()) == SOURCE_ATOM_COUNT
        assert not any(line.startswith("HETATM") for line in SOURCE.read_text().splitlines())
        assert {tuple(sorted((int(line[17:21]), int(line[31:35]))))
                for line in SOURCE.read_text().splitlines() if line.startswith("SSBOND")} == \
            SOURCE_DISULFIDES
        for path in (HOST, POLICY, ROOT / "out/host/wwwroot/index.html"):
            assert path.is_file(), path
        ARTIFACTS.mkdir(parents=True, exist_ok=True)

    def test_00_preflight_source_preparation_and_manual_position(self):
        """Bound the real source/pose prerequisites without starting Memgen."""
        self._preflight("unequal")

    def test_01_preflight_pure_dopc_and_manual_position(self):
        """Bound the separate pure-DOPC route before its costly provider run."""
        self._preflight("pure-dopc")

    def _preflight(self, route_key: str):
        folder = Path(tempfile.mkdtemp(prefix=f"preflight-1crn-{route_key}-", dir=ARTIFACTS))
        with patch.dict(os.environ, {"PIM_BROWSER_HOST": str(HOST),
                                     "PIM_AMBERTOOLS_HOME": str(ROOT / "out/ambertools26")}), \
                running_host_without_ppm(folder / "workspace", POLICY, ppm=None) as base, \
                sync_playwright() as playwright:
            browser = chromium(playwright)
            try:
                page = browser.new_page(viewport={"width": 1024, "height": 768})
                page.goto(base, wait_until="domcontentloaded")
                prepared = prepare_exact_crambin(page)
                adopted = choose_membrane_and_manual_pose(
                    page, prepared["protein"]["subjectId"], route_key)
                assert adopted["study"]["adoptedPlacementProposalId"] == \
                    adopted["placement"]["proposalId"]
                assert adopted["attempt"] is None and adopted["stages"] == []
                ppm = next(item for item in adopted["placementMethods"]
                           if item["method"] == "PPM")
                assert ppm["standing"] == "unavailable", ppm
                route = next(item for item in adopted["constructionRoutes"]
                             if item["policyId"] == ROUTE_ID)
                assert route["available"], route
                page.screenshot(path=str(folder / "preflight-adopted-manual-1024.png"),
                                full_page=True, animations="disabled")
                page.get_by_role("navigation", name="Research work areas").get_by_role(
                    "button", name="Preparation", exact=True).click()
                recipe = page.locator(".construction-route").filter(has_text=ROUTE_LABEL)
                expect(recipe).to_have_count(1)
                recipe_text = recipe.locator("div").first.bounding_box()
                recipe_action = recipe.get_by_role("button", name="Build and minimize")
                expect(recipe_action).to_be_enabled()
                action_box = recipe_action.bounding_box()
                assert recipe_text and action_box
                assert recipe_text["y"] + recipe_text["height"] <= action_box["y"] - 1, \
                    (recipe_text, action_box)
                page.screenshot(path=str(folder / "preflight-preparation-recipes-1024.png"),
                                full_page=True, animations="disabled")
            finally:
                try:
                    with urlopen(base + "/api/state", timeout=10) as response:
                        (folder / "preflight-account.json").write_text(json.dumps(
                            json.load(response), indent=2) + "\n")
                except Exception as failure:
                    (folder / "account-read-error.txt").write_text(str(failure) + "\n")
                browser.close()

    def test_one_unequal_leaflet_build_to_export_readback(self):
        """Real four-component unequal/sterol case with 1CRN and export read-back."""
        self._run_connected_route("unequal")

    def test_two_pure_dopc_build_to_export_readback(self):
        """Separate real pure-DOPC route with unrelated 1CRN and export read-back."""
        self._run_connected_route("pure-dopc")

    def _run_connected_route(self, route_key: str, source_case: str = "crambin",
                             prepare_source=None):
        expected = INTENDED_BY_ROUTE[route_key]
        source_spec = SOURCE_CASES[source_case]
        assert sha256(source_spec["path"].read_bytes()).hexdigest() == source_spec["sha256"]
        if prepare_source is None:
            assert source_case == "crambin"
            prepare_source = prepare_exact_crambin
        folder = Path(tempfile.mkdtemp(
            prefix=f"real-memgen-{source_case}-{route_key}-", dir=ARTIFACTS))
        workspace = folder / "workspace"
        download_path = folder / "downloaded-stage.zip"
        print(f"Connected Memgen route workspace: {folder}", flush=True)
        # Invoke the published Host directly so start-local.sh cannot discover
        # and re-enable PPM when this route is deliberately manual.
        with patch.dict(os.environ, {"PIM_BROWSER_HOST": str(HOST),
                                     "PIM_AMBERTOOLS_HOME": str(ROOT / "out/ambertools26")}), \
                running_host_without_ppm(workspace, POLICY, ppm=None) as base, \
                sync_playwright() as playwright:
            browser = chromium(playwright)
            try:
                page = browser.new_page(viewport={"width": 1672, "height": 941},
                                        accept_downloads=True)
                observe_rendered_execution_phases(page)
                page.goto(base, wait_until="domcontentloaded")
                prepared = prepare_source(page)
                adopted = choose_membrane_and_manual_pose(
                    page, prepared["protein"]["subjectId"], route_key)
                revision_id = adopted["study"]["id"]
                admitted = start_named_route(page, ROUTE_LABEL, ROUTE_ID)
                attempt_id = admitted["attempt"]["attemptId"]
                assert admitted["attempt"]["studyRevisionId"] == revision_id
                assert admitted["stages"] == []
                before_duplicate = current_state(page)["attempt"]
                assert before_duplicate["status"] in ("pending", "running"), (
                    f"Memgen stopped before duplicate-command check: "
                    f"{before_duplicate['status']} / {before_duplicate.get('failureCode')} / "
                    f"{before_duplicate.get('message')}")
                duplicate = page.request.post(base + "/api/commands", data=json.dumps({
                    "kind": "buildAndMinimize", "data": {"policyId": ROUTE_ID},
                    "expectedRevision": current_state(page)["revision"],
                }), headers={"Content-Type": "application/json", "Origin": base})
                assert duplicate.status == 422, duplicate.status
                assert current_state(page)["attempt"]["attemptId"] == attempt_id
                seen_phases = set()
                captured_phases = set()
                rendered_phase_events = []
                rendered_final_minimization = False
                reattachment = None
                deadline = time.monotonic() + 3300
                while time.monotonic() < deadline:
                    if page is None:
                        # Keep the verified reattachment, then leave Mol* closed
                        # during the long provider packing search. The browser's
                        # software renderer otherwise competes with Packmol.
                        with urlopen(base + "/api/state", timeout=10) as response:
                            current = json.load(response)
                        pending = current.get("attempt") or {}
                        assert pending.get("attemptId") == attempt_id, pending.get("attemptId")
                        assert pending.get("studyRevisionId") == revision_id, pending.get(
                            "studyRevisionId")
                        assert current["study"]["id"] == revision_id
                        if pending.get("status") == "running" and \
                                pending.get("phase") == "providerPacking" and \
                                pending.get("constructed") is None:
                            assert current["stages"] == []
                            time.sleep(1.0)
                            continue
                        page = browser.new_page(viewport={"width": 1672, "height": 941},
                                                accept_downloads=True)
                        observe_rendered_execution_phases(page)
                        page.goto(base, wait_until="domcontentloaded")
                        await_state(page, lambda state: (state.get("attempt") or {}).get(
                            "attemptId") == attempt_id,
                            "same Memgen attempt after packing", timeout=120)
                        page.get_by_role("navigation", name="Research work areas").get_by_role(
                            "button", name="Preparation", exact=True).click()
                        page.get_by_role("button", name="Review attempt", exact=True).click()
                    current = current_state(page)
                    attempt = current.get("attempt") or {}
                    phase = attempt.get("phase")
                    if phase:
                        seen_phases.add(phase)
                    if attempt.get("status") in ("pending", "running") and phase in (
                            "providerPopulation", "providerPacking", "providerCleanup",
                            "amberParameterization", "providerConditioningRestrained",
                            "providerConditioningUnrestrained"):
                        assert attempt.get("constructed") is None, attempt
                        assert current["stages"] == [], (phase, current["stages"])
                    if phase in ("providerPacking", "providerConditioningUnrestrained") and \
                            phase not in captured_phases:
                        width, height = (1024, 768) if phase == "providerPacking" else (820, 760)
                        page.set_viewport_size({"width": width, "height": height})
                        expect(page.locator(".stage-strip")).to_be_visible()
                        if phase == "providerPacking":
                            expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                                re.compile("provider packing", re.IGNORECASE), timeout=30000)
                            assert page.locator(".stage-card").count() == 0
                            expect(page.get_by_label("System stages and attempts")).to_contain_text(
                                "No minimized or equilibrated system stage has completed.")
                        else:
                            expect(page.locator(
                                '[aria-label="Observed preparation progress"] '
                                '.execution-progress-body.is-running p.help-text')).to_contain_text(
                                    "Current operation: provider Conditioning Unrestrained.",
                                    timeout=30000)
                        page.screenshot(path=str(folder / f"{phase}-{width}.png"),
                                        full_page=True, animations="disabled")
                        captured_phases.add(phase)
                        if phase == "providerPacking" and reattachment is None:
                            rendered_phase_events.extend(page.evaluate(
                                "window.__renderedMemgenPhaseEvents || []"))
                            rendered_final_minimization |= attempt_id in page.evaluate(
                                "window.__renderedFinalMinimizationAttemptIds || []")
                            page.close()
                            page = browser.new_page(viewport={"width": 1672, "height": 941},
                                                    accept_downloads=True)
                            observe_rendered_execution_phases(page)
                            page.goto(base, wait_until="domcontentloaded")
                            reopened = await_state(page, lambda state: state.get("attempt") is not None and
                                                   state["attempt"]["attemptId"] == attempt_id,
                                                   "same Memgen attempt after browser reattachment", timeout=120)
                            assert reopened["study"]["id"] == revision_id
                            assert len(list((workspace / "attempts").iterdir())) == 1
                            page.get_by_role("navigation", name="Research work areas").get_by_role(
                                "button", name="Preparation", exact=True).click()
                            page.get_by_role("button", name="Review attempt", exact=True).click()
                            attempt_review = page.get_by_label("Current preparation attempt")
                            expect(attempt_review).to_be_visible()
                            attempt_review.locator("details.execution-exact-identity summary").click()
                            expect(attempt_review).to_contain_text(attempt_id)
                            if reopened["attempt"]["status"] == "running" and \
                                    reopened["attempt"].get("constructed") is None:
                                assert reopened["stages"] == []
                                assert page.locator(".stage-card").count() == 0
                                if reopened["attempt"].get("phase") == "providerPacking":
                                    expect(page.get_by_label("Observed preparation progress")).to_contain_text(
                                        re.compile("provider packing", re.IGNORECASE), timeout=30000)
                            reattachment = {"attemptId": attempt_id,
                                            "status": reopened["attempt"]["status"],
                                            "phase": reopened["attempt"].get("phase"),
                                            "stageCount": len(reopened["stages"])}
                            rendered_phase_events.extend(page.evaluate(
                                "window.__renderedMemgenPhaseEvents || []"))
                            rendered_final_minimization |= attempt_id in page.evaluate(
                                "window.__renderedFinalMinimizationAttemptIds || []")
                            page.close()
                            page = None
                    if attempt.get("constructed") is not None or attempt.get("status") in (
                            "failed", "stopped", "resourceRefused", "unobserved"):
                        break
                    time.sleep(0.35)
                else:
                    raise AssertionError("Timed out before a checked Memgen handoff")
                assert attempt["attemptId"] == attempt_id
                assert attempt["constructed"] is not None, attempt
                rendered_phase_events.extend(page.evaluate(
                    "window.__renderedMemgenPhaseEvents || []"))
                rendered_final_minimization |= attempt_id in page.evaluate(
                    "window.__renderedFinalMinimizationAttemptIds || []")
                assert all(event["attemptId"] == attempt_id
                           for event in rendered_phase_events), rendered_phase_events
                rendered_phases = {event["phase"] for event in rendered_phase_events}
                assert "providerPacking" in seen_phases, seen_phases
                # Reattachment can occur during a brief provider phase. The
                # controlled browser phase test covers every label; here only
                # phases with a live page and an actual screenshot are required.
                required_rendered = {"provider Packing"}
                if "providerConditioningUnrestrained" in captured_phases:
                    required_rendered.add("provider Conditioning Unrestrained")
                assert required_rendered <= rendered_phases, (
                    "Provider phases missing from the actor view", sorted(rendered_phases))
                assert reattachment is not None and reattachment["attemptId"] == attempt_id
                assert not any(action["kind"] == "continueMinimization" for action in current["actions"])
                if attempt["status"] == "running":
                    assert current["stages"] == []
                checked_trials = [trial for trial in attempt["trials"]
                                  if trial["standing"] == "checked"]
                assert len(checked_trials) == 1, attempt["trials"]
                trial = checked_trials[0]
                assert trial["trialIndex"] >= 0
                assert {(item["physicalSide"], item["speciesId"]) for item in
                        trial["achievedLipidCounts"]} == set(expected)
                (folder / "observed-provider-phases.json").write_text(json.dumps({
                    "attemptId": attempt_id, "observedPhases": sorted(seen_phases),
                    "capturedPhases": sorted(captured_phases),
                    "renderedProviderPhases": sorted(rendered_phases),
                    "renderedPhaseEvents": rendered_phase_events,
                    "trials": attempt["trials"],
                    "browserReattachment": reattachment,
                }, indent=2) + "\n")
                print(f"Checked Memgen trial under attempt {attempt_id}; phases {sorted(seen_phases)}",
                      flush=True)
                minimization_or_stage = await_state(page, lambda state: (
                    ((state.get("attempt") or {}).get("status") == "running" and
                     (state.get("attempt") or {}).get("phase") == "finalMinimization") or
                    any(item["attemptId"] == attempt_id and item["kind"] == "Minimization" and
                        item["status"] == "completed" for item in state["stages"])),
                    "observable final minimization or completed stage", timeout=2400)
                if (minimization_or_stage["attempt"]["status"] == "running" and
                        minimization_or_stage["attempt"].get("phase") == "finalMinimization"):
                    seen_phases.add("finalMinimization")
                completed, stage = await_completed_stage(page, attempt_id, timeout=2400)
                rendered_final_minimization |= attempt_id in page.evaluate(
                    "window.__renderedFinalMinimizationAttemptIds || []")
                observed_phases = json.loads((folder / "observed-provider-phases.json").read_text())
                observed_phases["observedPhases"] = sorted(seen_phases)
                observed_phases["finalMinimizationRenderedWhileRunning"] = \
                    rendered_final_minimization
                (folder / "observed-provider-phases.json").write_text(
                    json.dumps(observed_phases, indent=2) + "\n")
                stage_id = stage["stageId"]
                assessment = stage["assessment"]
                assert stage["studyRevisionId"] == revision_id
                assert stage["constructed"]["subjectId"] == attempt["constructed"]["subjectId"]
                assert stage["observation"]["termination"] == "converged"
                assert assessment["checkStanding"] in STANDINGS
                assert assessment["currentlyApplicable"]
                page.locator(f'.stage-card[title^="Stage {stage_id} · "]').click()
                inspected = await_state(page, lambda state: (state.get("inspection") or {}).get(
                    "subjectId") == stage_id, "exact completed Memgen stage inspection")
                assert inspected["inspection"]["studyRevisionId"] == revision_id
                for width, height in ((1672, 941), (1024, 768), (820, 760)):
                    capture_review(page, folder, "export-ready", width, height,
                                   stage_id, assessment["id"])
                with page.expect_download(timeout=600000) as transfer:
                    page.get_by_role("button", name="Export with status").click()
                transfer.value.save_as(download_path)
                delivered = await_state(page, lambda state: (stage_in(state, stage_id).get(
                    "export") or {}).get("status") == "verified", "same-stage Memgen export")
                exported = stage_in(delivered, stage_id)["export"]
                assert exported["stageId"] == stage_id
                assert exported["assessmentId"] == assessment["id"]
                assert exported["sha256"] == sha256(download_path.read_bytes()).hexdigest()
                assert exported["byteLength"] == download_path.stat().st_size
                print(f"Completed and downloaded {stage_id}: {exported['sha256']}", flush=True)
            finally:
                try:
                    with urlopen(base + "/api/state", timeout=10) as response:
                        (folder / "final-account.json").write_text(json.dumps(
                            json.load(response), indent=2) + "\n")
                except Exception as failure:
                    (folder / "account-read-error.txt").write_text(str(failure) + "\n")
                browser.close()
        checked = subprocess.run([
            str(ROOT / "out/python/bin/python"), str(Path(__file__).resolve()),
            "--inspect-bundle", str(download_path), stage_id, attempt_id,
            revision_id, assessment["id"], assessment["checkStanding"], route_key,
            trial["trialId"], str(folder), source_case,
        ], cwd=ROOT, text=True, capture_output=True, timeout=600, check=False)
        self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
        evidence = json.loads(checked.stdout)
        assert evidence["bundleSha256"] == sha256(download_path.read_bytes()).hexdigest()
        (folder / "independent-bundle-readback.json").write_text(json.dumps(
            evidence, indent=2) + "\n")
        self.assertTrue(rendered_final_minimization,
            "The real OpenMM final minimization was absent from the running actor view")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--inspect-bundle":
    print(json.dumps(inspect_bundle(Path(sys.argv[2]), *sys.argv[3:10],
                                    Path(sys.argv[10]), sys.argv[11])))
elif __name__ == "__main__":
    unittest.main()
