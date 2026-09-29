"""Focused slice-1 checks of the real one-request worker and source observation."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src" / "ProteinInMembrane.Host"
PYTHON = os.environ.get("PIM_TEST_PYTHON", sys.executable)


def atom_line(serial: int, name: str, residue: str, chain: str, number: int,
              xyz: tuple[float, float, float], element: str) -> str:
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:4s} {residue:3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n")


def two_alanines(model_number: int, chain: str) -> str:
    atoms = [
        ("N", 1, (0.000, 0.000, 0.000), "N"),
        ("CA", 1, (1.460, 0.000, 0.000), "C"),
        ("C", 1, (2.020, 1.420, 0.000), "C"),
        ("O", 1, (1.340, 2.400, 0.000), "O"),
        ("CB", 1, (2.000, -0.770, -1.200), "C"),
        ("N", 2, (3.340, 1.520, 0.000), "N"),
        ("CA", 2, (4.020, 2.820, 0.000), "C"),
        ("C", 2, (5.530, 2.700, 0.000), "C"),
        ("O", 2, (6.100, 1.600, 0.000), "O"),
        ("CB", 2, (3.480, 3.600, -1.200), "C"),
        ("OXT", 2, (6.200, 3.720, 0.000), "O"),
    ]
    return (f"MODEL     {model_number:4d}\n" +
            "".join(atom_line(index, name, "ALA", chain, residue, xyz, element)
                    for index, (name, residue, xyz, element) in enumerate(atoms, 1)) +
            "TER\nENDMDL\n")


def invoke(directory: Path, operation: str, payload: dict, request_id: str = "slice-one"):
    request = {"requestId": request_id, "operation": operation,
               "workingDirectory": str(directory), "payload": payload}
    environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT), PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT))
    result = subprocess.run(
        [PYTHON, "-m", "ProteinInMembraneSystem.worker"],
        input=json.dumps(request) + "\n", text=True, capture_output=True,
        cwd=SOURCE_ROOT, env=environment, check=False, timeout=30)
    events = [json.loads(line) for line in result.stdout.splitlines()]
    return result, events


class WorkerExchangeTests(unittest.TestCase):
    def test_malformed_request_has_one_terminal_error(self):
        environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT),
                           PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT))
        result = subprocess.run(
            [PYTHON, "-m", "ProteinInMembraneSystem.worker"],
            input="{malformed\n", text=True, capture_output=True,
            cwd=SOURCE_ROOT, env=environment, check=False, timeout=30)
        events = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([event["kind"] for event in events], ["error"])
        self.assertEqual(events[0]["requestId"], "unidentified")
        self.assertEqual(events[0]["payload"]["failureCode"], "invalidRequest")

    def test_unsupported_request_is_correlated_and_terminal_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, events = invoke(Path(temporary), "unsupported_slice_one_operation", {})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([event["kind"] for event in events], ["progress", "error"])
        self.assertTrue(all(event["requestId"] == "slice-one" for event in events))
        self.assertEqual(events[-1]["payload"]["failureCode"], "unsupportedOperation")
        self.assertNotIn("qualification", events[-1]["payload"])

    def test_real_source_inspection_keeps_model_and_chain_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-models.pdb"
            source.write_text(two_alanines(1, "A") + two_alanines(2, "B") + "END\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            result, events = invoke(directory, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": digest,
                "maxAtoms": 25, "sourceKind": "upload", "prediction": None,
            })
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertEqual([event["kind"] for event in events], ["progress", "result"])
        observed = events[-1]["payload"]["observations"]
        self.assertEqual(observed["sourceFormat"], "pdb")
        self.assertEqual([model["index"] for model in observed["models"]], [0, 1])
        self.assertEqual([model["chains"][0]["name"] for model in observed["models"]], ["A", "B"])
        self.assertEqual([model["atomCount"] for model in observed["models"]], [11, 11])
        self.assertIsNone(observed["prediction"])

    def test_partner_display_names_preserve_codes_and_ignore_entity_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "named-partners.pdb"
            partner_atoms = (
                "HETATM   13  O   HOH B   3       9.000   4.000   0.000  1.00 20.00           O\n"
                "HETATM   14  FE  HEM B   4      10.000   4.000   0.000  1.00 20.00          FE\n"
                "HETATM   15  NA   NA B   5      11.000   4.000   0.000  1.00 20.00          NA\n"
                "HETATM   16  CL   CL B   6      12.000   4.000   0.000  1.00 20.00          CL\n"
                "HETATM   17  C1  GOL B   7      13.000   4.000   0.000  1.00 20.00           C\n"
            )
            source.write_text(two_alanines(1, "A").replace("ENDMDL\n", "") +
                              partner_atoms + "ENDMDL\nEND\n")
            result, events = invoke(directory, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "maxAtoms": 30, "sourceKind": "upload", "prediction": None,
            })
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            partners = events[-1]["payload"]["observations"]["models"][0]["partners"]
            self.assertEqual([(item["sourceId"], item["kind"], item["displayName"])
                              for item in partners], [
                ("0:B:3::HOH", "water", "Water"),
                ("0:B:4::HEM", "nonpolymer", "Heme"),
                ("0:B:5::NA", "ion", "Sodium ion"),
                ("0:B:6::CL", "ion", "Chloride ion"),
                ("0:B:7::GOL", "nonpolymer", "GOL"),
            ])

            # The authored mmCIF has an entity identifier, not a molecule name.
            # Make it numeric for an unknown compound without changing any address.
            fixture = ROOT / "tests" / "ProteinInMembraneSystem.AcceptanceTests" / \
                "SelectAndPrepareProtein" / "fixtures" / "actor-selection-identities.cif"
            cif_source = directory / "numeric-entity.cif"
            cif_source.write_text(fixture.read_text().replace("HEM!", "2").replace("HEM", "GOL"))
            result, events = invoke(directory, "inspect_source", {
                "sourcePath": str(cif_source),
                "sourceSha256": hashlib.sha256(cif_source.read_bytes()).hexdigest(),
                "maxAtoms": 30, "sourceKind": "upload", "prediction": None,
            })
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            for model in events[-1]["payload"]["observations"]["models"]:
                self.assertEqual({item["displayName"] for item in model["partners"]}, {"GOL"})
                self.assertTrue(all(item["label"] == "GOL" and item["sourceId"].endswith(":GOL")
                                    for item in model["partners"]))

    def test_real_source_reports_partner_and_requires_exact_alternate_location(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            first_cb = atom_line(5, "CB", "ALA", "A", 1, (2, -0.77, -1.2), "C")
            other_cb = atom_line(12, "CB", "ALA", "A", 1, (2, -0.77, 1.2), "C")
            alt_a = first_cb[:16] + "A" + first_cb[17:]
            alt_b = other_cb[:16] + "B" + other_cb[17:]
            partner = "HETATM   13  C1  GOL B   3       9.000   4.000   0.000  1.00 20.00           C\n"
            source = directory / "alternate-and-partner.pdb"
            source.write_text(two_alanines(1, "A").replace(first_cb, alt_a + alt_b)
                              .replace("ENDMDL\n", "") + partner + "ENDMDL\nEND\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            result, events = invoke(directory, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": digest,
                "maxAtoms": 25, "sourceKind": "upload", "prediction": None,
            })
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            model = events[-1]["payload"]["observations"]["models"][0]
            self.assertEqual(model["residues"][0]["alternateLocations"], ["A", "B"])
            self.assertEqual(model["partners"][0]["sourceId"], "0:B:3::GOL")
            self.assertEqual(model["partners"][0]["kind"], "nonpolymer")

            geometry_spec = json.loads((ROOT / "config" / "policies" /
                                        "protein-slice1.json").read_text())[
                                            "proteinStructuralPolicies"][0]["measurement"]
            selection = {
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "retainedPartners": [], "altlocChoices": [],
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                "geometrySpec": geometry_spec,
            }
            refused, refusal_events = invoke(directory, "inspect_preparation_changes", selection)
            self.assertNotEqual(refused.returncode, 0)
            self.assertEqual(refusal_events[-1]["payload"]["failureCode"], "approvalRequired")
            selection["altlocChoices"] = [{
                "residue": {"model": 0, "chain": "A", "residue": 1,
                            "insertionCode": "", "copyId": "A"},
                "altloc": "A",
            }]
            selected, selected_events = invoke(directory, "inspect_preparation_changes", selection)
            self.assertEqual(selected.returncode, 0, selected.stderr + selected.stdout)
            observed = selected_events[-1]["payload"]["observations"]
            self.assertEqual(observed["selectedAtomCount"], 11)
            self.assertEqual(observed["previewChains"], [{
                "sourceChain": "A", "copyId": "A", "previewChain": "A",
            }])
            self.assertEqual(observed["assessmentStanding"], "Observed")
            selection["retainedPartners"] = model["partners"]
            unsupported, unsupported_events = invoke(directory, "inspect_preparation_changes", selection)
            self.assertNotEqual(unsupported.returncode, 0)
            self.assertEqual(unsupported_events[-1]["payload"]["failureCode"], "unsupportedChemistry")
            selection["retainedPartners"] = [dict(model["partners"][0], insertionCode="A")]
            wrong_address, wrong_address_events = invoke(directory, "inspect_preparation_changes", selection)
            self.assertNotEqual(wrong_address.returncode, 0)
            self.assertEqual(wrong_address_events[-1]["payload"]["failureCode"], "invalidSelection")

    def test_missing_backbone_is_explicitly_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "missing-backbone.pdb"
            source.write_text(two_alanines(1, "A").replace(
                atom_line(4, "O", "ALA", "A", 1, (1.34, 2.4, 0), "O"), "") + "END\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            result, events = invoke(directory, "inspect_preparation_changes", {
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "retainedPartners": [], "altlocChoices": [],
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
            })
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(events[-1]["payload"]["failureCode"], "unsupportedBackbone")

    def test_noncanonical_polymer_is_not_silently_parameterized(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "noncanonical.pdb"
            source.write_text(two_alanines(1, "A").replace(
                "ALA A   1", "MSE A   1") + "END\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            result, events = invoke(directory, "inspect_preparation_changes", {
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "retainedPartners": [], "altlocChoices": [],
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
            })
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(events[-1]["payload"]["failureCode"], "unsupportedChemistry")

    def test_hash_mismatch_and_escaping_source_path_do_not_open_a_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            wrong_hash_result, wrong_hash_events = invoke(directory, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": "0" * 64, "maxAtoms": 25,
                "sourceKind": "upload", "prediction": None,
            })
            outside = directory.parent / "outside-slice-one.pdb"
            escape_result, escape_events = invoke(directory, "inspect_source", {
                "sourcePath": str(outside), "sourceSha256": "0" * 64, "maxAtoms": 25,
                "sourceKind": "upload", "prediction": None,
            })
        self.assertNotEqual(wrong_hash_result.returncode, 0)
        self.assertEqual(wrong_hash_events[-1]["payload"]["failureCode"], "inputMismatch")
        self.assertNotEqual(escape_result.returncode, 0)
        self.assertEqual(escape_events[-1]["payload"]["failureCode"], "invalidPath")
        self.assertEqual(sum(event["kind"] == "error" for event in wrong_hash_events), 1)
        self.assertEqual(sum(event["kind"] == "error" for event in escape_events), 1)

    def test_real_preparation_preserves_selected_heavy_atoms_and_observes_geometry(self):
        catalogue = json.loads((ROOT / "config" / "policies" / "protein-slice1.json").read_text())
        geometry_spec = catalogue["proteinStructuralPolicies"][0]["measurement"]
        force_field_asset = dict(catalogue["proteinChemicalStates"][0]["forceFieldFiles"][0])
        installed_asset = (ROOT / "config" / "policies" / force_field_asset["path"]).resolve()
        self.assertEqual(hashlib.sha256(installed_asset.read_bytes()).hexdigest(),
                         force_field_asset["sha256"])

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines.pdb"
            source.write_text(two_alanines(1, "A") + "END\n")
            source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            staged_asset = directory / "protein.ff19SB.xml"
            shutil.copyfile(installed_asset, staged_asset)
            force_field_asset["path"] = str(staged_asset)
            result, events = invoke(directory, "prepare_protein", {
                "studyRevisionId": "known-two-alanine-revision",
                "nominalPh": 7.0,
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                "sourcePath": str(source), "sourceSha256": source_hash,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "altlocChoices": [], "retainedPartners": [],
                "heavyAtomApprovals": [], "approvedDisulfides": [], "residueVariants": [],
                "forceFieldFiles": [force_field_asset], "geometrySpec": geometry_spec,
            })
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual([event["kind"] for event in events],
                             ["progress", "progress", "progress", "result"])
            terminal = events[-1]["payload"]
            self.assertEqual(terminal["studyRevisionId"], "known-two-alanine-revision")
            observed = terminal["observations"]
            # The authored source has 2 ALA residues, 11 heavy atoms, and no
            # alternate sites or disulfides. Terminal hydrogens add 12 atoms.
            self.assertEqual((observed["sourceAtomCount"], observed["preparedAtomCount"],
                              observed["retainedResidueCount"]), (11, 23, 2))
            self.assertEqual(observed["addedHeavyAtoms"], [])
            self.assertEqual(observed["actualDisulfides"], [])
            self.assertEqual(observed["geometryWarnings"], [])
            self.assertEqual({item["kind"] for item in observed["geometry"]["kinds"]},
                             {"covalentBond", "chainContinuity"})
            self.assertTrue(all(item["standing"] == "Observed" and
                                item["eligibleCount"] == item["measuredCount"] > 0
                                for item in observed["geometry"]["kinds"]))
            artifacts = {item["role"]: item for item in terminal["artifacts"]}
            self.assertEqual(set(artifacts),
                             {"preparedPdb", "preparedBondGraph", "correspondenceJson"})
            for artifact in artifacts.values():
                self.assertEqual(hashlib.sha256(Path(artifact["path"]).read_bytes()).hexdigest(),
                                 artifact["sha256"])
            correspondence = json.loads(Path(artifacts["correspondenceJson"]["path"]).read_text())
            self.assertTrue(correspondence["complete"])
            self.assertEqual(correspondence["sourceId"], source_hash)
            self.assertEqual(correspondence["resultId"], artifacts["preparedPdb"]["sha256"])
            self.assertEqual(len(correspondence["atoms"]), 23)
            self.assertEqual(sum(item["sourceAtomId"] is not None
                                 for item in correspondence["atoms"]), 11)
            bond_graph = json.loads(Path(artifacts["preparedBondGraph"]["path"]).read_text())
            self.assertEqual(len(bond_graph["atoms"]), 23)
            self.assertEqual(len(bond_graph["residues"]), 2)
            def atom_identity(index):
                atom = bond_graph["atoms"][index]
                residue = bond_graph["residues"][atom["residueIndex"]]
                return residue["id"], atom["name"]
            peptide_bonds = [bond for bond in bond_graph["bonds"] if
                             {atom_identity(index) for index in bond["atomIndices"]} ==
                             {("1", "C"), ("2", "N")}]
            self.assertEqual(len(peptide_bonds), 1)

    def test_real_heavy_atom_repair_needs_exact_approval_and_records_correspondence(self):
        catalogue = json.loads((ROOT / "config" / "policies" / "protein-slice1.json").read_text())
        asset = dict(catalogue["proteinChemicalStates"][0]["forceFieldFiles"][0])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "missing-sidechain.pdb"
            source.write_text(two_alanines(1, "A").replace(
                atom_line(5, "CB", "ALA", "A", 1, (2, -0.77, -1.2), "C"), "") + "END\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            staged_asset = directory / "protein.ff19SB.xml"
            shutil.copyfile((ROOT / "config" / "policies" / asset["path"]).resolve(),
                            staged_asset)
            asset["path"] = str(staged_asset)
            selected_residue = {"model": 0, "chain": "A", "residue": 1,
                                "insertionCode": "", "copyId": "A"}
            payload = {
                "studyRevisionId": "repair-revision", "nominalPh": 7.0,
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "altlocChoices": [], "retainedPartners": [],
                "approvedHeavyAtoms": [], "approvedDisulfides": [],
                "residueVariants": [], "forceFieldFiles": [asset],
                "geometrySpec": catalogue["proteinStructuralPolicies"][0]["measurement"],
            }
            inspected, inspection_events = invoke(directory, "inspect_preparation_changes", payload)
            self.assertEqual(inspected.returncode, 0, inspected.stderr + inspected.stdout)
            self.assertEqual(inspection_events[-1]["payload"]["observations"][
                "missingNonbackboneHeavyAtoms"], [{
                    "residue": selected_residue, "atomName": "CB",
                }])
            unapproved, unapproved_events = invoke(directory, "prepare_protein", payload)
            self.assertNotEqual(unapproved.returncode, 0)
            self.assertEqual(unapproved_events[-1]["payload"]["failureCode"], "approvalRequired")
            self.assertEqual(sum(event["kind"] == "error" for event in unapproved_events), 1)

            payload["approvedHeavyAtoms"] = [{
                "residue": selected_residue, "atomName": "CB", "decisionId": "approval-1",
            }]
            prepared, prepared_events = invoke(directory, "prepare_protein", payload)
            self.assertEqual(prepared.returncode, 0, prepared.stderr + prepared.stdout)
            terminal = prepared_events[-1]["payload"]
            self.assertEqual(terminal["observations"]["addedHeavyAtoms"], [{
                "residue": selected_residue, "atomName": "CB",
            }])
            mapping_path = next(Path(item["path"]) for item in terminal["artifacts"]
                                if item["role"] == "correspondenceJson")
            mapping = json.loads(mapping_path.read_text())
            added = [atom for atom in mapping["atoms"]
                     if atom["approvedChangeId"] == "approval-1"]
            self.assertEqual(len(added), 1)
            self.assertEqual(added[0]["role"], "generated")
            self.assertEqual(added[0]["element"], "C")

    def test_real_variable_residue_requires_and_reports_exact_chemical_state(self):
        catalogue = json.loads((ROOT / "config" / "policies" / "protein-slice1.json").read_text())
        asset = dict(catalogue["proteinChemicalStates"][0]["forceFieldFiles"][0])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "cysteine-and-alanine.pdb"
            first_second_residue_atom = atom_line(6, "N", "ALA", "A", 2,
                                                  (3.34, 1.52, 0), "N")
            sulfur = atom_line(12, "SG", "CYS", "A", 1,
                               (3.8, -0.7, -1.2), "S")
            source.write_text(two_alanines(1, "A").replace(
                "ALA A   1", "CYS A   1").replace(
                    first_second_residue_atom,
                    sulfur + first_second_residue_atom) + "END\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            staged_asset = directory / "protein.ff19SB.xml"
            shutil.copyfile((ROOT / "config" / "policies" / asset["path"]).resolve(),
                            staged_asset)
            asset["path"] = str(staged_asset)
            selected_residue = {"model": 0, "chain": "A", "residue": 1,
                                "insertionCode": "", "copyId": "A"}
            payload = {
                "studyRevisionId": "explicit-cysteine-state", "nominalPh": 7.0,
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "altlocChoices": [], "retainedPartners": [],
                "approvedHeavyAtoms": [], "approvedDisulfides": [],
                "residueVariants": [], "forceFieldFiles": [asset],
                "geometrySpec": catalogue["proteinStructuralPolicies"][0]["measurement"],
            }
            unresolved, unresolved_events = invoke(directory, "prepare_protein", payload)
            self.assertNotEqual(unresolved.returncode, 0)
            self.assertEqual(unresolved_events[-1]["payload"]["failureCode"],
                             "chemicalStateUnresolved")
            payload["residueVariants"] = [{
                "residue": selected_residue, "variant": "CYS", "decisionId": None,
            }]
            prepared, prepared_events = invoke(directory, "prepare_protein", payload)
            self.assertEqual(prepared.returncode, 0, prepared.stderr + prepared.stdout)
            observed = prepared_events[-1]["payload"]["observations"]
            self.assertEqual(observed["actualResidueVariants"], [{
                "residue": selected_residue, "variant": "CYS", "decisionId": None,
            }])
            self.assertEqual(observed["actualDisulfides"], [])
            self.assertEqual(observed["geometryWarnings"], [])

    def test_real_histidine_choice_reaches_provider_with_corresponding_result(self):
        """6QWR model 1 is an identified provider fixture, not a new biological claim."""
        catalogue = json.loads((ROOT / "config" / "policies" / "protein-slice1.json").read_text())
        policy = catalogue["proteinChemicalStates"][0]
        asset = dict(policy["forceFieldFiles"][0])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "6QWR.pdb"
            shutil.copyfile(ROOT / "config" / "policies" / "source-assets" / "6QWR.pdb", source)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            staged_asset = directory / "protein.ff19SB.xml"
            shutil.copyfile((ROOT / "config" / "policies" / asset["path"]).resolve(),
                            staged_asset)
            asset["path"] = str(staged_asset)
            def address(number, insertion_code=""):
                return {"model": 0, "chain": "A", "residue": number,
                        "insertionCode": insertion_code, "copyId": "A"}
            variants = []
            seen = set()
            for line in source.read_text().splitlines():
                if line.startswith("ENDMDL"):
                    break
                if not line.startswith("ATOM"):
                    continue
                name = line[17:20].strip()
                number = int(line[22:26])
                insertion_code = line[26].strip()
                key = number, insertion_code
                if key in seen:
                    continue
                seen.add(key)
                if name in policy["defaultVariants"]:
                    variants.append({"residue": address(number, insertion_code),
                                     "variant": policy["defaultVariants"][name],
                                     "decisionId": None})
                elif name == "HIS":
                    variants.append({"residue": address(number, insertion_code),
                                     "variant": "HID", "decisionId": "fixture-HIS-108"})
            self.assertEqual([item for item in variants if item["variant"] == "HID"], [{
                "residue": address(108), "variant": "HID", "decisionId": "fixture-HIS-108"}])
            result, events = invoke(directory, "prepare_protein", {
                "studyRevisionId": "real-histidine-review", "nominalPh": 7.0,
                "disulfideCandidateMaxSgDistanceAngstrom": 2.5,
                "sourcePath": str(source), "sourceSha256": digest,
                "modelIndex": 0, "assemblyId": None,
                "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
                "altlocChoices": [], "retainedPartners": [],
                "approvedHeavyAtoms": [{"residue": address(211), "atomName": "OXT",
                                        "decisionId": "fixture-OXT-211"}],
                "approvedDisulfides": [], "residueVariants": variants,
                "forceFieldFiles": [asset],
                "geometrySpec": catalogue["proteinStructuralPolicies"][0]["measurement"],
            }, request_id="real-histidine-review")
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            terminal = events[-1]["payload"]
            self.assertEqual([event["kind"] for event in events],
                             ["progress", "progress", "progress", "result"])
            self.assertEqual(terminal["studyRevisionId"], "real-histidine-review")
            observed = terminal["observations"]
            # The provider reports actual chemistry; the root retains decision IDs.
            self.assertIn({"residue": address(108), "variant": "HID",
                           "decisionId": None}, observed["actualResidueVariants"])
            self.assertEqual(observed["addedHeavyAtoms"],
                             [{"residue": address(211), "atomName": "OXT"}])
            self.assertEqual(observed["preparedAtomCount"], 3205)
            self.assertEqual(observed["geometryWarnings"], [])
            artifacts = {item["role"]: item for item in terminal["artifacts"]}
            self.assertEqual(set(artifacts),
                             {"preparedPdb", "preparedBondGraph", "correspondenceJson"})
            for item in artifacts.values():
                self.assertEqual(hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest(),
                                 item["sha256"])
            correspondence = json.loads(Path(artifacts["correspondenceJson"]["path"]).read_text())
            self.assertTrue(correspondence["complete"])
            self.assertEqual(correspondence["sourceId"], digest)
            self.assertEqual(correspondence["resultId"], artifacts["preparedPdb"]["sha256"])


if __name__ == "__main__":
    unittest.main()
