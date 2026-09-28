"""Real bounded PDBFixer/OpenMM recommendation observations on small protein fixtures."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from openmm.app import ForceField, Modeller, PDBFile
from pdbfixer import PDBFixer

from test_worker_exchange import ROOT, atom_line, invoke, two_alanines


CATALOGUE = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())


def fixture_payload(directory: Path, source: Path) -> dict:
    chemical = CATALOGUE["proteinChemicalStates"][0]
    asset = dict(chemical["forceFieldFiles"][0])
    installed = (ROOT / "config/policies" / asset["path"]).resolve()
    staged = directory / "protein.ff19SB.xml"
    shutil.copyfile(installed, staged)
    asset["path"] = str(staged)
    return {"studyRevisionId": "identified-recommendation-fixture",
            "sourcePath": str(source), "sourceSha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "modelIndex": 0, "assemblyId": None,
            "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
            "altlocChoices": [], "retainedPartners": [], "approvedDisulfides": [],
            "overrides": [], "nominalPh": 7.0, "seed": chemical["recommendation"]["seed"],
            "forceFieldFiles": [asset], "permittedVariants": chemical["permittedVariants"],
            "disulfideCandidateMaxSgDistanceAngstrom": chemical["disulfideCandidateMaxSgDistanceAngstrom"],
            "geometrySpec": CATALOGUE["proteinStructuralPolicies"][0]["measurement"]}


def prepared_histidines(path: Path) -> None:
    starting = path.with_name("starting-alanines.pdb")
    starting.write_text(two_alanines(1, "A") + "END\n")
    fixer = PDBFixer(filename=str(starting))
    fixer.applyMutations(["ALA-1-HIS", "ALA-2-HIS"], "A")
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with path.open("w") as stream:
        PDBFile.writeFile(fixer.topology, fixer.positions, stream, keepIds=True)


def internal_variant_fragment(path: Path, variant_name: str) -> None:
    """Use only the first 6QWR model's observed 106–110 peptide geometry."""
    first_model = (ROOT / "config/policies/source-assets/6QWR.pdb").read_text().split("ENDMDL")[0]
    lines = [line for line in first_model.splitlines() if line.startswith("ATOM") and
             line[21] == "A" and 106 <= int(line[22:26]) <= 110]
    starting = path.with_name("first-model-peptide.pdb")
    starting.write_text("\n".join(lines) + "\nTER\nEND\n")
    fixer = PDBFixer(filename=str(starting))
    fixer.applyMutations([f"HIS-108-{variant_name}"], "A")
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with path.open("w") as stream:
        PDBFile.writeFile(fixer.topology, fixer.positions, stream, keepIds=True)


class RecommendationWorkerTests(unittest.TestCase):
    def test_authorized_seed_reproduces_checked_candidate_with_heavy_repair(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "missing-sidechain.pdb"
            missing_cb = atom_line(5, "CB", "ALA", "A", 1, (2, -0.77, -1.2), "C")
            source.write_text(two_alanines(1, "A").replace(missing_cb, "") + "END\n")
            payload = fixture_payload(directory, source)
            inspected, events = invoke(directory, "inspect_preparation_changes", payload)
            self.assertEqual(inspected.returncode, 0, events[-1])
            payload["proposedHeavyAtoms"] = events[-1]["payload"]["observations"][
                "missingNonbackboneHeavyAtoms"]
            self.assertEqual(len(payload["proposedHeavyAtoms"]), 1)
            proposed, plan_events = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(proposed.returncode, 0, plan_events[-1])
            plan = plan_events[-1]["payload"]["observations"]
            payload["approvedHeavyAtoms"] = [{"residue": atom["residue"],
                                                "atomName": atom["atomName"],
                                                "decisionId": "exact-repair-approval"}
                                               for atom in payload["proposedHeavyAtoms"]]
            payload["residueVariants"] = [{"residue": choice["residue"],
                                           "variant": choice["variant"],
                                           "decisionId": "whole-plan-authorization"}
                                          for choice in plan["choices"]]
            payload["normalizeProteinHydrogens"] = True
            payload["planSeed"] = payload["seed"]
            actual, actual_events = invoke(directory, "prepare_protein", payload)
            self.assertEqual(actual.returncode, 0, actual_events[-1])
            self.assertEqual(next(item["sha256"] for item in actual_events[-1]["payload"]["artifacts"]
                                  if item["role"] == "preparedPdb"), plan["candidateSha256"])

    def test_reviewed_disulfide_pairs_are_preserved_in_joint_candidate(self):
        # RCSB 1CRN, downloaded 2026-09-28 from files.rcsb.org/download/1CRN.pdb.
        # The three observed sulfur pairs are provider observations, not an automatic approval.
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "1CRN-RCSB.pdb"
            shutil.copyfile(ROOT / "tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein" /
                            "fixtures/1CRN-RCSB.pdb", source)
            payload = fixture_payload(directory, source)
            inspected, events = invoke(directory, "inspect_preparation_changes", payload)
            self.assertEqual(inspected.returncode, 0, events[-1])
            pairs = events[-1]["payload"]["observations"]["possibleDisulfides"]
            self.assertEqual(len(pairs), 3)
            payload["proposedHeavyAtoms"] = []
            payload["approvedDisulfides"] = [{"first": pair["first"], "second": pair["second"],
                                               "decisionId": f"exact-bond-approval-{index}"}
                                              for index, pair in enumerate(pairs)]
            recommended, outcome = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(recommended.returncode, 0, outcome[-1])
            actual = outcome[-1]["payload"]["observations"]["candidateObservations"]
            expected = {(json.dumps(pair["first"], sort_keys=True),
                         json.dumps(pair["second"], sort_keys=True)) for pair in pairs}
            self.assertEqual({(json.dumps(pair["first"], sort_keys=True),
                               json.dumps(pair["second"], sort_keys=True))
                              for pair in actual["actualDisulfides"]}, expected)
            # A declined observed pair cannot be silently reintroduced by PDB loading.
            payload["approvedDisulfides"] = payload["approvedDisulfides"][1:]
            refused, refusal_events = invoke(directory, "recommend_preparation", payload)
            self.assertNotEqual(refused.returncode, 0)
            self.assertEqual(refusal_events[-1]["payload"]["failureCode"], "approvalRequired")

    def test_reviewed_conformer_is_kept_in_joint_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-alanines-alternate.pdb"
            original = atom_line(5, "CB", "ALA", "A", 1, (2, -0.77, -1.2), "C")
            shifted = atom_line(12, "CB", "ALA", "A", 1, (2.1, -0.77, -1.2), "C")
            source.write_text(two_alanines(1, "A").replace(
                original, original[:16] + "A" + original[17:] +
                shifted[:16] + "B" + shifted[17:]) + "END\n")
            payload = fixture_payload(directory, source)
            payload["altlocChoices"] = [{"residue": {"model": 0, "chain": "A", "copyId": "A",
                                                        "residue": 1, "insertionCode": ""},
                                        "altloc": "A", "decisionId": "exact-conformer-approval"}]
            payload["proposedHeavyAtoms"] = []
            result, events = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(result.returncode, 0, events[-1])
            self.assertEqual(events[-1]["payload"]["observations"]["candidateObservations"]
                             ["unresolvedAlternateLocations"], [])

    def test_internal_acidic_and_basic_overrides_use_exact_templates(self):
        for ordinary, alternative in (("ASP", "ASH"), ("GLU", "GLH"), ("LYS", "LYN")):
            with self.subTest(ordinary=ordinary, alternative=alternative), \
                 tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                source = directory / f"internal-{ordinary.lower()}.pdb"
                internal_variant_fragment(source, ordinary)
                payload = fixture_payload(directory, source)
                inspected, events = invoke(directory, "inspect_preparation_changes", payload)
                self.assertEqual(inspected.returncode, 0, events[-1])
                payload["proposedHeavyAtoms"] = events[-1]["payload"]["observations"][
                    "missingNonbackboneHeavyAtoms"]
                payload["overrides"] = [{"residue": {"model": 0, "chain": "A", "copyId": "A",
                                                     "residue": 108, "insertionCode": ""},
                                         "variant": alternative}]
                result, outcome = invoke(directory, "recommend_preparation", payload)
                self.assertEqual(result.returncode, 0, outcome[-1])
                plan = outcome[-1]["payload"]["observations"]
                self.assertEqual(next(item for item in plan["choices"] if
                                      item["residue"]["residue"] == 108)["variant"], alternative)
                self.assertTrue(plan["jointParameterizationObserved"])

    def test_method_variants_override_and_joint_candidate_are_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "two-histidines.pdb"
            prepared_histidines(source)
            source_bytes = source.read_bytes()
            payload = fixture_payload(directory, source)
            inspected, inspection_events = invoke(directory, "inspect_preparation_changes", payload)
            self.assertEqual(inspected.returncode, 0, inspection_events[-1])
            payload["proposedHeavyAtoms"] = inspection_events[-1]["payload"]["observations"][
                "missingNonbackboneHeavyAtoms"]
            first, first_events = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(first.returncode, 0, first_events[-1])
            first_plan = first_events[-1]["payload"]["observations"]
            self.assertEqual(len(first_plan["choices"]), 2)
            self.assertTrue(first_plan["jointParameterizationObserved"])
            self.assertEqual(source.read_bytes(), source_bytes)
            self.assertEqual(first_plan["candidateSha256"], next(item["sha256"] for item in
                first_events[-1]["payload"]["artifacts"] if item["role"] == "preparedPdb"))
            self.assertEqual({item["variant"] for item in first_plan["choices"]}, {"HID"})
            repeated, repeated_events = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(repeated.returncode, 0, repeated_events[-1])
            self.assertEqual(first_plan["planSha256"],
                             repeated_events[-1]["payload"]["observations"]["planSha256"])
            override = first_plan["choices"][0]
            for variant in ("HID", "HIE", "HIP"):
                payload["overrides"] = [{"residue": override["residue"], "variant": variant}]
                second, second_events = invoke(directory, "recommend_preparation", payload)
                self.assertEqual(second.returncode, 0, second_events[-1])
                second_plan = second_events[-1]["payload"]["observations"]
                self.assertNotEqual(first_plan["planSha256"], second_plan["planSha256"])
                self.assertEqual(second_plan["choices"][0]["variant"], variant)
                self.assertTrue(second_plan["choices"][0]["overridden"])
                self.assertEqual(second_plan["choices"][1]["variant"], "HID")
                self.assertEqual({(item["residue"]["residue"], item["variant"]) for item in
                                  second_plan["candidateObservations"]["actualResidueVariants"]},
                                 {(1, variant), (2, "HID")})
            self.assertEqual(source.read_bytes(), source_bytes)

    def test_source_protein_hydrogens_are_enumerated_and_regenerated(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source-with-hydrogens.pdb"
            naked = directory / "starting.pdb"
            naked.write_text(two_alanines(1, "A") + "END\n")
            asset = CATALOGUE["proteinChemicalStates"][0]["forceFieldFiles"][0]
            force_field = ForceField(str((ROOT / "config/policies" / asset["path"]).resolve()))
            loaded = PDBFile(str(naked))
            modeller = Modeller(loaded.topology, loaded.positions)
            modeller.addHydrogens(force_field, pH=7.0)
            with source.open("w") as stream:
                PDBFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
            source_bytes = source.read_bytes()
            payload = fixture_payload(directory, source)
            payload["proposedHeavyAtoms"] = []
            result, events = invoke(directory, "recommend_preparation", payload)
            self.assertEqual(result.returncode, 0, events[-1])
            observations = events[-1]["payload"]["observations"]
            self.assertGreater(len(observations["removedSourceHydrogens"]), 0)
            self.assertEqual({json.dumps(item, sort_keys=True) for item in
                              observations["candidateObservations"]["removedSourceHydrogens"]},
                             {json.dumps(item, sort_keys=True) for item in
                              observations["removedSourceHydrogens"]})
            correspondence = json.loads(Path(next(item["path"] for item in
                events[-1]["payload"]["artifacts"] if item["role"] == "correspondenceJson")).read_text())
            source_hydrogen_names = {item["atomName"] for item in observations["removedSourceHydrogens"]}
            regenerated = [item for item in correspondence["atoms"] if item["element"] == "H" and
                           item["resultAtomId"].split(":")[-1] in source_hydrogen_names]
            self.assertTrue(regenerated)
            self.assertTrue(all(item["role"] == "generated" and item["sourceAtomId"] is None
                                for item in regenerated))
            self.assertEqual(source.read_bytes(), source_bytes)


if __name__ == "__main__":
    unittest.main()
