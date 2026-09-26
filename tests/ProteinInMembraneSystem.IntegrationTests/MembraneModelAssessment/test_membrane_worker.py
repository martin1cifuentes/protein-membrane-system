"""Slice-2 observations through the real one-request OpenMM worker process."""

from __future__ import annotations

import hashlib
import importlib.metadata
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
OPENMM_DATA = Path(importlib.metadata.distribution("openmm").locate_file("openmm/app/data"))
EXPECTED_SOURCE_HASHES = {
    "POPC": "a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29",
    "lipid21": "4ed8dc7485a3df14f55e1d9d12fb2a99813495293ba07c2d910fe0f66185e23a",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_molecule(source: Path, output: Path, species: str) -> tuple[int, int]:
    """Preserve a provider sample's first exact molecule, including bonds and coordinates."""
    from openmm import unit
    from openmm.app import PDBFile, Topology

    sample = PDBFile(str(source))
    first = next(sample.topology.residues())
    topology = Topology()
    residue = topology.addResidue(species, topology.addChain("A"), id="1")
    copied = {atom: topology.addAtom(atom.name, atom.element, residue)
              for atom in first.atoms()}
    for left, right in sample.topology.bonds():
        if left.residue == first and right.residue == first:
            topology.addBond(copied[left], copied[right])
    positions = unit.Quantity(
        [sample.positions[atom.index].value_in_unit(unit.nanometer)
         for atom in first.atoms()], unit.nanometer)
    with output.open("w", encoding="ascii") as stream:
        PDBFile.writeFile(topology, positions, stream, keepIds=True)
    round_trip = PDBFile(str(output))
    return len(list(round_trip.topology.atoms())), len(list(round_trip.topology.bonds()))


def invoke(directory: Path, payload: dict, request_id: str = "membrane-slice-two"):
    request = {"requestId": request_id, "operation": "assess_membrane",
               "workingDirectory": str(directory), "payload": payload}
    environment = dict(os.environ, PYTHONPATH=str(SOURCE_ROOT),
                       PIM_OWNER_SOURCE_ROOT=str(SOURCE_ROOT))
    result = subprocess.run(
        [PYTHON, "-m", "ProteinInMembraneSystem.worker"],
        input=json.dumps(request) + "\n", text=True, capture_output=True,
        cwd=SOURCE_ROOT, env=environment, check=False, timeout=60)
    return result, [json.loads(line) for line in result.stdout.splitlines()]


def assets(directory: Path, names: tuple[str, ...]) -> dict:
    source_ff = OPENMM_DATA / "amber19" / "lipid21.xml"
    if "POPC" in names:
        assert digest(OPENMM_DATA / "POPC.pdb") == EXPECTED_SOURCE_HASHES["POPC"]
    assert digest(source_ff) == EXPECTED_SOURCE_HASHES["lipid21"]
    force_field = directory / "lipid21.xml"
    shutil.copyfile(source_ff, force_field)
    force_field_hash = digest(force_field)
    representations = []
    counts = {"POPC": 134, "POPE": 125}
    for name in names:
        coordinate = directory / f"{name}-one-molecule.pdb"
        actual_count, bonds = one_molecule(OPENMM_DATA / f"{name}.pdb", coordinate, name)
        assert actual_count == counts[name]
        assert bonds == counts[name] - 1
        representations.append({
            "speciesId": name,
            "chemistryId": f"OpenMM-8.6.0-Lipid21-{name}",
            "category": "lipid",
            "templatePath": str(force_field), "templateSha256": force_field_hash,
            "coordinateTemplatePath": str(coordinate),
            "coordinateTemplateSha256": digest(coordinate),
            "atomCount": actual_count, "netChargeElementary": 0,
            # Positive values satisfy the worker's schema guards. They are
            # test inputs, not claimed area/volume policy evidence.
            "areaPerMoleculeAngstromSquared": 1.0,
            "volumeAngstromCubed": 1.0,
            "headAtomIndices": [1, 20],
            "forceFieldFamily": "Lipid21", "forceFieldVersion": "8.6.0",
            "limitations": [],
            "stereoChecks": [
                {"kind": "tetrahedral", "atomNames": ["O21", "C1", "C3", "HS"],
                 "expected": "negative"},
                {"kind": "alkene", "atomNames": ["C28", "C29", "C210", "C211"],
                 "expected": "cis"},
            ],
        })
    return {
        "studyRevisionId": "revision-2", "membraneModelId": "model-2",
        "upper": {"physicalSide": "upper", "fractions": [{"speciesId": names[0], "fraction": 1.0}]},
        "lower": {"physicalSide": "lower", "fractions": [{"speciesId": names[0], "fraction": 1.0}]},
        "speciesRepresentations": representations,
        "forceFieldFiles": [{"id": "openmm-8.6.0-lipid21", "version": "8.6.0",
                             "family": "Lipid21", "path": str(force_field),
                             "sha256": force_field_hash}],
    }


def catalogue_selection(directory: Path, names: tuple[str, ...]) -> dict:
    catalogue_path = ROOT / "config" / "policies" / "protein-membrane-slice2.json"
    catalogue = json.loads(catalogue_path.read_text(encoding="utf-8"))
    by_species = {item["speciesId"]: item for item in catalogue["lipids"]}
    source_base = catalogue_path.parent
    selected = []
    for species in names:
        item = dict(by_species[species])
        source_coordinate = (source_base / item["coordinateTemplatePath"]).resolve()
        assert digest(source_coordinate) == item["coordinateTemplateSha256"]
        staged_coordinate = directory / f"{species}{source_coordinate.suffix}"
        shutil.copyfile(source_coordinate, staged_coordinate)
        item["coordinateTemplatePath"] = str(staged_coordinate)
        selected.append(item)
    source_force_field = (source_base / selected[0]["templatePath"]).resolve()
    assert all(source_force_field == (source_base / item["templatePath"]).resolve()
               for item in selected)
    assert all(digest(source_force_field) == item["templateSha256"] for item in selected)
    staged_force_field = directory / "lipid21.xml"
    shutil.copyfile(source_force_field, staged_force_field)
    for item in selected:
        item["templatePath"] = str(staged_force_field)
    return {
        "studyRevisionId": "catalogue-revision", "membraneModelId": "catalogue-model",
        "upper": {"physicalSide": "upper", "fractions": [
            {"speciesId": names[0], "fraction": 1.0},
        ]},
        "lower": {"physicalSide": "lower", "fractions": [
            {"speciesId": names[0], "fraction": 1.0},
        ]},
        "speciesRepresentations": selected,
        "forceFieldFiles": [{"id": "openmm-8.6.0-lipid21", "version": "8.6.0",
                             "family": "Lipid21", "path": str(staged_force_field),
                             "sha256": digest(staged_force_field)}],
    }


def catalogue_mixture(directory: Path) -> dict:
    payload = catalogue_selection(directory, ("POPC", "CHL1"))
    payload["studyRevisionId"] = "cholesterol-revision"
    payload["membraneModelId"] = "cholesterol-model"
    payload["upper"]["fractions"] = [
        {"speciesId": "POPC", "fraction": 0.60},
        {"speciesId": "CHL1", "fraction": 0.40},
    ]
    payload["lower"]["fractions"] = [
        {"speciesId": "POPC", "fraction": 0.40},
        {"speciesId": "CHL1", "fraction": 0.60},
    ]
    return payload


class MembraneWorkerTests(unittest.TestCase):
    def test_every_qualified_catalogue_species_has_exact_real_worker_support(self):
        expected_counts = {
            "DLPC": 106, "DLPE": 97, "DMPC": 118, "DPPC": 130,
            "DOPC": 138, "POPC": 134, "POPE": 125, "CHL1": 74,
        }
        catalogue = json.loads((ROOT / "config" / "policies" /
                                "protein-membrane-slice2.json").read_text(encoding="utf-8"))
        self.assertEqual({item["speciesId"] for item in catalogue["lipids"]},
                         set(expected_counts))
        for species, atom_count in expected_counts.items():
            with self.subTest(species=species), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                payload = catalogue_selection(directory, (species,))
                request_id = f"catalogue-{species}"
                result, events = invoke(directory, payload, request_id)

                self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                self.assertEqual([event["kind"] for event in events],
                                 ["progress", "progress", "result"])
                self.assertTrue(all(event["requestId"] == request_id for event in events))
                terminal = events[-1]["payload"]
                self.assertEqual(terminal["standing"], "observed")
                self.assertEqual(terminal["studyRevisionId"], "catalogue-revision")
                self.assertEqual(terminal["artifacts"], [])
                self.assertEqual(terminal["observations"], {
                    "species": [{"speciesId": species,
                                 "chemistryId": payload["speciesRepresentations"][0]["chemistryId"],
                                 "coordinateAtomCount": atom_count,
                                 "parameterAtomCount": atom_count,
                                 "atomIdentityAndBondMatch": True, "warnings": []}],
                    "combinationWarnings": [], "combinedParameterizationObserved": True,
                })
                self.assertNotIn("qualification", terminal)

    def test_exact_pure_popc_is_observed_with_correlated_terminal_and_no_qualification(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            result, events = invoke(directory, payload)

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertEqual([event["kind"] for event in events],
                         ["progress", "progress", "result"])
        self.assertTrue(all(event["requestId"] == "membrane-slice-two" for event in events))
        terminal = events[-1]["payload"]
        self.assertEqual(terminal["studyRevisionId"], "revision-2")
        self.assertEqual(terminal["artifacts"], [])
        observed = terminal["observations"]
        self.assertEqual(observed["species"], [{
            "speciesId": "POPC", "chemistryId": "OpenMM-8.6.0-Lipid21-POPC",
            "coordinateAtomCount": 134, "parameterAtomCount": 134,
            "atomIdentityAndBondMatch": True, "warnings": [],
        }])
        self.assertTrue(observed["combinedParameterizationObserved"])
        self.assertEqual(observed["combinationWarnings"], [])
        self.assertNotIn("qualification", terminal)

    def test_unknown_construction_footprint_does_not_block_membrane_local_template_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            representation = payload["speciesRepresentations"][0]
            representation["areaPerMoleculeAngstromSquared"] = 0.0
            representation["volumeAngstromCubed"] = 0.0
            result, events = invoke(directory, payload, "membrane-local-without-footprint")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        observed = events[-1]["payload"]["observations"]
        self.assertTrue(observed["species"][0]["atomIdentityAndBondMatch"])
        self.assertTrue(observed["combinedParameterizationObserved"])
        self.assertEqual(observed["combinationWarnings"], [])

    def test_mixed_species_with_distinct_leaflet_inputs_are_parameterized_and_checked_separately(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC", "POPE"))
            payload["upper"]["fractions"] = [
                {"speciesId": "POPC", "fraction": 0.25},
                {"speciesId": "POPE", "fraction": 0.75},
            ]
            payload["lower"]["fractions"] = [
                {"speciesId": "POPC", "fraction": 0.60},
                {"speciesId": "POPE", "fraction": 0.40},
            ]
            result, events = invoke(directory, payload, "asymmetric-mixture")

            payload["lower"]["fractions"][1]["fraction"] = 0.30
            incoherent, incoherent_events = invoke(directory, payload, "incoherent-lower-leaflet")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertTrue(all(event["requestId"] == "asymmetric-mixture" for event in events))
        observed = events[-1]["payload"]["observations"]
        self.assertEqual([(item["speciesId"], item["coordinateAtomCount"],
                           item["parameterAtomCount"]) for item in observed["species"]],
                         [("POPC", 134, 134), ("POPE", 125, 125)])
        self.assertTrue(all(item["atomIdentityAndBondMatch"] and not item["warnings"]
                            for item in observed["species"]))
        self.assertTrue(observed["combinedParameterizationObserved"])
        self.assertEqual(observed["combinationWarnings"], [])
        self.assertEqual(incoherent.returncode, 0, incoherent.stderr + incoherent.stdout)
        self.assertTrue(any("lower leaflet fractions" in warning
                            for warning in incoherent_events[-1]["payload"]["observations"]["combinationWarnings"]))

    def test_qualified_catalogue_popc_cholesterol_mixture_has_exact_local_representation(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = catalogue_mixture(directory)
            result, events = invoke(directory, payload, "popc-cholesterol-mixture")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertEqual([event["kind"] for event in events], ["progress", "progress", "result"])
        self.assertTrue(all(event["requestId"] == "popc-cholesterol-mixture" for event in events))
        terminal = events[-1]["payload"]
        self.assertEqual(terminal["studyRevisionId"], "cholesterol-revision")
        observed = terminal["observations"]
        self.assertEqual([(item["speciesId"], item["coordinateAtomCount"],
                           item["parameterAtomCount"]) for item in observed["species"]],
                         [("POPC", 134, 134), ("CHL1", 74, 74)])
        self.assertTrue(all(item["atomIdentityAndBondMatch"] and not item["warnings"]
                            for item in observed["species"]))
        self.assertTrue(observed["combinedParameterizationObserved"])
        self.assertEqual(observed["combinationWarnings"], [])
        self.assertNotIn("qualification", terminal)

    def test_digest_mismatch_is_a_correlated_error_without_observations(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            payload["speciesRepresentations"][0]["coordinateTemplateSha256"] = "0" * 64
            result, events = invoke(directory, payload, "tampered-coordinate-hash")

        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([event["kind"] for event in events], ["progress", "error"])
        self.assertTrue(all(event["requestId"] == "tampered-coordinate-hash" for event in events))
        self.assertEqual(events[-1]["payload"]["failureCode"], "inputMismatch")
        self.assertNotIn("observations", events[-1]["payload"])

    def test_tampered_force_field_and_missing_catalogue_assets_are_refused(self):
        for damage, expected_code in (("tampered-force-field", "inputMismatch"),
                                      ("missing-force-field", "missingInput"),
                                      ("missing-coordinate", "missingInput")):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                payload = catalogue_selection(directory, ("POPC",))
                force_field = Path(payload["forceFieldFiles"][0]["path"])
                coordinate = Path(payload["speciesRepresentations"][0]["coordinateTemplatePath"])
                if damage == "tampered-force-field":
                    with force_field.open("ab") as stream:
                        stream.write(b"\n<!-- altered after catalogue identification -->\n")
                elif damage == "missing-force-field":
                    force_field.unlink()
                else:
                    coordinate.unlink()
                result, events = invoke(directory, payload, damage)

                self.assertNotEqual(result.returncode, 0)
                self.assertEqual([event["kind"] for event in events], ["progress", "error"])
                self.assertTrue(all(event["requestId"] == damage for event in events))
                self.assertEqual(events[-1]["payload"]["standing"], "failed")
                self.assertEqual(events[-1]["payload"]["failureCode"], expected_code)
                self.assertNotIn("observations", events[-1]["payload"])

    def test_declared_net_charge_must_match_parameterized_species(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = catalogue_selection(directory, ("POPC",))
            payload["speciesRepresentations"][0]["netChargeElementary"] = 1.0
            result, events = invoke(directory, payload, "wrong-net-charge")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertEqual([event["kind"] for event in events],
                         ["progress", "progress", "result"])
        self.assertTrue(all(event["requestId"] == "wrong-net-charge" for event in events))
        observed = events[-1]["payload"]["observations"]
        self.assertFalse(observed["species"][0]["atomIdentityAndBondMatch"])
        self.assertTrue(any("net charge" in warning.lower()
                            for warning in observed["species"][0]["warnings"]))
        self.assertFalse(observed["combinedParameterizationObserved"])

    def test_exact_hash_but_missing_atom_cannot_match_complete_parameter_template(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            representation = payload["speciesRepresentations"][0]
            coordinate = Path(representation["coordinateTemplatePath"])
            text = coordinate.read_text(encoding="ascii")
            coordinate.write_text("\n".join(line for line in text.splitlines()
                                            if not (line.startswith("HETATM") and
                                                    line[12:16].strip() == "H16Z")) + "\n",
                                  encoding="ascii")
            representation["coordinateTemplateSha256"] = digest(coordinate)
            result, events = invoke(directory, payload, "altered-atom-graph")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        observed = events[-1]["payload"]["observations"]
        self.assertFalse(observed["species"][0]["atomIdentityAndBondMatch"])
        self.assertTrue(observed["species"][0]["warnings"])
        self.assertFalse(observed["combinedParameterizationObserved"])
        self.assertTrue(observed["combinationWarnings"])

    def test_mirrored_lipid_coordinates_do_not_establish_the_named_stereochemistry(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            representation = payload["speciesRepresentations"][0]
            coordinate = Path(representation["coordinateTemplatePath"])
            lines = []
            for line in coordinate.read_text(encoding="ascii").splitlines():
                if line.startswith("HETATM"):
                    line = line[:30] + f"{-float(line[30:38]):8.3f}" + line[38:]
                lines.append(line)
            coordinate.write_text("\n".join(lines) + "\n", encoding="ascii")
            representation["coordinateTemplateSha256"] = digest(coordinate)
            result, events = invoke(directory, payload, "mirrored-configuration")

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        observed = events[-1]["payload"]["observations"]
        self.assertFalse(observed["species"][0]["atomIdentityAndBondMatch"],
                         "An exact atom graph with reflected handedness is a different molecular configuration.")
        self.assertTrue(observed["species"][0]["warnings"])
        self.assertFalse(observed["combinedParameterizationObserved"])

    def test_missing_or_unknown_stereo_descriptor_cannot_observe_exact_species(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = assets(directory, ("POPC",))
            representation = payload["speciesRepresentations"][0]
            representation["stereoChecks"] = []
            missing, missing_events = invoke(directory, payload, "missing-stereo-check")
            representation["stereoChecks"] = [{
                "kind": "unrecognized", "atomNames": ["O21", "C1", "C3", "HS"],
                "expected": "negative",
            }]
            unknown, unknown_events = invoke(directory, payload, "unknown-stereo-check")

        for result, events in ((missing, missing_events), (unknown, unknown_events)):
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            observed = events[-1]["payload"]["observations"]
            self.assertFalse(observed["species"][0]["atomIdentityAndBondMatch"])
            self.assertTrue(observed["species"][0]["warnings"])
            self.assertFalse(observed["combinedParameterizationObserved"])


if __name__ == "__main__":
    unittest.main()
