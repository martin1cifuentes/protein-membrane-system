"""Real preparation and placement crossings for selected ordinary water and ions.

The resulting prepared/oriented artifacts can seed separate full Memgen probes.
This test does not claim that a full provider construction has run.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

from openmm import NonbondedForce, unit
from openmm.app import ForceField, PDBFile


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/ProteinInMembraneSystem.IntegrationTests/ProteinPreparation"))
from test_worker_exchange import atom_line, invoke, two_alanines  # noqa: E402


CATALOGUE = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
REVISION = "retained-ordinary-members-ala2"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retained_source(extra_sodium: int = 0, *, unsupported_partner: bool = False) -> str:
    """Ala2 with selected complete TIP3P water and ion source chemistry.

    Additional Na residues occupy a separated 5 x 2 x 2 aqueous grid.  They
    are input molecules, never an adjustment to Memgen's charge argument.
    """
    if extra_sodium not in (0, 4, 20):
        raise ValueError("Only the named base and charge-branch fixtures are defined")
    members = [
        ("O", "HOH", 50, (12.000, 7.000, 27.000), "O"),
        ("H1", "HOH", 50, (12.957, 7.000, 27.000), "H"),
        ("H2", "HOH", 50, (11.760, 7.927, 27.000), "H"),
        ("NA", "NA", 51, (11.000, -6.000, 27.000), "NA"),
        ("CL", "CL", 52, (-6.000, 6.000, -27.000), "CL"),
    ]
    sodium_grid = [(x, y, z) for z in (-27.0, 27.0)
                   for y in (-12.0, 12.0) for x in (-12.0, -6.0, 0.0, 6.0, 12.0)]
    members.extend(("NA", "NA", 53 + index, point, "NA")
                   for index, point in enumerate(sodium_grid[:extra_sodium]))
    if unsupported_partner:
        members.append(("C1", "GOL", 80, (18.0, 0.0, 0.0), "C"))
    partner_lines = "".join(
        atom_line(index, name, residue, "B", number, xyz, element).replace(
            "ATOM  ", "HETATM", 1)
        for index, (name, residue, number, xyz, element) in enumerate(members, 12))
    return two_alanines(1, "A").replace("ENDMDL\n", partner_lines + "TER\nENDMDL\n") + "END\n"


def staged_preparation_payload(directory: Path, source: Path) -> tuple[dict, list[Path]]:
    chemical = CATALOGUE["proteinChemicalStates"][0]
    assets = []
    staged = []
    for item in chemical["forceFieldFiles"]:
        installed = (ROOT / "config/policies" / item["path"]).resolve()
        if not installed.is_file() or digest(installed) != item["sha256"]:
            raise AssertionError(f"Pinned preparation asset is unavailable: {item['id']}")
        copied = directory / installed.name
        shutil.copyfile(installed, copied)
        assets.append(dict(item, path=str(copied)))
        staged.append(copied)
    return ({
        "studyRevisionId": REVISION,
        "sourcePath": str(source), "sourceSha256": digest(source),
        "modelIndex": 0, "assemblyId": None,
        "chainSelections": [{"sourceChain": "A", "copyId": "A"}],
        "altlocChoices": [], "retainedPartners": [], "approvedDisulfides": [],
        "overrides": [], "nominalPh": 7.0,
        "seed": chemical["recommendation"]["seed"],
        "forceFieldFiles": assets, "permittedVariants": chemical["permittedVariants"],
        "disulfideCandidateMaxSgDistanceAngstrom":
            chemical["disulfideCandidateMaxSgDistanceAngstrom"],
        "geometrySpec": CATALOGUE["proteinStructuralPolicies"][0]["measurement"],
    }, staged)


def checked_result(test: unittest.TestCase, directory: Path, operation: str,
                   payload: dict) -> dict:
    process, events = invoke(directory, operation, payload, f"retained-{operation}")
    test.assertEqual(process.returncode, 0, process.stderr + process.stdout)
    test.assertTrue(events, process.stderr)
    test.assertEqual(events[-1]["kind"], "result", events[-1])
    result = events[-1]["payload"]
    for artifact in result.get("artifacts", []):
        test.assertEqual(digest(Path(artifact["path"])), artifact["sha256"])
    return result


def artifact_path(result: dict, role: str) -> Path:
    matches = [Path(item["path"]) for item in result["artifacts"] if item["role"] == role]
    if len(matches) != 1:
        raise AssertionError(f"Expected one {role} artifact, observed {len(matches)}")
    return matches[0]


def atom_identity(topology) -> list[tuple[str, str, str, str, str]]:
    return [(atom.residue.chain.id, atom.residue.id, atom.residue.name,
             atom.name, atom.element.symbol) for atom in topology.atoms()]


def pdb_coordinates(path: Path) -> list[tuple[float, float, float]]:
    return [(float(line[30:38]), float(line[38:46]), float(line[46:54]))
            for line in path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))]


def bond_indices(topology) -> set[tuple[int, int]]:
    return {tuple(sorted((first.index, second.index))) for first, second in topology.bonds()}


def charges(topology, force_field_paths: list[Path]) -> list[float]:
    system = ForceField(*(str(path) for path in force_field_paths)).createSystem(topology)
    nonbonded = [system.getForce(index) for index in range(system.getNumForces())
                 if isinstance(system.getForce(index), NonbondedForce)]
    if len(nonbonded) != 1:
        raise AssertionError("Expected exactly one parameterized nonbonded force")
    return [nonbonded[0].getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
            for index in range(system.getNumParticles())]


class RetainedMemberCrossingTests(unittest.TestCase):
    def test_real_preparation_and_manual_placement_keep_exact_water_and_ions(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "ala2-retained-source.pdb"
            source.write_text(retained_source(), encoding="ascii")
            source_sha = digest(source)
            payload, ff_paths = staged_preparation_payload(work, source)

            inspected = checked_result(self, work, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": source_sha,
                "maxAtoms": 1000, "sourceKind": "upload", "prediction": None,
            })
            partners = inspected["observations"]["models"][0]["partners"]
            self.assertEqual({(item["sourceId"], item["kind"]) for item in partners}, {
                ("0:B:50::HOH", "water"), ("0:B:51::NA", "ion"),
                ("0:B:52::CL", "ion"),
            })
            payload["retainedPartners"] = partners
            preview = checked_result(self, work, "inspect_preparation_changes", payload)
            proposed = preview["observations"]["missingNonbackboneHeavyAtoms"]
            self.assertEqual(proposed, [])
            payload["proposedHeavyAtoms"] = proposed
            proposed_result = checked_result(self, work, "recommend_preparation", payload)
            plan = proposed_result["observations"]
            self.assertEqual(plan["proposedHeavyAtoms"], [])
            self.assertEqual(digest(source), source_sha)

            authorized = dict(payload,
                              approvedHeavyAtoms=[], normalizeProteinHydrogens=True,
                              planSeed=payload["seed"],
                              residueVariants=[{"residue": item["residue"],
                                                "variant": item["variant"],
                                                "decisionId": "accepted-fixture-plan"}
                                               for item in plan["choices"]])
            prepared_result = checked_result(self, work, "prepare_protein", authorized)
            prepared = artifact_path(prepared_result, "preparedPdb")
            graph = artifact_path(prepared_result, "preparedBondGraph")
            correspondence = artifact_path(prepared_result, "correspondenceJson")
            self.assertEqual(digest(prepared), plan["candidateSha256"])
            prepared_topology = PDBFile(str(prepared)).topology
            prepared_identity = atom_identity(prepared_topology)
            self.assertEqual(prepared_result["observations"]["preparedAtomCount"],
                             len(prepared_identity))

            mapping = json.loads(correspondence.read_text())
            self.assertTrue(mapping["complete"])
            self.assertEqual((mapping["sourceId"], mapping["resultId"]),
                             (source_sha, digest(prepared)))
            self.assertEqual(len(mapping["atoms"]), len(prepared_identity))
            retained = [item for item in mapping["atoms"]
                        if item["moleculeRole"] in {"water", "ion"}]
            self.assertEqual({(item["sourceResidue"]["chain"],
                               item["sourceResidue"]["residue"], item["moleculeRole"])
                              for item in retained},
                             {("B", 50, "water"), ("B", 51, "ion"), ("B", 52, "ion")})
            self.assertEqual(len(retained), 5)
            self.assertTrue(all(item["role"] == "source" and item["sourceAtomId"]
                                for item in retained))

            sidecar = json.loads(graph.read_text())
            graph_bonds = {tuple(bond["atomIndices"]) for bond in sidecar["bonds"]}
            self.assertEqual(graph_bonds, bond_indices(prepared_topology))
            water_atoms = {index: atom.name for index, atom in enumerate(prepared_topology.atoms())
                           if atom.residue.name == "HOH"}
            self.assertEqual(set(water_atoms.values()), {"O", "H1", "H2"})
            water_bonds = {tuple(sorted((water_atoms[first], water_atoms[second])))
                           for first, second in graph_bonds
                           if first in water_atoms and second in water_atoms}
            self.assertEqual(water_bonds, {("H1", "O"), ("H2", "O")})
            self.assertEqual(sum(first in water_atoms or second in water_atoms
                                 for first, second in graph_bonds), 2)
            ion_atoms = {atom.index: atom for atom in prepared_topology.atoms()
                         if atom.residue.name in {"NA", "CL"}}
            self.assertEqual({(atom.residue.name, atom.name) for atom in ion_atoms.values()},
                             {("NA", "NA"), ("CL", "CL")})
            self.assertFalse(any(first in ion_atoms or second in ion_atoms
                                 for first, second in graph_bonds))
            prepared_charges = charges(prepared_topology, ff_paths)
            self.assertAlmostEqual(prepared_charges[next(index for index, atom in ion_atoms.items()
                                                         if atom.residue.name == "NA")], 1.0, places=6)
            self.assertAlmostEqual(prepared_charges[next(index for index, atom in ion_atoms.items()
                                                         if atom.residue.name == "CL")], -1.0, places=6)
            self.assertAlmostEqual(sum(prepared_charges), 0.0, places=5)

            placed = checked_result(self, work, "place_manual", {
                "preparedPdbPath": str(prepared), "preparedSha256": digest(prepared),
                "preparedProteinId": digest(prepared),
                "preparedAtomCount": len(prepared_identity), "maximumAtomCount": 120000,
                "startingPosition": "center",
                "offsetXAngstrom": 0.0, "offsetYAngstrom": 0.0,
                "offsetZAngstrom": 0.0,
                "rotationXDegrees": 0.0, "rotationYDegrees": 0.0,
                "rotationZDegrees": 0.0, "leafletEnvelopeAngstrom": 23.0,
            })
            oriented = artifact_path(placed, "orientedPdb")
            oriented_topology = PDBFile(str(oriented)).topology
            self.assertEqual(atom_identity(oriented_topology), prepared_identity)
            self.assertEqual(bond_indices(oriented_topology), graph_bonds)
            self.assertEqual(placed["observations"]["orientedAtomCount"],
                             len(prepared_identity))
            for before, after in zip(prepared_charges, charges(oriented_topology, ff_paths)):
                self.assertAlmostEqual(after, before, places=7)
            prepared_xyz, oriented_xyz = pdb_coordinates(prepared), pdb_coordinates(oriented)
            shift = tuple(oriented_xyz[0][axis] - prepared_xyz[0][axis] for axis in range(3))
            self.assertEqual(len(prepared_xyz), len(oriented_xyz))
            for before, after in zip(prepared_xyz, oriented_xyz):
                for axis in range(3):
                    self.assertAlmostEqual(after[axis] - before[axis], shift[axis], places=3)
            self.assertEqual(digest(source), source_sha)

    def test_unsupported_selected_partner_is_refused_before_preparation(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "ala2-unsupported-partner.pdb"
            source.write_text(retained_source(unsupported_partner=True), encoding="ascii")
            payload, _ = staged_preparation_payload(work, source)
            inspected = checked_result(self, work, "inspect_source", {
                "sourcePath": str(source), "sourceSha256": digest(source),
                "maxAtoms": 1000, "sourceKind": "upload", "prediction": None,
            })
            payload["retainedPartners"] = inspected["observations"]["models"][0]["partners"]
            process, events = invoke(work, "inspect_preparation_changes", payload,
                                     "retained-unsupported")
            self.assertNotEqual(process.returncode, 0)
            self.assertEqual(events[-1]["kind"], "error")
            self.assertEqual(events[-1]["payload"]["failureCode"], "unsupportedChemistry")

    def test_charge_branch_sources_add_real_distinct_sodium_residues(self):
        for extra in (0, 4, 20):
            with self.subTest(extra_sodium=extra):
                lines = retained_source(extra).splitlines()
                sodium = [line for line in lines if line.startswith("HETATM")
                          and line[17:20].strip() == "NA"]
                self.assertEqual(len(sodium), 1 + extra)
                self.assertEqual(len({(line[21], line[22:26]) for line in sodium}), 1 + extra)
                self.assertEqual(sum(line.startswith("HETATM") for line in lines), 5 + extra)


if __name__ == "__main__":
    unittest.main()
