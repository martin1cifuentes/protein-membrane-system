"""Small artifact-level proof for the read-only final-State inspector."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from openmm import (Context, CustomExternalForce, HarmonicBondForce, Platform,
                    System, Vec3, VerletIntegrator, XmlSerializer, unit)
from openmm.app import PDBxFile, Topology, element

from inspect_final_state import FILES, InspectionError, _contacts, inspect, main


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SavedFinalStateInspectorCases(unittest.TestCase):
    def test_hydrogen_only_short_contact_is_reported_without_heavy_overlap(self):
        # Protein C-H and lipid C-H are distinct bonded molecules. Their H-H
        # minimum is 0.2 A, while their only heavy-heavy pair is 3.0 A.
        xyz = np.asarray(((0, 0, 0), (1, 0, 0), (3, 0, 0), (1.2, 0, 0)), dtype=float)
        observed = _contacts(xyz, np.asarray((40, 40, 40), dtype=float),
                             np.asarray(("protein", "protein", "lipid", "lipid")),
                             np.asarray(("C", "H", "C", "H")),
                             np.asarray((0, 0, 2, 2)), 6.0)
        all_pair = observed["allAtom"]["protein|lipid"]
        whole = observed["allAtom"]["wholeSystem"]
        heavy = observed["heavy"]["protein|lipid"]
        self.assertAlmostEqual(0.2, all_pair["nearestAllAtomAngstrom"])
        self.assertAlmostEqual(0.2, whole["nearestAllAtomAngstrom"])
        self.assertGreater(all_pair["allAtomPairsBelow1_5Angstrom"], 0)
        self.assertAlmostEqual(3.0, heavy["nearestHeavyAngstrom"])
        self.assertEqual(0, heavy["heavyPairsBelow1_5Angstrom"])

    def test_synthetic_saved_artifacts_project_force_and_track_contact_classes(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            topology = Topology()
            topology.setUnitCellDimensions(Vec3(40, 40, 80) * unit.angstrom)
            atoms = []
            labels = [
                ("C", "protein", "backbone", None, "PRO", "A"),
                ("H", "protein", "body", None, "PRO", "A"),
                ("P", "lipid", "head", "upper", "DMP", "B"),
                ("C", "lipid", "tail", "upper", "DMP", "B"),
                ("P", "lipid", "head", "lower", "DMP", "C"),
                ("O", "water", "body", None, "HOH", "D"),
                ("O", "water", "body", None, "HOH", "E"),
                ("Na", "ion", "body", None, "NA", "F"),
                ("Na", "ion", "body", None, "NA", "G"),
            ]
            chains = {}
            residues = {}
            for index, (symbol, _, _, _, residue_name, chain_id) in enumerate(labels):
                if chain_id not in chains:
                    chains[chain_id] = topology.addChain(chain_id)
                    residues[chain_id] = topology.addResidue(residue_name, chains[chain_id], id="1")
                atoms.append(topology.addAtom(symbol, element.get_by_symbol(symbol),
                                              residues[chain_id], id=str(index + 1)))
            topology.addBond(atoms[0], atoms[1])
            topology.addBond(atoms[2], atoms[3])
            first = [Vec3(*point) for point in (
                (0.1, 0, 4), (1.1, 0, 4),  # constrained 1 nm protein bond
                (0.1, 0, 5), (0.3, 0, 4), (0.1, 0, 3),
                (0.1, 0, 4), (3.99, 0, 4), (3.95, 0, 4), (0.03, 0, 4))]
            last = list(first)
            last[6] = Vec3(3.0, 0, 4)  # remove water-water seam pair
            last[8] = Vec3(2.0, 0, 4)  # remove ion-ion seam pair

            system = System()
            for _ in labels:
                system.addParticle(12.0)
            system.setDefaultPeriodicBoxVectors(*(Vec3(*point) * unit.nanometer
                                                  for point in ((4, 0, 0), (0, 4, 0), (0, 0, 8))))
            system.addConstraint(0, 1, 1.0 * unit.nanometer)
            normal = HarmonicBondForce()
            normal.addBond(0, 1, 0.8 * unit.nanometer,
                           1000 * unit.kilojoule_per_mole / unit.nanometer**2)
            system.addForce(normal)
            transverse = CustomExternalForce("10*y")
            transverse.addParticle(0, [])
            system.addForce(transverse)
            integrator = VerletIntegrator(0.001 * unit.picoseconds)
            context = Context(system, integrator, Platform.getPlatformByName("Reference"))
            context.setPositions(first * unit.nanometer)
            initial = context.getState(getPositions=True, getForces=True, getEnergy=True)
            context.setPositions(last * unit.nanometer)
            final = context.getState(getPositions=True, getForces=True, getEnergy=True)
            del context, integrator

            with (work / FILES["topologyCif"]).open("w") as stream:
                PDBxFile.writeFile(topology, initial.getPositions(), stream, keepIds=True)
            with (work / FILES["minimizedCif"]).open("w") as stream:
                PDBxFile.writeFile(topology, final.getPositions(), stream, keepIds=True)
            (work / FILES["systemXml"]).write_text(XmlSerializer.serialize(system))
            (work / FILES["stateXml"]).write_text(XmlSerializer.serialize(initial))
            (work / FILES["minimizedStateXml"]).write_text(XmlSerializer.serialize(final))
            sidecar = {
                "formatVersion": 1,
                "chains": [{"id": chain.id} for chain in topology.chains()],
                "residues": [{"chainIndex": residue.chain.index, "id": residue.id,
                              "name": residue.name, "insertionCode": residue.insertionCode}
                             for residue in topology.residues()],
                "atoms": [{"residueIndex": atom.residue.index, "id": atom.id,
                           "name": atom.name, "element": atom.element.symbol}
                          for atom in topology.atoms()],
                "bonds": [{"atomIndices": [0, 1], "order": None},
                          {"atomIndices": [2, 3], "order": None}],
                "boxVectorsAngstrom": [[40.0, 0, 0], [0, 40.0, 0], [0, 0, 80.0]],
            }
            (work / FILES["topologyJson"]).write_text(json.dumps(sidecar))
            correspondence = {
                "complete": True, "sourceId": "a" * 64,
                "resultId": digest(work / FILES["topologyCif"]),
                "atoms": [{"resultAtomIndex": index, "resultAtomId": f"fixture:{index}",
                           "moleculeRole": role, "atomRole": atom_role,
                           "physicalSide": side, "element": symbol,
                           "generatedSpeciesId": ("DMPC" if role == "lipid" else
                                                  "HOH" if role == "water" else
                                                  "NA" if role == "ion" else None),
                           "generatedComponentRole": ("lipid" if role == "lipid" else
                                                      "water" if role == "water" else
                                                      "positiveIon" if role == "ion" else None)}
                          for index, (symbol, role, atom_role, side, _, _) in enumerate(labels)],
            }
            (work / FILES["correspondenceJson"]).write_text(json.dumps(correspondence))
            expected = {"atomCount": len(labels), "sourceId": "a" * 64,
                        "sha256ByRole": {role: digest(work / filename)
                                         for role, filename in FILES.items()}}
            manifest = work / "expected.json"
            report_path = work / "inspection.json"
            manifest.write_text(json.dumps(expected))
            self.assertEqual(0, main(["--constructed-dir", str(work),
                                      "--expected", str(manifest),
                                      "--output", str(report_path)]))
            self.assertEqual("observed", json.loads(report_path.read_text())["status"])

            report = inspect(work, work, expected)
            self.assertEqual("observed", report["status"])
            self.assertEqual(9, report["atomCount"])
            self.assertEqual(2, report["population"]["proteinAtomCount"])
            self.assertEqual([{"speciesId": "DMPC", "physicalSide": "lower", "count": 1},
                              {"speciesId": "DMPC", "physicalSide": "upper", "count": 1}],
                             report["population"]["lipidCounts"])
            self.assertEqual((2, 2, 0),
                             (report["population"]["waterCount"],
                              report["population"]["sodiumCount"],
                              report["population"]["chlorideCount"]))
            self.assertEqual(1, report["finalForce"]["constraintCount"])
            self.assertAlmostEqual(math.sqrt(100 / 9),
                                   report["finalForce"]["tangentPerParticleRmsKjMolNm"], places=5)
            self.assertGreater(report["finalForce"]["rawComponentRmsKjMolNm"], 50)
            self.assertLess(report["finalForce"]["maximumRelativeConstraintError"], 1e-8)
            initial_contacts = report["stages"]["initial"]["heavyPeriodicContacts"]
            final_contacts = report["stages"]["final"]["heavyPeriodicContacts"]
            initial_all = report["stages"]["initial"]["allAtomPeriodicContacts"]
            final_all = report["stages"]["final"]["allAtomPeriodicContacts"]
            self.assertEqual(1, initial_contacts["water|water"]["periodicImageHeavyPairsBelow2_2Angstrom"])
            self.assertEqual(0, final_contacts["water|water"]["heavyPairsWithinRadius"])
            self.assertEqual(1, initial_contacts["ion|ion"]["periodicImageHeavyPairsBelow2_2Angstrom"])
            self.assertEqual(0, final_contacts["ion|ion"]["heavyPairsWithinRadius"])
            self.assertEqual(1, initial_all["water|water"]["periodicImageAllAtomPairsBelow2_2Angstrom"])
            self.assertEqual(0, final_all["water|water"]["allAtomPairsWithinRadius"])
            self.assertLess(report["allAtomContactChange"]["ion|ion"]["changeInAllAtomPairsBelow1_5Angstrom"], 0)
            self.assertIsNotNone(initial_all["wholeSystem"]["nearestAllAtomAngstrom"])
            self.assertGreater(initial_contacts["water|ion"]["heavyPairsWithinRadius"], 0)
            self.assertGreater(initial_contacts["protein|lipid"]["heavyPairsWithinRadius"], 0)
            for stage in ("initial", "final"):
                organization = report["stages"][stage]["relativeOrganization"]
                self.assertAlmostEqual(20, organization["leafletHeadSeparationAngstrom"])
                self.assertAlmostEqual(0, organization["proteinBilayerMidplaneOffsetAngstrom"])

            wrong = json.loads(json.dumps(expected))
            wrong["sha256ByRole"]["minimizedStateXml"] = "0" * 64
            with self.assertRaisesRegex(InspectionError, "SHA-256 mismatch for minimizedStateXml"):
                inspect(work, work, wrong)
            bad_correspondence = dict(correspondence, complete=False)
            (work / FILES["correspondenceJson"]).write_text(json.dumps(bad_correspondence))
            wrong["sha256ByRole"]["minimizedStateXml"] = expected["sha256ByRole"]["minimizedStateXml"]
            wrong["sha256ByRole"]["correspondenceJson"] = digest(work / FILES["correspondenceJson"])
            with self.assertRaisesRegex(InspectionError, "correspondence is incomplete"):
                inspect(work, work, wrong)
            wrong_side = json.loads(json.dumps(correspondence))
            wrong_side["atoms"][4]["physicalSide"] = "upper"
            (work / FILES["correspondenceJson"]).write_text(json.dumps(wrong_side))
            wrong["sha256ByRole"]["correspondenceJson"] = digest(work / FILES["correspondenceJson"])
            with self.assertRaisesRegex(InspectionError, "both DMPC leaflets"):
                inspect(work, work, wrong)


if __name__ == "__main__":
    unittest.main()
