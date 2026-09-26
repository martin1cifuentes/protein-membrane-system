"""Red stage-observer proofs for exact optional-stage molecular readback.

The five-atom system exercises the real OpenMM stage observer. It makes no
scientific qualification claim. Deliberately rehashed CIF/sidecar mutations
must prevent the worker from reporting a completed-stage observation.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from openmm import Context, NonbondedForce, System, Vec3, VerletIntegrator, XmlSerializer, unit
from openmm.app import PDBxFile, Topology, element

from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import observe_stage
from ProteinInMembraneSystem.worker.exchange import WorkError
from ProteinInMembraneSystem.worker.parameterized_structure import _topology_data


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ExactEquilibratedStageReadback(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="pim-equilibration-readback-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        topology = Topology()
        chain = topology.addChain("A")
        atoms = (("ALA", "CA", element.carbon, "protein", "backbone", None),
                 ("DMP", "P", element.phosphorus, "lipid", "head", "upper"),
                 ("DMP", "P", element.phosphorus, "lipid", "head", "lower"),
                 ("HOH", "O", element.oxygen, "water", "waterAtom", None),
                 ("NA", "NA", element.sodium, "ion", "ionAtom", None))
        correspondence = []
        for index, (residue_name, atom_name, chemical_element, role, atom_role, side) in enumerate(atoms):
            residue = topology.addResidue(residue_name, chain, id=str(index + 1))
            topology.addAtom(atom_name, chemical_element, residue, id=str(index + 1))
            correspondence.append({
                "resultAtomIndex": index, "moleculeRole": role, "atomRole": atom_role,
                "physicalSide": side, "element": chemical_element.symbol,
                "sourceResidue": {"model": 0, "chain": "A", "residue": index + 1,
                                  "insertionCode": "", "copyId": "A"},
            })
        cell = tuple(Vec3(*(3.0 if component == axis else 0.0 for component in range(3)))
                     for axis in range(3))
        topology.setPeriodicBoxVectors(cell * unit.nanometer)
        self.positions = [Vec3(*point) for point in
                          ((0.3, 0.3, 1.5), (0.3, 0.3, 2.0), (0.3, 0.3, 1.0),
                           (1.1, 1.1, 1.5), (2.0, 2.0, 2.0))] * unit.nanometer
        system = System()
        nonbonded = NonbondedForce()
        nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(0.8 * unit.nanometer)
        for _, _, chemical_element, *_ in atoms:
            system.addParticle(chemical_element.mass)
            nonbonded.addParticle(0 * unit.elementary_charge, 0.3 * unit.nanometer,
                                  0.1 * unit.kilojoule_per_mole)
        system.addForce(nonbonded)
        system.setDefaultPeriodicBoxVectors(*[vector * unit.nanometer for vector in cell])
        integrator = VerletIntegrator(0.001 * unit.picoseconds)
        context = Context(system, integrator)
        context.setPositions(self.positions)
        context.setVelocitiesToTemperature(100 * unit.kelvin, 37)
        state = context.getState(getPositions=True, getVelocities=True,
                                 getForces=True, getEnergy=True)
        del context, integrator

        self.topology = topology
        self.cif = self.work / "equilibrated.cif"
        self.sidecar = self.work / "equilibrated-topology.json"
        self.system_xml = self.work / "equilibrated-system.xml"
        self.state_xml = self.work / "equilibrated-state.xml"
        self.correspondence_json = self.work / "correspondence.json"
        self._write_cif(self.cif, self.positions)
        self.sidecar.write_text(json.dumps(_topology_data(topology)), encoding="utf-8")
        self.system_xml.write_text(XmlSerializer.serialize(system), encoding="utf-8")
        self.state_xml.write_text(XmlSerializer.serialize(state), encoding="utf-8")
        self.correspondence_json.write_text(json.dumps({"complete": True, "atoms": correspondence}),
                                            encoding="utf-8")
        radii = {"C": 1.7, "P": 1.8, "O": 1.52, "Na": 2.27}
        self.payload = {
            "topologyCifPath": str(self.cif), "topologyCifSha256": digest(self.cif),
            "topologyJsonPath": str(self.sidecar), "topologyJsonSha256": digest(self.sidecar),
            "systemXmlPath": str(self.system_xml), "systemXmlSha256": digest(self.system_xml),
            "stateXmlPath": str(self.state_xml), "stateXmlSha256": digest(self.state_xml),
            "correspondencePath": str(self.correspondence_json),
            "correspondenceSha256": digest(self.correspondence_json),
            "stageKind": "equilibration",
            "localObservationSpec": {
                "contactRolePairs": [{"firstMoleculeRole": "protein", "secondMoleculeRole": "lipid"}],
                "atomRadiusByElementAngstrom": radii,
                "contactSearchRadiusAngstrom": 6.0, "maximumReportedPairs": 16,
                "usePeriodicBoundary": True,
                "requiredMetricNames": ["minimumIntermolecularDistanceAngstrom",
                                        "minimumIntermolecularHeavyAtomDistanceAngstrom",
                                        "upperLipidHeadMeanZAngstrom", "lowerLipidHeadMeanZAngstrom",
                                        "leafletHeadSeparationAngstrom",
                                        "proteinBilayerMidplaneOffsetAngstrom"],
            },
            "stageProteinGeometrySpec": {
                "requiredKinds": ["covalentBond", "chainContinuity", "nonbondedDistance"],
                "atomRadiusByElementAngstrom": radii,
                "neighborSearchRadiusAngstrom": 4.0, "excludedBondHops": 3,
                "maximumReportedPairs": 16,
            },
        }

    def _write_cif(self, path: Path, positions: object) -> None:
        with path.open("w", encoding="utf-8") as stream:
            PDBxFile.writeFile(self.topology, positions, stream, keepIds=True)

    def test_rehashed_coordinate_mutation_is_refused_before_stage_observation(self) -> None:
        baseline = observe_stage(self.work, self.payload, lambda *_: None)
        self.assertTrue(baseline["observations"]["atomOrderMatched"])
        self.assertTrue(baseline["observations"]["bondsMatched"])

        changed = self.work / "mutated-equilibrated.cif"
        positions = list(self.positions.value_in_unit(unit.nanometer))
        original = positions[0]
        positions[0] = Vec3(original.x + 0.002, original.y, original.z)
        self._write_cif(changed, positions * unit.nanometer)
        with self.assertRaises(WorkError) as mismatch:
            observe_stage(self.work,
                          dict(self.payload, topologyCifPath=str(changed),
                               topologyCifSha256=digest(changed)), lambda *_: None)
        self.assertEqual("inputMismatch", mismatch.exception.code)

    def test_rehashed_sidecar_cell_mutation_is_refused_before_stage_observation(self) -> None:
        baseline = observe_stage(self.work, self.payload, lambda *_: None)
        self.assertTrue(baseline["observations"]["atomOrderMatched"])
        self.assertTrue(baseline["observations"]["bondsMatched"])

        changed = self.work / "mutated-equilibrated-topology.json"
        data = json.loads(self.sidecar.read_text(encoding="utf-8"))
        data["boxVectorsAngstrom"][0][0] += 0.02
        changed.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaises(WorkError) as mismatch:
            observe_stage(self.work,
                          dict(self.payload, topologyJsonPath=str(changed),
                               topologyJsonSha256=digest(changed)), lambda *_: None)
        self.assertEqual("inputMismatch", mismatch.exception.code)

    def test_rehashed_sidecar_bond_mutation_is_refused_before_stage_observation(self) -> None:
        changed = self.work / "mutated-equilibrated-bonds.json"
        data = json.loads(self.sidecar.read_text(encoding="utf-8"))
        data["bonds"].append({"atomIndices": [0, 1], "order": 1})
        changed.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaises(WorkError) as mismatch:
            observe_stage(self.work,
                          dict(self.payload, topologyJsonPath=str(changed),
                               topologyJsonSha256=digest(changed)), lambda *_: None)
        self.assertEqual("inputMismatch", mismatch.exception.code)


if __name__ == "__main__":
    unittest.main()
