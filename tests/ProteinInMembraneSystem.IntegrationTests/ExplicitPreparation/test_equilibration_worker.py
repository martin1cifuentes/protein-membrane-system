"""Focused OpenMM worker proofs for optional-stage measurements and artifacts.

These synthetic systems verify mechanics; they do not qualify an equilibration
protocol or a protein/lipid class for actor use.
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

import numpy as np
from openmm import (Context, MonteCarloMembraneBarostat, NonbondedForce, System,
                    Vec3, VerletIntegrator, XmlSerializer, unit)
from openmm.app import PDBxFile, Topology, element

from ProteinInMembraneSystem.ExplicitPreparation.OptionalEquilibrationProcedure.worker.equilibration import (
    _advance_temperature_control, _count_within_distance, _initial_temperature,
    _protein_restraint_indices, equilibrate,
)
from ProteinInMembraneSystem.CompletedStageExport.worker.export_verification import verify_export
from ProteinInMembraneSystem.worker.exchange import WorkError
from ProteinInMembraneSystem.worker.parameterized_structure import _topology_data


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class OptionalEquilibrationWorkerMechanics(unittest.TestCase):
    def test_declared_temperature_ramp_and_protein_restraint_targets_are_exact(self) -> None:
        class IntegratorRecorder:
            def __init__(self) -> None:
                self.events: list[tuple[str, float | int]] = []

            def setTemperature(self, temperature) -> None:
                self.events.append(("temperature", temperature.value_in_unit(unit.kelvin)))

            def step(self, steps: int) -> None:
                self.events.append(("steps", steps))

        ramp = {"steps": 3, "temperatureKelvin": 303.0,
                "initialTemperatureKelvin": 100.0, "pressureBar": None}
        self.assertEqual(100.0, _initial_temperature(ramp))
        integrator = IntegratorRecorder()
        _advance_temperature_control(integrator, 100.0, 303.0, 3, 0, 2)
        _advance_temperature_control(integrator, 100.0, 303.0, 3, 2, 1)
        self.assertEqual([100.0, 201.5, 303.0],
                         [value for name, value in integrator.events if name == "temperature"])
        self.assertEqual([1, 1, 1], [value for name, value in integrator.events if name == "steps"])
        fixed = IntegratorRecorder()
        _advance_temperature_control(fixed, None, 303.0, 3, 0, 3)
        self.assertEqual([("steps", 3)], fixed.events)
        with self.assertRaises(WorkError):
            _initial_temperature(dict(ramp, pressureBar=1.0))

        topology = Topology()
        chain = topology.addChain("A")
        residue = topology.addResidue("ALA", chain)
        for name, chemical_element in (("N", element.nitrogen), ("CA", element.carbon),
                                       ("C", element.carbon), ("H", element.hydrogen),
                                       ("CB", element.carbon)):
            topology.addAtom(name, chemical_element, residue)
        atoms = list(topology.atoms())
        backbone = [0, 1, 2, 3]
        heavy = [0, 1, 2, 4]
        self.assertEqual(backbone, _protein_restraint_indices("backbone", atoms, backbone, heavy))
        self.assertEqual(heavy, _protein_restraint_indices("protein-heavy", atoms, backbone, heavy))
        self.assertEqual([0, 1, 2], _protein_restraint_indices("backbone-heavy", atoms, backbone, heavy))
        self.assertEqual([1], _protein_restraint_indices("protein-ca", atoms, backbone, heavy))
        with self.assertRaises(WorkError):
            _protein_restraint_indices("protein-ca", atoms, [0, 2, 3], heavy)

    def test_periodic_contact_count_matches_full_scan_for_orthorhombic_and_tilted_cells(self) -> None:
        rng = np.random.default_rng(4831)
        for box in (np.diag([82.0, 63.0, 111.0]),
                    np.array([[82.0, 0.0, 0.0], [9.0, 63.0, 0.0],
                              [3.0, 7.0, 111.0]])):
            fractional = rng.random((320, 3))
            fractional[:4] = np.array([[0.998, 0.51, 0.51], [0.002, 0.51, 0.51],
                                       [0.44, 0.998, 0.51], [0.44, 0.002, 0.51]])
            first = list(range(0, 160))
            reference = list(range(160, 320)) + [1, 3]
            inverse = np.linalg.inv(box)
            for cutoff in (1.5, 4.0, 12.0, 120.0):
                observed = _count_within_distance(fractional, box, inverse, first, reference, cutoff)
                expected = 0
                for index in first:
                    delta = fractional[reference] - fractional[index]
                    delta -= np.rint(delta)
                    wrapped = delta @ box
                    expected += int(np.any(np.sum(wrapped * wrapped, axis=1) <= cutoff * cutoff))
                self.assertEqual(expected, observed, (box.tolist(), cutoff))

    def test_one_declared_window_writes_stage_cell_in_own_topology_sidecar(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pim-equilibration-worker-") as raw:
            work = Path(raw)
            topology = Topology()
            chain = topology.addChain("A")
            atom_specs = (("ALA", "CA", element.carbon),
                          ("DMP", "P", element.phosphorus),
                          ("DMP", "P", element.phosphorus),
                          ("HOH", "O", element.oxygen),
                          ("NA", "NA", element.sodium))
            for ordinal, (residue_name, atom_name, chemical_element) in enumerate(atom_specs, 1):
                residue = topology.addResidue(residue_name, chain, id=str(ordinal))
                topology.addAtom(atom_name, chemical_element, residue, id=str(ordinal))
            box = tuple(Vec3(value, 0, 0) if axis == 0 else
                        Vec3(0, value, 0) if axis == 1 else Vec3(0, 0, value)
                        for axis, value in enumerate((3.0, 3.0, 3.0)))
            topology.setPeriodicBoxVectors(box * unit.nanometer)
            positions = [Vec3(*point) for point in ((0.3, 0.3, 1.5), (1.0, 1.0, 2.0),
                                                    (1.0, 1.0, 1.0), (0.9, 0.3, 1.5),
                                                    (2.0, 2.0, 2.0))] * unit.nanometer
            system = System()
            nonbonded = NonbondedForce()
            nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
            nonbonded.setCutoffDistance(0.8 * unit.nanometer)
            for chemical_element in (item[2] for item in atom_specs):
                system.addParticle(chemical_element.mass)
                nonbonded.addParticle(0 * unit.elementary_charge, 0.3 * unit.nanometer,
                                      0.1 * unit.kilojoule_per_mole)
            system.addForce(nonbonded)
            system.setDefaultPeriodicBoxVectors(*[vector * unit.nanometer for vector in box])
            integrator = VerletIntegrator(0.001 * unit.picoseconds)
            context = Context(system, integrator)
            context.setPositions(positions)
            context.setVelocitiesToTemperature(100 * unit.kelvin, 31)
            state = context.getState(getPositions=True, getVelocities=True,
                                     getForces=True, getEnergy=True)
            del context, integrator

            cif = work / "minimized.cif"
            sidecar = work / "minimized-topology.json"
            system_xml = work / "minimized-system.xml"
            state_xml = work / "minimized-state.xml"
            with cif.open("w", encoding="utf-8") as stream:
                PDBxFile.writeFile(topology, positions, stream, keepIds=True)
            sidecar.write_text(json.dumps(_topology_data(topology)), encoding="utf-8")
            system_xml.write_text(XmlSerializer.serialize(system), encoding="utf-8")
            state_xml.write_text(XmlSerializer.serialize(state), encoding="utf-8")

            observables = [
                ("temperature", "K", "system", "temperature", "none", "none", "none", None, [], [], []),
                ("leafletSeparation", "angstrom", "membraneOrganization", "centroidSeparationZ",
                 "upper-lipid-head", "lower-lipid-head", "none", None, [1], [2], []),
                ("proteinOffset", "angstrom", "proteinPlacement", "proteinMidplaneOffset",
                 "protein-backbone", "upper-lipid-head", "lower-lipid-head", None, [0], [1], [2]),
                ("hydration", "count", "hydrationIons", "countWithinDistance",
                 "water-oxygen", "protein-or-lipid-heavy", "none", 5.0, [3], [0, 1, 2], []),
                ("ions", "count", "hydrationIons", "countWithinDistance",
                 "all-ions", "protein-or-lipid-heavy", "none", 5.0, [4], [0, 1, 2], []),
            ]
            declarations = []
            resolved = []
            for (name, measurement_unit, category, method, first, second, third, cutoff,
                 indices, comparison, third_indices) in observables:
                common = {"name": name, "unit": measurement_unit, "scope": "synthetic",
                          "category": category, "method": method,
                          "distanceCutoffAngstrom": cutoff}
                declarations.append(dict(common, atomSelector=first, comparisonSelector=second,
                                         thirdSelector=third))
                resolved.append(dict(common, atomIndices=indices,
                                     comparisonAtomIndices=comparison,
                                     thirdAtomIndices=third_indices))
            control = {"name": "unrestrained-observation", "steps": 3,
                       "reportIntervalSteps": 1, "timestepPicoseconds": 0.001,
                       "temperatureKelvin": 303.0, "frictionPerPicosecond": 1.0,
                       "pressureBar": 1000.0, "pressureMode": "xyisotropiczfree",
                       "barostatFrequencySteps": 1, "surfaceTensionBarNm": 0.0,
                       "proteinRestraintKjMolNm2": 0.0, "lipidRestraintKjMolNm2": 0.0}
            protocol = {"id": "synthetic-worker-mechanics", "randomSeed": 41,
                        "stages": [control], "extensionWindow": dict(control, name="extension"),
                        "maximumExtensions": 0, "maximumSampleCount": 3,
                        "maximumFrameBytes": 3 * 5 * 24,
                        "requiredObservations": [row["name"] for row in declarations],
                        "observables": declarations,
                        "sufficiencyRules": [
                            {"observableName": row["name"], "blockSizeSamples": 1,
                             "minimumEffectiveBlocks": 3,
                             "maximumAbsoluteFirstVsLastBlockMeanDifference": 1e9,
                             "maximumAbsoluteLagOneBlockCorrelation": 0.99}
                            for row in declarations],
                        "comparisonBasis": "synthetic-mechanics-only"}
            payload = {"topologyCifPath": str(cif), "topologyCifSha256": digest(cif),
                       "topologyJsonPath": str(sidecar), "topologyJsonSha256": digest(sidecar),
                       "systemXmlPath": str(system_xml), "systemXmlSha256": digest(system_xml),
                       "minimizedStateXmlPath": str(state_xml),
                       "minimizedStateXmlSha256": digest(state_xml),
                       "studyRevisionId": "revision-1", "attemptId": "attempt-1",
                       "stageId": "a" * 32,
                       "sourceMinimizedStageId": "minimized-1",
                       "protocolSha256": "a" * 64,
                       "validatedPolicyId": protocol["id"],
                       "proteinBackboneAtomIndices": [0], "proteinHeavyAtomIndices": [0],
                       "lipidHeavyAtomIndices": [1, 2],
                       "resolvedObservables": resolved, "protocol": protocol}
            too_small = dict(protocol, maximumFrameBytes=3 * 5 * 24 - 1)
            with self.assertRaises(WorkError) as budget:
                equilibrate(work, dict(payload, stageId="8" * 32, protocol=too_small),
                            lambda *_: None)
            self.assertEqual("resourceRefusal", budget.exception.code)
            self.assertFalse((work / ("equilibration-" + "8" * 32)).exists())
            progress = []
            result = equilibrate(work, payload, lambda kind, item: progress.append((kind, item)))
            artifacts = {item["role"]: Path(item["path"]) for item in result["artifacts"]}
            self.assertTrue(all(path.parent == work / ("equilibration-" + "a" * 32)
                                for path in artifacts.values()))
            self.assertEqual({"equilibratedCif", "equilibratedTopologyJson",
                              "equilibratedSystemXml", "equilibratedStateXml",
                              "sampledPositionsF64", "sampledFramesManifest"}, set(artifacts))
            for item in result["artifacts"]:
                self.assertEqual(item["sha256"], digest(Path(item["path"])))
            frame_manifest = json.loads(artifacts["sampledFramesManifest"].read_text())
            self.assertEqual("protein-in-membrane.equilibration-frames.v1",
                             frame_manifest["schemaVersion"])
            self.assertEqual("float64-le-xyz-angstrom", frame_manifest["encoding"])
            self.assertEqual(("revision-1", "attempt-1", "a" * 32, "minimized-1", "a" * 64),
                             tuple(frame_manifest[key] for key in
                                   ("studyRevisionId", "attemptId", "stageId",
                                    "sourceMinimizedStageId", "protocolSha256")))
            self.assertEqual(digest(sidecar), frame_manifest["sourceTopologySha256"])
            self.assertEqual(digest(system_xml), frame_manifest["sourceSystemSha256"])
            self.assertEqual(digest(artifacts["equilibratedTopologyJson"]),
                             frame_manifest["finalTopologySha256"])
            self.assertEqual(digest(artifacts["equilibratedSystemXml"]),
                             frame_manifest["finalSystemSha256"])
            self.assertEqual(3, frame_manifest["frameCount"])
            self.assertEqual(3 * 5 * 24, frame_manifest["positionsByteLength"])
            frame_data = artifacts["sampledPositionsF64"].read_bytes()
            self.assertEqual(digest(artifacts["sampledPositionsF64"]),
                             frame_manifest["positionsSha256"])
            frames = np.frombuffer(frame_data, dtype="<f8").reshape(3, 5, 3)
            for index, item in enumerate(frame_manifest["frames"]):
                self.assertEqual(("unrestrained-observation", index + 1, index * 5 * 24, 5 * 24),
                                 (item["windowName"], item["step"], item["offsetBytes"],
                                  item["lengthBytes"]))
                self.assertEqual(hashlib.sha256(frame_data[index * 5 * 24:(index + 1) * 5 * 24])
                                 .hexdigest(), item["positionsSha256"])
            final_state = XmlSerializer.deserialize(artifacts["equilibratedStateXml"].read_text())
            final_system = XmlSerializer.deserialize(artifacts["equilibratedSystemXml"].read_text())
            final_box = final_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
            np.testing.assert_allclose(frames[-1], final_state.getPositions(asNumpy=True)
                                       .value_in_unit(unit.angstrom), atol=1e-10)
            np.testing.assert_allclose(frame_manifest["frames"][-1]["boxVectorsAngstrom"],
                                       final_box, atol=1e-10)
            system_box = np.asarray([vector.value_in_unit(unit.angstrom)
                                     for vector in final_system.getDefaultPeriodicBoxVectors()])
            final_sidecar = json.loads(artifacts["equilibratedTopologyJson"].read_text())
            self.assertFalse(np.allclose(_topology_data(topology)["boxVectorsAngstrom"], final_box,
                                         atol=1e-8), "The controlled barostat must exercise a changed cell")
            self.assertTrue(np.allclose(system_box, final_box, atol=1e-8))
            self.assertTrue(np.allclose(final_sidecar["boxVectorsAngstrom"], final_box, atol=1e-8))
            self.assertEqual(_topology_data(topology)["atoms"], final_sidecar["atoms"])
            self.assertEqual(3, len(result["observations"]["samples"]))
            self.assertEqual("completed", result["observations"]["termination"])
            self.assertTrue(result["observations"]["finalObservationUnrestrained"])
            self.assertEqual(1000.0, result["observations"]["windows"][0]["pressureBar"])
            self.assertTrue(any(kind == "equilibrationProgress" for kind, _ in progress))
            export_payload = {
                "topologyCifPath": str(artifacts["equilibratedCif"]),
                "topologyCifSha256": digest(artifacts["equilibratedCif"]),
                "topologyJsonPath": str(artifacts["equilibratedTopologyJson"]),
                "topologyJsonSha256": digest(artifacts["equilibratedTopologyJson"]),
                "systemXmlPath": str(artifacts["equilibratedSystemXml"]),
                "systemXmlSha256": digest(artifacts["equilibratedSystemXml"]),
                "stateXmlPath": str(artifacts["equilibratedStateXml"]),
                "stateXmlSha256": digest(artifacts["equilibratedStateXml"]),
                "coordinateReadBackToleranceAngstrom": 0.01,
                "cellLengthReadBackToleranceAngstrom": 0.01,
                "cellAngleReadBackToleranceDegrees": 0.1,
            }
            verified = verify_export(work, export_payload, lambda *_: None)
            self.assertTrue(verified["observations"]["readBackMatched"],
                            verified["observations"]["warnings"])
            old_topology = dict(export_payload, topologyJsonPath=str(sidecar),
                                topologyJsonSha256=digest(sidecar))
            mismatched = verify_export(work, old_topology, lambda *_: None)
            self.assertFalse(mismatched["observations"]["readBackMatched"])
            self.assertFalse(mismatched["observations"]["cellMatched"])
            first_hashes = {role: digest(path) for role, path in artifacts.items()}
            second_result = equilibrate(work, dict(payload, stageId="b" * 32), lambda *_: None)
            self.assertTrue(all(Path(item["path"]).parent == work / ("equilibration-" + "b" * 32)
                                for item in second_result["artifacts"]))
            self.assertEqual(first_hashes, {role: digest(path) for role, path in artifacts.items()})
            with self.assertRaises(WorkError) as duplicate:
                equilibrate(work, payload, lambda *_: None)
            self.assertEqual("resultConflict", duplicate.exception.code)
            self.assertEqual(first_hashes, {role: digest(path) for role, path in artifacts.items()})
            with self.assertRaises(WorkError) as invalid:
                equilibrate(work, dict(payload, stageId="../escape"), lambda *_: None)
            self.assertEqual("invalidRequest", invalid.exception.code)

            # The cross-language NVT wire spelling is "none": null pressure,
            # frequency and tension mean no pressure-coupling force at all.
            nvt_control = dict(control, pressureBar=None, pressureMode="none",
                               barostatFrequencySteps=None, surfaceTensionBarNm=None)
            nvt_protocol = dict(protocol, id="synthetic-nvt-worker-mechanics",
                                stages=[nvt_control],
                                extensionWindow=dict(nvt_control, name="extension"))
            nvt_payload = dict(payload, stageId="c" * 32,
                               validatedPolicyId=nvt_protocol["id"], protocol=nvt_protocol)
            nvt_result = equilibrate(work, nvt_payload, lambda *_: None)
            self.assertIsNone(nvt_result["observations"]["windows"][0]["pressureBar"])
            nvt_system_path = next(Path(item["path"]) for item in nvt_result["artifacts"]
                                   if item["role"] == "equilibratedSystemXml")
            nvt_system = XmlSerializer.deserialize(nvt_system_path.read_text())
            self.assertFalse(any(isinstance(nvt_system.getForce(index), MonteCarloMembraneBarostat)
                                 for index in range(nvt_system.getNumForces())))
            ramp_control = dict(nvt_control, name="heat", steps=3,
                                initialTemperatureKelvin=100.0,
                                proteinRestraintSelector="protein-heavy",
                                proteinRestraintKjMolNm2=2092.0,
                                lipidRestraintKjMolNm2=2092.0)
            ca_control = dict(nvt_control, name="ca-release", steps=3,
                              proteinRestraintSelector="protein-ca",
                              proteinRestraintKjMolNm2=418.4)
            ramp_protocol = dict(nvt_protocol, id="synthetic-ramp-worker-mechanics",
                                 stages=[ramp_control, ca_control, nvt_control],
                                 maximumSampleCount=9, maximumFrameBytes=9 * 5 * 24)
            ramp_progress = []
            ramp_result = equilibrate(work, dict(payload, stageId="f" * 32,
                validatedPolicyId=ramp_protocol["id"], protocol=ramp_protocol),
                lambda kind, item: ramp_progress.append((kind, item)))
            self.assertEqual(["heat", "ca-release", "unrestrained-observation"],
                             [window["name"] for window in ramp_result["observations"]["windows"]])
            self.assertEqual(9, len(ramp_result["observations"]["samples"]))
            self.assertEqual(42, next(item["velocityInitializationSeed"] for kind, item in ramp_progress
                                      if kind == "equilibrationSeedsUsed" and item["window"] == "heat"))
            legacy_hold = dict(nvt_control, name="legacy-backbone-hold",
                               proteinRestraintKjMolNm2=1.0)
            legacy_protocol = dict(nvt_protocol, id="synthetic-legacy-fixed-temperature",
                                   stages=[legacy_hold, nvt_control], maximumSampleCount=6,
                                   maximumFrameBytes=6 * 5 * 24)
            legacy_payload = dict(payload, stageId="9" * 32,
                                  validatedPolicyId=legacy_protocol["id"], protocol=legacy_protocol)
            legacy_payload.pop("proteinHeavyAtomIndices")
            legacy_result = equilibrate(work, legacy_payload, lambda *_: None)
            self.assertEqual(["legacy-backbone-hold", "unrestrained-observation"],
                             [window["name"] for window in legacy_result["observations"]["windows"]])
            refused_control = dict(nvt_control, barostatFrequencySteps=1)
            refused_protocol = dict(nvt_protocol, stages=[refused_control])
            with self.assertRaises(WorkError) as mixed_controls:
                equilibrate(work, dict(nvt_payload, stageId="d" * 32,
                                       protocol=refused_protocol), lambda *_: None)
            self.assertEqual("invalidProtocol", mixed_controls.exception.code)
            mismatched_mode = dict(nvt_control, pressureMode="constantVolume")
            with self.assertRaises(WorkError) as vocabulary_mismatch:
                equilibrate(work, dict(nvt_payload, stageId="e" * 32,
                                       protocol=dict(nvt_protocol, stages=[mismatched_mode],
                                                     extensionWindow=dict(mismatched_mode,
                                                                          name="extension"))),
                            lambda *_: None)
            self.assertEqual("invalidProtocol", vocabulary_mismatch.exception.code)
            stray_tension = dict(nvt_control, surfaceTensionBarNm=0.0)
            with self.assertRaises(WorkError) as nvt_tension_mismatch:
                equilibrate(work, dict(nvt_payload, stageId="f" * 32,
                                       protocol=dict(nvt_protocol, stages=[stray_tension],
                                                     extensionWindow=dict(stray_tension,
                                                                          name="extension"))),
                            lambda *_: None)
            self.assertEqual("invalidProtocol", nvt_tension_mismatch.exception.code)


if __name__ == "__main__":
    unittest.main()
