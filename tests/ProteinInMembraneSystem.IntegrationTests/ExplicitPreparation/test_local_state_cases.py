"""Independently calculable Slice 4 local-state observation cases.

These fixtures exercise the production observer on real OpenMM topology and
coordinate objects. Each expected distance follows directly from the small
coordinate arrangement, without using another production measurement helper.
"""

from __future__ import annotations

import copy
import math
from pathlib import Path
import sys
import unittest

from openmm import Vec3, unit
from openmm.app import Topology, element


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from ProteinInMembraneSystem.ExplicitPreparation.worker.local_state_observations import (  # noqa: E402
    observe_local_state,
)
from ProteinInMembraneSystem.worker.exchange import WorkError  # noqa: E402


def stage(atoms, *, bonds=(), cell=None):
    """Build one atom per residue; only declared bonds join molecules."""
    topology = Topology()
    entries = []
    mapping = []
    positions = []
    for index, (symbol, role, atom_role, side, xyz) in enumerate(atoms):
        chain = topology.addChain(str(index))
        residue = topology.addResidue("FIX", chain, id=str(index + 1))
        entries.append(topology.addAtom(symbol, element.get_by_symbol(symbol), residue))
        mapping.append({"resultAtomIndex": index, "moleculeRole": role,
                        "atomRole": atom_role, "physicalSide": side,
                        "element": symbol})
        positions.append(Vec3(*xyz))
    for first, second in bonds:
        topology.addBond(entries[first], entries[second])
    if cell is not None:
        topology.setUnitCellDimensions(Vec3(*cell) * unit.angstrom)
    return topology, positions * unit.angstrom, {"complete": True, "atoms": mapping}


def specification(*, metrics, pairs, radius=12.0, periodic=True, radii=None):
    return {"requiredMetricNames": list(metrics),
            "contactSearchRadiusAngstrom": radius, "maximumReportedPairs": 32,
            "usePeriodicBoundary": periodic,
            "contactRolePairs": [{"firstMoleculeRole": first, "secondMoleculeRole": second}
                                 for first, second in pairs],
            "atomRadiusByElementAngstrom": radii or
            {"C": 1.7, "H": 1.2, "O": 1.52, "P": 1.8, "Na": 1.0}}


def measures(result):
    return {item["name"]: item["value"] for item in result["measurements"]}


def pairs(result):
    return {(item["firstMoleculeRole"], item["secondMoleculeRole"]): item
            for item in result["rolePairMeasurements"]}


class LocalStateObservationCases(unittest.TestCase):
    def test_ion_and_water_seams_are_attributed_separately_from_protein_lipid_contact(self):
        # At x=0.3 and x=29.7/29.9 in a 30 A cell, minimum-image distances
        # are 0.6/0.4 A. The water-ion seam is still closer at 0.2 A.
        atoms = [
            ("C", "protein", "backbone", None, (0.3, 0, 30)),
            ("O", "water", "body", None, (29.7, 0, 30)),
            ("Na", "ion", "body", None, (29.9, 0, 30)),
            ("P", "lipid", "head", "upper", (2.3, 0, 40)),
            ("P", "lipid", "head", "lower", (2.3, 0, 20)),
        ]
        top, xyz, mapped = stage(atoms, cell=(30, 30, 60))
        spec = specification(metrics=("minimumIntermolecularDistanceAngstrom",
                                      "minimumIntermolecularHeavyAtomDistanceAngstrom",
                                      "leafletHeadSeparationAngstrom",
                                      "proteinBilayerMidplaneOffsetAngstrom"),
                             pairs=(("protein", "water"), ("protein", "ion"),
                                    ("water", "ion"), ("protein", "lipid")))
        observed = observe_local_state(top, xyz, mapped, spec)
        self.assertEqual("Observed", observed["standing"])
        self.assertEqual({"protein|water", "protein|ion", "water|ion", "protein|lipid"},
                         set(observed["coveredRolePairs"]))
        values = measures(observed)
        self.assertAlmostEqual(0.2, values["minimumIntermolecularDistanceAngstrom"])
        self.assertAlmostEqual(0.2, values["minimumIntermolecularHeavyAtomDistanceAngstrom"])
        self.assertAlmostEqual(20, values["leafletHeadSeparationAngstrom"])
        self.assertAlmostEqual(0, values["proteinBilayerMidplaneOffsetAngstrom"])
        by_role = pairs(observed)
        self.assertAlmostEqual(0.6, by_role[("protein", "water")]["minimumDistanceAngstrom"])
        self.assertAlmostEqual(0.4, by_role[("protein", "ion")]["minimumDistanceAngstrom"])
        self.assertAlmostEqual(0.2, by_role[("water", "ion")]["minimumDistanceAngstrom"])
        self.assertEqual(2, by_role[("protein", "lipid")]["pairsWithinSearchRadius"])
        self.assertAlmostEqual(math.sqrt(104), by_role[("protein", "lipid")]["minimumDistanceAngstrom"])

        # A common translation, including wrapping all coordinates back into
        # the cell, cannot change contact distances or relative organization.
        shifted = [(s, r, a, side,
                    ((point[0] + 10) % 30, (point[1] + 7) % 30,
                     (point[2] + 25) % 60))
                   for s, r, a, side, point in atoms]
        top2, xyz2, mapped2 = stage(shifted, cell=(30, 30, 60))
        translated = observe_local_state(top2, xyz2, mapped2, spec)
        self.assertEqual("Observed", translated["standing"])
        for metric, expected in values.items():
            self.assertAlmostEqual(expected, measures(translated)[metric])
        for key, expected in by_role.items():
            actual = pairs(translated)[key]
            self.assertEqual(expected["pairsWithinSearchRadius"], actual["pairsWithinSearchRadius"])
            self.assertAlmostEqual(expected["minimumDistanceAngstrom"], actual["minimumDistanceAngstrom"])

    def test_hydrogen_close_contact_is_distinct_from_heavy_overlap(self):
        atoms = [
            ("C", "protein", "body", None, (0.0, 0, 0)),
            ("H", "protein", "body", None, (1.0, 0, 0)),
            ("O", "water", "body", None, (1.8, 0, 0)),
            ("H", "water", "body", None, (1.2, 0, 0)),
        ]
        spec = specification(metrics=("minimumIntermolecularDistanceAngstrom",
                                      "minimumIntermolecularHeavyAtomDistanceAngstrom"),
                             pairs=(("protein", "water"),), radius=3, periodic=False)
        top, xyz, mapped = stage(atoms, bonds=((0, 1), (2, 3)))
        observed = observe_local_state(top, xyz, mapped, spec)
        self.assertEqual("Observed", observed["standing"])
        self.assertAlmostEqual(0.2, measures(observed)["minimumIntermolecularDistanceAngstrom"])
        self.assertAlmostEqual(1.8, measures(observed)["minimumIntermolecularHeavyAtomDistanceAngstrom"])
        self.assertGreaterEqual(measures(observed)["minimumIntermolecularHeavyAtomDistanceAngstrom"], 1.5)

        # Move only the water's heavy atom: this becomes a distinct severe
        # heavy contact even though hydrogen diagnostics were already short.
        heavy_overlap = list(atoms)
        heavy_overlap[2] = ("O", "water", "body", None, (1.4, 0, 0))
        top2, xyz2, mapped2 = stage(heavy_overlap, bonds=((0, 1), (2, 3)))
        changed = observe_local_state(top2, xyz2, mapped2, spec)
        self.assertEqual("Observed", changed["standing"])
        self.assertAlmostEqual(1.4, measures(changed)["minimumIntermolecularHeavyAtomDistanceAngstrom"])
        self.assertLess(measures(changed)["minimumIntermolecularHeavyAtomDistanceAngstrom"], 1.5)

    def test_crossed_tail_does_not_relabel_head_but_swapped_leaflet_identity_reverses_separation(self):
        atoms = [
            ("C", "protein", "backbone", None, (0, 0, 30)),
            ("P", "lipid", "head", "upper", (5, 0, 40)),
            ("C", "lipid", "tail", "upper", (5, 0, 18)),  # below midplane
            ("P", "lipid", "head", "lower", (5, 0, 20)),
            ("C", "lipid", "tail", "lower", (5, 0, 30)),
        ]
        spec = specification(metrics=("leafletHeadSeparationAngstrom",
                                      "proteinBilayerMidplaneOffsetAngstrom"),
                             pairs=(("protein", "lipid"),), radius=12)
        top, xyz, mapped = stage(atoms, bonds=((1, 2), (3, 4)), cell=(40, 40, 80))
        observed = observe_local_state(top, xyz, mapped, spec)
        self.assertEqual("Observed", observed["standing"])
        self.assertAlmostEqual(20, measures(observed)["leafletHeadSeparationAngstrom"])
        self.assertAlmostEqual(0, measures(observed)["proteinBilayerMidplaneOffsetAngstrom"])

        switched = copy.deepcopy(mapped)
        for entry in switched["atoms"]:
            if entry["moleculeRole"] == "lipid":
                entry["physicalSide"] = "lower" if entry["physicalSide"] == "upper" else "upper"
        inverted = observe_local_state(top, xyz, switched, spec)
        self.assertEqual("Observed", inverted["standing"])
        self.assertAlmostEqual(-20, measures(inverted)["leafletHeadSeparationAngstrom"])

    def test_missing_cell_mapping_coordinates_and_coverage_are_unavailable(self):
        atoms = [("C", "protein", "backbone", None, (0, 0, 30)),
                 ("P", "lipid", "head", "upper", (5, 0, 40)),
                 ("P", "lipid", "head", "lower", (5, 0, 20))]
        top, xyz, mapped = stage(atoms, cell=(40, 40, 80))
        spec = specification(metrics=("leafletHeadSeparationAngstrom",),
                             pairs=(("protein", "lipid"),), radius=12)
        self.assertEqual("Observed", observe_local_state(top, xyz, copy.deepcopy(mapped), spec)["standing"])

        cases = []
        no_cell, _, _ = stage(atoms)
        cases.append(("missing periodic topology", no_cell, xyz, mapped, spec, "lacks a cell"))
        too_short, _, _ = stage(atoms, cell=(20, 20, 80))
        cases.append(("periodic search cannot cover cell", too_short, xyz, mapped, spec,
                      "too small"))
        incomplete = copy.deepcopy(mapped)
        incomplete["complete"] = False
        cases.append(("incomplete correspondence", top, xyz, incomplete, spec, "incomplete"))
        duplicate = copy.deepcopy(mapped)
        duplicate["atoms"][2]["resultAtomIndex"] = 1
        cases.append(("duplicate atom index", top, xyz, duplicate, spec, "exactly once"))
        changed_element = copy.deepcopy(mapped)
        changed_element["atoms"][1]["element"] = "C"
        cases.append(("topology identity mismatch", top, xyz, changed_element, spec, "differs"))
        cases.append(("missing coordinate", top, xyz[:-1], mapped, spec, "absent"))
        nonfinite = [Vec3(0, 0, 30), Vec3(5, 0, float("nan")), Vec3(5, 0, 20)] * unit.angstrom
        cases.append(("nonfinite coordinate", top, nonfinite, mapped, spec, "nonfinite"))
        missing_radius = copy.deepcopy(spec)
        missing_radius["atomRadiusByElementAngstrom"].pop("P")
        cases.append(("uncovered element", top, xyz, mapped, missing_radius, "no declared radius"))
        missing_role = copy.deepcopy(spec)
        missing_role["contactRolePairs"] = [{"firstMoleculeRole": "protein", "secondMoleculeRole": "ion"}]
        cases.append(("uncovered role", top, xyz, mapped, missing_role, "no identified atoms"))
        missing_head = copy.deepcopy(mapped)
        missing_head["atoms"][2]["atomRole"] = "tail"
        cases.append(("missing lower head group", top, xyz, missing_head, spec,
                      "required protein-backbone or leaflet-head"))
        for name, case_top, case_xyz, case_map, case_spec, reason in cases:
            with self.subTest(name=name):
                result = observe_local_state(case_top, case_xyz, copy.deepcopy(case_map), case_spec)
                self.assertEqual("Unavailable", result["standing"])
                self.assertIn(reason, result["unavailableReason"])
                self.assertEqual([], result["measurements"])

    def test_ambiguous_backbone_image_empty_search_and_resource_bound_do_not_look_favorable(self):
        ambiguous = [("C", "protein", "backbone", None, (0, 0, 0)),
                     ("C", "protein", "backbone", None, (0, 0, 30)),
                     ("P", "lipid", "head", "upper", (5, 0, 40)),
                     ("P", "lipid", "head", "lower", (5, 0, 20))]
        top, xyz, mapped = stage(ambiguous, cell=(40, 40, 60))
        spec = specification(metrics=("proteinBilayerMidplaneOffsetAngstrom",),
                             pairs=(("protein", "lipid"),), radius=12)
        result = observe_local_state(top, xyz, mapped, spec)
        self.assertEqual("Unavailable", result["standing"])
        self.assertIn("ambiguous", result["unavailableReason"])

        separated = [("C", "protein", "body", None, (0, 0, 0)),
                     ("O", "water", "body", None, (20, 0, 0))]
        top2, xyz2, mapped2 = stage(separated, cell=(50, 50, 50))
        short = specification(metrics=("minimumIntermolecularDistanceAngstrom",),
                              pairs=(("protein", "water"),), radius=4)
        result = observe_local_state(top2, xyz2, mapped2, short)
        self.assertEqual("Unavailable", result["standing"])
        self.assertIn("No intermolecular atom pair", result["unavailableReason"])
        self.assertEqual([], result["measurements"])

        limited = copy.deepcopy(short)
        limited["maximumReportedPairs"] = 10_001
        with self.assertRaises(WorkError) as caught:
            observe_local_state(top2, xyz2, mapped2, limited)
        self.assertEqual("resourceRefused", caught.exception.code)


if __name__ == "__main__":
    unittest.main()
