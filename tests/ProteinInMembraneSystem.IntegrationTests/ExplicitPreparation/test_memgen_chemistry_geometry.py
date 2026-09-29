"""Controlled negative mutations for the independent Memgen artifact oracle."""

import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np

from inspect_memgen_chemistry_geometry import (
    InspectionError,
    component_image_contacts,
    graph_map,
    periodic_contacts,
    stereo_value,
)
from inspect_memgen_route import (
    raw_aqueous_region_account,
    recipe_blocks,
    retained_input_matches_intermediate,
    retained_provenance_matches,
)


class ChemistryGeometryMutationTests(unittest.TestCase):
    def test_retained_ion_pdb_abbreviation_is_narrow_and_final_copy_binding_is_exact(self):
        source = [{"name": "CA", "residue": "ALA", "chain": "A", "residueId": "1",
                   "insertionCode": "", "element": "C", "xyz": (1., 2., 3.)},
                  {"name": "NA", "residue": "NA", "chain": "B", "residueId": "2",
                   "insertionCode": "", "element": "NA", "xyz": (4., 5., 6.)},
                  {"name": "CL", "residue": "CL", "chain": "B", "residueId": "3",
                   "insertionCode": "", "element": "CL", "xyz": (7., 8., 9.)}]
        intermediate = copy.deepcopy(source)
        intermediate[1]["element"] = "N"
        intermediate[2]["element"] = "C"
        self.assertTrue(retained_input_matches_intermediate(source, intermediate))
        for index, field, changed in ((0, "element", "N"), (1, "name", "NX"),
                                      (1, "residueId", "4"), (2, "xyz", (7., 8., 8.))):
            with self.subTest(index=index, field=field):
                altered = copy.deepcopy(intermediate)
                altered[index][field] = changed
                self.assertFalse(retained_input_matches_intermediate(source, altered))
        prepared = {"sourceAtomId": "B:3:CL", "role": "source",
                    "sourceResidue": {"chain": "B", "residue": 3, "copyId": "B"},
                    "approvedChangeId": None, "moleculeRole": "ion",
                    "atomRole": "partnerAtom", "element": "Cl"}
        final = dict(prepared, generatedComponentRole=None, generatedSpeciesId=None)
        self.assertTrue(retained_provenance_matches(prepared, final))
        for field, changed in (("element", "C"), ("sourceAtomId", "B:4:CL"),
                               ("sourceResidue", {"chain": "A", "residue": 3, "copyId": "A"}),
                               ("generatedComponentRole", "negativeIon")):
            with self.subTest(field=field):
                altered = copy.deepcopy(final)
                altered[field] = changed
                self.assertFalse(retained_provenance_matches(prepared, altered))

    def test_raw_aqueous_planes_and_inventory_reject_mutations_and_measure_deviation(self):
        recipe = """pbc -30. -30. -45. 30. 30. 45.
structure PROT0.pdb
 number 1
end structure
structure WAT.pdb
 number 1
 below plane 0. 0. 1. -23.0
end structure
structure Na+.pdb
 number 1
 above plane 0. 0. 1. 23.0
end structure
structure Cl-.pdb
 number 1
 below plane 0. 0. 1. -23.0
end structure
"""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "packmol.inp"
            path.write_text(recipe, encoding="ascii")
            blocks, _ = recipe_blocks(path)
            molecules = [
                {"blockIndex": 1, "species": "WAT", "side": "lower", "atoms": [
                    {"name": "O", "xyz": (0., 0., -24.)},
                    {"name": "H1", "xyz": (0., 0., -23.1)},
                    {"name": "H2", "xyz": (0., 0., -23.01)}]},
                {"blockIndex": 2, "species": "Na+", "side": "upper", "atoms": [
                    {"name": "Na+", "xyz": (0., 0., 23.1)}]},
                {"blockIndex": 3, "species": "Cl-", "side": "lower", "atoms": [
                    {"name": "Cl-", "xyz": (0., 0., -23.1)}]},
            ]
            accepted, _, measured = raw_aqueous_region_account(blocks, molecules)
            self.assertTrue(accepted)
            self.assertEqual(measured["examinedAtomCount"], 5)
            self.assertEqual(measured["insideEnvelopeAtomCount"], 0)
            misplaced = copy.deepcopy(molecules)
            misplaced[0]["atoms"][2]["xyz"] = (0., 0., -22.99)
            accepted, _, measured = raw_aqueous_region_account(blocks, misplaced)
            self.assertTrue(accepted)
            self.assertEqual(measured["insideEnvelopeAtomCount"], 1)
            self.assertAlmostEqual(measured["maximumInsideDepthAngstrom"], 0.01)
            self.assertEqual(measured["insideAtomsBySpeciesSide"], {"WAT:lower": 1})
            wrong_plane = copy.deepcopy(blocks)
            wrong_plane[1]["planes"][0]["height"] = -22.
            self.assertFalse(raw_aqueous_region_account(wrong_plane, molecules)[0])
            wrong_direction = copy.deepcopy(blocks)
            wrong_direction[2]["planes"][0]["direction"] = "below"
            self.assertFalse(raw_aqueous_region_account(wrong_direction, molecules)[0])
            extra_plane = copy.deepcopy(blocks)
            extra_plane[1]["planes"].append({"direction": "above", "normal": (0., 0., 1.),
                                             "height": 23.})
            self.assertFalse(raw_aqueous_region_account(extra_plane, molecules)[0])
            missing_ion = molecules[:-1]
            self.assertFalse(raw_aqueous_region_account(blocks, missing_ion)[0])

    def test_bond_graph_mutation_is_rejected(self):
        reference = (["C", "C", "C", "O"], ["A", "B", "C", "D"],
                     [{1}, {0, 2}, {1, 3}, {2}])
        changed = (["C", "C", "C", "O"], ["A", "B", "C", "D"],
                   [{1, 2, 3}, {0}, {0}, {0}])
        with self.assertRaises(InspectionError):
            graph_map(reference, changed)

    def test_tetrahedral_inversion_and_alkene_trans_mutations(self):
        tetra = np.asarray([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
                            [-1., -1., -1.]])
        unchanged = stereo_value("tetrahedral", tetra)
        inverted = tetra.copy()
        inverted[:, 0] *= -1
        self.assertGreater(abs(unchanged), 1.)
        self.assertAlmostEqual(inverted_value := stereo_value("tetrahedral", inverted),
                               -unchanged)
        self.assertLess(unchanged * inverted_value, 0.)
        cis = np.asarray([[-1., 1., 0.], [0., 0., 0.], [1., 0., 0.],
                          [2., 1., 0.]])
        trans = cis.copy()
        trans[3, 1] = -1.
        self.assertGreater(stereo_value("alkene", cis), 0.5)
        self.assertLess(stereo_value("alkene", trans), 0.5)

    def test_face_edge_corner_intercomponent_contacts_are_detected(self):
        cell = np.asarray([10., 10., 10.])
        for far, expected in [([9.8, 0.2, 0.2], "face"),
                              ([9.8, 9.8, 0.2], "edge"),
                              ([9.8, 9.8, 9.8], "corner")]:
            with self.subTest(expected=expected):
                screen = periodic_contacts(np.asarray([[0.2, 0.2, 0.2], far]),
                                           cell, ["C", "C"], ["protein", "lipid"],
                                           [0, 1])
                self.assertEqual(len(screen["violations"]), 1)
                self.assertEqual(screen["violations"][0]["imageKind"], expected)
                self.assertLess(screen["violations"][0]["distanceAngstrom"], 1.5)

    def test_same_component_self_image_contact_is_detected(self):
        # Every adjacent bond is 0.96 A, yet the spanning chain nearly meets
        # its own periodic image at x=10 A. Intercomponent screening skips it.
        points = np.asarray([[0.96 * i, 2., 2.] for i in range(11)])
        edges = [{j for j in (i - 1, i + 1) if 0 <= j < 11} for i in range(11)]
        cell = np.asarray([10., 10., 10.])
        self.assertEqual(periodic_contacts(points, cell, ["C"] * 11,
                                           ["protein"] * 11, [0] * 11)["violations"], [])
        images = component_image_contacts(points, cell, edges, ["C"] * 11,
                                          ["protein"] * 11, [list(range(11))])
        self.assertTrue(any(item["imageKind"] == "face" and
                            item["distanceAngstrom"] < 1.5 for item in images))


if __name__ == "__main__":
    unittest.main()
