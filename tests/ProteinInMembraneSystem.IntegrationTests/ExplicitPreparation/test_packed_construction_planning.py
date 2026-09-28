"""Focused constrained-packing inputs and real small Packmol exchange.

The small 20-lipid case checks identity transport and constraints.  It is not
evidence of a complete hydrated, parameterized or minimized product route.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from ProteinInMembraneSystem.ExplicitPreparation.worker.packed_construction import (
    _check_protein_images, _head_and_tails, _packmol_text, _protein_clearance_bounds,
    _read_xyz, _template_block, _write_xyz, allocate_counts,
)
from ProteinInMembraneSystem.worker.exchange import WorkError


class PackedConstructionPlanningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from openmm.app import PDBxFile

        catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
        cls.representations = {item["speciesId"]: item for item in catalogue["lipids"]}
        cls.references = {
            species: PDBxFile(str(ROOT / "config/policies" /
                                  cls.representations[species]["coordinateTemplatePath"]))
            for species in ("DLPC", "DMPC", "CHL1")}

    def test_largest_remainder_retains_every_positive_species_without_silent_rounding(self):
        self.assertEqual({"DLPC": 24, "DMPC": 23},
                         allocate_counts([("DLPC", 0.5), ("DMPC", 0.5)], 47))
        self.assertEqual({"DLPC": 1, "DMPC": 46},
                         allocate_counts([("DLPC", 0.02), ("DMPC", 0.98)], 47))
        with self.assertRaises(WorkError):
            allocate_counts([("DLPC", 0.001), ("DMPC", 0.001), ("POPE", 0.998)], 47)

    def test_cholesterol_hydroxyl_index_and_bond_mapping(self):
        raw = self.representations["CHL1"]
        reference = self.references["CHL1"]
        head, tails, span = _head_and_tails(raw, reference)
        self.assertEqual(28, head)
        self.assertEqual("O3", list(reference.topology.atoms())[head - 1].name)
        self.assertEqual(3, len(tails))
        self.assertGreater(span, 15)
        with self.assertRaises(WorkError):
            _head_and_tails(dict(raw, headAtomIndices=[27]), reference)
        with self.assertRaises(WorkError):
            _head_and_tails(dict(raw, areaPerMoleculeAngstromSquared=0), reference)

    def test_clearance_uses_whole_construct_and_refuses_resource_limit(self):
        points = [(-5.0, -4.0, -2.0), (7.0, 8.0, 3.0)]
        xmin, ymin, zmin, xmax, ymax, zmax = _protein_clearance_bounds(
            points, [17.0, 21.0], 1, 180, 1)
        self.assertGreaterEqual(xmin * -1 - 5, 10)
        self.assertGreaterEqual(xmax - 7, 10)
        self.assertGreaterEqual(ymax - 8, 10)
        self.assertGreaterEqual(8 - ymin, 10)
        self.assertGreaterEqual(zmax - zmin, 20)
        shifted = [(x + 19, y - 7, z) for x, y, z in points]
        sxmin, symin, _, sxmax, symax, _ = _protein_clearance_bounds(
            shifted, [17.0, 21.0], 1, 180, 1)
        self.assertAlmostEqual(sxmin, -sxmax)
        self.assertAlmostEqual(symin, -symax)
        self.assertGreaterEqual(sxmax - max(x for x, _, _ in shifted), 10)
        self.assertGreaterEqual(min(x for x, _, _ in shifted) - sxmin, 10)
        self.assertGreaterEqual(symax - max(y for _, y, _ in shifted), 10)
        self.assertGreaterEqual(min(y for _, y, _ in shifted) - symin, 10)
        self.assertEqual(26, _check_protein_images(points,
                         (xmin, ymin, zmin, xmax, ymax, zmax), 20))
        with self.assertRaises(WorkError):
            _protein_clearance_bounds(points, [21.0], 1, 40, 1)

    def test_periodic_protein_images_check_faces_edges_and_corners(self):
        box = (-10, -10, -10, 10, 10, 10)
        for name, first, second in (
                ("face", (-9.8, 0, 0), (9.8, 0, 0)),
                ("edge", (-9.8, -9.8, 0), (9.8, 9.8, 0)),
                ("corner", (-9.8, -9.8, -9.8), (9.8, 9.8, 9.8))):
            with self.subTest(name=name), self.assertRaisesRegex(WorkError, name):
                _check_protein_images([first, second], box, 1.5)
        self.assertEqual(26, _check_protein_images([(-1, -1, -1), (1, 1, 1)],
                                                   box, 1.5))

    def test_small_real_packmol_preserves_blocks_and_head_tail_sides(self):
        from openmm import unit

        executable = ROOT / "out/reconciled-design/packmol/venv311/bin/packmol"
        if not executable.is_file():
            raise unittest.SkipTest("The identified Packmol 21.2.3 executable is not installed")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            paths = {}
            for species in ("DLPC", "DMPC"):
                reference = self.references[species]
                points = [tuple(float(value) for value in position.value_in_unit(unit.angstrom))
                          for position in reference.positions]
                center = [sum(point[axis] for point in points) / len(points) for axis in range(3)]
                path = directory / f"{species}.xyz"
                symbols = tuple(atom.element.symbol for atom in reference.topology.atoms())
                _write_xyz(path, symbols,
                           [tuple(point[axis] - center[axis] for axis in range(3))
                            for point in points], species)
                paths[species] = path
            blocks = [_template_block(species, self.references[species],
                                      self.representations[species], side, 5)
                      for side in ("upper", "lower") for species in ("DLPC", "DMPC")]
            recipe, expected = _packmol_text(directory, blocks, paths, [],
                                             (-25, -25, -60, 25, 25, 60), 20260927)
            input_path = directory / "input.inp"
            input_path.write_text(recipe)
            result = subprocess.run([str(executable), "-i", str(input_path)],
                                    cwd=directory, text=True, capture_output=True,
                                    timeout=60, check=False)
            self.assertEqual(0, result.returncode, result.stdout[-800:])
            self.assertIn("Success!", result.stdout)
            points = _read_xyz(directory / "packed.xyz", expected)
            self.assertEqual(2240, len(points))
            offset = 0
            for block in blocks:
                for _ in range(block.count):
                    molecule = points[offset:offset + len(block.symbols)]
                    offset += len(block.symbols)
                    head_z = molecule[block.head_index - 1][2]
                    tails_z = [molecule[index - 1][2] for index in block.tail_indices]
                    tail_z = sum(tails_z) / len(tails_z)
                    self.assertGreater((head_z - tail_z) * (1 if block.side == "upper" else -1), 0)
                    sign = 1 if block.side == "upper" else -1
                    self.assertGreaterEqual(head_z * sign + 0.02,
                                            max(8, block.head_tail_span - 1))
                    self.assertTrue(all(z * sign <= 4 + 0.02 for z in tails_z))
            self.assertEqual(offset, len(points))


if __name__ == "__main__":
    unittest.main()
