"""Exact native-reference headgroup frames, including inverted patch orientation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "ProteinInMembrane.Host"))
from ProteinInMembraneSystem.PlacementAssessment.worker.placement_assessment import _headgroup_boundary  # noqa: E402
from ProteinInMembraneSystem.worker.exchange import WorkError  # noqa: E402


class ManualFrameTests(unittest.TestCase):
    def test_all_identified_pure_templates_supply_an_orientation_independent_frame(self):
        catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            extents = {}
            for representation in catalogue["lipids"]:
                species = representation["speciesId"]
                if species == "CHL1":
                    continue
                source = ROOT / "config/policies" / representation["coordinateTemplatePath"]
                target = work / f"{species}.cif"
                shutil.copyfile(source, target)
                descriptor = dict(representation, coordinateTemplatePath=target.name)
                self.assertEqual(hashlib.sha256(target.read_bytes()).hexdigest(),
                                 descriptor["coordinateTemplateSha256"])
                extents[species] = _headgroup_boundary(work, [descriptor])
                self.assertGreater(extents[species], 5.0)
                self.assertLess(extents[species], 80.0)
            self.assertGreater(extents["DPPC"], 18.0,
                               "The native DPPC reference stores its phosphate below its tails")

    def test_head_index_must_identify_the_actual_phosphate(self):
        catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
        representation = next(item for item in catalogue["lipids"] if item["speciesId"] == "DPPC")
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = ROOT / "config/policies" / representation["coordinateTemplatePath"]
            target = work / "DPPC.cif"
            shutil.copyfile(source, target)
            incorrect = dict(representation, coordinateTemplatePath=target.name,
                             headAtomIndices=[representation["headAtomIndices"][0] + 1])
            with self.assertRaises(WorkError):
                _headgroup_boundary(work, [incorrect])


if __name__ == "__main__":
    unittest.main()
