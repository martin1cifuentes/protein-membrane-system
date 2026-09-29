"""The adopted construction frame determines manual upper and lower positions."""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "ProteinInMembrane.Host"))
from ProteinInMembraneSystem.PlacementAssessment.worker.placement_assessment import place_manual  # noqa: E402
from ProteinInMembraneSystem.worker.exchange import WorkError  # noqa: E402

SOURCE = Path(__file__).parent / "fixtures" / "6qwr-prepared.pdb"


def _heavy_z(path: Path) -> list[float]:
    return [float(line[46:54]) for line in path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip() not in {"H", "D"}]


class ManualFrameTests(unittest.TestCase):
    def _place(self, position: str, envelope: float) -> tuple[Path, dict]:
        self.temporary = tempfile.TemporaryDirectory()
        work = Path(self.temporary.name)
        prepared = work / "prepared.pdb"
        shutil.copyfile(SOURCE, prepared)
        payload = {
            "preparedPdbPath": str(prepared),
            "preparedSha256": hashlib.sha256(prepared.read_bytes()).hexdigest(),
            "preparedProteinId": "6qwr-controlled-frame",
            "preparedAtomCount": sum(line.startswith(("ATOM  ", "HETATM"))
                                     for line in prepared.read_text().splitlines()),
            "maximumAtomCount": 120000,
            "startingPosition": position,
            "offsetXAngstrom": 0.0, "offsetYAngstrom": 0.0, "offsetZAngstrom": 0.0,
            "rotationXDegrees": 0.0, "rotationYDegrees": 0.0, "rotationZDegrees": 0.0,
            "leafletEnvelopeAngstrom": envelope,
        }
        result = place_manual(work, payload, lambda *_: None)
        return work / "manual-placement.pdb", result

    def tearDown(self) -> None:
        if hasattr(self, "temporary"):
            self.temporary.cleanup()

    def test_selected_23_angstrom_envelope_anchors_both_physical_sides(self) -> None:
        upper, result = self._place("upper", 23.0)
        self.assertAlmostEqual(23.0, result["observations"]["headgroupBoundaryAngstrom"])
        self.assertAlmostEqual(25.0, min(_heavy_z(upper)), places=2)
        self.temporary.cleanup()
        del self.temporary
        lower, result = self._place("lower", 23.0)
        self.assertAlmostEqual(23.0, result["observations"]["headgroupBoundaryAngstrom"])
        self.assertAlmostEqual(-25.0, max(_heavy_z(lower)), places=2)

    def test_absent_frame_cannot_position_a_construct(self) -> None:
        with self.assertRaises(WorkError) as caught:
            self._place("upper", 0.0)
        self.assertEqual("missingFrame", caught.exception.code)


if __name__ == "__main__":
    unittest.main()
