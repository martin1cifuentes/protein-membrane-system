"""Exact pinned POPC patch and construction-boundary checks without assembly."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src/ProteinInMembrane.Host"
sys.path.insert(0, str(SOURCE_ROOT))

from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import (
    _native_atom_name, _native_lipid_side, _native_molecule_bonds,
    _native_molecule_matches, _native_provider_identity,
    _native_verified_custom_popc_patch,
)
from ProteinInMembraneSystem.worker.exchange import WorkError


class NativePopcMechanics(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from openmm import unit
        from openmm.app import PDBFile, PDBxFile

        catalogue = json.loads((ROOT / "config/policies/protein-membrane-slice4.json").read_text())
        cls.lipid = next(item for item in catalogue["lipids"] if item["speciesId"] == "POPC")
        cls.patch_path = ROOT / "out/python/lib/python3.11/site-packages/openmm/app/data/POPC.pdb"
        cls.custom_path = ROOT / "config/policies/membrane-templates/POPC-OpenMM86-63x63.pdb"
        cls.reference_path = ROOT / "config/policies" / cls.lipid["coordinateTemplatePath"]
        if not cls.patch_path.is_file():
            raise unittest.SkipTest("The identified OpenMM POPC patch is not installed")
        for path, expected in ((cls.patch_path,
                                "a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29"),
                               (cls.reference_path, cls.lipid["coordinateTemplateSha256"])):
            assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        cls.patch = PDBFile(str(cls.patch_path))
        cls.reference = PDBxFile(str(cls.reference_path))
        cls.lipids = [residue for residue in cls.patch.topology.residues() if residue.name == "POP"]
        cls.canonical = next(cls.reference.topology.residues())
        cls.expected_bonds = _native_molecule_bonds(cls.reference.topology, cls.canonical)
        cls.points = [tuple(float(value) for value in point.value_in_unit(unit.angstrom))
                      for point in cls.patch.positions]
        cls.center = cls.patch.topology.getUnitCellDimensions().value_in_unit(unit.angstrom)[2] / 2

    def matched(self, residue, positions=None):
        return _native_molecule_matches(
            residue, self.reference, "POPC", positions or self.patch.positions,
            self.lipid["stereoChecks"], _native_molecule_bonds(self.patch.topology, residue),
            self.expected_bonds)

    def test_exact_native_popc_reference_cis_geometry_and_physical_leaflets(self) -> None:
        self.assertEqual(128, len(self.lipids))
        self.assertEqual(134, len(list(self.canonical.atoms())))
        sides = Counter()
        mismatches = []
        for residue in self.lipids:
            try:
                atoms = self.matched(residue)
            except WorkError as error:
                mismatches.append((residue.id, error.code))
                continue
            sides[_native_lipid_side(atoms, self.points, self.center, "POPC")] += 1
        # This installed patch has one actual inverted sn-2 glycerol center.
        # The mechanism must report it rather than treating source-patch
        # availability as qualification of every retained POPC molecule.
        self.assertEqual([("109", "providerMismatch")], mismatches)
        self.assertEqual(Counter({"upper": 64, "lower": 63}), sides)

    def test_changed_cis_alkene_and_wrong_leaflet_orientation_refuse(self) -> None:
        from openmm import Vec3, unit

        residue = self.lipids[0]
        atoms = self.matched(residue)
        by_name = {_native_atom_name(atom.name, "POPC"): atom for atom in atoms}
        self.assertEqual({"C218", "C316"},
                         {name for name in by_name if name in {"C218", "C316"}})
        second = self.points[by_name["C29"].index]
        third = self.points[by_name["C210"].index]
        fourth = self.points[by_name["C211"].index]
        axis = [third[i] - second[i] for i in range(3)]
        right = [fourth[i] - third[i] for i in range(3)]
        axial = sum(right[i] * axis[i] for i in range(3)) / sum(value * value for value in axis)
        reflected = tuple(third[i] + 2 * axial * axis[i] - right[i] for i in range(3))
        moved = list(self.patch.positions)
        moved[by_name["C211"].index] = Vec3(*reflected) * unit.angstrom
        with self.assertRaises(WorkError) as alkene:
            self.matched(residue, moved)
        self.assertEqual("providerMismatch", alkene.exception.code)

        wrong_side = list(self.points)
        head_z = self.points[by_name["P"].index][2]
        for name in ("C218", "C316"):
            x, y, _ = wrong_side[by_name[name].index]
            wrong_side[by_name[name].index] = (x, y, head_z + 10)
        with self.assertRaises(WorkError) as leaflet:
            _native_lipid_side(atoms, wrong_side, self.center, "POPC")
        self.assertEqual("providerMismatch", leaflet.exception.code)

    def test_provider_requires_exact_installed_popc_resource(self) -> None:
        import openmm

        payload = {"providerName": "OpenMM Modeller.addMembrane",
                   "providerVersion": openmm.version.full_version,
                   "lipidTypeArgument": "POPC", "nativePatchPath": str(self.patch_path),
                   "nativePatchSha256": "a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29"}
        observed_path, observed_sha, _ = _native_provider_identity(payload)
        self.assertEqual(self.patch_path.resolve(), observed_path)
        self.assertEqual(payload["nativePatchSha256"], observed_sha)
        for changed, reason in ((dict(payload, nativePatchSha256="0" * 64), "inputMismatch"),
                                (dict(payload, nativePatchPath=str(self.reference_path)), "providerMismatch"),
                                (dict(payload, lipidTypeArgument="POPE"), "unsupportedPolicy")):
            with self.subTest(changed=changed), self.assertRaises(WorkError) as refused:
                _native_provider_identity(changed)
            self.assertEqual(reason, refused.exception.code)

    def test_identified_custom_patch_is_exact_balanced_deletion_only(self) -> None:
        from openmm import Vec3, unit
        from openmm.app import PDBFile

        self.assertEqual("92a19c3470605d44ca779772f306ad693afdd11f1f78176655a3f1d2ee3a77b4",
                         hashlib.sha256(self.custom_path.read_bytes()).hexdigest())
        custom = PDBFile(str(self.custom_path))
        _native_verified_custom_popc_patch(self.patch, custom, self.reference,
                                           self.lipid["stereoChecks"], ["62", "109"])
        self.assertEqual(Counter({"POP": 126, "HOH": 5120}),
                         Counter(residue.name for residue in custom.topology.residues()))
        shifted = list(custom.positions)
        shifted[0] = shifted[0] + Vec3(1, 0, 0) * unit.angstrom
        with self.assertRaises(WorkError) as altered:
            _native_verified_custom_popc_patch(self.patch,
                SimpleNamespace(topology=custom.topology, positions=shifted),
                self.reference, self.lipid["stereoChecks"], ["62", "109"])
        self.assertEqual("providerMismatch", altered.exception.code)

        changed_cell = PDBFile(str(self.custom_path))
        vectors = changed_cell.topology.getPeriodicBoxVectors()
        changed_cell.topology.setPeriodicBoxVectors(
            (vectors[0], vectors[1], vectors[2] + Vec3(1, 0, 0) * unit.angstrom))
        with self.assertRaises(WorkError) as angle:
            _native_verified_custom_popc_patch(self.patch, changed_cell,
                self.reference, self.lipid["stereoChecks"], ["62", "109"])
        self.assertEqual("providerMismatch", angle.exception.code)

    def test_custom_provider_binds_installed_source_and_separate_derived_digest(self) -> None:
        import openmm

        payload = {"providerName": "OpenMM Modeller.addMembrane",
                   "providerVersion": openmm.version.full_version,
                   "lipidTypeArgument": "POPC", "nativePatchMode": "popc-62-109-deletion",
                   "nativeSourcePatchPath": str(self.patch_path),
                   "nativeSourcePatchSha256": "a35daa948562a67c142ea5ba38c1e18cc7a50bc1afc559338b2565d23fd38a29",
                   "nativePatchPath": str(self.custom_path),
                   "nativePatchSha256": "92a19c3470605d44ca779772f306ad693afdd11f1f78176655a3f1d2ee3a77b4",
                   "removedNativeLipidResidueIds": ["62", "109"]}
        observed_path, observed_sha, _ = _native_provider_identity(payload)
        self.assertEqual(self.custom_path.resolve(), observed_path)
        self.assertEqual(payload["nativePatchSha256"], observed_sha)
        for changed, reason in ((dict(payload, nativeSourcePatchSha256="0" * 64), "inputMismatch"),
                                (dict(payload, nativeSourcePatchPath=str(self.custom_path)), "providerMismatch"),
                                (dict(payload, nativePatchSha256="0" * 64), "inputMismatch"),
                                (dict(payload, nativePatchMode="other"), "unsupportedPolicy"),
                                (dict(payload, removedNativeLipidResidueIds=["109"]), "unsupportedPolicy")):
            with self.subTest(changed=changed), self.assertRaises(WorkError) as refused:
                _native_provider_identity(changed)
            self.assertEqual(reason, refused.exception.code)


if __name__ == "__main__":
    unittest.main()
