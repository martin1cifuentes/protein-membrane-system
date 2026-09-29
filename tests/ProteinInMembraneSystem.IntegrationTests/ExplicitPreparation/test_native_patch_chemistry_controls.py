"""Independent positive and deliberately altered controls for native chemistry checks.

These tests establish detection of the reported source-patch defects, not
physical qualification of any proposed replacement membrane patch.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import (
    _native_atom_name, _native_molecule_bonds, _native_molecule_matches,
    _native_provider_identity,
)
from ProteinInMembraneSystem.worker.exchange import WorkError


class NativePatchChemistryControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from openmm.app import PDBFile, PDBxFile

        catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
        cls.representations = {item["speciesId"]: item for item in catalogue["lipids"]}
        cls.patch = {}
        cls.reference = {}
        for species in ("DOPC", "DPPC", "POPC"):
            path = ROOT / "out/python/lib/python3.11/site-packages/openmm/app/data" / f"{species}.pdb"
            if not path.is_file():
                raise unittest.SkipTest("The exact OpenMM native patches are not installed")
            cls.patch[species] = PDBFile(str(path))
            cls.reference[species] = PDBxFile(str(ROOT / "config/policies" /
                                                cls.representations[species]["coordinateTemplatePath"]))

    def check(self, species, residue, positions=None):
        reference = self.reference[species]
        expected = _native_molecule_bonds(reference.topology, next(reference.topology.residues()))
        return _native_molecule_matches(
            residue, reference, species, positions or self.patch[species].positions,
            self.representations[species]["stereoChecks"],
            _native_molecule_bonds(self.patch[species].topology, residue), expected)

    def test_source_defect_ids_and_known_correct_molecules(self):
        for species, residue_name, expected_bad in (
                ("DOPC", "DOP", {"4", "90"}),
                ("DPPC", "DPP", {"6", "70"}),
                ("POPC", "POP", {"109"})):
            with self.subTest(species=species):
                residues = [r for r in self.patch[species].topology.residues()
                            if r.name == residue_name]
                self.assertEqual(128, len(residues))
                bad = set()
                for residue in residues:
                    try:
                        self.check(species, residue)
                    except WorkError as error:
                        self.assertEqual("providerMismatch", error.code)
                        bad.add(residue.id)
                self.assertEqual(expected_bad, bad)
                # The single-molecule reference itself supplies a separate
                # positive atom/bond/stereo control for the same descriptor.
                reference = self.reference[species]
                exact = next(reference.topology.residues())
                self.assertEqual(len(list(exact.atoms())),
                                 len(_native_molecule_matches(
                                     exact, reference, species, reference.positions,
                                     self.representations[species]["stereoChecks"],
                                     _native_molecule_bonds(reference.topology, exact),
                                     _native_molecule_bonds(reference.topology, exact))))

    def test_deliberate_glycerol_inversion_is_rejected(self):
        from openmm import Vec3, unit

        for species, residue_name in (("DPPC", "DPP"), ("POPC", "POP")):
            residue = next(r for r in self.patch[species].topology.residues()
                           if r.name == residue_name)
            atoms = list(residue.atoms())
            by_name = {_native_atom_name(atom.name, species): atom.index for atom in atoms}
            self.check(species, residue)
            moved = list(self.patch[species].positions)
            left, right = by_name["C1"], by_name["C3"]
            moved[left], moved[right] = moved[right], moved[left]
            with self.subTest(species=species), self.assertRaises(WorkError) as refused:
                self.check(species, residue, moved)
            self.assertEqual("providerMismatch", refused.exception.code)

    def test_deliberate_cis_to_trans_change_is_rejected(self):
        from openmm import Vec3, unit

        species = "DOPC"
        residue = next(r for r in self.patch[species].topology.residues() if r.name == "DOP")
        atoms = list(residue.atoms())
        by_name = {_native_atom_name(atom.name, species): atom.index for atom in atoms}
        self.check(species, residue)
        # Reflect the fourth atom of one declared cis descriptor across the
        # central bond axis.  The worker must reject the observed trans form.
        descriptor = self.representations[species]["stereoChecks"][1]["atomNames"]
        a, b, c, d = (by_name[name] for name in descriptor)
        positions = self.patch[species].positions.value_in_unit(unit.angstrom)
        second, third, fourth = positions[b], positions[c], positions[d]
        axis = [third[i] - second[i] for i in range(3)]
        right = [fourth[i] - third[i] for i in range(3)]
        axial = sum(right[i] * axis[i] for i in range(3)) / sum(value * value for value in axis)
        reflected = [third[i] + 2 * axial * axis[i] - right[i] for i in range(3)]
        moved = list(self.patch[species].positions)
        moved[d] = Vec3(*reflected) * unit.angstrom
        with self.assertRaises(WorkError) as refused:
            self.check(species, residue, moved)
        self.assertEqual("providerMismatch", refused.exception.code)

    def test_superseded_deletion_derivatives_are_not_selectable(self):
        import hashlib
        import openmm

        cases = (("DOPC", "DOP", ["4", "90"],
                  "05bc09ae7c7ac6cbdb221868fede000f310b86999105fda1cfbbb262135c7f10",
                  "f1fe152d3c9bf22ac5eec70ff3d207d7506816ee03dcc910766b2f29d7fbddee"),
                 ("DPPC", "DPP", ["6", "70"],
                  "6020f03c08d38b8469c3074da826fdf979a1ea627354459a092e27558af503d6",
                  "c961627244da03a88aad66124fb8e877f84fa750e80b4f430f20137a975e1926"))
        for species, residue_name, removed, source_sha, derived_sha in cases:
            with self.subTest(species=species):
                source_path = ROOT / "out/python/lib/python3.11/site-packages/openmm/app/data" / f"{species}.pdb"
                derived_path = ROOT / "config/policies/membrane-templates" / f"{species}-OpenMM86-63x63.pdb"
                self.assertEqual(source_sha, hashlib.sha256(source_path.read_bytes()).hexdigest())
                self.assertEqual(derived_sha, hashlib.sha256(derived_path.read_bytes()).hexdigest())
                payload = {"providerName": "OpenMM Modeller.addMembrane",
                           "providerVersion": openmm.version.full_version,
                           "lipidTypeArgument": species,
                           "nativePatchMode": "balanced-defect-deletion",
                           "nativeSourcePatchPath": str(source_path),
                           "nativeSourcePatchSha256": source_sha,
                           "nativePatchPath": str(derived_path),
                           "nativePatchSha256": derived_sha,
                           "removedNativeLipidResidueIds": removed}
                with self.assertRaises(WorkError) as refused:
                    _native_provider_identity(payload)
                self.assertEqual("unsupportedPolicy", refused.exception.code)


if __name__ == "__main__":
    unittest.main()
