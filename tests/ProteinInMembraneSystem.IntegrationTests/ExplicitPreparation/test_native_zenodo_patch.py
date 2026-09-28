"""Exact CC-BY POPC source conversion and worker provider binding."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import openmm
from openmm.app import PDBxFile


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from ProteinInMembraneSystem.ExplicitPreparation.worker.construction import (  # noqa: E402
    _native_provider_identity, _native_verified_mapped_popc_patch,
)
from ProteinInMembraneSystem.worker.exchange import WorkError  # noqa: E402


ASSETS = ROOT / "config/policies/membrane-templates"
MANIFEST = json.loads((ASSETS / "POPC-Lipid21-Zenodo-64x64-map.json").read_text())


class MappedZenodoPatchTests(unittest.TestCase):
    def test_exact_source_mapping_and_every_lipid_match_selected_chemistry(self):
        source = ASSETS / MANIFEST["sourceFile"]
        mapped = ASSETS / MANIFEST["mappedFile"]
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), MANIFEST["sourceSha256"])
        self.assertEqual(hashlib.sha256(mapped.read_bytes()).hexdigest(), MANIFEST["mappedSha256"])
        catalogue = json.loads((ROOT / "config/policies/protein-membrane-current.json").read_text())
        lipid = next(item for item in catalogue["lipids"] if item["speciesId"] == "POPC")
        reference = PDBxFile(str(ROOT / "config/policies" / lipid["coordinateTemplatePath"]))
        _native_verified_mapped_popc_patch(PDBxFile(str(mapped)), reference,
                                           lipid["stereoChecks"])
        with tempfile.TemporaryDirectory() as temporary:
            reproduction = Path(temporary) / "reproduced.cif"
            result = subprocess.run([sys.executable,
                                     str(ROOT / "scripts/map-zenodo-popc-patch.py"),
                                     str(reproduction)], capture_output=True, text=True,
                                    timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(hashlib.sha256(reproduction.read_bytes()).hexdigest(),
                             MANIFEST["mappedSha256"])

    def test_provider_binds_both_pinned_source_and_mapped_bytes(self):
        payload = {"providerName": "OpenMM Modeller.addMembrane",
                   "providerVersion": openmm.version.full_version,
                   "lipidTypeArgument": "POPC",
                   "nativePatchMode": "mapped-lipid21-zenodo-popc",
                   "nativeSourcePatchPath": str(ASSETS / MANIFEST["sourceFile"]),
                   "nativeSourcePatchSha256": MANIFEST["sourceSha256"],
                   "nativePatchPath": str(ASSETS / MANIFEST["mappedFile"]),
                   "nativePatchSha256": MANIFEST["mappedSha256"]}
        path, digest, version = _native_provider_identity(payload)
        self.assertEqual(path, (ASSETS / MANIFEST["mappedFile"]).resolve())
        self.assertEqual(digest, MANIFEST["mappedSha256"])
        self.assertEqual(version, openmm.version.full_version)
        for changed in (dict(payload, nativeSourcePatchSha256="0" * 64),
                        dict(payload, nativePatchSha256="0" * 64),
                        dict(payload, lipidTypeArgument="DPPC")):
            with self.subTest(changed=changed), self.assertRaises(WorkError) as refused:
                _native_provider_identity(changed)
            self.assertEqual(refused.exception.code, "providerMismatch")


if __name__ == "__main__":
    unittest.main()
