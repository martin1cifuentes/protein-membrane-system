"""Source-identified AlphaFold DB confidence and PAE-axis observations."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from ProteinInMembraneSystem.worker.exchange import WorkError, artifact, require_mapping, require_text, sha256, verify_sha256, work_path
from ProteinInMembraneSystem.worker.pae_matrix import _pae

def _unavailable(prediction: dict[str, Any], coordinate_hash: str, confidence: list[dict[str, Any]], reason: str):
    return {
        "recordId": prediction["recordId"], "coordinateSha256": coordinate_hash,
        "localConfidence": confidence, "paeStanding": "Unavailable", "paeReason": reason,
        "paeAxisResidueCount": None, "paeMappingPath": None, "paeMappingSha256": None,
        "limitations": [],
    }, []


def inspect_prediction(directory: Path, payload: dict[str, Any], structure: Any,
                       source_hash: str, peptide_subchains: set[str]):
    """Observe pLDDT/PAE only for an exactly identified predicted source."""
    if payload.get("sourceKind") != "alphafold":
        return None, []
    prediction = require_mapping(payload.get("prediction"), "prediction")
    record_id = require_text(prediction.get("recordId"), "prediction.recordId")
    if str(prediction.get("coordinateSha256", "")).lower() != source_hash.lower():
        raise WorkError("inputMismatch", "Prediction identity does not match staged coordinates")

    # The AFDB mmCIF B_iso field contains pLDDT. It is not read as pLDDT for
    # experimental or arbitrary uploaded structures, where B factors mean other things.
    from .protein_preparation import _address, _residue_kind

    confidence: list[dict[str, Any]] = []
    axis: list[tuple[int, dict[str, Any]]] = []
    if len(structure) != 1:
        for model_index, model in enumerate(structure):
            for chain in model:
                for residue in chain:
                    if _residue_kind(residue, peptide_subchains) == "protein":
                        confidence.append({"residue": _address(chain.name, residue, model_index),
                                           "pLddt": None, "standing": "Unavailable",
                                           "reason": "Multiple coordinate models prevent a unique prediction-confidence mapping."})
        return _unavailable(prediction, source_hash, confidence,
                            "The identified prediction has multiple coordinate models; PAE axes cannot be assigned uniquely.")
    for chain in structure[0]:
        for residue in chain:
            if _residue_kind(residue, peptide_subchains) != "protein":
                continue
            address = _address(chain.name, residue)
            values = [float(atom.b_iso) for atom in residue if atom.element.name not in {"H", "D"}]
            if values and all(math.isfinite(value) and 0 <= value <= 100 for value in values) and max(values)-min(values) <= 0.02:
                confidence.append({"residue": address, "pLddt": sum(values)/len(values),
                                   "standing": "Observed", "reason": None})
            else:
                confidence.append({"residue": address, "pLddt": None, "standing": "Unavailable",
                                   "reason": "Residue B-factor fields do not provide one consistent 0–100 confidence value."})
            label_seq = getattr(residue, "label_seq", None)
            if isinstance(label_seq, int) and label_seq > 0:
                axis.append((label_seq, address))

    if prediction.get("paeStanding") != "Available":
        return _unavailable(prediction, source_hash, confidence,
                            prediction.get("paeReason") or "Identified PAE data were not acquired.")
    try:
        pae_path = work_path(directory, prediction.get("paePath"), "prediction.paePath")
        pae_hash = verify_sha256(pae_path, require_text(prediction.get("paeSha256"), "prediction.paeSha256"),
                                 "prediction.paeSha256")
        matrix = _pae(pae_path)
    except (WorkError, ValueError, json.JSONDecodeError) as exc:
        return _unavailable(prediction, source_hash, confidence, f"Identified PAE could not be interpreted: {exc}")
    if len(axis) != len(confidence) or len(axis) != len(matrix):
        return _unavailable(prediction, source_hash, confidence,
                            "PAE axis length does not match uniquely numbered protein residues.")
    axis.sort(key=lambda entry: entry[0])
    labels = [entry[0] for entry in axis]
    if len(set(labels)) != len(labels) or labels != list(range(labels[0], labels[0] + len(labels))):
        return _unavailable(prediction, source_hash, confidence,
                            "Protein label_seq identifiers do not give one contiguous PAE axis.")
    sequence_start = prediction.get("sequenceStart")
    sequence_end = prediction.get("sequenceEnd")
    if isinstance(sequence_start, int) and isinstance(sequence_end, int) and sequence_end-sequence_start+1 != len(axis):
        return _unavailable(prediction, source_hash, confidence,
                            "Identified prediction span does not match the PAE axis length.")

    mapping_path = directory / "pae-axis-mapping.json"
    mapping_path.write_text(json.dumps({"recordId": record_id, "coordinateSha256": source_hash,
                                        "paeSha256": pae_hash, "axisResidues": [address for _, address in axis],
                                        "method": "unique contiguous mmCIF label_seq order"},
                                       separators=(",", ":")), encoding="utf-8")
    return ({"recordId": record_id, "coordinateSha256": source_hash,
             "localConfidence": confidence, "paeStanding": "Observed", "paeReason": None,
             "paeAxisResidueCount": len(axis), "paeMappingPath": str(mapping_path),
             "paeMappingSha256": sha256(mapping_path), "limitations": []},
            [artifact(directory, mapping_path, "paeAxisMapping")])
