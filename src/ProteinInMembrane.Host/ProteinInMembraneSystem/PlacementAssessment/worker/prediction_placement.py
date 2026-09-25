"""Directional PAE observations for a proposed protein placement."""

from __future__ import annotations

import heapq
import json
from pathlib import Path
from typing import Any

from ProteinInMembraneSystem.worker.exchange import WorkError, require_integer, require_mapping, require_text, verify_sha256, work_path
from ProteinInMembraneSystem.worker.pae_matrix import _pae

def summarize_prediction_evidence(directory: Path, payload: dict[str, Any], progress: Any) -> dict[str, Any]:
    pae_path = work_path(directory, payload.get("paePath"), "paePath")
    pae_hash = verify_sha256(pae_path, require_text(payload.get("paeSha256"), "paeSha256"), "paeSha256")
    mapping_path = work_path(directory, payload.get("paeMappingPath"), "paeMappingPath")
    verify_sha256(mapping_path, require_text(payload.get("paeMappingSha256"), "paeMappingSha256"),
                  "paeMappingSha256")
    mapping = require_mapping(json.loads(mapping_path.read_text(encoding="utf-8")), "PAE mapping")
    record_id = require_text(payload.get("recordId"), "recordId")
    if mapping.get("recordId") != record_id or str(mapping.get("coordinateSha256", "")).lower() != str(payload.get("coordinateSha256", "")).lower() or str(mapping.get("paeSha256", "")).lower() != pae_hash.lower():
        raise WorkError("inputMismatch", "PAE axis mapping does not identify this prediction and coordinate source")
    matrix = _pae(pae_path)
    addresses = mapping.get("axisResidues")
    if not isinstance(addresses, list) or len(addresses) != len(matrix):
        raise WorkError("invalidEvidence", "PAE axis mapping and matrix lengths disagree")
    index = {json.dumps(address, sort_keys=True): ordinal for ordinal, address in enumerate(addresses)}
    if len(index) != len(addresses):
        raise WorkError("invalidEvidence", "PAE axis mapping repeats a residue identity")
    first = payload.get("firstRegion")
    second = payload.get("secondRegion")
    maximum_pairs = require_integer(payload.get("maximumReportedPairs"), "maximumReportedPairs", 1)
    if not isinstance(first, list) or not first or not isinstance(second, list) or not second or maximum_pairs > 1000:
        raise WorkError("invalidRequest", "PAE regions must be nonempty and reported pairs bounded by 1000")
    if (len({json.dumps(item, sort_keys=True) for item in first}) != len(first) or
            len({json.dumps(item, sort_keys=True) for item in second}) != len(second)):
        raise WorkError("invalidRequest", "A PAE region repeats an exact residue identity")
    if len(first) * len(second) > 2_000_000:
        raise WorkError("resourceRefused", "Requested PAE region-pair summary exceeds its bounded comparison size")
    first_indices = [index.get(json.dumps(item, sort_keys=True)) for item in first]
    second_indices = [index.get(json.dumps(item, sort_keys=True)) for item in second]
    if any(value is None for value in first_indices + second_indices):
        observations = {"recordId": record_id, "standing": "Unavailable",
                        "reason": "A requested residue is absent from the exact PAE axis mapping.",
                        "firstAlignedOnSecond": None, "secondAlignedOnFirst": None}
    else:
        def directional(first_aligned_on_second: bool):
            count = 0
            total = 0.0
            minimum = None
            maximum = None
            high: list[tuple[float, int, int]] = []
            for a, i in enumerate(first_indices):
                for b, j in enumerate(second_indices):
                    # AFDB convention: matrix[alignment residue][target residue].
                    # FirstAlignedOnSecond aligns on j; SecondAlignedOnFirst aligns on i.
                    # An asymmetric 2x2 fixture with matrix[0][1] != matrix[1][0]
                    # must distinguish these directions in Implementation tests.
                    value = float(matrix[j][i] if first_aligned_on_second else matrix[i][j])
                    count += 1
                    total += value
                    minimum = value if minimum is None else min(minimum, value)
                    maximum = value if maximum is None else max(maximum, value)
                    entry = (value, a, b)
                    if len(high) < maximum_pairs:
                        heapq.heappush(high, entry)
                    elif entry > high[0]:
                        heapq.heapreplace(high, entry)
            return {"possiblePairCount": len(first) * len(second), "validPairCount": count,
                    "minimumAngstrom": minimum,
                    "maximumAngstrom": maximum,
                    "meanAngstrom": total / count,
                    "highestErrorPairs": [{"first": first[a], "second": second[b],
                                           "predictedAlignedErrorAngstrom": value}
                                          for value, a, b in sorted(high, reverse=True)]}
        observations = {"recordId": record_id, "standing": "Observed", "reason": None,
                        "firstAlignedOnSecond": directional(True),
                        "secondAlignedOnFirst": directional(False)}
    return {"artifacts": [], "observations": observations,
            "provider": {"name": "AlphaFold DB PAE JSON", "version": "identified source artifact"}}
