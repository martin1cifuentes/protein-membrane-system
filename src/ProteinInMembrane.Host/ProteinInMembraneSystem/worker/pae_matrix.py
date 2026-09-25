"""Bounded parsing of an identified PAE matrix for source and placement observations."""

from __future__ import annotations

import json
import math
from pathlib import Path

from ProteinInMembraneSystem.worker.exchange import WorkError

_MAX_PAE_AXIS = 8_000
_MAX_PAE_BYTES = 100_000_000

def _pae(path: Path) -> list[list[float]]:
    if path.stat().st_size > _MAX_PAE_BYTES:
        raise WorkError("resourceRefused", "Identified PAE exceeds the bounded evidence size")
    with path.open("r", encoding="utf-8") as stream:
        data = json.load(stream)
    if not isinstance(data, list) or len(data) != 1 or not isinstance(data[0], dict):
        raise WorkError("invalidEvidence", "PAE JSON is not the identified AlphaFold DB matrix format")
    matrix = data[0].get("predicted_aligned_error")
    if not isinstance(matrix, list) or not matrix or len(matrix) > _MAX_PAE_AXIS:
        raise WorkError("invalidEvidence", "PAE matrix axis is absent or exceeds the bounded size")
    n = len(matrix)
    maximum = data[0].get("max_predicted_aligned_error")
    if not isinstance(maximum, (int, float)) or not math.isfinite(float(maximum)) or maximum <= 0:
        raise WorkError("invalidEvidence", "PAE maximum is not a finite positive value")
    for row in matrix:
        if not isinstance(row, list) or len(row) != n:
            raise WorkError("invalidEvidence", "PAE matrix is not square")
        if any(not isinstance(value, (int, float)) or not math.isfinite(float(value)) or
               value < 0 or value > maximum for value in row):
            raise WorkError("invalidEvidence", "PAE matrix contains a nonfinite or out-of-range value")
    return matrix
