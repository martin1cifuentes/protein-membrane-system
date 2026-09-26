"""One-request JSON crossing between the host and the local scientific worker.

Large structures are attempt-scoped files. Standard output is reserved for
newline-delimited JSON events; providers' output is captured in files.

Wire request: {requestId, operation, workingDirectory, payload}. All payloads
are lower-camel serialization of the corresponding C# *Payload record in
ProteinInMembraneSystem.Contracts.cs. Source coordinates, force-field XML, coordinate templates,
System/State XML and generated molecular files are attempt-scoped under
workingDirectory; the host stages and hashes them before requesting work.
Only the explicit PPM executable and the verified OpenMM-installed membrane
patch may be read outside that directory. Operations are inspect_source,
inspect_preparation_changes, assess_membrane, prepare_protein, place_ppm,
adjust_placement, measure_placement, summarize_prediction_evidence,
construct_system, minimize, equilibrate, observe_stage and verify_export. Each
process reads exactly one request and writes progress events
followed by one terminal result or error event. The terminal payload echoes
requestId and any study/attempt/stage IDs, with observed facts and artifact
hashes or a distinct failure code; it never asserts C# product qualification.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable


@dataclass(slots=True)
class WorkError(Exception):
    code: str
    message: str
    details: dict[str, Any] | None = None

    def __str__(self) -> str:
        return self.message


def require_mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise WorkError("invalidRequest", f"{name} must be an object")
    return value


def require_text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkError("invalidRequest", f"{name} must be nonempty text")
    return value


def require_integer(value: Any, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise WorkError("invalidRequest", f"{name} must be an integer >= {minimum}")
    return value


def require_number(value: Any, name: str, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkError("invalidRequest", f"{name} must be a number")
    result = float(value)
    if not __import__("math").isfinite(result) or (minimum is not None and result < minimum):
        raise WorkError("invalidRequest", f"{name} is outside its finite range")
    return result


def work_path(directory: Path, value: Any, name: str, *, must_exist: bool = True) -> Path:
    """Resolve only files staged in the request's own working directory."""
    raw = require_text(value, name)
    path = (directory / raw).resolve()
    if path == directory or directory not in path.parents:
        raise WorkError("invalidPath", f"{name} must stay within the work directory")
    if must_exist and not path.is_file():
        raise WorkError("missingInput", f"{name} is unavailable", {"path": raw})
    return path


def verify_sha256(path: Path, expected: Any, name: str) -> str:
    actual = sha256(path)
    if expected is not None and actual.lower() != require_text(expected, name).lower():
        raise WorkError("inputMismatch", f"{name} does not match the staged file")
    return actual


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(directory: Path, path: Path, role: str) -> dict[str, str]:
    resolved = path.resolve()
    if directory not in resolved.parents or not resolved.is_file():
        raise WorkError("unobservedOutput", f"{role} artifact was not produced in the work directory")
    return {"role": role, "path": str(resolved), "sha256": sha256(resolved)}


def emit(request_id: str, kind: str, payload: dict[str, Any]) -> None:
    sys.stdout.write(json.dumps({"requestId": request_id, "kind": kind, "payload": payload}, allow_nan=False) + "\n")
    sys.stdout.flush()


def _dispatch(operation: str) -> Callable[[Path, dict[str, Any], Callable[[str, dict[str, Any]], None]], dict[str, Any]]:
    owner_modules = {
        "inspect_source": "ProteinPreparation.worker.protein_preparation",
        "inspect_preparation_changes": "ProteinPreparation.worker.protein_preparation",
        "prepare_protein": "ProteinPreparation.worker.protein_preparation",
        "assess_membrane": "MembraneModelAssessment.worker.membrane_model",
        "place_ppm": "PlacementAssessment.worker.placement_assessment",
        "adjust_placement": "PlacementAssessment.worker.placement_assessment",
        "measure_placement": "PlacementAssessment.worker.placement_assessment",
        "summarize_prediction_evidence": "PlacementAssessment.worker.prediction_placement",
        "construct_system": "ExplicitPreparation.worker.construction",
        "observe_stage": "ExplicitPreparation.worker.construction",
        "minimize": "ExplicitPreparation.Minimization.worker.minimization",
        "equilibrate": "ExplicitPreparation.OptionalEquilibrationProcedure.worker.equilibration",
        "verify_export": "CompletedStageExport.worker.export_verification",
    }
    owner_module = owner_modules.get(operation)
    if owner_module is None:
        raise WorkError("unsupportedOperation", f"Unsupported worker operation: {operation}")
    source_root_value = os.environ.get("PIM_OWNER_SOURCE_ROOT")
    if not source_root_value:
        raise WorkError("missingSource", "The owner-local Python source root was not supplied")
    source_root = Path(source_root_value).resolve()
    if not (source_root / "ProteinInMembraneSystem").is_dir():
        raise WorkError("missingSource", "The owner-local Python source tree is unavailable")
    sys.path.insert(0, str(source_root))
    module = importlib.import_module(f"ProteinInMembraneSystem.{owner_module}")
    return getattr(module, operation)


def main() -> int:
    request_id = "unidentified"
    identity: dict[str, Any] = {}
    try:
        line = sys.stdin.readline()
        if not line:
            raise WorkError("invalidRequest", "No work request was supplied")
        try:
            request = require_mapping(json.loads(line), "request")
        except (ValueError, TypeError) as exc:
            raise WorkError("invalidRequest", f"Invalid request JSON: {exc}") from exc
        request_id = require_text(request.get("requestId"), "requestId")
        operation = require_text(request.get("operation"), "operation")
        directory = Path(require_text(request.get("workingDirectory"), "workingDirectory")).resolve()
        if not directory.is_dir() or directory == Path(directory.anchor):
            raise WorkError("invalidPath", "workingDirectory must be an existing attempt directory")
        payload = require_mapping(request.get("payload"), "payload")
        identity = {key: payload[key] for key in ("studyRevisionId", "attemptId", "stageId") if key in payload}

        def progress(stage: str, detail: dict[str, Any] | None = None) -> None:
            emit(request_id, "progress", {**identity, "operation": operation, "stage": stage, "detail": detail or {}})

        progress("started")
        result = _dispatch(operation)(directory, payload, progress)
        emit(request_id, "result", {"requestId": request_id, **identity, "standing": "observed", **result})
        return 0
    except WorkError as exc:
        emit(request_id, "error", {"requestId": request_id, **identity, "standing": "failed", "failureCode": exc.code,
                                    "failureMessage": exc.message, "details": exc.details or {}})
        return 2
    except ImportError as exc:
        emit(request_id, "error", {"requestId": request_id, **identity, "standing": "failed",
                                   "failureCode": "dependencyUnavailable", "failureMessage": str(exc)})
        return 3
    except Exception as exc:
        emit(request_id, "error", {"requestId": request_id, **identity, "standing": "failed",
                                   "failureCode": "mechanicsFailed", "failureMessage": str(exc)})
        return 4
