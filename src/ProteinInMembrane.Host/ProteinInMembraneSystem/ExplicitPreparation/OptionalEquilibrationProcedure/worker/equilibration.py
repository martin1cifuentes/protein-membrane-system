"""Owner-local staged equilibration and adequacy observations.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import json
import hashlib
import math
import re
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (WorkError, artifact, sha256, require_integer, require_mapping,
                                        require_number, require_text)
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _provider, _load_stage, _rms_force, _stage_measurements, _final_state, _topology_data,
    _write_stage_coordinates)


def _index_selection(value: Any, name: str, system: Any) -> list[int]:
    if not isinstance(value, list):
        raise WorkError("invalidRequest", f"{name} must be an explicit atom-index array")
    indices = [require_integer(item, name) for item in value]
    if len(set(indices)) != len(indices) or any(item >= system.getNumParticles() for item in indices):
        raise WorkError("invalidRequest", f"{name} contains duplicate or out-of-range atom indices")
    return indices


def _restraint(system: Any, indices: list[int], constant: float, reference: Any, name: str) -> None:
    from openmm import CustomExternalForce, unit

    if constant == 0:
        return
    if not indices:
        raise WorkError("invalidProtocol", f"{name} requests restraints without identified target atoms")
    force = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.setName(name)
    for parameter in ("k", "x0", "y0", "z0"):
        force.addPerParticleParameter(parameter)
    positions = reference.value_in_unit(unit.nanometer)
    for index in indices:
        x, y, z = positions[index]
        force.addParticle(index, [constant, x, y, z])
    system.addForce(force)


def _protein_restraint_indices(selector: str, atoms: list[Any], backbone: list[int],
                               protein_heavy: list[int]) -> list[int]:
    heavy = set(protein_heavy)
    backbone_heavy = [index for index in backbone if index in heavy]
    targets = {
        "backbone": backbone,
        "protein-heavy": protein_heavy,
        "backbone-heavy": backbone_heavy,
        "protein-ca": [index for index in backbone_heavy
                       if atoms[index].name == "CA" and atoms[index].element.symbol == "C"],
    }.get(selector)
    if targets is None:
        raise WorkError("invalidProtocol", "Protein restraint selector is not declared by this worker")
    if not targets:
        raise WorkError("invalidProtocol", f"Protein restraint selector {selector} has no exact target atoms")
    return targets


def _initial_temperature(control: dict[str, Any]) -> float | None:
    start_raw = control.get("initialTemperatureKelvin")
    if start_raw is None:
        return None
    start = require_number(start_raw, "initialTemperatureKelvin", 0.000001)
    target = require_number(control.get("temperatureKelvin"), "temperatureKelvin", 0.000001)
    steps = require_integer(control.get("steps"), "stage steps", 2)
    if start == target or control.get("pressureBar") is not None:
        raise WorkError("invalidProtocol", "A temperature ramp needs distinct endpoints and constant volume")
    return start


def _advance_temperature_control(integrator: Any, initial: float | None, target: float,
                                 total_steps: int, completed: int, chunk: int) -> None:
    if initial is None:
        integrator.step(chunk)
        return
    from openmm import unit

    # The first and last integration steps use the exact declared endpoints.
    for step in range(completed, completed + chunk):
        kelvin = initial + (target - initial) * step / (total_steps - 1)
        integrator.setTemperature(kelvin * unit.kelvin)
        integrator.step(1)


def _barostat(system: Any, control: dict[str, Any]):
    from openmm import MonteCarloMembraneBarostat, unit

    pressure = control.get("pressureBar")
    if pressure is None:
        if (str(control.get("pressureMode", "none")).lower() != "none" or
                control.get("barostatFrequencySteps") is not None or
                control.get("surfaceTensionBarNm") is not None):
            raise WorkError("invalidProtocol", "An NVT stage must not request pressure-coupling controls")
        return None
    pressure = require_number(pressure, "pressureBar", 0)
    frequency = require_integer(control.get("barostatFrequencySteps"), "barostatFrequencySteps", 1)
    tension = require_number(control.get("surfaceTensionBarNm"), "surfaceTensionBarNm")
    mode = require_text(control.get("pressureMode"), "pressureMode").lower()
    modes = {
        "xyisotropiczfree": (MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree),
        "xyisotropiczfixed": (MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFixed),
        "xyanisotropiczfree": (MonteCarloMembraneBarostat.XYAnisotropic, MonteCarloMembraneBarostat.ZFree),
    }
    if mode not in modes:
        raise WorkError("unsupportedProtocol", "Pressure mode is not one of the explicit membrane-barostat modes")
    xy, z = modes[mode]
    temperature = require_number(control.get("temperatureKelvin"), "stage temperatureKelvin", 0.000001)
    barostat = MonteCarloMembraneBarostat(pressure * unit.bar, tension * unit.bar * unit.nanometer,
                                          temperature * unit.kelvin, xy, z, frequency)
    system.addForce(barostat)
    return barostat


def _instant_temperature(state: Any, system: Any) -> float | None:
    from openmm import CMMotionRemover, unit

    massive = sum(system.getParticleMass(index).value_in_unit(unit.dalton) > 0
                  for index in range(system.getNumParticles()))
    dof = 3 * massive - system.getNumConstraints()
    if any(isinstance(system.getForce(index), CMMotionRemover) for index in range(system.getNumForces())):
        dof -= 3
    if dof <= 0:
        return None
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    value = 2 * kinetic / (dof * 0.00831446261815324)
    return value if math.isfinite(value) else None


def _equilibration_observables(protocol: dict[str, Any], payload: dict[str, Any],
                              system: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate policy-selected measurements; never infer scientific thresholds here."""
    raw_observables = protocol.get("observables")
    resolved_observables = payload.get("resolvedObservables")
    raw_rules = protocol.get("sufficiencyRules")
    if (not isinstance(raw_observables, list) or not raw_observables or
            not isinstance(resolved_observables, list) or
            len(resolved_observables) != len(raw_observables) or not isinstance(raw_rules, list)):
        raise WorkError("missingPolicy", "Equilibration needs declared observables and sufficiency rules")
    method_units = {"potentialEnergy": "kJ/mol", "kineticEnergy": "kJ/mol",
                    "rmsForce": "kJ/mol/nm", "temperature": "K",
                    "cellAreaXY": "angstrom^2", "cellHeight": "angstrom",
                    "centroidZ": "angstrom", "centroidSeparationZ": "angstrom",
                    "proteinMidplaneOffset": "angstrom",
                    "countWithinDistance": "count"}
    methods = {"potentialEnergy", "kineticEnergy", "rmsForce", "temperature", "cellAreaXY", "cellHeight",
               "centroidZ", "centroidSeparationZ", "proteinMidplaneOffset", "countWithinDistance"}
    categories = {"system", "membraneOrganization", "proteinPlacement", "hydrationIons"}
    observables = []
    names = set()
    selectors = {"none", "protein-backbone", "upper-lipid-head", "lower-lipid-head",
                 "all-lipid-head", "water-oxygen", "all-ions", "protein-or-lipid-heavy"}
    for raw, resolved in zip(raw_observables, resolved_observables):
        declaration = require_mapping(raw, "declared equilibration observable")
        item = require_mapping(resolved, "resolved equilibration observable")
        for field in ("name", "unit", "scope", "category", "method", "distanceCutoffAngstrom"):
            if declaration.get(field) != item.get(field):
                raise WorkError("inputMismatch", f"Resolved observable disagrees with policy field {field}")
        first_selector = require_text(declaration.get("atomSelector"), "observable atomSelector")
        second_selector = require_text(declaration.get("comparisonSelector"), "observable comparisonSelector")
        third_selector = require_text(declaration.get("thirdSelector"), "observable thirdSelector")
        if first_selector not in selectors or second_selector not in selectors or third_selector not in selectors:
            raise WorkError("unsupportedProtocol", "An observable selector is not supported")
        name = require_text(item.get("name"), "observable name")
        method = require_text(item.get("method"), "observable method")
        category = require_text(item.get("category"), "observable category")
        if name in names or method not in methods or category not in categories:
            raise WorkError("unsupportedProtocol", "An observable is duplicated or has an unsupported method/category")
        names.add(name)
        if require_text(item.get("unit"), "observable unit") != method_units[method]:
            raise WorkError("invalidProtocol", f"{name} unit does not match its measurement method")
        require_text(item.get("scope"), "observable scope")
        indices = _index_selection(item.get("atomIndices", []), f"{name} atomIndices", system)
        comparison = _index_selection(item.get("comparisonAtomIndices", []), f"{name} comparisonAtomIndices", system)
        third = _index_selection(item.get("thirdAtomIndices", []), f"{name} thirdAtomIndices", system)
        if (first_selector == "none") != (len(indices) == 0) or \
                (second_selector == "none") != (len(comparison) == 0) or \
                (third_selector == "none") != (len(third) == 0):
            raise WorkError("inputMismatch", f"{name} resolved indices disagree with declared selectors")
        if method in {"centroidZ", "centroidSeparationZ", "proteinMidplaneOffset", "countWithinDistance"} and not indices:
            raise WorkError("invalidProtocol", f"{name} requires exact atom indices")
        if method in {"centroidSeparationZ", "proteinMidplaneOffset", "countWithinDistance"} and not comparison:
            raise WorkError("invalidProtocol", f"{name} requires exact comparison atom indices")
        if method == "proteinMidplaneOffset" and not third:
            raise WorkError("invalidProtocol", f"{name} requires exact lower-leaflet atom indices")
        if method != "proteinMidplaneOffset" and third:
            raise WorkError("invalidProtocol", f"{name} has an unused third atom selection")
        if method not in {"centroidSeparationZ", "proteinMidplaneOffset", "countWithinDistance"} and comparison:
            raise WorkError("invalidProtocol", f"{name} has an unused comparison selection")
        if method not in {"centroidZ", "centroidSeparationZ", "proteinMidplaneOffset", "countWithinDistance"} and indices:
            raise WorkError("invalidProtocol", f"{name} has an unused atom selection")
        cutoff = item.get("distanceCutoffAngstrom")
        if method == "countWithinDistance":
            require_number(cutoff, f"{name} distanceCutoffAngstrom", 0.000001)
        elif cutoff is not None:
            raise WorkError("invalidProtocol", f"{name} has an unused distance cutoff")
        observables.append(item)
    if {item["category"] for item in observables} != categories:
        raise WorkError("missingPolicy", "An applicable policy needs all four observed condition categories")
    rules = []
    ruled = set()
    for raw in raw_rules:
        rule = require_mapping(raw, "equilibration sufficiency rule")
        name = require_text(rule.get("observableName"), "rule observableName")
        if name not in names or name in ruled:
            raise WorkError("invalidProtocol", "A sufficiency rule names no unique declared observable")
        ruled.add(name)
        require_integer(rule.get("blockSizeSamples"), f"{name} blockSizeSamples", 1)
        require_integer(rule.get("minimumEffectiveBlocks"), f"{name} minimumEffectiveBlocks", 1)
        require_number(rule.get("maximumAbsoluteFirstVsLastBlockMeanDifference"),
                       f"{name} maximumAbsoluteFirstVsLastBlockMeanDifference", 0)
        correlation = require_number(rule.get("maximumAbsoluteLagOneBlockCorrelation"),
                                     f"{name} maximumAbsoluteLagOneBlockCorrelation", 0)
        if correlation >= 1:
            raise WorkError("invalidProtocol", "Lag-one correlation bound must be below one")
        rules.append(rule)
    if ruled != names:
        raise WorkError("missingPolicy", "Every declared observable needs an explicit sufficiency rule")
    return observables, rules


def _count_within_distance(fractional_positions: Any, box: Any, inverse_box: Any,
                           first_indices: list[int], reference_indices: list[int],
                           cutoff: float) -> int:
    """Count first-group atoms near a second group in the periodic cell.

    Fractional bins are at least as wide as the real-space cutoff projected
    onto each reciprocal-cell axis. This limits each query to adjacent bins,
    including periodic seams, while retaining the worker's fractional
    minimum-image calculation.
    """
    import numpy as np

    reciprocal_bounds = cutoff * np.linalg.norm(inverse_box, axis=0)
    divisions = np.asarray([max(1, min(512, int(math.floor(1.0 / bound))))
                            for bound in reciprocal_bounds], dtype=np.int64)
    first = np.asarray(first_indices, dtype=np.int64)
    reference = np.asarray(reference_indices, dtype=np.int64)
    first_bins = np.minimum((fractional_positions[first] * divisions).astype(np.int64), divisions - 1)
    reference_bins = np.minimum((fractional_positions[reference] * divisions).astype(np.int64), divisions - 1)

    def keys(addresses: Any) -> Any:
        return (addresses[:, 0] * divisions[1] + addresses[:, 1]) * divisions[2] + addresses[:, 2]

    order = np.argsort(keys(reference_bins), kind="stable")
    sorted_keys = keys(reference_bins)[order]
    sorted_reference = reference[order]
    maximum_bin_population = int(np.unique(sorted_keys, return_counts=True)[1].max())
    batch_size = max(1, min(2048, 1_000_000 // maximum_bin_population))
    matched = np.zeros(len(first), dtype=bool)
    offsets = [(-1, 0, 1) if size > 2 else (0, 1) if size == 2 else (0,)
               for size in divisions]
    cutoff_squared = cutoff * cutoff
    for dx in offsets[0]:
        for dy in offsets[1]:
            for dz in offsets[2]:
                active = np.flatnonzero(~matched)
                if not len(active):
                    return len(first)
                # Bound temporary pair arrays if a broad cutoff covers most
                # of the cell; the usual local search remains vectorized.
                for begin in range(0, len(active), batch_size):
                    batch = active[begin:begin + batch_size]
                    shifted = (first_bins[batch] + (dx, dy, dz)) % divisions
                    query_keys = keys(shifted)
                    starts = np.searchsorted(sorted_keys, query_keys, side="left")
                    ends = np.searchsorted(sorted_keys, query_keys, side="right")
                    lengths = ends - starts
                    total = int(lengths.sum())
                    if not total:
                        continue
                    repeated = np.repeat(np.arange(len(batch)), lengths)
                    starts_repeated = np.repeat(starts, lengths)
                    within = np.arange(total) - np.repeat(np.cumsum(lengths) - lengths, lengths)
                    candidates = sorted_reference[starts_repeated + within]
                    pairs = batch[repeated]
                    delta = fractional_positions[candidates] - fractional_positions[first[pairs]]
                    delta -= np.rint(delta)
                    wrapped = delta @ box
                    hits = np.einsum("ij,ij->i", wrapped, wrapped) <= cutoff_squared
                    np.logical_or.at(matched, pairs, hits)
    return int(matched.sum())


def _observable_values(state: Any, system: Any, observables: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Measure only bounded, named scalars from an actual OpenMM state."""
    import numpy as np
    from openmm import unit

    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    area = float(np.linalg.norm(np.cross(box[0], box[1])))
    if not math.isfinite(area) or area <= 0:
        raise WorkError("nonfiniteObservation", "Periodic lateral area is unavailable")
    height = abs(float(np.linalg.det(box))) / area
    inverse_box = np.linalg.inv(box)
    fractional_positions = (positions @ inverse_box) % 1.0

    def periodic_centroid_normal(atom_indices: list[int]) -> float:
        # A circular mean along the periodic membrane normal avoids a false
        # jump when an identified group crosses the cell origin. This is a
        # coordinate measurement, not a protein-placement support judgment.
        fractions = fractional_positions[atom_indices, 2]
        angles = 2 * math.pi * fractions
        sine = float(np.mean(np.sin(angles)))
        cosine = float(np.mean(np.cos(angles)))
        if math.hypot(sine, cosine) <= 1e-12:
            raise WorkError("unobservedOutput", "A periodic group has no attributable normal centroid")
        return (math.atan2(sine, cosine) / (2 * math.pi)) % 1.0

    values = []
    for item in observables:
        method = item["method"]
        indices = item.get("atomIndices", [])
        comparison = item.get("comparisonAtomIndices", [])
        third = item.get("thirdAtomIndices", [])
        if method == "potentialEnergy":
            value = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        elif method == "kineticEnergy":
            value = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        elif method == "rmsForce":
            value = _rms_force(state, system.getNumParticles())
        elif method == "temperature":
            value = _instant_temperature(state, system)
        elif method == "cellAreaXY":
            value = area
        elif method == "cellHeight":
            value = height
        elif method == "centroidZ":
            value = periodic_centroid_normal(indices) * height
        elif method == "centroidSeparationZ":
            fractional_difference = periodic_centroid_normal(indices) - periodic_centroid_normal(comparison)
            value = ((fractional_difference + 0.5) % 1.0 - 0.5) * height
        elif method == "proteinMidplaneOffset":
            upper = periodic_centroid_normal(comparison)
            lower = periodic_centroid_normal(third)
            leaflet_separation = (upper - lower + 0.5) % 1.0 - 0.5
            midpoint = (lower + leaflet_separation / 2) % 1.0
            protein_difference = periodic_centroid_normal(indices) - midpoint
            value = ((protein_difference + 0.5) % 1.0 - 0.5) * height
        elif method == "countWithinDistance":
            value = _count_within_distance(fractional_positions, box, inverse_box,
                                           indices, comparison,
                                           float(item["distanceCutoffAngstrom"]))
        else:
            raise WorkError("unsupportedProtocol", f"Unsupported observable method {method}")
        if value is None or not math.isfinite(float(value)):
            raise WorkError("nonfiniteObservation", f"Observable {item['name']} is unavailable or nonfinite")
        values.append({"name": item["name"], "value": float(value),
                       "unit": item["unit"], "scope": item["scope"]})
    return values


def _assess_equilibration_samples(samples: list[dict[str, Any]], rules: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], bool]:
    """Fixed-block, lag-one observation adequacy; not a scientific equilibrium verdict."""
    assessments = []
    for rule in rules:
        name = rule["observableName"]
        trace = [measurement["value"] for sample in samples
                 for measurement in sample["measurements"] if measurement["name"] == name]
        width = rule["blockSizeSamples"]
        blocks = [sum(trace[index:index + width]) / width
                  for index in range(0, len(trace) - width + 1, width)]
        # At least three complete blocks are needed to observe a lag-one
        # correlation. Nullable statistics distinguish absence from zero.
        effective = 0.0
        difference = correlation = None
        sufficient = False
        if len(blocks) >= 3:
            mean = sum(blocks) / len(blocks)
            variance = sum((value - mean) ** 2 for value in blocks)
            if variance == 0:
                correlation = 0.0
            else:
                correlation = sum((blocks[index] - mean) * (blocks[index + 1] - mean)
                                  for index in range(len(blocks) - 1)) / variance
                correlation = max(-1.0, min(1.0, correlation))
            positive = max(0.0, correlation)
            effective = len(blocks) * (1 - positive) / (1 + positive)
            difference = blocks[-1] - blocks[0]
            sufficient = (effective >= rule["minimumEffectiveBlocks"] and
                          abs(difference) <= rule["maximumAbsoluteFirstVsLastBlockMeanDifference"] and
                          abs(correlation) <= rule["maximumAbsoluteLagOneBlockCorrelation"])
        assessments.append({"observableName": name, "sampleCount": len(trace),
                            "effectiveBlockCount": effective,
                            "firstVsLastBlockMeanDifference": difference,
                            "lagOneBlockCorrelation": correlation,
                            "sufficient": sufficient})
    return assessments, all(item["sufficient"] for item in assessments)


class _FrameSeriesWriter:
    """Stream exact sampled positions; hold only one frame and bounded metadata."""

    def __init__(self, stage_directory: Path, atom_count: int, maximum_bytes: int):
        self.atom_count = atom_count
        self.maximum_bytes = maximum_bytes
        self.path = stage_directory / "sampled-positions.f64"
        self.manifest_path = stage_directory / "sampled-frames.json"
        self.path.touch(exist_ok=False)
        self.digest = hashlib.sha256()
        self.byte_length = 0
        self.frames: list[dict[str, Any]] = []

    def append(self, state: Any, window: str, step: int) -> None:
        import numpy as np
        from openmm import unit

        positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
                               dtype="<f8")
        box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom),
                         dtype=float)
        if (positions.shape != (self.atom_count, 3) or not np.isfinite(positions).all() or
                box.shape != (3, 3) or not np.isfinite(box).all() or
                not math.isfinite(float(np.linalg.det(box))) or np.linalg.det(box) <= 0):
            raise WorkError("nonfiniteObservation", "Sampled positions or periodic cell are unavailable")
        frame = positions.tobytes(order="C")
        if self.byte_length + len(frame) > self.maximum_bytes:
            raise WorkError("resourceRefusal", "Sampled periodic frames exceed their declared byte budget")
        offset = self.byte_length
        # Each bounded frame is closed before the next integration chunk. A
        # stopped or invalid one-request worker cannot leave an open stream.
        with self.path.open("ab") as stream:
            stream.write(frame)
        self.digest.update(frame)
        self.byte_length += len(frame)
        self.frames.append({"windowName": window, "step": step,
                            "offsetBytes": offset, "lengthBytes": len(frame),
                            "positionsSha256": hashlib.sha256(frame).hexdigest(),
                            "boxVectorsAngstrom": box.tolist()})

    def finish(self, payload: dict[str, Any], final_paths: dict[str, Path]) -> None:
        protocol = payload["protocol"]
        manifest = {
            "schemaVersion": "protein-in-membrane.equilibration-frames.v1",
            "encoding": "float64-le-xyz-angstrom",
            "studyRevisionId": payload["studyRevisionId"],
            "attemptId": payload["attemptId"],
            "stageId": payload["stageId"],
            "sourceMinimizedStageId": payload["sourceMinimizedStageId"],
            "validatedPolicyId": payload["validatedPolicyId"],
            "protocolSha256": payload["protocolSha256"],
            "sourceCoordinateSha256": payload["topologyCifSha256"],
            "sourceTopologySha256": payload["topologyJsonSha256"],
            "sourceSystemSha256": payload["systemXmlSha256"],
            "sourceStateSha256": payload["minimizedStateXmlSha256"],
            "finalCoordinateSha256": sha256(final_paths["cif"]),
            "finalTopologySha256": sha256(final_paths["topology"]),
            "finalSystemSha256": sha256(final_paths["system"]),
            "finalStateSha256": sha256(final_paths["state"]),
            "atomCount": self.atom_count,
            "frameCount": len(self.frames),
            "maximumFrameBytes": protocol["maximumFrameBytes"],
            "positionsFile": self.path.name,
            "positionsSha256": self.digest.hexdigest(),
            "positionsByteLength": self.byte_length,
            "frames": self.frames,
        }
        manifest_bytes = json.dumps(manifest, separators=(",", ":"), allow_nan=False).encode("utf-8")
        if len(manifest_bytes) > 64 * 1024 * 1024:
            raise WorkError("resourceRefusal", "Sampled periodic frame manifest exceeds its fixed size cap")
        self.manifest_path.write_bytes(manifest_bytes)


def equilibrate(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm import LangevinMiddleIntegrator, System, XmlSerializer, unit

    stage_id = require_text(payload.get("stageId"), "stageId")
    if re.fullmatch(r"[0-9a-f]{32}", stage_id) is None:
        raise WorkError("invalidRequest", "stageId must be an identified 32-character lowercase hex ID")
    require_text(payload.get("studyRevisionId"), "studyRevisionId")
    require_text(payload.get("attemptId"), "attemptId")
    require_text(payload.get("sourceMinimizedStageId"), "sourceMinimizedStageId")
    protocol_sha256 = require_text(payload.get("protocolSha256"), "protocolSha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", protocol_sha256) is None:
        raise WorkError("invalidRequest", "protocolSha256 must identify the exact declared procedure")
    pdb, original_system, state = _load_stage(directory, payload, "minimizedStateXmlPath")
    validated = require_text(payload.get("validatedPolicyId"), "validatedPolicyId")
    protocol = require_mapping(payload.get("protocol"), "protocol")
    if validated != require_text(protocol.get("id"), "protocol id") or "candidate" in validated.lower():
        raise WorkError("policyNotValidated", "An exact validated, non-candidate policy ID is required")
    stages = protocol.get("stages")
    if not isinstance(stages, list) or not stages:
        raise WorkError("invalidProtocol", "Equilibration requires a nonempty finite stage sequence")
    extension = require_mapping(protocol.get("extensionWindow"), "extensionWindow")
    maximum_extensions = require_integer(protocol.get("maximumExtensions"), "maximumExtensions")
    maximum_samples = require_integer(protocol.get("maximumSampleCount"), "maximumSampleCount", 1)
    observables, rules = _equilibration_observables(protocol, payload, original_system)
    required = protocol.get("requiredObservations")
    if not isinstance(required, list) or not required or any(not isinstance(name, str) for name in required):
        raise WorkError("missingPolicy", "Required observations must be declared")
    if not set(required).issubset({item["name"] for item in observables}):
        raise WorkError("missingPolicy", "Required observations lack declared sampled measurements")
    backbone = _index_selection(payload.get("proteinBackboneAtomIndices"), "proteinBackboneAtomIndices", original_system)
    protein_heavy = _index_selection(payload.get("proteinHeavyAtomIndices", []), "proteinHeavyAtomIndices", original_system)
    lipid = _index_selection(payload.get("lipidHeavyAtomIndices"), "lipidHeavyAtomIndices", original_system)
    topology_atoms = list(pdb.topology.atoms())
    if (set(backbone) & set(lipid) or set(protein_heavy) & set(lipid) or
            any(topology_atoms[index].element.symbol == "H" for index in protein_heavy)):
        raise WorkError("invalidProtocol", "Restraint target groups overlap or include a hydrogen as heavy")
    random_seed = require_integer(protocol.get("randomSeed"), "protocol randomSeed", 1)
    all_controls = [require_mapping(raw, "equilibration stage") for raw in stages]
    extension_name = require_text(extension.get("name"), "extension name")
    stage_names = [require_text(control.get("name"), "stage name") for control in all_controls]
    if (len(extension_name) > 128 or any(len(name) > 128 for name in stage_names) or
            len(set(stage_names)) != len(stage_names) or any(
            name == extension_name or name.startswith(extension_name + "-") for name in stage_names)):
        raise WorkError("invalidProtocol", "Equilibration windows need distinct, noncolliding names")
    if not all_controls[-1].get("proteinRestraintKjMolNm2") == 0 or not all_controls[-1].get("lipidRestraintKjMolNm2") == 0:
        raise WorkError("invalidProtocol", "Final observation stage must be unrestrained")
    if extension.get("proteinRestraintKjMolNm2") != 0 or extension.get("lipidRestraintKjMolNm2") != 0:
        raise WorkError("invalidProtocol", "Observation extensions must be unrestrained")
    if all_controls[-1].get("initialTemperatureKelvin") is not None or extension.get("initialTemperatureKelvin") is not None:
        raise WorkError("invalidProtocol", "Final observation and extension windows must use a fixed temperature")
    for field in ("timestepPicoseconds", "temperatureKelvin", "pressureBar", "pressureMode",
                  "barostatFrequencySteps", "surfaceTensionBarNm", "frictionPerPicosecond"):
        if extension.get(field) != all_controls[-1].get(field):
            raise WorkError("invalidProtocol", f"Observation extensions change final-window condition {field}")
    for control in [*all_controls, extension]:
        _initial_temperature(control)
        if control.get("proteinRestraintSelector", "backbone") not in (
                "backbone", "protein-heavy", "backbone-heavy", "protein-ca"):
            raise WorkError("invalidProtocol", "Unknown protein restraint selector")
        _barostat(System(), control)
    planned_samples = sum(math.ceil(require_integer(control.get("steps"), "stage steps", 1) /
                                    require_integer(control.get("reportIntervalSteps"), "reportIntervalSteps", 1))
                          for control in all_controls)
    planned_samples += maximum_extensions * math.ceil(
        require_integer(extension.get("steps"), "extension steps", 1) /
        require_integer(extension.get("reportIntervalSteps"), "extension reportIntervalSteps", 1))
    if planned_samples > maximum_samples:
        raise WorkError("resourceRefusal", "Declared observation plan exceeds its maximum sample count")
    frame_budget = require_integer(protocol.get("maximumFrameBytes"), "maximumFrameBytes", 1)
    if maximum_samples > 100_000 or frame_budget > 16 * 1024 * 1024 * 1024 or \
            planned_samples * original_system.getNumParticles() * 24 > frame_budget:
        raise WorkError("resourceRefusal", "Declared periodic frame series exceeds its bounded byte budget")
    stage_directory = directory / f"equilibration-{stage_id}"
    try:
        stage_directory.mkdir()
    except FileExistsError as exc:
        raise WorkError("resultConflict", "This exact optional stage already has an output directory") from exc
    frame_writer = _FrameSeriesWriter(stage_directory, original_system.getNumParticles(), frame_budget)
    previous = state
    windows = []
    samples = []
    adequacy_samples = []
    assessments = []
    adequacy = False
    extensions_performed = 0
    final_system = None
    final_unrestrained = False
    stage_index = 0
    while True:
        if stage_index < len(all_controls):
            control = all_controls[stage_index]
            name = require_text(control.get("name"), "stage name")
        else:
            control = extension
            name = f"{extension_name}-{extensions_performed}"
        steps = require_integer(control.get("steps"), "stage steps", 1)
        report = require_integer(control.get("reportIntervalSteps"), "reportIntervalSteps", 1)
        timestep = require_number(control.get("timestepPicoseconds"), "timestepPicoseconds", 0.000001)
        temperature = require_number(control.get("temperatureKelvin"), "stage temperatureKelvin", 0.000001)
        initial_temperature = _initial_temperature(control)
        friction = require_number(control.get("frictionPerPicosecond"), "frictionPerPicosecond", 0.000001)
        protein_k = require_number(control.get("proteinRestraintKjMolNm2"), "proteinRestraintKjMolNm2", 0)
        lipid_k = require_number(control.get("lipidRestraintKjMolNm2"), "lipidRestraintKjMolNm2", 0)
        system = XmlSerializer.deserialize(XmlSerializer.serialize(original_system))
        if protein_k:
            target = _protein_restraint_indices(control.get("proteinRestraintSelector", "backbone"),
                                                topology_atoms, backbone, protein_heavy)
            _restraint(system, target, protein_k, previous.getPositions(), f"protein-restraint-{name}")
        if lipid_k:
            _restraint(system, lipid, lipid_k, previous.getPositions(), f"lipid-restraint-{name}")
        barostat = _barostat(system, control)
        if barostat is not None:
            barostat.setRandomNumberSeed(random_seed + stage_index * 3 + 2)
        integrator = LangevinMiddleIntegrator((initial_temperature if initial_temperature is not None
                                              else temperature) * unit.kelvin,
                                               friction / unit.picosecond,
                                               timestep * unit.picoseconds)
        integrator.setRandomNumberSeed(random_seed + stage_index * 3)
        from openmm import Context

        context = Context(system, integrator)
        context.setPeriodicBoxVectors(*previous.getPeriodicBoxVectors())
        context.setPositions(previous.getPositions())
        velocities = previous.getVelocities()
        if velocities is not None:
            context.setVelocities(velocities)
        initialized_velocities = stage_index == 0 and (initial_temperature is not None or velocities is None or
            all(sum(float(v) ** 2 for v in velocity.value_in_unit(unit.nanometer / unit.picosecond)) < 1e-20
                for velocity in velocities))
        if initialized_velocities:
            context.setVelocitiesToTemperature((initial_temperature if initial_temperature is not None
                                                else temperature) * unit.kelvin, random_seed + stage_index * 3 + 1)
        completed = 0
        while completed < steps:
            chunk = min(report, steps - completed)
            _advance_temperature_control(integrator, initial_temperature, temperature,
                                         steps, completed, chunk)
            completed += chunk
            sampled = _final_state(context)
            sample = {"windowName": name, "step": completed,
                      "measurements": _observable_values(sampled, system, observables)}
            frame_writer.append(sampled, name, completed)
            samples.append(sample)
            if stage_index >= len(all_controls) - 1:
                adequacy_samples.append(sample)
            progress("equilibrationProgress", {"window": name, "completedSteps": completed,
                                                "requestedSteps": steps})
        observed = _final_state(context)
        measurements = _stage_measurements(observed, system.getNumParticles())
        windows.append({"name": name, "requestedSteps": steps, "completedSteps": completed,
                        "temperatureKelvin": _instant_temperature(observed, system),
                        # This is the provider's configured pressure target,
                        # not an instantaneous pressure measurement.
                        "pressureBar": (barostat.getDefaultPressure().value_in_unit(unit.bar)
                                        if barostat is not None else None),
                        "measurements": measurements,
                        "warnings": []})
        progress("equilibrationSeedsUsed", {"window": name,
                                            "integratorSeed": random_seed + stage_index * 3,
                                            "velocityInitializationSeed": random_seed + stage_index * 3 + 1 if initialized_velocities else None,
                                            "barostatSeed": random_seed + stage_index * 3 + 2 if barostat is not None else None})
        previous = observed
        final_system = system
        final_unrestrained = protein_k == 0 and lipid_k == 0
        del context, integrator
        if stage_index >= len(all_controls) - 1:
            assessments, adequacy = _assess_equilibration_samples(adequacy_samples, rules)
        stage_index += 1
        if stage_index < len(all_controls):
            continue
        if adequacy or extensions_performed >= maximum_extensions:
            break
        extensions_performed += 1
    if final_system is None:
        raise WorkError("unobservedOutput", "No equilibration window completed")
    final_xml = stage_directory / "equilibrated-state.xml"
    system_xml = stage_directory / "equilibrated-system.xml"
    final_cif = stage_directory / "equilibrated-coordinates.cif"
    final_topology = stage_directory / "equilibrated-topology.json"
    final_system.setDefaultPeriodicBoxVectors(*previous.getPeriodicBoxVectors())
    final_xml.write_text(XmlSerializer.serialize(previous), encoding="utf-8")
    system_xml.write_text(XmlSerializer.serialize(final_system), encoding="utf-8")
    _write_stage_coordinates(final_cif, pdb.topology, previous.getPositions(), previous.getPeriodicBoxVectors())
    final_topology.write_text(json.dumps(_topology_data(pdb.topology), separators=(",", ":")),
                              encoding="utf-8")
    frame_writer.finish(payload, {"cif": final_cif, "topology": final_topology,
                                  "system": system_xml, "state": final_xml})
    progress("equilibrationWindowsObserved", {"count": len(windows)})
    return {"artifacts": [artifact(directory, final_xml, "equilibratedStateXml"),
                          artifact(directory, system_xml, "equilibratedSystemXml"),
                          artifact(directory, final_cif, "equilibratedCif"),
                          artifact(directory, final_topology, "equilibratedTopologyJson"),
                          artifact(directory, frame_writer.path, "sampledPositionsF64"),
                          artifact(directory, frame_writer.manifest_path, "sampledFramesManifest")],
            "observations": {"windows": windows, "samples": samples,
                             "observationAssessments": assessments,
                             "extensionsPerformed": extensions_performed,
                             "observationAdequacy": "adequate" if adequacy else "insufficientAtBound",
                             "finalAtomCount": original_system.getNumParticles(),
                             "termination": "completed",
                             "finalObservationUnrestrained": final_unrestrained,
                             "numericalWarnings": []},
            "provider": _provider("OpenMM staged dynamics")}
