"""Owner-local OpenMM minimization mechanics.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import WorkError, artifact, require_integer, require_number
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _provider, _load_stage, _state_context, _final_state, _write_stage_coordinates)


_TANGENT_FORCE_METHOD = "constraint-tangent-per-particle-v1"


def _constrained_tangent_force(system: Any, state: Any) -> tuple[float, float, float]:
    """Observe physical force on the final HBonds constraint tangent space.

    OpenMM's L-BFGS objective temporarily adds harmonic constraint forces and
    projects the final positions back onto constraints.  Raw final State forces
    therefore include legitimate constraint-normal reactions.  This Euclidean
    projection uses the *actual* System constraint Jacobian; it observes the
    remaining physical force without guessing the provider's private stop code.
    See OpenMM 8.6 ReferenceMinimize.cpp evaluate() and final applyConstraints().
    """
    import numpy as np
    from openmm import unit

    count = system.getNumParticles()
    positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
    physical = np.asarray(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer), dtype=float)
    if count < 1 or positions.shape != (count, 3) or physical.shape != (count, 3) or not (
            np.isfinite(positions).all() and np.isfinite(physical).all()):
        raise WorkError("nonfiniteObservation", "Final positions or forces are absent or nonfinite")
    if system.getNumConstraints() < 1:
        raise WorkError("unsupportedPolicy", "Selected HBonds treatment has no observed constraints")
    for index in range(count):
        mass = system.getParticleMass(index).value_in_unit(unit.dalton)
        if not math.isfinite(mass) or mass <= 0:
            raise WorkError("unsupportedPolicy", "Final constrained-force observation requires positive atom masses")

    parent = list(range(count))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    constraints = []
    maximum_relative_error = 0.0
    seen_pairs = set()
    for ordinal in range(system.getNumConstraints()):
        first, second, target = system.getConstraintParameters(ordinal)
        first, second = int(first), int(second)
        desired = target.value_in_unit(unit.nanometer)
        pair = tuple(sorted((first, second)))
        if (first < 0 or first >= count or second < 0 or second >= count or
                first == second or pair in seen_pairs or
                not math.isfinite(desired) or desired <= 0):
            raise WorkError("unsupportedPolicy", "Invalid or duplicated final System constraint")
        seen_pairs.add(pair)
        delta = positions[first] - positions[second]
        observed = float(np.linalg.norm(delta))
        if not math.isfinite(observed) or observed <= 0:
            raise WorkError("nonfiniteObservation", "A final constrained distance is invalid")
        maximum_relative_error = max(maximum_relative_error, abs(observed - desired) / desired)
        constraints.append((first, second, delta / observed))
        a, b = root(first), root(second)
        if a != b:
            parent[b] = a
    groups: dict[int, dict[str, Any]] = {}
    for first, second, direction in constraints:
        key = root(first)
        group = groups.setdefault(key, {"atoms": set(), "constraints": []})
        group["atoms"].update((first, second))
        group["constraints"].append((first, second, direction))
    tangent = physical.copy()
    for group in groups.values():
        atoms = sorted(group["atoms"])
        rows = group["constraints"]
        if len(atoms) > 16 or len(rows) > 16:
            raise WorkError("unsupportedPolicy", "Constraint component exceeds the bounded HBonds route")
        index = {atom: local for local, atom in enumerate(atoms)}
        jacobian = np.zeros((len(rows), 3 * len(atoms)), dtype=float)
        for row, (first, second, direction) in enumerate(rows):
            jacobian[row, 3*index[first]:3*index[first]+3] = direction
            jacobian[row, 3*index[second]:3*index[second]+3] = -direction
        force = physical[atoms].reshape(-1)
        try:
            normal_coefficients, _, rank, _ = np.linalg.lstsq(jacobian.T, force, rcond=1e-12)
        except np.linalg.LinAlgError as exc:
            raise WorkError("unobservedConvergence", "Constraint-force projection is numerically unavailable") from exc
        if rank != len(rows):
            raise WorkError("unobservedConvergence", "Final constraint Jacobian is rank deficient")
        remainder = force - jacobian.T @ normal_coefficients
        if np.max(np.abs(jacobian @ remainder)) > 1e-6 * max(1.0, float(np.max(np.abs(force)))):
            raise WorkError("unobservedConvergence", "Constraint-force projection residual is unresolved")
        tangent[atoms] = remainder.reshape((len(atoms), 3))
    tangent_rms = math.sqrt(float(np.sum(tangent * tangent)) / count)
    raw_component_rms = math.sqrt(float(np.mean(physical * physical)))
    if not all(math.isfinite(value) and value >= 0 for value in
               (tangent_rms, raw_component_rms, maximum_relative_error)):
        raise WorkError("nonfiniteObservation", "Final constrained-force observation is nonfinite")
    return tangent_rms, raw_component_rms, maximum_relative_error


def minimize(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm import LocalEnergyMinimizer, VerletIntegrator, XmlSerializer, unit

    pdb, system, state = _load_stage(directory, payload, "stateXmlPath")
    maximum = require_integer(payload.get("maxIterations"), "maxIterations", 1)
    target = require_number(payload.get("rmsForceTargetKjMolNm"), "rmsForceTargetKjMolNm", 0.000001)
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = _state_context(system, state, integrator)
    initial = context.getState(getEnergy=True)
    initial_energy = initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(initial_energy):
        raise WorkError("nonfiniteObservation", "Pre-minimization potential energy is nonfinite")
    progress("minimizationStarted", {"maxIterations": maximum})
    LocalEnergyMinimizer.minimize(context, target, maximum)
    final = _final_state(context)
    final_energy = final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    final_force, raw_force, maximum_constraint_error = _constrained_tangent_force(system, final)
    applied_constraint_tolerance = integrator.getConstraintTolerance()
    if not math.isfinite(final_energy):
        raise WorkError("nonfiniteObservation", "Final minimized potential energy is nonfinite")
    final_xml = directory / "minimized-state.xml"
    final_cif = directory / "minimized-coordinates.cif"
    final_xml.write_text(XmlSerializer.serialize(final), encoding="utf-8")
    _write_stage_coordinates(final_cif, pdb.topology, final.getPositions(), final.getPeriodicBoxVectors())
    unrestrained = not any(force.__class__.__name__ in {"CustomExternalForce", "CustomCentroidBondForce"}
                           for force in (system.getForce(index) for index in range(system.getNumForces())))
    del context, integrator
    final_criterion_observed = (maximum_constraint_error <= applied_constraint_tolerance and
                                final_force <= target)
    progress("minimizationObserved", {"finalRmsForceKjMolNm": final_force,
                                      "forceMeasurementMethod": _TANGENT_FORCE_METHOD,
                                      "maximumRelativeConstraintError": maximum_constraint_error})
    return {"artifacts": [artifact(directory, final_xml, "minimizedStateXml"),
                          artifact(directory, final_cif, "minimizedCif")],
            "observations": {"initialPotentialEnergyKjMol": initial_energy,
                             "finalPotentialEnergyKjMol": final_energy,
                             "finalRmsForceKjMolNm": final_force,
                             "iterations": None,
                             "termination": "converged" if final_criterion_observed else "unknown",
                             "finalTreatmentUnrestrained": unrestrained,
                             "finalAtomCount": system.getNumParticles(),
                             "finalRawRmsForceKjMolNm": raw_force,
                             "maximumRelativeConstraintError": maximum_constraint_error,
                             "appliedConstraintTolerance": applied_constraint_tolerance,
                             "finalRmsForceMethod": _TANGENT_FORCE_METHOD,
                             "numericalWarnings": [] if final_criterion_observed else
                             ["The final constraint-tangent RMS-force target or constraint tolerance was not observed; the provider's private stop reason is unknown."]},
            "provider": _provider("OpenMM LocalEnergyMinimizer")}
