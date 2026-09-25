"""Owner-local OpenMM minimization mechanics.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import WorkError, artifact, require_integer, require_number
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _provider, _load_stage, _rms_force, _state_context, _final_state, _write_stage_coordinates)


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
    final_force = _rms_force(final, system.getNumParticles())
    if not math.isfinite(final_energy):
        raise WorkError("nonfiniteObservation", "Final minimized potential energy is nonfinite")
    final_xml = directory / "minimized-state.xml"
    final_cif = directory / "minimized-coordinates.cif"
    final_xml.write_text(XmlSerializer.serialize(final), encoding="utf-8")
    _write_stage_coordinates(final_cif, pdb.topology, final.getPositions(), final.getPeriodicBoxVectors())
    unrestrained = not any(force.__class__.__name__ in {"CustomExternalForce", "CustomCentroidBondForce"}
                           for force in (system.getForce(index) for index in range(system.getNumForces())))
    del context, integrator
    progress("minimizationObserved", {"finalRmsForceKjMolNm": final_force})
    return {"artifacts": [artifact(directory, final_xml, "minimizedStateXml"),
                          artifact(directory, final_cif, "minimizedCif")],
            "observations": {"initialPotentialEnergyKjMol": initial_energy,
                             "finalPotentialEnergyKjMol": final_energy,
                             "finalRmsForceKjMolNm": final_force,
                             "iterations": None,
                             "termination": "converged" if final_force <= target else "max_iterations",
                             "finalTreatmentUnrestrained": unrestrained,
                             "finalAtomCount": system.getNumParticles(),
                             "numericalWarnings": [] if final_force <= target else
                             ["The declared final RMS-force target was not observed."]},
            "provider": _provider("OpenMM LocalEnergyMinimizer")}
