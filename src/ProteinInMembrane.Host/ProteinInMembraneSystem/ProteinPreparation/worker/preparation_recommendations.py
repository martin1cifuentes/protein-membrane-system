"""Bounded, disposable protein starting-state plan using the selected OpenMM method.

This operation observes a jointly parameterizable candidate.  It does not
authorize repairs, promote a prepared protein, or infer an optimal protonation
state.  The host binds the returned plan to its study and applies it only after
one explicit actor authorization.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.ProteinPreparation.worker.protein_preparation import (
    _ALLOWED_VARIANTS,
    _address_from_key,
    _apply_altlocs,
    _check_membership,
    _establish_approved_disulfides,
    _fix_heavy_atoms,
    _icode,
    _normalize_address,
    _open_structure,
    _approved_sulfur_pairs,
    _selected_structure,
    prepare_protein,
)
from ProteinInMembraneSystem.worker.exchange import (
    WorkError,
    require_integer,
    require_mapping,
    require_number,
    require_text,
    sha256,
    verify_sha256,
    work_path,
)
from ProteinInMembraneSystem.worker.parameterized_structure import _force_field_files


def _resolved_variant(residue: Any, provider_value: str | None) -> str | None:
    """Resolve OpenMM's None for ordinary states from actual bonded hydrogen names."""
    name = residue.name
    if name not in _ALLOWED_VARIANTS:
        if provider_value is not None:
            raise WorkError("unsupportedChemistry", "The provider chose a variant outside the selected catalogue")
        return None
    hydrogens = {atom.name for atom in residue.atoms() if atom.element.symbol == "H"}
    if name == "ASP":
        resolved = "ASH" if "HD2" in hydrogens else "ASP"
    elif name == "GLU":
        resolved = "GLH" if "HE2" in hydrogens else "GLU"
    elif name == "LYS":
        resolved = "LYS" if "HZ3" in hydrogens else "LYN"
    elif name == "CYS":
        resolved = "CYS" if "HG" in hydrogens else "CYX"
    else:
        nd1, ne2 = "HD1" in hydrogens, "HE2" in hydrogens
        resolved = "HIP" if nd1 and ne2 else "HID" if nd1 else "HIE" if ne2 else "HIN"
    if resolved not in _ALLOWED_VARIANTS[name] or provider_value not in (None, resolved):
        raise WorkError("unsupportedChemistry", "The provider's named state and actual hydrogen pattern disagree",
                        {"residue": residue.id, "providerVariant": provider_value, "resolvedVariant": resolved})
    return resolved


def recommend_preparation(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm.app import ForceField, Modeller

    source_path = work_path(directory, payload.get("sourcePath"), "sourcePath")
    source_digest = verify_sha256(source_path, payload.get("sourceSha256"), "sourceSha256")
    model_index = require_integer(payload.get("modelIndex"), "modelIndex")
    nominal_ph = require_number(payload.get("nominalPh"), "nominalPh", 0)
    if nominal_ph > 14:
        raise WorkError("invalidRequest", "nominalPh exceeds the bounded OpenMM range")
    seed = require_integer(payload.get("seed"), "seed")
    if seed > 2**32 - 1:
        raise WorkError("invalidRequest", "seed exceeds the bounded 32-bit range")
    import random
    import numpy as np

    random.seed(seed)
    np.random.seed(seed)
    allowed = {require_text(value, "permitted variant") for value in payload.get("permittedVariants", [])}
    if not allowed:
        raise WorkError("invalidRequest", "The identified chemical-state variants are absent")
    ff_files = payload.get("forceFieldFiles")
    ff_paths = _force_field_files(directory, ff_files)
    source = _open_structure(source_path)
    selected, chain_map = _selected_structure(source, payload)
    _apply_altlocs(selected, payload, chain_map)
    _check_membership(selected, payload, chain_map)
    import gemmi

    removed_source_hydrogens = []
    for chain in selected[0]:
        for residue in chain:
            if residue.entity_type != gemmi.EntityType.Polymer:
                continue
            for index in range(len(residue) - 1, -1, -1):
                atom = residue[index]
                if atom.element.name in {"H", "D"}:
                    key = (chain.name, int(residue.seqid.num), _icode(residue.seqid.icode), atom.name)
                    removed_source_hydrogens.append({"residue": _address_from_key(key[:3], chain_map, model_index),
                                                     "atomName": key[3]})
                    del residue[index]
    selected_path = directory / "recommendation-selected-no-protein-h.pdb"
    selected.write_pdb(str(selected_path))
    progress("recommendationSelected", {"removedSourceHydrogenCount": len(removed_source_hydrogens)})

    approvals = []
    for raw in payload.get("proposedHeavyAtoms", []):
        item = require_mapping(raw, "proposed heavy atom")
        approvals.append({"residue": item.get("residue"),
                          "atomName": require_text(item.get("atomName"), "atomName"),
                          "decisionId": "recommendation-plan-only"})
    trial = dict(payload, approvedHeavyAtoms=approvals)
    fixer, added_heavy, _ = _fix_heavy_atoms(selected_path, trial, chain_map)
    approved_disulfides = _approved_sulfur_pairs(payload, chain_map)
    _establish_approved_disulfides(fixer.topology, fixer.positions, approved_disulfides,
                                   require_number(payload.get("disulfideCandidateMaxSgDistanceAngstrom"),
                                                  "disulfideCandidateMaxSgDistanceAngstrom", 0.000001))
    overrides = {}
    for raw in payload.get("overrides", []):
        item = require_mapping(raw, "state override")
        key = _normalize_address(item.get("residue"), chain_map, model_index)
        variant = require_text(item.get("variant"), "override variant")
        if key in overrides or variant not in allowed:
            raise WorkError("invalidSelection", "An override is repeated or outside the identified catalogue")
        overrides[key] = variant
    requested = []
    seen = set()
    for residue in fixer.topology.residues():
        key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
        variant = overrides.get(key)
        if variant is not None:
            if residue.name not in _ALLOWED_VARIANTS or variant not in _ALLOWED_VARIANTS[residue.name]:
                raise WorkError("invalidSelection", "The override is not an option for its exact residue")
            seen.add(key)
        requested.append(variant)
    if seen != set(overrides):
        raise WorkError("invalidSelection", "An override does not match the selected molecular construct")

    random.seed(seed)
    np.random.seed(seed)
    modeller = Modeller(fixer.topology, fixer.positions)
    force_field = ForceField(*ff_paths)
    returned = modeller.addHydrogens(force_field, pH=nominal_ph, variants=requested)
    if len(returned) != len(list(modeller.topology.residues())):
        raise WorkError("providerMismatch", "The recommendation method did not account for every residue")
    choices = []
    for residue, provider_variant in zip(modeller.topology.residues(), returned):
        variant = _resolved_variant(residue, provider_variant)
        if variant is None:
            continue
        if variant not in allowed:
            raise WorkError("unsupportedChemistry", "A returned variant is outside the selected catalogue",
                            {"residue": residue.id, "variant": variant})
        key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
        choices.append({"residue": _address_from_key(key, chain_map, model_index),
                        "variant": variant, "overridden": key in overrides})
    if len(choices) != len({json.dumps(item["residue"], sort_keys=True) for item in choices}):
        raise WorkError("providerMismatch", "The method repeated a variable residue address")
    # Full-system template matching checks the whole proposed combination; no
    # sitewise list or successful addHydrogens call is treated as joint support.
    force_field.createSystem(modeller.topology)
    progress("recommendationJointCheck", {"stateChoiceCount": len(choices),
                                           "heavyAtomCount": len(added_heavy)})

    explicit = [{"residue": item["residue"], "variant": item["variant"],
                 "decisionId": "recommendation-plan-only"} for item in choices]
    final_trial = dict(trial, residueVariants=explicit, normalizeProteinHydrogens=True, planSeed=seed)
    random.seed(seed)
    np.random.seed(seed)
    checked = prepare_protein(directory, final_trial, progress)
    observed = checked["observations"]
    expected = {(json.dumps(item["residue"], sort_keys=True), item["variant"]) for item in choices}
    actual = {(json.dumps(item["residue"], sort_keys=True), item["variant"])
              for item in observed["actualResidueVariants"]}
    if actual != expected or len(actual) != len(choices):
        raise WorkError("providerMismatch", "Joint candidate states differ from method-resolved plan")
    if sha256(source_path) != source_digest:
        raise WorkError("inputMismatch", "The source changed during recommendation")
    candidate = next(item for item in checked["artifacts"] if item["role"] == "preparedPdb")
    identity = {
        "sourceSha256": source_digest, "modelIndex": model_index,
        "assemblyId": payload.get("assemblyId"), "chainSelections": payload.get("chainSelections"),
        "altlocChoices": payload.get("altlocChoices"), "retainedPartners": payload.get("retainedPartners"),
        "approvedDisulfides": payload.get("approvedDisulfides"),
        "nominalPh": nominal_ph, "seed": seed,
        "forceFieldFiles": [{"sha256": require_text(item.get("sha256"), "force-field digest")}
                            for item in ff_files],
        "proposedHeavyAtoms": sorted(payload.get("proposedHeavyAtoms", []),
                                     key=lambda item: json.dumps(item, sort_keys=True)),
        "choices": choices, "candidateSha256": candidate["sha256"],
    }
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"artifacts": checked["artifacts"],
            "observations": {"planSha256": digest, "candidateSha256": candidate["sha256"],
                             "choices": choices,
                             "proposedHeavyAtoms": [{"residue": _address_from_key(key[:3], chain_map, model_index),
                                                     "atomName": key[3]} for key in sorted(added_heavy)],
                             "removedSourceHydrogens": sorted(removed_source_hydrogens,
                                                              key=lambda item: json.dumps(item, sort_keys=True)),
                             "candidateAtomCount": observed["preparedAtomCount"],
                             "candidateObservations": observed,
                             "method": "OpenMM Modeller.addHydrogens starting-state rule",
                             "methodVersion": importlib.metadata.version("openmm"),
                             "nominalPh": nominal_ph, "seed": seed,
                             "jointParameterizationObserved": True},
            "provider": {"name": "PDBFixer/OpenMM starting-state recommendation",
                         "version": "/".join(importlib.metadata.version(name)
                                             for name in ("pdbfixer", "openmm"))}}
