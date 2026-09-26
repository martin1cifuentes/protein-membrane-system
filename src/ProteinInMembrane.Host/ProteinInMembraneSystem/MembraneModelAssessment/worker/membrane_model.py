"""Owner-local membrane representation and exact parameter-template observations.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (WorkError, require_integer, require_mapping,
                                        require_number, require_text, verify_sha256, work_path)
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _force_field_files, _provider, _copy_molecule,
    _coordinate_file, _nonbonded_charge)


def _template_bonds(template: Any) -> set[tuple[str, str]]:
    return {tuple(sorted((template.atoms[a].name, template.atoms[b].name))) for a, b in template.bonds}


def _residue_bonds(topology: Any, residue: Any) -> set[tuple[str, str]]:
    return {tuple(sorted((a.name, b.name))) for a, b in topology.bonds()
            if a.residue == residue and b.residue == residue}


def _exact_template_match(force_field: Any, topology: Any) -> tuple[bool, int, list[str]]:
    warnings: list[str] = []
    residues = list(topology.residues())
    try:
        templates = force_field.getMatchingTemplates(topology)
    except ValueError as exc:
        return False, 0, [f"Parameter template match failed: {exc}"]
    if len(templates) != len(residues):
        return False, 0, ["Parameter template count differs from coordinate residue count."]
    parameter_atoms = 0
    matched = True
    for residue, template in zip(residues, templates):
        coordinate_atoms = {(atom.name, atom.element.symbol if atom.element else "") for atom in residue.atoms()}
        reference_atoms = {(atom.name, atom.element.symbol if atom.element else "") for atom in template.atoms}
        parameter_atoms += len(reference_atoms)
        if coordinate_atoms != reference_atoms:
            matched = False
            warnings.append(f"Residue {residue.name} coordinate atom names/elements differ from {template.name}.")
        if _residue_bonds(topology, residue) != _template_bonds(template):
            matched = False
            warnings.append(f"Residue {residue.name} coordinate bonds differ from {template.name}.")
    return matched, parameter_atoms, warnings


def _stereo_geometry(pdb: Any, checks: Any) -> list[str]:
    """Observe declared local configuration, without deciding membrane support."""
    from openmm import unit

    if not isinstance(checks, list) or not checks:
        return ["No identified stereochemical geometry checks accompany this molecular representation."]
    atoms = list(pdb.topology.atoms())
    by_name = {atom.name: atom for atom in atoms}
    if len(by_name) != len(atoms):
        return ["Coordinate atom names are not unique for stereochemical correspondence."]
    neighbors: dict[str, set[str]] = {atom.name: set() for atom in atoms}
    for first, second in pdb.topology.bonds():
        neighbors[first.name].add(second.name)
        neighbors[second.name].add(first.name)

    def position(name: str) -> tuple[float, float, float]:
        values = pdb.positions[by_name[name].index].value_in_unit(unit.angstrom)
        return tuple(float(value) for value in values)

    def subtract(first: tuple[float, ...], second: tuple[float, ...]) -> tuple[float, ...]:
        return tuple(first[index] - second[index] for index in range(3))

    def dot(first: tuple[float, ...], second: tuple[float, ...]) -> float:
        return sum(first[index] * second[index] for index in range(3))

    def cross(first: tuple[float, ...], second: tuple[float, ...]) -> tuple[float, ...]:
        return (first[1] * second[2] - first[2] * second[1],
                first[2] * second[0] - first[0] * second[2],
                first[0] * second[1] - first[1] * second[0])

    warnings: list[str] = []
    for index, raw in enumerate(checks, 1):
        if not isinstance(raw, dict):
            warnings.append(f"Stereochemical check {index} has no identified descriptor.")
            continue
        names = raw.get("atomNames")
        kind = raw.get("kind")
        expected = raw.get("expected")
        if (not isinstance(names, list) or len(names) != 4 or
                any(not isinstance(name, str) or name not in by_name for name in names) or
                len(set(names)) != 4):
            warnings.append(f"Stereochemical check {index} does not identify four distinct coordinate atoms.")
            continue
        first, second, third, fourth = (position(name) for name in names)
        label = ", ".join(names)
        if kind == "tetrahedral" and expected in ("positive", "negative"):
            common_centers = set.intersection(*(neighbors[name] for name in names))
            if len(common_centers) != 1:
                warnings.append(f"Tetrahedral check {label} does not identify one bonded center.")
                continue
            signed_volume = dot(subtract(first, fourth),
                                cross(subtract(second, fourth), subtract(third, fourth)))
            if (not math.isfinite(signed_volume) or
                    (signed_volume <= 1.0 if expected == "positive" else signed_volume >= -1.0)):
                warnings.append(f"Tetrahedral check {label} does not match the declared {expected} handedness.")
        elif kind == "alkene" and expected in ("cis", "trans"):
            if names[1] not in neighbors[names[0]] or names[2] not in neighbors[names[1]] or \
                    names[3] not in neighbors[names[2]]:
                warnings.append(f"Alkene check {label} does not follow one bonded four-atom path.")
                continue
            axis = subtract(third, second)
            axis_length_squared = dot(axis, axis)
            if not math.isfinite(axis_length_squared) or axis_length_squared <= 1e-8:
                warnings.append(f"Alkene check {label} has an unresolved central bond geometry.")
                continue
            left = subtract(first, second)
            right = subtract(fourth, third)
            left_projected = tuple(left[i] - dot(left, axis) / axis_length_squared * axis[i]
                                   for i in range(3))
            right_projected = tuple(right[i] - dot(right, axis) / axis_length_squared * axis[i]
                                    for i in range(3))
            denominator = math.sqrt(dot(left_projected, left_projected) *
                                    dot(right_projected, right_projected))
            cosine = dot(left_projected, right_projected) / denominator if denominator > 1e-8 else math.nan
            if not math.isfinite(cosine) or (cosine <= 0.5 if expected == "cis" else cosine >= -0.5):
                warnings.append(f"Alkene check {label} does not match the declared {expected} configuration.")
        else:
            warnings.append(f"Stereochemical check {index} has an unsupported kind or expected configuration.")
    return warnings


def assess_membrane(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm.app import ForceField, Topology

    representations = payload.get("speciesRepresentations")
    if not isinstance(representations, list) or not representations:
        raise WorkError("invalidRequest", "speciesRepresentations must contain the proposed molecule representations")
    assets = payload.get("forceFieldFiles")
    if not isinstance(assets, list) or not assets:
        raise WorkError("missingPolicy", "forceFieldFiles must identify staged hashed assets")
    paths = _force_field_files(directory, assets)
    asset_by_path = {str(work_path(directory, require_mapping(asset, "force-field asset").get("path"),
                                   "force-field path")): asset for asset in assets}
    coordinates = []
    for item in representations:
        item = require_mapping(item, "molecular representation")
        template = work_path(directory, item.get("templatePath"), "templatePath")
        template_sha = require_text(item.get("templateSha256"), "templateSha256")
        verify_sha256(template, template_sha, "templateSha256")
        selected_asset = asset_by_path.get(str(template))
        if (selected_asset is None or
                require_text(selected_asset.get("sha256"), "force-field sha256").lower() != template_sha.lower() or
                require_text(selected_asset.get("version"), "force-field version") !=
                require_text(item.get("forceFieldVersion"), "species forceFieldVersion") or
                require_text(selected_asset.get("family"), "force-field family") !=
                require_text(item.get("forceFieldFamily"), "species forceFieldFamily")):
            raise WorkError("inputMismatch", "A species has no exact staged force-field asset identity")
        coordinate = work_path(directory, item.get("coordinateTemplatePath"), "coordinateTemplatePath")
        verify_sha256(coordinate, item.get("coordinateTemplateSha256"), "coordinateTemplateSha256")
        coordinates.append(coordinate)
    force_field = ForceField(*paths)
    species = []
    combined = Topology()
    for item, coordinate in zip(representations, coordinates):
        pdb = _coordinate_file(coordinate)
        atom_count = len(list(pdb.topology.atoms()))
        match, parameter_count, warnings = _exact_template_match(force_field, pdb.topology)
        stereo_warnings = _stereo_geometry(pdb, item.get("stereoChecks"))
        if stereo_warnings:
            match = False
            warnings.extend(stereo_warnings)
        heads = item.get("headAtomIndices")
        if not isinstance(heads, list) or (item.get("category") == "lipid" and not heads) or any(
                not isinstance(index, int) or index < 1 or index > atom_count for index in heads):
            match = False
            warnings.append("Head-atom indices do not identify atoms in the coordinate template.")
        if atom_count != require_integer(item.get("atomCount"), "declared atomCount", 1):
            match = False
            warnings.append("Declared atom count differs from coordinate template atom count.")
        if len(list(pdb.topology.residues())) != 1:
            match = False
            warnings.append("A species coordinate template must represent one complete molecule/residue.")
        if match:
            try:
                actual_charge = _nonbonded_charge(force_field.createSystem(pdb.topology))
            except ValueError as exc:
                match = False
                warnings.append(f"Single-species parameter assignment failed: {exc}")
            else:
                declared_charge = require_number(item.get("netChargeElementary"), "netChargeElementary")
                if abs(actual_charge - declared_charge) > 1e-5:
                    match = False
                    warnings.append("Declared net charge differs from the parameterized coordinate template.")
        if match:
            _copy_molecule(combined, pdb.topology)
        species.append({"speciesId": require_text(item.get("speciesId"), "speciesId"),
                        "chemistryId": require_text(item.get("chemistryId"), "chemistryId"),
                        "coordinateAtomCount": atom_count, "parameterAtomCount": parameter_count,
                        "atomIdentityAndBondMatch": match, "warnings": warnings})
    combined_observed: bool | None = None
    combination_warnings: list[str] = []
    known_species = {item["speciesId"] for item in species}
    for side in ("upper", "lower"):
        composition = require_mapping(payload.get(side), f"{side} leaflet composition")
        fractions = composition.get("fractions")
        if not isinstance(fractions, list) or not fractions:
            combination_warnings.append(f"{side} leaflet has no represented lipid fractions.")
            continue
        names = [require_text(item.get("speciesId"), "lipid speciesId") for item in fractions]
        if len(set(names)) != len(names) or any(name not in known_species for name in names):
            combination_warnings.append(f"{side} leaflet references a duplicate or unrepresented lipid species.")
        total = sum(require_number(item.get("fraction"), "lipid fraction", 0) for item in fractions)
        if abs(total - 1.0) > 1e-6:
            combination_warnings.append(f"{side} leaflet fractions sum to {total:.8f} rather than one.")
    if all(item["atomIdentityAndBondMatch"] for item in species):
        try:
            force_field.createSystem(combined)
            combined_observed = True
        except ValueError as exc:
            combined_observed = False
            combination_warnings.append(f"Combined force-field parameterization failed: {exc}")
    else:
        combination_warnings.append("Combination parameterization was not attempted because a species failed exact template correspondence.")
    progress("templateFactsObserved", {"speciesCount": len(species)})
    return {"artifacts": [],
            "observations": {"species": species, "combinationWarnings": combination_warnings,
                             "combinedParameterizationObserved": combined_observed},
            "provider": _provider("OpenMM ForceField")}
