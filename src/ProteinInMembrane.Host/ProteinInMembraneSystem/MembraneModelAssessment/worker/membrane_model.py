"""Owner-local membrane representation and exact parameter-template observations.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

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
        heads = item.get("headAtomIndices")
        if not isinstance(heads, list) or (item.get("category") == "lipid" and not heads) or any(
                not isinstance(index, int) or index < 1 or index > atom_count for index in heads):
            match = False
            warnings.append("Head-atom indices do not identify atoms in the coordinate template.")
        for label, field in (("area per molecule", "areaPerMoleculeAngstromSquared"),
                             ("molecular volume", "volumeAngstromCubed")):
            if item.get("category") == "lipid":
                require_number(item.get(field), label, 0.000001)
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
