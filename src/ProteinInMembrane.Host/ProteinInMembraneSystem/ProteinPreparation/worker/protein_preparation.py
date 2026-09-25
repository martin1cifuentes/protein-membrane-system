"""Structural reading and bounded candidate protein mechanics.

Gemmi interprets source models and selected assembly copies. PDBFixer is used
only for approved non-backbone heavy atoms. Modeller places hydrogens only
after explicit policy variants have been supplied for variable residues.
Neither a generated coordinate file nor parameter matching is a product
assessment of protein suitability.
"""

from __future__ import annotations

import json
import importlib.metadata
import math
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import WorkError, artifact, require_integer, require_mapping, require_number, require_text, sha256, verify_sha256, work_path
from ProteinInMembraneSystem.worker.protein_identity import _icode, _address_from_key, _atom_address

_CANONICAL = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
_BACKBONE = frozenset({"N", "CA", "C", "O"})
_VARIABLE = {"ASP", "GLU", "CYS", "HIS", "LYS"}
_ALLOWED_VARIANTS = {
    "ASP": {"ASP", "ASH"},
    "GLU": {"GLU", "GLH"},
    "CYS": {"CYS", "CYX"},
    "HIS": {"HID", "HIE", "HIP"},
    "LYS": {"LYS", "LYN"},
}


def _address(chain: str, residue: Any, model: int = 0, copy_id: str = "") -> dict[str, Any]:
    return {
        "model": model,
        "chain": chain,
        "residue": int(residue.seqid.num),
        "insertionCode": _icode(residue.seqid.icode),
        "copyId": copy_id,
    }


def _key(address: dict[str, Any]) -> tuple[str, int, str]:
    return (
        require_text(address.get("chain"), "residue address chain"),
        require_integer(address.get("residue"), "residue address number", -999999),
        _icode(address.get("insertionCode")),
    )


def _open_structure(path: Path):
    import gemmi

    structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Detect)
    if len(structure) == 0:
        raise WorkError("invalidStructure", "Structural source contains no coordinate model")
    structure.setup_entities()
    return structure


def _assembly_copy_origin(copy_name: str, source_names: list[str]) -> str:
    """Find the unique source chain behind a full Gemmi output chain ID."""
    matches = [name for name in source_names
               if copy_name.startswith(name) and copy_name[len(name):].isdigit()]
    if len(matches) != 1:
        raise WorkError("ambiguousAssembly", "Assembly-generated chain name does not identify one source chain copy",
                        {"copy": copy_name, "candidateSourceChains": matches})
    return matches[0]


def _residue_kind(residue: Any, peptide_subchains: set[str]) -> str:
    """Classify by Gemmi entity semantics, never by ATOM/HETATM alone."""
    import gemmi

    if residue.is_water() or residue.entity_type == gemmi.EntityType.Water:
        return "solvent"
    if residue.entity_type == gemmi.EntityType.NonPolymer:
        return "heterogen"
    if residue.entity_type == gemmi.EntityType.Polymer:
        if residue.subchain in peptide_subchains or gemmi.find_tabulated_residue(residue.name).is_amino_acid():
            return "protein"
    return "unknown"


def inspect_source(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    import gemmi

    path = work_path(directory, payload.get("sourcePath"), "sourcePath")
    verify_sha256(path, payload.get("sourceSha256"), "sourceSha256")
    max_atoms = require_integer(payload.get("maxAtoms", 2_000_000), "maxAtoms", 1)
    structure = _open_structure(path)
    peptide_subchains = {subchain for entity in structure.entities
                         if entity.entity_type == gemmi.EntityType.Polymer and
                         gemmi.sequence_kind(entity.polymer_type) == gemmi.ResidueKind.AA
                         for subchain in entity.subchains}
    model_summaries: list[dict[str, Any]] = []
    total_atoms = 0
    for index, model in enumerate(structure):
        chains = []
        partners = []
        residues = []
        for chain in model:
            residue_count = len(chain)
            atom_count = sum(len(residue) for residue in chain)
            total_atoms += atom_count
            chains.append({
                "name": chain.name,
                "residueCount": residue_count,
                "atomCount": atom_count,
                "residues": [
                    {"address": _address(chain.name, residue, index), "name": residue.name,
                     "alternateLocations": sorted({str(atom.altloc) for atom in residue if atom.has_altloc()})}
                    for residue in chain
                ],
            })
            for residue in chain:
                residue_atoms = {atom.name for atom in residue if atom.element.name not in {"H", "D"}}
                address = _address(chain.name, residue, index)
                residue_kind = _residue_kind(residue, peptide_subchains)
                residues.append({"address": _address(chain.name, residue, index), "name": residue.name,
                                 "residueKind": residue_kind,
                                 "backboneHeavyAtomsComplete": _BACKBONE.issubset(residue_atoms),
                                 "alternateLocations": sorted({str(atom.altloc) for atom in residue if atom.has_altloc()}),
                                 "missingNonbackboneHeavyAtomNames": [],
                                 "possibleDisulfidePartner": None,
                                 "missingAtomAssessmentStanding": "Unavailable",
                                 "disulfideAssessmentStanding": "Unavailable",
                                 "limitations": []})
                if residue_kind == "heterogen":
                    partners.append({"sourceId": f"{index}:{chain.name}:{residue.seqid.num}:{_icode(residue.seqid.icode)}:{residue.name}",
                                     "label": residue.name, "kind": "nonpolymer",
                                     "atomCount": len(residue), "chain": chain.name,
                                     "residue": int(residue.seqid.num)})
        assemblies = []
        for assembly in structure.assemblies:
            expanded = gemmi.make_assembly(assembly, model, gemmi.HowToNameCopiedChain.AddNumber)
            copies = []
            source_names = [chain.name for chain in model]
            for copy in expanded:
                origin = _assembly_copy_origin(copy.name, source_names)
                copies.append({"sourceChain": origin, "copyId": copy.name})
            assemblies.append({"name": assembly.name, "chainCopies": copies})
        model_summaries.append({"index": index, "chains": chains,
                                "assemblies": assemblies, "partners": partners, "residues": residues,
                                "atomCount": sum(c["atomCount"] for c in chains)})
    if total_atoms > max_atoms:
        raise WorkError("resourceRefused", "Structural source exceeds the declared atom inspection limit", {"atomCount": total_atoms})
    source_format = "mmcif" if path.name.lower().endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")) else (
        "pdb" if path.name.lower().endswith((".pdb", ".pdb.gz", ".ent", ".ent.gz")) else "coordinate")
    from .prediction_source import inspect_prediction
    prediction, prediction_artifacts = inspect_prediction(directory, payload, structure,
                                                          sha256(path), peptide_subchains)
    return {
        "artifacts": prediction_artifacts,
        "observations": {
            "sourceFormat": source_format,
            "models": model_summaries,
            "prediction": prediction,
        },
        "provider": {"name": "Gemmi", "version": __import__("gemmi").__version__},
    }


def _selected_structure(source: Any, payload: dict[str, Any]):
    import gemmi

    model_index = require_integer(payload.get("modelIndex"), "modelIndex")
    if model_index >= len(source):
        raise WorkError("invalidSelection", "Selected coordinate model is absent")
    source_model = source[model_index]
    assembly_id = payload.get("assemblyId")
    if assembly_id:
        assembly = next((item for item in source.assemblies if item.name == assembly_id), None)
        if assembly is None:
            raise WorkError("invalidSelection", "Selected biological assembly is absent")
        model = gemmi.make_assembly(assembly, source_model, gemmi.HowToNameCopiedChain.AddNumber)
    else:
        model = source_model
    raw_selected = payload.get("chainSelections")
    if not isinstance(raw_selected, list) or not raw_selected:
        raise WorkError("invalidSelection", "chainSelections must identify at least one exact chain instance")
    selected = []
    for raw in raw_selected:
        item = require_mapping(raw, "chain selection")
        source_chain = require_text(item.get("sourceChain"), "sourceChain")
        copy_id = require_text(item.get("copyId"), "copyId")
        if assembly_id and _assembly_copy_origin(copy_id, [chain.name for chain in source_model]) != source_chain:
            raise WorkError("invalidSelection", "An assembly copy ID must be its exact output chain name")
        if not assembly_id and copy_id != source_chain:
            raise WorkError("invalidSelection", "A direct source chain's copy identity must equal its chain name")
        selected_name = copy_id if assembly_id else source_chain
        selected.append((selected_name, source_chain, copy_id))
    if len({name for name, _, _ in selected}) != len(selected):
        raise WorkError("invalidSelection", "chainSelections contains duplicate or colliding instances")
    available = {chain.name for chain in model}
    if not {name for name, _, _ in selected}.issubset(available):
        raise WorkError("invalidSelection", "A selected chain instance is absent", {"availableChains": sorted(available)})
    for name, _, _ in selected:
        if sum(chain.name == name for chain in model) != 1:
            raise WorkError("ambiguousAssembly", "Selected chain instance does not identify exactly one structure chain",
                            {"selectedChain": name})

    chosen = gemmi.Structure()
    chosen.name = source.name
    chosen.cell = source.cell
    selected_model = gemmi.Model("1")
    chain_map = []
    names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(selected) > len(names):
        raise WorkError("unsupportedRepresentation", "Selected chain count exceeds the bounded PDB exchange representation")
    for short_name, (original_name, source_chain, copy_id) in zip(names, selected):
        copied = selected_model.add_chain(model[original_name])
        copied.name = short_name
        if assembly_id:
            observed_origin = _assembly_copy_origin(original_name, [chain.name for chain in source_model])
            if observed_origin != source_chain or original_name != copy_id:
                raise WorkError("invalidSelection", "Requested assembly copy conflicts with inspected source identity")
        chain_map.append({"preparedChain": short_name, "selectedChain": original_name,
                          "sourceChain": source_chain, "copyId": copy_id})
    chosen.add_model(selected_model)
    return chosen, chain_map


def _normalize_address(address: Any, chain_map: list[dict[str, str]], model_index: int) -> tuple[str, int, str]:
    item = require_mapping(address, "residue address")
    if require_integer(item.get("model"), "residue address model") != model_index:
        raise WorkError("invalidSelection", "A residue decision addresses a different source model")
    name = require_text(item.get("chain"), "residue address chain")
    copy_id = str(item.get("copyId") or "")
    matches = [entry["preparedChain"] for entry in chain_map
               if (name == entry["selectedChain"] and (not copy_id or copy_id == entry["copyId"]))
               or (name == entry["sourceChain"] and copy_id == entry["copyId"])]
    if len(matches) != 1:
        raise WorkError("invalidSelection", "Residue decision does not identify exactly one selected chain copy")
    return matches[0], require_integer(item.get("residue"), "residue address number", -999999), _icode(item.get("insertionCode"))


def _apply_altlocs(structure: Any, payload: dict[str, Any], chain_map: list[dict[str, str]]) -> None:
    choices: dict[tuple[str, int, str], str] = {}
    for entry in payload.get("altlocChoices", []):
        item = require_mapping(entry, "altloc choice")
        key = _normalize_address(item.get("residue"), chain_map, payload["modelIndex"])
        if key in choices:
            raise WorkError("invalidSelection", "A residue has conflicting alternate-location choices")
        choices[key] = require_text(item.get("altloc"), "altloc")
    seen = set()
    for chain in structure[0]:
        for residue in chain:
            key = _key(_address(chain.name, residue))
            options = {str(atom.altloc) for atom in residue if atom.has_altloc()}
            if not options:
                continue
            choice = choices.get(key)
            if choice not in options:
                raise WorkError("approvalRequired", "Alternate conformation lacks an exact approved choice", {"residue": _address(chain.name, residue), "alternatives": sorted(options)})
            seen.add(key)
            for index in range(len(residue) - 1, -1, -1):
                atom = residue[index]
                if atom.has_altloc() and str(atom.altloc) != choice:
                    del residue[index]
            for atom in residue:
                if atom.has_altloc():
                    atom.altloc = "\0"
    if set(choices) != seen:
        raise WorkError("invalidSelection", "An alternate-location choice does not address an affected selected residue")


def _check_membership(structure: Any, payload: dict[str, Any], chain_map: list[dict[str, str]]) -> None:
    import gemmi

    retained = {}
    for entry in payload.get("retainedPartners", []):
        item = require_mapping(entry, "retained partner")
        identifier = require_text(item.get("sourceId"), "partner sourceId")
        parts = identifier.split(":", 4)
        if len(parts) != 5 or require_integer(int(parts[0]), "partner model") != payload["modelIndex"]:
            raise WorkError("invalidSelection", "A retained partner does not identify the chosen source model")
        chain_name = require_text(item.get("chain"), "partner chain")
        residue_number = require_integer(item.get("residue"), "partner residue", -999999)
        if chain_name != parts[1] or residue_number != int(parts[2]) or item.get("label") != parts[4]:
            raise WorkError("invalidSelection", "Retained partner fields disagree with its source identity")
        retained[identifier] = (chain_name, residue_number, parts[3], parts[4])
    seen_retained = set()
    by_prepared_chain = {entry["preparedChain"]: entry for entry in chain_map}
    for chain in structure[0]:
        for index in range(len(chain) - 1, -1, -1):
            residue = chain[index]
            address = _address(chain.name, residue)
            if residue.is_water() or residue.entity_type == gemmi.EntityType.Water:
                del chain[index]
                continue
            if residue.entity_type == gemmi.EntityType.Polymer:
                if residue.name not in _CANONICAL:
                    raise WorkError("unsupportedChemistry", "A selected protein residue is noncanonical", {"residue": address, "name": residue.name})
                atoms = {atom.name for atom in residue if atom.element.name != "H"}
                if not _BACKBONE.issubset(atoms):
                    raise WorkError("unsupportedBackbone", "A selected retained backbone is incomplete", {"residue": address, "missing": sorted(_BACKBONE - atoms)})
                # Normalize the selected canonical polymer to PDB ATOM records
                # for the bounded PDBFixer exchange, preserving source identity
                # separately in the correspondence.
                residue.het_flag = "A"
                continue
            if residue.entity_type != gemmi.EntityType.NonPolymer:
                raise WorkError("unresolvedStructure", "A selected residue's molecular role cannot be established",
                                {"residue": address, "name": residue.name})
            source_chain = by_prepared_chain[chain.name]["sourceChain"]
            matches = [identifier for identifier, parts in retained.items()
                       if parts == (source_chain, int(residue.seqid.num), _icode(residue.seqid.icode), residue.name)]
            if matches:
                seen_retained.update(matches)
                residue.het_flag = "H"
            else:
                del chain[index]
    if set(retained) != seen_retained:
        raise WorkError("invalidSelection", "A retained partner does not match the selected structure")


def _approved_atom_keys(payload: dict[str, Any], chain_map: list[dict[str, str]]) -> dict[tuple[str, int, str, str], str]:
    approved = {}
    for entry in payload.get("approvedHeavyAtoms", []):
        item = require_mapping(entry, "approved heavy atom")
        key = _normalize_address(item.get("residue"), chain_map, payload["modelIndex"])
        name = require_text(item.get("atomName"), "approved atom name")
        full = (*key, name)
        if full in approved:
            raise WorkError("invalidSelection", "Duplicate approved atom decision")
        approved[full] = require_text(item.get("decisionId"), "decisionId")
    return approved


def _fix_heavy_atoms(pdb_path: Path, payload: dict[str, Any], chain_map: list[dict[str, str]]):
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb_path))
    # No missing-residue, mutation, or nonstandard-residue operation is called.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    approved = _approved_atom_keys(payload, chain_map)
    needed = set()
    for mapping in (fixer.missingAtoms, fixer.missingTerminals):
        for residue, names in mapping.items():
            address = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
            for missing_atom in names:
                name = require_text(getattr(missing_atom, "name", missing_atom), "missing atom name")
                if name in _BACKBONE:
                    raise WorkError("unsupportedBackbone", "Heavy-atom completion would create a missing retained backbone atom", {"address": address, "atom": name})
                needed.add((*address, name))
    if needed != set(approved):
        raise WorkError("approvalRequired", "Approved heavy-atom changes do not equal the observed missing atoms",
                        {"needed": sorted(needed), "approved": sorted(approved)})
    if needed:
        fixer.addMissingAtoms()
    return fixer, needed, approved


def _explicit_variants(topology: Any, payload: dict[str, Any], chain_map: list[dict[str, str]]) -> list[str | None]:
    selections = {}
    for entry in payload.get("residueVariants", []):
        item = require_mapping(entry, "residue variant")
        key = _normalize_address(item.get("residue"), chain_map, payload["modelIndex"])
        if key in selections:
            raise WorkError("invalidSelection", "A residue has conflicting selected variants")
        selections[key] = require_text(item.get("variant"), "variant")
    output: list[str | None] = []
    seen = set()
    for residue in topology.residues():
        key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
        if residue.name in _VARIABLE:
            variant = selections.get(key)
            if variant not in _ALLOWED_VARIANTS[residue.name]:
                raise WorkError("chemicalStateUnresolved", "A chemically variable residue lacks a supported explicit variant",
                                {"chain": key[0], "residue": key[1], "name": residue.name})
            output.append(variant)
            seen.add(key)
        else:
            if key in selections:
                raise WorkError("unsupportedChemistry", "An unsupported residue-state override was requested",
                                {"chain": key[0], "residue": key[1], "name": residue.name})
            output.append(None)
    if set(selections) != seen:
        raise WorkError("invalidSelection", "A selected chemical state does not match the prepared topology")
    return output


def _sulfur_pair_keys(topology: Any) -> set[tuple[tuple[str, int, str], tuple[str, int, str]]]:
    pairs = set()
    for first, second in topology.bonds():
        if first.name != "SG" or second.name != "SG":
            continue
        if first.residue.name not in {"CYS", "CYX"} or second.residue.name not in {"CYS", "CYX"}:
            continue
        a = (first.residue.chain.id, int(first.residue.id), _icode(first.residue.insertionCode))
        b = (second.residue.chain.id, int(second.residue.id), _icode(second.residue.insertionCode))
        if a == b:
            raise WorkError("providerMismatch", "A disulfide bond joins sulfur atoms in one residue")
        pairs.add(tuple(sorted((a, b))))
    return pairs


def _approved_sulfur_pairs(payload: dict[str, Any], chain_map: list[dict[str, str]]) -> set[tuple[tuple[str, int, str], tuple[str, int, str]]]:
    pairs = set()
    occupied = set()
    for raw in payload.get("approvedDisulfides", []):
        choice = require_mapping(raw, "approved disulfide")
        require_text(choice.get("decisionId"), "disulfide decisionId")
        a = _normalize_address(choice.get("first"), chain_map, payload["modelIndex"])
        b = _normalize_address(choice.get("second"), chain_map, payload["modelIndex"])
        if a == b or a in occupied or b in occupied:
            raise WorkError("invalidSelection", "Approved disulfide pairs overlap or self-link")
        occupied.update((a, b))
        pairs.add(tuple(sorted((a, b))))
    return pairs


def _establish_approved_disulfides(topology: Any, positions: Any,
                                    approved: set[tuple[tuple[str, int, str], tuple[str, int, str]]],
                                    maximum_distance_angstrom: float) -> None:
    from openmm import unit

    actual = _sulfur_pair_keys(topology)
    if actual - approved:
        raise WorkError("approvalRequired", "The provider inferred a disulfide that was not approved",
                        {"unexpectedPairs": [list(pair) for pair in sorted(actual - approved)]})
    sites = {}
    for atom in topology.atoms():
        if atom.name == "SG" and atom.residue.name in {"CYS", "CYX"}:
            key = (atom.residue.chain.id, int(atom.residue.id), _icode(atom.residue.insertionCode))
            if key in sites:
                raise WorkError("correspondenceFailed", "A cysteine has multiple selected sulfur atoms")
            sites[key] = atom
    xyz = positions.value_in_unit(unit.angstrom)
    for a, b in approved - actual:
        if a not in sites or b not in sites:
            raise WorkError("invalidSelection", "An approved disulfide lacks an exact sulfur atom")
        distance = math.dist(xyz[sites[a].index], xyz[sites[b].index])
        if not math.isfinite(distance) or distance > maximum_distance_angstrom:
            raise WorkError("invalidSelection", "Approved disulfide pair no longer matches the inspected sulfur proximity",
                            {"distanceAngstrom": distance})
        topology.addBond(sites[a], sites[b])
    if _sulfur_pair_keys(topology) != approved:
        raise WorkError("providerMismatch", "Prepared disulfide bonds differ from the approved pair set")


def inspect_preparation_changes(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from pdbfixer import PDBFixer
    from ProteinInMembraneSystem.worker.geometry_observations import observe_protein_geometry, unavailable_geometry

    source_path = work_path(directory, payload.get("sourcePath"), "sourcePath")
    verify_sha256(source_path, payload.get("sourceSha256"), "sourceSha256")
    source = _open_structure(source_path)
    selected, chain_map = _selected_structure(source, payload)
    _apply_altlocs(selected, payload, chain_map)
    _check_membership(selected, payload, chain_map)
    model_index = require_integer(payload.get("modelIndex"), "modelIndex")
    disulfide_distance = require_number(payload.get("disulfideCandidateMaxSgDistanceAngstrom"),
                                        "disulfideCandidateMaxSgDistanceAngstrom", 0.000001)
    selected_path = directory / "preparation-change-selection.pdb"
    selected.write_pdb(str(selected_path))
    protein_only = directory / "preparation-change-protein-only.pdb"
    with protein_only.open("w", encoding="utf-8") as stream:
        for line in selected_path.read_text(encoding="utf-8").splitlines():
            if line.startswith(("ATOM  ", "TER   ", "CRYST1")):
                stream.write(line + "\n")
        stream.write("END\n")
    selected_residues = {(chain.name, int(residue.seqid.num), _icode(residue.seqid.icode)): residue.name
                         for chain in selected[0] for residue in chain if residue.het_flag == "A"}
    missing = []
    standing = "Observed"
    limitations = []
    geometry = unavailable_geometry(payload.get("geometrySpec"),
                                    "Selected-source geometry has not been measured.")
    try:
        fixer = PDBFixer(filename=str(protein_only))
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer_residues = {(residue.chain.id, int(residue.id), _icode(residue.insertionCode)): residue.name
                          for residue in fixer.topology.residues()}
        if fixer_residues != selected_residues:
            raise WorkError("correspondenceFailed", "PDBFixer residue identity differs from the exact selected protein")
        for mapping in (fixer.missingAtoms, fixer.missingTerminals):
            for residue, names in mapping.items():
                key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
                for missing_atom in names:
                    name = require_text(getattr(missing_atom, "name", missing_atom), "missing atom name")
                    if name in _BACKBONE:
                        raise WorkError("unsupportedBackbone", "A required retained backbone atom is missing",
                                        {"residue": _address_from_key(key, chain_map, model_index), "atomName": name})
                    missing.append(_atom_address((*key, name), chain_map, model_index))
        geometry = observe_protein_geometry(fixer.topology, fixer.positions,
                                            payload.get("geometrySpec"), chain_map, model_index)
    except WorkError:
        raise
    except (ValueError, KeyError, RuntimeError) as exc:
        standing = "Unavailable"
        missing = []
        limitations.append(f"Missing-heavy-atom provider could not assess the selected candidate: {exc}")
    sulfur = []
    for chain in selected[0]:
        for residue in chain:
            if residue.name != "CYS":
                continue
            atoms = [atom for atom in residue if atom.name == "SG" and not atom.has_altloc()]
            if len(atoms) == 1:
                key = (chain.name, int(residue.seqid.num), _icode(residue.seqid.icode))
                sulfur.append((_address_from_key(key, chain_map, model_index), atoms[0].pos))
    disulfides = []
    for first in range(len(sulfur)):
        for second in range(first + 1, len(sulfur)):
            a, b = sulfur[first][1], sulfur[second][1]
            distance = math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))
            if distance <= disulfide_distance:
                disulfides.append({"first": sulfur[first][0], "second": sulfur[second][0],
                                   "distanceAngstrom": distance})
    progress("preparationChangesObserved", {"missingAtomCount": len(missing), "candidateDisulfides": len(disulfides)})
    return {"artifacts": [artifact(directory, selected_path, "selectedProteinPreview")],
            "observations": {"selectedAtomCount": sum(len(residue) for chain in selected[0] for residue in chain),
                             "previewChains": [{"sourceChain": entry["sourceChain"],
                                                "copyId": entry["copyId"],
                                                "previewChain": entry["preparedChain"]}
                                               for entry in chain_map],
                             "missingNonbackboneHeavyAtoms": missing,
                             "possibleDisulfides": disulfides,
                             "assessmentStanding": standing,
                             "limitations": limitations,
                             "geometry": geometry},
            "provider": {"name": "Gemmi/PDBFixer candidate inspection",
                         "version": "/".join(importlib.metadata.version(name) for name in ("gemmi", "pdbfixer"))}}


def prepare_protein(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm.app import ForceField, Modeller, PDBFile
    from ProteinInMembraneSystem.worker.geometry_observations import observe_protein_geometry

    source_path = work_path(directory, payload.get("sourcePath"), "sourcePath")
    verify_sha256(source_path, payload.get("sourceSha256"), "sourceSha256")
    structure = _open_structure(source_path)
    selected, chain_map = _selected_structure(structure, payload)
    _apply_altlocs(selected, payload, chain_map)
    _check_membership(selected, payload, chain_map)
    import gemmi

    retained_partner_residues = {(chain.name, int(residue.seqid.num), _icode(residue.seqid.icode))
                                 for chain in selected[0] for residue in chain
                                 if residue.entity_type == gemmi.EntityType.NonPolymer}
    source_atom_count = sum(len(residue) for chain in selected[0] for residue in chain)
    selected_atoms = {(chain.name, int(residue.seqid.num), _icode(residue.seqid.icode), atom.name)
                      for chain in selected[0] for residue in chain for atom in residue}
    selected_hydrogens = {(chain.name, int(residue.seqid.num), _icode(residue.seqid.icode), atom.name)
                          for chain in selected[0] for residue in chain for atom in residue
                          if atom.element.name in {"H", "D"}}
    selected_path = directory / "selected-protein.pdb"
    selected.write_pdb(str(selected_path))
    progress("selectedStructure", {"chainCount": len(chain_map)})

    fixer, added_heavy, approved = _fix_heavy_atoms(selected_path, payload, chain_map)
    approved_disulfides = _approved_sulfur_pairs(payload, chain_map)
    disulfide_distance = require_number(payload.get("disulfideCandidateMaxSgDistanceAngstrom"),
                                        "disulfideCandidateMaxSgDistanceAngstrom", 0.000001)
    _establish_approved_disulfides(fixer.topology, fixer.positions, approved_disulfides, disulfide_distance)
    requested_variants = _explicit_variants(fixer.topology, payload, chain_map)
    from ProteinInMembraneSystem.worker.parameterized_structure import _force_field_files

    force_field_paths = _force_field_files(directory, payload.get("forceFieldFiles"))
    force_field = ForceField(*force_field_paths)
    before_hydrogen = {(r.chain.id, int(r.id), _icode(r.insertionCode), a.name) for r in fixer.topology.residues() for a in r.atoms() if a.element.symbol == "H"}
    modeller = Modeller(fixer.topology, fixer.positions)
    nominal_ph = require_number(payload.get("nominalPh"), "nominalPh")
    actual_variants = modeller.addHydrogens(force_field, pH=nominal_ph, variants=requested_variants)
    if _sulfur_pair_keys(modeller.topology) != approved_disulfides:
        raise WorkError("providerMismatch", "Hydrogen addition changed the exact approved disulfide pair set")
    geometry = observe_protein_geometry(modeller.topology, modeller.positions,
                                        payload.get("geometrySpec"), chain_map, payload["modelIndex"])
    for expected, actual, residue in zip(requested_variants, actual_variants, modeller.topology.residues()):
        if expected is not None and expected != actual:
            raise WorkError("providerMismatch", "Hydrogen provider returned a different residue variant", {"residue": residue.id, "expected": expected, "actual": actual})
    after_hydrogen = {(r.chain.id, int(r.id), _icode(r.insertionCode), a.name) for r in modeller.topology.residues() for a in r.atoms() if a.element.symbol == "H"}
    # Template matching is a mechanics check, not an assessment that the
    # selected chemical states or structural evidence are scientifically sound.
    force_field.createSystem(modeller.topology)
    prepared_path = directory / "prepared-protein.pdb"
    with prepared_path.open("w", encoding="utf-8") as stream:
        PDBFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
    from ProteinInMembraneSystem.worker.parameterized_structure import _full_atom_sequence, _topology_data

    prepared_bond_graph = directory / "prepared-bond-graph.json"
    prepared_bond_graph.write_text(json.dumps(_topology_data(modeller.topology), separators=(",", ":")),
                                   encoding="utf-8")
    if _full_atom_sequence(PDBFile(str(prepared_path)).topology) != _full_atom_sequence(modeller.topology):
        raise WorkError("correspondenceFailed", "Prepared PDB atom identities differ from the approved bonded topology")
    prepared_atoms: list[tuple[str, int, str, str]] = []
    correspondence_atoms: list[dict[str, Any]] = []
    for atom in modeller.topology.atoms():
        residue = atom.residue
        key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode), atom.name)
        prepared_atoms.append(key)
        address = _atom_address(key, chain_map, payload["modelIndex"])
        result_id = f"{atom.index}:{key[0]}:{key[1]}:{key[2]}:{key[3]}"
        source_id = None
        role = "generated"
        approval_id = approved.get(key)
        if key in selected_atoms:
            source_id = f"{address['residue']['model']}:{address['residue']['chain']}:{address['residue']['copyId']}:{key[1]}:{key[2]}:{key[3]}"
            role = "source"
            approval_id = None
        elif atom.element.symbol != "H" and key not in added_heavy:
            raise WorkError("correspondenceFailed", "An unexplained heavy atom appeared in the prepared result",
                            {"atom": result_id})
        molecule_role = "retainedPartner" if key[:3] in retained_partner_residues else "protein"
        atom_role = "backbone" if molecule_role == "protein" and atom.name in _BACKBONE else (
            "sidechain" if molecule_role == "protein" else "partnerAtom")
        correspondence_atoms.append({"resultAtomIndex": atom.index, "resultAtomId": result_id,
                                     "sourceAtomId": source_id, "role": role,
                                     "sourceResidue": address["residue"],
                                     "moleculeRole": molecule_role, "atomRole": atom_role,
                                     "element": atom.element.symbol if atom.element else "",
                                     "approvedChangeId": approval_id})
    if len(prepared_atoms) != len(set(prepared_atoms)):
        raise WorkError("correspondenceFailed", "Prepared atom names do not uniquely identify each retained residue atom")
    lost_heavy = (selected_atoms - selected_hydrogens) - set(prepared_atoms)
    if lost_heavy:
        raise WorkError("correspondenceFailed", "Protein preparation lost selected heavy atoms",
                        {"missing": [list(key) for key in sorted(lost_heavy)]})
    removed_hydrogens = selected_hydrogens - set(prepared_atoms)
    added_hydrogens = after_hydrogen - selected_hydrogens
    mapping_path = directory / "protein-correspondence.json"
    mapping_path.write_text(json.dumps({
        "sourceId": sha256(source_path), "resultId": sha256(prepared_path),
        "atoms": correspondence_atoms, "complete": len(correspondence_atoms) == len(prepared_atoms),
    }, indent=2), encoding="utf-8")
    variant_observations = []
    for residue, actual in zip(modeller.topology.residues(), actual_variants):
        if residue.name in _VARIABLE or actual is not None:
            key = (residue.chain.id, int(residue.id), _icode(residue.insertionCode))
            variant_observations.append({"residue": _address_from_key(key, chain_map, payload["modelIndex"]),
                                         "variant": actual or residue.name, "decisionId": None})
    disulfide_observations = [
        {"first": _address_from_key(a, chain_map, payload["modelIndex"]),
         "second": _address_from_key(b, chain_map, payload["modelIndex"])}
        for a, b in sorted(_sulfur_pair_keys(modeller.topology))
    ]
    progress("candidatePrepared")
    return {
        "artifacts": [artifact(directory, prepared_path, "preparedPdb"),
                      artifact(directory, prepared_bond_graph, "preparedBondGraph"),
                      artifact(directory, mapping_path, "correspondenceJson")],
        "observations": {
            "sourceAtomCount": source_atom_count,
            "preparedAtomCount": len(prepared_atoms),
            "retainedResidueCount": len(list(modeller.topology.residues())),
            "missingBackboneResidues": [], "noncanonicalResidues": [],
            "unresolvedAlternateLocations": [], "unparameterizedResidues": [],
            "addedHeavyAtoms": [_atom_address(key, chain_map, payload["modelIndex"]) for key in sorted(added_heavy)],
            "addedHydrogens": [_atom_address(key, chain_map, payload["modelIndex"]) for key in sorted(added_hydrogens)],
            "removedSourceHydrogens": [_atom_address(key, chain_map, payload["modelIndex"]) for key in sorted(removed_hydrogens)],
            "actualResidueVariants": variant_observations,
            "actualDisulfides": disulfide_observations,
            "geometryWarnings": geometry["limitations"],
            "geometry": geometry,
            "correspondedResultAtomCount": len(correspondence_atoms),
        },
        "provider": {"name": "Gemmi/PDBFixer/OpenMM",
                     "version": "/".join(importlib.metadata.version(name) for name in ("gemmi", "pdbfixer", "openmm"))},
    }
