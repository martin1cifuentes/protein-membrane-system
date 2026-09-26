"""Run the local PPM 2.0 orientation mechanic and preserve its raw evidence.

PPM's implicit-membrane optimum is only a placement proposal. The host owns
whether the proposed topology, orientation, and explicit membrane are credible.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

import numpy as np

from ProteinInMembraneSystem.worker.exchange import (WorkError, artifact, require_mapping, require_number,
                       require_text, verify_sha256, work_path)


def _atom_signature(line: str) -> tuple[str, str, str, str, str, str]:
    return (line[12:16].strip(), line[16:17].strip(), line[17:20].strip(),
            line[21:22].strip(), line[22:27].strip(), line[76:78].strip())


def _atoms(path: Path) -> list[tuple[str, str, str, str, str, str]]:
    return [_atom_signature(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() != "DUM"]


def _position(line: str, failure_code: str, label: str) -> np.ndarray:
    try:
        position = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    except ValueError as exc:
        raise WorkError(failure_code, f"{label} has an unreadable coordinate") from exc
    if not bool(np.all(np.isfinite(position))):
        raise WorkError(failure_code, f"{label} has a nonfinite coordinate")
    return position


def _element(line: str, failure_code: str, label: str) -> str:
    element = line[76:78].strip().upper()
    if not element:
        raise WorkError(failure_code, f"{label} has no explicit element identity")
    return element


def _source_atoms(lines: list[str]) -> list[str]:
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise WorkError("invalidStructure", "Prepared PDB must describe one coordinate model")
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms or any(line[17:20].strip() == "DUM" for line in atoms):
        raise WorkError("invalidStructure", "Prepared PDB needs protein atoms without PPM plane markers")
    signatures = [_atom_signature(line) for line in atoms]
    if len(set(signatures)) != len(signatures):
        raise WorkError("invalidStructure", "Prepared PDB has ambiguous repeated atom identities")
    for line in atoms:
        _position(line, "invalidStructure", "Prepared PDB atom")
        _element(line, "invalidStructure", "Prepared PDB atom")
    return atoms


def _parse_output(path: Path) -> tuple[list[str], list[str], float, float]:
    """Read only PPM's first oriented model; DUM planes are observations, not atoms."""
    protein: list[str] = []
    markers: dict[str, list[float]] = {"N": [], "O": []}
    marker_ids: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ENDMDL", "END   ", "END\n")):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        _position(line, "invalidProviderOutput", "PPM output atom")
        if line[17:20].strip() == "DUM":
            kind = line[12:16].strip()
            if kind not in markers:
                raise WorkError("invalidProviderOutput", "PPM has an unidentified plane marker")
            marker_ids.append(f"{kind}:{line[22:27].strip()}")
            markers[kind].append(float(line[46:54]))
        else:
            _element(line, "invalidProviderOutput", "PPM output atom")
            protein.append(line)
    if not protein:
        raise WorkError("invalidProviderOutput", "PPM output has no oriented protein coordinates")
    if not markers["N"] or len(markers["N"]) != len(markers["O"]):
        raise WorkError("invalidProviderOutput", "PPM did not expose paired membrane plane markers")
    lower, upper = markers["N"][0], markers["O"][0]
    if lower >= upper or any(abs(z - lower) > 0.001 for z in markers["N"]) or any(
            abs(z - upper) > 0.001 for z in markers["O"]):
        raise WorkError("invalidProviderOutput", "PPM membrane plane markers disagree")
    return protein, marker_ids, lower, upper


def _finite_legend_value(output: str, label: str) -> float | None:
    # PPM output formats have varied; only a plainly labelled scalar is read.
    match = re.search(rf"(?im)(?:^\s*|\b){label}\s*[:=]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", output)
    if match is None:
        return None
    value = float(match.group(1))
    return value if math.isfinite(value) else None


def _fit_provider_transform(source: list[str], provider: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Fit one proper rigid transform, allowing only PDB's 0.001 Å rounding error."""
    source_heavy = [line for line in source if _element(line, "invalidStructure", "Prepared atom") not in {"H", "D"}]
    provider_heavy = [line for line in provider if _element(line, "invalidProviderOutput", "PPM atom") not in {"H", "D"}]
    if [_atom_signature(line) for line in source_heavy] != [_atom_signature(line) for line in provider_heavy]:
        raise WorkError("providerMismatch", "PPM did not preserve ordered prepared heavy-atom identities",
                        {"inputHeavyAtoms": len(source_heavy), "outputHeavyAtoms": len(provider_heavy)})
    if len(provider) not in {len(provider_heavy), len(source)}:
        raise WorkError("providerMismatch", "PPM returned only part of the prepared hydrogen/deuterium set")
    if len(provider) == len(source) and [_atom_signature(line) for line in provider] != [
            _atom_signature(line) for line in source]:
        raise WorkError("providerMismatch", "PPM returned hydrogen/deuterium identities out of correspondence")
    if len(source_heavy) < 4:
        raise WorkError("providerMismatch", "Too few corresponding heavy atoms establish a rigid transform")
    before = np.stack([_position(line, "invalidStructure", "Prepared atom") for line in source_heavy])
    after = np.stack([_position(line, "invalidProviderOutput", "PPM atom") for line in provider_heavy])
    left_center, right_center = before.mean(axis=0), after.mean(axis=0)
    left, right = before - left_center, after - right_center
    u, singular, vt = np.linalg.svd(left.T @ right)
    if singular[-1] <= 1e-6:
        raise WorkError("providerMismatch", "Corresponding heavy atoms do not establish a three-dimensional transform")
    correction = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    translation = right_center - left_center @ rotation
    residual = np.linalg.norm(before @ rotation + translation - after, axis=1)
    # Two independently rounded PDB coordinate sets give observed residuals
    # below 0.001 Å for both distributed 1RSY and prepared 6QWR probes.
    if not bool(np.all(np.isfinite(residual))) or float(residual.max()) > 0.005:
        raise WorkError("providerMismatch", "PPM coordinates are not one proper rigid transform of the prepared atoms",
                        {"maxHeavyAtomResidualAngstrom": float(residual.max())})
    if len(provider) == len(source):
        full = np.stack([_position(line, "invalidProviderOutput", "PPM atom") for line in provider])
        original = np.stack([_position(line, "invalidStructure", "Prepared atom") for line in source])
        if float(np.linalg.norm(original @ rotation + translation - full, axis=1).max()) > 0.005:
            raise WorkError("providerMismatch", "PPM hydrogen coordinates disagree with its heavy-atom transform")
    return rotation, translation


def _write_oriented(source_lines: list[str], target: Path,
                    rotation: np.ndarray, translation: np.ndarray) -> int:
    count = 0
    rendered = []
    for line in source_lines:
        if line.startswith("CRYST1"):
            # This isolated, newly oriented protein has no established periodic cell.
            continue
        if line.startswith(("ATOM  ", "HETATM")):
            position = _position(line, "invalidStructure", "Prepared atom") @ rotation + translation
            fields = [f"{value:8.3f}" for value in position]
            if any(len(field) != 8 for field in fields):
                raise WorkError("unsupportedRepresentation", "Oriented atom exceeds the PDB coordinate field")
            rendered.append(f"{line[:30]}{''.join(fields)}{line[54:]}")
            count += 1
        else:
            rendered.append(line)
    with target.open("x", encoding="utf-8") as stream:
        stream.writelines(rendered)
    return count


_PPM_TIMEOUT_SECONDS = 180


def place_ppm(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    source = work_path(directory, payload.get("preparedPdbPath"), "preparedPdbPath")
    source_hash = require_text(payload.get("preparedSha256"), "preparedSha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", source_hash) is None:
        raise WorkError("invalidRequest", "Prepared coordinate identity must be an exact SHA-256 digest")
    verify_sha256(source, source_hash, "preparedSha256")
    executable = Path(require_text(payload.get("ppmExecutablePath"), "ppmExecutablePath")).resolve()
    if not executable.is_file():
        raise WorkError("dependencyUnavailable", "The selected local PPM executable is not installed")
    if not os.access(executable, os.X_OK):
        raise WorkError("dependencyUnavailable", "The selected local PPM executable cannot be run")
    version = require_text(payload.get("ppmVersion"), "ppmVersion")
    executable_hash = require_text(payload.get("ppmExecutableSha256"), "ppmExecutableSha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", executable_hash) is None:
        raise WorkError("invalidRequest", "PPM executable identity must be an exact SHA-256 digest")
    verify_sha256(executable, executable_hash, "ppmExecutableSha256")
    library = Path(require_text(payload.get("ppmResidueLibraryPath"), "ppmResidueLibraryPath")).resolve()
    if not library.is_file():
        raise WorkError("dependencyUnavailable", "The selected PPM residue library is not installed")
    library_hash = require_text(payload.get("ppmResidueLibrarySha256"), "ppmResidueLibrarySha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", library_hash) is None:
        raise WorkError("invalidRequest", "PPM residue-library identity must be an exact SHA-256 digest")
    verify_sha256(library, library_hash, "ppmResidueLibrarySha256")
    topology = require_text(payload.get("topologyKind"), "topologyKind").strip().lower()
    if topology not in {"membrane-spanning", "one-surface-associated"}:
        raise WorkError("unsupportedTopology", "Requested placement class is not an established PPM candidate route")
    nterminal_side = require_text(payload.get("ppmNterminalSide"), "ppmNterminalSide").strip().lower()
    if nterminal_side not in {"in", "out"}:
        raise WorkError("invalidSelection", "PPM requires an explicit in/out first-subunit N-terminal side")
    ppm_dir = directory / "ppm"
    ppm_dir.mkdir(exist_ok=False)
    staged = ppm_dir / "protein.pdb"
    shutil.copyfile(source, staged)
    verify_sha256(staged, source_hash, "preparedSha256")
    staged_library = ppm_dir / "res.lib"
    shutil.copyfile(library, staged_library)
    verify_sha256(staged_library, library_hash, "ppmResidueLibrarySha256")
    staged_executable = ppm_dir / "immers"
    shutil.copy2(executable, staged_executable)
    verify_sha256(staged_executable, executable_hash, "ppmExecutableSha256")
    source_lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    before = _source_atoms(source_lines)
    # opm.f reads (i2,1x,a3,1x,a80) from stdin; the three-character
    # topology entry is right-padded rather than inferred from a protein type.
    heterogens = any(line.startswith("HETATM") and line[17:20].strip() not in {"HOH", "WAT"}
                     for line in staged.read_text(encoding="utf-8").splitlines())
    record = f"{1 if heterogens else 0:2d} {nterminal_side:<3} {staged.name:<80}\n"
    progress("orientationProviderStarted", {"topologyKind": topology, "ppmNterminalSide": nterminal_side})
    try:
        completed = subprocess.run([str(staged_executable)], cwd=ppm_dir, input=record,
                                   text=True, capture_output=True, check=False,
                                   timeout=_PPM_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired as exc:
        raise WorkError("providerTimeout", "PPM exceeded the local orientation time limit",
                        {"timeoutSeconds": _PPM_TIMEOUT_SECONDS}) from exc
    stdout = ppm_dir / "ppm-stdout.txt"
    stderr = ppm_dir / "ppm-stderr.txt"
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise WorkError("providerFailed", "PPM did not complete the orientation calculation",
                        {"exitCode": completed.returncode, "stdoutPath": str(stdout), "stderrPath": str(stderr)})
    raw = ppm_dir / "proteinout.pdb"
    if not raw.is_file():
        raise WorkError("unobservedOutput", "PPM returned without producing its expected oriented PDB")
    provider_atoms, marker_ids, lower, upper = _parse_output(raw)
    rotation, translation = _fit_provider_transform(before, provider_atoms)
    oriented = directory / "oriented-protein.pdb"
    atom_count = _write_oriented(source_lines, oriented, rotation, translation)
    if atom_count != len(before) or _atoms(oriented) != _atoms(source):
        raise WorkError("providerMismatch", "The reconstructed orientation lost prepared atom identities")
    midplane = (lower + upper) / 2
    thickness = upper - lower
    warnings = ["PPM's implicit symmetric DOPC membrane is an orientation approximation, not the selected explicit lipid mixture.",
                "The executable version was supplied by the host and was not established by this PPM invocation."]
    tilt = _finite_legend_value(completed.stdout, "Tilt angle")
    if tilt is None:
        tilt = _finite_legend_value(completed.stdout, "tilt")
    progress("orientationCandidateObserved", {"alignedAtomCount": atom_count})
    return {
        "artifacts": [artifact(directory, oriented, "orientedPdb"),
                      artifact(directory, raw, "ppmRawOutput"), artifact(directory, stdout, "ppmStdout"),
                      artifact(directory, stderr, "ppmStderr")],
        "observations": {"midplaneAngstrom": midplane, "thicknessAngstrom": thickness,
                         "tiltDegrees": tilt, "alignedSourceAtomCount": atom_count,
                         "numericalOutput": ([{"name": "PPM reported tilt", "value": tilt,
                                               "unit": "degree", "scope": "local implicit DOPC calculation"}]
                                             if tilt is not None else []), "planeMarkerIds": marker_ids,
                         "interpretationWarnings": warnings,
                         "assumedMembrane": "PPM 2.0 implicit symmetric DOPC membrane"},
        "provider": {"name": "PPM 2.0", "version": f"host-declared {version}; binary-sha256={executable_hash.lower()}; res.lib-sha256={library_hash.lower()}"},
    }


def _rotated(x: float, y: float, z: float, a: float, b: float, c: float) -> tuple[float, float, float]:
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    y, z = ca * y - sa * z, sa * y + ca * z
    x, z = cb * x + sb * z, -sb * x + cb * z
    return cc * x - sc * y, sc * x + cc * y, z


def adjust_placement(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    source = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    verify_sha256(source, require_text(payload.get("orientedPdbSha256"), "orientedPdbSha256"), "orientedPdbSha256")
    require_text(payload.get("sourceProposalId"), "sourceProposalId")
    depth = require_number(payload.get("depthShiftAngstrom"), "depthShiftAngstrom")
    x_degrees = require_number(payload.get("tiltAboutXDegrees"), "tiltAboutXDegrees")
    y_degrees = require_number(payload.get("tiltAboutYDegrees"), "tiltAboutYDegrees")
    normal_degrees = require_number(payload.get("rotationAboutNormalDegrees"), "rotationAboutNormalDegrees")
    require_text(payload.get("rationale"), "rationale")
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    atom_lines = _source_atoms(lines)
    positions = [tuple(_position(line, "invalidStructure", "Placement atom")) for line in atom_lines]
    centroid = tuple(sum(point[axis] for point in positions) / len(positions) for axis in range(3))
    angles = tuple(math.radians(value) for value in (x_degrees, y_degrees, normal_degrees))
    adjusted = []
    for position in positions:
        relative = tuple(position[axis] - centroid[axis] for axis in range(3))
        rx, ry, rz = _rotated(*relative, *angles)
        result = (rx + centroid[0], ry + centroid[1], rz + centroid[2] + depth)
        if not all(math.isfinite(value) and len(f"{value:8.3f}") == 8 for value in result):
            raise WorkError("unsupportedRepresentation", "Adjusted coordinate exceeds the PDB coordinate field")
        adjusted.append(result)
    output = directory / "adjusted-placement.pdb"
    if output.resolve() == source.resolve() or output.exists():
        raise WorkError("invalidPath", "Adjustment output would overwrite an existing source or earlier observation")
    iterator = iter(adjusted)
    with output.open("w", encoding="utf-8") as stream:
        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")):
                x, y, z = next(iterator)
                stream.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            else:
                stream.write(line)
    if _atoms(source) != _atoms(output):
        raise WorkError("correspondenceFailed", "Placement adjustment altered atom membership or order")
    progress("adjustedPlacementObserved", {"atomCount": len(atom_lines)})
    return {"artifacts": [artifact(directory, output, "adjustedPdb")],
            "observations": {"sourceAtomCount": len(atom_lines), "adjustedAtomCount": len(adjusted),
                             "appliedDepthShiftAngstrom": depth, "appliedTiltAboutXDegrees": x_degrees,
                             "appliedTiltAboutYDegrees": y_degrees,
                             "appliedRotationAboutNormalDegrees": normal_degrees,
                             "geometryWarnings": []},
            "provider": {"name": "Scientific worker rigid transform", "version": "1"}}


def measure_placement(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    source = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    verify_sha256(source, require_text(payload.get("orientedPdbSha256"), "orientedPdbSha256"), "orientedPdbSha256")
    require_text(payload.get("proposalId"), "proposalId")
    midplane = require_number(payload.get("membraneMidplaneAngstrom"), "membraneMidplaneAngstrom")
    lower = require_number(payload.get("coreLowerZAngstrom"), "coreLowerZAngstrom")
    upper = require_number(payload.get("coreUpperZAngstrom"), "coreUpperZAngstrom")
    if not lower < midplane < upper:
        raise WorkError("invalidRequest", "Selected membrane core must surround its declared midplane")
    source_addresses = payload.get("residueAddressesInOrder")
    if not isinstance(source_addresses, list) or not source_addresses:
        raise WorkError("missingCorrespondence", "Placement measurement requires ordered source-residue identities")
    expected_chains = payload.get("outputChainIdsInOrder")
    if not isinstance(expected_chains, list) or len(expected_chains) != len(source_addresses):
        raise WorkError("missingCorrespondence", "Placement measurement requires each prepared output-chain identity")
    expected_atoms = payload.get("expectedResultAtomIdsInOrder")
    if not isinstance(expected_atoms, list) or not expected_atoms:
        raise WorkError("missingCorrespondence", "Placement measurement requires every prepared result-atom identity")
    groups: list[tuple[tuple[str, str, str, str], list[tuple[str, float]]]] = []
    total = above = within = below = 0
    observed_atom_ids = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() == "DUM":
            raise WorkError("invalidStructure", "PPM plane markers cannot be measured as protein atoms")
        z = float(_position(line, "invalidStructure", "Oriented PDB atom")[2])
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip())
        try:
            numeric_residue = int(key[1])
        except ValueError as exc:
            raise WorkError("invalidStructure", "Oriented PDB has a nonnumeric residue identifier") from exc
        observed_atom_ids.append(f"{total}:{key[0]}:{numeric_residue}:{key[2]}:{line[12:16].strip()}")
        if not groups or groups[-1][0] != key:
            groups.append((key, []))
        groups[-1][1].append((line[12:16].strip(), z))
        total += 1
    if len(groups) != len(source_addresses):
        raise WorkError("correspondenceFailed", "Oriented coordinate residues differ from ordered source identities")
    if total == 0:
        raise WorkError("invalidStructure", "Oriented placement contains no protein atoms")
    if observed_atom_ids != expected_atoms:
        raise WorkError("correspondenceFailed", "Oriented placement atoms differ from the prepared result identities",
                        {"expectedAtoms": len(expected_atoms), "observedAtoms": total})
    if len({key for key, _ in groups}) != len(groups):
        raise WorkError("correspondenceFailed", "Oriented placement repeats a residue out of source order")
    observations = []
    source_to_output: dict[tuple[str, str], str] = {}
    output_to_source: dict[str, tuple[str, str]] = {}
    for (chain, residue_id, insertion, name), atoms, raw_address, expected_chain in (
            (key, atoms, address, output_chain) for (key, atoms), address, output_chain
            in zip(groups, source_addresses, expected_chains)):
        address = require_mapping(raw_address, "source residue address")
        try:
            numeric_residue = int(residue_id)
        except ValueError as exc:
            raise WorkError("invalidStructure", "Oriented PDB has a nonnumeric residue identifier") from exc
        if numeric_residue != address.get("residue") or insertion != str(address.get("insertionCode") or ""):
            raise WorkError("correspondenceFailed", "Source residue number/insertion disagrees with oriented coordinate order")
        for field in ("model", "chain", "copyId"):
            if field not in address:
                raise WorkError("missingCorrespondence", f"Source residue address lacks {field}")
        source_chain = require_text(address["chain"], "source residue chain")
        copy_id = require_text(address["copyId"], "source residue copyId")
        if require_text(expected_chain, "prepared output chain") != chain:
            raise WorkError("correspondenceFailed", "Prepared output-chain mapping disagrees with oriented coordinates")
        source_key = (source_chain, copy_id)
        if (source_key in source_to_output and source_to_output[source_key] != chain) or (
                chain in output_to_source and output_to_source[chain] != source_key):
            raise WorkError("correspondenceFailed", "Source chain/copy and prepared-chain mapping conflicts")
        source_to_output[source_key] = chain
        output_to_source[chain] = source_key
        zs = [z for _, z in atoms]
        in_core = sum(lower <= z <= upper for z in zs)
        above_core = sum(z > upper for z in zs)
        below_core = sum(z < lower for z in zs)
        within += in_core
        above += above_core
        below += below_core
        observations.append({"address": address, "name": name, "atomCount": len(atoms),
                             "outputChainId": chain, "outputResidueId": residue_id,
                             "outputInsertionCode": insertion,
                             "minZAngstrom": min(zs), "maxZAngstrom": max(zs),
                             "meanZAngstrom": sum(zs) / len(zs), "atomsWithinCore": in_core,
                             "backboneAtomsWithinCore": sum(atom in {"N", "CA", "C", "O"} and lower <= z <= upper
                                                            for atom, z in atoms),
                             "atomsAboveCore": above_core, "atomsBelowCore": below_core})
    progress("placementGeometryObserved", {"residueCount": len(groups), "atomCount": total})
    return {"artifacts": [],
            "observations": {"atomCount": total, "residues": observations,
                             "atomsWithinCore": within, "atomsAboveCore": above, "atomsBelowCore": below,
                             "proteinZMinAngstrom": min(item["minZAngstrom"] for item in observations),
                             "proteinZMaxAngstrom": max(item["maxZAngstrom"] for item in observations),
                             "limitations": []},
            "provider": {"name": "Scientific worker placement geometry", "version": "1"}}
