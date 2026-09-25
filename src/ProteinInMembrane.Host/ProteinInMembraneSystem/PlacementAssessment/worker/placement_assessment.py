"""Run the local PPM 2.0 orientation mechanic and preserve its raw evidence.

PPM's implicit-membrane optimum is only a placement proposal. The host owns
whether the proposed topology, orientation, and explicit membrane are credible.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (WorkError, artifact, require_mapping, require_number,
                       require_text, verify_sha256, work_path)


def _atom_signature(line: str) -> tuple[str, str, str, str, str, str]:
    return (line[12:16].strip(), line[16:17].strip(), line[17:20].strip(),
            line[21:22].strip(), line[22:27].strip(), line[76:78].strip())


def _atoms(path: Path) -> list[tuple[str, str, str, str, str, str]]:
    return [_atom_signature(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() != "DUM"]


def _parse_output(path: Path, oriented: Path) -> tuple[list[str], list[float]]:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    first_model: list[str] = []
    marker_ids: list[str] = []
    marker_z: list[float] = []
    for line in lines:
        if line.startswith("ENDMDL") or line.startswith("END   "):
            break
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "DUM":
            marker_ids.append(f"{line[12:16].strip()}:{line[22:27].strip()}")
            try:
                marker_z.append(float(line[46:54]))
            except ValueError as exc:
                raise WorkError("invalidProviderOutput", "PPM plane marker has no finite z-coordinate") from exc
            continue
        if line.startswith(("ATOM  ", "HETATM", "TER   ", "CRYST1")):
            first_model.append(line)
    if not first_model or not any(line.startswith("ATOM  ") for line in first_model):
        raise WorkError("invalidProviderOutput", "PPM output has no oriented protein coordinates")
    with oriented.open("w", encoding="utf-8") as stream:
        stream.writelines(first_model)
        stream.write("END\n")
    return marker_ids, marker_z


def _finite_legend_value(output: str, label: str) -> float | None:
    # PPM output formats have varied; only a plainly labelled scalar is read.
    match = re.search(rf"(?im)^\s*{label}\s*[:=]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\b", output)
    if match is None:
        return None
    value = float(match.group(1))
    return value if math.isfinite(value) else None


def place_ppm(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    source = work_path(directory, payload.get("preparedPdbPath"), "preparedPdbPath")
    verify_sha256(source, payload.get("preparedSha256"), "preparedSha256")
    executable = Path(require_text(payload.get("ppmExecutablePath"), "ppmExecutablePath")).resolve()
    if not executable.is_file():
        raise WorkError("dependencyUnavailable", "The selected local PPM executable is not installed")
    version = require_text(payload.get("ppmVersion"), "ppmVersion")
    executable_hash = require_text(payload.get("ppmExecutableSha256"), "ppmExecutableSha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", executable_hash) is None:
        raise WorkError("invalidRequest", "PPM executable identity must be an exact SHA-256 digest")
    verify_sha256(executable, executable_hash, "ppmExecutableSha256")
    topology = require_text(payload.get("topologyKind"), "topologyKind").strip().lower()
    if topology not in {"membrane-spanning", "one-surface-associated"}:
        raise WorkError("unsupportedTopology", "Requested placement class is not an established PPM candidate route")
    nterminal_side = require_text(payload.get("ppmNterminalSide"), "ppmNterminalSide").strip().lower()
    if nterminal_side not in {"in", "out"}:
        raise WorkError("invalidSelection", "PPM requires an explicit in/out first-subunit N-terminal side")
    ppm_dir = directory / "ppm"
    ppm_dir.mkdir(exist_ok=True)
    staged = ppm_dir / "protein.pdb"
    shutil.copyfile(source, staged)
    # opm.f reads (i2,1x,a3,1x,a80) from stdin; the three-character
    # topology entry is right-padded rather than inferred from a protein type.
    heterogens = any(line.startswith("HETATM") and line[17:20].strip() not in {"HOH", "WAT"}
                     for line in staged.read_text(encoding="utf-8").splitlines())
    record = f"{1 if heterogens else 0:2d} {nterminal_side:<3} {staged.name:<80}\n"
    progress("orientationProviderStarted", {"topologyKind": topology, "ppmNterminalSide": nterminal_side})
    completed = subprocess.run([str(executable)], cwd=ppm_dir, input=record,
                               text=True, capture_output=True, check=False)
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
    oriented = directory / "oriented-protein.pdb"
    marker_ids, marker_z = _parse_output(raw, oriented)
    before = _atoms(source)
    after = _atoms(oriented)
    if len(before) != len(after) or before != after:
        raise WorkError("providerMismatch", "PPM oriented structure does not preserve the selected atom identities",
                        {"inputAtoms": len(before), "outputAtoms": len(after)})
    distinct_z = sorted(set(marker_z))
    midplane = (distinct_z[0] + distinct_z[-1]) / 2 if len(distinct_z) >= 2 else None
    thickness = distinct_z[-1] - distinct_z[0] if len(distinct_z) >= 2 else None
    warnings = ["PPM's implicit membrane is an orientation approximation, not the selected explicit lipid mixture.",
                "The executable version was supplied by the host and was not established by this PPM invocation."]
    if not marker_ids:
        warnings.append("PPM did not expose parseable plane markers in the oriented model.")
    tilt = _finite_legend_value(completed.stdout, "tilt")
    progress("orientationCandidateObserved", {"alignedAtomCount": len(after)})
    return {
        "artifacts": [artifact(directory, oriented, "orientedPdb"),
                      artifact(directory, raw, "ppmRawOutput"), artifact(directory, stdout, "ppmStdout"),
                      artifact(directory, stderr, "ppmStderr")],
        "observations": {"midplaneAngstrom": midplane, "thicknessAngstrom": thickness,
                         "tiltDegrees": tilt, "alignedSourceAtomCount": len(after),
                         "numericalOutput": [], "planeMarkerIds": marker_ids,
                         "interpretationWarnings": warnings,
                         "assumedMembrane": "PPM 2.0 implicit symmetric membrane model"},
        "provider": {"name": "PPM 2.0", "version": f"host-declared {version}; binary-sha256={executable_hash.lower()}"},
    }


def _rotated(x: float, y: float, z: float, a: float, b: float, c: float) -> tuple[float, float, float]:
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    y, z = ca * y - sa * z, sa * y + ca * z
    x, z = cb * x + sb * z, -sb * x + cb * z
    return cc * x - sc * y, sc * x + cc * y, z


def adjust_placement(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    source = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    require_text(payload.get("sourceProposalId"), "sourceProposalId")
    depth = require_number(payload.get("depthShiftAngstrom"), "depthShiftAngstrom")
    x_degrees = require_number(payload.get("tiltAboutXDegrees"), "tiltAboutXDegrees")
    y_degrees = require_number(payload.get("tiltAboutYDegrees"), "tiltAboutYDegrees")
    normal_degrees = require_number(payload.get("rotationAboutNormalDegrees"), "rotationAboutNormalDegrees")
    require_text(payload.get("rationale"), "rationale")
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atom_lines:
        raise WorkError("invalidStructure", "Placement adjustment requires a nonempty molecule-only PDB")
    if any(line[17:20].strip() == "DUM" for line in atom_lines):
        raise WorkError("invalidStructure", "PPM plane markers cannot be transformed as protein atoms")
    positions = []
    for line in atom_lines:
        try:
            positions.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError as exc:
            raise WorkError("invalidStructure", "Placement PDB contains an unreadable coordinate") from exc
    centroid = tuple(sum(point[axis] for point in positions) / len(positions) for axis in range(3))
    angles = tuple(math.radians(value) for value in (x_degrees, y_degrees, normal_degrees))
    adjusted = []
    for position in positions:
        relative = tuple(position[axis] - centroid[axis] for axis in range(3))
        rx, ry, rz = _rotated(*relative, *angles)
        result = (rx + centroid[0], ry + centroid[1], rz + centroid[2] + depth)
        if not all(math.isfinite(value) and abs(value) < 9999.999 for value in result):
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
    require_text(payload.get("proposalId"), "proposalId")
    midplane = require_number(payload.get("membraneMidplaneAngstrom"), "membraneMidplaneAngstrom")
    lower = require_number(payload.get("coreLowerZAngstrom"), "coreLowerZAngstrom")
    upper = require_number(payload.get("coreUpperZAngstrom"), "coreUpperZAngstrom")
    if not lower < midplane < upper:
        raise WorkError("invalidRequest", "Selected membrane core must surround its declared midplane")
    source_addresses = payload.get("residueAddressesInOrder")
    if not isinstance(source_addresses, list) or not source_addresses:
        raise WorkError("missingCorrespondence", "Placement measurement requires ordered source-residue identities")
    groups: list[tuple[tuple[str, str, str, str], list[tuple[str, float]]]] = []
    total = above = within = below = 0
    for line in source.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() == "DUM":
            raise WorkError("invalidStructure", "PPM plane markers cannot be measured as protein atoms")
        try:
            z = float(line[46:54])
        except ValueError as exc:
            raise WorkError("invalidStructure", "Oriented PDB has an unreadable z-coordinate") from exc
        if not math.isfinite(z):
            raise WorkError("invalidStructure", "Oriented PDB has a nonfinite z-coordinate")
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip())
        if not groups or groups[-1][0] != key:
            groups.append((key, []))
        groups[-1][1].append((line[12:16].strip(), z))
        total += 1
    if len(groups) != len(source_addresses):
        raise WorkError("correspondenceFailed", "Oriented coordinate residues differ from ordered source identities")
    if total == 0:
        raise WorkError("invalidStructure", "Oriented placement contains no protein atoms")
    observations = []
    for (chain, residue_id, insertion, name), atoms, raw_address in (
            (key, atoms, address) for (key, atoms), address in zip(groups, source_addresses)):
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
