"""Constrained mixed-leaflet starting coordinates from exact molecular templates.

Packmol chooses initial rigid-body positions.  It does not determine molecular
identity, bonds, parameters, scientific suitability, or completed preparation.
Those identities come from the already selected templates and checked protein.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
import math
from pathlib import Path
import re
import subprocess
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import WorkError


@dataclass(frozen=True)
class PackingBlock:
    species: str
    side: str
    count: int
    head_index: int
    tail_indices: tuple[int, ...]
    head_tail_span: float
    symbols: tuple[str, ...]
    reference: Any


@dataclass(frozen=True)
class PackedBuild:
    topology: Any
    positions: Any
    proposed_counts: tuple[tuple[str, str, int], ...]
    cell: tuple[float, float, float]
    bounds: tuple[float, float, float, float, float, float]
    seed: int
    packing_attempts: int
    packmol_output_sha256: str
    packmol_recipe_sha256: str
    water_sites_before_ion_replacement: int


def allocate_counts(fractions: list[tuple[str, float]], total: int) -> dict[str, int]:
    """Largest remainders with a positive copy and <1-copy error per species."""
    if (not fractions or total < len(fractions) or
            any(not math.isfinite(f) or f <= 0 for _, f in fractions) or
            len({name for name, _ in fractions}) != len(fractions) or
            abs(sum(f for _, f in fractions) - 1.0) > 1e-9):
        raise WorkError("invalidRequest", "Leaflet fractions or requested population are incoherent")
    targets = {name: total * f for name, f in fractions}
    counts = {name: int(math.floor(value)) for name, value in targets.items()}
    remaining = total - sum(counts.values())
    for name in sorted(counts, key=lambda item: (-(targets[item] - counts[item]), item))[:remaining]:
        counts[name] += 1
    if any(counts[name] < 1 or abs(counts[name] - targets[name]) >= 1.0 - 1e-10
           for name in counts):
        raise WorkError("resourceRefused", "Finite leaflet counts cannot retain every chosen species within one molecule")
    return counts


def _head_and_tails(representation: dict[str, Any], reference: Any) -> tuple[int, tuple[int, ...], float]:
    from openmm import unit

    atoms = list(reference.topology.atoms())
    names = [atom.name for atom in atoms]
    species = representation["speciesId"]
    declared = representation.get("headAtomIndices")
    if species == "CHL1":
        if (declared != [28] or len(atoms) < 28 or names[27] != "O3"):
            raise WorkError("invalidTemplate", "CHL1's one-based hydroxyl anchor does not resolve to O3")
        neighbors = {second.name if first.name == "O3" else first.name
                     for first, second in reference.topology.bonds()
                     if first.name == "O3" or second.name == "O3"}
        if neighbors != {"C3", "H3'"}:
            raise WorkError("invalidTemplate", "CHL1's O3 hydroxyl bond identity changed")
        if representation.get("areaPerMoleculeAngstromSquared") != 38.2:
            raise WorkError("unsupportedPolicy", "CHL1's starting footprint is not the identified estimate")
        head = 27
        tails = tuple(names.index(name) for name in ("C25", "C26", "C27") if name in names)
        if len(tails) != 3:
            raise WorkError("invalidTemplate", "CHL1 lacks its identified nonpolar terminal group")
    else:
        if (not isinstance(declared, list) or len(declared) != 1 or
                not isinstance(declared[0], int) or
                declared[0] < 1 or declared[0] > len(atoms) or
                names[declared[0] - 1] != "P"):
            raise WorkError("invalidTemplate", f"{species} lacks its exact phosphate anchor")
        head = declared[0] - 1
        tails = []
        for chain in ("2", "3"):
            candidates = [(int(match.group(1)), index) for index, name in enumerate(names)
                          if (match := re.fullmatch(rf"C{chain}(\d+)", name))]
            if not candidates:
                raise WorkError("invalidTemplate", f"{species} lacks a terminal acyl anchor")
            tails.append(max(candidates)[1])
        tails = tuple(tails)
    points = [tuple(float(value) for value in position.value_in_unit(unit.angstrom))
              for position in reference.positions]
    head_point = points[head]
    tail_point = tuple(sum(points[index][axis] for index in tails) / len(tails)
                       for axis in range(3))
    span = math.dist(head_point, tail_point)
    if not math.isfinite(span) or span < 5 or span > 40:
        raise WorkError("invalidTemplate", f"{species} has no finite polar-to-tail extent")
    return head + 1, tuple(index + 1 for index in tails), span


def _write_xyz(path: Path, symbols: tuple[str, ...], points: list[tuple[float, float, float]],
               label: str) -> None:
    if len(symbols) != len(points) or not symbols:
        raise WorkError("invalidTemplate", f"{label} XYZ atom identities and coordinates differ")
    path.write_text(f"{len(symbols)}\n{label}\n" + "".join(
        f"{symbol} {point[0]:.8f} {point[1]:.8f} {point[2]:.8f}\n"
        for symbol, point in zip(symbols, points)), encoding="ascii")


def _read_xyz(path: Path, expected: list[str]) -> list[tuple[float, float, float]]:
    if not path.is_file():
        raise WorkError("providerMismatch", "Packmol did not return the identified XYZ candidate")
    lines = path.read_text(encoding="ascii").splitlines()
    try:
        declared = int(lines[0].strip())
    except (IndexError, ValueError) as error:
        raise WorkError("providerMismatch", "Packmol XYZ has no valid atom count") from error
    if declared != len(expected) or len(lines) != declared + 2:
        raise WorkError("providerMismatch", "Packmol XYZ changed the expected atom population")
    result = []
    for index, (line, symbol) in enumerate(zip(lines[2:], expected)):
        fields = line.split()
        if len(fields) != 4 or fields[0].upper() != symbol.upper():
            raise WorkError("providerMismatch", f"Packmol XYZ changed atom identity or order at {index}")
        try:
            point = tuple(float(value) for value in fields[1:])
        except ValueError as error:
            raise WorkError("providerMismatch", "Packmol XYZ has an invalid coordinate") from error
        if not all(math.isfinite(value) for value in point):
            raise WorkError("providerMismatch", "Packmol XYZ has a nonfinite coordinate")
        result.append(point)
    return result


def _packmol_version(executable: Path, expected_hash: str, version: str) -> None:
    if (not executable.is_file() or
            hashlib.sha256(executable.read_bytes()).hexdigest().lower() != expected_hash.lower()):
        raise WorkError("providerMismatch", "The identified Packmol executable is absent or changed")
    try:
        result = subprocess.run([str(executable), "-h"], capture_output=True, text=True,
                                timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise WorkError("providerUnavailable", "The identified Packmol executable cannot run") from error
    if f"Version {version}" not in result.stdout:
        raise WorkError("providerMismatch", "Packmol does not report the pinned provider version")


def _packmol_text(root: Path, blocks: list[PackingBlock], input_paths: dict[str, Path],
                  fixed: list[tuple[Path, tuple[str, ...], list[tuple[float, float, float]]]],
                  bounds: tuple[float, float, float, float, float, float], seed: int) -> tuple[str, list[str]]:
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    expected = []
    parts = ["tolerance 2.0", "filetype xyz", f"output {root / 'packed.xyz'}",
             f"pbc {xmin:.6f} {ymin:.6f} {zmin:.6f} {xmax:.6f} {ymax:.6f} {zmax:.6f}",
             f"seed {seed}"]
    for path, symbols, _ in fixed:
        parts.extend([f"structure {path}", "  number 1", "  fixed 0.0 0.0 0.0 0.0 0.0 0.0",
                      "end structure"])
        expected.extend(symbols)
    for block in blocks:
        if block.count <= 0:
            continue
        sign = 1 if block.side == "upper" else -1
        # Membrane-specific Packmol recipes use one-sided head/tail planes.
        # Two-sided windows for both anchors overconstrain a dense rigid lipid
        # bilayer and prevent tails from interdigitating near its midplane.
        head_bound = max(8.0, block.head_tail_span - 1.0)
        tail_bound = 4.0
        core_half = block.head_tail_span + 11.0
        if core_half + 15.0 > min(zmax, -zmin):
            raise WorkError("resourceRefused", "The cell leaves no aqueous layer beyond the lipid heads")
        parts.extend([f"structure {input_paths[block.species]}", f"  number {block.count}",
                      f"  atoms {block.head_index}",
                      f"    {'over' if sign > 0 else 'below'} plane 0. 0. 1. "
                      f"{sign * head_bound:.6f}",
                      "  end atoms", f"  atoms {' '.join(map(str, block.tail_indices))}",
                      f"    {'below' if sign > 0 else 'over'} plane 0. 0. 1. "
                      f"{sign * tail_bound:.6f}",
                      "  end atoms", "end structure"])
        expected.extend(block.symbols * block.count)
    return "\n".join(parts) + "\n", expected


def _fractions(raw: Any, side: str, representations: dict[str, dict[str, Any]]) -> list[tuple[str, float]]:
    if not isinstance(raw, dict) or raw.get("physicalSide") != side or not isinstance(raw.get("fractions"), list):
        raise WorkError("invalidRequest", f"The {side} leaflet has no identified composition")
    result = []
    for item in raw["fractions"]:
        if not isinstance(item, dict) or item.get("speciesId") not in representations:
            raise WorkError("invalidRequest", f"The {side} leaflet names an unavailable template")
        value = item.get("fraction")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise WorkError("invalidRequest", f"The {side} leaflet has an invalid fraction")
        result.append((item["speciesId"], float(value)))
    if (not result or len({name for name, _ in result}) != len(result) or
            any((not math.isfinite(value) or value <= 0) for _, value in result) or
            abs(sum(value for _, value in result) - 1) > 1e-9):
        raise WorkError("invalidRequest", f"The {side} leaflet fractions are incoherent")
    return result


def _template_block(species: str, reference: Any, representation: dict[str, Any],
                    side: str, count: int) -> PackingBlock:
    head, tails, span = _head_and_tails(representation, reference)
    symbols = tuple(atom.element.symbol for atom in reference.topology.atoms())
    if any(not symbol for symbol in symbols):
        raise WorkError("invalidTemplate", f"{species} lacks an exact element sequence")
    return PackingBlock(species, side, count, head, tails, span, symbols, reference)


def _protein_clearance_bounds(protein_points: list[tuple[float, float, float]],
                              spans: list[float], scale: float, maximum_cell: float,
                              minimum_padding: float) -> tuple[float, float, float, float, float, float]:
    from .construction import _coordinate_extent

    p = _coordinate_extent(protein_points)
    # The box remains centered on the membrane frame.  Centering it on the
    # protein would erase an intentional in-plane placement offset.
    packing_boundary_inset = 0.0  # Packmol PBC permits atoms at the face, Å
    pair_exclusion = 2.0  # Packmol tolerance, Å
    numerical_tolerance = 0.1  # accumulated XYZ rounding and contact tolerance, Å
    self_image_half_gap = minimum_padding * 10.0  # declared nm to Å
    face_clearance = max(self_image_half_gap,
                         pair_exclusion + packing_boundary_inset + numerical_tolerance)
    lateral_half = [max(abs(p[0]), abs(p[1])) + face_clearance,
                    max(abs(p[2]), abs(p[3])) + face_clearance]
    lateral = [max(2 * half, 2 * minimum_padding * 10.0 + numerical_tolerance) * scale
               for half in lateral_half]
    core_half = max(spans) + 11.0
    half_height = max(core_half + 16.0, abs(p[4]) + face_clearance,
                      abs(p[5]) + face_clearance) * scale
    if max(*lateral, 2 * half_height) > maximum_cell:
        raise WorkError("resourceRefused", "Protein clearance and aqueous layers exceed the declared cell bound")
    bounds = (-lateral[0] / 2, -lateral[1] / 2, -half_height,
              lateral[0] / 2, lateral[1] / 2, half_height)
    if any(not math.isfinite(value) for value in bounds):
        raise WorkError("invalidGeometry", "The retained protein has no finite packing bounds")
    return bounds


def _check_protein_images(points: list[tuple[float, float, float]],
                          bounds: tuple[float, float, float, float, float, float],
                          minimum_gap: float) -> int:
    """Check all 26 translated constructs, including edge and corner images."""
    if not points or not math.isfinite(minimum_gap) or minimum_gap <= 0:
        raise WorkError("invalidGeometry", "Protein-image check lacks finite input")
    lengths = tuple(bounds[axis + 3] - bounds[axis] for axis in range(3))
    if any(not math.isfinite(length) or length <= 0 for length in lengths):
        raise WorkError("invalidGeometry", "Protein-image check lacks a finite periodic cell")
    p = (tuple(min(point[axis] for point in points) for axis in range(3)),
         tuple(max(point[axis] for point in points) for axis in range(3)))
    cell_size = minimum_gap
    cells: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
    for point in points:
        key = tuple(math.floor(value / cell_size) for value in point)
        cells.setdefault(key, []).append(point)
    checked = 0
    for image in itertools.product((-1, 0, 1), repeat=3):
        if image == (0, 0, 0):
            continue
        checked += 1
        shift = tuple(image[axis] * lengths[axis] for axis in range(3))
        lower_bound = math.sqrt(sum(max(0.0, p[0][axis] - p[1][axis] - shift[axis],
                                        p[0][axis] + shift[axis] - p[1][axis]) ** 2
                                    for axis in range(3)))
        if lower_bound >= minimum_gap - 1e-6:
            continue
        for point in points:
            moved = tuple(point[axis] + shift[axis] for axis in range(3))
            key = tuple(math.floor(value / cell_size) for value in moved)
            for delta in itertools.product((-1, 0, 1), repeat=3):
                neighbor = tuple(key[axis] + delta[axis] for axis in range(3))
                if any(math.dist(moved, existing) < minimum_gap - 1e-6
                       for existing in cells.get(neighbor, ())):
                    crossing = ("face", "edge", "corner")[sum(value != 0 for value in image) - 1]
                    raise WorkError("periodicClash", f"The retained construct clashes with its {crossing} periodic image")
    return checked


def _pack_molecules(directory: Path, packing: dict[str, Any],
                    representations: dict[str, dict[str, Any]], references: dict[str, Any],
                    protein_symbols: tuple[str, ...],
                    protein_points: list[tuple[float, float, float]],
                    maximum_cell: float, maximum_atoms: int,
                    minimum_padding: float,
                    progress: Callable) -> tuple[Any, Any, tuple[tuple[str, str, int], ...],
                                                  tuple[float, float, float, float, float, float],
                                                  int, int, str, str]:
    from openmm import Vec3, unit
    from openmm.app import Modeller
    from ProteinInMembraneSystem.worker.exchange import require_integer, require_text

    executable = Path(require_text(packing.get("executablePath"), "packing executablePath")).resolve()
    executable_sha = require_text(packing.get("executableSha256"), "packing executableSha256")
    version = require_text(packing.get("version"), "packing version")
    if version != "21.2.3" or executable_sha.lower() != (
            "746db4014630b5df54cf210146c9d83cbc3eed806372905b58fd666ba7c89f5c"):
        raise WorkError("unsupportedPolicy", "The selected constrained-packing recipe is not verified")
    _packmol_version(executable, executable_sha, version)
    seed = require_integer(packing.get("seed"), "packing seed", 1)
    max_attempts = require_integer(packing.get("maximumAttempts"), "packing maximumAttempts", 1)
    seconds = require_integer(packing.get("secondsPerAttempt"), "packing secondsPerAttempt", 1)
    if max_attempts > 3 or seconds > 300 or seed > 2_000_000_000:
        raise WorkError("unsupportedPolicy", "Packing resource settings exceed the identified recipe")
    upper = _fractions(packing.get("upper"), "upper", representations)
    lower = _fractions(packing.get("lower"), "lower", representations)
    selected = sorted({name for name, _ in upper + lower})
    if len(selected) < 2 and upper == lower and len(upper) == 1:
        raise WorkError("unsupportedPolicy", "A pure symmetric phospholipid needs its native recipe")
    templates = {name: _template_block(name, references[name], representations[name], "upper", 1)
                 for name in selected}
    input_paths = {}
    for species in selected:
        reference = references[species]
        positions = [tuple(float(value) for value in point.value_in_unit(unit.angstrom))
                     for point in reference.positions]
        center = tuple(sum(point[axis] for point in positions) / len(positions) for axis in range(3))
        path = directory / f"packmol-{species}.xyz"
        _write_xyz(path, templates[species].symbols,
                   [tuple(point[axis] - center[axis] for axis in range(3)) for point in positions],
                   species)
        input_paths[species] = path
    protein_xyz = directory / "packmol-retained-construct.xyz"
    _write_xyz(protein_xyz, protein_symbols, protein_points, "retained construct")
    spans = [templates[name].head_tail_span for name in selected]
    failures = []
    for attempt in range(max_attempts):
        scale = 1.0 + 0.12 * attempt
        bounds = _protein_clearance_bounds(protein_points, spans, scale, maximum_cell,
                                           minimum_padding)
        _check_protein_images(protein_points, bounds, 2 * minimum_padding * 10.0)
        lx, ly = bounds[3] - bounds[0], bounds[4] - bounds[1]
        area = lx * ly
        core_atoms = [point for point in protein_points if abs(point[2]) <= max(spans) + 5]
        if core_atoms:
            x_span = max(point[0] for point in core_atoms) - min(point[0] for point in core_atoms) + 4
            y_span = max(point[1] for point in core_atoms) - min(point[1] for point in core_atoms) + 4
            occupied = min(0.8 * area, x_span * y_span)
        else:
            occupied = 0.0
        available = area - occupied
        blocks = []
        try:
            for side, fractions in (("upper", upper), ("lower", lower)):
                mean_footprint = sum(fraction * float(representations[name]["areaPerMoleculeAngstromSquared"])
                                     for name, fraction in fractions)
                if not math.isfinite(mean_footprint) or mean_footprint <= 0:
                    raise WorkError("unsupportedPolicy", "A selected lipid lacks a recorded starting footprint")
                total = max(len(fractions), int(math.floor(available / mean_footprint + 0.5)))
                counts = allocate_counts(fractions, total)
                for name, _ in fractions:
                    blocks.append(_template_block(name, references[name], representations[name],
                                                  side, counts[name]))
        except WorkError as error:
            failures.append(str(error))
            continue
        predicted = len(protein_points) + sum(block.count * len(block.symbols) for block in blocks)
        if predicted + 3 * int(area * (2 * (bounds[5] - max(spans) - 5)) / 30) > maximum_atoms:
            raise WorkError("resourceRefused", "Proposed packing and aqueous sites exceed the atom bound")
        recipe, expected = _packmol_text(directory, blocks, input_paths,
                                         [(protein_xyz, protein_symbols, protein_points)],
                                         bounds, seed + attempt)
        recipe_path = directory / "packmol-input.inp"
        output_path = directory / "packed.xyz"
        output_path.unlink(missing_ok=True)
        recipe_path.write_text(recipe, encoding="ascii")
        progress("constrainedPackingStarted", {"attempt": attempt + 1, "seed": seed + attempt,
                                               "cellAngstrom": [lx, ly, bounds[5] - bounds[2]],
                                               "lipidCounts": [{"speciesId": b.species,
                                                                "physicalSide": b.side,
                                                                "count": b.count} for b in blocks]})
        try:
            result = subprocess.run([str(executable), "-i", str(recipe_path)],
                                    cwd=directory, capture_output=True, text=True,
                                    timeout=seconds, check=False)
        except subprocess.TimeoutExpired:
            failures.append(f"Packmol attempt {attempt + 1} exceeded {seconds} seconds")
            continue
        if result.returncode != 0 or "Success!" not in result.stdout:
            failures.append(f"Packmol attempt {attempt + 1} did not converge: " +
                            (result.stdout + result.stderr)[-400:])
            continue
        points = _read_xyz(output_path, expected)
        if any(math.dist(a, b) > 0.0002 for a, b in zip(points[:len(protein_points)], protein_points)):
            raise WorkError("providerMismatch", "Packmol moved the fixed retained construct")
        from openmm.app import PDBFile
        # PDBFile is used only for the already verified, unchanged protein
        # topology; no generated packed molecule is round-tripped through PDB.
        oriented_path = Path(require_text(packing.get("orientedPdbPath"), "orientedPdbPath"))
        oriented = PDBFile(str(oriented_path))
        topology = Modeller(oriented.topology, [Vec3(*point) for point in protein_points] * unit.angstrom)
        offset = len(protein_points)
        for index, block in enumerate(blocks):
            for copy in range(block.count):
                atom_points = points[offset:offset + len(block.symbols)]
                offset += len(block.symbols)
                topology.add(block.reference.topology,
                             [Vec3(*point) for point in atom_points] * unit.angstrom)
                list(topology.topology.chains())[-1].id = f"L{index:02d}{copy:04d}"
        if offset != len(points):
            raise WorkError("providerMismatch", "Packmol returned an unaccounted molecular block")
        topology.topology.setUnitCellDimensions(
            Vec3(lx, ly, bounds[5] - bounds[2]) * unit.angstrom)
        proposed = tuple((block.species, block.side, block.count) for block in blocks)
        progress("constrainedPackingReturned", {"atomCount": topology.topology.getNumAtoms(),
                                                "lipidCounts": [{"speciesId": s, "physicalSide": side,
                                                                 "count": count} for s, side, count in proposed]})
        return (topology, blocks, proposed, bounds, seed + attempt, attempt + 1,
                hashlib.sha256(output_path.read_bytes()).hexdigest(),
                hashlib.sha256(recipe_path.read_bytes()).hexdigest())
    raise WorkError("providerFailed", "Bounded Packmol attempts did not establish a candidate: " +
                    "; ".join(failures[-3:]))
