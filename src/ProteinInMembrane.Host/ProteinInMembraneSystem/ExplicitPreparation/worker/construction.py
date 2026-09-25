"""Owner-local construction, packing, and constructed-stage observations.

The C# owning boundary decides qualification; this module returns observations.
"""

from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (WorkError, artifact, require_integer, require_mapping,
                                        require_number, require_text, sha256, verify_sha256,
                                        work_path)
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _force_field_files, _provider, _copy_molecule, _coordinate_file, _topology_data,
    _read_topology_data, _nonbonded_charge, _full_atom_sequence,
    _bond_indices, _system_bonds_match, _load_stage, _stage_measurements,
    _state_context, _final_state)

def _coordinate_extent(points: list[tuple[float, float, float]]) -> tuple[float, float, float, float, float, float]:
    if not points:
        raise WorkError("invalidStructure", "No atoms can define the requested geometry")
    return (min(p[0] for p in points), max(p[0] for p in points),
            min(p[1] for p in points), max(p[1] for p in points),
            min(p[2] for p in points), max(p[2] for p in points))


def _write_xyz(path: Path, topology: Any, positions: Any) -> None:
    from openmm import unit

    atoms = list(topology.atoms())
    coordinates = positions.value_in_unit(unit.angstrom)
    if len(atoms) != len(coordinates):
        raise WorkError("inputMismatch", "Coordinate template and topology have different atom counts")
    lines = [str(len(atoms)), "Verified molecular template"]
    for atom, coordinate in zip(atoms, coordinates):
        if atom.element is None or any(not math.isfinite(float(value)) for value in coordinate):
            raise WorkError("invalidTemplate", "An XYZ template requires identified elements and finite coordinates")
        lines.append(f"{atom.element.symbol} {coordinate[0]:.8f} {coordinate[1]:.8f} {coordinate[2]:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_xyz(path: Path, expected_elements: list[str]) -> list[tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        declared = int(lines[0].strip())
    except (IndexError, ValueError) as exc:
        raise WorkError("providerMismatch", "Packmol XYZ has no valid atom-count header") from exc
    if declared != len(expected_elements) or len(lines) != declared + 2:
        raise WorkError("providerMismatch", "Packmol XYZ atom count differs from exact template expansion")
    coordinates = []
    for index, (line, element) in enumerate(zip(lines[2:], expected_elements)):
        fields = line.split()
        if len(fields) != 4 or fields[0].upper() != element.upper():
            raise WorkError("providerMismatch", f"Packmol XYZ atom {index} does not preserve template element/order")
        try:
            point = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise WorkError("providerMismatch", "Packmol XYZ has an invalid coordinate") from exc
        if not all(math.isfinite(value) for value in point):
            raise WorkError("providerMismatch", "Packmol XYZ has a nonfinite coordinate")
        coordinates.append(point)
    return coordinates


def _signed_volume_at(points: list[tuple[float, float, float]], indices: tuple[int, int, int, int]) -> float:
    a, b, c, d = (points[index] for index in indices)
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    w = tuple(d[i] - a[i] for i in range(3))
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) -
            u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0]))


def _signed_tetrahedral_volume(points: list[tuple[float, float, float]]) -> tuple[tuple[int, int, int, int], float] | None:
    if len(points) < 4:
        return None
    for last in range(3, len(points)):
        indices = (0, 1, 2, last)
        signed = _signed_volume_at(points, indices)
        if abs(signed) > _XYZ_CHIRAL_WITNESS_MIN_ANGSTROM_CUBED:
            return indices, signed
    return None


def _disc_union_area(atoms: list[tuple[float, float, float, float]], z0: float, z1: float,
                     clearance: float, resolution: float) -> float:
    selected = [(x, y, radius + clearance) for x, y, z, radius in atoms if z0 <= z <= z1]
    if not selected:
        return 0.0
    x_min = min(x - radius for x, y, radius in selected)
    x_max = max(x + radius for x, y, radius in selected)
    y_min = min(y - radius for x, y, radius in selected)
    y_max = max(y + radius for x, y, radius in selected)
    nx = math.ceil((x_max - x_min) / resolution)
    ny = math.ceil((y_max - y_min) / resolution)
    if nx * ny > 2_000_000:
        raise WorkError("resourceRefused", "Requested projected-area grid exceeds the bounded two-million-cell calculation")
    covered: set[tuple[int, int]] = set()
    for x, y, radius in selected:
        left = max(0, math.floor((x - radius - x_min) / resolution))
        right = min(nx - 1, math.ceil((x + radius - x_min) / resolution))
        bottom = max(0, math.floor((y - radius - y_min) / resolution))
        top = min(ny - 1, math.ceil((y + radius - y_min) / resolution))
        for ix in range(left, right + 1):
            xx = x_min + (ix + 0.5) * resolution
            for iy in range(bottom, top + 1):
                yy = y_min + (iy + 0.5) * resolution
                if (xx - x) ** 2 + (yy - y) ** 2 <= radius ** 2:
                    covered.add((ix, iy))
    return len(covered) * resolution * resolution


def _aqueous_sphere_union_volume(atoms: list[tuple[float, float, float, float]],
                                 origin: tuple[float, float, float], extent: tuple[float, float, float],
                                 excluded_z: list[tuple[float, float]], clearance: float,
                                 resolution: float) -> float:
    grid = [math.ceil(length / resolution) for length in extent]
    if grid[0] * grid[1] * grid[2] > 2_000_000:
        raise WorkError("resourceRefused", "Requested aqueous-volume grid exceeds the bounded two-million-voxel calculation")
    covered: set[tuple[int, int, int]] = set()
    for x, y, z, radius in atoms:
        radius += clearance
        limits = [(max(0, math.floor((center - radius - origin[axis]) / resolution)),
                   min(grid[axis] - 1, math.ceil((center + radius - origin[axis]) / resolution)))
                  for axis, center in enumerate((x, y, z))]
        for ix in range(limits[0][0], limits[0][1] + 1):
            xx = origin[0] + (ix + 0.5) * resolution
            for iy in range(limits[1][0], limits[1][1] + 1):
                yy = origin[1] + (iy + 0.5) * resolution
                lateral = (xx - x) ** 2 + (yy - y) ** 2
                if lateral > radius ** 2:
                    continue
                for iz in range(limits[2][0], limits[2][1] + 1):
                    zz = origin[2] + (iz + 0.5) * resolution
                    if any(lower <= zz <= upper for lower, upper in excluded_z):
                        continue
                    if lateral + (zz - z) ** 2 <= radius ** 2:
                        covered.add((ix, iy, iz))
    return len(covered) * resolution ** 3


def measure_construction_inputs(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm.app import ForceField, PDBFile
    from openmm import unit

    protein_path = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    protein_graph_path = work_path(directory, payload.get("preparedBondGraphPath"), "preparedBondGraphPath")
    verify_sha256(protein_graph_path, require_text(payload.get("preparedBondGraphSha256"),
                                                   "preparedBondGraphSha256"), "preparedBondGraphSha256")
    pdb = PDBFile(str(protein_path))
    protein_bonds = _read_topology_data(protein_graph_path)
    if _full_atom_sequence(pdb.topology) != _full_atom_sequence(protein_bonds):
        raise WorkError("inputMismatch", "Oriented protein differs from approved bonded preparation")
    ff = ForceField(*_force_field_files(directory, payload.get("forceFieldFiles")))
    system = ff.createSystem(protein_bonds)
    charge = _nonbonded_charge(system)
    radii = require_mapping(payload.get("atomRadiusByElementAngstrom"), "atomRadiusByElementAngstrom")
    clearance = require_number(payload.get("projectionClearanceAngstrom"), "projectionClearanceAngstrom", 0)
    resolution = require_number(payload.get("gridResolutionAngstrom"), "gridResolutionAngstrom", 0.000001)
    regions = payload.get("leafletRegions")
    if not isinstance(regions, list) or len(regions) != 2:
        raise WorkError("invalidRequest", "Exactly upper and lower leaflet projection regions are required")
    region_by_side = {require_text(item.get("physicalSide"), "physicalSide").lower(): item for item in regions}
    if set(region_by_side) != {"upper", "lower"}:
        raise WorkError("invalidRequest", "Leaflet projection regions must identify upper and lower sides")
    points = []
    radius_atoms = []
    for atom, position in zip(pdb.topology.atoms(), pdb.positions):
        x, y, z = position.value_in_unit(unit.angstrom)
        element = atom.element.symbol if atom.element else ""
        radius = require_number(radii.get(element), f"radius for {element}", 0)
        points.append((x, y, z))
        radius_atoms.append((x, y, z, radius))
    x0, x1, y0, y1, z0, z1 = _coordinate_extent(points)
    areas = {}
    midplane = require_number(payload.get("membraneMidplaneAngstrom"), "membraneMidplaneAngstrom")
    region_bounds = {}
    for side, item in region_by_side.items():
        minimum = require_number(item.get("zMinAngstrom"), "leaflet zMinAngstrom")
        maximum = require_number(item.get("zMaxAngstrom"), "leaflet zMaxAngstrom")
        if minimum >= maximum or (side == "upper" and minimum < midplane) or (side == "lower" and maximum > midplane):
            raise WorkError("invalidRequest", "Leaflet projection region does not correspond to the declared membrane side")
        areas[side] = _disc_union_area(radius_atoms, minimum, maximum, clearance, resolution)
        region_bounds[side] = (minimum, maximum)
    lateral = require_number(payload.get("lateralClearanceAngstrom"), "lateralClearanceAngstrom", 0)
    water = require_number(payload.get("waterMarginAngstrom"), "waterMarginAngstrom", 0)
    head = require_number(payload.get("headRegionThicknessAngstrom"), "headRegionThicknessAngstrom", 0)
    volume_resolution = require_number(payload.get("volumeGridResolutionAngstrom"), "volumeGridResolutionAngstrom", 0.000001)
    candidate_origin = (x0 - lateral, y0 - lateral,
                        min(z0, region_bounds["lower"][0] - head) - water)
    candidate_max = (x1 + lateral, y1 + lateral,
                     max(z1, region_bounds["upper"][1] + head) + water)
    candidate_extent = tuple(candidate_max[index] - candidate_origin[index] for index in range(3))
    excluded = [(region_bounds["lower"][0] - head, region_bounds["lower"][1]),
                (region_bounds["upper"][0], region_bounds["upper"][1] + head)]
    aqueous_volume = _aqueous_sphere_union_volume(radius_atoms, candidate_origin, candidate_extent,
                                                    excluded, clearance, volume_resolution)
    progress("constructionInputsObserved")
    return {"artifacts": [],
            "observations": {"proteinXMinAngstrom": x0, "proteinXMaxAngstrom": x1,
                             "proteinYMinAngstrom": y0, "proteinYMaxAngstrom": y1,
                             "proteinZMinAngstrom": z0, "proteinZMaxAngstrom": z1,
                             "upperOccludedAreaAngstromSquared": areas["upper"],
                             "lowerOccludedAreaAngstromSquared": areas["lower"],
                             "proteinAqueousOccludedVolumeAngstromCubed": aqueous_volume,
                             "candidateCellOriginAngstrom": list(candidate_origin),
                             "candidateCellAngstrom": list(candidate_extent),
                             "netChargeElementary": charge,
                             "approximationWarnings": ["Projected occlusion is a declared-grid disc union and aqueous exclusion is a declared-grid sphere union, not exact lipid or water accessible volume."]},
            "provider": _provider("OpenMM parameterization and declared-grid geometry")}


def _region(value: Any, cell: list[float], origin: list[float]) -> tuple[float, float, float, float, float, float]:
    item = require_mapping(value, "component regionAngstrom")
    bounds = tuple(require_number(item.get(key), key) for key in ("x0", "y0", "z0", "x1", "y1", "z1"))
    if not (origin[0] <= bounds[0] < bounds[3] <= origin[0] + cell[0] and
            origin[1] <= bounds[1] < bounds[4] <= origin[1] + cell[1] and
            origin[2] <= bounds[2] < bounds[5] <= origin[2] + cell[2]):
        raise WorkError("invalidRequest", "Component region must be positive and inside the declared cell")
    return bounds


def _packmol_component_script(name: str, count: int, bounds: tuple[float, ...],
                              heads: list[int], head_bounds: tuple[float, ...] | None) -> list[str]:
    lines = [f"structure {name}", f"  number {count}",
             "  inside box " + " ".join(f"{number:.5f}" for number in bounds)]
    if heads:
        if head_bounds is None:
            raise WorkError("invalidRequest", "Lipids require an explicit policy-derived head region")
        lines += ["  atoms " + " ".join(str(value) for value in heads),
                  "    inside box " + " ".join(f"{number:.5f}" for number in head_bounds),
                  "  end atoms"]
    lines.append("end structure")
    return lines


def _bounded_executable(value: Any, name: str) -> Path:
    executable = Path(require_text(value, name)).resolve()
    if not executable.is_file():
        raise WorkError("dependencyUnavailable", f"{name} does not identify an installed local executable")
    return executable


def _atom_sequence(topology: Any) -> list[tuple[str, str]]:
    return [(atom.name, atom.element.symbol if atom.element else "") for atom in topology.atoms()]


def construct_system(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm import CMMotionRemover, Context, NonbondedForce, Vec3, VerletIntegrator, XmlSerializer, unit
    from openmm.app import ForceField, HBonds, PME, PDBFile, PDBxFile, Topology

    protein_path = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    prepared_path = work_path(directory, payload.get("preparedPdbPath"), "preparedPdbPath")
    prepared_sha = require_text(payload.get("preparedPdbSha256"), "preparedPdbSha256")
    verify_sha256(prepared_path, prepared_sha, "preparedPdbSha256")
    prepared_mapping = require_mapping(payload.get("preparedCorrespondence"), "preparedCorrespondence")
    protein_graph_path = work_path(directory, payload.get("preparedBondGraphPath"), "preparedBondGraphPath")
    verify_sha256(protein_graph_path, require_text(payload.get("preparedBondGraphSha256"),
                                                   "preparedBondGraphSha256"), "preparedBondGraphSha256")
    executable = _bounded_executable(payload.get("packmolExecutablePath"), "packmolExecutablePath")
    version = require_text(payload.get("packmolVersion"), "packmolVersion")
    executable_hash = require_text(payload.get("packmolExecutableSha256"), "packmolExecutableSha256")
    if re.fullmatch(r"[0-9a-fA-F]{64}", executable_hash) is None:
        raise WorkError("invalidRequest", "Packmol executable identity must be an exact SHA-256 digest")
    verify_sha256(executable, executable_hash, "packmolExecutableSha256")
    cell_input = payload.get("cellAngstrom")
    if not isinstance(cell_input, list) or len(cell_input) != 3:
        raise WorkError("invalidRequest", "cellAngstrom must contain three orthorhombic dimensions")
    cell = [require_number(value, "cell dimension", 20) for value in cell_input]
    origin_input = payload.get("cellOriginAngstrom")
    if not isinstance(origin_input, list) or len(origin_input) != 3:
        raise WorkError("invalidRequest", "cellOriginAngstrom must contain three coordinates")
    origin = [require_number(value, "cell origin") for value in origin_input]
    tolerance = require_number(payload.get("toleranceAngstrom"), "toleranceAngstrom", 0.1)
    components = payload.get("components")
    if not isinstance(components, list) or not components:
        raise WorkError("invalidRequest", "At least one explicit molecular component is required")
    protein = PDBFile(str(protein_path))
    prepared_protein = PDBFile(str(prepared_path))
    protein_bonds = _read_topology_data(protein_graph_path)
    if _full_atom_sequence(protein.topology) != _full_atom_sequence(protein_bonds):
        raise WorkError("inputMismatch", "Oriented protein atom identities differ from approved prepared bond graph")
    if _full_atom_sequence(protein.topology) != _full_atom_sequence(prepared_protein.topology):
        raise WorkError("inputMismatch", "Oriented protein atom identities differ from the prepared source")
    mapped_atoms = prepared_mapping.get("atoms")
    if (prepared_mapping.get("complete") is not True or
            prepared_mapping.get("resultId") != prepared_sha or
            not isinstance(mapped_atoms, list) or
            len(mapped_atoms) != len(list(protein.topology.atoms()))):
        raise WorkError("inputMismatch", "Prepared correspondence is incomplete or addresses another protein")
    if sorted(require_integer(require_mapping(item, "prepared atom mapping").get("resultAtomIndex"),
                              "resultAtomIndex") for item in mapped_atoms) != list(range(len(mapped_atoms))):
        raise WorkError("inputMismatch", "Prepared correspondence does not address each ordered atom exactly once")
    mapped_atoms.sort(key=lambda item: item["resultAtomIndex"])
    for index, (mapped, atom) in enumerate(zip(mapped_atoms, prepared_protein.topology.atoms())):
        if (mapped.get("moleculeRole") not in {"protein", "retainedPartner"} or
                mapped.get("element") != (atom.element.symbol if atom.element else "") or
                mapped.get("role") not in {"source", "generated"} or
                (mapped.get("role") == "source" and not mapped.get("sourceAtomId")) or
                (mapped.get("role") == "generated" and mapped.get("sourceAtomId") is not None)):
            raise WorkError("inputMismatch", f"Prepared correspondence atom {index} conflicts with its source")
    protein_positions = [tuple(position.value_in_unit(unit.angstrom)) for position in protein.positions]
    x0, x1, y0, y1, z0, z1 = _coordinate_extent(protein_positions)
    if not (origin[0] < x0 and x1 < origin[0] + cell[0] and
            origin[1] < y0 and y1 < origin[1] + cell[1] and
            origin[2] < z0 and z1 < origin[2] + cell[2]):
        raise WorkError("invalidGeometry", "Oriented protein does not fit inside the derived periodic cell")
    pack_dir = directory / "packing"
    pack_dir.mkdir(exist_ok=True)
    _write_xyz(pack_dir / "protein.xyz", protein.topology, protein.positions)
    expected = Topology()
    _copy_molecule(expected, protein_bonds)
    expected_sequence = _atom_sequence(protein.topology)
    role_by_index = [(item["moleculeRole"], require_text(item.get("atomRole"), "prepared atomRole"))
                     for item in mapped_atoms]
    side_by_index: list[str | None] = [None] * len(role_by_index)
    provenance_by_index = [{"sourceAtomId": item.get("sourceAtomId"),
                            "sourceResidue": item.get("sourceResidue"),
                            "approvedChangeId": item.get("approvedChangeId"),
                            "role": item["role"], "generatedSpeciesId": None,
                            "generatedComponentRole": None}
                           for item in mapped_atoms]
    result_atom_ids = ["protein:" + json.dumps([index, atom.name, atom.residue.name], separators=(",", ":"))
                       for index, atom in enumerate(protein.topology.atoms())]
    script = [f"tolerance {tolerance:.5f}", "filetype xyz", "output packed.xyz",
              "pbc " + " ".join(f"{value:.5f}" for value in (*origin, *(origin[i] + cell[i] for i in range(3)))),
              "structure protein.xyz", "  number 1",
              "  fixed 0.0 0.0 0.0 0.0 0.0 0.0",
              "end structure"]
    species_counts = []
    blocks: list[tuple[int, int, str, tuple[float, ...], tuple[float, ...] | None,
                       list[int], list[tuple[float, float, float]]]] = []
    protein_atom_count = len(expected_sequence)
    water_count = positive_count = negative_count = 0
    for index, raw in enumerate(components):
        item = require_mapping(raw, "construction component")
        count = require_integer(item.get("count"), "component count")
        if count == 0:
            continue
        role = require_text(item.get("role"), "component role")
        side = require_text(item.get("physicalSide"), "component physicalSide")
        if role not in {"lipid", "water", "positiveIon", "negativeIon"} or side not in {"upper", "lower"}:
            raise WorkError("invalidRequest", "Each component needs a declared chemical role and physical side")
        species_id = require_text(item.get("speciesId"), "component speciesId")
        template_path = work_path(directory, item.get("templateCoordinatePath"), "component templateCoordinatePath")
        verify_sha256(template_path,
                      require_text(item.get("templateCoordinateSha256"), "component templateCoordinateSha256"),
                      "component templateCoordinateSha256")
        template = _coordinate_file(template_path)
        atoms = list(template.topology.atoms())
        if len(list(template.topology.residues())) != 1 or not atoms:
            raise WorkError("invalidTemplate", "Each component coordinate template must be one complete molecule")
        heads = item.get("headAtomIndices")
        if not isinstance(heads, list) or any(not isinstance(number, int) or number < 1 or number > len(atoms) for number in heads):
            raise WorkError("invalidRequest", "headAtomIndices must be one-based template atom indices")
        if role == "lipid" and not heads:
            raise WorkError("invalidTemplate", "A lipid requires explicit head-atom placement indices")
        bounds = _region(item.get("regionAngstrom"), cell, origin)
        raw_head_bounds = item.get("headRegionAngstrom")
        head_bounds = _region(raw_head_bounds, cell, origin) if raw_head_bounds is not None else None
        if role == "lipid":
            if head_bounds is None:
                raise WorkError("invalidRequest", "Upper/lower lipid packing requires a declared head region")
            if any(head_bounds[i] < bounds[i] or head_bounds[i+3] > bounds[i+3] for i in range(3)):
                raise WorkError("invalidRequest", "Head region must be within the molecule region")
        elif head_bounds is not None:
            raise WorkError("invalidRequest", "Only a lipid may declare a head region")
        staged_name = f"component-{index:03d}.xyz"
        _write_xyz(pack_dir / staged_name, template.topology, template.positions)
        script += _packmol_component_script(staged_name, count, bounds, heads, head_bounds)
        blocks.append((len(expected_sequence), count, side, bounds, head_bounds, heads,
                       [tuple(point.value_in_unit(unit.angstrom)) for point in template.positions]))
        for copy_index in range(count):
            _copy_molecule(expected, template.topology)
            expected_sequence.extend(_atom_sequence(template.topology))
            molecule_role = "ion" if role in {"positiveIon", "negativeIon"} else role
            role_by_index.extend((molecule_role, "head" if atom.index + 1 in heads else "body")
                                 for atom in atoms)
            side_by_index.extend(side for _ in atoms)
            provenance_by_index.extend({"sourceAtomId": None, "sourceResidue": None,
                                        "approvedChangeId": None, "role": "generated",
                                        "generatedSpeciesId": species_id,
                                        "generatedComponentRole": role} for _ in atoms)
            result_atom_ids.extend("component:" + json.dumps(
                [index, side, role, species_id, copy_index, atom.index], separators=(",", ":"))
                                   for atom in atoms)
        species_counts.append({"physicalSide": side, "role": role,
                               "speciesId": species_id, "count": count})
        if role == "water":
            water_count += count
        elif role == "positiveIon":
            positive_count += count
        elif role == "negativeIon":
            negative_count += count
    recipe = pack_dir / "packmol.inp"
    recipe.write_text("\n".join(script) + "\n", encoding="utf-8")
    progress("packingProviderStarted", {"expectedAtomCount": len(expected_sequence)})
    packed = pack_dir / "packed.xyz"
    maximum_attempts = require_integer(payload.get("maximumPackingAttempts"), "maximumPackingAttempts", 1)
    stdout = stderr = None
    for attempt_index in range(maximum_attempts):
        if packed.exists():
            packed.unlink()
        completed = subprocess.run([str(executable)], cwd=pack_dir, input=recipe.read_text(encoding="utf-8"),
                                   capture_output=True, text=True, check=False)
        stdout = pack_dir / f"packmol-{attempt_index + 1:03d}-stdout.txt"
        stderr = pack_dir / f"packmol-{attempt_index + 1:03d}-stderr.txt"
        stdout.write_text(completed.stdout, encoding="utf-8")
        stderr.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode == 0 and "Success!" in completed.stdout and packed.is_file():
            break
        progress("packingAttemptUnresolved", {"attempt": attempt_index + 1,
                                               "maximumAttempts": maximum_attempts})
    else:
        raise WorkError("providerFailed", "Packmol did not establish a completed packed configuration",
                        {"maximumAttempts": maximum_attempts, "stdoutPath": str(stdout),
                         "stderrPath": str(stderr)})
    coordinates = _read_xyz(packed, [element for _, element in expected_sequence])
    if any(math.dist(actual, original) > _XYZ_FIXED_ROUNDTRIP_ANGSTROM for actual, original in
           zip(coordinates[:protein_atom_count], protein_positions)):
        raise WorkError("providerMismatch", "Packmol moved or changed the fixed oriented protein")
    for start, count, side, bounds, head_bounds, heads, template_points in blocks:
        width = len(template_points)
        source_span = [math.dist(template_points[i], template_points[j]) for i in range(width) for j in range(i)]
        tetrahedron = _signed_tetrahedral_volume(template_points)
        for copy_index in range(count):
            placed = coordinates[start + copy_index * width:start + (copy_index + 1) * width]
            if len(placed) != width or any(not all(bounds[axis] - _XYZ_REGION_ROUNDTRIP_ANGSTROM <= point[axis] <=
                                                    bounds[axis + 3] + _XYZ_REGION_ROUNDTRIP_ANGSTROM
                                                    for axis in range(3)) for point in placed):
                raise WorkError("providerMismatch", "A packed molecule lies outside its declared region")
            if head_bounds is not None:
                if any(not all(head_bounds[axis] - _XYZ_REGION_ROUNDTRIP_ANGSTROM <= placed[index - 1][axis] <=
                                   head_bounds[axis + 3] + _XYZ_REGION_ROUNDTRIP_ANGSTROM
                                   for axis in range(3)) for index in heads):
                    raise WorkError("providerMismatch", "A lipid head lies outside its declared leaflet head region")
                body = [point[2] for index, point in enumerate(placed, 1) if index not in heads]
                head_mean_z = sum(placed[index - 1][2] for index in heads) / len(heads)
                if not body or (side == "upper" and
                                head_mean_z <= sum(body) / len(body)):
                    raise WorkError("providerMismatch", "Upper-leaflet lipid head/body orientation was not established")
                if side == "lower" and head_mean_z >= sum(body) / len(body):
                    raise WorkError("providerMismatch", "Lower-leaflet lipid head/body orientation was not established")
            placed_span = [math.dist(placed[i], placed[j]) for i in range(width) for j in range(i)]
            if any(abs(a - b) > _XYZ_RIGID_DISTANCE_ROUNDTRIP_ANGSTROM for a, b in zip(source_span, placed_span)):
                raise WorkError("providerMismatch", "A packed molecule lost rigid template geometry")
            if tetrahedron is not None:
                indices, signed = tetrahedron
                result_signed = _signed_volume_at(placed, indices)
                if abs(result_signed - signed) > max(_XYZ_CHIRAL_VOLUME_ROUNDTRIP_ANGSTROM_CUBED,
                                                     abs(signed) * _XYZ_CHIRAL_VOLUME_ROUNDTRIP_FRACTION):
                    raise WorkError("providerMismatch", "A packed molecule lost template handedness")
    vectors = (Vec3(cell[0] / 10, 0, 0) * unit.nanometer,
               Vec3(0, cell[1] / 10, 0) * unit.nanometer,
               Vec3(0, 0, cell[2] / 10) * unit.nanometer)
    expected.setPeriodicBoxVectors(vectors)
    positions = [Vec3(*coordinate) * unit.angstrom for coordinate in coordinates]
    settings = require_mapping(payload.get("systemSettings"), "systemSettings")
    method_name = require_text(settings.get("nonbondedMethod"), "nonbondedMethod")
    constraint_name = require_text(settings.get("constraints"), "constraints")
    if method_name != "PME" or constraint_name != "HBonds":
        raise WorkError("unsupportedPolicy", "Only the declared PME/HBonds all-atom route is materialized")
    nonbonded_method = {"PME": PME}[method_name]
    constraints = {"HBonds": HBonds}[constraint_name]
    cutoff_nanometers = require_number(settings.get("nonbondedCutoffNanometers"),
                                        "nonbondedCutoffNanometers", 0.000001)
    if cutoff_nanometers >= min(cell) / 20:
        raise WorkError("invalidGeometry", "Nonbonded cutoff must be less than half the shortest periodic cell dimension")
    rigid_water = settings.get("rigidWater")
    if not isinstance(rigid_water, bool):
        raise WorkError("invalidRequest", "rigidWater must be an explicit policy boolean")
    ewald_error = require_number(settings.get("ewaldErrorTolerance"), "ewaldErrorTolerance", 0.000000001)
    switch_raw = settings.get("switchDistanceNanometers")
    switch_nanometers = (None if switch_raw is None else
                          require_number(switch_raw, "switchDistanceNanometers", 0.000001))
    if switch_nanometers is not None and switch_nanometers >= cutoff_nanometers:
        raise WorkError("invalidRequest", "Lennard-Jones switching must begin before the nonbonded cutoff")
    dispersion = settings.get("useDispersionCorrection")
    remove_cm_motion = settings.get("removeCMMotion")
    if not isinstance(dispersion, bool) or not isinstance(remove_cm_motion, bool):
        raise WorkError("invalidRequest", "Dispersion and center-of-mass controls must be explicit policy booleans")
    hydrogen_raw = settings.get("hydrogenMassDaltons")
    hydrogen_mass = (None if hydrogen_raw is None else
                     require_number(hydrogen_raw, "hydrogenMassDaltons", 0.000001) * unit.dalton)
    ff = ForceField(*_force_field_files(directory, payload.get("forceFieldFiles")))
    system = ff.createSystem(expected, nonbondedMethod=nonbonded_method,
                             nonbondedCutoff=cutoff_nanometers * unit.nanometer,
                             constraints=constraints, rigidWater=rigid_water,
                             ewaldErrorTolerance=ewald_error,
                             switchDistance=(None if switch_nanometers is None else
                                             switch_nanometers * unit.nanometer),
                             removeCMMotion=remove_cm_motion, hydrogenMass=hydrogen_mass)
    nonbonded_forces = [force for force in (system.getForce(index)
                          for index in range(system.getNumForces())) if isinstance(force, NonbondedForce)]
    if len(nonbonded_forces) != 1:
        raise WorkError("providerMismatch", "Selected PME policy did not produce exactly one nonbonded force")
    nonbonded = nonbonded_forces[0]
    nonbonded.setUseDispersionCorrection(dispersion)
    if (nonbonded.getNonbondedMethod() != NonbondedForce.PME or
        abs(nonbonded.getCutoffDistance().value_in_unit(unit.nanometer) - cutoff_nanometers) > 1e-10 or
        abs(nonbonded.getEwaldErrorTolerance() - ewald_error) > 1e-12 or
        nonbonded.getUseSwitchingFunction() != (switch_nanometers is not None) or
        (switch_nanometers is not None and
         abs(nonbonded.getSwitchingDistance().value_in_unit(unit.nanometer) - switch_nanometers) > 1e-10) or
        nonbonded.getUseDispersionCorrection() != dispersion or
        sum(isinstance(system.getForce(index), CMMotionRemover)
            for index in range(system.getNumForces())) != int(remove_cm_motion)):
        raise WorkError("providerMismatch", "Realized OpenMM system controls differ from the declared policy")
    constrained_pairs = {tuple(sorted(system.getConstraintParameters(index)[:2]))
                         for index in range(system.getNumConstraints())}
    bonded_hydrogens = set()
    for first, second in expected.bonds():
        if first.element.symbol == "H" or second.element.symbol == "H":
            if tuple(sorted((first.index, second.index))) not in constrained_pairs:
                raise WorkError("providerMismatch", "HBonds policy left a hydrogen bond unconstrained")
            if first.element.symbol == "H" and second.element.symbol != "H":
                bonded_hydrogens.add(first.index)
            if second.element.symbol == "H" and first.element.symbol != "H":
                bonded_hydrogens.add(second.index)
    if rigid_water:
        for residue in expected.residues():
            atoms = list(residue.atoms())
            if not atoms or any(role_by_index[atom.index][0] != "water" for atom in atoms):
                continue
            hydrogens = [atom.index for atom in atoms if atom.element.symbol == "H"]
            if any(tuple(sorted((first, second))) not in constrained_pairs
                   for first in hydrogens for second in hydrogens if first < second):
                raise WorkError("providerMismatch", "Rigid-water policy left a water hydrogen pair unconstrained")
    if hydrogen_mass is not None:
        target_hydrogen_mass = hydrogen_mass.value_in_unit(unit.dalton)
        for atom in expected.atoms():
            if atom.index in bonded_hydrogens and role_by_index[atom.index][0] != "water":
                if abs(system.getParticleMass(atom.index).value_in_unit(unit.dalton) - target_hydrogen_mass) > 1e-6:
                    raise WorkError("providerMismatch", "Hydrogen mass repartitioning differs from the declared policy")
    charge = _nonbonded_charge(system)
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator)
    context.setPositions(positions)
    initial = context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
    energy = initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        raise WorkError("nonfiniteObservation", "Initial parameterized system has nonfinite potential energy")
    topology_path = directory / "constructed-topology.cif"
    with topology_path.open("w", encoding="utf-8") as stream:
        PDBxFile.writeFile(expected, positions, stream, keepIds=True)
    topology_json = directory / "constructed-topology.json"
    topology_json.write_text(json.dumps(_topology_data(expected), separators=(",", ":")), encoding="utf-8")
    readback = PDBxFile(str(topology_path))
    if _full_atom_sequence(readback.topology) != _full_atom_sequence(expected) or len(readback.positions) != len(coordinates):
        raise WorkError("providerMismatch", "Constructed mmCIF does not preserve atom order and identity")
    system_path = directory / "constructed-system.xml"
    state_path = directory / "constructed-state.xml"
    system_path.write_text(XmlSerializer.serialize(system), encoding="utf-8")
    state_path.write_text(XmlSerializer.serialize(initial), encoding="utf-8")
    correspondence = directory / "constructed-correspondence.json"
    if (len(provenance_by_index) != len(expected_sequence) or
            len(result_atom_ids) != len(expected_sequence) or
            len(set(result_atom_ids)) != len(result_atom_ids)):
        raise WorkError("correspondenceFailed", "Constructed atom provenance is incomplete or nonunique")
    correspondence_atoms = []
    for index, atom in enumerate(expected.atoms()):
        molecule_role, atom_role = role_by_index[index]
        provenance = provenance_by_index[index]
        correspondence_atoms.append({
            "resultAtomIndex": index, "resultAtomId": result_atom_ids[index],
            "sourceAtomId": provenance["sourceAtomId"], "role": provenance["role"],
            "moleculeRole": molecule_role, "atomRole": atom_role,
            "physicalSide": side_by_index[index],
            "element": atom.element.symbol if atom.element else "",
            "sourceResidue": provenance["sourceResidue"],
            "approvedChangeId": provenance["approvedChangeId"],
            "generatedSpeciesId": provenance["generatedSpeciesId"],
            "generatedComponentRole": provenance["generatedComponentRole"]})
    correspondence.write_text(json.dumps({"sourceId": sha256(protein_path),
                                          "resultId": sha256(topology_path),
                                          "atoms": correspondence_atoms,
                                          "complete": len(correspondence_atoms) == len(expected_sequence)}, indent=2),
                              encoding="utf-8")
    from .local_state_observations import observe_local_state

    local_state = observe_local_state(expected, initial.getPositions(),
                                      {"complete": True, "atoms": correspondence_atoms},
                                      payload.get("localObservationSpec"))
    del context, integrator
    progress("constructedCoordinatesObserved", {"atomCount": len(expected_sequence)})
    return {"artifacts": [artifact(directory, topology_path, "topologyCif"),
                          artifact(directory, topology_json, "topologyJson"),
                          artifact(directory, packed, "packedXyz"),
                          artifact(directory, system_path, "systemXml"),
                          artifact(directory, state_path, "stateXml"),
                          artifact(directory, correspondence, "correspondenceJson"),
                          artifact(directory, recipe, "packmolRecipe"),
                          artifact(directory, stdout, "packmolStdout"), artifact(directory, stderr, "packmolStderr")],
            "observations": {"atomCount": len(expected_sequence), "speciesCounts": species_counts,
                             "actualCellAngstrom": cell, "waterCount": water_count,
                             "positiveIonCount": positive_count, "negativeIonCount": negative_count,
                             "netChargeElementary": charge,
                             "initialPotentialEnergyKjMol": energy,
                             "correspondedResultAtomCount": len(correspondence_atoms),
                             "contactWarnings": local_state["limitations"],
                             "geometryWarnings": local_state["limitations"],
                             "parameterWarnings": [], "localState": local_state},
            "provider": _provider("Packmol/OpenMM",
                                  f"host-declared {version}; binary-sha256={executable_hash.lower()} / OpenMM {_provider('OpenMM')['version']}")}


def observe_stage(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    from openmm import VerletIntegrator, unit
    from .local_state_observations import observe_local_state
    from ProteinInMembraneSystem.worker.geometry_observations import observe_stage_protein_geometry

    pdb, system, state = _load_stage(directory, payload, "stateXmlPath")
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = _state_context(system, state, integrator)
    observed = _final_state(context)
    measurements = _stage_measurements(observed, system.getNumParticles())
    sidecar = _read_topology_data(work_path(directory, payload.get("topologyJsonPath"), "topologyJsonPath"))
    atom_order_matched = _full_atom_sequence(pdb.topology) == _full_atom_sequence(sidecar)
    bonds_matched = (_bond_indices(pdb.topology) == _bond_indices(sidecar) and
                     _system_bonds_match(pdb.topology, system))
    correspondence_path = work_path(directory, payload.get("correspondencePath"), "correspondencePath")
    verify_sha256(correspondence_path, require_text(payload.get("correspondenceSha256"),
                                                    "correspondenceSha256"), "correspondenceSha256")
    correspondence = require_mapping(json.loads(correspondence_path.read_text(encoding="utf-8")),
                                     "stage correspondence")
    pdb.topology.setPeriodicBoxVectors(observed.getPeriodicBoxVectors())
    local_state = observe_local_state(pdb.topology, observed.getPositions(), correspondence,
                                      payload.get("localObservationSpec"))
    protein_geometry = observe_stage_protein_geometry(
        pdb.topology, observed.getPositions(), correspondence,
        payload.get("stageProteinGeometrySpec"))
    del context, integrator
    progress("stageMeasurementsObserved", {"stageKind": payload.get("stageKind")})
    return {"artifacts": [],
            "observations": {"atomCount": system.getNumParticles(), "measurements": measurements,
                             "atomOrderMatched": atom_order_matched,
                             "bondsMatched": bonds_matched,
                             "contactWarnings": local_state["limitations"],
                             "structuralWarnings": local_state["limitations"] + protein_geometry["limitations"],
                             "numericalWarnings": [], "localState": local_state,
                             "proteinGeometry": protein_geometry},
            "provider": _provider("OpenMM stage observation")}
