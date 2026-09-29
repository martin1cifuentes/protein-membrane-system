"""Independent read-back of one completed Memgen construction worker exchange.

Usage: out/python/bin/python inspect_memgen_route.py CASE_DIRECTORY [MINIMIZATION_DIRECTORY]

This inspector reads provider and Amber files directly. With a complete direct
minimization exchange, it also checks the final OpenMM worker files. Worker
evidence does not establish host completion, assessment, browser delivery,
export, or biological suitability.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path: Path, *, identity_only: bool = False) -> list[dict]:
    result = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            result.append({"name": line[12:16].strip(), "residue": line[17:20].strip(),
                           "chain": line[21:22].strip(), "residueId": line[22:26].strip(),
                           "insertionCode": line[26:27].strip(),
                           "element": line[76:78].strip().upper(),
                           "xyz": None if identity_only else
                                  tuple(float(line[i:i + 8]) for i in (30, 38, 46))})
    return result


def signature(atom: dict) -> tuple:
    return (atom["name"], atom["residue"], atom["chain"], atom["residueId"],
            atom["insertionCode"], atom["element"], atom["xyz"])


def retained_input_matches_intermediate(source: list[dict], intermediate: list[dict]) -> bool:
    """Allow only Memgen's observed Na/Cl PDB element-column abbreviation.

    The provider writes retained NA/CL as N/C in PROT and packed PDBs. Atom
    names, residue and copy labels, coordinates, and all other elements must
    remain exact. Amber atomic numbers and bound source provenance are checked
    separately at the final handoff.
    """
    abbreviations = {("NA", "NA", "NA", "N"), ("CL", "CL", "CL", "C")}
    return len(source) == len(intermediate) and all(
        signature(before)[:5] == signature(after)[:5] and
        before["xyz"] == after["xyz"] and
        (before["element"] == after["element"] or
         (before["residue"], before["name"], before["element"],
          after["element"]) in abbreviations)
        for before, after in zip(source, intermediate))


def retained_provenance_matches(prepared: dict, final: dict) -> bool:
    """Bind a reordered final Amber atom to its selected source/copy identity."""
    fields = ("sourceAtomId", "role", "sourceResidue", "approvedChangeId",
              "moleculeRole", "atomRole", "element", "physicalSide")
    return (final.get("generatedComponentRole") is None and
            final.get("generatedSpeciesId") is None and
            all(prepared.get(field) == final.get(field) for field in fields))


def components(topology) -> list[list[int]]:
    parent = list(range(topology.getNumAtoms()))

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in topology.bonds():
        parent[root(second.index)] = root(first.index)
    groups = defaultdict(list)
    for index in range(len(parent)):
        groups[root(index)].append(index)
    return list(groups.values())


def recipe_blocks(path: Path) -> tuple[list[dict], list[float]]:
    text = path.read_text(encoding="ascii")
    box = re.findall(r"^pbc\s+([-+\d. ]+)\s*$", text, re.MULTILINE)
    if len(box) != 1 or len(box[0].split()) != 6:
        raise ValueError("Packmol input lacks one exact periodic box")
    bounds = [float(value) for value in box[0].split()]
    blocks = []
    for name, body in re.findall(r"^structure\s+(\S+)\s*\n(.*?)^end structure\s*$",
                                 text, re.MULTILINE | re.DOTALL):
        count = re.search(r"^\s*number\s+(\d+)\s*$", body, re.MULTILINE)
        if count is None:
            raise ValueError(f"Packmol block {name} lacks a population")
        planes = []
        for line in body.splitlines():
            words = line.split()
            if len(words) >= 2 and words[0] in {"below", "above", "over"} and words[1] == "plane":
                if len(words) != 6:
                    raise ValueError(f"Packmol block {name} has an incomplete plane")
                normal = tuple(float(value) for value in words[2:5])
                height = float(words[5])
                if not all(math.isfinite(value) for value in (*normal, height)):
                    raise ValueError(f"Packmol block {name} has a nonfinite plane")
                planes.append({"direction": words[0], "normal": normal, "height": height})
        side = None
        z_plane = next((plane for plane in planes if plane["normal"] == (0.0, 0.0, 1.0)), None)
        if z_plane:
            height = z_plane["height"]
            side = "lower" if height < 0 else "upper" if height > 0 else None
        blocks.append({"template": name, "count": int(count.group(1)),
                       "side": side, "planes": planes})
    if not blocks:
        raise ValueError("Packmol input lacks molecular blocks")
    return blocks, [bounds[i + 3] - bounds[i] for i in range(3)]


def raw_aqueous_region_account(blocks: list[dict], molecules: list[dict],
                               envelope: float = 23.0) -> tuple[bool, str, dict]:
    """Bind selected planes and measure every raw solvent/ion atom against them.

    Packmol can stop at its loop cap with small plane deviations. The 23 A
    plane is its initial-placement control, not a new all-atom acceptance floor.
    The conditioned handoff has its own declared chemistry and contact checks.
    """
    aqueous = {"WAT", "Na+", "Cl-"}
    measurements = {"examinedAtomCount": 0, "insideEnvelopeAtomCount": 0,
                    "maximumInsideDepthAngstrom": 0.0, "insideAtomsBySpeciesSide": {}}
    if not math.isfinite(envelope) or envelope != 23.0:
        return False, "The selected outer leaflet envelope is not 23 A", measurements
    inside = Counter()
    block_count = 0
    for number, block in enumerate(blocks):
        if Path(block["template"]).stem not in aqueous:
            continue
        block_count += 1
        side = block["side"]
        planes = block["planes"]
        expected_direction = "below" if side == "lower" else "above" if side == "upper" else None
        expected_height = -envelope if side == "lower" else envelope
        if (expected_direction is None or len(planes) != 1 or
                planes[0] != {"direction": expected_direction,
                              "normal": (0.0, 0.0, 1.0), "height": expected_height}):
            return False, f"{block['template']} block {number} lacks its exact aqueous plane", measurements
        members = [item for item in molecules if item["blockIndex"] == number]
        if len(members) != block["count"]:
            return False, f"{block['template']} block {number} changed its molecule count", measurements
        for molecule in members:
            if molecule["species"] != Path(block["template"]).stem or molecule["side"] != side:
                return False, f"{block['template']} block {number} changed its source identity or side", measurements
            for atom in molecule["atoms"]:
                z = atom["xyz"][2]
                if not math.isfinite(z):
                    return False, f"{block['template']} block {number} has a nonfinite atom coordinate", measurements
                measurements["examinedAtomCount"] += 1
                depth = z + envelope if side == "lower" else envelope - z
                if depth > 0:
                    measurements["insideEnvelopeAtomCount"] += 1
                    measurements["maximumInsideDepthAngstrom"] = max(
                        measurements["maximumInsideDepthAngstrom"], depth)
                    inside[(molecule["species"], side)] += 1
    if block_count == 0 or any(item["species"] in aqueous and
                               (item["blockIndex"] < 0 or item["blockIndex"] >= len(blocks) or
                                Path(blocks[item["blockIndex"]]["template"]).stem not in aqueous)
                               for item in molecules):
        return False, "A raw generated aqueous molecule lacks a selected recipe block", measurements
    measurements["insideAtomsBySpeciesSide"] = {f"{species}:{side}": count
                                                  for (species, side), count in sorted(inside.items())}
    return True, "", measurements


def inspect_minimized(directory: Path, case: Path, construction: dict,
                      report: dict, check) -> None:
    """Read a separate final OpenMM worker result using the construction IDs."""
    import gemmi
    import numpy as np
    from openmm import Context, HarmonicBondForce, Platform, VerletIntegrator, XmlSerializer, unit

    request = json.loads((directory / "request.json").read_text())
    payload = request["payload"]
    terminal = [json.loads(line) for line in (directory / "events.jsonl").read_text().splitlines()
                if json.loads(line).get("kind") in {"result", "error"}]
    if len(terminal) != 1 or terminal[0]["kind"] != "result":
        check("minimizationTerminal", False)
        return
    result = terminal[0]["payload"]
    observed = result.get("observations") or {}
    check("minimizationCorrelation", result.get("standing") == "observed" and
          result.get("requestId") == request.get("requestId") and
          result.get("studyRevisionId") == construction.get("studyRevisionId") ==
              payload.get("studyRevisionId") and
          result.get("attemptId") == construction.get("attemptId") == payload.get("attemptId") and
          result.get("stageId") == payload.get("stageId") and
          payload.get("maxIterations") == 20000 and
          payload.get("rmsForceTargetKjMolNm") == 10.0)
    original = {item["role"]: item for item in construction["artifacts"]}
    handoff = json.loads((directory / "construction-handoff.json").read_text())
    five = ("topologyCif", "topologyJson", "systemXml", "stateXml", "correspondenceJson")
    check("constructionHandoffIdentity", all(handoff["sha256ByRole"].get(role) ==
          original[role]["sha256"] and Path(original[role]["path"]).is_file() and
          digest(Path(original[role]["path"])) == original[role]["sha256"]
          for role in five) and handoff.get("requestId") == construction.get("requestId") and
          handoff.get("attemptId") == construction.get("attemptId") and
          handoff.get("studyRevisionId") == construction.get("studyRevisionId"))
    keys = ("topologyCif", "topologyJson", "systemXml", "stateXml")
    check("minimizationConsumedExactHandoff", all(
          digest(Path(payload[role + "Path"])) == payload[role + "Sha256"] ==
          original[role]["sha256"] for role in keys))
    artifacts = {item["role"]: item for item in result.get("artifacts", [])}
    check("sevenHashedMolecularArtifacts", len(artifacts) == 2 and
          set(artifacts) == {"minimizedStateXml", "minimizedCif"} and all(
          Path(item["path"]).is_file() and digest(Path(item["path"])) == item["sha256"]
          for item in artifacts.values()))
    if not {"minimizedStateXml", "minimizedCif"}.issubset(artifacts):
        return
    topology = json.loads(Path(payload["topologyJsonPath"]).read_text())
    system = XmlSerializer.deserialize(Path(payload["systemXmlPath"]).read_text())
    initial = XmlSerializer.deserialize(Path(payload["stateXmlPath"]).read_text())
    state = XmlSerializer.deserialize(Path(artifacts["minimizedStateXml"]["path"]).read_text())
    cif = gemmi.make_structure_from_block(gemmi.cif.read_file(
        artifacts["minimizedCif"]["path"]).sole_block())
    cif_entries = [(chain, residue, atom) for model in cif for chain in model
                   for residue in chain for atom in residue]
    cif_atoms = [entry[2] for entry in cif_entries]
    count = system.getNumParticles()
    check("oneOrderedFinalMolecularModel", len(cif) == 1 and
          count == len(cif_atoms) == len(topology["atoms"]) == observed.get("finalAtomCount") and
          all(atom.name == topology["atoms"][i]["name"] and
              atom.element.name.upper() == topology["atoms"][i]["element"].upper() and
              residue.name == topology["residues"][topology["atoms"][i]["residueIndex"]]["name"] and
              str(residue.seqid.num) == topology["residues"][topology["atoms"][i]["residueIndex"]]["id"] and
              chain.name == topology["chains"][topology["residues"][topology["atoms"][i]["residueIndex"]]["chainIndex"]]["id"]
              for i, (chain, residue, atom) in enumerate(cif_entries)))
    positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
    cif_positions = np.asarray([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in cif_atoms])
    difference = float(np.max(np.linalg.norm(positions - cif_positions, axis=1)))
    check("finalCifMatchesStatePositions", math.isfinite(difference) and difference <= 0.01)
    initial_positions = np.asarray(initial.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
    check("minimizationUsedConstructedState", len(initial_positions) == count and
          float(np.max(np.linalg.norm(initial_positions - positions, axis=1))) > 0)
    vectors = state.getPeriodicBoxVectors().value_in_unit(unit.angstrom)
    source_vectors = system.getDefaultPeriodicBoxVectors()
    deviation = max(abs(vectors[i][j] - source_vectors[i][j].value_in_unit(unit.angstrom))
                    for i in range(3) for j in range(3))
    check("finalCellMatchesSystemAndCif", deviation <= 0.01 and
          all(abs(getattr(cif.cell, axis) - vectors[i][i]) <= 0.01
              for i, axis in enumerate("abc")))
    bonded = {tuple(sorted(item["atomIndices"])) for item in topology["bonds"]}
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    terms = {tuple(sorted(force.getBondParameters(i)[:2])) for force in forces
             if isinstance(force, HarmonicBondForce) for i in range(force.getNumBonds())}
    constraints = {tuple(sorted(system.getConstraintParameters(i)[:2]))
                   for i in range(system.getNumConstraints())}
    extra = constraints - bonded
    check("finalBondAndConstraintInventory", bonded.issubset(terms | constraints) and
          not (terms - bonded) and all(
          topology["atoms"][a]["residueIndex"] == topology["atoms"][b]["residueIndex"] and
          topology["atoms"][a]["element"] == topology["atoms"][b]["element"] == "H"
          for a, b in extra))

    # Re-evaluate the saved System at the final positions; a serialized State
    # force vector is evidence only after this independent calculation agrees.
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "1"})
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    evaluated = context.getState(getForces=True, getEnergy=True)
    serialized_force = np.asarray(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    force = np.asarray(evaluated.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    force_difference = force - serialized_force
    maximum_force_difference = float(np.max(np.abs(force_difference)))
    rms_force_difference = math.sqrt(float(np.mean(force_difference * force_difference)))
    energy_difference = abs((evaluated.getPotentialEnergy() - state.getPotentialEnergy())
                            .value_in_unit(unit.kilojoule_per_mole))
    check("freshSystemForcesMatchSerializedState", maximum_force_difference <= 0.05 and
          rms_force_difference <= 0.01 and energy_difference <= 0.1)
    del context, integrator

    # Independently project freshly evaluated forces by solving the
    # constraint-normal Gram system, rather than using worker State forces.
    pos_nm = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    parent = list(range(count))
    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    rows = []
    maximum_error = 0.0
    for index in range(system.getNumConstraints()):
        first, second, distance = system.getConstraintParameters(index)
        delta = pos_nm[first] - pos_nm[second]
        norm = float(np.linalg.norm(delta))
        target = distance.value_in_unit(unit.nanometer)
        maximum_error = max(maximum_error, abs(norm - target) / target)
        rows.append((first, second, delta / norm))
        parent[root(second)] = root(first)
    groups = defaultdict(list)
    for row in rows:
        groups[root(row[0])].append(row)
    projected = force.copy()
    for members in groups.values():
        atoms = sorted({index for a, b, _ in members for index in (a, b)})
        local = {atom: i for i, atom in enumerate(atoms)}
        jacobian = np.zeros((len(members), 3 * len(atoms)))
        for i, (a, b, direction) in enumerate(members):
            jacobian[i, 3 * local[a]:3 * local[a] + 3] = direction
            jacobian[i, 3 * local[b]:3 * local[b] + 3] = -direction
        flat = force[atoms].reshape(-1)
        multipliers = np.linalg.solve(jacobian @ jacobian.T, jacobian @ flat)
        projected[atoms] = (flat - jacobian.T @ multipliers).reshape((len(atoms), 3))
    tangent_rms = math.sqrt(float(np.sum(projected * projected)) / count)
    raw_rms = math.sqrt(float(np.mean(force * force)))
    check("independentFinalForceAndConstraint", math.isfinite(tangent_rms) and
          tangent_rms <= 10.0 and maximum_error <= observed["appliedConstraintTolerance"] and
          abs(tangent_rms - observed["finalRmsForceKjMolNm"]) < 1e-3 and
          abs(raw_rms - observed["finalRawRmsForceKjMolNm"]) < 1e-3 and
          abs(maximum_error - observed["maximumRelativeConstraintError"]) < 1e-7 and
          observed["termination"] == "converged" and
          observed["finalTreatmentUnrestrained"] is True)
    report["measurements"]["finalMinimization"] = {
        "stageId": result["stageId"], "atomCount": count,
        "maximumCifStateDeviationAngstrom": difference,
        "constraintTangentRmsKjMolNm": tangent_rms,
        "rawComponentRmsKjMolNm": raw_rms,
        "maximumRelativeConstraintError": maximum_error,
        "maximumSerializedForceDifferenceKjMolNm": maximum_force_difference,
        "rmsSerializedForceDifferenceKjMolNm": rms_force_difference,
        "serializedEnergyDifferenceKjMol": energy_difference,
        "topologyBonds": len(bonded), "systemConstraints": len(constraints)}


def inspect(case: Path, minimization: Path | None = None) -> dict:
    from openmm import CMAPTorsionForce, NonbondedForce, XmlSerializer, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PDBxFile, PME
    from openmm.app.internal.amber_file_parser import PrmtopLoader

    case = case.resolve()
    request = json.loads((case / "request.json").read_text())
    events = [json.loads(line) for line in (case / "events.jsonl").read_text().splitlines()]
    terminal = [event for event in events if event.get("kind") in {"result", "error"}]
    report = {"case": str(case), "requestId": request.get("requestId"),
              "standing": "incomplete", "checks": {}, "measurements": {},
              "unavailableCoverage": [], "failures": []}
    manifest_path = case.parent / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        recorded = manifest.get("sourceVersionsSha256", {})
        current_source = Path(__file__).resolve().parents[3] / \
            "src/ProteinInMembrane.Host/ProteinInMembraneSystem/ExplicitPreparation/worker/packmol_memgen.py"
        bound_source = recorded.get(
            "src/ProteinInMembrane.Host/ProteinInMembraneSystem/ExplicitPreparation/worker/packmol_memgen.py")
        report["evidenceLabel"] = manifest.get("kind", "fixture with unspecified scientific scope")
        report["sourceRelation"] = ("currentSource" if bound_source and
            current_source.is_file() and digest(current_source) == bound_source else "priorOrUnboundSource")

    def check(label: str, valid: bool, detail=None):
        report["checks"][label] = bool(valid)
        if not valid:
            report["failures"].append(label if detail is None else f"{label}: {detail}")

    payload = request["payload"]
    if manifest_path.is_file():
        cases = [item for item in manifest.get("cases", [])
                 if Path(item["requestPath"]).resolve() == (case / "request.json")]
        check("oneExternallyBoundIntendedCase", len(cases) == 1)
        if len(cases) == 1:
            intended = cases[0]
            check("externalRequestDigestAndTrial", digest(case / "request.json") ==
                  intended["requestSha256"] and payload.get("trialId") == intended["trialId"])
            check("externalLeafletIntent", all(
                  [[part["speciesId"], part["fraction"]] for part in
                   payload["target" + side.capitalize()]["fractions"]] == intended[side]
                  for side in ("lower", "upper")))
            check("externalPreparedAndPlacementInputs", all(
                  payload[key + "Sha256"] == intended.get(label, manifest.get(label)) and
                  digest(Path(payload[key + "Path"])) == intended.get(label, manifest.get(label))
                  for key, label in (("preparedPdb", "preparedSha256"),
                                     ("orientedPdb", "orientedSha256"),
                                     ("preparedBondGraph", "bondGraphSha256"))))
            source_path = intended.get("sourcePath", manifest.get("sourcePath"))
            source_sha = intended.get("sourceSha256", manifest.get("sourceSha256"))
            check("externalSourceDigest", bool(source_path and source_sha) and
                  Path(source_path).is_file() and digest(Path(source_path)) == source_sha)
    else:
        check("oneExternallyBoundIntendedCase", False, "No source-versioned case manifest")

    if len(terminal) != 1 or terminal[0]["kind"] != "result":
        check("oneObservedTerminalResult", False, "No single successful terminal event")
        report["standing"] = "failed"
        return report
    result = terminal[0]["payload"]
    observed = result.get("observations") or {}
    check("correlatedObservedResult", result.get("standing") == "observed" and
          result.get("requestId") == request.get("requestId") and
          result.get("studyRevisionId") == payload.get("studyRevisionId") and
          result.get("attemptId") == payload.get("attemptId") and
          observed.get("trial", {}).get("trialId") == payload.get("trialId"))
    inputs = [(key, Path(payload[key + "Path"]), payload[key + "Sha256"])
              for key in ("preparedPdb", "orientedPdb", "preparedBondGraph")]
    inputs.extend((item["id"], Path(item["path"]), item["sha256"]) for item in
                  payload["providerAssets"] + payload["forceFieldFiles"])
    for representation in payload["selectedSpeciesRepresentations"] + [
            payload["water"], payload["sodium"], payload["chloride"]]:
        for name in ("template", "coordinateTemplate"):
            inputs.append((representation["speciesId"] + ":" + name,
                           Path(representation[name + "Path"]),
                           representation[name + "Sha256"]))
    for name, path, expected in inputs:
        check(f"boundInput:{name}", path.is_file() and digest(path) == expected)
    artifacts = {entry["role"]: entry for entry in result.get("artifacts", [])}
    check("uniqueArtifactRoles", len(artifacts) == len(result.get("artifacts", [])))
    for role, entry in artifacts.items():
        path = Path(entry["path"])
        check(f"artifact:{role}", path.is_file() and digest(path) == entry.get("sha256"))
    for entry in observed.get("providerIntermediates", []):
        check(f"intermediateTopLevel:{entry['role']}", artifacts.get(entry["role"]) == entry)
    required = {"providerInput", "providerProtein", "providerPackmolInput",
                "providerPacked", "providerLeapInput", "providerParameterized",
                "amberTopology", "amberCoordinates", "providerRestrainedRestart",
                "amberFinalRestart", "providerCommand", "providerOptions", "providerLog",
                "systemXml", "stateXml", "correspondenceJson",
                "preparedChargeLeapTopology", "providerRestrainedInput",
                "providerUnrestrainedInput"}
    check("methodStageArtifacts", required.issubset(artifacts))
    if not required.issubset(artifacts) or report["failures"]:
        report["standing"] = "failed"
        return report
    path = lambda role: Path(artifacts[role]["path"])
    trial = path("providerInput").parent
    source = records(path("providerInput"))
    prot = records(path("providerProtein"))
    raw = records(path("providerPacked"))
    check("sourceToPROT", retained_input_matches_intermediate(source, prot))
    check("sourceInRawPacking", retained_input_matches_intermediate(
          source, raw[:len(source)]))

    command = json.loads(path("providerCommand").read_text())
    options = json.loads(path("providerOptions").read_text())
    def argument(flag):
        return (command[command.index(flag) + 1]
                if command.count(flag) == 1 and command.index(flag) + 1 < len(command)
                else None)
    lower = payload["targetLower"]["fractions"]
    upper = payload["targetUpper"]["fractions"]
    lipids = ":".join(item["speciesId"] for item in lower) + "//" + \
             ":".join(item["speciesId"] for item in upper)
    ratios = ":".join(format(item["fraction"], ".17g") for item in lower) + "//" + \
             ":".join(format(item["fraction"], ".17g") for item in upper)
    flags = {"--preoriented", "--notprotonate", "--nottrim", "--noxy_cen",
             "--keep", "--pbc", "--tight_box", "--salt", "--parametrize", "--minimize"}
    forbidden = {"--random", "--vol", "--salt_override", "--nocounter", "--charmm"}
    exact = {"--pdb": "construct.pdb", "--lipids": lipids, "--ratio": ratios,
             "--leaflet": "23", "--dist": str(payload["lateralPaddingAngstrom"]),
             "--dist_wat": str(payload["aqueousPaddingAngstrom"]),
             "--salt_c": "Na+", "--salt_a": "Cl-", "--saltcon": "0.15",
             "--ffprot": "ff19SB", "--fflip": "lipid21", "--ffwat": "tip3p",
             "--engine": "sander", "--nloop": "20", "--nloop_all": "100",
             "--maxit": "20", "--tolerance": "2", "--sd_steps": "250", "--cg_steps": "250"}
    check("selectedMethodCommand", flags.issubset(command) and not forbidden.intersection(command)
          and all(argument(flag) is not None and
                  float(argument(flag)) == float(value) if flag in
                  {"--leaflet", "--dist", "--dist_wat", "--saltcon", "--nloop",
                   "--nloop_all", "--maxit", "--tolerance", "--sd_steps", "--cg_steps"}
                  else argument(flag) == value for flag, value in exact.items()))
    check("providerRecordedOptions", options.get("lipids") == [lipids] and
          options.get("ratio") == [ratios] and options.get("ffwat") == "tip3p" and
          options.get("engine") == "sander" and options.get("pbc") is True and
          options.get("tight_box") is True and options.get("salt") is True)
    provider_utils = next(Path(item["path"]) for item in payload["providerAssets"]
                          if Path(item["path"]).name == "utils.py")
    table = next(ast.literal_eval(node.value) for node in ast.parse(
        provider_utils.read_text(encoding="utf-8")).body if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "charged"
                for target in node.targets))
    track = None
    residue_estimate = 0.0
    # The fixed-width number is read independently from the consumed file.
    with path("providerInput").open(encoding="ascii") as stream:
        for line in stream:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name, number = line[17:20].strip(), line[22:26].strip()
            if name in table and number != track:
                residue_estimate += table[name]
                track = number
    preflight = PrmtopLoader(str(path("preparedChargeLeapTopology")))
    preflight_charge = sum(preflight.getCharges())
    charged_delta = preflight_charge - residue_estimate
    supplied_delta = float(argument("--charge_pdb_delta") or 0)
    check("preparedChargeAndProviderCorrection", math.isfinite(preflight_charge) and
          abs(preflight_charge - observed["conditions"]["preparedFormalChargeElementary"]) < 1e-6 and
          abs(residue_estimate - observed["conditions"]["providerResidueNameChargeElementary"]) < 1e-6 and
          abs(charged_delta - supplied_delta) < 1e-4)
    restrained = path("providerRestrainedInput").read_text()
    unrestrained = path("providerUnrestrainedInput").read_text() if "providerUnrestrainedInput" in artifacts else ""
    check("conditioningControls", all(value in restrained for value in
          ("ntr = 1", "restraint_wt = 10", "maxcyc = 500", "ncyc = 250")) and
          all(value in unrestrained for value in ("maxcyc = 500", "ncyc = 250")) and
          "ntr = 1" not in unrestrained)

    blocks, planned_cell = recipe_blocks(path("providerPackmolInput"))
    cursor = 0
    molecules = []
    template_identity = True
    for block_index, block in enumerate(blocks):
        # Templates are read only for atom inventory. The pinned Memgen
        # archive's DLPE.pdb contains a malformed coordinate for C112; the
        # packed and parameterized molecular coordinates remain strict reads.
        template = records(trial / block["template"], identity_only=True)
        for _ in range(block["count"]):
            segment = raw[cursor:cursor + len(template)]
            if len(segment) != len(template):
                break
            template_identity &= [(atom["name"], atom["element"]) for atom in segment] == [
                (atom["name"], atom["element"]) for atom in template]
            molecules.append({"species": Path(block["template"]).stem,
                              "side": block["side"], "blockIndex": block_index,
                              "atoms": segment})
            cursor += len(template)
    check("packedMoleculeInventory", cursor == len(raw) and len(molecules) ==
          sum(block["count"] for block in blocks) and molecules and
          retained_input_matches_intermediate(source, molecules[0]["atoms"]))
    check("packedTemplateAtomIdentity", template_identity)
    aqueous_ready, aqueous_detail, raw_aqueous = raw_aqueous_region_account(blocks, molecules)
    check("rawGeneratedAqueousPlaneAndInventory", aqueous_ready, aqueous_detail)
    report["measurements"]["rawGeneratedAqueousPlaneDeviation"] = raw_aqueous
    raw_packmol_log = path("providerPackmolLog").read_text(encoding="utf-8", errors="replace")
    report["measurements"]["rawPackmolSearchStoppedAtLoopCap"] = (
        "STOP: Maximum number of GENCAN loops achieved." in raw_packmol_log)
    raw_counts = Counter((item["side"], item["species"]) for item in molecules
                         if item["species"] in {x["speciesId"] for x in lower + upper})
    post = records(path("providerPostCleanup")) if "providerPostCleanup" in artifacts else raw
    remaining = Counter(signature(atom) for atom in post)
    kept = Counter()
    removed = Counter()
    for molecule in molecules:
        atoms = Counter(signature(atom) for atom in molecule["atoms"])
        if all(remaining[key] >= number for key, number in atoms.items()):
            remaining.subtract(atoms)
            kept[(molecule["side"], molecule["species"])] += 1
        elif all(remaining[key] == 0 for key in atoms):
            removed[(molecule["side"], molecule["species"])] += 1
        else:
            check("wholeMoleculeCleanup", False, molecule["species"])
    check("postCleanupMoleculeInventory", not +remaining)
    selected = {(side, item["speciesId"]): item["fraction"]
                for side, leaflet in (("lower", lower), ("upper", upper)) for item in leaflet}
    proposed = {(item["physicalSide"], item["speciesId"]): item["count"]
                for item in observed["trial"]["proposedLipidCounts"]}
    achieved = {(item["physicalSide"], item["speciesId"]): item["count"]
                for item in observed["trial"]["achievedLipidCounts"]}
    check("providerPopulationAccount", selected.keys() == proposed.keys() == achieved.keys() and
          all(raw_counts[key] == proposed[key] and kept[key] == achieved[key] and
              achieved[key] > 0 for key in selected))
    declared_removed = {(item["physicalSide"], item["speciesId"]): item["count"]
                        for item in observed["trial"]["cleanupRemovedLipidCounts"]}
    check("cleanupRemovalAccount", all(removed[key] == declared_removed.get(key, 0)
          for key in selected) and all(key in selected for key in declared_removed))
    report["measurements"]["packedLipids"] = {f"{a}:{b}": raw_counts[(a, b)] for a, b in selected}
    report["measurements"]["postCleanupLipids"] = {f"{a}:{b}": kept[(a, b)] for a, b in selected}
    if "CHL1" not in {species for _, species in selected}:
        report["unavailableCoverage"].append("Sterol cleanup branch not exercised by this case")
    else:
        provider_log = path("providerLog").read_text()
        removed_any = sum(removed.values()) > 0
        check("sterolPiercingDetectorObserved",
              "Potential piercing lipid tails will be searched" in provider_log and
              "Lipid piercing finder failed" not in provider_log and
              "Lipid piercing removal failed" not in provider_log)
        check("postSterolCleanupArtifact",
              ("providerPostCleanup" in artifacts) if removed_any else
              ("providerPostCleanup" not in artifacts and
               "No piercing lipid found!" in provider_log))

    final = AmberInpcrdFile(str(path("amberFinalRestart")))
    initial = AmberInpcrdFile(str(path("amberCoordinates")))
    amber = AmberPrmtopFile(str(path("amberTopology")), periodicBoxVectors=final.boxVectors)
    raw_top = PrmtopLoader(str(path("amberTopology")))
    fresh = amber.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
        constraints=HBonds, rigidWater=True, ewaldErrorTolerance=0.0005,
        removeCMMotion=True, hydrogenMass=None)
    saved = XmlSerializer.deserialize(path("systemXml").read_text())
    state = XmlSerializer.deserialize(path("stateXml").read_text())
    check("directAmberToOpenMMParameters", XmlSerializer.serialize(fresh) ==
          XmlSerializer.serialize(saved))
    check("amberOrderAndParticles", amber.topology.getNumAtoms() == len(final.positions) ==
          saved.getNumParticles() == observed.get("atomCount"))
    all_atoms = list(amber.topology.atoms())
    parameterized = records(path("providerParameterized"))
    pdb_residue_ordinals = []
    previous_residue = None
    ordinal = -1
    for atom in parameterized:
        residue_key = (atom["chain"], atom["residueId"], atom["insertionCode"], atom["residue"])
        if residue_key != previous_residue:
            ordinal += 1
            previous_residue = residue_key
        pdb_residue_ordinals.append(ordinal)
    # LEaP renames only the first amino-terminal H1 to H when it writes the
    # topology. The neighboring H2/H3 and N must prove this exact alias.
    n_terminal_alias = (len(parameterized) >= 4 and
        [atom["name"] for atom in parameterized[:4]] == ["N", "H1", "H2", "H3"] and
        all(pdb_residue_ordinals[i] == 0 for i in range(4)))
    def indexed_identity(index: int) -> bool:
        pdb = parameterized[index]
        amber_atom = all_atoms[index]
        alias = (n_terminal_alias and index == 1 and pdb["name"] == "H1" and
                 amber_atom.name == "H" and amber_atom.element.symbol == "H")
        raw_residue = raw_top._raw_data["RESIDUE_LABEL"][pdb_residue_ordinals[index]].strip()
        return (pdb["name"] == amber_atom.name or alias) and \
            pdb["element"].upper() == amber_atom.element.symbol.upper() and \
            pdb["residue"] == raw_residue and \
            amber_atom.residue.name == ("HOH" if raw_residue == "WAT" else raw_residue) and \
            pdb_residue_ordinals[index] == amber_atom.residue.index
    check("leapPdbToAmberIndexedIdentity", len(parameterized) == len(all_atoms) and
          all(indexed_identity(i) for i in range(len(all_atoms))))
    check("prmtopAtomicNumberImportConsistency", raw_top.has_atomic_number and
          len(raw_top._raw_data["ATOMIC_NUMBER"]) == len(all_atoms) and
          all(int(raw_top._raw_data["ATOMIC_NUMBER"][i]) == atom.element.atomic_number
              for i, atom in enumerate(all_atoms)))
    force = [saved.getForce(i) for i in range(saved.getNumForces())]
    nonbonded = [item for item in force if isinstance(item, NonbondedForce)]
    charge = sum(item.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                 for i in range(saved.getNumParticles()) for item in nonbonded)
    check("amberNetCharge", len(nonbonded) == 1 and math.isfinite(charge) and
          abs(charge) <= 1e-4 and abs(charge - observed["conditions"]["finalNetChargeElementary"]) < 1e-6)
    raw_cmap = int(raw_top._raw_data.get("CMAP_COUNT", ["0"])[0])
    imported_cmap = sum(item.getNumTorsions() for item in force if isinstance(item, CMAPTorsionForce))
    check("cmapTermCorrespondence", raw_cmap == imported_cmap)
    if raw_cmap == 0:
        report["unavailableCoverage"].append("No ff19SB CMAP term in this protein; positive CMAP preservation untested")
    state_positions = state.getPositions().value_in_unit(unit.angstrom)
    final_positions = final.positions.value_in_unit(unit.angstrom)
    max_position_delta = max(math.dist(x, y) for x, y in zip(state_positions, final_positions))
    check("conditionedRestartIsOpenMMState", max_position_delta < 1e-7)
    cell = [final.boxVectors[i][i].value_in_unit(unit.angstrom) for i in range(3)]
    check("actualPeriodicCell", all(abs(a - b) < 0.02 for a, b in
          zip(cell, observed["actualCellAngstrom"])) and
          all(abs(a - b) < 0.02 for a, b in zip(cell, planned_cell)))

    index_by_coordinate = defaultdict(list)
    for index, atom in enumerate(parameterized):
        index_by_coordinate[(atom["xyz"], atom["element"])].append(index)
    source_map = {}
    for index, atom in enumerate(source):
        matches = index_by_coordinate[(atom["xyz"], atom["element"])]
        if len(matches) == 1 and atom["name"] == parameterized[matches[0]]["name"]:
            source_map[index] = matches[0]
    prepared_mapping = payload["preparedCorrespondence"]
    prepared_atoms = prepared_mapping.get("atoms", [])
    prepared_by_index = {atom["resultAtomIndex"]: atom for atom in prepared_atoms}
    final_mapping = json.loads(path("correspondenceJson").read_text())
    final_atoms = final_mapping.get("atoms", [])
    final_by_index = {atom["resultAtomIndex"]: atom for atom in final_atoms}
    check("retainedSourceToAmberMap", len(parameterized) == len(all_atoms) and
          len(source_map) == len(set(source_map.values())) == len(source) and
          len(prepared_by_index) == len(prepared_atoms) == len(source) and
          len(final_by_index) == len(final_atoms) == len(all_atoms) and
          prepared_mapping.get("complete") is True and
          final_mapping.get("complete") is True and
          final_mapping.get("sourceId") == payload["orientedPdbSha256"] and
          final_mapping.get("resultId") == artifacts["topologyCif"]["sha256"] and
          all(indexed_identity(mapped_index) and
              source[source_index]["element"].upper() ==
                  all_atoms[mapped_index].element.symbol.upper() and
              prepared_by_index[source_index]["element"].upper() ==
                  source[source_index]["element"].upper() and
              final_by_index[mapped_index]["element"].upper() ==
                  all_atoms[mapped_index].element.symbol.upper() and
              retained_provenance_matches(prepared_by_index[source_index],
                                          final_by_index[mapped_index])
              for source_index, mapped_index in source_map.items()))
    if len(source_map) == len(source):
        graph = json.loads(Path(payload["preparedBondGraphPath"]).read_text())
        expected_bonds = {tuple(sorted((source_map[a], source_map[b])))
                          for item in graph["bonds"] for a, b in [item["atomIndices"]]}
        mapped = set(source_map.values())
        actual_bonds = {tuple(sorted((a.index, b.index))) for a, b in amber.topology.bonds()
                        if a.index in mapped and b.index in mapped}
        check("retainedPreparedBondGraph", expected_bonds == actual_bonds)
        initial_positions = initial.positions.value_in_unit(unit.angstrom)
        shift = [initial_positions[0][axis] - parameterized[0]["xyz"][axis]
                 for axis in range(3)]
        conditioned = [[source[i]["xyz"][axis] +
                        ((final_positions[source_map[i]][axis] - shift[axis] -
                          source[i]["xyz"][axis] + cell[axis] / 2) % cell[axis] -
                         cell[axis] / 2)
                        for axis in range(3)] for i in range(len(source))]
        gaps = [cell[axis] - (max(point[axis] for point in conditioned) -
                              min(point[axis] for point in conditioned)) for axis in range(3)]
        check("conditionedRetainedImageClearance", all(gap >= 0 and
              abs(gap - observed["proteinPeriodicImageGapsAngstrom"][axis]) < 1e-5 and
              gap + 0.02 >= 2 * (payload["aqueousPaddingAngstrom"] if axis == 2 else
                                 payload["lateralPaddingAngstrom"])
              for axis, gap in enumerate(gaps)))
        report["measurements"]["conditionedRetainedImageGapsAngstrom"] = gaps
    else:
        report["unavailableCoverage"].append("Retained bond graph could not be independently mapped")

    groups = components(amber.topology)
    references = {}
    for item in payload["selectedSpeciesRepresentations"]:
        top = PDBxFile(item["coordinateTemplatePath"]).topology
        formula = tuple(sorted(Counter(atom.element.symbol for atom in top.atoms()).items()))
        check(f"distinctSelectedFormula:{item['speciesId']}", formula not in references)
        references[formula] = item["speciesId"]
    lipid_counts = Counter()
    water = sodium = chloride = 0
    mapped = set(source_map.values())
    midplane = (10 * payload["membraneCenterZNanometers"] +
        (initial.positions[0][2] - parameterized[0]["xyz"][2] * unit.angstrom)
        .value_in_unit(unit.angstrom))
    check("conditionedOutputFrameMidplane", math.isfinite(midplane) and
          abs(midplane - observed["outputFrameMidplaneZAngstrom"]) < 1e-6)
    conditioned_aqueous = {"examinedAtomCount": 0, "insideEnvelopeAtomCount": 0,
                           "maximumInsideDepthAngstrom": 0.0, "insideAtomsBySpeciesSide": {}}
    conditioned_inside = Counter()
    for group in groups:
        if mapped.intersection(group):
            continue
        formula = Counter(all_atoms[i].element.symbol for i in group)
        species = references.get(tuple(sorted(formula.items())))
        aqueous_species = None
        if species:
            head = "O1" if species == "CHL1" else "P31"
            heads = [i for i in group if all_atoms[i].name == head]
            if len(heads) != 1:
                check("uniqueLipidHead", False, species)
                continue
            z = final_positions[heads[0]][2] - midplane
            z = (z + cell[2] / 2) % cell[2] - cell[2] / 2
            lipid_counts[("upper" if z > 0 else "lower", species)] += 1
        elif formula == Counter({"O": 1, "H": 2}):
            water += 1
            aqueous_species = "WAT"
        elif formula == Counter({"Na": 1}):
            sodium += 1
            aqueous_species = "Na+"
        elif formula == Counter({"Cl": 1}):
            chloride += 1
            aqueous_species = "Cl-"
        else:
            check("generatedMoleculeChemistry", False, str(formula))
        if aqueous_species is not None:
            for index in group:
                relative_z = (final_positions[index][2] - midplane + cell[2] / 2) % cell[2] - cell[2] / 2
                side = "upper" if relative_z > 0 else "lower"
                depth = 23.0 - abs(relative_z)
                conditioned_aqueous["examinedAtomCount"] += 1
                if depth > 0:
                    conditioned_aqueous["insideEnvelopeAtomCount"] += 1
                    conditioned_aqueous["maximumInsideDepthAngstrom"] = max(
                        conditioned_aqueous["maximumInsideDepthAngstrom"], float(depth))
                    conditioned_inside[(aqueous_species, side)] += 1
    conditioned_aqueous["insideAtomsBySpeciesSide"] = {
        f"{species}:{side}": count for (species, side), count in sorted(conditioned_inside.items())}
    report["measurements"]["conditionedGeneratedAqueousEnvelopeDeviation"] = conditioned_aqueous
    check("finalLipidFormulaAndLeaflet", all(lipid_counts[key] == achieved[key] for key in selected))
    conditions = observed["conditions"]
    check("finalWaterIonInventory", (water, sodium, chloride) ==
          (conditions["finalWaterCount"] - conditions["retainedWaterCount"],
           conditions["finalSodiumCount"] - conditions["retainedSodiumCount"],
           conditions["finalChlorideCount"] - conditions["retainedChlorideCount"]))
    raw_water = sum(block["count"] for block in blocks if block["template"] == "WAT.pdb")
    raw_sodium = sum(block["count"] for block in blocks if block["template"] == "Na+.pdb")
    raw_chloride = sum(block["count"] for block in blocks if block["template"] == "Cl-.pdb")
    check("providerGeneratedAqueousAccount", (raw_water, raw_sodium, raw_chloride) ==
          (conditions["providerGeneratedWaterCount"],
           conditions["providerGeneratedSodiumCount"],
           conditions["providerGeneratedChlorideCount"]))
    check("leapAqueousChangeAccount",
          water == raw_water - conditions["leapRemovedGeneratedWaterCount"] and
          sodium == raw_sodium + conditions["leapAddedSodiumCount"] -
              conditions["leapRemovedGeneratedSodiumCount"] and
          chloride == raw_chloride + conditions["leapAddedChlorideCount"] -
              conditions["leapRemovedGeneratedChlorideCount"])
    provider_log = path("providerLog").read_text()
    def logged_number(label: str) -> float | None:
        matches = re.findall(r"^\s*" + re.escape(label) + r"\s*=\s*([-+\d.]+)",
                             provider_log, re.MULTILINE)
        return float(matches[0]) if len(matches) == 1 else None
    volumes = {side: logged_number(side.capitalize() + " water box vol")
               for side in ("lower", "upper")}
    salt_numbers = {side: logged_number(side.capitalize() + " salt number")
                    for side in ("lower", "upper")}
    logged_ions = {(side, ion): logged_number(side.capitalize() + " " + sign + " charges")
                   for side in ("lower", "upper")
                   for ion, sign in (("Na+.pdb", "positive"), ("Cl-.pdb", "negative"))}
    regions = {region["side"]: region for region in conditions["aqueousRegions"]}
    planned_ions = Counter()
    for block in blocks:
        if block["template"] in {"Na+.pdb", "Cl-.pdb"}:
            planned_ions[(block["side"], block["template"])] += block["count"]
    check("providerSaltLogAndRegionalInventory", set(regions) == set(volumes) and
          all(volumes[side] is not None and salt_numbers[side] is not None and
              abs(regions[side]["estimatedVolumeAngstromCubed"] - volumes[side]) < 0.01 and
              math.floor(max(0.0, volumes[side] - 0.005 - 1e-8) *
                  0.15 * 6.02214086e23 / 1e27) <=
                  salt_numbers[side] == regions[side]["flooredNominalSaltCount"] <=
                  math.floor((volumes[side] + 0.005 + 1e-8) *
                      0.15 * 6.02214086e23 / 1e27) and
              all(logged_ions[(side, ion)] == planned_ions[(side, ion)] ==
                  regions[side]["providerGenerated" + species + "Count"]
                  for ion, species in (("Na+.pdb", "Sodium"), ("Cl-.pdb", "Chloride")))
              for side in volumes))
    # Memgen branches on its residue-name estimate plus the integer
    # --charge_pdb_delta. LEaP's formally integral charge can carry noise.
    provider_branch_charge = residue_estimate + int(round(supplied_delta))
    check("providerChargeCompensationBranch", all(value is not None for value in salt_numbers.values())
          and conditions["saltBranch"] == ("chargeCompensated" if
          abs(provider_branch_charge) / 2 < min(salt_numbers.values()) else "neutralizationOnly"))
    volume = sum(value for value in volumes.values() if value is not None)
    molar_factor = 6.02214076e23 / 1e27
    check("reportedSaltArithmeticFromProviderVolumes", all(value is not None for value in volumes.values())
          and volume > 0 and conditions["finalWaterCount"] > 0 and
          abs(conditions["sodiumAqueousMolar"] - conditions["finalSodiumCount"] /
              (molar_factor * volume)) < 1e-8 and
          abs(conditions["chlorideAqueousMolar"] - conditions["finalChlorideCount"] /
              (molar_factor * volume)) < 1e-8 and
          abs(conditions["sodiumFiniteWaterMolar"] - 55.4 * conditions["finalSodiumCount"] /
              conditions["finalWaterCount"]) < 1e-8 and
          abs(conditions["chlorideFiniteWaterMolar"] - 55.4 * conditions["finalChlorideCount"] /
              conditions["finalWaterCount"]) < 1e-8)
    report["unavailableCoverage"].append(
        "Provider aqueous-volume geometry is checked against its hashed log, not independently reintegrated")
    report["measurements"].update({"amberAtomCount": len(final.positions), "cellAngstrom": cell,
        "netChargeElementary": charge, "cmapTerms": raw_cmap,
        "finalLipids": {f"{a}:{b}": lipid_counts[(a, b)] for a, b in selected},
        "generatedWater": water, "generatedSodium": sodium, "generatedChloride": chloride,
        "maximumRestartStateDeviationAngstrom": max_position_delta})
    report["unavailableCoverage"].extend([
        "Full stereochemistry and periodic face/edge/corner contact oracle requires separate check",
    ])
    if minimization is not None or (case / "direct-minimization/construction-handoff.json").is_file():
        inspect_minimized((minimization or case / "direct-minimization").resolve(),
                          case, result, report, check)
        report["unavailableCoverage"].append(
            "Direct final worker result does not establish a host-owned completed stage, assessment or export")
    else:
        report["unavailableCoverage"].append(
            "A construct-only exchange cannot establish final minimization, assessment or export")
    report["standing"] = "failed" if report["failures"] else \
        "checksIncomplete" if report["unavailableCoverage"] else "checksPassed"
    return report


def main() -> int:
    if len(sys.argv) not in (2, 3):
        raise SystemExit("Usage: inspect_memgen_route.py CASE_DIRECTORY [MINIMIZATION_DIRECTORY]")
    try:
        report = inspect(Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) == 3 else None)
    except Exception as error:
        report = {"case": sys.argv[1], "standing": "failed", "checks": {},
                  "measurements": {}, "unavailableCoverage": [],
                  "failures": [f"inspector:{type(error).__name__}: {error}"]}
    print(json.dumps(report, indent=2, sort_keys=True))
    return 1 if report["standing"] == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
