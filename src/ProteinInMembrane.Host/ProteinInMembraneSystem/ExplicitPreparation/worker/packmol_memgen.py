"""Bounded crossing to the complete, locally installed PACKMOL-Memgen method.

This module translates exact product inputs and observes provider output.  The
provider owns population, packing, hydration, ions, LEaP and preliminary Amber
conditioning.  No product-side membrane-building substitute is implemented.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import tarfile
import time
from typing import Any, Callable

from ProteinInMembraneSystem.worker.exchange import (
    WorkError, artifact, require_integer, require_mapping, require_number,
    require_text, sha256, verify_sha256, work_path,
)
from ProteinInMembraneSystem.worker.parameterized_structure import (
    _bond_indices, _full_atom_sequence, _nonbonded_charge, _provider,
    _system_bonds_match, _topology_data,
)

_SPECIES = frozenset({"DLPC", "DLPE", "DMPC", "DOPC", "DPPC", "POPC", "POPE", "CHL1"})
_LIPID_HEAD = {name: "P31" for name in _SPECIES if name != "CHL1"} | {"CHL1": "O1"}
_SOURCE_HASHES = {
    "main.py": "7892dd3ad13028ffc572ed84085c8cfef8a323f1dab9b032bdad7fe8a55ceeea",
    "amber.py": "426299f36b88c70eece307ee0a9d4ea5d5a71b679ea12dd2463f0c5ca9a2fdd7",
    "utils.py": "c36d2784ec7957111059c887ea002ea68bd90bfb25377094fc4f0793b102f395",
}
_PROVIDER_CHARGED_RESIDUES = {
    "ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1, "HIP": 1,
    "Cl-": -1, "MG": 2, "Na+": 1, "CA": 2, "OHE": -0.308100,
    "A": -1, "A5": -0.3081, "A3": -0.6919, "DA": -1,
    "DA5": -0.3079, "DA3": -0.6921, "C": -1, "C5": -0.3081,
    "C3": -0.6919, "DC": -1, "DC5": -0.3079, "DC3": -0.6921,
    "G": -1, "G5": -0.3081, "G3": -0.6919, "DG": -1,
    "DG5": -0.3079, "DG3": -0.6921, "U": -1, "U5": -0.3081,
    "U3": -0.6919, "DT": -1, "DT5": -0.3079, "DT3": -0.6921,
    "PTR": -2, "SEP": -2, "TPO": -2, "Y1P": -1, "S1P": -1, "T1P": -1,
    "H1D": 0, "H2D": -1, "H1E": 0, "H2E": -1,
    "NME": 1, "ACE": -1,
}
_LEAP_FILES = frozenset({
    "parm19.dat", "frcmod.ff19SB", "lipid21.dat", "frcmod.tip3p",
    "frcmod.ions1lm_126_tip3p", "frcmod.ionsjc_tip3p",
    "frcmod.ions234lm_126_tip3p", "amino19.lib", "aminoct12.lib",
    "aminont12.lib", "atomic_ions.lib", "lipid21.lib", "solvents.lib",
})


def _asset_paths(payload: dict[str, Any], amberhome: Path) -> dict[str, Path]:
    raw = payload.get("providerAssets")
    if not isinstance(raw, list) or not raw:
        raise WorkError("providerUnavailable", "The complete Memgen route lacks identified provider assets")
    result: dict[str, Path] = {}
    for item in raw:
        asset = require_mapping(item, "provider asset")
        key = require_text(asset.get("id"), "provider asset id")
        if key in result:
            raise WorkError("invalidRequest", "A provider asset identity is repeated")
        require_text(asset.get("version"), "provider asset version")
        path = Path(require_text(asset.get("path"), "provider asset path")).resolve()
        if not path.is_file():
            raise WorkError("providerUnavailable", f"Identified provider asset {key} is unavailable")
        verify_sha256(path, asset.get("sha256"), f"provider asset {key}")
        result[key] = path
    packaged = amberhome / "lib/python3.12/site-packages/packmol_memgen"
    for name, expected in _SOURCE_HASHES.items():
        path = packaged / (name if name == "main.py" else f"lib/{name}")
        if not path.is_file() or sha256(path) != expected:
            raise WorkError("providerMismatch", f"Reviewed Memgen source {name} changed")
        if path not in result.values():
            raise WorkError("unsupportedPolicy", f"Memgen source {name} is absent from the bound asset list")
    for name in ("memgen.parm", "pdbs.tar.gz"):
        path = packaged / "data" / name
        if path not in result.values():
            raise WorkError("unsupportedPolicy", f"Memgen {name} is absent from the bound asset list")
    for name in ("packmol-memgen", "packmol", "tleap", "sander", "ambpdb"):
        path = amberhome / "bin" / name
        if not path.is_file() or path not in result.values():
            raise WorkError("providerUnavailable", f"Bound AmberTools executable {name} is unavailable")
    return result


def _settings(payload: dict[str, Any]) -> dict[str, Any]:
    settings = require_mapping(payload.get("memgen"), "Memgen construction settings")
    exact = {
        "engine": "sander", "proteinForceField": "ff19SB", "lipidForceField": "lipid21",
        "waterForceField": "tip3p", "lateralPaddingAngstrom": 15.0,
        "aqueousPaddingAngstrom": 17.5, "leafletEnvelopeAngstrom": 23.0,
        "localPackingLoops": 20, "totalPackingLoops": 100,
        "packingOptimizerIterations": 20, "packingToleranceAngstrom": 2.0,
        "conditioningSteepestDescentSteps": 250,
        "conditioningConjugateGradientSteps": 250,
        "conditioningRestraintKcalMolAngstromSquared": 10.0,
        "maximumGeometryRetries": 2,
    }
    for key, value in exact.items():
        if settings.get(key) != value:
            raise WorkError("unsupportedPolicy", f"Memgen {key} differs from the selected method")
    for key in ("preoriented", "doNotProtonate", "doNotTrim", "doNotCenterXy",
                "retainIntermediates", "useRatio", "usePbc", "useTightBox",
                "useSalt", "parameterize", "condition"):
        if settings.get(key) is not True:
            raise WorkError("unsupportedPolicy", f"Memgen {key} must be selected")
    if (payload.get("route") != "packmolMemgen" or
            payload.get("saltConvention") != "memgenChargeCompensated" or
            payload.get("providerName") != "PACKMOL-Memgen" or
            payload.get("providerVersion") != "2026.3.25" or
            payload.get("positiveIonArgument") != "Na+" or
            payload.get("negativeIonArgument") != "Cl-" or
            payload.get("ionicStrengthMolar") != 0.15 or
            payload.get("nativePatchPath") is not None or
            payload.get("nativePatchSha256") is not None):
        raise WorkError("unsupportedPolicy", "The selected general construction route changed")
    return settings


def _fractions(payload: dict[str, Any]) -> tuple[list[tuple[str, float]], list[tuple[str, float]]]:
    known = {item.get("speciesId"): item for item in payload.get("selectedSpeciesRepresentations", [])
             if isinstance(item, dict)}
    if not known or set(known) - _SPECIES or len(known) != len(payload.get("selectedSpeciesRepresentations", [])):
        raise WorkError("unsupportedPolicy", "Selected molecular representations are incomplete or repeated")
    def read(key: str, side: str) -> list[tuple[str, float]]:
        raw = require_mapping(payload.get(key), key)
        values = raw.get("fractions")
        if raw.get("physicalSide") != side or not isinstance(values, list) or not values:
            raise WorkError("invalidRequest", f"The {side} target leaflet is absent")
        result = []
        for item in values:
            fraction = require_mapping(item, "leaflet fraction")
            species = require_text(fraction.get("speciesId"), "species id")
            value = require_number(fraction.get("fraction"), "species fraction", 0)
            if value > 0:
                if species not in known or species in {name for name, _ in result}:
                    raise WorkError("unsupportedPolicy", "Leaflet has unavailable or repeated selected chemistry")
                result.append((species, value))
        if not result or abs(sum(value for _, value in result) - 1.0) > 1e-9:
            raise WorkError("invalidRequest", "Leaflet fractions do not sum to one")
        return result
    return read("targetLower", "lower"), read("targetUpper", "upper")


def _force_field_assets(payload: dict[str, Any], amberhome: Path) -> list[Path]:
    raw = payload.get("forceFieldFiles")
    if not isinstance(raw, list) or len(raw) != len(_LEAP_FILES):
        raise WorkError("unsupportedPolicy", "Memgen needs its exact LEaP libraries and parameter files")
    paths = []
    for item in raw:
        asset = require_mapping(item, "Amber force-field asset")
        path = Path(require_text(asset.get("path"), "Amber force-field path")).resolve()
        if path.name not in _LEAP_FILES or path.name in {p.name for p in paths} or not path.is_file():
            raise WorkError("providerUnavailable", "A selected LEaP asset is missing or repeated")
        if amberhome / "dat/leap" not in path.parents:
            raise WorkError("unsupportedPolicy", "A LEaP asset is outside the selected AmberTools runtime")
        verify_sha256(path, asset.get("sha256"), f"LEaP asset {path.name}")
        paths.append(path)
    if {p.name for p in paths} != _LEAP_FILES:
        raise WorkError("unsupportedPolicy", "LEaP asset binding omits selected chemistry")
    return paths


def _provider_residue_name_charge(input_path: Path, source: list[dict[str, Any]]) -> float:
    """Reproduce the pinned provider's residue estimate from its exact input.

    The selected input has unique sequential residue numbers, including for
    distinct source insertion codes.  Memgen tracks the PDB residue-number
    column when adding each charged residue's table value.
    """
    emitted = _pdb_records(input_path)
    if len(emitted) != len(source) or any(
            left["name"] != right["name"] or
            left["resname"] != right["resname"] or
            left["chain"] != right["chain"] or
            left["element"] != right["element"] or
            left["xyz"] != right["xyz"]
            for left, right in zip(emitted, source)):
        raise WorkError("inputMismatch", "Charge estimate does not match exact written provider input")
    charge = 0.0
    tracked_residue = None
    for record in emitted:
        name = record["resname"]
        if name in _PROVIDER_CHARGED_RESIDUES and record["resid"] != tracked_residue:
            charge += _PROVIDER_CHARGED_RESIDUES[name]
            tracked_residue = record["resid"]
    return float(charge)


def _source_charge(trial: Path, amberhome: Path, source: list[dict[str, Any]],
                   approved: Any, disulfides: list[tuple[int, int]]) -> tuple[
                       tuple[float, float, float], dict[int, int]]:
    """Measure the selected prepared construct under the same Amber libraries.

    The provider's residue-name estimator is deliberately separate from the
    charge in Amber's actual parameterized topology.  The difference is the
    provider's documented `--charge_pdb_delta` correction.
    """
    from openmm.app import AmberPrmtopFile
    from openmm.app.internal.amber_file_parser import PrmtopLoader

    provider_estimate = _provider_residue_name_charge(trial / "construct.pdb", source)
    lines = ["source leaprc.lipid21", "source leaprc.protein.ff19SB",
             "source leaprc.water.tip3p", "SYS = loadpdb construct.pdb"]
    lines.extend(f"bond SYS.{first}.SG SYS.{second}.SG" for first, second in disulfides)
    lines.extend(["savepdb SYS prepared-charge.pdb",
                  "saveAmberParm SYS prepared-charge.top prepared-charge.crd", "quit"])
    script = trial / "prepared-charge.leap.in"
    script.write_text("\n".join(lines) + "\n", encoding="ascii")
    environment = dict(os.environ, AMBERHOME=str(amberhome))
    environment["PATH"] = str(amberhome / "bin") + os.pathsep + environment.get("PATH", "")
    with (trial / "prepared-charge.leap.log").open("wb") as output:
        result = subprocess.run([str(amberhome / "bin/tleap"), "-f", script.name],
                                cwd=trial, env=environment, stdout=output,
                                stderr=subprocess.STDOUT, check=False)
    top = trial / "prepared-charge.top"
    if result.returncode or not top.is_file() or not top.stat().st_size:
        raise WorkError("unrepresentableInput", "Amber cannot parameterize the exact prepared chemistry")
    parsed = AmberPrmtopFile(str(top)).topology
    preflight = _pdb_records(trial / "prepared-charge.pdb")
    source_to_amber = _coordinate_correspondence(source, preflight)
    parsed_atoms = list(parsed.atoms())
    if (parsed.getNumAtoms() != approved.getNumAtoms() or len(source_to_amber) != len(source) or
            len(preflight) != len(source) or
            any(parsed_atoms[source_to_amber[i]].element != atom.element
                for i, atom in enumerate(approved.atoms())) or
            {tuple(sorted((source_to_amber[a.index], source_to_amber[b.index])))
             for a, b in approved.bonds()} != _bond_indices(parsed)):
        raise WorkError("unrepresentableInput", "Amber preflight changed prepared atoms or bonds")
    charge = float(sum(PrmtopLoader(str(top)).getCharges()))
    delta = charge - provider_estimate
    if (not math.isfinite(charge) or abs(charge - round(charge)) > 1e-4 or
            abs(delta - round(delta)) > 1e-4):
        raise WorkError("unrepresentableInput", "Prepared Amber charge cannot be represented by Memgen's integral correction")
    return (charge, provider_estimate, delta), source_to_amber


def _approved_disulfide_pairs(approved: Any) -> list[tuple[Any, Any]]:
    pairs = []
    occupied = set()
    for first, second in approved.bonds():
        if (first.name == second.name == "SG" and first.residue is not second.residue and
                first.residue.name in {"CYS", "CYX"} and second.residue.name in {"CYS", "CYX"}):
            if first.index in occupied or second.index in occupied:
                raise WorkError("inputMismatch", "An approved sulfur belongs to multiple disulfides")
            occupied.update((first.index, second.index))
            pairs.append((first, second))
    return pairs


def _amber_prepared_names(approved: Any) -> dict[int, tuple[str, str]]:
    """Translate only selected atom names/states that Amber spells differently.

    The retained graph and its explicit hydrogens, rather than a new
    protonation decision, determine each translation.  LEaP's ff19SB
    N-terminal template uses H1/H2/H3; OpenMM writes the same first hydrogen
    as H.  LEaP otherwise defaults an unsuffixed HIS to HIE, which is wrong
    for an already prepared ND1-protonated histidine.  OpenMM writes bonded
    cysteines as CYS; Amber's corresponding no-thiol-hydrogen state is CYX.
    """
    bonds = {tuple(sorted((first.index, second.index))) for first, second in approved.bonds()}
    disulfide_sulfurs = {atom.index for pair in _approved_disulfide_pairs(approved)
                         for atom in pair}
    result: dict[int, tuple[str, str]] = {}
    for residue in approved.residues():
        atoms = list(residue.atoms())
        names = {atom.name: atom for atom in atoms}
        if len(names) != len(atoms):
            raise WorkError("unrepresentableInput", "A prepared residue repeats an atom name")
        sulfur = names.get("SG")
        if sulfur is not None and sulfur.index in disulfide_sulfurs:
            if "HG" in names or any(
                    atom.element.symbol == "H" and
                    tuple(sorted((sulfur.index, atom.index))) in bonds for atom in atoms):
                raise WorkError("unrepresentableInput", "An approved disulfide retains a thiol hydrogen")
            for atom in atoms:
                result[atom.index] = (atom.name, "CYX")
        if residue.name in {"HIS", "HID", "HIE", "HIP"}:
            if not {"ND1", "NE2"}.issubset(names):
                raise WorkError("unrepresentableInput", "Prepared histidine lacks its imidazole nitrogens")
            for proton, nitrogen in (("HD1", "ND1"), ("HE2", "NE2")):
                if proton in names and tuple(sorted((names[proton].index,
                                                     names[nitrogen].index))) not in bonds:
                    raise WorkError("unrepresentableInput", "Histidine proton is not on its selected nitrogen")
            delta, epsilon = "HD1" in names, "HE2" in names
            if not (delta or epsilon):
                raise WorkError("unrepresentableInput", "Prepared histidine has no selected ring proton")
            state = "HIP" if delta and epsilon else ("HID" if delta else "HIE")
            if residue.name != "HIS" and residue.name != state:
                raise WorkError("unrepresentableInput", "Histidine label contradicts retained hydrogen state")
            for atom in atoms:
                result[atom.index] = (atom.name, state)
        if {"N", "H", "H2", "H3"}.issubset(names):
            if "H1" in names or any(tuple(sorted((names["N"].index,
                                                  names[hydrogen].index))) not in bonds
                                    for hydrogen in ("H", "H2", "H3")):
                raise WorkError("unrepresentableInput", "N-terminal hydrogens lack an exact Amber mapping")
            result[names["H"].index] = ("H1", result.get(names["H"].index,
                                                            ("H", residue.name))[1])
    return result


def _input_pdb(directory: Path, trial: Path, payload: dict[str, Any]) -> tuple[Path, list[dict[str, Any]], list[tuple[int, int]]]:
    source = work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")
    verify_sha256(source, payload.get("orientedPdbSha256"), "orientedPdbSha256")
    prepared = work_path(directory, payload.get("preparedPdbPath"), "preparedPdbPath")
    verify_sha256(prepared, payload.get("preparedPdbSha256"), "preparedPdbSha256")
    graph = work_path(directory, payload.get("preparedBondGraphPath"), "preparedBondGraphPath")
    verify_sha256(graph, payload.get("preparedBondGraphSha256"), "preparedBondGraphSha256")
    from openmm.app import PDBFile
    from .construction import _bond_indices, _full_atom_sequence, _read_topology_data
    oriented, original = PDBFile(str(source)), PDBFile(str(prepared))
    approved = _read_topology_data(graph)
    if (_full_atom_sequence(oriented.topology) != _full_atom_sequence(original.topology) or
            _full_atom_sequence(oriented.topology) != _full_atom_sequence(approved) or
            _bond_indices(oriented.topology) != _bond_indices(approved)):
        raise WorkError("inputMismatch", "Oriented construct does not preserve prepared molecular identity")
    correspondence = require_mapping(payload.get("preparedCorrespondence"), "prepared correspondence")
    mapped = correspondence.get("atoms")
    if (correspondence.get("complete") is not True or not isinstance(mapped, list) or
            len(mapped) != oriented.topology.getNumAtoms() or
            sorted(item.get("resultAtomIndex") for item in mapped) != list(range(len(mapped)))):
        raise WorkError("inputMismatch", "Prepared atom correspondence is incomplete")
    records = []
    lines = []
    previous_residue = None
    residue_index = 0
    leap_residue_by_atom = []
    serial_by_atom = []
    approved_disulfides = _approved_disulfide_pairs(approved)
    amber_names = _amber_prepared_names(approved)
    for line in source.read_text(encoding="ascii").splitlines():
        if line.startswith("CRYST1"):
            continue
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 54 or not line[6:11].strip().isdigit() or not line[22:26].strip().lstrip("-").isdigit():
                raise WorkError("unrepresentableInput", "Prepared construct exceeds fixed-width provider input")
            residue_key = (line[21:22], line[22:26], line[26:27])
            if residue_key != previous_residue:
                residue_index += 1
                previous_residue = residue_key
            if residue_index > 9999:
                raise WorkError("unrepresentableInput", "Provider PDB residue numbering would wrap")
            leap_residue_by_atom.append(residue_index)
            serial_by_atom.append(int(line[6:11]))
            name, resname = amber_names.get(len(records),
                                            (line[12:16].strip(), line[17:20].strip()))
            if name != line[12:16].strip():
                line = line[:12] + f" {name:<3}" + line[16:]
            if resname != line[17:20].strip():
                line = line[:17] + f"{resname:>3}" + line[20:]
            records.append({"name": line[12:16].strip(), "resname": line[17:20].strip(),
                            "resid": line[22:26].strip(), "chain": line[21:22],
                            "xyz": tuple(float(line[i:i + 8]) for i in (30, 38, 46)),
                            "element": line[76:78].strip().upper()})
            line = line[:22] + f"{residue_index:4d}" + line[26:]
        lines.append(line)
    if len(records) != len(mapped) or not all(math.isfinite(v) for item in records for v in item["xyz"]):
        raise WorkError("inputMismatch", "Provider input does not contain every finite prepared atom")
    if len(serial_by_atom) != len(set(serial_by_atom)):
        raise WorkError("unrepresentableInput", "Prepared PDB repeats an atom serial")
    duplicate_conect = {frozenset((serial_by_atom[first.index], serial_by_atom[second.index]))
                        for first, second in approved_disulfides}
    provider_lines = []
    for line in lines:
        if line.startswith("CONECT") and duplicate_conect:
            try:
                serials = [int(line[index:index + 5]) for index in range(6, len(line), 5)
                           if line[index:index + 5].strip()]
            except ValueError as error:
                raise WorkError("unrepresentableInput", "Prepared PDB has malformed bond serials") from error
            if not serials:
                raise WorkError("unrepresentableInput", "Prepared PDB has an empty bond record")
            retained = [serial for serial in serials[1:]
                        if frozenset((serials[0], serial)) not in duplicate_conect]
            if not retained:
                continue
            line = "CONECT" + f"{serials[0]:5d}" + "".join(f"{serial:5d}" for serial in retained)
        provider_lines.append(line)
    target = trial / "construct.pdb"
    target.write_text("\n".join(provider_lines) + "\n", encoding="ascii")
    disulfides = [(leap_residue_by_atom[first.index], leap_residue_by_atom[second.index])
                  for first, second in approved_disulfides]
    return target, records, disulfides


def _run_provider(trial: Path, payload: dict[str, Any], settings: dict[str, Any],
                  lower: list[tuple[str, float]], upper: list[tuple[str, float]],
                  amberhome: Path, disulfides: list[tuple[int, int]],
                  charges: tuple[float, float, float], progress: Callable) -> None:
    executable = amberhome / "bin/packmol-memgen"
    options = [str(executable), "--pdb", "construct.pdb", "--lipids",
               ":".join(species for species, _ in lower) + "//" +
               ":".join(species for species, _ in upper), "--ratio",
               ":".join(format(value, ".17g") for _, value in lower) + "//" +
               ":".join(format(value, ".17g") for _, value in upper),
               "--output", "system.pdb", "--log", "packmol-memgen.log",
               "--packlog", "packmol", "--preoriented", "--notprotonate", "--nottrim",
               "--noxy_cen", "--keep", "--pbc", "--tight_box", "--leaflet", "23",
               "--dist", format(require_number(payload.get("lateralPaddingAngstrom"), "lateralPaddingAngstrom"), ".17g"),
               "--dist_wat", format(require_number(payload.get("aqueousPaddingAngstrom"), "aqueousPaddingAngstrom"), ".17g"),
               "--salt", "--salt_c", "Na+", "--salt_a", "Cl-", "--saltcon", "0.15",
               "--ffprot", "ff19SB", "--fflip", "lipid21", "--ffwat", "tip3p",
               "--engine", "sander", "--parametrize", "--minimize",
               "--nloop", "20", "--nloop_all", "100", "--maxit", "20",
               "--tolerance", "2", "--sd_steps", "250", "--cg_steps", "250"]
    formal, estimated, delta = charges
    for key, value in (("preparedFormalChargeElementary", formal),
                       ("providerResidueNameChargeElementary", estimated),
                       ("chargePdbDeltaElementary", delta)):
        supplied = payload.get(key)
        if supplied is not None and abs(require_number(supplied, key) - value) > 1e-4:
            raise WorkError("inputMismatch", f"{key} differs from exact Amber charge observation")
    if abs(formal - estimated - delta) > 1e-4 or abs(delta - round(delta)) > 1e-4:
        raise WorkError("inputMismatch", "Memgen charge correction does not match selected chemistry")
    if abs(delta) > 1e-4:
        options.extend(["--charge_pdb_delta", str(int(round(delta)))])
    for first, second in disulfides:
        if any(not isinstance(value, int) or value < 1 or value > 9999 for value in (first, second)):
            raise WorkError("inputMismatch", "An authorized disulfide lacks a provider residue map")
        options.extend(["--leapline", f"bond SYS.{first}.SG SYS.{second}.SG"])
    (trial / "provider-command.json").write_text(json.dumps(options, indent=2) + "\n", encoding="utf-8")
    environment = dict(os.environ)
    environment["AMBERHOME"] = str(amberhome)
    environment["PATH"] = str(amberhome / "bin") + os.pathsep + environment.get("PATH", "")
    detail = {"trialId": payload.get("trialId"), "trialIndex": payload.get("trialIndex")}
    progress("providerPopulation", {**detail, "message": "Memgen is deriving the trial population and cell."})
    with (trial / "provider-stdout.log").open("wb") as stdout:
        process = subprocess.Popen(options, cwd=trial, env=environment,
                                   stdout=stdout, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            observed = set()
            # Packmol writes system.pdb while its optimizer is still running.
            # Memgen creates leap.in only after packing and postprocessing.
            stages = (("leap.in", "providerCleanup", False),
                      ("leap_.log", "amberParameterization", True),
                      ("min.out", "providerConditioningRestrained", True))
            while (returncode := process.poll()) is None:
                # The generated Packmol input exists before Packmol has begun
                # arranging molecules. Its own first optimizer loop is the
                # attributable start of packing.
                packlog = trial / "packmol.log"
                if "providerPacking" not in observed and packlog.is_file():
                    with packlog.open("r", encoding="utf-8", errors="replace") as log:
                        packing_started = any("Starting GENCAN loop:" in line for line in log)
                    if packing_started:
                        observed.add("providerPacking")
                        progress("providerPacking", {**detail,
                                                     "message": "PACKMOL began molecular packing optimization."})
                for filename, stage, nonempty in stages:
                    marker = trial / filename
                    if (stage not in observed and marker.is_file() and
                            (not nonempty or marker.stat().st_size > 0)):
                        observed.add(stage)
                        progress(stage, {**detail, "message": f"Memgen produced {filename}."})
                if ("providerConditioningUnrestrained" not in observed and
                        any(path.stat().st_size > 0 for path in trial.glob("*_min.out"))):
                    observed.add("providerConditioningUnrestrained")
                    progress("providerConditioningUnrestrained",
                             {**detail, "message": "Memgen produced its final unrestrained output."})
                time.sleep(0.5)
        except BaseException:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
            raise
    if returncode != 0:
        log = "\n".join((trial / name).read_text(encoding="utf-8", errors="replace")[-50000:]
                        for name in ("packmol-memgen.log", "provider-stdout.log") if (trial / name).is_file())
        if "concentration of ions required to neutralize" in log:
            code = "chargeBudgetRefusal"
        elif "ratio for lipid" in log and "too small" in log:
            code = "zeroRoundedSpecies"
        else:
            code = "providerFailed"
        details: dict[str, Any] = {}
        recipe = trial / "packmol.inp"
        if code == "providerFailed" and recipe.is_file():
            pbc = re.findall(r"^pbc\s+([-+\d. ]+)\s*$", recipe.read_text(), flags=re.MULTILINE)
            if len(pbc) == 1 and len(pbc[0].split()) == 6:
                bounds = list(map(float, pbc[0].split()))
                lengths = [bounds[i + 3] - bounds[i] for i in range(3)]
                source = _pdb_records(trial / "construct.pdb")
                spans = [max(item["xyz"][i] for item in source) -
                         min(item["xyz"][i] for item in source) for i in range(3)]
                required = [2 * payload["lateralPaddingAngstrom"]] * 2 + [
                    2 * payload["aqueousPaddingAngstrom"]]
                for axis, (length, span, clearance) in enumerate(zip(lengths, spans, required)):
                    if length - span < clearance:
                        code = "insufficientCellClearance"
                        details["axis"] = "aqueous" if axis == 2 else "lateral"
                        break
        raise WorkError(code, "Full Memgen workflow failed; inspect exact trial diagnostics", details)


def _pdb_records(path: Path) -> list[dict[str, Any]]:
    records = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise WorkError("malformedProviderArtifact", f"{path.name} has a truncated atom record")
        try:
            xyz = tuple(float(line[index:index + 8]) for index in (30, 38, 46))
        except ValueError as error:
            raise WorkError("malformedProviderArtifact", f"{path.name} has a nonnumeric coordinate") from error
        if not all(math.isfinite(value) for value in xyz):
            raise WorkError("nonfiniteObservation", f"{path.name} has a nonfinite coordinate")
        records.append({"name": line[12:16].strip(), "resname": line[17:20].strip(),
                        "resid": line[22:26].strip(), "chain": line[21:22],
                        "xyz": xyz, "element": line[76:78].strip().upper()})
    if not records:
        raise WorkError("malformedProviderArtifact", f"{path.name} contains no atoms")
    return records


def _same_records(before: list[dict[str, Any]], after: list[dict[str, Any]],
                  *, position_tolerance: float = 0.002) -> bool:
    def normalized(item: dict[str, Any]) -> str:
        if item["resname"] == "ILE" and item["name"] == "CD":
            return "CD1"
        return {"OT1": "O", "OT2": "OXT"}.get(item["name"], item["name"])
    return len(before) == len(after) and all(
        normalized(first) == second["name"] and
        first["resname"] == second["resname"] and
        math.dist(first["xyz"], second["xyz"]) <= position_tolerance
        for first, second in zip(before, after))


def _template_records(archive: Path, species: str) -> list[dict[str, Any]]:
    with tarfile.open(archive, "r:gz") as handle:
        try:
            member = handle.extractfile(f"{species}.pdb")
        except KeyError as error:
            raise WorkError("providerMismatch", f"Memgen lacks a {species} coordinate template") from error
        if member is None:
            raise WorkError("providerMismatch", f"Memgen cannot read its {species} coordinate template")
        raw = member.read().decode("ascii").splitlines()
    # The archive is immutable and verified. Parse without writing a second
    # product template or depending on a provider extraction side effect.
    result = []
    for line in raw:
        if line.startswith(("ATOM  ", "HETATM")):
            result.append({"name": line[12:16].strip(), "resname": line[17:20].strip()})
    if not result:
        raise WorkError("providerMismatch", f"Memgen {species} template has no atoms")
    return result


def _packmol_recipe(path: Path) -> list[dict[str, Any]]:
    """Read actual provider-authored populations and leaflet constraints."""
    blocks = []
    current: dict[str, Any] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        words = line.strip().split()
        if len(words) == 2 and words[0] == "structure":
            current = {"file": Path(words[1]).name, "count": None, "side": None}
        elif current is not None and len(words) == 2 and words[0] == "number":
            current["count"] = int(words[1])
        elif current is not None and len(words) > 1 and words[0] in {"below", "over", "above"} and words[1] == "plane":
            side = "lower" if words[0] == "below" else "upper"
            if current["side"] is not None and current["side"] != side and \
                    Path(current["file"]).stem not in _SPECIES:
                raise WorkError("malformedProviderArtifact", "A Packmol block crosses physical sides")
            if current["side"] is None:
                current["side"] = side
        elif words[:2] == ["end", "structure"]:
            if current is None or current["count"] is None or current["count"] < 0:
                raise WorkError("malformedProviderArtifact", "Packmol block lacks a population")
            blocks.append(current)
            current = None
    if current is not None or not blocks:
        raise WorkError("malformedProviderArtifact", "Packmol recipe is incomplete")
    return blocks


def _packed_molecules(records: list[dict[str, Any]], source: list[dict[str, Any]],
                      templates: dict[str, list[dict[str, Any]]],
                      recipe: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not _same_records(source, records[:len(source)]):
        raise WorkError("correspondenceFailed", "Packing changed the fixed retained construct")
    cursor = len(source)
    molecules = []
    for block in recipe:
        if block["file"] == "PROT0.pdb":
            if block["count"] != 1:
                raise WorkError("correspondenceFailed", "Provider did not fix one complete construct")
            continue
        species = Path(block["file"]).stem
        template = templates.get(species)
        if template is None:
            raise WorkError("correspondenceFailed", f"Packmol recipe introduces unselected {species}")
        for _ in range(block["count"]):
            part = records[cursor:cursor + len(template)]
            if len(part) != len(template) or any(
                    atom["name"] != expected["name"] or atom["resname"] != expected["resname"]
                    for atom, expected in zip(part, template)):
                raise WorkError("correspondenceFailed", f"Packed {species} differs from bound provider template")
            anchor = _LIPID_HEAD.get(species, template[0]["name"] if len(template) == 1 else "O")
            heads = [atom["xyz"][2] for atom in part if atom["name"] == anchor]
            if len(heads) != 1 or heads[0] == 0:
                raise WorkError("correspondenceFailed", f"{species} has no unique physical-side anchor")
            side = "upper" if heads[0] > 0 else "lower"
            if species in _SPECIES and block["side"] != side:
                raise WorkError("correspondenceFailed", f"Packed {species} moved to the wrong leaflet")
            molecules.append({"species": species, "side": side, "start": cursor,
                              "atoms": part, "headZ": heads[0]})
            cursor += len(part)
    if cursor != len(records):
        raise WorkError("correspondenceFailed", "Packed output has unaccounted atoms")
    return molecules


def _provider_artifacts(directory: Path, trial: Path, packed_base: str,
                        *, include_final: bool) -> list[dict[str, str]]:
    named = {
        "providerOptions": "packmol-memgen.json", "providerLog": "packmol-memgen.log",
        "providerCommand": "provider-command.json", "providerStdout": "provider-stdout.log",
        "providerInput": "construct.pdb",
        "preparedChargeLeapInput": "prepared-charge.leap.in",
        "preparedChargeLeapLog": "prepared-charge.leap.log",
        "preparedChargeLeapTopology": "prepared-charge.top",
        "preparedChargeLeapCoordinates": "prepared-charge.crd",
        "preparedChargePdb": "prepared-charge.pdb",
        "providerProtein": "PROT0.pdb", "providerPackmolInput": "packmol.inp",
        "providerPackmolLog": "packmol.log", "providerPacked": "system.pdb",
        "providerLeapInput": "leap.in", "providerLeapLog": "leap_.log",
        "providerParameterized": f"{packed_base}_lipid.pdb",
        "amberCoordinates": f"{packed_base}_lipid.crd",
        "amberTopology": f"{packed_base}_lipid.top",
        "providerRestrainedInput": "min.in", "providerRestrainedLog": "min.out",
        "providerRestrainedRestart": "min.restrt",
        "providerUnrestrainedInput": f"{packed_base}_min.in",
        "providerUnrestrainedLog": f"{packed_base}_min.out",
        "amberFinalRestart": f"{packed_base}_min.restrt",
    }
    if packed_base != "system":
        named["providerPostCleanup"] = f"{packed_base}.pdb"
    found = [artifact(directory, trial / name, role) for role, name in named.items()
             if (trial / name).is_file()]
    if include_final and any(not (trial / name).is_file() for role, name in named.items()
                             if role != "providerPostCleanup"):
        raise WorkError("unobservedOutput", "Memgen did not return every required method artifact")
    return found


def _coordinate_correspondence(before: list[dict[str, Any]],
                               after: list[dict[str, Any]], *, tolerance: float = 0.025
                               ) -> dict[int, int]:
    """Associate preconditioning PDB atoms by chemical name and coordinates.

    No nearest-neighbor identity is accepted beyond PDB rounding.  LEaP names
    retained HOH as WAT.  It can also place a replacement ion at a removed
    water atom's exact coordinate; different chemical names must leave that
    pair explicitly unmatched.
    """
    scale = 1 / tolerance
    buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for index, item in enumerate(after):
        buckets[tuple(math.floor(value * scale) for value in item["xyz"])].append(index)
    result = {}
    used = set()
    for index, item in enumerate(before):
        key = tuple(math.floor(value * scale) for value in item["xyz"])
        matches = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for candidate in buckets.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if (candidate not in used and
                                item["name"] == after[candidate]["name"] and
                                (item["resname"] == after[candidate]["resname"] or
                                 (item["resname"], after[candidate]["resname"]) == ("HOH", "WAT")) and
                                math.dist(item["xyz"], after[candidate]["xyz"]) <= tolerance):
                            matches.append(candidate)
        if len(matches) > 1:
            raise WorkError("correspondenceFailed", "Provider coordinates do not give unique atom identity")
        if matches:
            result[index] = matches[0]
            used.add(matches[0])
    return result


def _components(topology: Any) -> tuple[list[set[int]], dict[int, int]]:
    parent = list(range(topology.getNumAtoms()))
    def find(index: int) -> int:
        while index != parent[index]:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    for first, second in topology.bonds():
        a, b = find(first.index), find(second.index)
        if a != b:
            parent[b] = a
    groups: dict[int, set[int]] = defaultdict(set)
    for index in range(len(parent)):
        groups[find(index)].add(index)
    result = list(groups.values())
    return result, {atom: group for group, members in enumerate(result) for atom in members}


def _verified_molecular_graph(actual: Any, actual_indices: set[int],
                              reference: Any, species: str) -> dict[int, int]:
    from collections import Counter as Multiset
    left = list(actual.atoms())
    right = list(reference.atoms())
    if len(actual_indices) != len(right) or Multiset(left[i].element.symbol for i in actual_indices) != \
            Multiset(atom.element.symbol for atom in right):
        raise WorkError("chemistryMismatch", f"Parameterized {species} has a changed molecular formula")
    actual_edges = {tuple(sorted((a.index, b.index))) for a, b in actual.bonds()
                    if a.index in actual_indices and b.index in actual_indices}
    reference_edges = {tuple(sorted((a.index, b.index))) for a, b in reference.bonds()}
    actual_adjacency = {i: set() for i in actual_indices}
    reference_adjacency = {i: set() for i in range(len(right))}
    for a, b in actual_edges:
        actual_adjacency[a].add(b); actual_adjacency[b].add(a)
    for a, b in reference_edges:
        reference_adjacency[a].add(b); reference_adjacency[b].add(a)
    def colours(adjacency: dict[int, set[int]], atoms: list[Any], indices: set[int]) -> dict[int, str]:
        labels = {i: atoms[i].element.symbol + ":" + str(len(adjacency[i])) for i in indices}
        for _ in range(8):
            labels = {i: __import__("hashlib").sha256(
                (labels[i] + "|" + "|".join(sorted(labels[j] for j in adjacency[i]))).encode()).hexdigest()
                for i in indices}
        return labels
    ac = colours(actual_adjacency, left, actual_indices)
    rc = colours(reference_adjacency, right, set(range(len(right))))
    if Multiset(ac.values()) != Multiset(rc.values()) or len(actual_edges) != len(reference_edges):
        raise WorkError("chemistryMismatch", f"Parameterized {species} differs from selected bonded chemistry")
    classes: dict[str, list[int]] = defaultdict(list)
    for index, colour in ac.items():
        classes[colour].append(index)
    mapping = {}
    used = set()
    order = sorted(range(len(right)), key=lambda i: (len(classes[rc[i]]), -len(reference_adjacency[i])))
    def match(position: int) -> bool:
        if position == len(order):
            return True
        index = order[position]
        for candidate in classes[rc[index]]:
            if candidate in used or any((other in reference_adjacency[index]) !=
                                         (mapping[other] in actual_adjacency[candidate]) for other in mapping):
                continue
            mapping[index] = candidate; used.add(candidate)
            if match(position + 1):
                return True
            used.remove(candidate); del mapping[index]
        return False
    if not match(0):
        raise WorkError("chemistryMismatch", f"Parameterized {species} molecular graph changed")
    return mapping


def _verify_stereo(reference: Any, mapping: dict[int, int], positions: Any,
                   checks: Any, species: str) -> None:
    from openmm import unit
    from .construction import _native_alkene_cosine, _native_signed_volume

    atoms = list(reference.atoms())
    if not isinstance(checks, list) or not checks:
        raise WorkError("invalidTemplate", f"Selected {species} lacks its stereo checks")
    points = [tuple(float(value) for value in positions[mapping[i]].value_in_unit(unit.angstrom))
              for i in range(len(atoms))]
    names = {atom.name: atom.index for atom in atoms}
    bonds = _bond_indices(reference)
    for check in checks:
        if check.get("kind") == "tetrahedral":
            _native_signed_volume(points, names, check, species)
        elif check.get("kind") == "alkene":
            _native_alkene_cosine(points, names, bonds, check, species)
        else:
            raise WorkError("invalidTemplate", f"Selected {species} has an unsupported stereo descriptor")


def _log_number(log: str, label: str) -> float:
    found = re.findall(rf"{re.escape(label)}\s*=\s*([-+]?\d+(?:\.\d+)?)", log)
    if len(found) != 1:
        raise WorkError("unobservedCondition", f"Memgen did not record {label}")
    return float(found[0])


def _provider_salt_branch(residue_name_charge: float, charge_pdb_delta: float,
                          aqueous: list[dict[str, Any]]) -> str:
    """Report Memgen's branch using its residue estimate and integer correction.

    LEaP's measured formal charge can carry floating noise at an integer
    boundary. Memgen branches on its own corrected charge instead.
    """
    provider_charge = residue_name_charge + int(round(charge_pdb_delta))
    return ("chargeCompensated" if abs(provider_charge) / 2 <
            min(item["flooredNominalSaltCount"] for item in aqueous)
            else "neutralizationOnly")


def _logged_salt_count_possible(volume: float, count: int) -> bool:
    """Account for Memgen's 0.01 A^3 volume rounding before its salt floor."""
    factor = 0.15 * 6.02214086e23 / 1e27
    low = math.floor(max(0.0, volume - 0.005 - 1e-8) * factor)
    high = math.floor((volume + 0.005 + 1e-8) * factor)
    return low <= count <= high


def _trial_summary(payload: dict[str, Any], standing: str, proposed: list[dict[str, Any]],
                   actual: list[dict[str, Any]], removed: list[dict[str, Any]],
                   proposed_cell: list[float], actual_cell: list[float],
                   conditions: dict[str, Any] | None, failure: str | None = None,
                   message: str | None = None) -> dict[str, Any]:
    return {"trialId": require_text(payload.get("trialId"), "trialId"),
            "trialIndex": require_integer(payload.get("trialIndex"), "trialIndex"),
            "standing": standing,
            "lateralPaddingAngstrom": require_number(payload.get("lateralPaddingAngstrom"), "lateralPaddingAngstrom"),
            "aqueousPaddingAngstrom": require_number(payload.get("aqueousPaddingAngstrom"), "aqueousPaddingAngstrom"),
            "proposedLipidCounts": proposed, "achievedLipidCounts": actual,
            "cleanupRemovedLipidCounts": removed,
            "proposedCellAngstrom": proposed_cell, "actualCellAngstrom": actual_cell,
            "conditions": conditions, "failureCode": failure, "message": message}


def _periodic_offset(value: float, reference: float, length: float) -> float:
    difference = value - reference
    return (difference + length / 2) % length - length / 2


def _conditioned_retained_gaps(points: list[tuple[float, float, float]],
                               mapped_source: dict[int, int],
                               source: list[dict[str, Any]],
                               cell: list[float],
                               amber_shift: tuple[float, float, float]) -> list[float]:
    """Measure retained-construct image gaps in the conditioned Amber state.

    Map each final atom back to the selected source image before taking its
    extent.  This preserves the original construct's relation across a box
    face when the restart wraps one of its atoms into the opposite image.
    """
    if len(mapped_source) != len(source) or set(mapped_source) != set(range(len(source))):
        raise WorkError("correspondenceFailed", "Conditioned retained atom mapping is incomplete")
    if any(not math.isfinite(length) or length <= 0 for length in cell):
        raise WorkError("invalidGeometry", "Conditioned periodic cell is invalid")
    unwrapped: list[tuple[float, float, float]] = []
    for index, record in enumerate(source):
        final = mapped_source[index]
        if final < 0 or final >= len(points):
            raise WorkError("correspondenceFailed", "Conditioned retained atom index is invalid")
        coordinate = tuple(record["xyz"][axis] + _periodic_offset(
            points[final][axis] - amber_shift[axis], record["xyz"][axis], cell[axis])
            for axis in range(3))
        if not all(math.isfinite(value) for value in coordinate):
            raise WorkError("nonfiniteObservation", "Conditioned retained coordinate is nonfinite")
        unwrapped.append(coordinate)
    gaps = [cell[axis] - (max(point[axis] for point in unwrapped) -
                          min(point[axis] for point in unwrapped)) for axis in range(3)]
    for axis, value in enumerate(gaps):
        if value < 0:
            raise WorkError("insufficientCellClearance",
                            "Conditioned retained construct intersects its periodic image",
                            {"axis": "aqueous" if axis == 2 else "lateral"})
    return gaps


def _amber_local_observation_spec(payload: dict[str, Any], amber_shift_z: float) -> tuple[dict[str, Any], float]:
    """Express the selected PDB-frame midplane in the verified Amber frame."""
    selected = require_mapping(payload.get("localObservationSpec"), "localObservationSpec")
    pdb_midplane = 10 * require_number(payload.get("membraneCenterZNanometers"),
                                        "membraneCenterZNanometers")
    specified = require_number(selected.get("referenceMidplaneZAngstrom"),
                               "referenceMidplaneZAngstrom")
    if abs(specified - pdb_midplane) > 1e-8:
        raise WorkError("inputMismatch", "Selected local midplane differs from the prepared PDB frame")
    output_midplane = pdb_midplane + amber_shift_z
    return {**selected, "referenceMidplaneZAngstrom": output_midplane}, output_midplane


def _verify_amber_pdb_order(top: Path, initial_coordinates: Path,
                            pdb: list[dict[str, Any]], topology: Any) -> tuple[tuple[float, float, float], float]:
    """Bind LEaP's PDB and initial NetCDF order to raw and imported Amber atoms.

    OpenMM standardizes some Amber atom names (N-terminal H1 becomes H), so
    its displayed topology names cannot be compared literally with LEaP PDB.
    The prmtop's original names and atomic numbers are the identity source.
    """
    from openmm import unit
    from openmm.app import AmberInpcrdFile, PDBFile, element
    from openmm.app.internal.amber_file_parser import PrmtopLoader

    raw = PrmtopLoader(str(top))
    initial = AmberInpcrdFile(str(initial_coordinates))
    atoms = list(topology.atoms())
    if (not raw.has_atomic_number or len(pdb) != raw.getNumAtoms() or
            len(initial.positions) != len(pdb) or len(atoms) != len(pdb)):
        raise WorkError("correspondenceFailed", "LEaP Amber/PDB atom inventory is incomplete")
    PDBFile._loadNameReplacementTables()
    numbers = raw._raw_data["ATOMIC_NUMBER"]
    coordinates = initial.positions.value_in_unit(unit.angstrom)
    translation = [coordinates[0][axis] - pdb[0]["xyz"][axis] for axis in range(3)]
    maximum_residual = 0.0
    for index, (record, atom) in enumerate(zip(pdb, atoms)):
        raw_name = raw.getAtomName(index).strip()
        raw_residue = raw.getResidueLabel(index).strip()
        normalized_residue = PDBFile._residueNameReplacements.get(raw_residue, raw_residue)
        normalized_name = PDBFile._atomNameReplacements.get(normalized_residue, {}).get(
            raw_name, raw_name)
        symbol = element.Element.getByAtomicNumber(int(numbers[index])).symbol
        if (record["name"] != raw_name or record["resname"] != raw_residue or
                record["element"] != symbol.upper() or atom.name != normalized_name or
                atom.residue.name != normalized_residue or atom.element.symbol != symbol):
            raise WorkError("correspondenceFailed", "LEaP Amber/PDB atom identity or order differs")
        residual = math.dist(
            coordinates[index],
            [record["xyz"][axis] + translation[axis] for axis in range(3)])
        maximum_residual = max(maximum_residual, residual)
        if residual > 0.002:
            raise WorkError("correspondenceFailed", "LEaP Amber/PDB coordinate order differs")
    return (tuple(float(value) for value in translation), maximum_residual)


def _amber_import(top: Path, restart: Path, payload: dict[str, Any]) -> tuple[Any, Any, Any, dict[str, Any]]:
    """Compare actual Amber terms and box with the selected OpenMM import."""
    import hashlib
    from openmm import (CMAPTorsionForce, CMMotionRemover, HarmonicAngleForce,
                        HarmonicBondForce, NonbondedForce, PeriodicTorsionForce, unit)
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME
    from openmm.app.internal.amber_file_parser import PrmtopLoader

    raw = PrmtopLoader(str(top))
    coordinate = AmberInpcrdFile(str(restart))
    if coordinate.boxVectors is None or len(coordinate.positions) != raw.getNumAtoms():
        raise WorkError("unobservedOutput", "Final Amber restart lacks the topology's periodic atoms")
    amber = AmberPrmtopFile(str(top), periodicBoxVectors=coordinate.boxVectors)
    settings = require_mapping(payload.get("systemSettings"), "system settings")
    if (settings.get("nonbondedMethod") != "PME" or settings.get("constraints") != "HBonds" or
            settings.get("rigidWater") is not True or settings.get("nonbondedCutoffNanometers") != 1.0 or
            settings.get("ewaldErrorTolerance") != 0.0005 or
            settings.get("useDispersionCorrection") is not True or
            settings.get("switchDistanceNanometers") is not None or
            settings.get("hydrogenMassDaltons") is not None):
        raise WorkError("unsupportedPolicy", "Amber import settings differ from the accepted final method")
    remove_cm = settings.get("removeCMMotion")
    if not isinstance(remove_cm, bool):
        raise WorkError("invalidRequest", "CMMotion removal setting must be explicit")
    system = amber.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                                constraints=HBonds, rigidWater=True,
                                ewaldErrorTolerance=0.0005, removeCMMotion=remove_cm,
                                hydrogenMass=None)
    topology = amber.topology
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    kinds = sorted(type(force).__name__ for force in forces)
    allowed = {"CMAPTorsionForce", "CMMotionRemover", "HarmonicAngleForce",
               "HarmonicBondForce", "NonbondedForce", "PeriodicTorsionForce"}
    if (set(kinds) - allowed or kinds.count("NonbondedForce") != 1 or
            kinds.count("HarmonicAngleForce") != 1 or
            kinds.count("PeriodicTorsionForce") != 1 or
            kinds.count("HarmonicBondForce") != 1 or
            kinds.count("CMMotionRemover") != int(remove_cm)):
        raise WorkError("parameterMismatch", "Imported Amber force inventory is incomplete or unsupported")
    nonbonded = next(force for force in forces if isinstance(force, NonbondedForce))
    nonbonded.setUseDispersionCorrection(True)
    if (nonbonded.getNonbondedMethod() != NonbondedForce.PME or
            abs(nonbonded.getCutoffDistance().value_in_unit(unit.nanometer) - 1) > 1e-10 or
            abs(nonbonded.getEwaldErrorTolerance() - 0.0005) > 1e-12 or
            not nonbonded.getUseDispersionCorrection() or
            nonbonded.getUseSwitchingFunction()):
        raise WorkError("parameterMismatch", "Imported Amber nonbonded settings changed")
    n = raw.getNumAtoms()
    if (n != topology.getNumAtoms() or n != system.getNumParticles() or
            nonbonded.getNumParticles() != n):
        raise WorkError("parameterMismatch", "Amber atom order/count changed on import")
    for index, mass in enumerate(raw.getMasses()):
        if abs(system.getParticleMass(index).value_in_unit(unit.dalton) - mass) > 1e-6:
            raise WorkError("parameterMismatch", "Amber particle mass changed on import")
    raw_charges = raw.getCharges()
    raw_lj = raw.getNonbondTerms()
    converted = []
    for index, (charge, (rvdw, epsilon)) in enumerate(zip(raw_charges, raw_lj)):
        found_charge, found_sigma, found_epsilon = nonbonded.getParticleParameters(index)
        observed = (found_charge.value_in_unit(unit.elementary_charge),
                    found_sigma.value_in_unit(unit.nanometer),
                    found_epsilon.value_in_unit(unit.kilojoule_per_mole))
        expected = (charge, rvdw * 2 ** (5 / 6), epsilon)
        if any(abs(a - b) > 1e-6 * max(1.0, abs(a)) for a, b in zip(expected, observed)):
            raise WorkError("parameterMismatch", "Amber charge or Lennard-Jones term changed on import")
        converted.append(observed)
    amber_charge = float(sum(raw_charges))
    imported_charge = float(sum(item[0] for item in converted))
    if abs(imported_charge - amber_charge) > 1e-6:
        raise WorkError("parameterMismatch", "Amber net charge changed on import")
    # Compare all bonded terms against the independent prmtop arrays, including
    # terms replaced by HBond/rigid-water constraints.
    raw_bonds = raw.getBondsWithH() + raw.getBondsNoH()
    actual_bonds = {}
    for force in forces:
        if isinstance(force, HarmonicBondForce):
            for index in range(force.getNumBonds()):
                a, b, length, k = force.getBondParameters(index)
                actual_bonds[tuple(sorted((a, b)))] = (
                    length.value_in_unit(unit.nanometer),
                    k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2))
    constraints = {}
    for index in range(system.getNumConstraints()):
        a, b, length = system.getConstraintParameters(index)
        constraints[tuple(sorted((a, b)))] = length.value_in_unit(unit.nanometer)
    for a, b, k, distance in raw_bonds:
        pair = tuple(sorted((a, b)))
        if pair in constraints:
            if abs(constraints[pair] - distance) > 1e-8:
                raise WorkError("parameterMismatch", "Amber constrained bond length changed")
        elif pair not in actual_bonds or any(abs(x - y) > 1e-6 * max(1, abs(x))
                                              for x, y in zip((distance, 2 * k), actual_bonds[pair])):
            raise WorkError("parameterMismatch", "Amber bonded term changed on import")
    if not _system_bonds_match(topology, system):
        raise WorkError("parameterMismatch", "Imported Amber molecular bond graph changed")
    angle = next(force for force in forces if isinstance(force, HarmonicAngleForce))
    if angle.getNumAngles() != len(raw.getAngles()):
        raise WorkError("parameterMismatch", "Amber angle terms changed on import")
    for index, (a, b, c, k, theta) in enumerate(raw.getAngles()):
        aa, bb, cc, observed_theta, observed_k = angle.getAngleParameters(index)
        if ((aa, bb, cc) != (a, b, c) or
                abs(observed_theta.value_in_unit(unit.radian) - theta) > 1e-8 or
                abs(observed_k.value_in_unit(unit.kilojoule_per_mole / unit.radian ** 2) - 2 * k) > 1e-6):
            raise WorkError("parameterMismatch", "Amber angle value changed on import")
    torsion = next(force for force in forces if isinstance(force, PeriodicTorsionForce))
    if torsion.getNumTorsions() != len(raw.getDihedrals()):
        raise WorkError("parameterMismatch", "Amber torsion terms changed on import")
    for index, (a, b, c, d, k, phase, periodicity) in enumerate(raw.getDihedrals()):
        aa, bb, cc, dd, pp, observed_phase, observed_k = torsion.getTorsionParameters(index)
        if ((aa, bb, cc, dd, pp) != (a, b, c, d, periodicity) or
                abs(observed_phase.value_in_unit(unit.radian) - phase) > 1e-8 or
                abs(observed_k.value_in_unit(unit.kilojoule_per_mole) - k) > 1e-6):
            raise WorkError("parameterMismatch", "Amber torsion value changed on import")
    maps = [force for force in forces if isinstance(force, CMAPTorsionForce)]
    cmap_count = len(raw.getCMAPDihedrals()) if raw.getNumMaps() else 0
    if (len(maps) != int(bool(raw.getNumMaps())) or
            (maps[0].getNumTorsions() if maps else 0) != cmap_count):
        raise WorkError("parameterMismatch", "Amber CMAP torsions changed on import")
    if maps:
        for index, expected in enumerate(raw.getCMAPDihedrals()):
            if tuple(maps[0].getTorsionParameters(index)) != tuple(expected):
                raise WorkError("parameterMismatch", "Amber CMAP atom ordering changed")
        for index, resolution in enumerate(raw.getCMAPResolutions()):
            size, values = maps[0].getMapParameters(index)
            original = raw.getCMAPParameters(index + 1)
            expected = [original[size * ((j + size // 2) % size) + ((i + size // 2) % size)] * 4.184
                        for i in range(size) for j in range(size)]
            observed = [value.value_in_unit(unit.kilojoule_per_mole) for value in values]
            if size != int(resolution) or any(abs(a - b) > 1e-6 for a, b in zip(observed, expected)):
                raise WorkError("parameterMismatch", "Amber CMAP map values changed")
    expected_exceptions = {}
    for a, b, charge, rmin, epsilon, scee, scnb in raw.get14Interactions():
        expected_exceptions[tuple(sorted((a, b)))] = (charge / scee, rmin * 2 ** (-1 / 6), epsilon / scnb)
    for a, excluded in enumerate(raw.getExcludedAtoms()):
        for b in excluded:
            expected_exceptions.setdefault(tuple(sorted((a, b))), (0.0, 0.1, 0.0))
    actual_exceptions = {}
    for index in range(nonbonded.getNumExceptions()):
        a, b, charge, sigma, epsilon = nonbonded.getExceptionParameters(index)
        actual_exceptions[tuple(sorted((a, b)))] = (
            charge.value_in_unit(unit.elementary_charge ** 2),
            sigma.value_in_unit(unit.nanometer),
            epsilon.value_in_unit(unit.kilojoule_per_mole))
    if expected_exceptions.keys() != actual_exceptions.keys() or any(
            any(abs(a - b) > 1e-6 * max(1, abs(a)) for a, b in zip(values, actual_exceptions[key]))
            for key, values in expected_exceptions.items()):
        raise WorkError("parameterMismatch", "Amber exclusions or scaled 1-4 terms changed")
    raw_box = coordinate.boxVectors
    imported_box = system.getDefaultPeriodicBoxVectors()
    deviation = max(abs(raw_box[i][j].value_in_unit(unit.angstrom) -
                        imported_box[i][j].value_in_unit(unit.angstrom))
                    for i in range(3) for j in range(3))
    if deviation > 1e-5:
        raise WorkError("parameterMismatch", "Amber final restart cell changed on import")
    parameter_digest = hashlib.sha256(json.dumps({"chargesAndLj": converted,
        "bonds": sorted((list(pair), list(values)) for pair, values in actual_bonds.items()),
        "exclusions": sorted((list(pair), list(values)) for pair, values in actual_exceptions.items()),
        "forces": kinds}, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    account = {"prmtopSha256": sha256(top), "finalRestartSha256": sha256(restart),
               "amberAtomCount": n, "importedAtomCount": system.getNumParticles(),
               "amberNetChargeElementary": amber_charge,
               "importedNetChargeElementary": imported_charge,
               "maximumCellVectorDeviationAngstrom": deviation,
               "amberCmapTermCount": cmap_count,
               "importedCmapTermCount": maps[0].getNumTorsions() if maps else 0,
               "importedForceKinds": kinds, "atomOrderPreserved": True,
               "bondedTermsPreserved": True, "nonbondedTermsPreserved": True,
               "exclusionsPreserved": True, "unitsPreserved": True,
               "parameterComparisonPerformed": True,
               "parameterCorrespondenceSha256": parameter_digest}
    return topology, coordinate.positions, system, account


def _retained_parameter_correspondence(prepared_top: Path, final_top: Path,
                                       mapping: dict[int, int]) -> None:
    """Compare every retained Amber parameter term across source and build."""
    from openmm.app.internal.amber_file_parser import PrmtopLoader

    before = PrmtopLoader(str(prepared_top))
    after = PrmtopLoader(str(final_top))
    mapped_set = set(mapping.values())
    if before.getNumAtoms() != len(mapping) or len(mapped_set) != len(mapping):
        raise WorkError("parameterMismatch", "Retained Amber parameter map is incomplete")
    left_lj, right_lj = before.getNonbondTerms(), after.getNonbondTerms()
    for source, final in mapping.items():
        values = ((before.getCharges()[source], after.getCharges()[final]),
                  (before.getMasses()[source], after.getMasses()[final]),
                  *zip(left_lj[source], right_lj[final]))
        if any(abs(a - b) > 1e-6 * max(1, abs(a)) for a, b in values):
            raise WorkError("parameterMismatch", "A retained source atom's Amber parameters changed")
    def indexed_terms(terms: list[tuple], count: int, mapped: bool) -> Counter[tuple]:
        found = Counter()
        for term in terms:
            indices = term[:count]
            if mapped:
                indices = tuple(mapping[index] for index in indices)
            elif any(index not in mapped_set for index in indices):
                continue
            parameters = tuple(round(float(value), 7) for value in term[count:])
            found[(tuple(indices), parameters)] += 1
        return found
    for getter, count in (("getBondsWithH", 2), ("getBondsNoH", 2),
                          ("getAngles", 3), ("getDihedrals", 4),
                          ("getImpropers", 4), ("getUreyBradleys", 2)):
        source_terms = indexed_terms(getattr(before, getter)(), count, True)
        final_terms = indexed_terms(getattr(after, getter)(), count, False)
        if source_terms != final_terms:
            raise WorkError("parameterMismatch", f"Retained {getter} terms changed in the built system")
    source_14 = indexed_terms(before.get14Interactions(), 2, True)
    final_14 = indexed_terms(after.get14Interactions(), 2, False)
    if source_14 != final_14:
        raise WorkError("parameterMismatch", "Retained Amber 1-4 interactions changed")
    source_excluded = {tuple(sorted((mapping[i], mapping[j])))
                       for i, row in enumerate(before.getExcludedAtoms()) for j in row}
    final_excluded = {tuple(sorted((i, j))) for i, row in enumerate(after.getExcludedAtoms())
                      for j in row if i in mapped_set and j in mapped_set}
    if source_excluded != final_excluded:
        raise WorkError("parameterMismatch", "Retained Amber exclusions changed")
    def cmap_terms(loader: Any, remap: bool) -> Counter[tuple]:
        found = Counter()
        if not loader.getNumMaps():
            return found
        for term in loader.getCMAPDihedrals():
            indices = tuple(mapping[index] for index in term[1:]) if remap else tuple(term[1:])
            if not remap and any(index not in mapped_set for index in indices):
                continue
            values = loader.getCMAPParameters(term[0] + 1)
            found[(indices, tuple(round(value, 7) for value in values))] += 1
        return found
    if cmap_terms(before, True) != cmap_terms(after, False):
        raise WorkError("parameterMismatch", "Retained ff19SB CMAP terms changed")


def construct_memgen(directory: Path, payload: dict[str, Any], progress: Callable) -> dict[str, Any]:
    """Execute and observe one complete, exact Memgen trial."""
    from openmm import Context, NonbondedForce, VerletIntegrator, XmlSerializer, unit
    from openmm.app import PDBxFile
    from .construction import _read_topology_data
    from .local_state_observations import observe_local_state

    require_text(payload.get("studyRevisionId"), "studyRevisionId")
    require_text(payload.get("attemptId"), "attemptId")
    require_text(payload.get("trialId"), "trialId")
    require_integer(payload.get("trialIndex"), "trialIndex")
    settings = _settings(payload)
    lower, upper = _fractions(payload)
    targets = {(side, species): fraction for side, fractions in (("lower", lower), ("upper", upper))
               for species, fraction in fractions}
    runtime_value = os.environ.get("PIM_AMBERTOOLS_HOME")
    if not runtime_value:
        raise WorkError("providerUnavailable", "PIM_AMBERTOOLS_HOME does not identify the selected AmberTools runtime")
    amberhome = Path(runtime_value).resolve()
    if not amberhome.is_dir():
        raise WorkError("providerUnavailable", "The selected AmberTools26 runtime is unavailable")
    assets = _asset_paths(payload, amberhome)
    force_fields = _force_field_assets(payload, amberhome)
    maximum_atoms = require_integer(payload.get("maximumAtomCount"), "maximumAtomCount", 1)
    maximum_cell = require_number(payload.get("maximumCellDimensionAngstrom"), "maximumCellDimensionAngstrom", 1)
    if maximum_atoms != 120000 or maximum_cell != 180:
        raise WorkError("unsupportedPolicy", "Memgen construction resource policy changed")
    trial = directory / f"memgen-trial-{payload['trialIndex']}"
    if trial.exists():
        raise WorkError("invalidRequest", "An exact Memgen trial directory already exists")
    trial.mkdir()
    proposed: list[dict[str, Any]] = []
    achieved: list[dict[str, Any]] = []
    removed: list[dict[str, Any]] = []
    proposed_cell: list[float] = []
    actual_cell: list[float] = []
    conditions: dict[str, Any] | None = None
    packed_base = "system"
    try:
        input_path, source, disulfides = _input_pdb(directory, trial, payload)
        approved = _read_topology_data(work_path(directory, payload.get("preparedBondGraphPath"),
                                                 "preparedBondGraphPath"))
        prepared = sorted(payload["preparedCorrespondence"]["atoms"],
                          key=lambda item: item["resultAtomIndex"])
        charges, preflight_map = _source_charge(trial, amberhome, source, approved, disulfides)
        _run_provider(trial, payload, settings, lower, upper, amberhome, disulfides,
                      charges, progress)
        for path in list(assets.values()) + force_fields:
            # A changed runtime between authorization and completed output
            # cannot establish the selected provider/parameter provenance.
            expected = next((item["sha256"] for item in payload["providerAssets"]
                             if Path(item["path"]).resolve() == path), None)
            if expected is None:
                expected = next(item["sha256"] for item in payload["forceFieldFiles"]
                                if Path(item["path"]).resolve() == path)
            verify_sha256(path, expected, f"completed provider asset {path.name}")
        options_path = trial / "packmol-memgen.json"
        options = require_mapping(json.loads(options_path.read_text(encoding="utf-8")), "Memgen options")
        supplied_lipids = ":".join(name for name, _ in lower) + "//" + ":".join(name for name, _ in upper)
        supplied_ratios = ":".join(format(value, ".17g") for _, value in lower) + "//" + \
                          ":".join(format(value, ".17g") for _, value in upper)
        if (options.get("pdb") != [input_path.name] or options.get("lipids") != [supplied_lipids] or
                options.get("ratio") != [supplied_ratios] or options.get("saltcon") != 0.15 or
                options.get("ffprot") != "ff19SB" or options.get("fflip") != "lipid21" or
                options.get("ffwat") != "tip3p" or options.get("engine") != "sander" or
                any(options.get(key) is not True for key in ("preoriented", "notprotonate",
                    "keep", "pbc", "tight_box", "salt", "parametrize", "minimize")) or
                options.get("nottrim") is not False or options.get("noxy_cen") is not False or
                options.get("vol") is not False or
                options.get("salt_override") is not False or options.get("nocounter") is not False or
                options.get("random") is not False):
            raise WorkError("providerMismatch", "Memgen recorded options differ from selected invocation")
        provider_protein = _pdb_records(trial / "PROT0.pdb")
        if not _same_records(source, provider_protein):
            raise WorkError("correspondenceFailed", "Memgen changed the fixed prepared construct before packing")
        recipe = _packmol_recipe(trial / "packmol.inp")
        if not recipe or recipe[0] != {"file": "PROT0.pdb", "count": 1, "side": None}:
            raise WorkError("providerMismatch", "Packmol did not keep one fixed complete construct first")
        pbc = re.findall(r"^pbc\s+([-+\d. ]+)\s*$", (trial / "packmol.inp").read_text(encoding="utf-8"),
                         flags=re.MULTILINE)
        if len(pbc) != 1 or len(pbc[0].split()) != 6:
            raise WorkError("unobservedOutput", "Provider recipe lacks one measured periodic cell")
        bounds = list(map(float, pbc[0].split()))
        proposed_cell = [bounds[i + 3] - bounds[i] for i in range(3)]
        if any(not math.isfinite(value) or value <= 0 for value in proposed_cell):
            raise WorkError("invalidGeometry", "Provider proposed a nonfinite periodic cell")
        archive = amberhome / "lib/python3.12/site-packages/packmol_memgen/data/pdbs.tar.gz"
        templates = {species: _template_records(archive, species)
                     for species in _SPECIES | {"WAT", "Na+", "Cl-"}}
        proposed_counts: Counter[tuple[str, str]] = Counter()
        for block in recipe:
            species = Path(block["file"]).stem
            if species in _SPECIES:
                proposed_counts[(block["side"], species)] += block["count"]
        proposed = [{"physicalSide": side, "speciesId": species,
                     "count": proposed_counts[(side, species)], "intendedFraction": fraction}
                    for (side, species), fraction in sorted(targets.items())]
        if any(item["count"] == 0 for item in proposed):
            raise WorkError("zeroRoundedSpecies", "Provider rounded a requested leaflet species to zero")
        raw = _pdb_records(trial / "system.pdb")
        packed = _packed_molecules(raw, provider_protein, templates, recipe)
        log = (trial / "packmol-memgen.log").read_text(encoding="utf-8", errors="replace")
        if any(message in log for message in ("Lipid piercing finder failed", "Lipid piercing removal failed")):
            raise WorkError("incompleteCleanup", "Memgen could not establish sterol piercing cleanup")
        post_path = trial / "system_noclash.pdb"
        post = _pdb_records(post_path) if post_path.is_file() else raw
        packed_base = "system_noclash" if post_path.is_file() else "system"
        raw_to_post = _coordinate_correspondence(raw, post)
        if any(index not in raw_to_post for index in range(len(source))):
            raise WorkError("correspondenceFailed", "Provider cleanup removed a retained construct atom")
        retained_molecules = []
        removed_molecules = []
        for molecule in packed:
            indices = range(molecule["start"], molecule["start"] + len(molecule["atoms"]))
            present = [index in raw_to_post for index in indices]
            if all(present):
                retained_molecules.append(molecule)
            elif not any(present) and molecule["species"] in _SPECIES:
                removed_molecules.append(molecule)
            else:
                raise WorkError("correspondenceFailed", "Cleanup removed part of an unsupported molecule")
        achieved_counts = Counter((item["side"], item["species"]) for item in retained_molecules
                                  if item["species"] in _SPECIES)
        removed_counts = Counter((item["side"], item["species"]) for item in removed_molecules)
        achieved = [{"physicalSide": side, "speciesId": species, "count": achieved_counts[(side, species)],
                     "intendedFraction": fraction} for (side, species), fraction in sorted(targets.items())]
        removed = [{"physicalSide": side, "speciesId": species, "count": removed_counts[(side, species)],
                    "intendedFraction": fraction} for (side, species), fraction in sorted(targets.items())
                   if removed_counts[(side, species)]]
        if any(item["count"] == 0 for item in achieved):
            raise WorkError("zeroRoundedSpecies", "Cleanup removed every molecule of a requested species")
        leap_script = (trial / "leap.in").read_text(encoding="utf-8")
        if (f"loadpdb {packed_base}.pdb" not in leap_script or
                "saveAmberParmNetcdf" not in leap_script or "set SYS box" not in leap_script or
                "setBox SYS vdw" in leap_script):
            raise WorkError("providerMismatch", "LEaP did not parameterize the exact post-cleanup periodic result")
        parameterized = _pdb_records(trial / f"{packed_base}_lipid.pdb")
        post_to_leap = _coordinate_correspondence(post, parameterized)
        raw_to_leap = {raw_index: post_to_leap[post_index] for raw_index, post_index in raw_to_post.items()
                       if post_index in post_to_leap}
        if any(index not in raw_to_leap for index in range(len(source))):
            raise WorkError("correspondenceFailed", "LEaP removed a retained source atom")
        topology, positions, system, imported = _amber_import(
            trial / f"{packed_base}_lipid.top", trial / f"{packed_base}_min.restrt", payload)
        amber_shift, _ = _verify_amber_pdb_order(
            trial / f"{packed_base}_lipid.top",
            trial / f"{packed_base}_lipid.crd", parameterized, topology)
        atoms = list(topology.atoms())
        if len(atoms) > maximum_atoms:
            raise WorkError("resourceRefused", "Actual Memgen result exceeds the authorized atom limit")
        vectors = system.getDefaultPeriodicBoxVectors()
        actual_cell = [vectors[i][i].value_in_unit(unit.angstrom) for i in range(3)]
        if (any(abs(vectors[i][j].value_in_unit(unit.angstrom)) > 1e-5
                for i in range(3) for j in range(3) if i != j) or
                any(value <= 20 or value > maximum_cell for value in actual_cell) or
                any(abs(a - b) > 0.02 for a, b in zip(actual_cell, proposed_cell))):
            raise WorkError("invalidGeometry", "Amber final cell differs from the bounded Packmol cell")
        def centered_z(index: int) -> float:
            coordinate = positions[index][2].value_in_unit(unit.angstrom) - amber_shift[2]
            return _periodic_offset(coordinate, 0, actual_cell[2])

        mapped_source = {i: raw_to_leap[i] for i in range(len(source))}
        approved_bonds = {tuple(sorted((mapped_source[a.index], mapped_source[b.index])))
                          for a, b in approved.bonds()}
        actual_source_bonds = {tuple(sorted((a.index, b.index))) for a, b in topology.bonds()
                               if a.index in mapped_source.values() and b.index in mapped_source.values()}
        if approved_bonds != actual_source_bonds or any(
                atoms[mapped_source[i]].element != atom.element for i, atom in enumerate(approved.atoms())):
            raise WorkError("chemistryMismatch", "Amber changed retained prepared atoms or bonds")
        _retained_parameter_correspondence(trial / "prepared-charge.top",
                                           trial / f"{packed_base}_lipid.top",
                                           {preflight_map[i]: mapped_source[i] for i in mapped_source})
        components, component_of = _components(topology)
        references = {}
        representations = {}
        for representation in payload.get("selectedSpeciesRepresentations", []):
            if not isinstance(representation, dict) or representation.get("speciesId") not in _SPECIES:
                continue
            path = work_path(directory, representation.get("coordinateTemplatePath"),
                             f"{representation['speciesId']} coordinate template")
            verify_sha256(path, representation.get("coordinateTemplateSha256"),
                          f"{representation['speciesId']} coordinate template")
            references[representation["speciesId"]] = PDBxFile(str(path)).topology
            representations[representation["speciesId"]] = representation
        assignment: dict[int, tuple[str, str | None, str | None, str | None]] = {}
        source_roles = [item["moleculeRole"] for item in prepared]
        for index, final in mapped_source.items():
            assignment[final] = ("retained", source_roles[index], None, None)
        generated_counts: Counter[tuple[str, str, str]] = Counter()
        provider_generated = Counter(item["species"] for item in packed)
        removed_at_leap = Counter()
        for molecule in retained_molecules:
            indices = range(molecule["start"], molecule["start"] + len(molecule["atoms"]))
            finals = [raw_to_leap.get(i) for i in indices]
            species = molecule["species"]
            if any(item is None for item in finals):
                if all(item is None for item in finals) and species in {"WAT", "Na+", "Cl-"}:
                    removed_at_leap[species] += 1
                    continue
                raise WorkError("correspondenceFailed", f"LEaP lost part of {species}")
            if species in _SPECIES:
                component = components[component_of[finals[0]]]
                if any(component_of[index] != component_of[finals[0]] for index in finals):
                    raise WorkError("chemistryMismatch", f"LEaP split {species} molecular chemistry")
                mapping = _verified_molecular_graph(topology, component, references[species], species)
                _verify_stereo(references[species], mapping, positions,
                               representations[species].get("stereoChecks"), species)
                if component != set(finals):
                    raise WorkError("chemistryMismatch", f"LEaP added untracked atoms to {species}")
                head_name = _LIPID_HEAD[species]
                head_local = [i for i, record in enumerate(molecule["atoms"])
                              if record["name"] == head_name]
                if len(head_local) != 1:
                    raise WorkError("chemistryMismatch", f"{species} lacks its selected leaflet head")
                head_z = centered_z(finals[head_local[0]])
                if head_z == 0 or ("upper" if head_z > 0 else "lower") != molecule["side"]:
                    raise WorkError("chemistryMismatch", f"Conditioned {species} crossed its selected leaflet")
                role, product_species = "lipid", species
            elif species == "WAT":
                role, product_species = "water", require_text(payload["water"]["speciesId"], "water species")
            elif species == "Na+":
                role, product_species = "positiveIon", require_text(payload["sodium"]["speciesId"], "sodium species")
            elif species == "Cl-":
                role, product_species = "negativeIon", require_text(payload["chloride"]["speciesId"], "chloride species")
            else:
                raise WorkError("correspondenceFailed", "Provider returned an unselected generated molecule")
            if any(index in assignment for index in finals):
                raise WorkError("correspondenceFailed", "Provider mapped an atom to two molecules")
            final_z = centered_z(finals[0])
            final_side = molecule["side"] if role == "lipid" else ("upper" if final_z > 0 else "lower")
            for index in finals:
                assignment[index] = ("generated", role, product_species, final_side)
            generated_counts[(role, final_side, product_species)] += 1
        # LEaP's addionsrand can add residual neutralizers by substituting
        # waters.  Every unmatched atom must be one complete Na/Cl ion.
        leap_added = Counter()
        for component in components:
            if component & assignment.keys():
                if component - assignment.keys() and not component.issubset(mapped_source.values()):
                    raise WorkError("correspondenceFailed", "LEaP returned an unaccounted atom in a mapped molecule")
                continue
            if len(component) != 1:
                raise WorkError("correspondenceFailed", "LEaP added an unaccounted multiatom molecule")
            index = next(iter(component))
            symbol = atoms[index].element.symbol
            if symbol not in {"Na", "Cl"}:
                raise WorkError("correspondenceFailed", "LEaP added an unaccounted chemical species")
            role, product_species = (("positiveIon", payload["sodium"]["speciesId"])
                                     if symbol == "Na" else ("negativeIon", payload["chloride"]["speciesId"]))
            side = "upper" if centered_z(index) > 0 else "lower"
            assignment[index] = ("generated", role, product_species, side)
            generated_counts[(role, side, product_species)] += 1
            leap_added[symbol] += 1
        if len(assignment) != len(atoms):
            raise WorkError("correspondenceFailed", "Not every Amber atom has selected provenance")
        all_components = Counter()
        for component in components:
            if len(component) == 3 and Counter(atoms[i].element.symbol for i in component) == Counter({"O": 1, "H": 2}):
                all_components["WAT"] += 1
            elif len(component) == 1 and atoms[next(iter(component))].element.symbol in {"Na", "Cl"}:
                all_components[atoms[next(iter(component))].element.symbol] += 1
        retained_components = Counter()
        for component in _components(approved)[0]:
            elements = Counter(list(approved.atoms())[i].element.symbol for i in component)
            if elements == Counter({"O": 1, "H": 2}):
                retained_components["WAT"] += 1
            elif len(component) == 1 and next(iter(elements)) in {"Na", "Cl"}:
                retained_components[next(iter(elements))] += 1
        aqueous = []
        for side, label in (("lower", "Lower"), ("upper", "Upper")):
            volume = _log_number(log, f"{label} water box vol")
            floored = int(_log_number(log, f"{label} salt number"))
            if volume <= 0 or not _logged_salt_count_possible(volume, floored):
                raise WorkError("unobservedCondition", "Provider salt count differs from recorded aqueous volume")
            aqueous.append({"side": side, "estimatedVolumeAngstromCubed": volume,
                            "flooredNominalSaltCount": floored,
                            "providerGeneratedSodiumCount": sum(block["count"] for block in recipe
                                if block["file"] == "Na+.pdb" and block["side"] == side),
                            "providerGeneratedChlorideCount": sum(block["count"] for block in recipe
                                if block["file"] == "Cl-.pdb" and block["side"] == side)})
        volume = sum(item["estimatedVolumeAngstromCubed"] for item in aqueous)
        water_final = all_components["WAT"]
        sodium_final = all_components["Na"]
        chloride_final = all_components["Cl"]
        added_sodium = leap_added["Na"]
        added_chloride = leap_added["Cl"]
        if (water_final <= 0 or added_sodium < 0 or added_chloride < 0 or
                water_final != retained_components["WAT"] + provider_generated["WAT"] - removed_at_leap["WAT"] or
                sodium_final != retained_components["Na"] + provider_generated["Na+"] +
                    added_sodium - removed_at_leap["Na+"] or
                chloride_final != retained_components["Cl"] + provider_generated["Cl-"] +
                    added_chloride - removed_at_leap["Cl-"] or
                abs(imported["importedNetChargeElementary"]) > 1e-4):
            raise WorkError("unobservedCondition", "Amber water/ion or charge conservation is unaccounted")
        # Reporting uses the defined SI constant; the provider's older
        # 6.02214086 value above is retained only for its integer floor.
        molar_factor = 6.02214076e23 / 1e27
        conditions = {"retainedWaterCount": retained_components["WAT"],
                      "providerGeneratedWaterCount": provider_generated["WAT"], "finalWaterCount": water_final,
                      "retainedSodiumCount": retained_components["Na"],
                      "retainedChlorideCount": retained_components["Cl"],
                      "providerGeneratedSodiumCount": provider_generated["Na+"],
                      "providerGeneratedChlorideCount": provider_generated["Cl-"],
                      "leapAddedSodiumCount": added_sodium, "leapAddedChlorideCount": added_chloride,
                      "leapRemovedRetainedWaterCount": 0,
                      "leapRemovedGeneratedWaterCount": removed_at_leap["WAT"],
                      "leapRemovedRetainedSodiumCount": 0,
                      "leapRemovedGeneratedSodiumCount": removed_at_leap["Na+"],
                      "leapRemovedRetainedChlorideCount": 0,
                      "leapRemovedGeneratedChlorideCount": removed_at_leap["Cl-"],
                      "finalSodiumCount": sodium_final, "finalChlorideCount": chloride_final,
                      "preparedFormalChargeElementary": charges[0],
                      "providerResidueNameChargeElementary": charges[1],
                      "chargePdbDeltaElementary": charges[2],
                      "finalNetChargeElementary": imported["importedNetChargeElementary"],
                      "estimatedAqueousVolumeAngstromCubed": volume,
                      "sodiumAqueousMolar": sodium_final / (molar_factor * volume),
                      "chlorideAqueousMolar": chloride_final / (molar_factor * volume),
                      "sodiumFiniteWaterMolar": 55.4 * sodium_final / water_final,
                      "chlorideFiniteWaterMolar": 55.4 * chloride_final / water_final,
                      "saltBranch": _provider_salt_branch(charges[1], charges[2], aqueous),
                      "aqueousRegions": aqueous}
        if sum(item["providerGeneratedSodiumCount"] for item in aqueous) != provider_generated["Na+"] or \
                sum(item["providerGeneratedChlorideCount"] for item in aqueous) != provider_generated["Cl-"]:
            raise WorkError("unobservedCondition", "Provider ion placement differs from its recipe")
        points = [tuple(float(value) for value in point.value_in_unit(unit.angstrom)) for point in positions]
        gaps = _conditioned_retained_gaps(points, mapped_source, source,
                                         actual_cell, amber_shift)
        deviation = max(math.sqrt(sum(_periodic_offset(
            points[mapped_source[i]][axis] - amber_shift[axis], item["xyz"][axis], actual_cell[axis]) ** 2
            for axis in range(3))) for i, item in enumerate(source))
        for residue in topology.residues():
            if residue.insertionCode == " ":
                residue.insertionCode = ""
        topology_path = directory / f"constructed-{payload['trialId']}-topology.cif"
        with topology_path.open("w", encoding="utf-8") as stream:
            PDBxFile.writeFile(topology, positions, stream, keepIds=True)
        topology_json = directory / f"constructed-{payload['trialId']}-topology.json"
        topology_json.write_text(json.dumps(_topology_data(topology), separators=(",", ":")), encoding="utf-8")
        system_path = directory / f"constructed-{payload['trialId']}-system.xml"
        state_path = directory / f"constructed-{payload['trialId']}-state.xml"
        integrator = VerletIntegrator(0.001 * unit.picoseconds)
        context = Context(system, integrator)
        context.setPositions(positions)
        state = context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            raise WorkError("nonfiniteObservation", "Imported Amber final state has nonfinite energy")
        system_path.write_text(XmlSerializer.serialize(system), encoding="utf-8")
        state_path.write_text(XmlSerializer.serialize(state), encoding="utf-8")
        del context, integrator
        correspondence_atoms = []
        reverse_source = {final: index for index, final in mapped_source.items()}
        for atom in atoms:
            index = atom.index
            kind, role, species, assigned_side = assignment[index]
            if kind == "retained":
                source_atom = prepared[reverse_source[index]]
                provenance = {key: source_atom.get(key) for key in
                              ("sourceAtomId", "role", "sourceResidue", "approvedChangeId")}
                atom_role = source_atom["atomRole"]
                side_value = source_atom.get("physicalSide")
                generated_role = None
            else:
                provenance = {"sourceAtomId": None, "role": "generated",
                              "sourceResidue": None, "approvedChangeId": None}
                atom_role = "head" if role == "lipid" and atom.name in {"P31", "O1"} else "body"
                side_value = assigned_side
                generated_role = role
            correspondence_atoms.append({"resultAtomIndex": index,
                "resultAtomId": f"amber:{payload['trialId']}:{index}:{atom.name}",
                "moleculeRole": "ion" if role in {"positiveIon", "negativeIon"} else role,
                "atomRole": atom_role, "physicalSide": side_value,
                "element": atom.element.symbol,
                "generatedSpeciesId": species,
                "generatedComponentRole": generated_role, **provenance})
        correspondence_path = directory / f"constructed-{payload['trialId']}-correspondence.json"
        correspondence = {"sourceId": sha256(work_path(directory, payload.get("orientedPdbPath"), "orientedPdbPath")),
                          "resultId": sha256(topology_path), "atoms": correspondence_atoms,
                          "complete": True}
        correspondence_path.write_text(json.dumps(correspondence, separators=(",", ":")), encoding="utf-8")
        output_local_spec, output_midplane = _amber_local_observation_spec(payload, amber_shift[2])
        local_state = observe_local_state(topology, state.getPositions(), correspondence,
                                          output_local_spec)
        if local_state.get("standing") != "Observed":
            raise WorkError("unobservedLocalState", "Imported Amber state lacks required local observations")
        provider_artifacts = _provider_artifacts(directory, trial, packed_base, include_final=True)
        trial_account = _trial_summary(payload, "checked", proposed, achieved, removed,
                                       proposed_cell, actual_cell, conditions)
        species_counts = [{"role": role, "physicalSide": side, "speciesId": species, "count": count}
                          for (role, side, species), count in sorted(generated_counts.items()) if count]
        handoff = [artifact(directory, topology_path, "topologyCif"),
                   artifact(directory, topology_json, "topologyJson"),
                   artifact(directory, system_path, "systemXml"),
                   artifact(directory, state_path, "stateXml"),
                   artifact(directory, correspondence_path, "correspondenceJson")]
        return {"artifacts": handoff + provider_artifacts,
                "observations": {"atomCount": len(atoms), "speciesCounts": species_counts,
                    "actualCellAngstrom": actual_cell, "proteinPeriodicImageGapsAngstrom": gaps,
                    "waterCount": water_final, "positiveIonCount": sodium_final,
                    "negativeIonCount": chloride_final,
                    "netChargeElementary": imported["importedNetChargeElementary"],
                    "proteinNetChargeElementary": charges[0],
                    "initialPotentialEnergyKjMol": energy,
                    "correspondedResultAtomCount": len(correspondence_atoms),
                    "maximumProteinCoordinateDeviationAngstrom": deviation,
                    "outputFrameMidplaneZAngstrom": output_midplane,
                    "proteinIdentityAndBondsPreserved": True,
                    "nativePatchSha256": None, "nativePatchMode": None,
                    "nativeSourcePatchSha256": None,
                    "contactWarnings": [], "geometryWarnings": [], "parameterWarnings": [],
                    "localState": local_state, "trial": trial_account, "conditions": conditions,
                    "providerIntermediates": provider_artifacts, "amberImport": imported},
                "provider": {"name": "PACKMOL-Memgen", "version": "2026.3.25"}}
    except WorkError as error:
        recipe_path = trial / "packmol.inp"
        if recipe_path.is_file() and (not proposed or not proposed_cell):
            try:
                recipe = _packmol_recipe(recipe_path)
                if not proposed:
                    counts = Counter()
                    for block in recipe:
                        species = Path(block["file"]).stem
                        if species in _SPECIES:
                            counts[(block["side"], species)] += block["count"]
                    proposed = [{"physicalSide": side, "speciesId": species,
                                 "count": counts[(side, species)], "intendedFraction": fraction}
                                for (side, species), fraction in sorted(targets.items())]
                if not proposed_cell:
                    pbc = re.findall(r"^pbc\s+([-+\d. ]+)\s*$", recipe_path.read_text(),
                                     flags=re.MULTILINE)
                    if len(pbc) == 1 and len(pbc[0].split()) == 6:
                        bounds = list(map(float, pbc[0].split()))
                        proposed_cell = [bounds[i + 3] - bounds[i] for i in range(3)]
            except (WorkError, ValueError):
                pass
        partial = _provider_artifacts(directory, trial, packed_base, include_final=False)
        details = dict(error.details or {})
        details["trial"] = _trial_summary(payload, "failed", proposed, achieved, removed,
                                           proposed_cell, actual_cell, conditions,
                                           error.code, error.message)
        details["artifacts"] = partial
        raise WorkError(error.code, error.message, details) from error
    except Exception as error:
        # Import/parser defects and unexpected provider-file formats are failed
        # trials with attributable artifacts, never checked constructions.
        partial = _provider_artifacts(directory, trial, packed_base, include_final=False)
        message = f"The selected provider crossing failed: {type(error).__name__}: {error}"
        raise WorkError("providerMechanicsFailed", message,
                        {"trial": _trial_summary(payload, "failed", proposed, achieved, removed,
                                                  proposed_cell, actual_cell, conditions,
                                                  "providerMechanicsFailed", message),
                         "artifacts": partial}) from error
