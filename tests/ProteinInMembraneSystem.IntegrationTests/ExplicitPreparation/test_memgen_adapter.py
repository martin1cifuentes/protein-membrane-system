"""Focused complete-provider crossing checks; full Memgen routes run separately.

These controls exercise the real adapter and pinned local assets without
claiming that a replaced process establishes packing or final construction.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/ProteinInMembrane.Host"))

from ProteinInMembraneSystem.ExplicitPreparation.worker import packmol_memgen as adapter
from ProteinInMembraneSystem.worker.exchange import WorkError


AMBERHOME = ROOT / "out/ambertools26"
POLICY = ROOT / "config/policies/protein-membrane-current.json"
PINNED_MAIN = AMBERHOME / "lib/python3.12/site-packages/packmol_memgen/main.py"
PINNED_MAIN_SHA256 = "7892dd3ad13028ffc572ed84085c8cfef8a323f1dab9b032bdad7fe8a55ceeea"


def pinned_provider_statement(line: int) -> object:
    """Select one statement in the exact reviewed provider source for controlled cases."""
    if hashlib.sha256(PINNED_MAIN.read_bytes()).hexdigest() != PINNED_MAIN_SHA256:
        raise AssertionError("The reviewed Memgen source changed")
    matches = [node for node in ast.walk(ast.parse(PINNED_MAIN.read_text(encoding="utf-8")))
               if isinstance(node, (ast.If, ast.Assign)) and node.lineno == line]
    if len(matches) != 1:
        raise AssertionError(f"Expected one pinned provider branch at line {line}")
    return matches[0]


def run_pinned_salt_case(charge: int, lower_volume: float, upper_volume: float) -> dict:
    """Execute the provider's neutralization and salt blocks without packing."""
    variables = {"charge": charge, "solvent_vol_down": lower_volume,
                 "solvent_vol_up": upper_volume,
                 "solvent_vol_tot": lower_volume + upper_volume,
                 "avogadro": 6.02214086e23, "saltcon": 0.15,
                 "self": SimpleNamespace(salt=True, salt_c="Na+",
                                         ion_dict={"Na+": (None, 1)}, charge_imbalance=0),
                 "logger": mock.Mock(), "override_salt": False, "sys": sys,
                 "pos_up": 0, "pos_down": 0, "neg_up": 0, "neg_down": 0}
    statements = [pinned_provider_statement(line) for line in (1893, 1899, 1900, 1901)]
    exec(compile(ast.Module(body=statements, type_ignores=[]), str(PINNED_MAIN), "exec"), variables)
    return variables


def prepared_topology(histidine: str = "HIS", *, delta: bool = True,
                      epsilon: bool = False, terminal: bool = True):
    from openmm.app import Topology, element

    topology = Topology()
    chain = topology.addChain("A")
    methionine = topology.addResidue("MET", chain, "1")
    nitrogen = topology.addAtom("N", element.nitrogen, methionine)
    if terminal:
        for name in ("H", "H2", "H3"):
            topology.addBond(nitrogen, topology.addAtom(name, element.hydrogen, methionine))
    histidine_residue = topology.addResidue(histidine, chain, "2")
    nd1 = topology.addAtom("ND1", element.nitrogen, histidine_residue)
    ne2 = topology.addAtom("NE2", element.nitrogen, histidine_residue)
    if delta:
        topology.addBond(nd1, topology.addAtom("HD1", element.hydrogen, histidine_residue))
    if epsilon:
        topology.addBond(ne2, topology.addAtom("HE2", element.hydrogen, histidine_residue))
    return topology


class MemgenAdapterControls(unittest.TestCase):
    def test_leap_water_residue_alias_preserves_retained_atom_correspondence(self):
        before = [
            {"name": name, "resname": "HOH", "xyz": xyz}
            for name, xyz in (("O", (1.0, 2.0, 3.0)),
                              ("H1", (1.957, 2.0, 3.0)),
                              ("H2", (0.76, 2.927, 3.0)))
        ]
        after = [dict(atom, resname="WAT") for atom in before]
        self.assertEqual({0: 0, 1: 1, 2: 2},
                         adapter._coordinate_correspondence(before, after))

    def test_leap_replacement_ion_at_water_coordinate_keeps_identity_distinct(self):
        before = [
            {"name": "O", "resname": "WAT", "xyz": (1.0, 2.0, 3.0)},
            {"name": "H1", "resname": "WAT", "xyz": (1.5, 2.0, 3.0)},
            {"name": "H2", "resname": "WAT", "xyz": (2.0, 2.0, 3.0)},
        ]
        after = [{"name": "Cl-", "resname": "Cl-", "xyz": (2.0, 2.0, 3.0)}]
        self.assertEqual({}, adapter._coordinate_correspondence(before, after))

    def test_exact_prepared_hydrogens_choose_amber_names_without_changing_graph(self):
        for delta, epsilon, expected in ((True, False, "HID"), (False, True, "HIE"),
                                         (True, True, "HIP")):
            with self.subTest(state=expected):
                topology = prepared_topology(delta=delta, epsilon=epsilon)
                bonds_before = {(a.index, b.index) for a, b in topology.bonds()}
                translated = adapter._amber_prepared_names(topology)
                methionine_h = next(atom for atom in topology.atoms()
                                    if atom.residue.name == "MET" and atom.name == "H")
                self.assertEqual(("H1", "MET"), translated[methionine_h.index])
                self.assertEqual({expected}, {translated[atom.index][1] for atom in topology.atoms()
                                              if atom.residue.name == "HIS"})
                self.assertEqual(bonds_before, {(a.index, b.index) for a, b in topology.bonds()})

    def test_inconsistent_or_unselected_histidine_state_is_refused(self):
        for topology in (prepared_topology(delta=False, epsilon=False),
                         prepared_topology(histidine="HID", delta=False, epsilon=True)):
            with self.subTest(residue=next(list(topology.residues())[1].atoms()).name):
                with self.assertRaises(WorkError) as refused:
                    adapter._amber_prepared_names(topology)
                self.assertEqual("unrepresentableInput", refused.exception.code)

    def test_approved_disulfides_use_cyx_without_duplicate_conect(self):
        from openmm import Vec3, unit
        from openmm.app import PDBFile, Topology, element
        from ProteinInMembraneSystem.worker.parameterized_structure import _topology_data

        topology = Topology()
        chain = topology.addChain("A")
        first = topology.addResidue("CYS", chain, "1")
        second = topology.addResidue("CYS", chain, "2")
        ligand = topology.addResidue("LIG", chain, "3")
        sulfur_a = topology.addAtom("SG", element.sulfur, first)
        sulfur_b = topology.addAtom("SG", element.sulfur, second)
        carbon_a = topology.addAtom("C1", element.carbon, ligand)
        carbon_b = topology.addAtom("C2", element.carbon, ligand)
        topology.addBond(sulfur_a, sulfur_b)
        topology.addBond(carbon_a, carbon_b)
        positions = [Vec3(index * 0.2, 0, 0) for index in range(4)] * unit.nanometer
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            prepared = directory / "prepared.pdb"
            with prepared.open("w", encoding="ascii") as stream:
                PDBFile.writeFile(topology, positions, stream, keepIds=True)
            graph = directory / "prepared-bonds.json"
            graph.write_text(json.dumps(_topology_data(PDBFile(str(prepared)).topology)),
                             encoding="utf-8")
            digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            payload = {"orientedPdbPath": str(prepared), "orientedPdbSha256": digest(prepared),
                       "preparedPdbPath": str(prepared), "preparedPdbSha256": digest(prepared),
                       "preparedBondGraphPath": str(graph), "preparedBondGraphSha256": digest(graph),
                       "preparedCorrespondence": {"complete": True,
                           "atoms": [{"resultAtomIndex": index} for index in range(4)]}}
            output, records, bonds = adapter._input_pdb(directory, directory, payload)
            self.assertEqual([(1, 2)], bonds)
            self.assertEqual(["CYX", "CYX", "LIG", "LIG"],
                             [record["resname"] for record in records])
            conect = [line for line in output.read_text().splitlines()
                      if line.startswith("CONECT")]
            self.assertEqual(["CONECT    3    4", "CONECT    4    3"], conect)

        protonated = Topology()
        protonated_chain = protonated.addChain("A")
        protonated_first = protonated.addResidue("CYS", protonated_chain, "1")
        protonated_sulfur = protonated.addAtom("SG", element.sulfur, protonated_first)
        thiol_hydrogen = protonated.addAtom("HG", element.hydrogen, protonated_first)
        protonated_second = protonated.addResidue("CYS", protonated_chain, "2")
        other_sulfur = protonated.addAtom("SG", element.sulfur, protonated_second)
        protonated.addBond(protonated_sulfur, thiol_hydrogen)
        protonated.addBond(protonated_sulfur, other_sulfur)
        with self.assertRaises(WorkError) as refused:
            adapter._amber_prepared_names(protonated)
        self.assertEqual("unrepresentableInput", refused.exception.code)

    def test_residue_charge_estimate_is_exactly_the_pinned_provider_table(self):
        source = AMBERHOME / "lib/python3.12/site-packages/packmol_memgen/lib/utils.py"
        if not source.is_file():
            self.skipTest("Pinned Memgen source is unavailable")
        self.assertEqual("c36d2784ec7957111059c887ea002ea68bd90bfb25377094fc4f0793b102f395",
                         hashlib.sha256(source.read_bytes()).hexdigest())
        tree = ast.parse(source.read_text(encoding="utf-8"))
        declared = next(ast.literal_eval(node.value) for node in tree.body
                        if isinstance(node, ast.Assign) and any(
                            isinstance(target, ast.Name) and target.id == "charged"
                            for target in node.targets))
        self.assertEqual(46, len(declared))
        self.assertEqual(declared, adapter._PROVIDER_CHARGED_RESIDUES)

    def test_charge_estimate_uses_exact_reindexed_provider_input(self):
        # Distinct charged residues 10A/10B share the original PDB number.
        # Exercise the actual input crossing before estimating provider charge.
        from openmm.app import PDBFile
        from ProteinInMembraneSystem.worker.parameterized_structure import _topology_data

        def atom(serial, residue, insertion, x):
            return (f"ATOM  {serial:5d} {'CA':4s} {residue:>3s} A"
                    f"{10:4d}{insertion}   {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
                    f"{1.0:6.2f}{0.0:6.2f}          {'C':>2s}\n")

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            prepared = directory / "prepared.pdb"
            oriented = directory / "oriented.pdb"
            prepared.write_text(atom(1, "ASP", "A", 0.0) +
                                atom(2, "GLU", "B", 1.0) + "END\n", encoding="ascii")
            oriented.write_bytes(prepared.read_bytes())
            graph = directory / "prepared-bonds.json"
            graph.write_text(json.dumps(_topology_data(PDBFile(str(prepared)).topology)),
                             encoding="utf-8")
            digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            payload = {"orientedPdbPath": str(oriented), "orientedPdbSha256": digest(oriented),
                       "preparedPdbPath": str(prepared), "preparedPdbSha256": digest(prepared),
                       "preparedBondGraphPath": str(graph), "preparedBondGraphSha256": digest(graph),
                       "preparedCorrespondence": {"complete": True,
                           "atoms": [{"resultAtomIndex": 0}, {"resultAtomIndex": 1}]}}
            trial = directory / "trial"
            trial.mkdir()
            path, source, disulfides = adapter._input_pdb(directory, trial, payload)
            self.assertEqual([], disulfides)
            self.assertEqual(["10", "10"], [item["resid"] for item in source])
            self.assertEqual(["1", "2"], [item["resid"] for item in adapter._pdb_records(path)])
            self.assertEqual(-2.0, adapter._provider_residue_name_charge(path, source))
            changed = [dict(item) for item in source]
            changed[1]["name"] = "CB"
            with self.assertRaises(WorkError) as refused:
                adapter._provider_residue_name_charge(path, changed)
            self.assertEqual("inputMismatch", refused.exception.code)

    def test_salt_branch_uses_provider_integer_charge_at_the_floored_boundary(self):
        regions = [{"flooredNominalSaltCount": 2}, {"flooredNominalSaltCount": 3}]
        # LEaP may measure 3.9999999 e for a formally +4 e construct.
        # That measured value would choose the wrong side of the strict `<`.
        measured_formal_charge = 3.9999999
        self.assertLess(abs(measured_formal_charge) / 2, 2)
        self.assertEqual("neutralizationOnly",
                         adapter._provider_salt_branch(3.0, 1.0000001, regions))
        self.assertEqual("chargeCompensated",
                         adapter._provider_salt_branch(2.0, 1.0, regions))

    def test_logged_volume_rounding_does_not_relabel_provider_salt_floor(self):
        # A true volume just above this boundary floors to two ions, while
        # Memgen's two-decimal log displays a volume that floors to one.
        boundary = 2 / (0.15 * 6.02214086e23 / 1e27)
        reported = round(boundary + 0.001, 2)
        self.assertEqual(1, int(reported * 0.15 * 6.02214086e23 / 1e27))
        self.assertTrue(adapter._logged_salt_count_possible(reported, 2))
        self.assertTrue(adapter._logged_salt_count_possible(reported, 1))
        self.assertFalse(adapter._logged_salt_count_possible(reported, 3))

    def test_provider_command_uses_lower_then_upper_and_selected_controls(self):
        class FinishedProcess:
            def poll(self):
                return 0

        with tempfile.TemporaryDirectory() as temporary:
            trial = Path(temporary)
            payload = {"trialId": "controlled-command", "trialIndex": 0,
                       "lateralPaddingAngstrom": 15.0,
                       "aqueousPaddingAngstrom": 17.5}
            observed_progress = []
            with mock.patch.object(adapter.subprocess, "Popen", return_value=FinishedProcess()) as invocation:
                adapter._run_provider(trial, payload, {}, [("DOPC", 0.7), ("CHL1", 0.3)],
                                      [("DPPC", 1.0)], AMBERHOME, [(2, 8)],
                                      (-1.0, -2.0, 1.0),
                                      lambda stage, detail: observed_progress.append(stage))
            argv = invocation.call_args.args[0]
            self.assertEqual("DOPC:CHL1//DPPC", argv[argv.index("--lipids") + 1])
            self.assertEqual("0.69999999999999996:0.29999999999999999//1",
                             argv[argv.index("--ratio") + 1])
            self.assertEqual("1", argv[argv.index("--charge_pdb_delta") + 1])
            self.assertEqual("bond SYS.2.SG SYS.8.SG", argv[argv.index("--leapline") + 1])
            for flag in ("--preoriented", "--notprotonate", "--nottrim", "--noxy_cen",
                         "--keep", "--pbc", "--tight_box", "--salt", "--parametrize",
                         "--minimize"):
                self.assertIn(flag, argv)
            for flag in ("--random", "--vol", "--salt_override", "--nocounter"):
                self.assertNotIn(flag, argv)
            self.assertEqual(["providerPopulation"], observed_progress)
            self.assertEqual(argv, json.loads((trial / "provider-command.json").read_text()))

    def test_pinned_unequal_water_boxes_can_suppress_salt_without_high_charge_refusal(self):
        # Controlled provider-source case, not a completed Memgen route.  The
        # lower box floors to two nominal pairs and the upper to four; charge
        # four reaches the strict half-charge boundary only in the lower box.
        values = run_pinned_salt_case(charge=4, lower_volume=26016.81,
                                      upper_volume=53000.0)
        self.assertEqual((2, 4), (values["saltnum_down"], values["saltnum_up"]))
        self.assertLess(values["con_neg"], 0.15)
        self.assertEqual((0, 0, 1, 3), (values["pos_down"], values["pos_up"],
                                          values["neg_down"], values["neg_up"]))
        self.assertEqual(4, values["neg_down"] + values["neg_up"])
        self.assertLess(values["neg_down"] + values["neg_up"],
                        values["saltnum_down"] + values["saltnum_up"])
        self.assertEqual(2, abs(values["charge"]) / 2)
        self.assertEqual(2, min(values["saltnum_down"], values["saltnum_up"]))

    def test_pinned_high_charge_refusal_becomes_typed_adapter_failure(self):
        # The real provider block refuses its fixed 0.15 M input before
        # packing; a failed process returns this exact source-owned message.
        with self.assertRaises(SystemExit) as refused:
            run_pinned_salt_case(charge=40, lower_volume=26016.81,
                                 upper_volume=53000.0)
        self.assertIn("concentration of ions required to neutralize",
                      str(refused.exception))

        class FailedProcess:
            def poll(self):
                return 1

        with tempfile.TemporaryDirectory() as temporary:
            trial = Path(temporary)
            (trial / "packmol-memgen.log").write_text(str(refused.exception), encoding="utf-8")
            payload = {"trialId": "controlled-high-charge", "trialIndex": 0,
                       "lateralPaddingAngstrom": 15.0,
                       "aqueousPaddingAngstrom": 17.5}
            with mock.patch.object(adapter.subprocess, "Popen", return_value=FailedProcess()):
                with self.assertRaises(WorkError) as translated:
                    adapter._run_provider(trial, payload, {}, [("DOPC", 1.0)],
                                          [("DOPC", 1.0)], AMBERHOME, [],
                                          (40.0, 40.0, 0.0), lambda *_: None)
            self.assertEqual("chargeBudgetRefusal", translated.exception.code)
            command = json.loads((trial / "provider-command.json").read_text())
            for flag, selected in (("--saltcon", "0.15"), ("--dist", "15"),
                                   ("--dist_wat", "17.5")):
                self.assertEqual(selected, command[command.index(flag) + 1])
            self.assertNotIn("--salt_override", command)

    def test_pinned_zero_rounded_positive_species_retains_intended_ratio_and_fails(self):
        # One controlled 21-lipid leaflet makes a positive 0.001 CHL1
        # fraction round to zero by Memgen's own integer-population branch.
        branch = pinned_provider_statement(1706)
        variables = {"composition": {0: {1: {"DOPC": 0.999, "CHL1": 0.001}}},
                     "bilayer": 0, "leaflet": 1, "lipid": "CHL1",
                     "lipnum_dict": {0: {1: 21}}, "logger": mock.Mock(), "sys": sys}
        with self.assertRaises(SystemExit) as refused:
            exec(compile(ast.Module(body=[branch], type_ignores=[]),
                         str(PINNED_MAIN), "exec"), variables)
        self.assertEqual(1, refused.exception.code)
        provider_message = variables["logger"].error.call_args.args[0]
        self.assertIn("ratio for lipid CHL1 is too small", provider_message)

        class FailedProcess:
            def poll(self):
                return 1

        with tempfile.TemporaryDirectory() as temporary:
            trial = Path(temporary)
            (trial / "packmol-memgen.log").write_text(provider_message, encoding="utf-8")
            payload = {"trialId": "controlled-trace-chl1", "trialIndex": 0,
                       "lateralPaddingAngstrom": 15.0,
                       "aqueousPaddingAngstrom": 17.5}
            with mock.patch.object(adapter.subprocess, "Popen", return_value=FailedProcess()):
                with self.assertRaises(WorkError) as translated:
                    adapter._run_provider(trial, payload, {},
                                          [("DOPC", 0.999), ("CHL1", 0.001)],
                                          [("DOPC", 1.0)], AMBERHOME, [],
                                          (0.0, 0.0, 0.0), lambda *_: None)
            self.assertEqual("zeroRoundedSpecies", translated.exception.code)
            command = json.loads((trial / "provider-command.json").read_text())
            self.assertEqual("DOPC:CHL1//DOPC", command[command.index("--lipids") + 1])
            self.assertEqual("0.999:0.001//1", command[command.index("--ratio") + 1])
            self.assertNotIn("--random", command)

        # The same controlled population keeps an ordinary positive fraction.
        variables["composition"][0][1] = {"DOPC": 0.95, "CHL1": 0.05}
        variables["logger"].reset_mock()
        exec(compile(ast.Module(body=[branch], type_ignores=[]),
                     str(PINNED_MAIN), "exec"), variables)
        variables["logger"].error.assert_not_called()

    def test_prepacking_refusal_keeps_one_failed_trial_without_invented_counts(self):
        # Controlled transport fixture: Memgen has not written a Packmol
        # recipe, so no proposed or achieved population can be reported.
        from openmm.app import PDBFile
        from ProteinInMembraneSystem.worker.parameterized_structure import _topology_data

        for code in ("zeroRoundedSpecies", "chargeBudgetRefusal"):
            with self.subTest(failure=code), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                prepared = directory / "prepared.pdb"
                prepared.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000"
                                    "  1.00  0.00           C\nEND\n", encoding="ascii")
                oriented = directory / "oriented.pdb"
                oriented.write_bytes(prepared.read_bytes())
                graph = directory / "prepared-graph.json"
                graph.write_text(json.dumps(_topology_data(PDBFile(str(prepared)).topology)))
                digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
                payload = {"studyRevisionId": "controlled-revision", "attemptId": "controlled-attempt",
                           "trialId": "controlled-trial-0", "trialIndex": 0,
                           "lateralPaddingAngstrom": 15.0, "aqueousPaddingAngstrom": 17.5,
                           "maximumAtomCount": 120000, "maximumCellDimensionAngstrom": 180.0,
                           "preparedPdbPath": str(prepared), "preparedPdbSha256": digest(prepared),
                           "orientedPdbPath": str(oriented), "orientedPdbSha256": digest(oriented),
                           "preparedBondGraphPath": str(graph), "preparedBondGraphSha256": digest(graph),
                           "preparedCorrespondence": {"complete": True,
                                                      "atoms": [{"resultAtomIndex": 0}]},
                           "selectedSpeciesRepresentations": [
                               {"speciesId": "DOPC"}, {"speciesId": "CHL1"}],
                           "targetLower": {"physicalSide": "lower", "fractions": [
                               {"speciesId": "DOPC", "fraction": 0.999},
                               {"speciesId": "CHL1", "fraction": 0.001}]},
                           "targetUpper": {"physicalSide": "upper", "fractions": [
                               {"speciesId": "DOPC", "fraction": 1.0}]}}
                with mock.patch.dict(adapter.os.environ, {"PIM_AMBERTOOLS_HOME": str(AMBERHOME)}), \
                        mock.patch.object(adapter, "_settings", return_value={}), \
                        mock.patch.object(adapter, "_asset_paths", return_value={}), \
                        mock.patch.object(adapter, "_force_field_assets", return_value=[]), \
                        mock.patch.object(adapter, "_source_charge",
                                          return_value=((0.0, 0.0, 0.0), {})), \
                        mock.patch.object(adapter, "_run_provider",
                                          side_effect=WorkError(code, "Pinned provider refusal")):
                    with self.assertRaises(WorkError) as refused:
                        adapter.construct_memgen(directory, payload, lambda *_: None)
                self.assertEqual(code, refused.exception.code)
                trial = refused.exception.details["trial"]
                self.assertEqual(("controlled-trial-0", 0, "failed", code),
                                 (trial["trialId"], trial["trialIndex"],
                                  trial["standing"], trial["failureCode"]))
                self.assertEqual((15.0, 17.5),
                                 (trial["lateralPaddingAngstrom"],
                                  trial["aqueousPaddingAngstrom"]))
                self.assertEqual([], trial["proposedLipidCounts"])
                self.assertEqual([], trial["achievedLipidCounts"])
                self.assertIsNone(trial["conditions"])
                self.assertEqual({"providerInput"},
                                 {item["role"] for item in refused.exception.details["artifacts"]})

    def test_provider_phases_require_attributable_tool_output(self):
        class OneObservedPoll:
            def __init__(self):
                self.polls = 0

            def poll(self):
                self.polls += 1
                return None if self.polls == 1 else 0

        payload = {"trialId": "controlled-progress", "trialIndex": 0,
                   "lateralPaddingAngstrom": 15.0,
                   "aqueousPaddingAngstrom": 17.5}
        phases = ["providerPopulation", "providerPacking", "providerCleanup",
                  "amberParameterization", "providerConditioningRestrained",
                  "providerConditioningUnrestrained"]
        cases = (
            ({}, phases[:1]),
            ({"packmol.log": "PACKMOL input accepted\n"}, phases[:1]),
            ({"packmol.log": "Starting GENCAN loop: 0\n"}, phases[:2]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n"}, phases[:3]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": ""}, phases[:3]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": "LEaP started\n"}, phases[:4]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": "LEaP started\n",
              "min.in": "restrained input\n", "min.out": ""}, phases[:4]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": "LEaP started\n",
              "min.in": "restrained input\n", "min.out": "sander started\n"}, phases[:5]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": "LEaP started\n",
              "min.out": "sander started\n", "system_min.in": "unrestrained input\n",
              "system_min.out": ""}, phases[:5]),
            ({"packmol.log": "Starting GENCAN loop: 0\n",
              "leap.in": "postprocessed leap input\n", "leap_.log": "LEaP started\n",
              "min.out": "sander started\n", "system_min.in": "unrestrained input\n",
              "system_min.out": "unrestrained sander started\n"}, phases),
        )
        for files, expected in cases:
            with self.subTest(files=tuple(files)), tempfile.TemporaryDirectory() as temporary:
                trial = Path(temporary)
                (trial / "packmol.inp").write_text("packmol recipe\n")
                (trial / "system.pdb").write_text("in-progress packing output\n")
                for name, content in files.items():
                    (trial / name).write_text(content)
                stages = []
                with mock.patch.object(adapter.subprocess, "Popen", return_value=OneObservedPoll()), \
                        mock.patch.object(adapter.time, "sleep"):
                    adapter._run_provider(trial, payload, {}, [("DOPC", 1.0)],
                                          [("DOPC", 1.0)], AMBERHOME, [],
                                          (0.0, 0.0, 0.0),
                                          lambda stage, detail: stages.append(stage))
                self.assertEqual(expected, stages)

    def test_missing_or_unknown_named_construction_route_cannot_select_native(self):
        from ProteinInMembraneSystem.ExplicitPreparation.worker import construction

        for route in (None, "unknownLegacyRoute"):
            with self.subTest(route=route), tempfile.TemporaryDirectory() as temporary:
                payload = {"studyRevisionId": "exact-revision", "attemptId": "one-attempt"}
                if route is not None:
                    payload["route"] = route
                with self.assertRaises(WorkError) as refused:
                    construction.construct_system(Path(temporary), payload,
                                                  lambda _stage, _detail: None)
                self.assertEqual("unsupportedPolicy", refused.exception.code)

    def test_leap_raw_amber_order_accepts_only_proven_openmm_name_normalization(self):
        from openmm import Vec3, unit
        from openmm.app import Topology, element

        topology = Topology()
        residue = topology.addResidue("ALA", topology.addChain(), "1")
        topology.addAtom("N", element.nitrogen, residue)
        topology.addAtom("H", element.hydrogen, residue)  # OpenMM's standard H1 name

        class RawAmber:
            has_atomic_number = True
            _raw_data = {"ATOMIC_NUMBER": [7, 1]}

            def getNumAtoms(self):
                return 2

            def getAtomName(self, index):
                return ("N", "H1")[index]

            def getResidueLabel(self, index):
                return "ALA"

        class InitialCoordinates:
            positions = [Vec3(5, 5, 5), Vec3(6, 5, 5)] * unit.angstrom

        records = [{"name": "N", "resname": "ALA", "element": "N", "xyz": (0, 0, 0)},
                   {"name": "H1", "resname": "ALA", "element": "H", "xyz": (1, 0, 0)}]
        with mock.patch("openmm.app.internal.amber_file_parser.PrmtopLoader",
                        return_value=RawAmber()), \
                mock.patch("openmm.app.AmberInpcrdFile", return_value=InitialCoordinates()):
            shift, residual = adapter._verify_amber_pdb_order(
                Path("raw.top"), Path("initial.crd"), records, topology)
            self.assertEqual((5.0, 5.0, 5.0), shift)
            self.assertLess(residual, 1e-12)
            self.assertAlmostEqual(10.0, adapter._periodic_offset(50.0, 40.0, 80.0))
            self.assertAlmostEqual(-10.0, adapter._periodic_offset(30.0, 40.0, 80.0))
            self.assertAlmostEqual(10.0, adapter._periodic_offset(130.0, 40.0, 80.0))
            changed_name = [dict(item) for item in records]
            changed_name[1]["name"] = "H2"
            with self.assertRaises(WorkError) as refused:
                adapter._verify_amber_pdb_order(Path("raw.top"), Path("initial.crd"),
                                                changed_name, topology)
            self.assertEqual("correspondenceFailed", refused.exception.code)
            changed_position = [dict(item) for item in records]
            changed_position[1]["xyz"] = (1.01, 0, 0)
            with self.assertRaises(WorkError) as refused:
                adapter._verify_amber_pdb_order(Path("raw.top"), Path("initial.crd"),
                                                changed_position, topology)
            self.assertEqual("correspondenceFailed", refused.exception.code)

    def test_local_observation_uses_measured_amber_frame_without_mutating_policy(self):
        selected = {"referenceMidplaneZAngstrom": 0.0,
                    "requiredMetricNames": ["upperLipidHeadMeanZAngstrom"]}
        payload = {"membraneCenterZNanometers": 0.0, "localObservationSpec": selected}
        mapped, midplane = adapter._amber_local_observation_spec(payload, 42.5545)
        self.assertAlmostEqual(42.5545, midplane)
        self.assertAlmostEqual(midplane, mapped["referenceMidplaneZAngstrom"])
        self.assertEqual(0.0, selected["referenceMidplaneZAngstrom"])
        with self.assertRaises(WorkError) as refused:
            adapter._amber_local_observation_spec(
                {"membraneCenterZNanometers": 0.0,
                 "localObservationSpec": {"referenceMidplaneZAngstrom": 3.0}}, 42.5545)
        self.assertEqual("inputMismatch", refused.exception.code)

    def test_conditioned_gap_uses_mapped_final_atoms_and_unwraps_box_faces(self):
        source = [{"xyz": (-1.0, 0.0, 0.0)}, {"xyz": (1.0, 0.0, 0.0)}]
        # Final atoms are reordered, and LEaP translated their frame by 40 Å.
        points = [(42.0, 40.0, 40.0), (38.0, 40.0, 40.0)]
        gaps = adapter._conditioned_retained_gaps(points, {0: 1, 1: 0}, source,
                                                 [80.0, 80.0, 80.0], (40.0, 40.0, 40.0))
        self.assertEqual([76.0, 80.0, 80.0], gaps)
        # A conditioned atom that crosses the face stays in the same source image.
        across_face = [{"xyz": (36.0, 0.0, 0.0)}, {"xyz": (39.0, 0.0, 0.0)}]
        self.assertEqual(75.0, adapter._conditioned_retained_gaps(
            [(36.0, 0.0, 0.0), (-39.0, 0.0, 0.0)], {0: 0, 1: 1},
            across_face, [80.0, 80.0, 80.0], (0.0, 0.0, 0.0))[0])
        with self.assertRaises(WorkError) as refused:
            adapter._conditioned_retained_gaps(
                [(-39.5, 0.0, 0.0), (39.5, 0.0, 0.0)], {0: 0, 1: 1},
                [{"xyz": (39.0, 0.0, 0.0)}, {"xyz": (-39.0, 0.0, 0.0)}],
                [80.0, 80.0, 80.0], (0.0, 0.0, 0.0))
        self.assertEqual("insufficientCellClearance", refused.exception.code)
        self.assertEqual("lateral", refused.exception.details["axis"])

    def test_selected_lipid_graph_rejects_missing_bond(self):
        from openmm.app import PDBxFile, Topology

        catalogue = json.loads(POLICY.read_text(encoding="utf-8"))
        selected = next(item for item in catalogue["lipids"] if item["speciesId"] == "DOPC")
        reference = PDBxFile(str(POLICY.parent / selected["coordinateTemplatePath"]))
        topology = reference.topology
        indices = set(range(topology.getNumAtoms()))
        self.assertEqual(len(indices), len(adapter._verified_molecular_graph(
            topology, indices, topology, "DOPC")))
        changed = Topology()
        chains = {chain.index: changed.addChain(chain.id) for chain in topology.chains()}
        atoms = {}
        for residue in topology.residues():
            new_residue = changed.addResidue(residue.name, chains[residue.chain.index], residue.id)
            for atom in residue.atoms():
                atoms[atom.index] = changed.addAtom(atom.name, atom.element, new_residue)
        for index, (first, second) in enumerate(topology.bonds()):
            if index:
                changed.addBond(atoms[first.index], atoms[second.index])
        with self.assertRaises(WorkError) as refused:
            adapter._verified_molecular_graph(changed, indices, topology, "DOPC")
        self.assertEqual("chemistryMismatch", refused.exception.code)

    def test_partial_provider_artifacts_are_hashed_but_not_a_handoff(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            trial = directory / "trial"
            trial.mkdir()
            (trial / "construct.pdb").write_bytes(b"prepared input\n")
            (trial / "packmol-memgen.log").write_bytes(b"observed provider log\n")
            artifacts = adapter._provider_artifacts(directory, trial, "system", include_final=False)
            by_role = {item["role"]: item for item in artifacts}
            self.assertEqual({"providerInput", "providerLog"}, set(by_role))
            for item in artifacts:
                self.assertEqual(hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest(),
                                 item["sha256"])
            with self.assertRaises(WorkError) as refused:
                adapter._provider_artifacts(directory, trial, "system", include_final=True)
            self.assertEqual("unobservedOutput", refused.exception.code)


if __name__ == "__main__":
    unittest.main()
