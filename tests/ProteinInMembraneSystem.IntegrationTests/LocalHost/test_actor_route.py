"""Exercise slice 1 through the published loopback host and real Python worker.

Run after scripts/build-local.sh using Python 3.11 in an environment that can
bind loopback: python3.11 tests/ProteinInMembraneSystem.IntegrationTests/LocalHost/test_actor_route.py.
Add --remote to exercise the actual selected RCSB and AlphaFold providers.
"""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "ProteinInMembraneSystem.IntegrationTests" / "ProteinPreparation"))
from test_worker_exchange import atom_line, two_alanines  # noqa: E402


def send(base: str, path: str, method: str = "GET", data: bytes | None = None,
         content_type: str | None = None, origin: str | None = None):
    headers = {}
    if content_type:
        headers["Content-Type"] = content_type
    if origin:
        headers["Origin"] = origin
    request = Request(base + path, data=data, headers=headers, method=method)
    try:
        with urlopen(request, timeout=60) as response:
            return response.status, response.read(), response.headers
    except HTTPError as failure:
        return failure.code, failure.read(), failure.headers


def command(base: str, kind: str, data: dict, revision: int, origin: str | None = None):
    return send(base, "/api/commands", "POST",
                json.dumps({"kind": kind, "data": data, "expectedRevision": revision}).encode(),
                "application/json", origin or base)


def require(condition: bool, message: str):
    if not condition:
        raise AssertionError(message)


def upload_pdb(base: str, source_text: str, provenance: str,
               filename: str = "protein.pdb") -> str:
    boundary = "slice1-exact-protein-upload"
    upload_body = (
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"provenance\"\r\n\r\n"
        f"{provenance}\r\n"
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"{filename}\"\r\n"
        "Content-Type: chemical/x-pdb\r\n\r\n"
        + source_text + "\r\n"
        + f"--{boundary}--\r\n").encode()
    status, body, _ = send(base, "/api/uploads", "POST", upload_body,
                           f"multipart/form-data; boundary={boundary}", base)
    require(status == 200, f"Declared local upload failed: {body!r}")
    return json.loads(body)["uploadToken"]


def provider_bytes(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "ProteinInMembraneSystem/0.1"})
    with urlopen(request, timeout=60) as response:
        return response.read()


def staged_source_bytes(workspace: Path, known: set[Path]) -> tuple[bytes, set[Path]]:
    files = set((workspace / "intake").glob("*/*.cif"))
    new_files = files - known
    require(len(new_files) == 1, "One exact remote coordinate asset must be staged")
    return next(iter(new_files)).read_bytes(), files


def main() -> None:
    host = ROOT / "out" / "host" / "ProteinInMembrane.Host.dll"
    worker = ROOT / "out" / "python" / "bin" / "python"
    catalogue = ROOT / "config" / "policies" / "protein-slice1.json"
    for file in (host, worker, catalogue, ROOT / "out" / "host" / "wwwroot" / "index.html"):
        require(file.is_file(), f"Build or policy prerequisite is missing: {file}")
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        with socket.socket() as probe:
            probe.bind(("127.0.0.1", 0))
            port = probe.getsockname()[1]
        base = f"http://127.0.0.1:{port}"
        environment = dict(os.environ,
                           PIM_PORT=str(port),
                           PIM_WORKSPACE_ROOT=str(directory / "workspace"),
                           PIM_POLICY_CATALOGUE=str(catalogue),
                           PIM_PPM_EXECUTABLE=str(directory / "absent-ppm"))
        log_path = directory / "host.log"
        with log_path.open("wb") as log:
            process = subprocess.Popen([str(ROOT / "scripts" / "start-local.sh")], cwd=ROOT,
                                       env=environment, stdout=log, stderr=subprocess.STDOUT)
            try:
                for _ in range(100):
                    if process.poll() is not None:
                        raise AssertionError("Host exited during startup: " + log_path.read_text())
                    try:
                        status, body, _ = send(base, "/api/state")
                        if status == 200:
                            break
                    except URLError:
                        time.sleep(0.1)
                else:
                    raise AssertionError("Host did not become ready: " + log_path.read_text())
                initial = json.loads(body)
                launch_log = log_path.read_text()
                require("Opening local protein-in-membrane workspace" in launch_log and
                        "PPM 2.0 is unavailable" in launch_log,
                        "The supported launcher must start slice 1 without the placement tool")
                require(initial["revision"] == 0, "An isolated host must start with revision 0")
                require(initial["protein"] is None and initial["placement"] is None,
                        "No scientific result may exist before selecting a source")
                require(not any("policy catalogue" in notice["message"].lower()
                                for notice in initial["notices"]),
                        "The identified slice-1 policy catalogue must load")

                # Foreign Origin never reaches the product root.
                status, _, _ = command(base, "searchSource", {"query": "alanine"}, 0,
                                       origin="http://foreign.example")
                require(status == 403, "A foreign Origin must be refused")
                require(json.loads(send(base, "/api/state")[1])["revision"] == 0,
                        "A refused foreign command must not change the account")

                token = upload_pdb(base, two_alanines(1, "A") + "END\n",
                                   "experimental", "two-alanines.pdb")

                status, body, _ = command(base, "selectSource", {"uploadToken": token}, 0)
                require(status == 200, f"Exact uploaded source was not inspected: {body!r}")
                selected = json.loads(body)
                require(selected["study"]["selectedSourceId"].startswith("upload:"),
                        "Selected source must retain local upload identity")
                require(selected["study"]["uploadProvenance"] == "experimental",
                        "Researcher-declared origin must remain a declaration")
                require(selected["sourcePrediction"] is None,
                        "An uploaded B-factor must not become verified prediction confidence")
                require(len(selected["sourceModels"]) == 1 and
                        selected["sourceModels"][0]["atomCount"] == 11 and
                        selected["sourceModels"][0]["chains"][0]["name"] == "A",
                        "The exact authored model and chain must be observed")
                require(selected["protein"] is None,
                        "Coordinates alone must not establish a prepared protein")

                status, body, _ = command(base, "selectProteinModel", {
                    "modelIndex": 0, "biologicalAssemblyId": None,
                    "chains": [{"sourceChain": "A", "copyId": "A"}],
                    "partners": [], "alternateLocations": [],
                }, selected["revision"])
                require(status == 200, f"Selected protein was not assessed: {body!r}")
                assessed = json.loads(body)
                protein = assessed["protein"]
                require(protein is not None and protein["status"] == "assessed",
                        f"Expected a corresponding Assessed Prepared Protein, got {protein!r}")
                require(protein["atomCount"] == 23 and protein["geometry"]["standing"] == "Observed",
                        "The actual prepared candidate and its geometry must be observed")
                require(assessed["study"]["modelIndex"] == 0 and
                        assessed["study"]["chainIds"] == ["A"],
                        "The assessed protein must preserve exact researcher selection")
                require(assessed["placement"] is None and assessed["stages"] == [],
                        "Protein preparation must not imply placement or a completed stage")

                status, _, _ = command(base, "selectSource", {"uploadToken": token}, 0)
                require(status == 422, "A stale revision must not repeat a source decision")
                require(json.loads(send(base, "/api/state")[1])["protein"]["status"] == "assessed",
                        "A stale request must not replace the assessed protein")

                status, body, _ = command(base, "selectInspectionSubject",
                                          {"subjectId": protein["subjectId"]}, assessed["revision"])
                require(status == 200, f"The assessed subject is not inspectable: {body!r}")
                inspected = json.loads(body)["inspection"]
                require(inspected is not None and inspected["subjectId"] == protein["subjectId"],
                        "The inspection context must identify this assessed protein")
                require(bool(inspected["structureUrl"]),
                        "The inspected protein must have a linked molecular structure")
                structure_status, structure_bytes, _ = send(base, inspected["structureUrl"])
                require(structure_status == 200 and b"ALA" in structure_bytes,
                        "The official host must serve the prepared molecular artifact")
                structure_token = urlsplit(inspected["structureUrl"]).path.rsplit("/", 1)[1]
                atom_path = (f"/api/inspection/atoms/{protein['subjectId']}/"
                             f"{structure_token}/0")
                atom_status, atom_body, _ = send(base, atom_path)
                require(atom_status == 200, "The first visible prepared atom needs verified identity")
                selected_atom = json.loads(atom_body)
                require(selected_atom["subjectId"] == protein["subjectId"] and
                        selected_atom["studyRevisionId"] == inspected["studyRevisionId"] and
                        selected_atom["structureToken"] == structure_token and
                        selected_atom["atomSiteIndex"] == 0 and
                        selected_atom["atom"]["resultAtomIndex"] == 0 and
                        selected_atom["atom"]["moleculeRole"] == "protein",
                        "The atom pick must resolve the exact selected protein correspondence")
                wrong_atom_path = (f"/api/inspection/atoms/wrong-subject/"
                                   f"{structure_token}/0")
                require(send(base, wrong_atom_path)[0] == 404 and
                        send(base, atom_path + "999999")[0] == 404,
                        "An unrelated subject or impossible row must have no atom identity")
                prepared_paths = list((directory / "workspace" / "protein-preparation").glob(
                    "*/prepared-protein.pdb"))
                require(len(prepared_paths) == 1 and prepared_paths[0].read_bytes() == structure_bytes,
                        "The selected structure URL must serve the exact prepared artifact")
                prepared_paths[0].write_bytes(structure_bytes + b"REMARK changed after selection\n")
                changed_status, _, _ = send(base, inspected["structureUrl"])
                require(changed_status == 404,
                        "A selected structure URL must not serve changed bytes under its old identity")
                require(send(base, atom_path)[0] == 404,
                        "A changed coordinate artifact must invalidate its atom-pick identity")
                prepared_paths[0].write_bytes(structure_bytes)
                restored_status, restored_bytes, _ = send(base, inspected["structureUrl"])
                require(restored_status == 200 and restored_bytes == structure_bytes,
                        "The original selected structure bytes must remain addressable after restoration")
                print("PASS: loopback upload, exact selection, real preparation, inspection, stale revision and Origin refusal")

                # A second source makes model, partner and alternate-location
                # choices distinguishable through the official actor route.
                alternate_cb = atom_line(5, "CB", "ALA", "B", 1,
                                         (2, -0.77, -1.2), "C")
                other_cb = atom_line(12, "CB", "ALA", "B", 1,
                                     (2, -0.77, 1.2), "C")
                alternate_model = two_alanines(2, "B").replace(
                    alternate_cb, alternate_cb[:16] + "A" + alternate_cb[17:] +
                    other_cb[:16] + "B" + other_cb[17:])
                partner = "HETATM   13  C1  GOL C   3       9.000   4.000   0.000  1.00 20.00           C\n"
                alternate_model = alternate_model.replace("ENDMDL\n", partner + "ENDMDL\n")
                mixed_token = upload_pdb(base, two_alanines(1, "A") +
                                         alternate_model + "END\n", "unknown",
                                         "models-alternate-partner.pdb")
                status, body, _ = command(base, "selectSource",
                                          {"uploadToken": mixed_token},
                                          json.loads(send(base, "/api/state")[1])["revision"])
                require(status == 200, f"Multi-model source inspection failed: {body!r}")
                mixed = json.loads(body)
                require(mixed["protein"] is None and
                        mixed["study"]["uploadProvenance"] == "unknown" and
                        mixed["sourcePrediction"] is None,
                        "A new unknown-origin upload withdraws the earlier assessed protein")
                require([model["index"] for model in mixed["sourceModels"]] == [0, 1],
                        "Two distinguishable source models must be observed")
                alternate_residue = next(residue for residue in mixed["sourceModels"][1]["residues"]
                                         if residue["address"]["chain"] == "B" and
                                         residue["address"]["residue"] == 1)
                require(alternate_residue["alternateLocations"] == ["A", "B"],
                        "The second model must expose its alternate coordinates")
                observed_partner = mixed["sourceModels"][1]["partners"][0]
                require(observed_partner["sourceId"] == "1:C:3::GOL",
                        "The second model must expose the exact nonprotein partner")
                selection = {
                    "modelIndex": 1, "biologicalAssemblyId": None,
                    "chains": [{"sourceChain": "B", "copyId": "B"}],
                    "partners": [{"sourceId": observed_partner["sourceId"],
                                  "retain": False, "reason": "Exclude this observed additive"}],
                    "alternateLocations": [{"residue": alternate_residue["address"],
                                            "altloc": "A", "decisionId": ""}],
                }
                status, _, _ = command(base, "selectProteinModel",
                                       {**selection, "modelIndex": 9}, mixed["revision"])
                require(status == 422, "An unobserved source model must be refused")
                status, _, _ = command(base, "selectProteinModel",
                                       {**selection, "partners": []}, mixed["revision"])
                require(status == 422, "An undecided observed partner must be refused")
                status, body, _ = command(base, "selectProteinModel", selection,
                                          mixed["revision"])
                require(status == 200, f"Exact second-model selection failed: {body!r}")
                proposed = json.loads(body)
                require(proposed["study"]["modelIndex"] == 1 and
                        proposed["study"]["chainIds"] == ["B"],
                        "The adopted protein must retain the chosen second model and chain")
                changes = proposed["protein"]["changes"]
                require(len(changes) == 1 and changes[0]["kind"] == "alternateLocation" and
                        proposed["protein"]["status"] != "assessed",
                        "An alternate-location choice must remain a review proposal")
                proposal = changes[0]
                status, body, _ = command(base, "selectInspectionSubject",
                                          {"subjectId": proposal["id"]}, proposed["revision"])
                require(status == 200 and
                        json.loads(body)["inspection"]["subjectId"] == proposal["id"],
                        "The exact alternate-location proposal must be inspectable")
                status, body, _ = command(base, "approvePreparationChange", {
                    "proposalId": proposal["id"], "approve": True, "rationale": "",
                }, json.loads(body)["revision"])
                require(status == 200, f"Reviewed alternate choice was not accepted: {body!r}")
                selected_alternate = json.loads(body)
                require(selected_alternate["protein"]["status"] == "assessed" and
                        selected_alternate["placement"] is None and
                        selected_alternate["stages"] == [],
                        "Approved exact alternate coordinates must yield only protein preparation")
                print("PASS: exact second model, chain, partner exclusion and reviewed alternate location")

                missing_cb = "".join(line for line in two_alanines(1, "A").splitlines(
                    keepends=True) if not (line.startswith("ATOM") and
                        line[12:16].strip() == "CB" and line[22:26].strip() == "2")) + "END\n"
                repair_token = upload_pdb(base, missing_cb, "experimental",
                                          "missing-sidechain.pdb")
                status, body, _ = command(base, "selectSource",
                                          {"uploadToken": repair_token},
                                          selected_alternate["revision"])
                require(status == 200, f"Repair source inspection failed: {body!r}")
                repair_source = json.loads(body)
                status, body, _ = command(base, "selectProteinModel", {
                    "modelIndex": 0, "biologicalAssemblyId": None,
                    "chains": [{"sourceChain": "A", "copyId": "A"}],
                    "partners": [], "alternateLocations": [],
                }, repair_source["revision"])
                require(status == 200, f"Repair candidate selection failed: {body!r}")
                repair_proposed = json.loads(body)
                changes = repair_proposed["protein"]["changes"]
                require(len(changes) == 1 and changes[0]["kind"] == "heavyAtom" and
                        repair_proposed["protein"]["status"] != "assessed",
                        "The missing CB must be a consequential bounded proposal")
                change_id = changes[0]["id"]
                status, _, _ = command(base, "approvePreparationChange", {
                    "proposalId": change_id, "approve": True, "rationale": "",
                }, repair_proposed["revision"])
                require(status == 422,
                        "A heavy-atom proposal must not be approved before linked inspection")
                status, body, _ = command(base, "selectInspectionSubject",
                                          {"subjectId": change_id},
                                          repair_proposed["revision"])
                require(status == 200 and
                        json.loads(body)["inspection"]["subjectId"] == change_id,
                        "The required heavy-atom region and evidence must be inspectable")
                status, body, _ = command(base, "approvePreparationChange", {
                    "proposalId": change_id, "approve": True, "rationale": "",
                }, json.loads(body)["revision"])
                require(status == 200, f"Reviewed exact heavy repair failed: {body!r}")
                repaired = json.loads(body)
                require(repaired["protein"]["status"] == "assessed" and
                        repaired["protein"]["atomCount"] == 23 and
                        repaired["protein"]["geometry"]["standing"] == "Observed" and
                        repaired["placement"] is None,
                        "The actual approved heavy-atom candidate must be reassessed as protein only")
                print("PASS: reviewed heavy-atom approval and real reassessed candidate")

                if "--remote" in sys.argv[1:]:
                    intake_files: set[Path] = set()
                    status, body, _ = command(base, "selectSource", {
                        "sourceKind": "rcsb", "exactIdentifier": "1CRN",
                    }, json.loads(send(base, "/api/state")[1])["revision"])
                    require(status == 200, f"Exact RCSB source retrieval failed: {body!r}")
                    rcsb = json.loads(body)
                    require(rcsb["study"]["selectedSourceId"] == "rcsb:1CRN" and
                            rcsb["study"]["selectedSourceKind"] == "rcsb" and
                            len(rcsb["sourceModels"]) > 0 and
                            rcsb["sourceModels"][0]["atomCount"] > 0,
                            "RCSB direct reference must identify the retrieved structure")
                    staged_rcsb, intake_files = staged_source_bytes(
                        directory / "workspace", intake_files)
                    require(hashlib.sha256(staged_rcsb).digest() == hashlib.sha256(
                        provider_bytes("https://files.rcsb.org/download/1CRN.cif")).digest(),
                        "Staged RCSB coordinates must match independently fetched exact entry bytes")
                    require(rcsb["protein"] is None,
                        "Changing source must withdraw the earlier assessed protein")

                    status, body, _ = command(base, "selectSource", {
                        "sourceKind": "alphafold", "exactIdentifier": "AF-P69905-F1",
                    }, rcsb["revision"])
                    require(status == 200, f"Exact AlphaFold source retrieval failed: {body!r}")
                    predicted = json.loads(body)
                    prediction = predicted["sourcePrediction"]
                    require(predicted["study"]["selectedSourceId"] == "alphafold:AF-P69905-F1" and
                            predicted["study"]["selectedSourceKind"] == "alphafold" and
                            prediction is not None and prediction["recordId"] == "AF-P69905-F1",
                            "AlphaFold coordinates and evidence must retain the selected model identity")
                    records = json.loads(provider_bytes(
                        "https://alphafold.ebi.ac.uk/api/prediction/P69905"))
                    matches = [record for record in records if
                               (record.get("modelEntityId") or record.get("entryId")) ==
                               "AF-P69905-F1"]
                    require(len(matches) == 1,
                            "Provider metadata must identify exactly one selected AlphaFold record")
                    staged_prediction, intake_files = staged_source_bytes(
                        directory / "workspace", intake_files)
                    require(hashlib.sha256(staged_prediction).digest() == hashlib.sha256(
                        provider_bytes(matches[0]["cifUrl"])).digest(),
                        "Staged prediction coordinates must match the exact metadata record")
                    pae_files = list((directory / "workspace" / "intake").glob("*/pae-*.json"))
                    require(len(pae_files) == 1 and hashlib.sha256(
                        pae_files[0].read_bytes()).digest() == hashlib.sha256(
                            provider_bytes(matches[0]["paeDocUrl"])).digest(),
                        "Staged PAE must belong to the selected prediction record")
                    require(prediction["paeStanding"] == "Observed" and
                            prediction["paeAxisResidueCount"] ==
                            len(prediction["localConfidence"]) > 0,
                            "The selected prediction's PAE axis and residue confidence must correspond")
                    require(predicted["protein"] is None,
                            "A predicted coordinate source alone is not an assessed protein")
                    print("PASS: live direct RCSB and AlphaFold references without discovery search")
            finally:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)


if __name__ == "__main__":
    main()
