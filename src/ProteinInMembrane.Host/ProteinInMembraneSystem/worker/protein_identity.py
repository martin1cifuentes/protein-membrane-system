"""Source residue and atom identity mapping shared by preparation and stage measurements."""

from __future__ import annotations

from typing import Any

from ProteinInMembraneSystem.worker.exchange import WorkError


def _icode(value: Any) -> str:
    return str(value or "").replace("\0", "").strip()


def _address_from_key(key: tuple[str, int, str], chain_map: list[dict[str, str]], model: int) -> dict[str, Any]:
    match = next((entry for entry in chain_map if entry["preparedChain"] == key[0]), None)
    if match is None:
        raise WorkError("correspondenceFailed", "A result atom belongs to an unknown selected chain")
    return {"model": model, "chain": match["sourceChain"], "residue": key[1],
            "insertionCode": key[2], "copyId": match["copyId"]}


def _atom_address(key: tuple[str, int, str, str], chain_map: list[dict[str, str]], model: int) -> dict[str, Any]:
    return {"residue": _address_from_key(key[:3], chain_map, model), "atomName": key[3]}
