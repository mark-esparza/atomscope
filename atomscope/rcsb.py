"""RCSB Protein Data Bank client: structure files + entry metadata."""

from __future__ import annotations

import re

from .http_util import FetchError, fetch_json, fetch_text

_PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")


def normalize_pdb_id(pdb_id: str) -> str:
    pid = (pdb_id or "").strip().upper()
    if not _PDB_ID_RE.match(pid):
        raise FetchError(f"'{pdb_id}' is not a valid 4-character PDB ID")
    return pid


def fetch_structure(pdb_id: str) -> str:
    """Return the raw PDB-format text for an entry."""
    pid = normalize_pdb_id(pdb_id)
    return fetch_text(f"https://files.rcsb.org/download/{pid}.pdb")


def fetch_entry_metadata(pdb_id: str) -> dict:
    """Return a compact, UI-friendly metadata dict for a PDB entry."""
    pid = normalize_pdb_id(pdb_id)
    data = fetch_json(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")

    struct = data.get("struct") or {}
    entry_info = data.get("rcsb_entry_info") or {}
    accession = data.get("rcsb_accession_info") or {}
    exptl = data.get("exptl") or [{}]

    resolution = None
    res_list = entry_info.get("resolution_combined")
    if isinstance(res_list, list) and res_list:
        resolution = res_list[0]

    methods = [e.get("method") for e in exptl if e.get("method")]

    return {
        "pdb_id": pid,
        "title": struct.get("title"),
        "experimental_method": ", ".join(methods) if methods else None,
        "resolution_A": resolution,
        "deposited": accession.get("initial_release_date"),
        "polymer_entity_count": entry_info.get("polymer_entity_count"),
        "deposited_atom_count": entry_info.get("deposited_atom_count"),
        "deposited_model_count": entry_info.get("deposited_model_count"),
        "molecular_weight_kDa": entry_info.get("molecular_weight"),
        "nonpolymer_count": entry_info.get("nonpolymer_entity_count"),
    }


def search_by_name(query: str, limit: int = 10) -> list[dict]:
    """Full-text search the PDB, returning [{id, score}] ranked entries."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    import json
    import urllib.parse

    url = (
        "https://search.rcsb.org/rcsbsearch/v2/query?json="
        + urllib.parse.quote(json.dumps(payload))
    )
    data = fetch_json(url)
    results = []
    for item in data.get("result_set", []):
        results.append({"pdb_id": item.get("identifier"), "score": item.get("score")})
    return results
