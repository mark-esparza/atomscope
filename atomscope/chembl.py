"""Optional ChEMBL cross-reference: is this chemical a known bioactive/drug?"""

from __future__ import annotations

import urllib.parse

from .http_util import FetchError, fetch_json

_BASE = "https://www.ebi.ac.uk/chembl/api/data"

_PHASE = {
    0: "preclinical / research compound",
    1: "Phase 1 clinical",
    2: "Phase 2 clinical",
    3: "Phase 3 clinical",
    4: "approved drug",
}


def lookup_molecule(name: str) -> dict | None:
    """Best-effort ChEMBL lookup by name. Returns None if not found/unavailable."""
    try:
        q = urllib.parse.quote(name)
        url = f"{_BASE}/molecule/search?q={q}&format=json&limit=1"
        data = fetch_json(url)
    except FetchError:
        return None

    molecules = data.get("molecules") or []
    if not molecules:
        return None
    mol = molecules[0]
    max_phase = mol.get("max_phase")
    try:
        phase_num = int(float(max_phase)) if max_phase is not None else None
    except (TypeError, ValueError):
        phase_num = None

    return {
        "chembl_id": mol.get("molecule_chembl_id"),
        "pref_name": mol.get("pref_name"),
        "max_phase": phase_num,
        "development_status": _PHASE.get(phase_num, "unknown"),
        "url": (
            f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol.get('molecule_chembl_id')}/"
            if mol.get("molecule_chembl_id")
            else None
        ),
    }
