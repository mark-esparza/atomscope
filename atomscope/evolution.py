"""Evolutionary conservation analysis from a protein-family alignment.

Pipeline (dependency-free):
  1. UniProt accession -> Pfam family (InterPro API) -> family alignment.
  2. Per-column Shannon-entropy conservation + consensus residue.
  3. Extract the loaded structure's sequence (with residue numbering).
  4. Needleman-Wunsch align the target to the family consensus to map each
     structure residue onto an alignment column, giving per-residue conservation.
  5. Score each detected pocket by the mean conservation of its lining residues.

Conservation here is a family-MSA signal (how invariant a position is across
homologs), not a phylogenetic reconstruction. Research-only.
"""

from __future__ import annotations

import gzip
import math
import urllib.parse

from .http_util import FetchError, fetch_bytes, fetch_json

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
_AA = "ACDEFGHIKLMNPQRSTVWY"
_LOG20 = math.log2(20)


# ---------------- alignment acquisition ----------------

def pfam_for_uniprot(accession: str) -> dict | None:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"
    try:
        data = fetch_json(url)
    except FetchError:
        return None
    results = data.get("results") or []
    if not results:
        return None
    meta = results[0].get("metadata", {})
    return {"pfam": meta.get("accession"), "name": meta.get("name")}


def fetch_family_alignment(pfam_id: str, kind: str = "alignment:seed") -> list[tuple]:
    """Return [(name, aligned_seq)] for a Pfam family (Stockholm, gap-normalised)."""
    url = (
        f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/"
        f"?annotation={urllib.parse.quote(kind)}"
    )
    raw = fetch_bytes(url)
    try:
        text = gzip.decompress(raw).decode("utf-8", "replace")
    except (OSError, EOFError):
        text = raw.decode("utf-8", "replace")

    rows: dict[str, list] = {}
    order: list[str] = []
    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        name, seq = parts
        if name not in rows:
            rows[name] = []
            order.append(name)
        rows[name].append(seq.strip())
    aln = []
    for name in order:
        s = "".join(rows[name]).upper().replace(".", "-")
        aln.append((name, s))
    # Keep only rows that share the modal alignment width (guards against junk).
    if not aln:
        return []
    from collections import Counter

    width = Counter(len(s) for _, s in aln).most_common(1)[0][0]
    return [(n, s) for (n, s) in aln if len(s) == width]


# ---------------- conservation ----------------

def column_conservation(alignment: list[tuple]):
    """Return (conservation[], consensus[], occupancy[]) per alignment column."""
    if not alignment:
        return [], [], []
    width = len(alignment[0][1])
    n = len(alignment)
    conservation, consensus, occupancy = [], [], []
    for col in range(width):
        counts: dict[str, int] = {}
        nongap = 0
        for _, seq in alignment:
            c = seq[col]
            if c == "-" or c not in _AA:
                continue
            counts[c] = counts.get(c, 0) + 1
            nongap += 1
        if nongap == 0:
            conservation.append(0.0)
            consensus.append("-")
            occupancy.append(0.0)
            continue
        entropy = 0.0
        for c, k in counts.items():
            p = k / nongap
            entropy -= p * math.log2(p)
        cons = 1.0 - entropy / _LOG20
        conservation.append(max(0.0, min(cons, 1.0)))
        consensus.append(max(counts, key=counts.get))
        occupancy.append(nongap / n)
    return conservation, consensus, occupancy


# ---------------- target sequence from structure ----------------

def structure_sequence(structure):
    """Pick the longest protein chain; return (one_letter_seq, [(chain,res_seq,res_name)])."""
    by_chain: dict[str, dict] = {}
    for a in structure.protein_atoms:
        by_chain.setdefault(a.chain, {})
        if a.res_seq not in by_chain[a.chain]:
            by_chain[a.chain][a.res_seq] = a.res_name
    if not by_chain:
        return "", []
    chain = max(by_chain, key=lambda c: len(by_chain[c]))
    residues = sorted(by_chain[chain].items())  # (res_seq, res_name)
    seq = []
    keys = []
    for res_seq, res_name in residues:
        seq.append(AA3TO1.get(res_name, "X"))
        keys.append((chain, res_seq, res_name))
    return "".join(seq), keys


# ---------------- Needleman-Wunsch (global) ----------------

def _nw_align(a: str, b: str):
    """Global alignment; returns list of (i|None, j|None) index pairs."""
    GAP = -4
    n, m = len(a), len(b)
    # Score matrix rows kept as lists; traceback via pointer matrix.
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * GAP
    for j in range(1, m + 1):
        score[0][j] = j * GAP
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = score[i], score[i - 1]
        for j in range(1, m + 1):
            sub = 5 if ai == b[j - 1] else -1
            row[j] = max(prev[j - 1] + sub, prev[j] + GAP, row[j - 1] + GAP)
    # Traceback.
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        cur = score[i][j]
        sub = 5 if a[i - 1] == b[j - 1] else -1
        if cur == score[i - 1][j - 1] + sub:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif cur == score[i - 1][j] + GAP:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    while i > 0:
        pairs.append((i - 1, None)); i -= 1
    while j > 0:
        pairs.append((None, j - 1)); j -= 1
    pairs.reverse()
    return pairs


# ---------------- orchestration ----------------

def analyze(structure, uniprots: list[str]) -> dict | None:
    """Full conservation analysis for a loaded structure. None if unavailable."""
    fam = None
    acc_used = None
    for acc in uniprots or []:
        fam = pfam_for_uniprot(acc)
        if fam and fam.get("pfam"):
            acc_used = acc
            break
    if not fam or not fam.get("pfam"):
        return None

    try:
        alignment = fetch_family_alignment(fam["pfam"])
    except FetchError:
        return None
    if len(alignment) < 5:
        return None

    cons, consensus, occ = column_conservation(alignment)
    # Core columns = reasonably occupied positions.
    core_cols = [i for i in range(len(cons)) if occ[i] >= 0.5]
    core_consensus = "".join(consensus[i] for i in core_cols)
    core_cons = [cons[i] for i in core_cols]
    if not core_consensus:
        return None

    target_seq, keys = structure_sequence(structure)
    if not target_seq:
        return None

    pairs = _nw_align(target_seq, core_consensus)
    residues = []
    mapped = 0
    cons_by_key: dict[tuple, float] = {}
    for ti, cj in pairs:
        if ti is None:
            continue
        chain, res_seq, res_name = keys[ti]
        if cj is None:
            residues.append(
                {"chain": chain, "res_seq": res_seq, "res_name": res_name,
                 "conservation": None}
            )
            continue
        c = round(core_cons[cj], 3)
        mapped += 1
        cons_by_key[(chain, res_seq)] = c
        residues.append(
            {"chain": chain, "res_seq": res_seq, "res_name": res_name,
             "conservation": c, "consensus": core_consensus[cj]}
        )

    ranked = sorted(
        (r for r in residues if r["conservation"] is not None),
        key=lambda r: -r["conservation"],
    )
    top_conserved = [
        {"res": f"{r['res_name']}{r['res_seq']}", "chain": r["chain"],
         "conservation": r["conservation"]}
        for r in ranked[:15]
    ]

    return {
        "pfam": fam["pfam"],
        "family_name": fam.get("name"),
        "uniprot": acc_used,
        "n_sequences": len(alignment),
        "target_length": len(target_seq),
        "mapped_residues": mapped,
        "coverage": round(mapped / len(target_seq), 2) if target_seq else 0,
        "residues": residues,
        "top_conserved": top_conserved,
        "_cons_by_key": cons_by_key,  # internal: for pocket scoring
    }


def annotate_pockets(evo: dict, pockets: list) -> list:
    """Add mean lining-residue conservation + a label to each pocket."""
    cons_by_key = evo.get("_cons_by_key", {})
    out = []
    for p in pockets:
        vals = []
        for r in p.get("lining_residues", []):
            c = cons_by_key.get((r["chain"], r["res_seq"]))
            if c is not None:
                vals.append(c)
        if vals:
            mean = sum(vals) / len(vals)
            label = (
                "evolutionarily conserved" if mean >= 0.55
                else "moderately conserved" if mean >= 0.35
                else "variable"
            )
        else:
            mean = None
            label = "unmapped"
        out.append(
            {
                "index": p["index"],
                "tier": p.get("tier"),
                "volume_A3": p.get("volume_A3"),
                "mean_conservation": round(mean, 3) if mean is not None else None,
                "conserved_residues": len(vals),
                "label": label,
            }
        )
    return out
