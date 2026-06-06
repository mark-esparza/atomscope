"""Geometric pocket / cavity detection (LIGSITE algorithm, pure Python).

Lets AtomScope dock into *apo* structures that have no bound reference ligand.

Method (LIGSITE, Hendlich et al. 1997):
  1. Lay a cubic grid over the protein.
  2. Mark grid points inside the protein's van der Waals volume as "protein".
  3. For each free grid point, scan 7 directions (3 axes + 4 cube diagonals).
     A direction is a "protein-solvent-protein" (PSP) event if protein is found
     within range on BOTH sides. Points enclosed in >= threshold directions are
     pocket points.
  4. Cluster adjacent pocket points; rank clusters by size (cavity volume).

Output per pocket: center, volume estimate, enclosure, and lining residues.
Pure geometry, no dependencies. Approximate / research-only.
"""

from __future__ import annotations

import math
from collections import deque

from .pdbparse import Structure

_RVDW = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "H": 1.20, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "ZN": 1.39, "FE": 1.40, "MG": 1.73, "CA": 2.00, "NA": 2.27, "K": 2.75,
}
_DEFAULT_RVDW = 1.70

MARGIN = 5.0          # box padding around protein (A)
FILL = 0.8            # radius inflation to close grid-discretization voids (A)
SCAN_RANGE = 8        # max grid steps to look for an enclosing protein point
PSP_THRESHOLD = 5     # min enclosed directions (of 7) to call a pocket point
MIN_POCKET_POINTS = 25  # discard clusters smaller than this (~ volume floor)
MAX_GRID_POINTS = 500_000

# 7 LIGSITE scan directions (axes + cube diagonals), as grid-index steps.
_DIRS = [
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
]

# Residue classes for pocket physicochemistry (Gap 2: pocket vs ligandable vs druggable).
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
POLAR_RES = {"SER", "THR", "ASN", "GLN", "TYR", "GLY"}
CHARGED_RES = {"ASP", "GLU", "LYS", "ARG", "HIS"}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}


def _rvdw(el: str) -> float:
    return _RVDW.get(el.upper(), _DEFAULT_RVDW)


class _CellGrid:
    def __init__(self, atoms, cell):
        self.cell = cell
        self.b = {}
        for a in atoms:
            k = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
            self.b.setdefault(k, []).append(a)

    def near(self, x, y, z):
        cx, cy, cz = int(x // self.cell), int(y // self.cell), int(z // self.cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.b.get((cx + dx, cy + dy, cz + dz))
                    if bucket:
                        yield from bucket


def detect_pockets(structure: Structure, max_pockets: int = 6) -> list[dict]:
    atoms = [a for a in structure.protein_atoms if a.element != "H"]
    if len(atoms) < 20:
        return []

    xs = [a.x for a in atoms]
    ys = [a.y for a in atoms]
    zs = [a.z for a in atoms]
    minx, miny, minz = min(xs) - MARGIN, min(ys) - MARGIN, min(zs) - MARGIN
    maxx, maxy, maxz = max(xs) + MARGIN, max(ys) + MARGIN, max(zs) + MARGIN

    # Choose spacing so the grid stays within a sane size.
    spacing = 1.0
    while True:
        nx = int((maxx - minx) / spacing) + 1
        ny = int((maxy - miny) / spacing) + 1
        nz = int((maxz - minz) / spacing) + 1
        if nx * ny * nz <= MAX_GRID_POINTS:
            break
        spacing += 0.2

    cell = _CellGrid(atoms, cell=5.0)
    nynz = ny * nz
    total = nx * nynz
    protein = bytearray(total)   # 1 = inside protein volume
    near = bytearray(total)      # 1 = within 5 A of protein (limits scan region)

    for i in range(nx):
        px = minx + i * spacing
        for j in range(ny):
            py = miny + j * spacing
            base = i * nynz + j * nz
            for k in range(nz):
                pz = minz + k * spacing
                is_prot = 0
                is_near = 0
                for a in cell.near(px, py, pz):
                    dx = px - a.x
                    dy = py - a.y
                    dz = pz - a.z
                    d2 = dx * dx + dy * dy + dz * dz
                    if d2 <= 25.0:  # 5 A
                        is_near = 1
                        r = _rvdw(a.element) + FILL
                        if d2 <= r * r:
                            is_prot = 1
                            break
                idx = base + k
                protein[idx] = is_prot
                near[idx] = is_near

    def enclosed(i, j, k, di, dj, dk):
        ii, jj, kk = i, j, k
        for _ in range(SCAN_RANGE):
            ii += di
            jj += dj
            kk += dk
            if ii < 0 or jj < 0 or kk < 0 or ii >= nx or jj >= ny or kk >= nz:
                return False
            if protein[ii * nynz + jj * nz + kk]:
                return True
        return False

    pocket = bytearray(total)
    for i in range(nx):
        for j in range(ny):
            base = i * nynz + j * nz
            for k in range(nz):
                idx = base + k
                if protein[idx] or not near[idx]:
                    continue
                psp = 0
                for (di, dj, dk) in _DIRS:
                    if enclosed(i, j, k, di, dj, dk) and enclosed(i, j, k, -di, -dj, -dk):
                        psp += 1
                if psp >= PSP_THRESHOLD:
                    pocket[idx] = psp  # store depth (5..7) for enclosure scoring

    clusters = _cluster(pocket, nx, ny, nz, nynz, nz)
    clusters.sort(key=len, reverse=True)

    results = []
    voxel_vol = spacing ** 3
    for ci, pts in enumerate(clusters[:max_pockets]):
        if len(pts) < MIN_POCKET_POINTS:
            break
        sx = sy = sz = 0.0
        psp_sum = 0
        coords = []
        for (i, j, k) in pts:
            wx = minx + i * spacing
            wy = miny + j * spacing
            wz = minz + k * spacing
            sx += wx
            sy += wy
            sz += wz
            psp_sum += pocket[i * nynz + j * nz + k]
            coords.append((wx, wy, wz))
        n = len(pts)
        center = (sx / n, sy / n, sz / n)
        mean_psp = psp_sum / n
        lining = _lining_residues(coords, cell)
        volume = n * voxel_vol
        assessment = _assess(volume, mean_psp, lining)
        results.append(
            {
                "index": ci,
                "center": [round(c, 2) for c in center],
                "n_points": n,
                "volume_A3": round(volume, 1),
                "enclosure": round(mean_psp, 2),
                "lining_residues": lining,
                "lining_residue_count": len(lining),
                "score": assessment["druggability_score"],
                "tier": assessment["tier"],
                "subscores": assessment["subscores"],
                "composition": assessment["composition"],
            }
        )

    results.sort(key=lambda p: p["score"], reverse=True)
    for new_index, p in enumerate(results):
        p["index"] = new_index
    return results


def _cluster(pocket, nx, ny, nz, nynz, sz):
    seen = bytearray(len(pocket))
    clusters = []
    neighbors = [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
        if not (dx == 0 and dy == 0 and dz == 0)
    ]
    for i in range(nx):
        for j in range(ny):
            base = i * nynz + j * sz
            for k in range(nz):
                idx = base + k
                if not pocket[idx] or seen[idx]:
                    continue
                comp = []
                q = deque([(i, j, k)])
                seen[idx] = 1
                while q:
                    ci, cj, ck = q.popleft()
                    comp.append((ci, cj, ck))
                    for (dx, dy, dz) in neighbors:
                        ni, nj, nk = ci + dx, cj + dy, ck + dz
                        if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
                            nidx = ni * nynz + nj * sz + nk
                            if pocket[nidx] and not seen[nidx]:
                                seen[nidx] = 1
                                q.append((ni, nj, nk))
                clusters.append(comp)
    return clusters


def _lining_residues(coords, cell, cutoff=4.0):
    cut2 = cutoff * cutoff
    seen = {}
    # Subsample pocket points for speed on large cavities.
    step = max(1, len(coords) // 200)
    for (wx, wy, wz) in coords[::step]:
        for a in cell.near(wx, wy, wz):
            dx = wx - a.x
            dy = wy - a.y
            dz = wz - a.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= cut2:
                rid = f"{a.chain}/{a.res_name}{a.res_seq}"
                cur = seen.get(rid)
                if cur is None:
                    seen[rid] = {
                        "res_id": rid,
                        "res_name": a.res_name,
                        "res_seq": a.res_seq,
                        "chain": a.chain,
                        "_d2": d2,
                    }
                elif d2 < cur["_d2"]:
                    cur["_d2"] = d2
    # Most-contacting (closest) residues first, so the displayed subset is representative.
    out = sorted(seen.values(), key=lambda r: r["_d2"])
    for r in out:
        r.pop("_d2", None)
    return out


def _assess(volume, mean_psp, lining):
    """Classify a cavity as pocket / ligandable / druggable with sub-scores.

    Addresses the literature's Gap 2: a detected *pocket* is geometric; a
    *ligandable* pocket is large and enclosed enough to bind a small molecule;
    a *druggable* pocket additionally has favourable (hydrophobic, enclosed,
    not-too-polar) chemistry. Heuristic — proxies SiteMap/fpocket druggability,
    not a trained model.
    """
    total = len(lining) or 1
    hydrophobic = sum(1 for r in lining if r["res_name"] in HYDROPHOBIC_RES)
    polar = sum(1 for r in lining if r["res_name"] in POLAR_RES)
    charged = sum(1 for r in lining if r["res_name"] in CHARGED_RES)
    aromatic = sum(1 for r in lining if r["res_name"] in AROMATIC_RES)

    hydrophobicity = hydrophobic / total
    polarity = (polar + charged) / total
    aromaticity = aromatic / total

    vol_term = min(volume / 400.0, 1.0)               # saturates ~400 A^3
    encl_term = max(0.0, min((mean_psp - 4.0) / 3.0, 1.0))  # PSP 5..7 -> 0.33..1

    # Druggability index (0-100): size + enclosure + hydrophobicity, polarity penalty.
    drug = 0.30 * vol_term + 0.30 * encl_term + 0.30 * hydrophobicity - 0.15 * polarity
    drug = round(100 * max(0.0, min(drug + 0.15, 1.0)), 1)

    # Tier thresholds (heuristic).
    if volume >= 80 and encl_term >= 0.4:
        tier = "druggable" if (drug >= 55 and hydrophobicity >= 0.40) else "ligandable"
    else:
        tier = "pocket"

    return {
        "tier": tier,
        "druggability_score": drug,
        "subscores": {
            "volume": round(vol_term, 2),
            "enclosure": round(encl_term, 2),
            "hydrophobicity": round(hydrophobicity, 2),
            "polarity": round(polarity, 2),
            "aromaticity": round(aromaticity, 2),
        },
        "composition": {
            "hydrophobic": hydrophobic,
            "polar": polar,
            "charged": charged,
            "aromatic": aromatic,
            "total": len(lining),
        },
    }
