"""Atomic-level protein-ligand interaction profiler.

Given a parsed Structure and a chosen hetero Component (ligand, ion, or metal),
this computes heavy-atom contacts between the component and the protein and
classifies each into an interaction type using geometric criteria. This mirrors
the approach of tools like PLIP, simplified to run dependency-free on
experimental coordinates (no explicit hydrogens assumed).

Criteria are research heuristics, not force-field calculations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .pdbparse import Component, Structure

# Distance cutoffs in Angstroms.
HB_MAX = 3.6          # N/O .. N/O polar contact (heavy-atom)
HB_S_MAX = 3.9        # contacts involving sulfur
SALT_MAX = 4.0        # charged group .. charged group
HYDRO_MIN = 2.8       # avoid covalent-bond-like distances
HYDRO_MAX = 4.0       # C .. C hydrophobic
METAL_MAX = 2.9       # metal .. coordinating atom
ARO_CENTROID_MAX = 4.5  # aromatic ring centroid .. nearest ligand heavy atom

POLAR_ELEMENTS = {"N", "O"}

# Charged protein side-chain atoms.
PROT_POS_ATOMS = {
    ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"),
    ("LYS", "NZ"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}
PROT_NEG_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"),
    ("GLU", "OE1"), ("GLU", "OE2"),
}

# Aromatic ring atoms by residue, for centroid-based pi detection.
AROMATIC_RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

INTERACTION_TYPES = [
    "metal_coordination",
    "salt_bridge",
    "hydrogen_bond",
    "hydrophobic",
    "aromatic",
]


@dataclass
class _Grid:
    cell: float
    buckets: dict

    @classmethod
    def build(cls, atoms, cell: float):
        buckets: dict = {}
        for a in atoms:
            key = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
            buckets.setdefault(key, []).append(a)
        return cls(cell, buckets)

    def neighbors(self, x: float, y: float, z: float):
        cx, cy, cz = int(x // self.cell), int(y // self.cell), int(z // self.cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.buckets.get((cx + dx, cy + dy, cz + dz))
                    if bucket:
                        yield from bucket


def _dist(a, b) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _atom_json(a, include_res: bool):
    d = {
        "name": a.name,
        "element": a.element,
        "serial": a.serial,
        "xyz": [round(a.x, 3), round(a.y, 3), round(a.z, 3)],
    }
    if include_res:
        d.update(
            res_name=a.res_name,
            res_seq=a.res_seq,
            chain=a.chain,
            res_id=f"{a.chain}/{a.res_name}{a.res_seq}",
        )
    return d


def _classify(lig, prot, d: float, is_metal_comp: bool) -> str | None:
    le, pe = lig.element, prot.element
    pkey = (prot.res_name, prot.name)

    # Metal coordination (component is a metal ion).
    if is_metal_comp and le in {"NA", "K", "MG", "CA", "MN", "FE", "CO", "NI",
                                "CU", "ZN", "MO", "W", "CD", "HG", "PT", "AU",
                                "AG", "PB", "BA", "SR", "LI", "AL", "V", "CR"}:
        if pe in {"O", "N", "S"} and d <= METAL_MAX:
            return "metal_coordination"
        return None

    # Salt bridge (charged .. charged), heuristic on protonation.
    if d <= SALT_MAX:
        if pkey in PROT_NEG_ATOMS and le == "N":
            return "salt_bridge"
        if pkey in PROT_POS_ATOMS and le == "O":
            return "salt_bridge"

    # Hydrogen bond / polar contact (heavy atom).
    if le in POLAR_ELEMENTS and pe in POLAR_ELEMENTS and d <= HB_MAX:
        return "hydrogen_bond"
    if ("S" in (le, pe)) and (le in POLAR_ELEMENTS or pe in POLAR_ELEMENTS) and d <= HB_S_MAX:
        return "hydrogen_bond"

    # Hydrophobic carbon contact.
    if le == "C" and pe == "C" and HYDRO_MIN <= d <= HYDRO_MAX:
        return "hydrophobic"

    return None


def _aromatic_interactions(component: Component, structure: Structure):
    """Detect aromatic ring centroid .. ligand contacts (possible pi-stacking)."""
    lig_heavy = [a for a in component.atoms if a.element != "H"]
    if not lig_heavy:
        return []

    # Group protein atoms by residue to assemble ring atom sets.
    by_res: dict[tuple, list] = {}
    for a in structure.protein_atoms:
        if a.res_name in AROMATIC_RINGS:
            by_res.setdefault((a.chain, a.res_seq, a.res_name), []).append(a)

    results = []
    for (chain, res_seq, res_name), res_atoms in by_res.items():
        ring_names = set(AROMATIC_RINGS[res_name])
        ring_atoms = [a for a in res_atoms if a.name in ring_names]
        if len(ring_atoms) < 3:
            continue
        cx = sum(a.x for a in ring_atoms) / len(ring_atoms)
        cy = sum(a.y for a in ring_atoms) / len(ring_atoms)
        cz = sum(a.z for a in ring_atoms) / len(ring_atoms)

        nearest = None
        nearest_d = ARO_CENTROID_MAX
        for la in lig_heavy:
            dd = math.sqrt((la.x - cx) ** 2 + (la.y - cy) ** 2 + (la.z - cz) ** 2)
            if dd < nearest_d:
                nearest_d = dd
                nearest = la
        if nearest is not None:
            centroid_atom = type(ring_atoms[0])(
                serial=ring_atoms[0].serial, name="ring-centroid",
                res_name=res_name, chain=chain, res_seq=res_seq, icode="",
                x=cx, y=cy, z=cz, element="C", is_hetero=False,
            )
            results.append(
                {
                    "type": "aromatic",
                    "distance": round(nearest_d, 2),
                    "ligand_atom": _atom_json(nearest, include_res=False),
                    "protein_atom": _atom_json(centroid_atom, include_res=True),
                }
            )
    return results


def profile_component(structure: Structure, component: Component) -> dict:
    """Return the full atomic interaction profile for one hetero component."""
    is_metal_comp = component.kind in ("metal", "ion")
    grid = _Grid.build(structure.protein_atoms, cell=5.5)

    raw: list[dict] = []
    # For hydrophobic noise control we keep only the closest per (residue).
    hydro_best: dict[tuple, dict] = {}

    for lig in component.atoms:
        if lig.element == "H":
            continue
        for prot in grid.neighbors(lig.x, lig.y, lig.z):
            d = _dist(lig, prot)
            kind = _classify(lig, prot, d, is_metal_comp)
            if kind is None:
                continue
            record = {
                "type": kind,
                "distance": round(d, 2),
                "ligand_atom": _atom_json(lig, include_res=False),
                "protein_atom": _atom_json(prot, include_res=True),
            }
            if kind == "hydrophobic":
                rkey = (prot.chain, prot.res_seq)
                cur = hydro_best.get(rkey)
                if cur is None or d < cur["distance"]:
                    hydro_best[rkey] = record
            else:
                raw.append(record)

    raw.extend(hydro_best.values())
    raw.extend(_aromatic_interactions(component, structure))

    # Counts + per-residue summary.
    counts = {t: 0 for t in INTERACTION_TYPES}
    res_summary: dict[str, dict] = {}
    for r in raw:
        counts[r["type"]] += 1
        pa = r["protein_atom"]
        rid = pa["res_id"]
        entry = res_summary.setdefault(
            rid,
            {
                "res_id": rid,
                "res_name": pa["res_name"],
                "res_seq": pa["res_seq"],
                "chain": pa["chain"],
                "types": set(),
                "total": 0,
                "min_distance": r["distance"],
            },
        )
        entry["types"].add(r["type"])
        entry["total"] += 1
        entry["min_distance"] = min(entry["min_distance"], r["distance"])

    residues = []
    for e in res_summary.values():
        e["types"] = sorted(e["types"])
        residues.append(e)
    residues.sort(key=lambda e: (-e["total"], e["min_distance"]))

    raw.sort(key=lambda r: (r["type"], r["distance"]))

    return {
        "component": {
            "label": component.label,
            "res_name": component.res_name,
            "chain": component.chain,
            "res_seq": component.res_seq,
            "kind": component.kind,
            "atom_count": len(component.atoms),
        },
        "interactions": raw,
        "counts": counts,
        "interaction_total": len(raw),
        "contact_residues": residues,
        "contact_residue_count": len(residues),
    }
