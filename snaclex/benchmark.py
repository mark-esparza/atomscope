"""Redocking benchmark harness (pure stdlib, no external engine).

Self-docks each crystallographic ligand back into its own (apo) receptor and
measures how close the predicted pose lands to the experimental one — the
standard first-order check on a docking method, and the audit's clearest
"measure before you change it" recommendation for the docking path.

It runs on local PDB files so it works offline / in CI; pass 4-character PDB IDs
instead to fetch from RCSB. The *same* harness scales to the public benchmarks
the audit names — point it at a directory of PoseBusters / CrossDocked / PDBbind
structures:

    python -m snaclex.benchmark structures/*.pdb --out benchmark_results.json

RMSD here is the nearest-atom proxy from ``docking.rmsd_to_reference`` (no
symmetry correction), so read sub-2 Å as "pose recovered", not as a rigorous
superposition RMSD. Optional Vina/GNINA backends (Phase 4) can be benchmarked
head-to-head through the same ``summarize`` once those engines are installed.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import statistics
import sys

from . import docking, interactions, pdbparse, pockets, rcsb


def _ligand_atoms(component):
    return [
        {"element": a.element, "x": a.x, "y": a.y, "z": a.z}
        for a in component.atoms
        if a.element != "H"
    ]


def redock_case(structure, component, *, seeds=220, mc_steps=40, seed=0):
    """Self-dock one ligand back into its receptor; return RMSD vs the crystal pose."""
    heavy = [a for a in component.atoms if a.element != "H"]
    if len(heavy) < 3:
        return {"label": component.label, "rmsd": None, "skipped": "too few heavy atoms"}
    center = docking.component_center(component)
    grid = docking.build_grid(structure, center)
    pose = docking.dock_with_grid(
        grid, _ligand_atoms(component), center,
        seeds=seeds, mc_steps=mc_steps, seed=seed,
    )
    return {
        "label": component.label,
        "res_name": component.res_name,
        "n_heavy_atoms": len(heavy),
        "score": pose["score"],
        "rmsd": docking.rmsd_to_reference(pose, component),
    }


def pocket_recovery(pockets_found, center, tol_A=6.0):
    """Did pocket detection recover the known site? Nearest pocket to `center`."""
    best = None
    for p in pockets_found:
        d = math.dist(p["center"], center)
        if best is None or d < best[0]:
            best = (d, p["index"])
    if best is None:
        return {"recovered": False, "distance_A": None, "rank": None, "n_pockets": 0}
    dist, idx = best
    return {
        "recovered": dist <= tol_A,
        "distance_A": round(dist, 2),
        "rank": idx + 1,
        "n_pockets": len(pockets_found),
    }


def _has_clash(pose, structure, tol_A=2.0):
    """True if any pose heavy atom overlaps a protein heavy atom (< tol_A)."""
    cell = 5.0
    buckets = {}
    for a in structure.protein_atoms:
        if a.element == "H":
            continue
        k = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
        buckets.setdefault(k, []).append(a)
    tol2 = tol_A * tol_A
    for (x, y, z) in pose["pose_coords"]:
        cx, cy, cz = int(x // cell), int(y // cell), int(z // cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for a in buckets.get((cx + dx, cy + dy, cz + dz), ()):
                        if (x - a.x) ** 2 + (y - a.y) ** 2 + (z - a.z) ** 2 < tol2:
                            return True
    return False


def benchmark_case(structure, component, *, seeds=220, mc_steps=40, seed=0,
                   pocket_tol_A=6.0, clash_tol_A=2.0, rmsd_success_A=2.0):
    """Evaluate SnaCleX against one known protein–ligand case.

    Self-docks the crystallographic ligand back into its own site and reports the
    research-credibility metrics: pocket recovery, pose RMSD, fraction of the
    experimental contact residues recovered, and a physical-plausibility (no
    severe clash) pass/fail.
    """
    center = docking.component_center(component)
    ref_heavy = [a for a in component.atoms if a.element != "H"]
    lig_atoms = [{"element": a.element, "x": a.x, "y": a.y, "z": a.z} for a in ref_heavy]

    grid = docking.build_grid(structure, center)
    pose = docking.dock_with_grid(
        grid, lig_atoms, center, seeds=seeds, mc_steps=mc_steps, seed=seed
    )
    rmsd = docking.rmsd_to_reference(pose, component)

    rec = pocket_recovery(pockets.detect_pockets(structure), center, pocket_tol_A)

    ref_prof = interactions.profile_component(structure, component)
    pose_comp = docking.pose_to_component(pose, component.res_name)
    pose_prof = interactions.profile_component(structure, pose_comp)
    ref_res = {r["res_id"] for r in ref_prof["contact_residues"]}
    pose_res = {r["res_id"] for r in pose_prof["contact_residues"]}
    recovered = sorted(ref_res & pose_res)

    clash = _has_clash(pose, structure, clash_tol_A)

    return {
        "ligand": component.label,
        "ligand_res_name": component.res_name,
        "n_heavy_atoms": len(ref_heavy),
        "pocket": rec,
        "pose_rmsd_A": rmsd,
        "rmsd_success": rmsd is not None and rmsd <= rmsd_success_A,
        "rmsd_success_threshold_A": rmsd_success_A,
        "interactions_recovered": len(recovered),
        "interactions_total": len(ref_res),
        "recovered_residues": recovered,
        "reference_residues": sorted(ref_res),
        "physical_plausibility": "fail" if clash else "pass",
        "score": pose["score"],
        "search": pose["search"],
    }


def benchmark_structure(structure, *, min_heavy=6, **kw):
    """Redock every ligand component with at least ``min_heavy`` heavy atoms."""
    cases = []
    for comp in structure.ligand_components:
        heavy = [a for a in comp.atoms if a.element != "H"]
        if len(heavy) < min_heavy:
            continue
        cases.append(redock_case(structure, comp, **kw))
    return cases


def summarize(cases, success_A=2.0):
    """Aggregate redock cases into the standard pose-recovery metrics."""
    rmsds = [c["rmsd"] for c in cases if c.get("rmsd") is not None]
    n_scored = len(rmsds)
    success = sum(1 for r in rmsds if r <= success_A)
    return {
        "n_cases": len(cases),
        "n_scored": n_scored,
        "success_threshold_A": success_A,
        "top1_success_rate": round(success / n_scored, 3) if n_scored else None,
        "median_rmsd_A": round(statistics.median(rmsds), 2) if rmsds else None,
        "mean_rmsd_A": round(statistics.fmean(rmsds), 2) if rmsds else None,
    }


# ---------------- CLI ----------------

def _load_source(source):
    if source.lower().endswith((".pdb", ".ent", ".cif", ".mmcif")):
        with open(source, encoding="utf-8") as fh:
            text = fh.read()
        return os.path.basename(source), pdbparse.parse_structure(text)
    # Otherwise treat it as a 4-character PDB ID and fetch it.
    return rcsb.normalize_pdb_id(source), pdbparse.parse_structure(rcsb.fetch_structure(source))


def run(sources, *, seeds=220, min_heavy=6):
    cases = []
    for src in sources:
        name, structure = _load_source(src)
        for c in benchmark_structure(structure, min_heavy=min_heavy, seeds=seeds):
            c["structure"] = name
            cases.append(c)
    return {
        "tool": f"SnaCleX v{__import__('snaclex').__version__}",
        "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "method": "self-docking, nearest-atom RMSD",
        "summary": summarize(cases),
        "cases": cases,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="SnaCleX redocking benchmark")
    parser.add_argument("sources", nargs="+", help="local .pdb files or 4-char PDB IDs")
    parser.add_argument("--seeds", type=int, default=220, help="MC restarts per dock")
    parser.add_argument("--min-heavy", type=int, default=6, help="min ligand heavy atoms")
    parser.add_argument("--out", help="write full results JSON to this path")
    args = parser.parse_args(argv)

    results = run(args.sources, seeds=args.seeds, min_heavy=args.min_heavy)
    for c in results["cases"]:
        rmsd = "skip" if c.get("rmsd") is None else f"{c['rmsd']:.2f} A"
        print(f"  {c['structure']:<16} {c['label']:<14} rmsd={rmsd:<10} score={c.get('score')}")
    s = results["summary"]
    print(
        f"\n{s['n_scored']}/{s['n_cases']} scored · "
        f"top-1 ≤{s['success_threshold_A']}Å: {s['top1_success_rate']} · "
        f"median {s['median_rmsd_A']} Å · mean {s['mean_rmsd_A']} Å"
    )
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)
        print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
