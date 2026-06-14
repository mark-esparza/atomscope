"""Method-transparency / provenance blocks for each analysis.

The audit's single "highest-leverage product idea" is a benchmark-first mode:
every algorithm surface should expose its method family, version, key
parameters, a plain-language "what this score means / does not mean", and its
limitations. Docking and screening already emit a methods block (see
``server._methods_block``); this module adds matching blocks for pocket
detection and conservation so *every* analytical result carries provenance.

These blocks are static for a given build (they describe the algorithm, not a
particular structure), which keeps them trivially unit-testable.
"""

from __future__ import annotations

from . import __version__
from . import evolution, pockets

_RESEARCH_DISCLAIMER = (
    "Research-only. A heuristic computed from a single static structure (and, "
    "for conservation, a family alignment). Not an affinity, not a validated "
    "binding-site or functional-site call, not clinical guidance. Validate with "
    "orthogonal evidence."
)


def tool() -> str:
    return f"SnaCleX v{__version__}"


def pocket_methods() -> dict:
    """Provenance for the LIGSITE pocket-detection + druggability scoring."""
    return {
        "tool": tool(),
        "method": "LIGSITE geometric cavity detection (Hendlich et al., 1997)",
        "method_family": "geometric pocket detection",
        "parameters": {
            "psp_threshold": pockets.PSP_THRESHOLD,
            "scan_range_steps": pockets.SCAN_RANGE,
            "min_pocket_points": pockets.MIN_POCKET_POINTS,
            "fill_radius_A": pockets.FILL,
            "box_margin_A": pockets.MARGIN,
            "max_pockets_reported": 6,
        },
        "scoring": (
            "Druggability index (0-100) = 0.30·volume + 0.30·enclosure + "
            "0.30·hydrophobicity − 0.15·polarity (+0.15 offset), clamped to "
            "[0,100]. A transparent rule-of-thumb proxy for SiteMap/fpocket-style "
            "druggability — NOT a trained model."
        ),
        "interpretation": (
            "Tiers: 'pocket' = geometric cavity only; 'ligandable' = large and "
            "enclosed enough to bind a small molecule; 'druggable' = also has "
            "favourable (hydrophobic, enclosed) chemistry. A higher score means a "
            "more promising cavity in relative terms; it does NOT confirm a real "
            "binding site or predict affinity."
        ),
        "limitations": [
            "Pure geometry on a single conformation — no induced fit or "
            "cryptic-pocket dynamics.",
            "Grid-discretized: very shallow or very small cavities may be missed "
            "or merged.",
            "Druggability is a heuristic combination of volume, enclosure and "
            "residue chemistry, not a benchmarked classifier.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def evolution_methods() -> dict:
    """Provenance for the Pfam-alignment conservation analysis."""
    return {
        "tool": tool(),
        "method": "Pfam/InterPro family-alignment conservation",
        "method_family": "evolutionary conservation",
        "parameters": {
            "alignment_source": "Pfam family alignment via InterPro / EMBL-EBI",
            "max_alignment_sequences": evolution.MAX_ALIGN_SEQS,
            "conservation_metric": "per-column Shannon entropy, 1 − H/log2(20)",
            "mapping": "Needleman-Wunsch global alignment of the structure "
                       "sequence onto alignment columns",
        },
        "interpretation": (
            "Per-residue conservation runs 0 (variable) to 1 (invariant) from "
            "the family alignment column entropy; a pocket's conservation is the "
            "mean over its lining residues. High conservation suggests functional "
            "or structural importance — it does NOT by itself mean a residue is "
            "part of a ligand-binding site."
        ),
        "limitations": [
            "Depends on the protein mapping to a Pfam family with a usable "
            "alignment; unavailable otherwise.",
            "Shallow or skewed alignments give noisy per-column entropy.",
            "Conservation reflects evolutionary pressure broadly (fold, catalysis, "
            "interfaces), not binding specifically.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }
