"""Turn an interaction profile into a plain-language summary + hypotheses."""

from __future__ import annotations

_TYPE_LABEL = {
    "metal_coordination": "metal coordination bond(s)",
    "salt_bridge": "possible salt bridge(s)",
    "hydrogen_bond": "hydrogen bond / polar contact(s)",
    "hydrophobic": "hydrophobic contact(s)",
    "aromatic": "aromatic contact(s) (possible pi-stacking)",
}


def summarize(profile: dict, metadata: dict | None = None) -> dict:
    counts = profile.get("counts", {})
    comp = profile.get("component", {})
    residues = profile.get("contact_residues", [])

    parts = []
    for t, label in _TYPE_LABEL.items():
        n = counts.get(t, 0)
        if n:
            parts.append(f"{n} {label}")
    breakdown = "; ".join(parts) if parts else "no significant heavy-atom contacts"

    top_res = ", ".join(
        f"{r['res_name']}{r['res_seq']}" for r in residues[:6]
    ) or "none"

    title = (metadata or {}).get("title") or "this protein"
    pdb_id = (metadata or {}).get("pdb_id", "")

    summary_text = (
        f"In {pdb_id}, the component {comp.get('label')} ({comp.get('kind')}) "
        f"engages {profile.get('contact_residue_count', 0)} protein residues "
        f"through {breakdown}. The most-contacted residues are: {top_res}."
    )

    hypotheses = _hypotheses(counts, comp, residues)

    return {
        "summary": summary_text,
        "binding_site_residues": top_res,
        "hypotheses": hypotheses,
        "disclaimer": (
            "Research-only. Interactions are geometric heuristics computed from "
            "experimental coordinates (no explicit hydrogens, no energy "
            "minimization). Not for clinical or diagnostic use."
        ),
    }


def _hypotheses(counts, comp, residues) -> list[str]:
    h = []
    kind = comp.get("kind")
    label = comp.get("label")

    if counts.get("metal_coordination"):
        h.append(
            f"{label} appears to be coordinated by protein side chains, suggesting "
            f"a structural or catalytic metal site; substitutions at coordinating "
            f"residues would be predicted to disrupt binding."
        )
    if counts.get("salt_bridge"):
        h.append(
            f"Electrostatic (salt-bridge) contacts dominate part of the interface, "
            f"suggesting binding is sensitive to charge — pH shifts or charge-altering "
            f"mutations near the pocket may modulate affinity."
        )
    if counts.get("hydrogen_bond", 0) >= 3:
        h.append(
            f"A hydrogen-bond network ({counts['hydrogen_bond']} polar contacts) anchors "
            f"the ligand, implying specificity is driven by directional polar interactions "
            f"rather than shape alone."
        )
    if counts.get("hydrophobic", 0) >= 4:
        h.append(
            f"Extensive hydrophobic packing ({counts['hydrophobic']} contacts) suggests a "
            f"largely apolar pocket; analogs with added hydrophobic bulk could be tested "
            f"for improved fit."
        )
    if counts.get("aromatic"):
        h.append(
            f"Aromatic (possible pi-stacking) contacts are present, suggesting ring-system "
            f"orientation contributes to binding and could be probed by aromatic-ring edits."
        )
    if residues:
        rid = ", ".join(f"{r['res_name']}{r['res_seq']}" for r in residues[:4])
        h.append(
            f"The residues {rid} form the predicted binding hot spot and are the priority "
            f"set for mutagenesis or comparison across homologs."
        )
    if not h:
        h.append(
            f"No strong directional interactions were detected for {label}; binding (if any) "
            f"may be weak, surface-level, or require a different conformational state."
        )
    return h
