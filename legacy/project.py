import argparse
import copy
import random
import re

# Preferred codon table (single representative codon per amino acid)
PREFERRED_CODON_TABLE = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTT",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
    "X": "NNN",
    "*": "TAA",
}

# Optional codon alternatives for biological realism.
CODON_OPTIONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "X": ["NNN"],
    "*": ["TAA", "TAG", "TGA"],
}

AMINO_ACID_PROPERTIES = {
    "A": ("Nonpolar", "No charge"),
    "C": ("Polar", "No charge"),
    "D": ("Polar", "Negative"),
    "E": ("Polar", "Negative"),
    "F": ("Nonpolar", "No charge"),
    "G": ("Nonpolar", "No charge"),
    "H": ("Polar", "Positive"),
    "I": ("Nonpolar", "No charge"),
    "K": ("Polar", "Positive"),
    "L": ("Nonpolar", "No charge"),
    "M": ("Nonpolar", "No charge"),
    "N": ("Polar", "No charge"),
    "P": ("Nonpolar", "No charge"),
    "Q": ("Polar", "No charge"),
    "R": ("Polar", "Positive"),
    "S": ("Polar", "No charge"),
    "T": ("Polar", "No charge"),
    "V": ("Nonpolar", "No charge"),
    "W": ("Nonpolar", "No charge"),
    "Y": ("Polar", "No charge"),
    "X": ("Unknown", "Unknown"),
    "*": ("Stop", "Stop"),
}

AMINO_ACID_GROUPS = {
    "A": "nonpolar",
    "V": "nonpolar",
    "L": "nonpolar",
    "I": "nonpolar",
    "M": "nonpolar",
    "F": "nonpolar",
    "W": "nonpolar",
    "P": "nonpolar",
    "G": "special",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "Y": "polar",
    "C": "polar",
    "K": "positive",
    "R": "positive",
    "H": "positive",
    "D": "negative",
    "E": "negative",
    "X": "unknown",
    "*": "stop",
}

VARIANT_PATTERN = re.compile(r"^(?:P\.)?([A-Z*])(\d+)([A-Z*])$")

PROTEIN_KNOWLEDGE_BASE = {
    "OPTN": {
        "protein_name": "Optineurin",
        "gene_name": "OPTN",
        "uniprot_id": "Q96CV9",
        "length": 577,
        "functions": [
            "Autophagy adaptor involved in selective cargo clearance",
            "Regulation of innate immune and stress signaling",
            "Vesicle trafficking support",
        ],
        "known_diseases": ["glaucoma", "ALS", "Paget disease"],
        "major_mechanisms": [
            "autophagy dysfunction",
            "ubiquitin-binding disruption",
            "NF-kB signaling alteration",
            "mitophagy impairment",
            "protein-interaction disruption",
        ],
        "pathways": [
            "autophagy",
            "mitophagy",
            "ubiquitin-mediated signaling",
            "NF-kB signaling",
            "vesicle trafficking",
        ],
        "expression": ["retina", "brain", "immune cells"],
        "domains": [
            {"name": "Coiled-coil region 1", "start": 143, "end": 164},
            {"name": "LC3-interacting region", "start": 169, "end": 209},
            {"name": "Coiled-coil region 2", "start": 289, "end": 412},
            {"name": "UBAN ubiquitin-binding region", "start": 470, "end": 510},
            {"name": "Zinc-finger region", "start": 547, "end": 576},
        ],
        "interface_regions": [
            {
                "partner": "TBK1",
                "start": 34,
                "end": 122,
                "evidence": "experimental",
                "pathway": "autophagy and innate immunity",
                "disease_relevance": ["glaucoma", "ALS"],
            },
            {
                "partner": "LC3/GABARAP",
                "start": 169,
                "end": 209,
                "evidence": "experimental",
                "pathway": "autophagosome formation",
                "disease_relevance": ["glaucoma", "ALS"],
            },
            {
                "partner": "Ubiquitin chains",
                "start": 470,
                "end": 510,
                "evidence": "experimental",
                "pathway": "ubiquitin signaling",
                "disease_relevance": ["glaucoma", "ALS"],
            },
            {
                "partner": "RIPK1",
                "start": 360,
                "end": 425,
                "evidence": "predicted",
                "pathway": "inflammatory signaling",
                "disease_relevance": ["glaucoma"],
            },
            {
                "partner": "MYO6",
                "start": 404,
                "end": 445,
                "evidence": "curated",
                "pathway": "vesicle trafficking",
                "disease_relevance": ["retinal stress"],
            },
        ],
        "ptm_sites": [50, 177, 193, 249, 476, 513],
        "known_variants": [
            {
                "variant": "p.E50K",
                "position": 50,
                "classification": "pathogenic",
                "diseases": ["glaucoma"],
                "frequency": "rare",
            },
            {
                "variant": "p.M98K",
                "position": 98,
                "classification": "risk-associated",
                "diseases": ["glaucoma"],
                "frequency": "low",
            },
            {
                "variant": "p.H486R",
                "position": 486,
                "classification": "likely pathogenic",
                "diseases": ["ALS"],
                "frequency": "rare",
            },
            {
                "variant": "p.E478G",
                "position": 478,
                "classification": "pathogenic",
                "diseases": ["ALS"],
                "frequency": "rare",
            },
        ],
    },
    "TP53": {
        "protein_name": "Cellular tumor antigen p53",
        "gene_name": "TP53",
        "uniprot_id": "P04637",
        "length": 393,
        "functions": [
            "DNA damage response transcription factor",
            "Cell-cycle arrest and apoptosis regulation",
            "Tumor suppressor complex formation",
        ],
        "known_diseases": ["cancer", "Li-Fraumeni syndrome"],
        "major_mechanisms": [
            "DNA-binding domain disruption",
            "loss of tumor suppressor activity",
            "dominant-negative missense effects",
            "tetramerization defects",
        ],
        "pathways": [
            "DNA damage response",
            "cell-cycle checkpoint",
            "apoptosis",
            "p53 signaling",
        ],
        "expression": ["broad tissue expression"],
        "domains": [
            {"name": "Transactivation domain", "start": 1, "end": 61},
            {"name": "Proline-rich region", "start": 62, "end": 93},
            {"name": "DNA-binding domain", "start": 94, "end": 292},
            {"name": "Tetramerization domain", "start": 325, "end": 356},
        ],
        "interface_regions": [
            {
                "partner": "DNA response elements",
                "start": 117,
                "end": 286,
                "evidence": "experimental",
                "pathway": "transcriptional control",
                "disease_relevance": ["cancer"],
            },
            {
                "partner": "MDM2",
                "start": 13,
                "end": 29,
                "evidence": "experimental",
                "pathway": "protein turnover",
                "disease_relevance": ["cancer"],
            },
            {
                "partner": "p53 monomers",
                "start": 325,
                "end": 356,
                "evidence": "experimental",
                "pathway": "tetramer assembly",
                "disease_relevance": ["cancer"],
            },
        ],
        "ptm_sites": [15, 20, 46, 118, 149, 315, 392],
        "known_variants": [
            {
                "variant": "p.R175H",
                "position": 175,
                "classification": "pathogenic",
                "diseases": ["cancer"],
                "frequency": "rare",
            },
            {
                "variant": "p.R248Q",
                "position": 248,
                "classification": "pathogenic",
                "diseases": ["cancer"],
                "frequency": "rare",
            },
            {
                "variant": "p.R273H",
                "position": 273,
                "classification": "pathogenic",
                "diseases": ["cancer"],
                "frequency": "rare",
            },
            {
                "variant": "p.P72R",
                "position": 72,
                "classification": "population variant",
                "diseases": ["none"],
                "frequency": "common",
            },
        ],
    },
    "MYOC": {
        "protein_name": "Myocilin",
        "gene_name": "MYOC",
        "uniprot_id": "Q99972",
        "length": 504,
        "functions": [
            "Extracellular matrix associated protein in trabecular meshwork",
            "Stress-response involvement in ocular tissue",
        ],
        "known_diseases": ["primary open-angle glaucoma"],
        "major_mechanisms": [
            "protein misfolding and aggregation",
            "endoplasmic-reticulum stress",
            "secretory trafficking disruption",
        ],
        "pathways": [
            "protein quality control",
            "ER stress response",
            "trabecular meshwork homeostasis",
        ],
        "expression": ["eye", "trabecular meshwork"],
        "domains": [
            {"name": "N-terminal coiled-coil", "start": 117, "end": 185},
            {"name": "Olfactomedin domain", "start": 244, "end": 501},
        ],
        "interface_regions": [
            {
                "partner": "ECM-associated partners",
                "start": 260,
                "end": 500,
                "evidence": "predicted",
                "pathway": "extracellular protein interactions",
                "disease_relevance": ["glaucoma"],
            }
        ],
        "ptm_sites": [112, 236, 334, 408],
        "known_variants": [
            {
                "variant": "p.Q368*",
                "position": 368,
                "classification": "pathogenic",
                "diseases": ["glaucoma"],
                "frequency": "rare",
            },
            {
                "variant": "p.P370L",
                "position": 370,
                "classification": "pathogenic",
                "diseases": ["glaucoma"],
                "frequency": "rare",
            },
            {
                "variant": "p.Y437H",
                "position": 437,
                "classification": "pathogenic",
                "diseases": ["glaucoma"],
                "frequency": "rare",
            },
        ],
    },
    "BRCA1": {
        "protein_name": "Breast cancer type 1 susceptibility protein",
        "gene_name": "BRCA1",
        "uniprot_id": "P38398",
        "length": 1863,
        "functions": [
            "DNA double-strand break repair",
            "Genome integrity maintenance",
            "Homologous recombination support",
        ],
        "known_diseases": ["breast cancer", "ovarian cancer"],
        "major_mechanisms": [
            "DNA repair complex disruption",
            "loss of homologous recombination function",
            "protein interaction network instability",
        ],
        "pathways": [
            "homologous recombination",
            "DNA damage checkpoint",
            "ubiquitin-mediated DNA repair",
        ],
        "expression": ["broad tissue expression"],
        "domains": [
            {"name": "RING domain", "start": 24, "end": 64},
            {"name": "Coiled-coil region", "start": 1364, "end": 1437},
            {"name": "BRCT repeat 1", "start": 1646, "end": 1736},
            {"name": "BRCT repeat 2", "start": 1760, "end": 1855},
        ],
        "interface_regions": [
            {
                "partner": "BARD1",
                "start": 24,
                "end": 64,
                "evidence": "experimental",
                "pathway": "E3 ligase function",
                "disease_relevance": ["breast cancer", "ovarian cancer"],
            },
            {
                "partner": "PALB2",
                "start": 1364,
                "end": 1437,
                "evidence": "experimental",
                "pathway": "homologous recombination",
                "disease_relevance": ["breast cancer", "ovarian cancer"],
            },
            {
                "partner": "Abraxas complex",
                "start": 1646,
                "end": 1855,
                "evidence": "curated",
                "pathway": "DNA damage signaling",
                "disease_relevance": ["breast cancer"],
            },
        ],
        "ptm_sites": [42, 308, 988, 1423, 1524, 1775],
        "known_variants": [
            {
                "variant": "p.C61G",
                "position": 61,
                "classification": "pathogenic",
                "diseases": ["breast cancer", "ovarian cancer"],
                "frequency": "rare",
            },
            {
                "variant": "p.M1775R",
                "position": 1775,
                "classification": "pathogenic",
                "diseases": ["breast cancer", "ovarian cancer"],
                "frequency": "rare",
            },
            {
                "variant": "p.I157T",
                "position": 157,
                "classification": "risk-associated",
                "diseases": ["breast cancer"],
                "frequency": "low",
            },
        ],
    },
    "APOE": {
        "protein_name": "Apolipoprotein E",
        "gene_name": "APOE",
        "uniprot_id": "P02649",
        "length": 317,
        "functions": [
            "Lipid transport and receptor binding",
            "Neuronal lipid homeostasis",
            "Inflammatory response modulation",
        ],
        "known_diseases": ["Alzheimer disease", "cardiovascular risk"],
        "major_mechanisms": [
            "lipid transport balance shift",
            "amyloid handling differences",
            "neuronal resilience modulation",
        ],
        "pathways": [
            "lipoprotein metabolism",
            "neuronal lipid transport",
            "amyloid-related pathways",
        ],
        "expression": ["liver", "brain"],
        "domains": [
            {"name": "Receptor binding region", "start": 136, "end": 150},
            {"name": "Lipid-binding C-terminal region", "start": 244, "end": 272},
        ],
        "interface_regions": [
            {
                "partner": "LDL receptor family",
                "start": 136,
                "end": 150,
                "evidence": "experimental",
                "pathway": "lipoprotein uptake",
                "disease_relevance": ["Alzheimer disease", "cardiovascular risk"],
            },
            {
                "partner": "Lipoprotein particles",
                "start": 244,
                "end": 272,
                "evidence": "experimental",
                "pathway": "lipid transport",
                "disease_relevance": ["Alzheimer disease"],
            },
        ],
        "ptm_sites": [130, 151, 194, 263],
        "known_variants": [
            {
                "variant": "p.C130R",
                "position": 130,
                "classification": "risk-associated",
                "diseases": ["Alzheimer disease"],
                "frequency": "common",
            },
            {
                "variant": "p.C176R",
                "position": 176,
                "classification": "risk-associated",
                "diseases": ["Alzheimer disease"],
                "frequency": "common",
            },
        ],
    },
}

VALID_AMINO_ACIDS = set(PREFERRED_CODON_TABLE.keys())


def normalize_protein_sequence(protein_seq):
    """Normalize user input by stripping spaces and uppercasing."""
    if protein_seq is None:
        return ""
    return "".join(protein_seq.split()).upper()


def normalize_protein_identifier(protein_id):
    """Normalize protein identifiers such as gene symbols."""
    if protein_id is None:
        return ""
    return "".join(str(protein_id).split()).upper()


def validate_protein_sequence(protein_seq):
    """Return (True, None) if valid, else (False, [invalid_chars])."""
    invalid_chars = sorted({aa for aa in protein_seq if aa not in VALID_AMINO_ACIDS})
    if invalid_chars:
        return False, invalid_chars
    return True, None


def protein_to_dna(protein_seq, strategy="preferred"):
    """Convert protein sequence to DNA using either preferred or random codons."""
    protein_seq = normalize_protein_sequence(protein_seq)
    is_valid, invalid_chars = validate_protein_sequence(protein_seq)
    if not is_valid:
        invalid = ", ".join(invalid_chars)
        raise ValueError(f"Invalid amino acid character(s): {invalid}")

    if strategy == "preferred":
        return "".join(PREFERRED_CODON_TABLE[aa] for aa in protein_seq)

    if strategy == "random":
        return "".join(random.choice(CODON_OPTIONS[aa]) for aa in protein_seq)

    raise ValueError("Invalid strategy. Use 'preferred' or 'random'.")


def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")


def count_nucleotides(sequence):
    sequence = sequence.upper()
    return {nucleotide: sequence.count(nucleotide) for nucleotide in "ACGTU"}


def calculate_polarity_and_charge_of_protein(protein_seq):
    protein_seq = normalize_protein_sequence(protein_seq)
    is_valid, invalid_chars = validate_protein_sequence(protein_seq)
    if not is_valid:
        invalid = ", ".join(invalid_chars)
        raise ValueError(f"Invalid amino acid character(s): {invalid}")

    polarity_count = {"Nonpolar": 0, "Polar": 0, "Unknown": 0, "Stop": 0}
    charge_count = {
        "Positive": 0,
        "Negative": 0,
        "No charge": 0,
        "Unknown": 0,
        "Stop": 0,
    }

    for aa in protein_seq:
        polarity, charge = AMINO_ACID_PROPERTIES[aa]
        polarity_count[polarity] += 1
        charge_count[charge] += 1

    return polarity_count, charge_count


def parse_variant_label(variant_label):
    """Parse HGVS-like protein variant labels such as p.E50K."""
    if variant_label is None:
        raise ValueError("Variant label is required.")

    normalized = str(variant_label).strip().upper()
    if not normalized:
        raise ValueError("Variant label is required.")

    match = VARIANT_PATTERN.match(normalized)
    if not match:
        raise ValueError(
            f"Invalid variant format: {variant_label}. Use forms like p.E50K or E50K."
        )

    ref_aa, pos_text, alt_aa = match.groups()
    return {
        "label": f"p.{ref_aa}{pos_text}{alt_aa}",
        "reference": ref_aa,
        "position": int(pos_text),
        "alternate": alt_aa,
    }


def parse_variant_list(raw_variants):
    """Parse a variant list from a comma/semicolon separated string or list."""
    if raw_variants is None:
        return []

    if isinstance(raw_variants, str):
        items = re.split(r"[,;]", raw_variants)
    else:
        items = list(raw_variants)

    variants = []
    seen = set()
    for item in items:
        cleaned = str(item).strip()
        if not cleaned:
            continue
        key = cleaned.upper()
        if key in seen:
            continue
        seen.add(key)
        variants.append(cleaned)

    return variants


def get_protein_record(protein_id):
    """Return curated protein record, or a minimal heuristic fallback."""
    key = normalize_protein_identifier(protein_id)
    if key in PROTEIN_KNOWLEDGE_BASE:
        return copy.deepcopy(PROTEIN_KNOWLEDGE_BASE[key])

    return {
        "protein_name": key or "Unknown protein",
        "gene_name": key or "N/A",
        "uniprot_id": "N/A",
        "length": None,
        "functions": ["No curated functions available in the local dataset"],
        "known_diseases": [],
        "major_mechanisms": [
            "Insufficient curated data: use external databases for stronger interpretation"
        ],
        "pathways": [],
        "expression": [],
        "domains": [],
        "interface_regions": [],
        "ptm_sites": [],
        "known_variants": [],
        "data_source": "heuristic_fallback",
    }


def _find_matching_known_variant(protein_record, variant_label, position):
    normalized_label = variant_label.upper()
    for known_variant in protein_record.get("known_variants", []):
        if known_variant.get("variant", "").upper() == normalized_label:
            return known_variant

    for known_variant in protein_record.get("known_variants", []):
        if known_variant.get("position") == position:
            return known_variant

    return None


def _regions_containing_position(position, regions):
    return [region for region in regions if region["start"] <= position <= region["end"]]


def _minimum_distance_to_positions(position, other_positions):
    if not other_positions:
        return None
    return min(abs(position - other_pos) for other_pos in other_positions)


def _amino_acid_group(amino_acid):
    return AMINO_ACID_GROUPS.get(amino_acid, "unknown")


def _substitution_disruption_score(reference, alternate):
    if reference == alternate:
        return 0.0

    if alternate == "*":
        return 2.0

    if reference == "*":
        return 1.0

    score = 0.7
    if _amino_acid_group(reference) != _amino_acid_group(alternate):
        score += 0.6

    if reference in {"G", "P", "C"} or alternate in {"G", "P", "C"}:
        score += 0.3

    return min(score, 2.0)


def _classification_weight(classification):
    normalized = str(classification).lower()
    if normalized == "pathogenic":
        return 2.0
    if normalized == "likely pathogenic":
        return 1.7
    if normalized == "risk-associated":
        return 1.4
    if normalized in {"vus", "unknown"}:
        return 1.0
    if normalized in {"benign", "likely benign"}:
        return 0.3
    if normalized == "population variant":
        return 0.4
    return 0.8


def _score_pathway_relevance(protein_record, disease_focus):
    if not disease_focus:
        return 0.8

    focus_lower = disease_focus.lower()
    known_diseases = [d.lower() for d in protein_record.get("known_diseases", [])]

    if any(focus_lower in disease for disease in known_diseases):
        return 1.6

    pathway_text = " ".join(protein_record.get("pathways", [])).lower()
    mechanism_text = " ".join(protein_record.get("major_mechanisms", [])).lower()
    combined = f"{pathway_text} {mechanism_text}"

    keyword_map = {
        "glaucoma": ["autophagy", "retina", "stress"],
        "als": ["autophagy", "mitophagy", "neur"],
        "cancer": ["dna", "cell-cycle", "tumor", "repair"],
        "alzheimer": ["amyloid", "lipid", "neur"],
        "autoimmune": ["immune", "antigen", "nf-kb"],
    }

    for disease_keyword, keywords in keyword_map.items():
        if disease_keyword in focus_lower and any(keyword in combined for keyword in keywords):
            return 1.2

    return 0.8


def _score_network_centrality(interaction_regions):
    degree = len(interaction_regions)
    if degree >= 6:
        return 1.5
    if degree >= 4:
        return 1.2
    if degree >= 2:
        return 0.9
    if degree == 1:
        return 0.6
    return 0.2


def _score_rarity(frequency_text, classification_text):
    frequency = str(frequency_text).lower()
    classification = str(classification_text).lower()

    if frequency in {"ultra-rare", "very rare"}:
        return 1.6
    if frequency == "rare":
        return 1.4
    if frequency == "low":
        return 1.1
    if frequency == "common":
        return 0.3
    if classification in {"benign", "likely benign", "population variant"}:
        return 0.3
    return 0.8


def pmds_band(score):
    if score >= 9.0:
        return "High"
    if score >= 7.0:
        return "Moderate-High"
    if score >= 4.5:
        return "Moderate"
    return "Low"


def score_variant_for_mechanism(protein_record, variant_label, disease_focus=None):
    """Compute a PMDS-like score for one variant against one protein context."""
    parsed = parse_variant_label(variant_label)
    position = parsed["position"]

    protein_length = protein_record.get("length")
    if protein_length and position > protein_length:
        raise ValueError(
            f"Variant {parsed['label']} position exceeds "
            f"{protein_record['gene_name']} length ({protein_length})."
        )

    domains = _regions_containing_position(position, protein_record.get("domains", []))
    interfaces = _regions_containing_position(
        position, protein_record.get("interface_regions", [])
    )

    known_variant = _find_matching_known_variant(protein_record, parsed["label"], position)
    known_classification = known_variant.get("classification") if known_variant else "unknown"
    known_frequency = known_variant.get("frequency") if known_variant else "unknown"
    known_diseases = known_variant.get("diseases", []) if known_variant else []

    structural_score = _substitution_disruption_score(parsed["reference"], parsed["alternate"])
    if domains:
        structural_score += 0.4
    if interfaces:
        structural_score += 0.5

    ptm_sites = protein_record.get("ptm_sites", [])
    ptm_distance = _minimum_distance_to_positions(position, ptm_sites)
    if ptm_distance is not None and ptm_distance <= 2:
        structural_score += 0.2

    structural_score = min(structural_score, 2.0)

    if known_variant:
        disease_proximity = _classification_weight(known_classification)
        if disease_focus and any(disease_focus.lower() in d.lower() for d in known_diseases):
            disease_proximity = min(disease_proximity + 0.2, 2.0)
    else:
        pathogenic_positions = [
            variant["position"]
            for variant in protein_record.get("known_variants", [])
            if variant.get("classification", "").lower() in {"pathogenic", "likely pathogenic"}
        ]
        nearest_pathogenic = _minimum_distance_to_positions(position, pathogenic_positions)
        if nearest_pathogenic is None:
            disease_proximity = 0.6
        elif nearest_pathogenic <= 5:
            disease_proximity = 1.4
        elif nearest_pathogenic <= 12:
            disease_proximity = 1.0
        else:
            disease_proximity = 0.7

    interface_score = 0.2
    if interfaces:
        interface_score = 1.4
        if known_variant and known_classification.lower() in {
            "pathogenic",
            "likely pathogenic",
            "risk-associated",
        }:
            interface_score = 1.8
    else:
        distances_to_interface = [
            _minimum_distance_to_positions(position, range(region["start"], region["end"] + 1))
            for region in protein_record.get("interface_regions", [])
        ]
        distances_to_interface = [distance for distance in distances_to_interface if distance]
        if distances_to_interface and min(distances_to_interface) <= 8:
            interface_score = 0.9

    pathway_relevance = _score_pathway_relevance(protein_record, disease_focus)

    if domains:
        conservation_score = 1.3
    elif protein_record.get("domains"):
        conservation_score = 0.7
    else:
        conservation_score = 0.5

    rarity_score = _score_rarity(known_frequency, known_classification)

    if ptm_distance is None:
        ptm_score = 0.2
    elif ptm_distance <= 2:
        ptm_score = 1.2
    elif ptm_distance <= 5:
        ptm_score = 0.8
    else:
        ptm_score = 0.3

    centrality_score = _score_network_centrality(protein_record.get("interface_regions", []))

    components = {
        "structural_disruption": round(structural_score, 2),
        "disease_variant_proximity": round(disease_proximity, 2),
        "interaction_interface": round(interface_score, 2),
        "pathway_relevance": round(pathway_relevance, 2),
        "conservation": round(conservation_score, 2),
        "population_rarity": round(rarity_score, 2),
        "ptm_proximity": round(ptm_score, 2),
        "network_centrality": round(centrality_score, 2),
    }

    total_score = round(sum(components.values()), 2)
    sorted_component_names = sorted(components, key=components.get, reverse=True)

    return {
        "variant": parsed["label"],
        "position": position,
        "classification": known_classification,
        "frequency": known_frequency,
        "associated_diseases": known_diseases,
        "domains": [domain["name"] for domain in domains],
        "interfaces": [region["partner"] for region in interfaces],
        "pmds_score": total_score,
        "pmds_band": pmds_band(total_score),
        "components": components,
        "top_drivers": sorted_component_names[:3],
    }


def rank_interactions(protein_record, disease_focus=None):
    """Rank interaction partners by support and disease/pathway relevance."""
    focus_lower = disease_focus.lower() if disease_focus else None
    ranked = []

    evidence_weight = {
        "experimental": 2.0,
        "curated": 1.7,
        "predicted": 1.2,
        "text-mined": 0.8,
    }

    for interaction in protein_record.get("interface_regions", []):
        score = evidence_weight.get(interaction.get("evidence", "predicted"), 1.0)

        if focus_lower:
            disease_matches = [
                disease
                for disease in interaction.get("disease_relevance", [])
                if focus_lower in disease.lower()
            ]
            if disease_matches:
                score += 1.0

            pathway_text = interaction.get("pathway", "").lower()
            if focus_lower in pathway_text:
                score += 0.4

        entry = {
            "partner": interaction["partner"],
            "evidence": interaction.get("evidence", "predicted"),
            "pathway": interaction.get("pathway", ""),
            "disease_relevance": interaction.get("disease_relevance", []),
            "rank_score": round(score, 2),
        }
        ranked.append(entry)

    ranked.sort(key=lambda item: item["rank_score"], reverse=True)
    return ranked


def generate_hypotheses(protein_record, scored_variants, disease_focus=None):
    """Generate concise, testable research hypotheses."""
    hypotheses = []
    high_impact = [variant for variant in scored_variants if variant["pmds_score"] >= 7.0]
    interface_hits = [variant for variant in high_impact if variant["interfaces"]]

    if interface_hits:
        hypotheses.append(
            "High-impact variants are enriched near protein interaction surfaces, "
            "suggesting altered partner-binding may drive disease-relevant signaling changes."
        )

    domain_hit_count = sum(1 for variant in high_impact if variant["domains"])
    if high_impact and domain_hit_count >= max(1, len(high_impact) // 2):
        hypotheses.append(
            "Most high-impact variants map to structured domains, which supports a mechanism "
            "involving local structural destabilization of functionally constrained regions."
        )

    if len(protein_record.get("interface_regions", [])) >= 3:
        hypotheses.append(
            "The protein appears to act as a network connector; perturbations may affect "
            "multiple pathways through interaction rewiring instead of complete loss of function."
        )

    if disease_focus:
        hypotheses.append(
            f"In the context of {disease_focus}, variants with elevated PMDS should be "
            "prioritized for follow-up in pathway-specific assays and interaction readouts."
        )

    if not hypotheses:
        hypotheses.append(
            "Current data are limited; prioritize generating structure-aware variant annotations "
            "and experimentally validated interaction data for stronger disease-mechanism hypotheses."
        )

    return hypotheses[:3]


def analyze_struct_interact(protein_id, disease_focus=None, variants=None):
    """Run a StructInteract-style analysis for one protein."""
    protein_record = get_protein_record(protein_id)
    variant_labels = parse_variant_list(variants)

    if not variant_labels and protein_record.get("known_variants"):
        prioritized = []
        for known_variant in protein_record["known_variants"]:
            if known_variant.get("classification", "").lower() in {
                "pathogenic",
                "likely pathogenic",
                "risk-associated",
            }:
                prioritized.append(known_variant["variant"])
            if len(prioritized) == 3:
                break
        variant_labels = prioritized

    scored_variants = [
        score_variant_for_mechanism(protein_record, variant, disease_focus=disease_focus)
        for variant in variant_labels
    ]

    mechanisms = list(protein_record.get("major_mechanisms", []))
    interactions = rank_interactions(protein_record, disease_focus=disease_focus)
    hypotheses = generate_hypotheses(
        protein_record, scored_variants, disease_focus=disease_focus
    )

    limitations = [
        "This MVP uses a local curated dataset and heuristic scoring; "
        "it is research-only and not clinical guidance."
    ]

    if protein_record.get("data_source") == "heuristic_fallback":
        limitations.append(
            "No curated local record found for this protein; integrate "
            "UniProt/ClinVar/STRING data for deeper interpretation."
        )

    if not scored_variants:
        limitations.append(
            "No variants were scored. Provide --variants to compute PMDS-style assessments."
        )

    return {
        "protein": {
            "protein_name": protein_record.get("protein_name"),
            "gene_name": protein_record.get("gene_name"),
            "uniprot_id": protein_record.get("uniprot_id"),
            "length": protein_record.get("length"),
            "known_diseases": protein_record.get("known_diseases", []),
            "functions": protein_record.get("functions", []),
            "pathways": protein_record.get("pathways", []),
            "expression": protein_record.get("expression", []),
        },
        "disease_focus": disease_focus,
        "mechanisms": mechanisms,
        "interactions": interactions,
        "variant_assessments": scored_variants,
        "hypotheses": hypotheses,
        "limitations": limitations,
    }


def format_struct_interact_report(result):
    """Format a user-facing StructInteract report."""
    protein = result["protein"]
    disease_focus = result.get("disease_focus") or "Not provided"

    lines = [
        "StructInteract Report",
        "====================",
        f"Protein: {protein['gene_name']} ({protein['protein_name']})",
        f"UniProt: {protein['uniprot_id']}",
        f"Length: {protein['length'] if protein['length'] is not None else 'Unknown'} aa",
        f"Disease focus: {disease_focus}",
        "",
        "Known disease links:",
    ]

    if protein["known_diseases"]:
        for disease in protein["known_diseases"]:
            lines.append(f"- {disease}")
    else:
        lines.append("- No curated disease links in local dataset")

    lines.append("")
    lines.append("Mechanism map:")
    if result["mechanisms"]:
        for mechanism in result["mechanisms"]:
            lines.append(f"- {mechanism}")
    else:
        lines.append("- No mechanism annotations available")

    lines.append("")
    lines.append("Top interaction partners:")
    if result["interactions"]:
        for interaction in result["interactions"][:5]:
            lines.append(
                f"- {interaction['partner']} | evidence: {interaction['evidence']} | "
                f"pathway: {interaction['pathway']} | rank: {interaction['rank_score']}"
            )
    else:
        lines.append("- No curated interactions available")

    lines.append("")
    lines.append("Variant PMDS-style assessment:")
    if result["variant_assessments"]:
        for item in result["variant_assessments"]:
            lines.append(
                f"- {item['variant']} -> {item['pmds_band']} ({item['pmds_score']}) | "
                f"class: {item['classification']} | top drivers: {', '.join(item['top_drivers'])}"
            )
    else:
        lines.append("- No variants scored")

    lines.append("")
    lines.append("Top hypotheses:")
    for index, hypothesis in enumerate(result["hypotheses"], start=1):
        lines.append(f"{index}. {hypothesis}")

    lines.append("")
    lines.append("Limitations:")
    for limitation in result["limitations"]:
        lines.append(f"- {limitation}")

    return "\n".join(lines)


def print_all_analyses(protein_sequence, strategy):
    dna_seq = protein_to_dna(protein_sequence, strategy=strategy)
    rna_seq = dna_to_rna(dna_seq)
    dna_count = count_nucleotides(dna_seq)
    rna_count = count_nucleotides(rna_seq)
    polarity_count, charge_count = calculate_polarity_and_charge_of_protein(protein_sequence)

    print(f"DNA sequence ({strategy} codons): {dna_seq}")
    print(f"RNA sequence: {rna_seq}")
    print(f"Nucleotide counts (DNA): {dna_count}")
    print(f"Nucleotide counts (RNA): {rna_count}")
    print(f"Polarity count: {polarity_count}")
    print(f"Charge count: {charge_count}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze protein sequences or run a StructInteract-style "
            "disease/interactome analysis."
        )
    )
    parser.add_argument(
        "protein_sequence",
        nargs="?",
        type=str,
        help="Protein sequence for menu-based sequence analysis",
    )
    parser.add_argument(
        "--codon-strategy",
        choices=["preferred", "random"],
        default="preferred",
        help="Codon selection mode for protein->DNA conversion",
    )
    parser.add_argument(
        "--struct-interact-protein",
        type=str,
        help="Gene/protein symbol to analyze with StructInteract (example: OPTN, TP53)",
    )
    parser.add_argument(
        "--disease-focus",
        type=str,
        default=None,
        help="Optional disease focus to weight mechanisms and interactions",
    )
    parser.add_argument(
        "--variants",
        type=str,
        default=None,
        help="Comma-separated protein variants (example: p.E50K,p.R214W)",
    )

    args = parser.parse_args()

    if args.struct_interact_protein:
        analysis = analyze_struct_interact(
            args.struct_interact_protein,
            disease_focus=args.disease_focus,
            variants=args.variants,
        )
        print(format_struct_interact_report(analysis))
        return

    if args.protein_sequence is None:
        parser.error(
            "protein_sequence is required unless --struct-interact-protein is provided."
        )

    protein_sequence = normalize_protein_sequence(args.protein_sequence)

    is_valid, invalid_chars = validate_protein_sequence(protein_sequence)
    if not is_valid:
        invalid = ", ".join(invalid_chars)
        print(f"Invalid protein sequence. Unknown character(s): {invalid}")
        return

    while True:
        print("\nMenu:")
        print("1. Convert to DNA")
        print("2. Convert to RNA")
        print("3. Count Nucleotides")
        print("4. Calculate Polarity and Charge of Protein")
        print("5. Analyze All")
        print("6. Exit")

        choice = input("Enter your choice: ").strip()

        try:
            if choice == "1":
                dna_seq = protein_to_dna(protein_sequence, strategy=args.codon_strategy)
                print(f"DNA sequence ({args.codon_strategy} codons): {dna_seq}")
            elif choice == "2":
                dna_seq = protein_to_dna(protein_sequence, strategy=args.codon_strategy)
                rna_seq = dna_to_rna(dna_seq)
                print(f"RNA sequence: {rna_seq}")
            elif choice == "3":
                dna_seq = protein_to_dna(protein_sequence, strategy=args.codon_strategy)
                dna_count = count_nucleotides(dna_seq)
                print(f"Nucleotide counts (DNA): {dna_count}")
                rna_seq = dna_to_rna(dna_seq)
                rna_count = count_nucleotides(rna_seq)
                print(f"Nucleotide counts (RNA): {rna_count}")
            elif choice == "4":
                polarity_count, charge_count = calculate_polarity_and_charge_of_protein(
                    protein_sequence
                )
                print(f"Polarity count: {polarity_count}")
                print(f"Charge count: {charge_count}")
            elif choice == "5":
                print_all_analyses(protein_sequence, args.codon_strategy)
            elif choice == "6":
                print("Exiting...")
                break
            else:
                print("Invalid choice. Please enter a number between 1 and 6.")
        except ValueError as err:
            print(err)


if __name__ == "__main__":
    main()
