import pytest

from project import (
    analyze_struct_interact,
    calculate_polarity_and_charge_of_protein,
    count_nucleotides,
    dna_to_rna,
    format_struct_interact_report,
    parse_variant_label,
    parse_variant_list,
    protein_to_dna,
)


def test_protein_to_dna_preferred():
    result = protein_to_dna("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT")
    expected = (
        "ATGACTGAATATAAACTTGTTGTTGTTGGTGCTGGTGGTGTTGGTAAATCTGCTCTTACT"
        "ATTCAACTTATTCAAAATCATTTTGTTGATGAATATGATCCTACT"
    )
    assert result == expected


def test_protein_to_dna_random_valid_output_length_and_charset():
    protein = "MTEYK"
    result = protein_to_dna(protein, strategy="random")
    assert len(result) == len(protein) * 3
    assert set(result).issubset({"A", "C", "G", "T", "N"})


def test_protein_to_dna_empty_sequence():
    assert protein_to_dna("") == ""


def test_protein_to_dna_invalid_sequence_raises_value_error():
    with pytest.raises(ValueError, match="Invalid amino acid"):
        protein_to_dna("INVALIDZ")


def test_protein_to_dna_with_stop_codon():
    assert protein_to_dna("M*") == "ATGTAA"


def test_dna_to_rna():
    result = dna_to_rna(
        "ATGACTGAATATAAACTTGTTGTTGTTGGTGCTGGTGGTGTTGGTAAATCTGCTCTTACT"
        "ATTCAACTTATTCAAAATCATTTTGTTGATGAATATGATCCTACT"
    )
    expected = (
        "AUGACUGAAUAUAAACUUGUUGUUGUUGGUGCUGGUGGUGUUGGUAAAUCUGCUCUUACU"
        "AUUCAACUUAUUCAAAAUCAUUUUGUUGAUGAAUAUGAUCCUACU"
    )
    assert result == expected


def test_count_nucleotides_dna():
    result = count_nucleotides(
        "ATGACTGAATATAAACTTGTTGTTGTTGGTGCTGGTGGTGTTGGTAAATCTGCTCTTACT"
        "ATTCAACTTATTCAAAATCATTTTGTTGATGAATATGATCCTACT"
    )
    expected = {"A": 27, "C": 14, "G": 20, "T": 44, "U": 0}
    assert result == expected


def test_count_nucleotides_empty():
    result = count_nucleotides("")
    expected = {"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}
    assert result == expected


def test_count_nucleotides_lowercase_rna():
    result = count_nucleotides("auguu")
    expected = {"A": 1, "C": 0, "G": 1, "T": 0, "U": 3}
    assert result == expected


def test_calculate_polarity_and_charge_of_protein():
    polarity_count, charge_count = calculate_polarity_and_charge_of_protein(
        "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT"
    )
    expected_polarity = {"Nonpolar": 19, "Polar": 16, "Unknown": 0, "Stop": 0}
    expected_charge = {"Positive": 3, "Negative": 4, "No charge": 28, "Unknown": 0, "Stop": 0}

    assert polarity_count == expected_polarity
    assert charge_count == expected_charge


def test_calculate_polarity_and_charge_accepts_lowercase():
    polarity_count, charge_count = calculate_polarity_and_charge_of_protein("m*")
    assert polarity_count == {"Nonpolar": 1, "Polar": 0, "Unknown": 0, "Stop": 1}
    assert charge_count == {"Positive": 0, "Negative": 0, "No charge": 1, "Unknown": 0, "Stop": 1}


def test_calculate_polarity_and_charge_invalid_sequence_raises_value_error():
    with pytest.raises(ValueError, match="Invalid amino acid"):
        calculate_polarity_and_charge_of_protein("ABZ")


def test_parse_variant_label_accepts_lowercase_and_prefix():
    parsed = parse_variant_label("p.e50k")
    assert parsed == {
        "label": "p.E50K",
        "reference": "E",
        "position": 50,
        "alternate": "K",
    }


def test_parse_variant_label_invalid_format_raises():
    with pytest.raises(ValueError, match="Invalid variant format"):
        parse_variant_label("E50")


def test_parse_variant_list_strips_and_dedupes():
    result = parse_variant_list(" p.E50K, p.E50K ; p.R214W ")
    assert result == ["p.E50K", "p.R214W"]


def test_analyze_struct_interact_known_protein():
    result = analyze_struct_interact(
        "OPTN", disease_focus="glaucoma", variants=["p.E50K", "p.H486R"]
    )

    assert result["protein"]["gene_name"] == "OPTN"
    assert "autophagy dysfunction" in result["mechanisms"]
    assert len(result["interactions"]) > 0
    assert len(result["variant_assessments"]) == 2
    assert all("pmds_score" in item for item in result["variant_assessments"])
    assert len(result["hypotheses"]) >= 1


def test_analyze_struct_interact_unknown_protein_uses_fallback():
    result = analyze_struct_interact("PROTEIN_X", variants=["p.A10V"])
    assert result["protein"]["gene_name"] == "PROTEIN_X"
    assert any("No curated local record found" in msg for msg in result["limitations"])


def test_analyze_struct_interact_variant_position_bounds_check():
    with pytest.raises(ValueError, match="position exceeds TP53 length"):
        analyze_struct_interact("TP53", variants=["p.R500H"])


def test_format_struct_interact_report_includes_key_sections():
    result = analyze_struct_interact("TP53", disease_focus="cancer", variants=["p.R175H"])
    report = format_struct_interact_report(result)

    assert "StructInteract Report" in report
    assert "Protein: TP53" in report
    assert "Variant PMDS-style assessment:" in report
    assert "Top hypotheses:" in report
