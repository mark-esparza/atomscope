from web_app import _parse_form_data, _render_page, _sequence_report


def test_parse_form_data_decodes_values():
    payload = b"mode=struct&protein=OPTN&variants=p.E50K%2Cp.H486R"
    parsed = _parse_form_data(payload)
    assert parsed["mode"] == "struct"
    assert parsed["protein"] == "OPTN"
    assert parsed["variants"] == "p.E50K,p.H486R"


def test_sequence_report_contains_key_sections():
    report = _sequence_report("MTEYK", "preferred")
    assert "Sequence Console Report" in report
    assert "DNA sequence (preferred codons):" in report
    assert "RNA sequence:" in report
    assert "Polarity count:" in report


def test_render_page_includes_terminal_and_protein_image():
    html = _render_page(
        form_data={
            "mode": "struct",
            "protein": "OPTN",
            "disease_focus": "glaucoma",
            "variants": "p.E50K,p.H486R",
            "sequence": "MTEYK",
            "codon_strategy": "preferred",
        },
        command_preview="python project.py --struct-interact-protein OPTN",
        output_text="ok",
        error_text="",
    )
    assert "Integrated Protein Data Bank Lookup" in html
    assert "Search PDB" in html
    assert "PS C:\\Users\\Marka\\Downloads\\project&gt;" in html
    assert "/static/protein_3d.svg" in html
    assert "https://search.rcsb.org/rcsbsearch/v2/query" in html
