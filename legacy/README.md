# Protein Sequencer + StructInteract MVP

A Python command-line tool with:
- sequence utilities (protein -> DNA/RNA + composition analysis)
- a research-only StructInteract MVP for single-protein disease and interactome analysis

## Features

### Sequence utilities
- Convert protein sequence -> DNA sequence
- Convert DNA sequence -> RNA sequence
- Count nucleotide frequencies in DNA/RNA output
- Calculate protein polarity and charge counts
- Run all analyses at once from a menu
- Choose codon strategy for translation:
  - `preferred`: deterministic single codon per amino acid
  - `random`: randomly chooses from synonymous codons

### StructInteract MVP
- Analyze any input protein symbol with local curated context
- Map disease mechanisms, pathways, and interaction partners
- Parse and score protein variants with a PMDS-style composite score
- Highlight top likely disruption drivers per variant
- Generate concise research hypotheses
- Produce a text report for export/use in notes

## Requirements

- Python 3.8+
- `pytest` for running tests

## Installation

1. Open this project folder.
2. (Optional) Create and activate a virtual environment.
3. Install test dependency:

```bash
pip install pytest
```

## Usage

### 1) Sequence menu mode

Run the tool by passing a protein sequence:

```bash
python project.py MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT
```

Use random codon selection:

```bash
python project.py MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT --codon-strategy random
```

Menu options:

1. Convert to DNA
2. Convert to RNA
3. Count Nucleotides
4. Calculate Polarity and Charge of Protein
5. Analyze All
6. Exit

### 2) StructInteract mode

Run a disease/interactome report for a protein:

```bash
python project.py --struct-interact-protein OPTN --disease-focus glaucoma --variants p.E50K,p.R214W
```

Variant input rules:
- Accepts forms like `p.E50K` or `E50K`
- Comma or semicolon separated lists are supported
- Duplicate variants are ignored

If `--variants` is omitted for a curated protein, the tool auto-prioritizes a few known variants.

### 3) Website mode (terminal-style UI + 3D protein image)

Start the local web app:

```bash
python web_app.py
```

Then open:

```text
http://127.0.0.1:8000
```

This web version keeps a command-console look and includes a bundled local 3D-style protein image panel.
It also includes an integrated RCSB Protein Data Bank lookup bar so users can search proteins and load PDB IDs into analysis.
Selecting **Use ID in analysis** now auto-populates protein, disease focus, variants, and protein sequence fields from RCSB entry metadata.
If port `8000` is already in use:

```bash
python web_app.py --port 8010
```

## Input Rules

- Input is normalized to uppercase and whitespace is removed.
- Valid amino acids are:
  - `ACDEFGHIKLMNPQRSTVWY`
  - `X` (unknown)
  - `*` (stop)
- Invalid characters trigger a clear error message.

## Running Tests

```bash
pytest -q
```

## Biological Notes and Limitations

- Protein -> DNA is not uniquely defined in biology because many amino acids have multiple codons.
- `preferred` mode uses a single representative codon for reproducibility.
- `random` mode better reflects codon degeneracy, but results vary run-to-run.
- This project does not perform organism-specific codon optimization or full translation framework logic.
- StructInteract scoring is heuristic and designed for research hypothesis generation only.
- Curated protein knowledge is local to this repository and intentionally limited in scope.
- This project is not intended for clinical decision-making.

## Project Files

- `project.py`: main CLI and analysis functions
- `web_app.py`: local web server for terminal-style browser UI
- `static/protein_3d.svg`: bundled 3D-style protein image used by the web UI
- `test_project.py`: unit tests
- `README.md`: project documentation
