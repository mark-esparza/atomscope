import argparse
import html
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from project import (
    analyze_struct_interact,
    calculate_polarity_and_charge_of_protein,
    count_nucleotides,
    dna_to_rna,
    format_struct_interact_report,
    normalize_protein_sequence,
    protein_to_dna,
)

ROOT_DIR = Path(__file__).resolve().parent
STATIC_DIR = ROOT_DIR / "static"
HOST = "127.0.0.1"
PORT = 8000


def _as_text(value, default=""):
    if value is None:
        return default
    return str(value).strip()


def _parse_form_data(body_bytes):
    parsed = parse_qs(body_bytes.decode("utf-8"), keep_blank_values=True)
    return {key: values[0] if values else "" for key, values in parsed.items()}


def _sequence_report(sequence, strategy):
    sequence = normalize_protein_sequence(sequence)
    dna_seq = protein_to_dna(sequence, strategy=strategy)
    rna_seq = dna_to_rna(dna_seq)
    dna_count = count_nucleotides(dna_seq)
    rna_count = count_nucleotides(rna_seq)
    polarity_count, charge_count = calculate_polarity_and_charge_of_protein(sequence)

    return "\n".join(
        [
            "Sequence Console Report",
            "=======================",
            f"Input sequence length: {len(sequence)} aa",
            f"Codon strategy: {strategy}",
            "",
            f"DNA sequence ({strategy} codons): {dna_seq}",
            f"RNA sequence: {rna_seq}",
            f"Nucleotide counts (DNA): {dna_count}",
            f"Nucleotide counts (RNA): {rna_count}",
            f"Polarity count: {polarity_count}",
            f"Charge count: {charge_count}",
        ]
    )


def _render_page(form_data, command_preview, output_text, error_text):
    mode = form_data.get("mode", "struct")
    protein = form_data.get("protein", "OPTN")
    disease_focus = form_data.get("disease_focus", "glaucoma")
    variants = form_data.get("variants", "p.E50K,p.H486R")
    sequence = form_data.get("sequence", "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT")
    codon_strategy = form_data.get("codon_strategy", "preferred")

    output_block = html.escape(output_text)
    error_block = html.escape(error_text)
    command_block = html.escape(command_preview)

    # Lightweight command preview updates so the web flow still feels script-driven.
    js_state = {
        "struct": "python project.py --struct-interact-protein {protein} --disease-focus "
        "{disease_focus} --variants {variants}",
        "sequence": "python project.py {sequence} --codon-strategy {codon_strategy}",
    }

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>StructInteract | PDB-Inspired Workbench</title>
  <style>
    :root {{
      --bg: #eff4fb;
      --panel: #ffffff;
      --panel-soft: #f7faff;
      --line: #d0dff4;
      --text: #1f324f;
      --muted: #4f6788;
      --accent: #0b61c8;
      --accent-2: #1f8df2;
      --danger: #b63043;
      --terminal-bg: #07152e;
      --terminal-line: #254573;
      --terminal-text: #d9e9ff;
      --terminal-muted: #93afd5;
      --terminal-accent: #9ee8ff;
      --shadow: rgba(10, 37, 74, 0.15);
    }}

    * {{
      box-sizing: border-box;
    }}

    body {{
      margin: 0;
      min-height: 100vh;
      background:
        radial-gradient(circle at 100% 0%, #e6eefb 0%, transparent 42%),
        radial-gradient(circle at 0% 100%, #e9f1ff 0%, transparent 40%),
        linear-gradient(145deg, #f5f8fd 0%, #f3f7fd 46%, #f1f6fc 100%);
      color: var(--text);
      font-family: "Consolas", "Cascadia Mono", "Menlo", monospace;
      padding: 0;
    }}

    .topbar {{
      border-bottom: 1px solid var(--line);
      background: linear-gradient(180deg, #ffffff 0%, #f4f8ff 100%);
      box-shadow: 0 6px 20px var(--shadow);
    }}

    .topbar-inner {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 14px 18px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 14px;
      flex-wrap: wrap;
    }}

    .brand {{
      font-size: 16px;
      color: var(--accent);
      font-weight: 700;
      letter-spacing: 0.02em;
    }}

    .nav-tag {{
      border: 1px solid var(--line);
      color: var(--muted);
      background: #fff;
      padding: 6px 10px;
      border-radius: 999px;
      font-size: 12px;
    }}

    .page {{
      max-width: 1240px;
      margin: 18px auto 30px auto;
      padding: 0 18px;
      display: grid;
      gap: 14px;
    }}

    .card {{
      border: 1px solid var(--line);
      border-radius: 14px;
      background: var(--panel);
      box-shadow: 0 10px 24px var(--shadow);
      overflow: hidden;
    }}

    .card-head {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
      background: linear-gradient(180deg, #ffffff 0%, #f4f8ff 100%);
    }}

    .card-title {{
      color: var(--accent);
      font-size: 14px;
      margin: 0;
    }}

    .lookup-body {{
      padding: 14px;
      background: var(--panel-soft);
    }}

    .lookup-grid {{
      display: grid;
      grid-template-columns: minmax(220px, 1.2fr) auto auto;
      gap: 8px;
      align-items: center;
    }}

    .lookup-grid input {{
      width: 100%;
      border: 1px solid var(--line);
      background: #fff;
      color: var(--text);
      border-radius: 8px;
      padding: 10px 11px;
      font-family: inherit;
      font-size: 14px;
    }}

    .lookup-actions {{
      display: flex;
      gap: 8px;
    }}

    .btn {{
      border: 1px solid var(--accent);
      background: linear-gradient(180deg, #1e84ea 0%, #0b61c8 100%);
      color: #fff;
      border-radius: 8px;
      padding: 9px 13px;
      font-family: inherit;
      cursor: pointer;
    }}

    .btn:hover {{
      filter: brightness(1.05);
    }}

    .btn-ghost {{
      border: 1px solid var(--line);
      background: #fff;
      color: var(--accent);
      border-radius: 8px;
      padding: 9px 11px;
      font-family: inherit;
      cursor: pointer;
    }}

    .lookup-status {{
      color: var(--muted);
      font-size: 12px;
      margin-top: 8px;
      min-height: 18px;
    }}

    .lookup-status.error {{
      color: var(--danger);
    }}

    .pdb-results {{
      margin-top: 10px;
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(230px, 1fr));
      gap: 10px;
    }}

    .pdb-card {{
      border: 1px solid var(--line);
      background: #fff;
      border-radius: 10px;
      padding: 10px;
      display: grid;
      gap: 7px;
    }}

    .pdb-card-id {{
      font-size: 14px;
      color: var(--accent);
      font-weight: 700;
      text-decoration: none;
    }}

    .pdb-card-title {{
      font-size: 12px;
      color: var(--text);
      line-height: 1.35;
      min-height: 34px;
    }}

    .pdb-meta {{
      font-size: 11px;
      color: var(--muted);
      line-height: 1.35;
    }}

    .pdb-card button {{
      width: fit-content;
      border: 1px solid var(--line);
      background: #fff;
      color: var(--accent);
      border-radius: 7px;
      padding: 6px 9px;
      font-family: inherit;
      cursor: pointer;
    }}

    .layout {{
      display: grid;
      grid-template-columns: 1.42fr 1fr;
      gap: 14px;
    }}

    .shell {{
      border: 1px solid var(--terminal-line);
      border-radius: 14px;
      background: var(--terminal-bg);
      box-shadow: 0 18px 28px rgba(0, 0, 0, 0.24);
      overflow: hidden;
    }}

    .shell-header {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 10px 12px;
      border-bottom: 1px solid var(--terminal-line);
      background: linear-gradient(180deg, #13264e 0%, #102042 100%);
    }}

    .window-lights {{
      display: flex;
      gap: 8px;
    }}

    .window-lights span {{
      width: 10px;
      height: 10px;
      border-radius: 50%;
      display: inline-block;
    }}

    .window-lights span:nth-child(1) {{ background: #ff6b6b; }}
    .window-lights span:nth-child(2) {{ background: #f6c453; }}
    .window-lights span:nth-child(3) {{ background: #42d288; }}

    .shell-title {{
      color: var(--terminal-muted);
      font-size: 12px;
      letter-spacing: 0.03em;
    }}

    .panel {{
      padding: 14px;
    }}

    .prompt {{
      color: var(--terminal-accent);
      margin: 0 0 10px 0;
      word-break: break-word;
    }}

    form {{
      display: grid;
      gap: 10px;
      margin-bottom: 14px;
      padding-bottom: 14px;
      border-bottom: 1px dashed var(--terminal-line);
    }}

    .row {{
      display: grid;
      grid-template-columns: 160px 1fr;
      gap: 10px;
      align-items: center;
    }}

    label {{
      color: var(--terminal-muted);
      font-size: 13px;
    }}

    form input, form select {{
      width: 100%;
      border: 1px solid var(--terminal-line);
      background: #0b1a38;
      color: var(--terminal-text);
      border-radius: 7px;
      padding: 8px 10px;
      font-family: inherit;
      font-size: 14px;
    }}

    .run-button {{
      width: fit-content;
      border: 1px solid #4f8dcf;
      background: #153160;
      color: #d5eaff;
      border-radius: 7px;
      padding: 8px 12px;
      font-family: inherit;
      cursor: pointer;
    }}

    .run-button:hover {{
      background: #1b3d77;
    }}

    pre {{
      margin: 0;
      padding: 14px;
      border: 1px solid var(--terminal-line);
      border-radius: 10px;
      background: #060f21;
      color: var(--terminal-text);
      white-space: pre-wrap;
      word-break: break-word;
      line-height: 1.35;
      min-height: 280px;
      max-height: 62vh;
      overflow: auto;
    }}

    .error {{
      color: var(--danger);
      margin-bottom: 10px;
      white-space: pre-wrap;
    }}

    .visual {{
      border: 1px solid var(--line);
      border-radius: 14px;
      background: #fff;
      box-shadow: 0 10px 24px var(--shadow);
      overflow: hidden;
    }}

    .visual-head {{
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
      background: linear-gradient(180deg, #ffffff 0%, #f4f8ff 100%);
      color: var(--accent);
      font-size: 14px;
      font-weight: 700;
    }}

    .visual-body {{
      padding: 14px;
    }}

    .visual h2 {{
      margin: 0 0 6px 0;
      font-size: 15px;
      color: var(--text);
    }}

    .visual p {{
      margin: 0 0 14px 0;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.4;
    }}

    .image-frame {{
      border: 1px solid #c7daf4;
      border-radius: 12px;
      padding: 12px;
      background: linear-gradient(160deg, #f4f8ff 0%, #ecf3ff 100%);
      box-shadow: inset 0 0 12px rgba(45, 88, 145, 0.12);
    }}

    .image-frame img {{
      width: 100%;
      display: block;
      border-radius: 8px;
    }}

    .legend {{
      margin-top: 12px;
      font-size: 12px;
      color: var(--muted);
      line-height: 1.45;
    }}

    @media (max-width: 960px) {{
      body {{
        padding: 0;
      }}

      .topbar-inner {{
        padding: 12px;
      }}

      .page {{
        padding: 0 12px 16px 12px;
      }}

      .lookup-grid {{
        grid-template-columns: 1fr;
      }}

      .layout {{
        grid-template-columns: 1fr;
      }}

      .row {{
        grid-template-columns: 1fr;
        gap: 6px;
      }}
    }}
  </style>
</head>
<body>
  <header class="topbar">
    <div class="topbar-inner">
      <div class="brand">StructInteract Workbench</div>
      <div class="nav-tag">PDB-style research layout</div>
      <div class="nav-tag">Research-only interpretation</div>
    </div>
  </header>

  <main class="page">
    <section class="card">
      <div class="card-head">
        <h1 class="card-title">Integrated Protein Data Bank Lookup</h1>
      </div>
      <div class="lookup-body">
        <div class="lookup-grid">
          <input id="pdb-query" type="text" value="{html.escape(protein)}" placeholder="Search by protein name, gene, organism, or keyword (e.g., hemoglobin)">
          <div class="lookup-actions">
            <button class="btn" type="button" id="pdb-search-btn">Search PDB</button>
            <button class="btn-ghost" type="button" id="use-query-btn">Use Query in Analysis</button>
          </div>
          <div class="lookup-actions">
            <a class="btn-ghost" href="https://www.rcsb.org" target="_blank" rel="noopener noreferrer">Open RCSB</a>
          </div>
        </div>
        <div id="pdb-status" class="lookup-status">Ready for RCSB Protein Data Bank search.</div>
        <div id="pdb-results" class="pdb-results"></div>
      </div>
    </section>

    <div class="layout">
      <section class="shell">
        <div class="shell-header">
          <div class="window-lights" aria-hidden="true"><span></span><span></span><span></span></div>
          <div class="shell-title">StructInteract console bridge (command-script style)</div>
        </div>
        <div class="panel">
          <form method="post" action="/">
            <div class="row">
              <label for="mode">Mode</label>
              <select id="mode" name="mode">
                <option value="struct" {"selected" if mode == "struct" else ""}>StructInteract analysis</option>
                <option value="sequence" {"selected" if mode == "sequence" else ""}>Sequence console</option>
              </select>
            </div>

            <div class="row">
              <label for="protein">Protein or PDB ID</label>
              <input id="protein" name="protein" value="{html.escape(protein)}" placeholder="OPTN or 4HHB">
            </div>

            <div class="row">
              <label for="disease_focus">Disease focus</label>
              <input id="disease_focus" name="disease_focus" value="{html.escape(disease_focus)}" placeholder="glaucoma">
            </div>

            <div class="row">
              <label for="variants">Variants (comma)</label>
              <input id="variants" name="variants" value="{html.escape(variants)}" placeholder="p.E50K,p.H486R">
            </div>

            <div class="row">
              <label for="sequence">Protein sequence</label>
              <input id="sequence" name="sequence" value="{html.escape(sequence)}" placeholder="MTEYKLV...">
            </div>

            <div class="row">
              <label for="codon_strategy">Codon strategy</label>
              <select id="codon_strategy" name="codon_strategy">
                <option value="preferred" {"selected" if codon_strategy == "preferred" else ""}>preferred</option>
                <option value="random" {"selected" if codon_strategy == "random" else ""}>random</option>
              </select>
            </div>

            <button class="run-button" type="submit">Run</button>
          </form>

          <p class="prompt">PS C:\\Users\\Marka\\Downloads\\project&gt; <span id="cmd-preview">{command_block}</span></p>
          {"<div class='error'>" + error_block + "</div>" if error_text else ""}
          <pre>{output_block}</pre>
        </div>
      </section>

      <aside class="visual">
        <div class="visual-head">Protein 3D Context Panel</div>
        <div class="visual-body">
          <h2>Local structural render</h2>
          <p>
            A built-in 3D-style protein image stays visible while you run analysis and search PDB entries.
            Use PDB search results above to jump to structure pages quickly.
          </p>
          <div class="image-frame">
            <img src="/static/protein_3d.svg" alt="3D style protein structure illustration">
          </div>
          <div class="legend">
            RCSB lookup is integrated above via public API endpoints.
            Use “Use ID in analysis” on any result card to inject the selected structure ID into the command console.
          </div>
        </div>
      </aside>
    </div>
  </main>

  <script>
    const commandFormats = {json.dumps(js_state)};
    const modeEl = document.getElementById("mode");
    const proteinEl = document.getElementById("protein");
    const diseaseEl = document.getElementById("disease_focus");
    const variantsEl = document.getElementById("variants");
    const sequenceEl = document.getElementById("sequence");
    const codonEl = document.getElementById("codon_strategy");
    const previewEl = document.getElementById("cmd-preview");
    const pdbQueryEl = document.getElementById("pdb-query");
    const pdbStatusEl = document.getElementById("pdb-status");
    const pdbResultsEl = document.getElementById("pdb-results");
    const pdbSearchBtn = document.getElementById("pdb-search-btn");
    const useQueryBtn = document.getElementById("use-query-btn");

    function escapeHtml(value) {{
      return String(value)
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;");
    }}

    function setPdbStatus(message, isError = false) {{
      pdbStatusEl.textContent = message;
      pdbStatusEl.classList.toggle("error", isError);
    }}

    function normalizeVariantToken(rawToken) {{
      const cleaned = (rawToken || "").toUpperCase().replace(/[^A-Z0-9*]/g, "");
      const match = cleaned.match(/^([A-Z*])(\\d+)([A-Z*])$/);
      if (!match) {{
        return null;
      }}
      return `p.${{match[1]}}${{match[2]}}${{match[3]}}`;
    }}

    function extractVariantsFromText(rawText) {{
      const text = (rawText || "").toUpperCase();
      const candidates = text.match(/[A-Z*]\\d+[A-Z*]/g) || [];
      const variants = [];
      const seen = new Set();
      for (const token of candidates) {{
        const normalized = normalizeVariantToken(token);
        if (!normalized || seen.has(normalized)) {{
          continue;
        }}
        seen.add(normalized);
        variants.push(normalized);
      }}
      return variants.slice(0, 6);
    }}

    function sanitizeProteinSequence(rawSequence) {{
      const collapsed = (rawSequence || "").toUpperCase().replace(/[^A-Z]/g, "");
      if (!collapsed) {{
        return "";
      }}
      const allowed = new Set("ACDEFGHIKLMNPQRSTVWYX*".split(""));
      return [...collapsed]
        .map((char) => (allowed.has(char) ? char : "X"))
        .join("");
    }}

    function inferDiseaseFocusFromPdbText(rawText) {{
      const text = (rawText || "").toLowerCase();
      const mapping = [
        {{
          label: "cancer",
          keywords: ["cancer", "oncogene", "tumor", "carcinoma", "p53", "brca", "kinase inhibitor"],
        }},
        {{
          label: "alzheimer disease",
          keywords: ["alzheimer", "amyloid", "apoe", "tau", "neurodegeneration"],
        }},
        {{
          label: "glaucoma",
          keywords: ["glaucoma", "myocilin", "optineurin", "retina", "optic nerve"],
        }},
        {{
          label: "als",
          keywords: ["als", "amyotrophic", "sod1", "tbk1", "motor neuron"],
        }},
        {{
          label: "infectious disease",
          keywords: ["virus", "viral", "bacterial", "pathogen", "sars", "hiv", "influenza"],
        }},
        {{
          label: "metabolic disease",
          keywords: ["metabolic", "enzyme deficiency", "diabetes", "lipid metabolism"],
        }},
      ];

      for (const rule of mapping) {{
        if (rule.keywords.some((keyword) => text.includes(keyword))) {{
          return rule.label;
        }}
      }}
      return "structural biology";
    }}

    function pickBestProteinEntity(polymerEntities) {{
      const proteins = polymerEntities.filter((entity) => {{
        const polymerType = (
          entity?.entity_poly?.rcsb_entity_polymer_type ||
          entity?.entity_poly?.type ||
          ""
        ).toLowerCase();
        return polymerType.includes("protein") || polymerType.includes("polypeptide");
      }});

      const candidates = proteins.length ? proteins : polymerEntities;
      candidates.sort((left, right) => {{
        const leftLength = (left?.entity_poly?.pdbx_seq_one_letter_code_can || "").length;
        const rightLength = (right?.entity_poly?.pdbx_seq_one_letter_code_can || "").length;
        return rightLength - leftLength;
      }});
      return candidates[0] || null;
    }}

    function getUniprotId(entity) {{
      const directIds = entity?.rcsb_polymer_entity_container_identifiers?.uniprot_ids || [];
      if (directIds.length) {{
        return directIds[0];
      }}

      const refs =
        entity?.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers || [];
      for (const ref of refs) {{
        const dbName = String(ref?.database_name || "").toLowerCase();
        if (dbName.includes("uniprot")) {{
          return ref?.database_accession || "";
        }}
      }}
      return "";
    }}

    async function fetchPolymerEntity(entryId, entityId) {{
      try {{
        const response = await fetch(
          `https://data.rcsb.org/rest/v1/core/polymer_entity/${{entryId}}/${{entityId}}`
        );
        if (!response.ok) {{
          return null;
        }}
        return await response.json();
      }} catch (_error) {{
        return null;
      }}
    }}

    function updateFormFromPdbEntry(entry) {{
      modeEl.value = "struct";
      proteinEl.value = (entry.analysisProtein || entry.id || "").toUpperCase();
      diseaseEl.value = entry.diseaseFocus || "structural biology";
      variantsEl.value = (entry.variants || []).join(",");
      sequenceEl.value = entry.sequence || "";
      refreshPreview();

      const variantText = entry.variants?.length
        ? entry.variants.join(", ")
        : "no mutation annotations found";
      setPdbStatus(
        `Loaded ${{entry.id}} -> protein=${{proteinEl.value}}, focus=${{diseaseEl.value}}, variants=${{variantText}}, sequence length=${{sequenceEl.value.length}}`
      );
    }}

    function refreshPreview() {{
      const mode = modeEl.value;
      const template = commandFormats[mode];
      const rendered = template
        .replace("{protein}", proteinEl.value || "OPTN")
        .replace("{disease_focus}", diseaseEl.value || "glaucoma")
        .replace("{variants}", variantsEl.value || "p.E50K,p.H486R")
        .replace("{sequence}", sequenceEl.value || "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT")
        .replace("{codon_strategy}", codonEl.value || "preferred");
      previewEl.textContent = rendered;
    }}

    [modeEl, proteinEl, diseaseEl, variantsEl, sequenceEl, codonEl].forEach((el) => {{
      el.addEventListener("input", refreshPreview);
      el.addEventListener("change", refreshPreview);
    }});

    useQueryBtn.addEventListener("click", () => {{
      proteinEl.value = (pdbQueryEl.value || "").trim().toUpperCase() || "OPTN";
      refreshPreview();
      setPdbStatus("Query copied into analysis field. Choose a PDB result card to auto-fill all fields.");
    }});

    async function fetchEntrySummary(entryId) {{
      try {{
        const detailResp = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${{entryId}}`);
        if (!detailResp.ok) {{
          return {{
            id: entryId,
            title: "Unable to load entry detail.",
            method: "N/A",
            resolution: "N/A",
            releaseDate: "N/A",
            organism: "N/A",
            sequence: "",
            variants: [],
            diseaseFocus: "structural biology",
            analysisProtein: entryId,
          }};
        }}
        const detail = await detailResp.json();
        const polymerEntityIds =
          detail?.rcsb_entry_container_identifiers?.polymer_entity_ids || [];
        const polymerEntities = (
          await Promise.all(
            polymerEntityIds.map((entityId) => fetchPolymerEntity(entryId, entityId))
          )
        ).filter(Boolean);

        const bestEntity = pickBestProteinEntity(polymerEntities);
        const entityDescription =
          bestEntity?.rcsb_polymer_entity?.pdbx_description ||
          bestEntity?.entity_poly?.pdbx_seq_one_letter_code ||
          "";
        const organism =
          bestEntity?.rcsb_entity_source_organism?.[0]?.ncbi_scientific_name || "N/A";
        const sequence = sanitizeProteinSequence(
          bestEntity?.entity_poly?.pdbx_seq_one_letter_code_can ||
            bestEntity?.entity_poly?.pdbx_seq_one_letter_code ||
            ""
        );
        const uniprotId = getUniprotId(bestEntity);
        const explicitMutationText = bestEntity?.rcsb_polymer_entity?.pdbx_mutation || "";
        let variants = extractVariantsFromText(explicitMutationText);
        if (!variants.length) {{
          const titleAndDescription = [detail?.struct?.title || "", entityDescription || ""].join(" ");
          if (/(mutant|mutation|variant)/i.test(titleAndDescription)) {{
            variants = extractVariantsFromText(titleAndDescription);
          }}
        }}
        const diseaseFocus = inferDiseaseFocusFromPdbText(
          [detail?.struct?.title || "", entityDescription || "", organism || ""].join(" ")
        );

        const analysisProtein = uniprotId || entryId;
        const method = detail?.exptl?.[0]?.method || "N/A";
        const resolution =
          detail?.rcsb_entry_info?.resolution_combined?.[0] !== undefined
            ? `${{detail.rcsb_entry_info.resolution_combined[0]}} A`
            : "N/A";
        const releaseDate = detail?.rcsb_accession_info?.initial_release_date
          ? detail.rcsb_accession_info.initial_release_date.slice(0, 10)
          : "N/A";
        return {{
          id: entryId,
          title: detail?.struct?.title || "No title available",
          method,
          resolution,
          releaseDate,
          organism,
          sequence,
          variants,
          diseaseFocus,
          analysisProtein,
        }};
      }} catch (_err) {{
        return {{
          id: entryId,
          title: "Unable to load entry detail.",
          method: "N/A",
          resolution: "N/A",
          releaseDate: "N/A",
          organism: "N/A",
          sequence: "",
          variants: [],
          diseaseFocus: "structural biology",
          analysisProtein: entryId,
        }};
      }}
    }}

    function renderPdbCards(entries) {{
      pdbResultsEl.innerHTML = "";

      for (const entry of entries) {{
        const card = document.createElement("article");
        card.className = "pdb-card";

        const link = document.createElement("a");
        link.className = "pdb-card-id";
        link.href = `https://www.rcsb.org/structure/${{entry.id}}`;
        link.target = "_blank";
        link.rel = "noopener noreferrer";
        link.textContent = entry.id;

        const title = document.createElement("div");
        title.className = "pdb-card-title";
        title.textContent = entry.title;

        const meta = document.createElement("div");
        meta.className = "pdb-meta";
        meta.innerHTML =
          `Method: ${{escapeHtml(entry.method)}}<br>` +
          `Resolution: ${{escapeHtml(entry.resolution)}}<br>` +
          `Release: ${{escapeHtml(entry.releaseDate)}}<br>` +
          `Organism: ${{escapeHtml(entry.organism || "N/A")}}<br>` +
          `Focus hint: ${{escapeHtml(entry.diseaseFocus || "structural biology")}}`;

        const useButton = document.createElement("button");
        useButton.type = "button";
        useButton.textContent = "Use ID in analysis";
        useButton.addEventListener("click", () => {{
          updateFormFromPdbEntry(entry);
        }});

        card.appendChild(link);
        card.appendChild(title);
        card.appendChild(meta);
        card.appendChild(useButton);
        pdbResultsEl.appendChild(card);
      }}
    }}

    async function runPdbSearch() {{
      const query = (pdbQueryEl.value || "").trim();
      if (!query) {{
        setPdbStatus("Enter a protein name or keyword before searching.", true);
        return;
      }}

      setPdbStatus("Searching RCSB Protein Data Bank...");
      pdbResultsEl.innerHTML = "";

      const payload = {{
        query: {{
          type: "terminal",
          service: "full_text",
          parameters: {{
            value: query
          }}
        }},
        return_type: "entry",
        request_options: {{
          paginate: {{
            start: 0,
            rows: 8
          }}
        }}
      }};

      try {{
        const response = await fetch("https://search.rcsb.org/rcsbsearch/v2/query", {{
          method: "POST",
          headers: {{
            "Content-Type": "application/json"
          }},
          body: JSON.stringify(payload)
        }});

        if (!response.ok) {{
          throw new Error(`Search request failed with status ${{response.status}}`);
        }}

        const searchData = await response.json();
        const ids = (searchData?.result_set || [])
          .map((entry) => entry.identifier)
          .filter(Boolean)
          .slice(0, 6);

        if (!ids.length) {{
          setPdbStatus("No PDB entries found for that query.");
          return;
        }}

        const detailedEntries = await Promise.all(ids.map((id) => fetchEntrySummary(id)));
        renderPdbCards(detailedEntries);
        setPdbStatus(`Found ${{detailedEntries.length}} matching PDB entries.`);
      }} catch (error) {{
        setPdbStatus(`PDB lookup failed: ${{error.message}}`, true);
      }}
    }}

    pdbSearchBtn.addEventListener("click", runPdbSearch);
    pdbQueryEl.addEventListener("keydown", (event) => {{
      if (event.key === "Enter") {{
        event.preventDefault();
        runPdbSearch();
      }}
    }});

    refreshPreview();
    runPdbSearch();
  </script>
</body>
</html>
"""


class StructInteractWebHandler(BaseHTTPRequestHandler):
    server_version = "StructInteractWeb/1.0"

    def _send_html(self, content):
        payload = content.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _send_static(self, relative_path):
        safe_path = STATIC_DIR / relative_path
        safe_path = safe_path.resolve()
        if STATIC_DIR not in safe_path.parents and safe_path != STATIC_DIR:
            self.send_error(403, "Forbidden")
            return
        if not safe_path.is_file():
            self.send_error(404, "Not Found")
            return

        extension = safe_path.suffix.lower()
        content_type = {
            ".svg": "image/svg+xml",
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".webp": "image/webp",
        }.get(extension, "application/octet-stream")

        data = safe_path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path.startswith("/static/"):
            relative = parsed.path.replace("/static/", "", 1)
            self._send_static(relative)
            return

        default_form = {
            "mode": "struct",
            "protein": "OPTN",
            "disease_focus": "glaucoma",
            "variants": "p.E50K,p.H486R",
            "sequence": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT",
            "codon_strategy": "preferred",
        }
        preview = (
            "python project.py --struct-interact-protein OPTN "
            "--disease-focus glaucoma --variants p.E50K,p.H486R"
        )
        welcome = (
            "StructInteract Web Console ready.\n"
            "Submit the form above to run analysis.\n\n"
            "Tip: switch mode to Sequence console if you want DNA/RNA outputs."
        )
        self._send_html(_render_page(default_form, preview, welcome, ""))

    def do_POST(self):
        if self.path != "/":
            self.send_error(404, "Not Found")
            return

        content_length = int(self.headers.get("Content-Length", "0"))
        raw_body = self.rfile.read(content_length)
        form_data = _parse_form_data(raw_body)

        mode = _as_text(form_data.get("mode"), "struct")
        protein = _as_text(form_data.get("protein"), "OPTN")
        disease_focus = _as_text(form_data.get("disease_focus"))
        variants = _as_text(form_data.get("variants"))
        sequence = _as_text(form_data.get("sequence"))
        codon_strategy = _as_text(form_data.get("codon_strategy"), "preferred")

        if mode == "sequence":
            preview = (
                f"python project.py {sequence or 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT'} "
                f"--codon-strategy {codon_strategy}"
            )
        else:
            preview = (
                f"python project.py --struct-interact-protein {protein or 'OPTN'} "
                f"--disease-focus {disease_focus or 'glaucoma'} "
                f"--variants {variants or 'p.E50K,p.H486R'}"
            )

        try:
            if mode == "sequence":
                report = _sequence_report(sequence, codon_strategy)
            else:
                variants_input = variants if variants else None
                disease_input = disease_focus if disease_focus else None
                result = analyze_struct_interact(
                    protein_id=protein or "OPTN",
                    disease_focus=disease_input,
                    variants=variants_input,
                )
                report = format_struct_interact_report(result)

            self._send_html(_render_page(form_data, preview, report, ""))
        except Exception as exc:  # noqa: BLE001
            error_message = f"Execution error: {exc}"
            self._send_html(_render_page(form_data, preview, "", error_message))


def run_server(host=HOST, port=PORT):
    server = ThreadingHTTPServer((host, port), StructInteractWebHandler)
    print(f"Serving StructInteract Web Console at http://{host}:{port}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server.")
    finally:
        server.server_close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the StructInteract terminal-style local web app."
    )
    parser.add_argument("--host", default=HOST, help="Host to bind (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=PORT, help="Port to bind (default: 8000)")
    args = parser.parse_args()
    run_server(host=args.host, port=args.port)
