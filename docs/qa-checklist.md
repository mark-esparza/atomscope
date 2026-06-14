# Manual QA checklist

Run before each deploy until the automated browser coverage (`e2e/`) is broad
enough to rely on alone. Open the browser devtools **Console** and watch for any
red errors throughout — there should be none.

## Load & overview
- [ ] Enter `1HSG`, click **Go** → structure loads; Overview shows title/metadata.
- [ ] Enter a name (e.g. `insulin`) → search results list; clicking one loads it.
- [ ] **Upload PDB…** a local `.pdb` file → loads like an RCSB entry.
- [ ] Load a large assembly (e.g. `6BCX`) → chain picker appears; pick a chain → loads.

## Per-tab
- [ ] **3D Atomic Viewer** renders; rotate/zoom; color modes switch; toggles work.
- [ ] **Interactions** — pick a bound molecule → contacts + per-residue summary.
- [ ] **Pockets** — Detect pockets → ranked cavities (or a clear "none" message);
      method/interpretation card shows.
- [ ] **Evolution** — Analyze conservation → residues/pockets scored, or a clear
      "unavailable" message; method card shows.
- [ ] **Docking** — dock a chemical into a site → job runs (status updates while
      polling) → pose + interactions + methods render.
- [ ] **Batch screen** — a few chemicals → ranked table.
- [ ] **Chemical** — look up a drug → properties + druglikeness.

## Reports & exports
- [ ] **Compile full report** populates.
- [ ] Export **.txt**, **.json**, **.csv** (after a screen), **.pdb** (after a dock)
      each download and open correctly.

## Footers & pages
- [ ] Footer links open `/api.html` (renders endpoints), `/privacy.html`, `/terms.html`.
- [ ] `/api/docs`, `/api/version` return JSON.

## Accessibility (spot check)
- [ ] Tab through controls — focus outline is visible on each.
- [ ] Status messages are announced (aria-live region updates).
- [ ] Narrow the window (< 820 px) — layout stacks instead of crowding.
