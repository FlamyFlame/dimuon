# Task 08 — Systematic-uncertainty variation framework

**Roadmap:** `analysis_roadmap_2026_06.md` chain [9]
**Depends on:** task 05 (corrected nominal). **Blocked input from
user:** nominal muon WP choice (Q2.8); MCP Run 3 HI recommendations
(Q2.7). **Reviewer:** `/review-analysis-code`.

## Objective

Set up the machinery to repeat the corrected crossx/R_AA measurement
under systematic variations and produce per-source uncertainty tables,
following the Run 2 note pattern (kb `run2_dimuon_note.md`
§Systematics: vary → ratio to nominal → smooth fit → quadrature sum).

## Source list (adapted from Run 2; refine with user)

| Source | Variation | Run 2 size (widths) |
|---|---|---|
| Muon WP | medium ↔ tight | 0.5–2% |
| Reco eff | ±1σ (MCP or MC-stat from dummy/final ε_reco) | ~0.1% |
| Trig eff | ε^nc fit-parameter ±1σ; fit vs interpolation below 8 GeV | small |
| dR trig correction | ε_dR = 1 vs MC cross-term | new (Run 3 specific) |
| Lumi / T_AA | normalization only | 1.5–1.6% lumi |
| Signal-region cuts | vary minv window / pair pT threshold / dR cut | n/a (crossx-specific) |
| Event selection (PbPb) | nominal vs `_alt` banana cuts (already produced: `event_sel_cuts_pbpb_20YY_alt.root`) | new |
| Centrality calib | FCal scale factor variation | new |
| Purity | template-fit fraction (task 07) | <2% (Run 2) |

## Design

- Variations enter as config flags on the existing RDF classes (WP
  flag, eff-file path, cut values) — avoid code duplication; one
  "variation registry" script in `pipelines/` that loops variations,
  reruns RDF + crossx with suffixed outputs
  (`_syst_<source>_<updown>`), and computes ratios.
- Ratio + smoothing plotter modeled on Run 2 note lower-panel fits.
- Output: per-source ratio ROOT files + a summary markdown/LaTeX table
  generator for the note.

## Steps

1. Create implementation tracking doc (multi-cycle work) with the
   variation list agreed with the user.
2. Implement registry + 2–3 cheapest variations first (WP, alt event
   selection, eff off) to validate the machinery.
3. Add remaining sources as their inputs become real (replace dummies).
4. Update roadmap ledger.

## Verification

- "Efficiency off" variation reproduces the uncorrected crossx exactly.
- Nominal rerun through the registry is bit-identical to task 05 output.
