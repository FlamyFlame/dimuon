# Prompt: FCal Weighting Slides

Please read the reference document at
`~/Documents/physics/heavy-ion/dimuon/dimuon/Analysis/docs/fcal_scaling_pbpb.md`
and generate **three Beamer slides** (LaTeX) for a physics talk, as described below.

---

## Context

This is for a PbPb heavy-ion dimuon analysis at ATLAS (Run 3).
Centrality is determined from the forward calorimeter transverse energy (FCal ET).
PbPb 2024 and 2025 data do not yet have Glauber-calibrated FCal thresholds.
We reuse the 2023 thresholds, but must first align the 2024/25 FCal ET distributions
to the 2023 reference using a per-event FCal-ET-dependent weight.

Key numbers:
- 80% centrality boundary: FCal ET > 0.063208 TeV (= `FCal_ET_Bins_PbPb2023[79]`)
- Weight: `w(FCal) = h23_norm(FCal) / hyr_norm(FCal)`, capped at 5
- Stored per event as `fcal_corr_weight`; default 1.0 for 2023 (reference year)
- Applied as an event weight in all cross-section and RAA histograms

Plots (copy to your slide directory as needed):
- Before scaling: `plots/fcal_scaling/fcal_before_scaling.png`
- After scaling:  `plots/fcal_scaling/fcal_after_scaling.png`

The before plot shows two panels (2024 left, 2025 right) with raw FCal distributions
overlaid on 2023, normalised to area = 1 in the 0–80% region.

The after plot is a 2-column, 2-row layout:
- Upper row: weighted 2024/25 FCal distributions overlaid on 2023
- Lower row: the weight w(FCal ET) vs FCal ET for each year, with a dashed line at w=1

---

## Slide 1 — FCal ET Reweighting: Procedure

**Title**: FCal ET Reweighting for PbPb 2024/25 Centrality

Content (bullet points covering derivation and application):

**Motivation**
- PbPb 2024/25: no Glauber-calibrated FCal thresholds yet
- Reuse PbPb 2023 thresholds → must align FCal ET distributions

**Weight derivation** (`plotting_codes/fcal_scaling/derive_fcal_scaling.cxx`)
- Apply full 5-cut Run-3 event selection to raw skim (each year)
- Restrict to centrality 0–80%: FCal ET > 0.063208 TeV
- Normalise each year's histogram to area = 1
- Per-bin weight: $w(\text{FCal}) = h_{23}^\text{norm} / h_\text{yr}^\text{norm}$, capped at 5
- Save as `TGraph("g_fcal_weight")` in `fcal_weight_pbpb_20YY.root`

**Application in NTuple processing** (`NTupleProcessingCode/PbPbExtras.c`)
- For yr 24/25: load weight graph at initialisation
- Per event: $w_\text{FCal} = \texttt{g\_fcal\_weight.Eval(FCal\_Et\_raw)}$ if FCal > boundary; else 0
- Stored in `PairPbPbExtras::fcal_corr_weight` (float, default 1.0 for yr 23)
- Centrality recalculated with 2023 thresholds after weight is set

**Application in histogram filling** (`RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx`)
- `weight_for_RAA = T_AA_factor(centrality) × fcal_corr_weight`
- Applied to all cross-section and RAA histograms; trigger-efficiency histograms unweighted

---

## Slide 2 — FCal ET Distribution Before Weighting

**Title**: FCal ET Distribution Before Weighting

Include the plot `fcal_before_scaling.png`.

Caption / talking points:
- Two panels: left = PbPb 2024 vs 2023, right = PbPb 2025 vs 2023
- Distributions normalised to unit area over FCal ET > 0.063 TeV (0–80% centrality)
- Trigger bias from single-muon trigger (mu4) selects preferentially central events
- 2024/25 FCal distributions differ in shape from 2023 → direct use of 2023 thresholds
  would give inconsistent centrality fractions
- Reweighting corrects this shape difference

---

## Slide 3 — FCal ET Distribution After Weighting

**Title**: FCal ET Distribution After Weighting

Include the plot `fcal_after_scaling.png`.

Caption / talking points:
- 2-column layout: left = 2024, right = 2025
- Upper panels: weighted 2024/25 FCal distribution (red) overlaid on 2023 (black)
  — good agreement after reweighting
- Lower panels: per-event weight $w(\text{FCal ET})$ vs FCal ET
  — weight near 1 over most of the range; deviations at high FCal (few-event tails, capped at 5)
- Weights close to unity confirm the shape difference between years is small
- 2023 FCal thresholds can be applied to reweighted 2024/25 events self-consistently

---

## Format instructions

- Use standard Beamer with `\usetheme{Madrid}` or similar clean theme
- Physics notation: FCal $E_\mathrm{T}$, centrality percentile
- Keep each slide to ~5–6 bullet points; use sub-bullets sparingly
- For the plot slides (2 and 3): one `\includegraphics` filling ~70% of the frame width,
  with a short caption below and 2–3 key talking points as bullets alongside or below
- No \maketitle or table of contents needed — just the three frames
