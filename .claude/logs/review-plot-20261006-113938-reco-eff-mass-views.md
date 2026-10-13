# Plot Review Log
**Task**: mass-axis pair reco-eff views (all OS), step 8 of pair_reco_eff_origin_test_and_views.md
**Log file**: review-plot-20261006-113938-reco-eff-mass-views.md
**Started**: 2026-10-06T11:39:38-04:00
**Status**: APPROVED at iteration 4
**Iterations completed**: 4
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 3 (1 WARNING, 2 INFO)
**Details**:
1. WARNING [PHYSICS-INTERPRETATION C1/C2] R2 "spread ≲3 %" understated: real 3.8–5.1 %; top-mass-bin dip ≈4 % (~10σ, all pT bins, both WPs) and forward (2.0–2.2) mass slope −17/−18 % unnamed. Fix: restate with numbers; name the dip (reco-mass migration out of window, hypothesis) and the forward slope; confirm cheaply or mark hypothesis.
2. INFO pT labels %.4g (consistent with approved views). 3. INFO 2.0–2.2 / 105–150 GeV panel low stats (covered by errors).
**Numerical verification**: all MATCH (edges, axes, under/overflow, Tight pT1 eta-integrated / |eta|<1 / 2.0–2.2 values, Medium spot-check).
**Amendments**: R2 rewritten with the real spreads; top-bin dip quantified against an upper bound from h_minv_zoomin_vs_truth_minv_zoomin_single_b_pass_tight (≲12 % top, ≲1 % bottom, observed 4 % / 1 %) → consistent with reco-mass migration; forward slope stated as HYPOTHESIS. Plots unchanged (numbers were correct).

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 5 (1 WARNING, 4 INFO)
**Details**:
1. WARNING lower-edge "≲1 %" labelled as bound but is not (reco<1.05 undercounts); Gaussian-smearing model gives σ≈0.020 GeV (1.1 GeV), ≈0.040–0.045 (2.9 GeV) → top loss ≈3.5–3.9 % vs 4 % observed, bottom ≈1.8 %.
2. INFO significance ≈8–9σ, six populated pT bins (not all). 3. INFO response matrix is single-b; note origin-independence of resolution. 4. INFO INDEX "nearly flat" wording. 5. INFO forward slope includes the 4 % migration; first→third ≈6.6 %.
**Numerical verification**: R2 numbers all MATCH; migration fractions reproduced; 12 % upper bound valid.
**Amendments**: R2 dip bullet rewritten (upper bound kept for the top edge only; model estimate for both edges, labelled as model; 8–9σ, six pT bins; single-b response note); forward-slope bullet notes the migration share (first→third ≈7 %); INDEX wording updated.

## Iteration 3
**Reviewer verdict**: FAIL
**Issues found**: 4 (1 WARNING, 3 INFO)
**Details**:
1. WARNING "(σ_m/m ≈ const)" contradicts σ_m 0.040–0.045 at 2.9 (const would give 0.054); value comes from the top-edge rows; exact uniform-bin inversion 0.041–0.049 → loss 3.6–4.3 %.
2. INFO name both bin widths (0.075 response, 0.455 analysis). 3. INFO "rigorous" → "upper bound (uniform truth density)". 4. INFO INDEX wording ambiguity re forward slope.
**Numerical verification**: σ inversion 0.0203 ✓; losses 3.51–3.95 %, 1.78 % ✓; bound 11.8 % ✓; slopes ✓.
**Amendments**: R2 model bullet rewritten (σ source = top-edge rows, exact inversion 0.041–0.049, σ_m/m 1.8 % → 1.4–1.7 % not const, both widths named, loss 3.6–4.3 %); bound bullet relabelled; INDEX wording clarified.

## Iteration 4
**Reviewer verdict**: PASS
**Issues found**: 3 INFO (INDEX clause separation; "match" wording for bottom bin; show bound formula) — all applied after the verdict.
**Numerical verification**: σ inversions 0.041–0.049 / 0.0203 reproduced numerically; losses 3.6–4.3 % / 1.75 %; bound 11.8 %; slopes ✓.

**Status**: APPROVED at iteration 4
**Summary**: Mass-axis views correct; R2 interpretation now quantitative and honestly labelled (top-bin dip = reco-mass migration, model 3.6–4.3 % vs 4 % observed; forward slope = hypothesis).
