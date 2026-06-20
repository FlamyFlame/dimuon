# Investigation: pp24 trigger-efficiency high-pT jump in crossx / R_AA

**Mode:** Investigation. **Created:** 2026-06-20. **Status:** ROOT CAUSE FOUND
(reviewed `/review-investigation` PASS iter 1, log
`.claude/logs/review-investigation-20260620-050718-pp-trig-eff-highpt-jump.md`).
Fix not yet applied — awaiting user decision on approach.

## Objective
Explain the large discontinuity / sudden jump UP in the pp24
trigger-efficiency-corrected cross-section (and R_AA, which drops/spikes at the
same pair pT) around pair pT ~50–60 GeV. Visible in
`plots/sanity_check_crossx/PP_2024_reco_eff_stages_pair_pt_in_eta.png` (red
"Reco×trig" curve shoots ~1–2 orders of magnitude above raw/reco-only at high pT,
huge error bars). No PbPb centrality shows it.

## Context
- pp crossx weight: `RDFBasedHistFillingPP.cxx:400-462` — 2mu4 PRODUCT weight
  `effcy_pair = ε₁·ε₂`, `w_trig = 1/effcy_pair`, **no floor/cap**.
- PbPb crossx weight: `RDFBasedHistFillingPbPb.cxx:933-985` — single-mu4 UNION
  weight `effcy_pair = ε₁+ε₂−ε₁ε₂`, `w_trig = 1/effcy_pair`.
- Single-muon ε lookup: `RDFBasedHistFillingData.cxx:596-625`
  (`EvaluateSingleMuonEffcyPtFitted`): `val = TF1::Eval(pt); if(val<0.01)val=0.01;`
- Fit files: `pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root`
  (pp); `pbpb_20YY/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root`
  (PbPb). Fits made separately (P2 trig-eff pipeline).

## Accumulated Findings (append-only)
1. **The turn-on TF1s return EXACTLY 0 above their fit range [4,60] GeV on
   read-back.** They are compiled-function TF1s (no analytic TFormula —
   `GetExpFormula()`==""). When such a TF1 is written to a ROOT file, only the
   sampled grid over [xmin,xmax] persists; `TF1::Eval(pt)` returns 0 for pt
   outside [4,60]. Verified pp TF1 `..._1_60_TO_2_00_divided`: Eval(60)=0.9893,
   Eval(62.5)=0.0000, Eval(70)=Eval(100)=0.000.
2. **Lookup floors ε to 0.01** (`RDFBasedHistFillingData.cxx:606`): the returned 0
   becomes 0.01. So every muon with pt>60 GeV gets ε=0.01 instead of its true
   ~0.97–0.99 plateau.
3. **pp 2mu4 product weight amplifies it, unfloored.** `w_trig = 1/(ε₁·ε₂)`: one
   >60 GeV muon → w_trig≈×100; both → ×10⁴. (Contrast w_reco, which IS floored at
   0.05 → ≤×20.) High pair_pt has tiny raw yield (0.015–1/bin) so a few
   mis-weighted pairs dominate → the jump.
4. **Data confirmation (pp).** trig×reco correction factor (corr/raw) per pair_pt
   bin from `histograms_real_pairs_pp_2024_2mu4_nominal.root`
   (`h2d_crossx_..._w_signal_cuts` / `..._no_trig_corr`, ProjectionX): smooth
   6.14→1.49 up to [48.7,58.3]; then **12.96 @[58.3,69.8], 73.7 @[69.8,83.6],
   80.8 @[83.6,100.2], 88.1 @[100.2,120]**. Onset bin straddles the 60 GeV
   single-muon fit edge.
5. **PbPb is immune via the union weight, NOT better fits.** PbPb fits ALSO return
   0 beyond [4,60] (Eval(70)=Eval(100)=0). But the union `ε₁+ε₂−ε₁ε₂` is bounded
   below by max(ε₁,ε₂): a floored leading muon (0.01) + an in-range subleading
   (~0.98) → denominator ≈0.98 → w_trig≈1.02. PbPb 2023 trig-corrected op yield
   falls smoothly 5762→…→1.86 with NO jump.

## Ruled Out (append-only)
- **Wrong file path / histogram name when applying the correction** — paths and
  hist names are correct (verified the code opens the right fit + hist files and
  reads the right histograms).
- **Obsolete/wrong 2D gap-region efficiency fallback** — not the cause; the gap
  2D fallback (`RDFBasedHistFillingData.cxx:612-619`) is a separate branch and the
  jump is driven by the fitted-TF1 branch returning 0 out of range.
- **Fit failure / wild fitted function** — fits are fine WITHIN [4,60]; the issue
  is read-back extrapolation returning 0, not a bad fit.
- **Breakdown of the no-correlation assumption (hypothesis 3)** — NOT the cause.
  The artifact is a sharp ×9 step localized exactly at the 60 GeV fit edge, not a
  smooth dR-correlation bias. mu4_mu4noL1 also requires both muons (product-type
  weight) using the SAME out-of-range fits → the jump would persist. The proposed
  full-chain mu4_mu4noL1 rerun is therefore unnecessary to explain THIS artifact
  (it can still be done later to study the correlation assumption itself, but it
  would not remove this jump).

## Recommended fix (follow-up — `/review-analysis-code`, NOT this investigation)
Make the lookup not propagate the out-of-range 0 as ε=0.01. Preferred (a):
**clamp pt into the TF1 range before Eval** (`pt → min(pt, xmax)`), so high-pT
muons get the plateau ε (~0.98) — physically correct since the trigger turn-on is
flat above threshold; fixes pp and PbPb consistently. Alternatives: (b) re-save
fits as analytic TFormula TF1s so Eval extrapolates the plateau; (c) floor/cap
w_trig (symptom only). After the fix: recompile RDF, refill pp (and PbPb) crossx,
replot sanity_check_crossx + crossx + R_AA, confirm the jump is gone.

## Fix applied & verified (2026-06-20) — DONE
7. **Fix (approach a) implemented and reviewed.**
   `RDFBasedHistFillingData.cxx::EvaluateSingleMuonEffcyPtFitted` (fitted-TF1
   branch): added `GetRange(xmin,xmax); pt_eval = clamp(pt, xmin, xmax);
   val = Eval(pt_eval)` before the existing `[0.01,1.0]` guard. So muons above the
   [4,60] turn-on range now get the plateau ε (~0.98) instead of 0→0.01. 2D gap
   fallback and the −1.0 return are unchanged; pp/PbPb weight formulas untouched;
   no cap added to w_trig (cause fixed, not symptom). Shared lookup → fixes pp and
   PbPb consistently. `/review-analysis-code` PASS iter 1 (log
   review-analysis-code-20260620-053030-trig-eff-plateau-clamp.md).
8. **Refilled** pp24 crossx + PbPb 23/24/25 (recompiled RDF, rc=0 each).
9. **Verification (corr/raw ratio per pair_pt bin, BEFORE → AFTER):**
   [58.3,69.8] 12.96→1.32; [69.8,83.6] 73.74→1.32; [83.6,100.2] 80.84→1.13;
   [100.2,120] 88.14→1.62. Low/mid-pt unchanged. Jump ELIMINATED. (Reviewer
   independently re-extracted: all MATCH.)
10. **Replotted & reviewed** (`/review-plot` PASS iter 1, log
    review-plot-20260620-054254-trig-clamp-crossx-raa-replot.md): PP stage plot
    red now tracks blue/black smoothly (no jump); pp24 crossx smooth falling
    spectrum; PbPb stage/crossx unchanged; R_AA vs pT/η/ctr smooth, no
    drop-to-zero discontinuity, PLACEHOLDER label present.

## Latest Stage
*(CLOSED 2026-06-20 — root cause found, fix applied (plateau-clamp in the shared
single-muon ε lookup), pp + PbPb refilled, all crossx/R_AA/sanity plots
regenerated and reviewed. pp high-pT trigger artifact eliminated; PbPb unchanged.
Changes uncommitted, awaiting user go-ahead to commit. Note: the underlying fits
remain compiled-function TF1s saved with only the [4,60] grid; if the fitting code
is ever changed to persist analytic TFormulas, the clamp becomes a harmless no-op
above the range. The mu4_mu4noL1 no-correlation study (roadmap, separate question)
is independent of this fix.)*
