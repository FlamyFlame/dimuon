# Plot Review Log
**Task**: Upfront hadronic/fake background reduction — |d0| discrimination (T3) + Δp/p distribution (T4) plot sets (doc E tf_upfront_bkg_reduction.md).
**Log file**: review-plot-20260706-200246-upfront-bkg-d0-dpop.md
**Started**: 2026-07-06T20:02:46-04:00
**Status**: APPROVED at iteration 2
**Iterations completed**: 2
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 1 WARNING + 4 INFO (0 CRITICAL)
**Details**:
1. [LAYOUT] WARNING — |d0| Mode-1 legend TLegend(0.16,0.68,0.48,0.90) overlaps the peak curves (mode1_shape_ptint_{signal,overlay} + ptbinned). Fix: relocate to a clear region (dpop top-right is the model).
2. [CONSISTENCY] INFO — "fake" color differs: d0 kGray+2 vs dpop kOrange+7.
3. [METHOD] INFO — Mode-1 norm source differs (d0 weighted, dpop unweighted); shape ~identical.
4. [REPRODUCIBILITY] INFO — dpop run_dpop.log shows PN/DPP_CUT redefinition (both .C define them) → same-session plotting failed; PNGs regenerated separately.
5. [AXES] INFO — dpop x-range [-0.3,0.6] vs spec ~[-0.6,0.6]; negative region empty, no physics lost.
**Physics (C1/C2/C3)**: ALL PASS — reviewer independently re-derived overlay shapes (hadronic weighted median 0.019 mm = smallest; real HF/prompt broad ~0.075) confirming no color/label bug. C3 RUN2-CROSSCHECK UNVERIFIED (neutral, per-muon provenance diagnostic, no Run2 analog).
**Numerical verification**: ALL MATCH (signal+overlay |d0| medians, real-muon 92% HF, Δp/p means, hadronic survival vs pT — every value reproduced from run_d0_final.log / run_dpop.log).
**Executor fixes applied (this iteration)**: (a) issue 1 WARNING — moved BOTH d0 Mode-1 & Mode-2 legends to (0.64,0.60,0.94,0.87) top-right + subtitles to (0.64,0.91), the empty region where d0 curves fall off at large |d0|. (b) issue 2 INFO — d0 fake color kGray+2(921)→kOrange+7(807) to match dpop. Re-ran d0 macro (all 8 PNGs). Issues 3/4/5 (INFO) left as documented; non-blocking. dpop set unchanged (its plots already PASS-clean).

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0 CRITICAL, 0 WARNING, 2 INFO (non-blocking)
- INFO: pt-binned plateau tails graze the top-right legend in a couple of busy panels (inherent to plateau shape; readable; optional bottom-left move) — no action.
- INFO: real HF/prompt render in slightly different red/blue shades between sets (both correct family; Δp/p set accepted iter 1) — no action.
**Numerical verification**: ALL MATCH (|d0| medians signal+overlay, Δp/p means, hadronic survival vs pT).
**Fixes confirmed**: |d0| Mode-1 & Mode-2 legends relocated top-right (no peak overlap); fake color harmonized (both sets orange). Physics C1/C2 PASS; C3 RUN2-CROSSCHECK UNVERIFIED (neutral).

**Status**: APPROVED at iteration 2
**Summary**: Both plot sets (|d0| discrimination T3 + Δp/p distribution T4) PASS — physics correct, numbers verified, labels/scales compliant, legend overlap fixed.
