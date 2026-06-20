# Analysis Code Review Log
**Task**: Fix |q·eta| mislabels (FIX A, comments/doc) + fullsim/overlay pair_eta->per-muon q*eta signal cut (FIX B)
**Log file**: review-analysis-code-20260619-235854-qeta-permuon-cut-fix.md
**Started**: 2026-06-20T03:58:54Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1 — executor work

### FIX A (comments/doc, no logic change)
- RDFBasedHistFilling/CommonEffcyConfig.h:15  "|q×eta| < 2.2" -> "-2.4 <= q*eta < 2.2 (one-sided)"
- RDFBasedHistFilling/CommonEffcyConfig.h:41  same
- RDFBasedHistFilling/CommonEffcyConfig.h:66  "per-muon |q×eta| < 2.2" -> "per-muon q*eta < 2.2 (one-sided, -2.4 <= q*eta < 2.2)"
- docs/tracking/raa_from_rdf_crossx.md §3a (line 48) "|q·η|<2.2 per muon" -> "−2.4≤q·η<2.2 per muon [one-sided, NO abs ...]"
- NOT touched (per task): RDFBasedHistFillingData.cxx:685 / .h:155 "|q*eta|<1.05" (legit symmetric barrel/endcap reco-eff split).

### FIX B (pair_eta -> per-muon q*eta signal cut)
Changed both truth ("truth_pair_eta < 2.2") and reco ("pair_eta < 2.2") to per-muon:
  truth: m1.truth_charge*m1.truth_eta < 2.2 && m2.truth_charge*m2.truth_eta < 2.2
  reco : m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2
in:
- RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsim.cxx:129 (truth), :131 (reco)
- RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsimOverlay.cxx:69 (truth), :71 (reco)
- RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsim.cxx:185 (truth), :186 (reco)
Matches reference data cut RDFBasedHistFillingPbPb.cxx:920 / PP.cxx:382 and truth refs PythiaTruth.cxx:412 / PowhegTruth.cxx:158.

### Flagged, NOT changed (reference/retired, per task)
- SingleBAnalysis/SingleBAnalysisBase.cxx:25  pair_eta < 2.2 (legacy, being retired)
- RDFBasedHistFilling/FillHistogramsCrossx_PP_clean.cxx:8  pair_eta < 2.2 (stale "clean" copy, not in active pipeline)

### Compilation (ACLiC, isolated sessions, source ~/setup.sh)
- RDFBasedHistFillingPythiaFullsim.cxx+        OK  -> (RDFBasedHistFillingPythiaFullsim) @0x...  (only -Wsign-compare warnings in ParamsSet.h, pre-existing)
- RDFBasedHistFillingPythiaFullsimOverlay.cxx+ OK  -> (RDFBasedHistFillingPythiaFullsimOverlay) @0x...
- RDFBasedHistFillingPowhegFullsim.cxx+        OK  -> (RDFBasedHistFillingPowhegFullsim) @0x...
No `error:` lines. Per-muon truth_charge/truth_eta and charge/eta members resolve against the struct. Samples NOT rerun (per task: compile + correctness only).

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: None found. FIX B reco/truth cuts character-for-character match references (PbPb.cxx:920/PP.cxx:382 reco; PythiaTruth.cxx:412/PowhegTruth.cxx:158 truth); one-sided <2.2, both muons, no abs; ternary & other terms intact. FIX A labels one-sided, binning arrays unchanged, |q*eta|<1.05 left untouched. Scope clean.
**Numerical verification**: No numbers to verify (cut-definition + comment change validated by compilation).

**Status**: APPROVED at iteration 1
**Summary**: Per-muon q·η<2.2 one-sided cut applied to all 3 fullsim/overlay fillers (was pair_eta); |q×eta| mislabels corrected in CommonEffcyConfig.h + tracking doc. All ACLiC-clean, reviewer PASS iter 1.
