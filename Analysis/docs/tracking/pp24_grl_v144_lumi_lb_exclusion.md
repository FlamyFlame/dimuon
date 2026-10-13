# pp24 GRL v144-pro42-02: luminosity 395.836 pb⁻¹ + run 488456 LB 84–220 exclusion

**Created:** 2026-10-12 · **Mode:** Implementation

## Objective
Move the pp24 (2mu4) normalisation to the new pp-ref GRL `SkimCode/xmls/physics_ppref2024_all_good.xml`
(`DetStatus-v144-pro42-02`), keeping numerator and denominator on the same events.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) pp24 NTuple processing drops run 488456 LB 84–220 (config header + `PPExtras` hook); (2) `PPBaseClass.h` pp24 2mu4 = 1/395.836; (3) `/review-analysis-code` PASS; (4) docs that quote 400.412 / 410.815 updated; (5) GRL rename `physics_2024ppRef_25ns.xml` → `…_legend_old_lumi_calib.xml` committed with `TrigRates_CA.py` repointed to `physics_ppref2024_all_good.xml`; (6) committed + pushed. The pp24 rerun (NTuple → hist filling → crossx) is the USER's.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure
1. **Motivation.** σ = N / (ε · L). L must be the luminosity of exactly the LBs that contribute to N.
2. **Equation.** pp24 crossx factor = 1/L_int [pb⁻¹], `PPBaseClass::CrossxFactorMap()[{24,"2mu4"}]`.
3. **Method.** The new GRL differs from the old (v132-pro42-01, used by the skim) only in run 488456: it starts at LB 221 instead of 84
   (verified, all other 7 runs identical LB ranges). New lumitable `IntNotes/data/luminosity/pp_2024/lumitable_pp_24_HLT_2mu4_L12MU3V.csv`:
   `Prescale Corrected` Total = 395.836 pb⁻¹ (old 400.412; run 488456 70.592 → 65.139). The skim has the old LBs, so the NTuple stage
   must DROP run 488456 LB 84–220 before any cut (same mechanism as `PbPbBadLBRanges`, `PbPbEventSelConfig.h`).
4. **Negative constraints.** Do not change any cut, binning or weight; no other run is touched. Numerator and denominator move together: never set 395.836 without the LB drop.

## Latest Stage
**2026-10-12 — DONE (review PASS).** Awaiting the user: cluster build of `PPExtras`, then the pp24 rerun (NTuple -> hist filling -> crossx/efficiency/plots).

## Progress Log

### 2026-10-12 — Implementation DONE
User rulings: (a) add the LB exclusion + 395.836; (b) rename old GRL to `physics_2024ppRef_25ns_legend_old_lumi_calib.xml`, repoint skim to the new GRL.
Files: NEW `NTupleProcessingCode/PPEventSelConfig.h` (`PPBadLBRanges(24)` = {488456, 84-220}); `PPExtras.{h,c}` (`evsel_bad_lbs_`, `PassEventSelExtra()`, picked up by `CallPassEventSel`, runs at `DimuonDataAlgCoreT.c:682` before any cut);
`PPBaseClass.h` 2mu4 = 1/395.836; `SkimCode/scripts/TrigRates_CA.py` pp2024 GRL -> `physics_ppref2024_all_good.xml`; docs (`pp_pipelines.md`, `plots_normalization_summary.md`, IntNotes `analysis_metadata.md`).
Tests: `PPEventSelConfig.h` standalone (LB 83/221 kept, 84/220 dropped, other runs/years untouched); GRL diff old vs new = only run 488456; lumitable run-488456 Good LBs 974 == new GRL. `/review-analysis-code` PASS it.1 (log `.claude/logs/review-analysis-code-*-pp24-grl-lb-exclusion.md`). PPExtras not compilable locally (old ROOT, C++20) -> cluster build pending.
**Flag (reviewer INFO, for the user / lumi group):** the new lumitable ALSO shifts the 7 unchanged runs by a uniform +0.266 % in `Prescale Corrected` (+0.877 pb^-1 total), i.e. a luminosity recalibration on top of the LB drop (-5.453 pb^-1); net 400.412 -> 395.836 (-1.14 %). Confirm the recalibration is intended and cite it.
Stale 400.412 remains in `IntNotes/tex/datasets.tex:70,82` (note-sync step) and in historical/closed tracking docs.
