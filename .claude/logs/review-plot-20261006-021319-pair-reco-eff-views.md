# Plot Review Log
**Task**: Review pair reco-eff origin test + views (plot_pp24_fullsim_pair_reco_eff_views.cxx; docs/tracking/pair_reco_eff_origin_test_and_views.md)
**Log file**: review-plot-20261006-021319-pair-reco-eff-views.md
**Started**: 2026-10-06T02:13:19-04:00
**Status**: APPROVED at iteration 2
**Iterations completed**: 2
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 7 (1 WARNING, 6 INFO)
**Details**:
1. [STATISTICS] WARNING: origin-test significance overstated — driven by low-N_eff cells; with N_eff>=100 in both samples Tight 84.8/88 (p=0.58), Medium 78.4/88 (0.76); N_eff>=1000 Tight 37.7/35, Medium 32.7/35. Fig B leading pull partly eta-composition (3.9σ→2.9σ reweighted). Fix: report N_eff-thresholded tests; restate R1; sample choice = user's.
2. INFO: ratio point off-scale (2.0<eta<2.2, last pT bin, Neff~1-2); ratio values not in table.
3. INFO: C1 edge drops at kinematic ΔR_max explained (muon near 4.5 GeV threshold).
4. INFO: empty ΔR bins at high pT = kinematic bound.
5. INFO: applied/ flat; out of scope.
6. INFO: "single-b" not defined on plot (established label).
7. INFO: Wilson width at eps=1 symmetric (bar above 1).
**Numerical verification**: axes MATCH; inclusive eps MATCH; fig B chi2 MATCH (Tight+Medium); fig A MATCH; max pull MATCH.
**Amendments**: N_eff>=100 / >=1000 3D-cell chi2 added to table+stdout (reproduce reviewer exactly); ratio values written to table; R1 restated in tracking doc (no cell-level significance; user decision); views also generated for _op (all_os/).

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 2 INFO (ratio 1.0000+-0.0000 where complement empty — correct; R1 "±1 %" wording acceptable)
**Numerical verification**: N_eff-thresholded chi2 MATCH; R1 numbers MATCH; independent PyROOT spot checks MATCH (0.58024, 0.58687, 0.29182).

**Status**: APPROVED at iteration 2
**Summary**: Origin test + views (single_b and all_os, Tight+Medium) pass; origin-test statistics restated with N_eff thresholds (no cell-level origin dependence in populated cells).
