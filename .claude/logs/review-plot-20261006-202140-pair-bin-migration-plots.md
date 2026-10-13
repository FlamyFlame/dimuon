# Plot Review Log
**Task**: Review NEW bin-migration plot set (plot_pp24_fullsim_pair_bin_migration.cxx, Tight + Medium) under pp24_reco_effcy_plots/<wp>/signal_region/all_os/bin_migration/. Physics Procedure = pair_bin_migration_bin_by_bin_feasibility.md §2–§3b.
**Log file**: review-plot-20261006-202140-pair-bin-migration-plots.md
**Started**: 2026-10-06T20:21:40-04:00
**Status**: APPROVED
**Iterations completed**: 3
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 9 (1 CRITICAL, 4 WARNING, 4 INFO)
**Details**:
1. [CRITICAL, label 16] bare "OS" in header on all 30 PNGs.
2. [WARNING, R1] cell-panel y-range and compare x-range driven by N_eff=2 cells; tails invisible on linear compare plot.
3. [WARNING] cells at exactly 0/1 drawn with zero error (normal approx).
4. [WARNING, label 17] 3D panel titles "9.0–12.8 GeV" lack variable name.
5. [WARNING, label 19] o_j undefined on plot.
6–8 [INFO] detached subscripts; "-2.2–-2.0" labels; "0.000" response text.
9 [INFO] report notes: ΔR top bin holds ~0.05 % of T; m bin 4 feed-in = J/ψ tail; doc "Medium within 0.002" wording.
**Numerical verification**: all MATCH; M/TR/R/Rout bookkeeping exact; C = p/(eps*s) holds.
**Amendments**: (1) "truth opposite-sign pairs"; (2) cell panels draw only cells with denominator N_eff >= 10 (kNeffDraw), y-range from drawn central values; compare plot -> CUMULATIVE N[T] fraction below x on log y (documented in code: tails 1e-4..1e-2); (3) 68% Wilson intervals with N_eff for cell points; (4) panel titles "9.0 < p_T^pair < 12.8 GeV", legend "1.080 < m_mumu < 1.535 GeV"; (5) legend "outside-range fraction = N_reco(truth out of range)/N_reco"; (6) dropped index subscripts (efficiency/stability/purity); (7) labels "[#minus2.2, #minus2.0)"; (8) text format .3g; (9) doc R1 wording updated. Tight + Medium regenerated.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 4 (1 WARNING, 3 INFO)
**Details**:
1. [WARNING R1] compare-plot y-floor 1e-5 hid m-cell F(0.5) ≈ 5e-6–9e-6 and ΔR stability F(0.6)=7.3e-6, contradicting the "nothing hidden" comment.
2. [INFO] dashed ΔR line gaps on steep rises. 3. [INFO] "N_{reco}" subscript gap. 4. [INFO] tails < 1e-4 come from N_eff ~1–12 cells (interpretation note).
**Numerical verification**: all MATCH (pT bin, ΔR top bin, 3D means, compare F(0.5) for ΔR 1.68e-4).
**Amendments**: y-floor = decade below the smallest nonzero F of either curve (now 1e-6); ΔR line style 7; legend N^{reco}. Regenerated Tight + Medium.

## Iteration 3
**Reviewer verdict**: PASS
**Issues found**: 2 INFO (N^{reco} kerning in 1D legend; mass legend 5 px from canvas edge) — left as is.
**Numerical verification**: compare-plot F(0.5)/min-F recomputed; all curves inside 1e-6..2 frame; values MATCH.

**Status**: APPROVED at iteration 3
**Summary**: Bin-migration plot set (30 PNG + 2 txt, Tight + Medium) approved after label (opposite-sign, physical panel ranges, o definition), statistics display (Wilson, N_eff>=10 draw cut) and cumulative-comparison fixes.
