# Analysis Code Review Log
**Task**: D7 numerator change (user 2026-10-13): reco-binned eps^reco_j = N[R_j]/N[T_j], R_j = matched && WP && reco cuts && truth in T (sum_i M_ij); feed-in excluded. Files: reco_binned_cells.h, build_pp24_fullsim_pair_reco_eff_reco_binned.C, fill_pp24_fullsim_pair_eta_smoothness_closure.C, fill_pp24_data_pair_eta_smoothness.C, plot_pp24_pair_eta_smoothness.cxx, plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx, Utilities/PairRecoEffRecoBinnedEvaluator.h. Constraint: crossx inputs unchanged.
**Log file**: review-analysis-code-20261013-164813-d7-numerator-truth-in-t.md
**Started**: 2026-10-13T16:48:13-04:00
**Status**: APPROVED
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 4 (all INFO)
**Details**:
1. [DOCUMENTATION] syst macro header: R_k not subset of T_k "because of feed-in" → now bin migration. INFO. FIXED (comment only).
2. [PHYSICS] data smoothness trig_reco_rb divides the raw data yield (still containing F_j) by the D7 eps → high by ~1/(1−o) ≈ +1.3 % (up to ~5 % top mass bin); common scale, shape conclusions unaffected; documented in tracking doc §2. INFO.
3. [CODE] fill_pp24_fullsim_pair_reco_sf_cells.C still fills R_sf incl. feed-in — consistent with ValueSF, not a crossx input; must change together with builder+ValueSF at propagation. INFO.
4. [NUMERICAL] eps = 1.000 ± 0.000 in two NeffT ≤ 2 cells (A/B/D gives zero error when A = all); pre-existing. INFO.
**Numerical verification**: all MATCH (inclusive 0.584097 Tight = TR/T exactly; Medium 0.635128; N[F] 0.0528358, 1.3245 %; readback 0.565195; eps>1 3 cells; syst central +0.0455(32); SF product / placeholder mtimes unchanged 02:00/02:04; product marker = new Value()).

**Status**: APPROVED at iteration 1
**Summary**: D7 numerator (truth in T, Σ_i M_ij) implemented consistently in the cell formula, builder, closure filler and syst inputs; crossx SF path byte-identical and rebuild-guarded.
