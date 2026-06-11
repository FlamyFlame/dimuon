# Plot Review Log
**Task**: Reorganize reconstruction-efficiency plot output into per-category subdirectories (signed/single_b_op_compr/2d/distr for Pythia; 1d/single_b_op_compr/2d for Powheg) under each (dataset, mode, WP) reco_effcy_plots dir; ranged/ and det_resp dirs untouched. 3 files: PythiaFullsimRecoEffPlotter.cxx, PowhegFullsimRecoEffPlotter.cxx, plot_reco_distr_singleb_vs_op_pp24.C. Compile both plotters with ACLiC. Bulk move of existing files handled separately by user.
**Log file**: review-plot-20260610-210732-reco-effcy-subdir-reorg.md
**Started**: 2026-06-10T21:07:32-04:00
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**:
None found. Routing verified for both plotters (Pythia: signed/2d/2d/single_b_op_compr; Powheg: 1d/2d/2d/single_b_op_compr). mkdir before every changed SaveAs. ranged/ and det_resp dirs untouched. distr macro repointed to pp24_reco_effcy_plots/medium/distr/. Both plotters compiled with ACLiC, .so produced, no errors.
**Numerical verification**:
No numbers to verify.

**Status**: APPROVED at iteration 1
**Summary**: Reco-effcy plot output routed into per-category subdirs (signed/single_b_op_compr/2d/distr for Pythia; 1d/single_b_op_compr/2d for Powheg); ranged/ and det_resp untouched; both plotters compile.
