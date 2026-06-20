# Plot Review Log
**Task**: Revert the PbPb cross-section plotter (plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx) from "differential_crossx" back to "TAA_weighted", undoing the plotter-only hunks of commits f8e4697 and 9c44533. Recompile (ACLiC) and rerun the plotter to regenerate the TAA_weighted/ plots, removing the now-orphaned differential_crossx/ dirs. Physics requirement: the TAA-weighted plots must be the same quantity that goes into R_AA.
**Log file**: review-plot-20260618-150015-pbpb-crossx-revert-taa-weighted.md
**Started**: 2026-06-18T15:00:15-04:00
**Status**: APPROVED
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1 — Executor work
**Code change**: Reverted `plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx`
to its pre-`f8e4697` state via `git checkout f8e4697^ -- <file>`. Result is
byte-identical to `f384c42` (verified `git diff --stat` empty). The working-tree
diff vs HEAD is exactly the reverse of the plotter hunks of f8e4697 + 9c44533:
TAA_weighted/ dir restored, _w_signal_cuts_ histos read, "[pb]" labels restored,
differential_crossx/ + _dsigma histos removed (both pt_120 main block and pt_150 block).
**Not touched** (intentional, out of scope of the plotting revert): the RDF
dsigma additions and the unrelated MC-data reco change in RDFBasedHistFillingPP.cxx
from f8e4697 — only the plotter was reverted.
**Run**: `root -l -b -q 'plot_single_b_crossx_pbpb.cxx(true)'` (sourced setup.sh).
Produced TAA_weighted/ in:
  - pbpb_23_24_25_combined/TAA_weighted/ — 30 PNGs (6 ctr × 5 plots), none 0-byte
  - pbpb_23_24_25_combined_pt_150/TAA_weighted/ — 12 PNGs (6 ctr × 2 plots), none 0-byte
**Cleanup blocked**: `rm -rf` of the orphaned differential_crossx/ dirs was DENIED
by the permission policy (twice). Left in place; will flag to user.
**Physics consistency (point 1)**: TAA plots read _w_signal_cuts_ histos filled with
`weight_for_RAA_trig_corr` (RDFBasedHistFillingPbPb.cxx:1061-1100), identical weight
to the RAA-input histo h3d_op_crossx_w_signal_cuts_vs_centr_... (line 1006-1009) that
RAA_plotting.cxx:31 consumes. Same quantity → consistent with R_AA.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 2 (both INFO)
**Details**:
1. [CLEANUP] INFO — orphaned differential_crossx/ dirs remain (rm -rf permission-blocked); manual removal when convenient; does not affect TAA_weighted output.
2. [LABEL] INFO — 2D/subplot labels read "[pb]" though quantity is T_AA-weighted R_AA-input yield; inherited verbatim from f384c42 (byte-identical), not a regression, out of scope.
**Numerical verification**: ALL MATCH — git diff vs f384c42 EMPTY; 30+12 PNGs none 0-byte; plotter reads _w_signal_cuts_ (no _dsigma); per-centrality _w_signal_cuts_ histos AND RAA-input 3D histo both filled with weight_for_RAA_trig_corr (= weight_for_RAA·w_reco·w_trig); RAA_plotting.cxx:31 consumes the matching histo.

**Status**: APPROVED at iteration 1
**Summary**: PbPb crossx plotter reverted byte-identical to the prior TAA_weighted version (f384c42); TAA_weighted plots regenerated (30 pt_120 + 12 pt_150). TAA-weighted plotted quantity uses the identical weight as the R_AA input — consistent with R_AA.
