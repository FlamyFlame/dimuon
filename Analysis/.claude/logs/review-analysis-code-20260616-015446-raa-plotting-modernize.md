# Analysis Code Review Log
**Task**: Modernize RAA_plotting.cxx to read RDF crossx outputs (combined PbPb + pp24), OS-SS subtraction, reco-corrected; case 6, index-wise ratio, 15-bin pT, cluster paths.
**Log file**: review-analysis-code-20260616-015446-raa-plotting-modernize.md
**Started**: 2026-06-16T05:54:46Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 1 WARNING (+ 2 INFO)
**Details**:
1. [CORRECTNESS/WARNING] Mode 3 (RAA vs centrality) ratio loop off-by-one (ibin=0..N-1: hits underflow bin 0, misses last real bin N), no error propagation (error not scaled by 1/pp_int), no zero-guard. Legacy branch not rewritten with modes 1/2.
2. [STYLE/INFO] duplicate Clone name + un-deleted clones (macro leaks; acceptable).
3. [INFO] combined-year naive sum normalization — pre-existing, flagged, not-a-bug.
**Amend**: rewrote mode-3 loop to iterate real bins 1..N, propagate error (er=ea/pp_int), guard pp_int!=0; removed noisy cout. Reran clean (0 warnings/errors), 3 PNGs. INFO items not addressed (not required).

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0 CRITICAL/WARNING (4 INFO: hcrossx_pp_proj not nullptr-init [always assigned before use]; mode-3 drops pp-integral error [by design, pp high-stats]; leftover diagnostic cout; signed/unsigned compares — all harmless)
**Details**: Mode-3 fix verified (real bins 1..N, er=ea/pp_int exact for constant denom, zero+bounds guards). Modes 1/2 ratio correct. OS−SS for pp+PbPb, year-combine clone/Add/close correct, no double-free, T_AA not double-counted, sign convention correct.

**Status**: APPROVED at iteration 2
**Summary**: RAA_plotting.cxx modernized (RDF inputs, OS−SS, combined years, reco-corrected, index-wise ratio, mode-3 off-by-one fixed); reviewed PASS.
