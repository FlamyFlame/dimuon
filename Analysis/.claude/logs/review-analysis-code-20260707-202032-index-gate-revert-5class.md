# Analysis Code Review Log
**Task**: Review the index-gate REVERT + new "real (untraced HIJING)" 5th class across the d0/Δp/p
provenance macros, and the restore-from-backup of the bkg_mc_provenance classifiers. Physics-correctness
against the authoritative §3 provenance definition (HIJING muons ARE real; hadronic = prob>0.5 &
(|id|!=13 || IsPrimary==0); NO index/barcode gate on real-vs-hadronic).
**Log file**: review-analysis-code-20260707-202032-index-gate-revert-5class.md
**Started**: 2026-07-08T00:20:32Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 2 (both INFO)
**Details**:
1. [DOCS/INFO] Stale "4-way"/"p=0..3" header/section comments (code correct, NC=5/NPROV=5). → FIXED post-review (d0 lines 3/67/250/292 → "5-way"; dpop_dist lines 3/22 → "5-way"/"p=0..4").
2. [INFO] plot_dpop_dist mode-1 normalizes the weighted "w" histos while the dpop_dist header describes mode-1 as unweighted — pre-existing, physically defensible (unit-area shape), not introduced by this change. No action.
**Numerical verification**: N/A (code-logic review only).
**What passed**: index gate removed from real-vs-hadronic in both classifiers; parent-map [0,lim) cutoff kept for HF/prompt sub-split ONLY (never gates real); 5-class reorder coherent everywhere (incl. dpop hadronic-survival table → index 3, combined "real (all)" = 0+1+2, array sizes 4→5/NPROV); untraced routing no double-count; MC AMI weights preserved (C6); bkg_mc_provenance.C + fill_weighted_fullsim.C carry NO index-gate remnant (C5 provenance OK).

**Status**: APPROVED at iteration 1
**Summary**: Index-gate revert + "real (untraced HIJING)" 5th class verified physically correct against §3; restored bkg classifiers clean. Zero CRITICAL/WARNING.
