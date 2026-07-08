# Plot Review Log
**Task**: Review regenerated d0/Δp/p/bkg_mc_provenance provenance study plots after the index-gate REVERT + new "real (untraced HIJING)" 5th class. HIJING muons ARE real.
**Log file**: review-plot-20260707-203156-index-gate-revert-5class.md
**Started**: 2026-07-08T00:31:56Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 3 (all INFO)
**Details**:
1. [LABEL/INFO] bkg make_plots.C unweighted stack labels real bin "prompt" (pre-existing mislabel; weighted variant already uses "real, truth-matched μ"). Out of scope for this revert; d0/dpop 5-way correctly use "real HF".
2. [COSMETIC/INFO] Minor legend/subtitle overlap on bkg survives_quality_cuts_weighted / provenance_stack_weighted SS panels — text readable, non-blocking (bkg plotter code unchanged by this revert, restored from backup).
3. [INFO] Signal pp24 draws the empty "real (untraced HIJING)" as a flat orange floor line + legend entry — JUDGED acceptable/honest (documents pp24 has no HIJING), no suppression needed.
**Numerical verification**: all MATCH — bkg overlay real/had/fake 96.57/3.05/0.38%; d0 real total 75659 == bkg single-mu real; pp24 untraced-HIJING=0; dpop untraced N=2031 > d0 1724 (no Δp/p cut); overlay real fraction high (HIJING=real).

**Status**: APPROVED at iteration 1
**Summary**: 5-class provenance plots (index-gate reverted, "real (untraced HIJING)" added; bkg restored) verified correct — 0 CRITICAL/0 WARNING. Physics sanity + all cross-checks pass.
