# Plot Review Log
**Task**: CERTIFY final |d0|+Δp/p set (index-gate reclassification + per-event MC weights on both modes + recolor + combined 3-way + new sample labels/filenames + Δp/p cut-line label). Must PASS with C1-C6 incl. C5 provenance + C6 MC-weight.
**Log file**: review-plot-20260707-120911-upfront-d0-dpop-certify.md
**Started**: 2026-07-07T12:09:11-04:00
**Status**: APPROVED at iteration 2
**Iterations completed**: 2
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL (3 WARNING + 2 INFO; ALL 32 plots certified correct/presentation-grade)
- WARNING 1-3: SUMMARY.md files stale (pre-index-gate overlay numbers + wrong "hadronic smallest |d0|" claim). d0 overlay real-prompt 7891→6167, hadronic 2390→4114, total 75659→73935, 91.7% HF, 0 unresolved; real-prompt (0.017) is smallest, not hadronic (0.032). dpop overlay realprompt 8671→6640 (91.8% HF), hadronic→7685.
- INFO 4-5: d0 top subtitle grazes top frame; dpop ptbinned sample label leading char grazes pad-left. (non-blocking)
- PASS on plots: C5 index-gate correct, C6 MC-weights correct (both modes weighted), recolor, combined, labels, cut-line, physics, all numbers-vs-logs match.
**Executor fixes**: refreshed BOTH SUMMARYs to the v8/v3 log numbers + corrected the "smallest |d0|" claim (real-prompt smallest). Fixed the Δp/p cut label per user (was buried low in the THStack → moved to the top strip, rotated, above the stack, reads in full). dpop re-plotted. Plots otherwise unchanged (already certified).

## Iteration 2
**Reviewer verdict**: PASS (0 CRITICAL, 0 WARNING; 2 non-blocking INFO: d0 top subtitle grazes frame, dpop 2x2 sample-label leading char grazes pad-left — both legible)
**Verified**: both SUMMARYs now match run_d0_v8.log / run_dpop_v3.log exactly (d0 overlay real-prompt 6167/0.017 smallest, hadronic 4114/0.032, 91.7% HF, 0 unresolved; dpop 91.8% HF, hadronic 7685); "hadronic smallest" claim corrected to real-prompt; Δp/p cut label readable at top of all dpop panels; no regression; C1-C6 all PASS.
**Status**: APPROVED at iteration 2
**Summary**: Final |d0|+Δp/p set CERTIFIED presentation-ready — index-gate provenance, per-event MC weights on both modes, recolor, combined 3-way, descriptive sample labels/filenames, Δp/p cut label fixed.
