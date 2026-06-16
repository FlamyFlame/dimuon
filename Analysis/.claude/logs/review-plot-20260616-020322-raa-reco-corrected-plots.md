# Plot Review Log
**Task**: Reco-corrected R_AA plots (vs pair pT, pair eta, centrality), PbPb 23+24+25 combined / pp24, OS−SS, from modernized RAA_plotting.cxx.
**Log file**: review-plot-20260616-020322-raa-reco-corrected-plots.md
**Started**: 2026-06-16T06:03:22Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 3 WARNING
**Details**: (1) legend/info block overlapped data; (2) "_Centrality: 10-20%" stray leading underscore; (3) y-axis title "R_AA" never applied.
**Amend**: moved legend to compact upper-right (0.52,0.60,0.90,0.90) + TextSize 0.030; removed underscore in ctr_labels; set y-title "R_{AA}" + explicit x-title per mode on first-drawn hist; shortened long info lines (were clipping at right pad edge). Reran clean (0 warnings), 3 PNGs.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 1 WARNING
**Details**: eta panels — info/legend block still overprinted data because the eta R_AA fills the full vertical range (no clear top region), unlike pT/ctr.
**Amend**: increased y-headroom ymax = ylim*1.7 (was *1.1) so data occupies lower ~60% and the upper-right legend clears it (all modes). Reran clean; eta plot verified — legend now above data.

## Iteration 3
**Reviewer verdict**: FAIL
**Issues found**: 1 WARNING
**Details**: eta right subplot — a low-stat acceptance-edge bin (|eta|~2.4) had a giant error bar spiking into the upper-right info text.
**Amend**: (a) restrict eta display to [-2.3,2.3] (drop sparse edge bins); (b) move legend to upper-CENTRE (0.30-0.70, 0.64-0.90) so it clears the extreme-edge error-bar spikes (which sit at far left/right) while the 1.7x headroom keeps it above the data bulk for all modes. Verified eta + pT panels: no text/data overlap.

## Iteration 4
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: None found. Legend/info clear of data in all panels/subplots; y-title R_AA, x-titles correct; clean centrality labels; eta restricted to [-2.3,2.3]; R_AA>=0. R_AA>1 magnitude is the documented placeholder-normalization behavior.

**Status**: APPROVED at iteration 4
**Summary**: Reco-corrected R_AA plots (vs pair pT/eta/centrality, PbPb 23+24+25 combined / pp24, OS−SS) reviewed PASS after legend/axis/headroom fixes.
