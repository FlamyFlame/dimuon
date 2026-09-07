# Plot Review Log
**Task**: Verify (and fix if needed) the linear-y-scale remake of DrawPairEtaIntegratedTwoSeries in Analysis/plotting_codes/single_b_analysis/SingleBCrossxPlotterBase.cxx. Change ONLY the pair-eta-dependence diagnostic plot (pp24 + PbPb 23+24+25 combined) from log-y to linear-y per the newly adopted convention (non-momentum, non-log-binned variables default to linear y — .claude/commands/review-plot.md R4). Specific concern: fixed-NDC legend box (tuned for the old log-scale ceil_pad=6.0 headroom) may now overlap the tallest data points under the new linear max*1.35 headroom, especially for PbPb whose pair-eta shape peaks sharply near eta=0.
**Log file**: review-plot-20260907-223044-mupt45-eta-linear-yscale.md
**Started**: 2026-09-07T22:30:44Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 2
**Details**:
1. [STALE OUTPUT / FRESHNESS] CRITICAL — reviewed PNGs were rendered before a later
   source edit (LeftMargin 0.13->0.17, and an (invalid, non-compiling) SetLabelFormat
   attempt) landed; reviewer could not see the current code's effect. Superseded by
   events below (SetLabelFormat did not compile; fixed via TGaxis::SetMaxDigits instead,
   see Step 6 in the tracking doc).
2. [LEGEND / FRAME OVERLAP] CRITICAL — fixed-NDC TLegend(0.15,0.78,0.97,0.94) extended
   above the frame's top border (default top margin 0.9), pixel-confirmed on both PNGs;
   pre-existing, independent of log/linear.

**Numerical verification**: all 4 counts + both percentages MATCH (pp24 705404/541427,
PbPb 275622/180446).

## Amendment (out-of-loop, user directive)
User sent "no review loop needed for these" before this reviewer's verdict arrived
(it had been running since before that message). Per that instruction, no further
reviewer subagent was spawned. Issue 2 (the only substantive one; Issue 1 was a
staleness artifact of ongoing work, not a defect) was fixed directly by the
orchestrator: headroom raised to max*1.4, legend moved to (0.15,0.73,0.97,0.90) —
inside the frame, clear of the tallest point. Recompiled, reran both driver macros,
visually re-verified both PNGs (no overlap, consistent x10^p headers via
TGaxis::SetMaxDigits(3)). See docs/tracking/muon_pt45_cut_diagnostic.md Step 6 for
the full record.

**Status**: CLOSED WITHOUT FORMAL RE-REVIEW (user directive) at 2026-09-07T22:52:00Z.
Fix applied and independently self-verified by the orchestrator (pixel/visual
inspection of regenerated PNGs); not re-submitted to a reviewer subagent.

## Execution (done before this log was opened)
- Edited `Analysis/plotting_codes/single_b_analysis/SingleBCrossxPlotterBase.cxx`,
  `DrawPairEtaIntegratedTwoSeries`: removed `c.SetLogy()`; replaced the log-only
  `ApplyCommonLogYRange({hpa,hpb}, ceil_pad=6.0)` headroom call with a linear
  computation (`eta2s_ymax` = max bin content over hpa/hpb; `SetMinimum(0.0)`,
  `SetMaximum(eta2s_ymax*1.35)`), with an inline comment explaining the switch and
  why the log-specific utility no longer applies.
- No other method touched (`DrawPairPtIntegratedTwoSeries`, `DrawPairPtByEtaTwoSeries`
  keep `SetLogx()+SetLogy()` — pair pT is log-binned/power-law, unaffected by this
  change).
- Recompiled via ACLiC (`root -l -b -q 'plot_pp24_muon_pt45_diagnostic.cxx+()'` and
  `'plot_pbpb_muon_pt45_diagnostic.cxx+()'`) and reran both driver macros.
- Regenerated PNGs (all 3 per dataset regenerate on each run; only
  `pair_eta_dependence_mupt_45_vs_40.png` actually changed pixel content):
  - `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp24/muon_pt45_diagnostic/pair_eta_dependence_mupt_45_vs_40.png`
  - `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pbpb_23_24_25_combined/muon_pt45_diagnostic/pair_eta_dependence_mupt_45_vs_40.png`
- Counts unchanged from the aae1445 commit (confirms no selection/histogram logic
  touched, only the draw/axis code): pp24 705404 (pT>4) -> 541427 (pT>4.5) pairs
  (-23.2458%); PbPb 23+24+25 combined 275622 -> 180446 pairs (-34.5314%). Printed by
  the macros themselves (`muon_pt45_diagnostic_summary.txt` in each output dir,
  byte-identical to the committed version).
