# Muon reconstructed pT > 4.5 GeV cut — statistics-vs-systematics diagnostic

**Type:** Implementation (diagnostic plot set + temporary code).

## Objective

The muon q·η spectra (`muon_gap_cuts_acceptance.md` F14,
`plot_muon_q_eta_pt_dependence.cxx`) show that muons with reconstructed
pT < 4.5 GeV have a markedly different q·η shape than higher-pT muons,
indicating higher sensitivity of their reconstruction efficiency to detector
gaps/edges. We are considering raising the muon reconstructed-pT cut from the
current 4 GeV to 4.5 GeV. Before that decision, we need to see the impact on
single-b signal statistics (RAW pair counts, no efficiency corrections),
since the power-law pT spectrum means a small pT-threshold change can cost a
disproportionate fraction of pairs.

**Done =**
1. Three PNGs for **pp24 data** in
   `plots/single_b_analysis/pp24/muon_pt45_diagnostic/`:
   (a) raw OS signal-region pair counts vs pair η (pair pT integrated),
   (b) vs pair pT (pair η integrated), (c) vs pair pT in the 9 pair-η
   panels (subplots) — each overlaying two series/markers: muon pT>4 GeV
   (current) and muon pT>4.5 GeV (candidate).
2. The same three PNGs for **Pb+Pb 2023+2024+2025 combined data** in
   `plots/single_b_analysis/pbpb_23_24_25_combined/muon_pt45_diagnostic/`.
3. Total percentage decrease in signal-region raw pair counts going from the
   4 GeV to the 4.5 GeV cut, reported separately for pp24 and for Pb+Pb
   23+24+25 combined.
4. All axes use the SAME pair-pT / pair-η binning as the crossx plots
   (`ParamsSet::pT_bins_120`, `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`)
   — no new binning invented (`.claude/CLAUDE.md` §Binnings).
5. Code clearly marks the 4.5 GeV filter as TEMPORARY/diagnostic, additive
   (does not touch the nominal 4 GeV signal selection or any existing
   histogram), and states where the permanent change would live if adopted.

**Stop-and-ask** = ANY physics-results-bending ambiguity (use judgment; if
unsure whether it's blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

1. **Motivation.** §Objective above. The reconstruction-efficiency systematic
   from detector gaps/edges is worse for pT<4.5 GeV muons; raising the cut
   trades signal statistics for a cleaner (less gap-sensitive) efficiency
   correction. The final call is the user's; this doc only quantifies the
   statistics side.
2. **Top-level quantities.** For each dataset (pp24; PbPb 23+24+25 combined),
   for each of the two muon-pT thresholds `pT_cut in {4.0, 4.5} GeV`:
   `N(pT_cut) = sum over OS pairs passing [current signal region] AND
   [muon1.pt > pT_cut] AND [muon2.pt > pT_cut]`, unweighted (raw count, no
   trigger/reco efficiency correction, no luminosity scaling, no background
   subtraction). `% decrease = 100 * (N(4.0) - N(4.5)) / N(4.0)`.
3. **Current signal region used as-is per dataset** (they currently differ —
   `signal_selection_change_impact.md` §0 — this diagnostic does NOT
   harmonize them, only asks "what happens to whichever selection is
   currently nominal for each dataset"):
   - **pp24:** `minv in (1.08,2.9) && pair_pt>8 && both muons pass
     ParamsSet::single_mu_fiducial_gap_cuts` (detector-gap fiducial cut),
     Tight WP, OS (`RDFBasedHistFillingPP.cxx:584-587,600`).
   - **PbPb (all 3 years):** `minv in (1.08,2.9) && pair_pt>8 &&
     m1.charge*m1.eta<2.2 && m2.charge*m2.eta<2.2` (older per-muon cut, not
     yet migrated to the fiducial-gap cut), Tight WP, OS
     (`RDFBasedHistFillingPbPb.cxx:977`).
   - Both already require reconstructed muon pT>4 GeV, but that cut is
     baked into **NTuple processing**, not the RDF signal region:
     `NTupleProcessingCode/DimuonDataAlgCoreT.c:596`
     (`if (m1.pt < 4 || m2.pt < 4) return false;`), shared by pp and PbPb
     data ntuple processing.
4. **Negative constraints.**
   - This is a pure DATA raw-count study. No MC, no AMI weights, no
     trigger/reco efficiency correction, no isospin weight — none of those
     are touched or relevant here.
   - Do NOT change, widen, or re-derive the minv / pair-pT / gap-cut parts
     of either dataset's current signal region. Only the muon-pT threshold
     is varied, additively, on top of the existing selection.
   - Do NOT retype the pair-pT or pair-η bin edges anywhere; read them from
     `ParamsSet::pT_bins_120` / `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`
     via the existing `SingleBCrossxPlotterBase` machinery
     (`.claude/CLAUDE.md` §Binnings).

## Design Decisions

- **Where the pT>4.5 filter is applied (temporary):** additively in the RDF
  hist-filling stage (`RDFBasedHistFillingPP.cxx` / `RDFBasedHistFillingPbPb.cxx`),
  as `.Filter("m1.pt > 4.5 && m2.pt > 4.5")` on top of the existing
  `df_single_b_crossx` / `df_crossx_ctr` node, producing ONE new,
  clearly-named, clearly-commented histogram per dataset
  (`..._diag_mupt45` suffix). The pre-existing pT>4 counts histogram
  (`h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts` for pp24;
  `h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>_counts` per
  centrality bin for PbPb, already `weight`=1 for data) is read UNCHANGED
  as the "pT>4 (current)" series. **Why RDF stage, not NTuple processing:**
  user's explicit instruction, for speed at this diagnosis stage — avoids a
  Condor rerun of NTuple processing. **If the 4.5 GeV cut is adopted:** the
  permanent change goes in `NTupleProcessingCode/DimuonDataAlgCoreT.c:596`
  (data) and the analogous truth-pT cuts in
  `NTupleProcessingCode/PythiaFullSimExtras.c:404` /
  `NTupleProcessingCode/PowhegFullSimExtras.c:167` (MC), ALL data and MC
  results must be rerun (full Condor + RDF + plotting chain), and the
  temporary `_diag_mupt45` blocks added here must be DELETED.
- **PbPb centrality handling:** PbPb has no centrality-inclusive counts
  histogram on disk; the new plotting code sums the per-centrality-bin 2D
  histograms (both the existing pT>4 counts and the new pT>4.5 diagnostic
  counts) the same way `SingleBCrossxPlotterPbPbCombined::GetHistObject`
  already sums "_counts"-suffixed histograms across the 3 years (simple sum,
  not luminosity-weighted — a raw count is a raw count).
- **Plotting:** reuse `SingleBCrossxPlotterBase` (both pp24 and PbPb-combined
  plotters derive from it) — add three new two-series overlay Draw methods
  (pair-η integrated 1D, pair-pT integrated 1D, pair-pT-in-pair-η-panels)
  modeled on the existing `DrawPairPtByEtaWithDrLines` (N-series overlay
  pattern) and `DrawPairPtByEta` (panel projection logic), rather than
  duplicating panel/binning code.

## Implementation Plan

1. RDF: add temporary `_diag_mupt45` counts histogram to
   `RDFBasedHistFillingPP.cxx` (per §Design Decisions) — per §3 pp24. →
   `/review-analysis-code`.
2. RDF: add temporary `_diag_mupt45` counts histogram to
   `RDFBasedHistFillingPbPb.cxx` (per centrality bin) — per §3 PbPb. →
   `/review-analysis-code`.
3. Recompile (ACLiC `.L ...cxx+`) and rerun crossx hist filling for pp24 and
   PbPb 23/24/25 (existing `run_crossx_hist_filling_*.sh` scripts — NTuple
   processing untouched, no Condor needed).
4. Plotting: add 3 new two-series overlay Draw methods to
   `SingleBCrossxPlotterBase.cxx` + two small driver macros (pp24, PbPb
   combined) producing the 6 PNGs + a CSV with the percentage decrease. →
   `/review-plot`.
5. Run the driver macros, verify outputs, compute and record the
   percentage decrease for both datasets.

## Progress Log

**Step 1 (pp24), DONE.** Added a temporary additive `_diag_mupt45` block to
`RDFBasedHistFillingPP.cxx` right after the existing
`h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts` booking: filters
`df_single_b_crossx_weighted` with `m1.pt > 4.5 && m2.pt > 4.5` and books an
unweighted `h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45` on
the same axes. ACLiC-compiled clean, ran via
`RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh`. `/review-analysis-code`
PASS at iteration 1 (log: `.claude/logs/review-analysis-code-20260906-221348-mupt45-diag-hists.md`).
**Result: pp24 705404 (pT>4) -> 541427 (pT>4.5) pairs.**

**Step 2 (PbPb), DONE but redesigned mid-step.** First tried the same additive
in-class approach in `RDFBasedHistFillingPbPb.cxx` (both from the post-weight
`df_crossx_ctr` and, after that failed, from the pre-weight `df_single_b_crossx`).
Both hit a PRE-EXISTING, unrelated bug: `EvaluateSingleMuonEffcyPtFitted: no
fitted turn-on ... on the CONTIGUOUS coarse q*eta binning` thrown mid-event-loop
and silently swallowed by ROOT (exits 0, near-empty output) -- the already-
tracked ACTIVE issue in `pp_trig_eff_highpt_jump.md` (PbPb has not migrated to
the pp24 fiducial-gap cut, so some muons fall outside the fitted turn-on
range). Because RDF's RunGraphs shares ONE event loop across every result tied
to a source tree, the throw kills ALL histograms from that loop, including
ones that never read the poisoned trigger-efficiency column -- so no in-class
additive block could work. **Asked the user** (AskUserQuestion); approved
"bypass it, but be careful to ensure the counts results are correct, and mark
the bug as a future TO-DO in tracking doc." `RDFBasedHistFillingPbPb.cxx` was
reverted to pristine (`git diff` clean, verified). Instead wrote a NEW
standalone macro `RDFBasedHistFilling/fill_pbpb_muon_pt45_diag_counts.cxx`
that opens a fresh `TChain("muon_pair_tree_sign2")` (OS) directly over the 3
years' `muon_pairs_pbpb_20YY_single_mu4_mindR_0_02.root` ntuple-processing
output, entirely decoupled from `RDFBasedHistFillingPbPb`'s poisoned trigger
machinery, and mirrors PbPb's EXACT current signal region + Tight WP verbatim
from `RDFBasedHistFillingPbPb.cxx:977,987`
(`pair_pass_tight && minv>1.08 && minv<2.9 && pair_pt>8 && m1.charge*m1.eta<2.2
&& m2.charge*m2.eta<2.2`), booking two unweighted counts histograms (pT>4
baseline, pT>4.5 diagnostic) on the SAME axes as pp24's
(`ParamsSet::pT_bins_120` x `ParamsSet::N_PAIR_ETA_CROSSX_BINS`). Combining
the 3 years is a plain TChain union -- correct for an unweighted count (no
luminosity weighting needed), matching the existing simple-sum-for-`_counts`
rule in `SingleBCrossxPlotterPbPbCombined::GetHistObject`. ACLiC-compiled
clean, ran successfully with no errors. `/review-analysis-code` PASS at
iteration 1 (same log, covered both pp24 and PbPb files together).
**Result: PbPb 23+24+25 combined 275622 (pT>4) -> 180446 (pT>4.5) pairs.**

**Step 4 (plotting), DONE.** Added three new two-series-overlay Draw methods
to `SingleBCrossxPlotterBase.cxx` (`DrawPairPtByEtaTwoSeries`,
`DrawPairEtaIntegratedTwoSeries`, `DrawPairPtIntegratedTwoSeries`), modeled on
the existing `DrawPairPtByEtaWithDrLines`/`DrawPairPtByEta` panel-projection
logic -- no panel/binning code duplicated. Two new driver macros
(`plot_pp24_muon_pt45_diagnostic.cxx`, `plot_pbpb_muon_pt45_diagnostic.cxx`)
produce the 6 requested PNGs + 2 summary txt files (percentage decrease).
**Self-review caught and fixed 3 real plotting defects before formal review**
(all verified by rerendering against the real ROOT output, not guessed):
(1) legend text overlapping data markers in the 9-panel plot (fixed:
geometry now mirrors the proven `DrawPairPtByEtaWithDrLines` two-box layout,
plus shortened in-panel labels -- the full labels stay on the two single-panel
PNGs); (2) legend text clipping past the frame edge on both single-panel PNGs
(fixed: wider boxes + smaller font + `ceil_pad` headroom, per
`Utilities/CommonLogYRange.h`); (3) a reproducible ROOT rendering bug: a
`TLegend` "lep"-option entry referencing a REAL histogram object, drawn
alongside a text-only `TLegend`/entries on a log-y pad, rendered as a
full-frame-height spike instead of a small icon (isolated with a minimal
standalone repro against the actual PbPb output; confirmed via independent
bin-content/error dump that the underlying data was always clean -- this was
purely a rendering artifact). Root cause not fully pinned down, but the fix
is verified: use `"p"` (marker-only) legend entries instead of `"lep"`, which
eliminates it entirely; applied to all three new methods. All 6 PNGs
re-rendered clean after each fix and visually verified.
**Percentage decrease: pp24 -23.2458%, PbPb 23+24+25 combined -34.5314%.**

**Step 4 (plotting), review DONE.** `/review-plot` PASS at iteration 1 (log:
`.claude/logs/review-plot-20260906-230220-mupt45-diag-plots.md`). Independent
bin-by-bin check confirmed 0/720 (pp24) and 0/720 (PbPb) bins where the
pT>4.5 series exceeds the pT>4 series -- the physically-required monotonicity
holds everywhere, not just in the totals. All 4 histogram integrals + both
percentages independently reproduced via a fresh ACLiC rebuild + rerun.

## Remaining Work

- **TODO (out of this diagnostic's scope, carried forward for whoever picks
  it up):** PbPb crossx trigger-weighted hist filling
  (`RDFBasedHistFillingPbPb::FillHistogramsCrossx`) is currently BROKEN --
  `EvaluateSingleMuonEffcyPtFitted` throws for some muons because PbPb's
  signal region has not migrated to the pp24 fiducial-gap cut. This
  pre-existing bug is owned by `pp_trig_eff_highpt_jump.md` (ACTIVE, blocked
  on a user decision), not this doc; flagged here per the user's explicit
  instruction when approving the bypass.
- **Human decision pending** on whether to adopt the 4.5 GeV cut. If adopted:
  the permanent change goes in `NTupleProcessingCode/DimuonDataAlgCoreT.c:596`
  (data) + `PythiaFullSimExtras.c`/`PowhegFullSimExtras.c` truth_pt (MC), ALL
  data and MC results must be rerun, and every `_diag_mupt45`/temporary block
  added here (both RDF files' comments, `fill_pbpb_muon_pt45_diag_counts.cxx`,
  the two driver macros, the three new Draw methods in
  `SingleBCrossxPlotterBase.cxx`) must be DELETED.

**Step 6 (2026-09-07), styling + rebinning follow-ups on the pair-eta and pair-pT plots.**
Counts/selection UNCHANGED throughout (verified after every rerun: pp24
705404->541427, Pb+Pb 275622->180446, byte-identical to Step 1-4) -- these are
draw-code-only changes.
1. **Pair-eta-dependence plots -> LINEAR y** (`pair_eta_dependence_mupt_45_vs_40.png`,
   both datasets): pair_eta is not momentum-like and is not log-binned, so per the
   newly adopted linear-y-default rule (`.claude/CLAUDE.md`-adjacent
   `.claude/commands/review-plot.md` R4, added this session; `feedback_log_scale_plots`
   memory corrected accordingly) this plot now defaults to linear y instead of log.
   `SingleBCrossxPlotterBase::DrawPairEtaIntegratedTwoSeries`: removed `SetLogy()`
   and the log-only `ApplyCommonLogYRange` headroom, replaced with a plain linear
   `SetMinimum(0)` / `SetMaximum(max*1.4)`. The pair-pT plots (pair pT is
   power-law/log-binned) correctly KEPT log y -- untouched.
   - An independent `/review-plot` reviewer subagent (spawned before the user said
     further review loops weren't needed for this thread) caught a real pre-existing
     defect surfaced by the rescale: the fixed-NDC `TLegend(0.15,0.78,0.97,0.94)`
     extended above the frame's top border (top margin default 0.9), so its top row
     was visually bisected by the axis line in BOTH datasets, independent of
     log/linear. Fixed: headroom raised to `max*1.4` and the legend moved to
     `(0.15,0.73,0.97,0.90)`, comfortably inside the frame and above the tallest
     point (verified visually post-fix, both datasets).
2. **Scientific-notation y-axis labels** (all datasets, the linear pair-eta plots):
   ROOT's `TAxis` has no per-tick label-format hook (unlike `TF1`/`TGraph`); the only
   available "N x 10^p" typesetting is `TGaxis::SetMaxDigits` (global/static), which
   by default (5 digits) fired for pp24 (6-digit range, "x10^3" header) but not for
   Pb+Pb (5-digit range, bare integers) -- an inconsistent look between the two
   datasets. Forced `TGaxis::SetMaxDigits(3)` (saved/restored around the Draw call)
   so both now show the header consistently (pp24 "x10^6", Pb+Pb "x10^3").
3. **User-requested rebinning VARIANT, log bins from 9 GeV** (2026-09-07): additive,
   NEW PNGs for the two pair-pT-dependence plots per dataset --
   `pair_pt_dependence_mupt_45_vs_40_logbins_from9.png` and
   `pair_pt_in_eta_subplots_mupt_45_vs_40_logbins_from9.png` -- log-spaced bins
   starting at 9 GeV instead of 8, SAME MAX (120 GeV) and SAME bin count (15) as
   `ParamsSet::pT_bins_120`. Opt-in and suffixed per `.claude/CLAUDE.md` Binnings
   item 4 -- does NOT touch `pT_bins_120` and does NOT overwrite the existing
   8-GeV-start PNGs. New additive histograms (`_logbins_from9` /
   `_diag_mupt45_logbins_from9` suffixes) added alongside the existing
   `_diag_mupt45` ones in `RDFBasedHistFillingPP.cxx` (pp24, full crossx pipeline
   rerun via `run_crossx_hist_filling_pp24.sh` to materialize them) and
   `fill_pbpb_muon_pt45_diag_counts.cxx` (Pb+Pb standalone macro, rerun). pair-eta-
   dependence is unaffected by the pair-pT axis and was not remade for this variant.
   User explicitly said no `/review-plot` loop was needed for this rebinning
   follow-up or the scientific-notation change; verified directly (binning,
   counts, rendering) by the orchestrator instead.

## Latest Stage

**DONE 2026-09-07.** All deliverables complete:
- Original 6 PNGs (3 per dataset x 2 datasets) + 2 summary txt files, both
  `/review-analysis-code` and `/review-plot` PASSED at iteration 1 (2026-09-06).
- Step 6 follow-ups (2026-09-07): pair-eta plots remade linear-y with a fixed
  legend/frame overlap and consistent scientific-notation y-labels; 4 NEW
  additive PNGs (2 per dataset) for a user-requested log-bins-from-9-GeV variant
  of the pair-pT plots, verified directly (no formal review loop, per explicit
  user instruction).
- pp24: 705404 (muon pT>4, current) -> 541427 (muon pT>4.5, candidate) raw
  signal-region OS pairs, **-23.2458%**.
- Pb+Pb 2023+2024+2025 combined: 275622 -> 180446 pairs, **-34.5314%**.
- Nothing left except the human decision on whether to adopt the 4.5 GeV cut
  (see Remaining Work) -- not an action item for this doc.
