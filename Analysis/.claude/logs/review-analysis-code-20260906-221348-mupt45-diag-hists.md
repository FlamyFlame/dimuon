# Analysis Code Review Log
**Task**: Review temporary diagnostic `_diag_mupt45` counts-histogram additions to RDFBasedHistFillingPP.cxx and RDFBasedHistFillingPbPb.cxx (muon reconstructed pT>4.5 GeV cut statistics diagnostic; see docs/tracking/muon_pt45_cut_diagnostic.md).
**Log file**: review-analysis-code-20260906-221348-mupt45-diag-hists.md
**Started**: 2026-09-06 22:13:48
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Step 1: Execute — summary

Files modified/created:
1. `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` — additive `_diag_mupt45` block right after
   the existing `h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts` booking. Filters
   `df_single_b_crossx_weighted` with `"m1.pt > 4.5 && m2.pt > 4.5"`, books an unweighted
   `h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45` on the same axes as the nominal
   counts histogram. No existing code touched.
2. `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` — **tried and reverted.** An analogous
   additive block inside the per-centrality loop (deriving first from the post-weight
   `df_crossx_ctr`, then from the pre-weight `df_single_b_crossx`) both hit a PRE-EXISTING,
   unrelated bug: `EvaluateSingleMuonEffcyPtFitted: no fitted turn-on ... on the CONTIGUOUS coarse
   q*eta binning` thrown mid-event-loop and swallowed by ROOT (exits 0, near-empty 398-byte output
   file). This is the already-tracked ACTIVE issue in `docs/tracking/pp_trig_eff_highpt_jump.md`
   (PbPb has not migrated to the pp24 fiducial-gap cut). Because RDF's RunGraphs shares ONE event
   loop across every result tied to a source tree, the throw kills ALL histograms from that loop —
   including ones (like mine) that never read the poisoned column — so no additive block inside
   `FillHistogramsCrossx` could work regardless of which pre-existing node it derived from. File
   is back to its pristine, unmodified state (`git diff` clean) after reverting both attempts.
3. `RDFBasedHistFilling/fill_pbpb_muon_pt45_diag_counts.cxx` — **NEW standalone macro** (asked the
   user first via AskUserQuestion; approved: "bypass it, but be careful to ensure the counts
   results are correct, and mark the bug as a future TO-DO in tracking doc"). Opens a fresh
   `TChain("muon_pair_tree_sign2")` (OS) over the 3 years' `muon_pairs_pbpb_20YY_single_mu4_mindR_0_02.root`
   combined ntuple-processing output directly, entirely decoupled from
   `RDFBasedHistFillingPbPb`/its poisoned trigger-efficiency machinery, and mirrors PbPb's EXACT
   current signal region + Tight-WP filter verbatim from `RDFBasedHistFillingPbPb.cxx:977,987`
   (`pair_pass_tight && minv>1.08 && minv<2.9 && pair_pt>8 && m1.charge*m1.eta<2.2 &&
   m2.charge*m2.eta<2.2`), then books two unweighted counts histograms (muon pT>4 baseline, muon
   pT>4.5 diagnostic) on the same axes as pp24's (`ParamsSet::pT_bins_120` x
   `ParamsSet::N_PAIR_ETA_CROSSX_BINS`/`PAIR_ETA_CROSSX_MIN/MAX`). Combining years is a plain
   TChain union — correct for an unweighted count (no luminosity weighting needed), matching the
   "simple sum for `_counts` histograms" rule already used by
   `SingleBCrossxPlotterPbPbCombined::GetHistObject`.

Compilation/runtime status:
- `RDFBasedHistFillingPP.cxx`: ACLiC-compiled clean (only pre-existing unrelated `-Wsign-compare`
  warnings). Ran via `run_crossx_hist_filling_pp24.sh` — event loop completed successfully, no
  errors, new histogram present and populated.
- `RDFBasedHistFillingPbPb.cxx`: reverted to pristine; not part of the final diff.
- `fill_pbpb_muon_pt45_diag_counts.cxx`: ACLiC-compiled clean (same pre-existing warnings only).
  Ran via `root -l -b -q 'fill_pbpb_muon_pt45_diag_counts.cxx+()'` — completed with no errors,
  wrote both histograms.

Numbers computed:
- **pp24** (from `histograms_real_pairs_pp_2024_2mu4_nominal.root`): nominal (pT>4) = 705404,
  diagnostic (pT>4.5) = 541427 → **23.2458 % decrease**.
- **PbPb 23+24+25 combined** (from
  `histograms_pbpb_23_24_25_muon_pt45_diag_counts.root`): baseline (pT>4) = 275622, diagnostic
  (pT>4.5) = 180446 → **34.5314 % decrease**.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: None found. All checks (additive-only diff, sign convention, exact signal-cut mirroring, binning, file paths, no-trigger-column-touched, comment accuracy) passed.
**Numerical verification**: All 4 histogram integrals and both percentages independently reproduced from the actual output ROOT files — exact MATCH. `Integral(0,-1,0,-1)` (incl. under/overflow) confirmed as the correct convention (a handful of pairs per histogram fall outside axis ranges and would be silently dropped by plain `Integral()`; difference <0.01%, does not change either headline percentage). Magnitude/sign of both percentage decreases (pp24 23.2%, PbPb 34.5%) judged physically plausible for a steeply falling dimuon-pT spectrum under a compounding both-muon pT-floor increase.

**Status**: APPROVED at iteration 1
**Summary**: Both the pp24 additive RDF block and the PbPb standalone-macro bypass are correct, additive/non-invasive, and produce independently-verified, physically plausible raw pair-count numbers (pp24 705404→541427, -23.2458%; PbPb 23+24+25 275622→180446, -34.5314%).
