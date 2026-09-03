# pp24 Pythia FullSim: projected statistics for a pT-hat production request

Mode: **IMPLEMENTATION**. Opened 2026-09-03.

## Objective

The user is deciding whether to request additional MC production for the two highest pT-hat
slices of the Pythia8 pp24-conditions FullSim FULL sample (pTHat 70-125 GeV and pTHat 125-300
GeV), from the current 320k events each to 1.2M events each. Two deliverables:

**D1.** Keep the existing pair-pT spectrum plot (markers + error bars, one series per pT-hat
slice, FULL sample) UNCHANGED, and produce a NEW companion plot showing what the spectrum would
look like if the two highest slices had 1.2M events instead of their current statistics.

**D2.** Two NEW CSV tables (same-sign, opposite-sign) of pair statistics on the (8 pair-pT bins) x
(3 sign-independent |pair eta| bins: 0-1, 1-2, 2-2.4) grid, with the two highest pT-hat slices'
contribution to each cell PROJECTED to 1.2M events. Per the user: use the statistics histograms
(absolute values), not numbers read off a PNG.

## Physics Procedure

### 1. Motivation

More MC events in a pT-hat slice do not change the slice's expected contribution to a
cross-section-weighted quantity (the per-pair weight `w = sigma_slice*genFiltEff*r_isospin/N_slice`
scales down as 1/N while the number of pairs scales up as N, so `sum(w)` is an unbiased estimator
of the slice's cross section independent of N) -- but the STATISTICAL UNCERTAINTY on that estimate
does shrink as 1/sqrt(N). "Projected statistics" therefore means: central values (cross section)
UNCHANGED, raw pair counts scaled by N_target/N_current, and sumw2 (the stat-error building block)
scaled by N_current/N_target. This is the SAME mechanism already implemented and used in this repo
for the TEST-to-FULL forecast (`plot_pythia_fullsim_kn_pt_crossx.cxx::ForecastScale`,
`plot_stat_error_forecast`): `SetBinError(err/sqrt(scale_factor))` with central values untouched.
Here the roles are FULL-to-a-LARGER-FULL, i.e. `scale_factor = N_target/N_current` applied to ONLY
the two affected slices, not gated by `g_is_test_sample`.

### 2. Current statistics (measured, not assumed)

Confirmed from the FULL-sample production log
(`/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/trigeff_rerun_tight_20260720_221052.log`,
lines 976 and 1009 -- the `InitInputFullsim` / `ProcessDataHook` event counts, i.e. the
authoritative N_beam the per-pair weight is built from, not an eyeballed number):
```
Fullsim pTH70_125  beam=pp N=319999/319999
Fullsim pTH125_300 beam=pp N=319999/319999
```
Both slices: **N_current = 319999** events (matches the user's "320k"). Target N_target =
1,200,000 for both. Scale factor `sf = N_target / N_current = 1200000/319999 = 3.75001...`,
identical for both slices to 5 digits.

### 3. D1 -- projected pair-pT spectrum

Source: `plotting_codes/pythia_plotting_codes/plot_pythia_fullsim_kn_pt_crossx.cxx`, the per-kn
marker panel of `plot_impl` (kn index 4 = pTHat 70-125, kn index 5 = pTHat 125-300), on the FULL
sample (`is_test_sample=false`, file `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root`,
trees `muon_pair_tree_kin{4,5}_sign2` for the from_same_b truth/reco single-b selection this macro
already applies). For kn4 and kn5 ONLY: after filling the differential histogram exactly as
`plot_impl` does (unchanged bin content), divide each bin's error by `sqrt(sf)` (per §1/§2). kn0-3
untouched (scale factor 1 -- they are not part of the request). New output files, ORIGINAL PNGs
untouched (user instruction: "keep that original plot").

### 4. D2 -- projected SS/OS pair-pT x |pair-eta| statistics CSVs

Source: `plotting_codes/trig_effcy/mc_based/write_mc_pair_statistics_tables.cxx`, which reads the
Step-3 statistics histograms `h_mc_{paircount,pairsumw,pairsumw2}_vs_pt_eta_{ss,os}` booked by
`RDFBasedHistFilling/FillMCTrigEffHists.cxx` on the Step-3 DENOMINATOR node (mc_trigger_efficiency.md
3.0/3.3 selection: generic muon cuts at the nominal WP, pT/eta fiducial, truth fiducial, forward
low-pT veto, q*eta gap cut; NO trigger requirement, NO signal-selection cut -- see that file's own
`header()` for the literal text). These are on `ParamsSet::pair_pt_coarse_bins` (8 log bins, 8-150
GeV) x `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` (9 bins, [-2.4,2.4]).

**Per-slice breakdown (needed to isolate kn4/kn5).** The Step-3 selection is currently only ever
applied to the SIGN-MERGED trees `muon_pair_tree_sign{1,2}` (all six pT-hat slices summed). Per
CLAUDE.md NTuple-Processing Provenance, prefer ntuple-processing OUTPUT over a standalone
re-derivation: the SAME mc_trig pair file
(`muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root`) ALREADY contains the
per-pT-hat-slice trees `muon_pair_tree_kin{0..5}_sign{1,2}` (`PythiaAlgCoreT::fill_kn_trees_fullsim`
-- `FillMuonPairTreePythia` fills BOTH the merged and the per-kn tree for every pair; confirmed by
direct entry-count comparison, byte-identical to the non-mc_trig kin-split file: kin4_sign1=68630,
kin4_sign2=231486, kin5_sign1=77889, kin5_sign2=239754). So a new opt-in `book_kn_stats` flag was
added to `FillMCTrigEffHists.cxx` that applies the EXACT SAME Step-3 selection to
`muon_pair_tree_kin{4,5}_sign{1,2}` and books the SAME three histograms with a `_kn{4,5}` suffix,
additively, in the SAME Step-3 output file. Byte-unchanged for every existing call site
(default false).

**Projection per cell (pair pT bin, pair eta bin, sign):**
```
count_proj = count_all_slices - count_kn4 - count_kn5 + sf*(count_kn4 + count_kn5)
sumw_proj  = sumw_all_slices                                    (UNCHANGED, per §1)
sumw2_proj = sumw2_all_slices - sumw2_kn4 - sumw2_kn5 + (sumw2_kn4 + sumw2_kn5)/sf
```
i.e. remove the two slices' CURRENT contribution and put back their PROJECTED one, leaving kn0-3
untouched -- never touching the merged all-slice histogram directly for those two subtractions
would double-count.

**Pair-eta merge.** The 9-bin coarse axis is merged into the 3 sign-independent |eta^pair| bins
0-1, 1-2, 2-2.4 with `MakeDrEtaGroups` (`dr_correction_cell_groups.h`) -- the SAME utility just
added 2026-09-03 for the dR-correction fit cells (`nocorr_etamerge` mode,
`docs/tracking/mc_trigeff_dr_binning_approaches.md`), reused here rather than re-implemented, per
CLAUDE.md 'never re-invent a binning'. The interior boundaries {1.0, 2.0} are existing edges of
`pair_eta_proj_ranges_coarse_incl_gap` on both sides of 0, so the fold is well defined.

### 5. Negative constraints

1. Do **NOT** touch the nominal (un-suffixed) `plot_pythia_fullsim_kn_pt_crossx.cxx` outputs.
2. Do **NOT** apply the sf factor to `sumw` (§1) -- only to the error-carrying quantities
   (bin error for D1, sumw2 for D2). A larger central value would misrepresent the projection as
   "more cross section", not "more precision".
3. Do **NOT** modify the STATISTICS BOOKKEEPING block already in `FillMCTrigEffHists.cxx`
   (all-slice histograms) -- the new per-kn block is additive only (`book_kn_stats`, default
   false).
4. Do **NOT** invent pair-eta group boundaries -- reuse `MakeDrEtaGroups`.

## Implementation Plan

1. **D2-prep** `RDFBasedHistFilling/FillMCTrigEffHists.cxx`: add opt-in `book_kn_stats` parameter;
   additive per-kn statistics histograms for kn4/kn5 (SS+OS), Step-3 selection, byte-identical to
   the existing STATISTICS BOOKKEEPING block. Per §4.
2. **D2-prep rerun** `FillMCTrigEffHists("pp_full", /*do_step3=*/true, /*tight=*/true, ...,
   book_kn_stats=true)` -- rewrites `mc_trig_eff_hists_pp24_full_step3.root` with the existing
   histograms unchanged plus the 12 new `_kn4`/`_kn5` keys.
3. **D2** New macro: projected SS/OS pair-pT x |eta| statistics CSVs, per §4's cell formula, using
   `MakeDrEtaGroups` for the eta merge.
4. **D1** New function in `plot_pythia_fullsim_kn_pt_crossx.cxx` (or a sibling macro): projected
   marker-panel plot, kn4/kn5 error bars scaled by `1/sqrt(sf)`, new output file names, per §3.
5. Reviews: `/review-analysis-code` (steps 1-2), `/review-plot` (steps 3-4). Docs + commit.

## Design Decisions

**D1. Central values are NEVER rescaled, only statistical-error quantities (bin error / sumw2).**
Per §1: N-scaling is an unbiased-estimator statement about `sum(weight)`, not about the counts
themselves being "more cross section". Confirmed against the existing repo convention
(`ForecastScale`/`plot_stat_error_forecast`), which does exactly this for the TEST->FULL forecast.

**D2. The per-slice breakdown reads the EXISTING per-kn trees already present in the mc_trig pair
file, rather than re-opening the raw NTUP farm or re-deriving the Step-3 selection standalone.**
Per CLAUDE.md NTuple-Processing Provenance rule #1 (prefer ntuple-processing OUTPUT); confirmed the
kin-indexed trees exist in `..._mc_trig_full.root` with entry counts byte-identical to the
non-mc_trig kin-split file used by the existing kn-overlay plot.

**D3. The pair-eta merge reuses `MakeDrEtaGroups` verbatim rather than a new grouping function.**
It is the SAME sign-independent 0-1/1-2/2-2.4 |eta^pair| fold the user adopted today for the
dR-correction fit cells (`nocorr_etamerge`); CLAUDE.md forbids re-inventing a binning that already
exists.

## Progress Log

- 2026-09-03 Step 0: doc created after triage of INDEX.md (`pp24_stats_and_powheg_fullsim_compr.md`
  read for the sibling W1 statistics-CSV precedent; `mc_trigeff_dr_binning_approaches.md` /
  `mc_trigger_efficiency.md` cross-referenced for the dR-correction pair-eta-merge precedent
  reused here) and code reconnaissance: `write_mc_pair_statistics_tables.cxx`,
  `plot_pythia_fullsim_kn_pt_crossx.cxx`, `FillMCTrigEffHists.cxx` (Step-3 selection + STATISTICS
  BOOKKEEPING block), `PythiaAlgCoreT.c` (confirmed `fill_kn_trees_fullsim` fills BOTH the merged
  and per-kn trees in one output file), `dr_correction_cell_groups.h` /
  `dr_correction_sample_cfg.h` (confirmed today's `nocorr_etamerge` |eta| fold and its exact
  0-1/1-2/2-2.4 boundaries match the user's request). Measured current pT-hat-slice statistics
  from the production log: N=319999 events each for pTHat 70-125 and 125-300 (matches "320k").

- 2026-09-03 Step 1 DONE. `FillMCTrigEffHists.cxx`: added `book_kn_stats` parameter (default
  false) + validation (requires `do_step3`); additive per-kn (kn4, kn5) x per-sign STATISTICS
  BOOKKEEPING block, identical Step-3 selection/aliases to the existing all-slice block, writing
  `h_mc_{paircount,pairsumw,pairsumw2}_vs_pt_eta_{ss,os}_kn{4,5}` (12 new keys) into the same
  Step-3 output file. Existing block untouched.

- 2026-09-03 Step 2 DONE. Recompiled and reran
  `FillMCTrigEffHists("pp_full", true, true, false, false, false, false, true)`. Wrote 30 TH1D +
  18 TH2D + 30 TH3D to `mc_trig_eff_hists_pp24_full_step3.root` (12 more TH2D than nominal: the 6
  new `_kn4`/`_kn5` histogram pairs x {count,sumw,sumw2}). Verified directly: e.g.
  `h_mc_paircount_vs_pt_eta_os_kn4` integral 125610 (of 231486 raw kn4 OS pairs before the Step-3
  selection -- consistent), `_kn5` 136380; `_ss_kn4` 30492, `_ss_kn5` 39641. Existing all-slice
  histograms unchanged (`h_mc_paircount_vs_pt_eta_os` integral 2.14321e6, matches the pre-existing
  reference in `pp24_stats_and_powheg_fullsim_compr.md`-adjacent runs).

- 2026-09-03 Step 3 DONE (D2). New macro
  `plotting_codes/trig_effcy/mc_based/write_mc_pair_statistics_tables_projected.cxx`, new shared
  header `Utilities/PtHatKn45ProjectedStats.h` (kNCurrent=319999, kNTarget=1.2e6, kSf, cited to the
  production log). Reuses `MakeDrEtaGroups` for the |eta| merge and `GetDrCorrSample` for paths.
  Output: `same_sign_pt_vs_abs_eta_counts_projected.csv` (current total 318649 -> projected
  511515.6, x1.605) and `opposite_sign_pt_vs_abs_eta_counts_projected.csv` (current 2143213 ->
  projected 2863688.6, x1.336) in `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/`.
  Both 8 pair-pT columns x 3 `|eta^pair|` rows (0-1, 1-2, 2-2.4), per-cell provenance formula in
  the CSV header.

- 2026-09-03 Step 4 DONE (D1). New `plot_impl_projected` + entry point
  `plot_pythia_fullsim_kn_pt_crossx_projected(use_tight_wp)` appended to
  `plot_pythia_fullsim_kn_pt_crossx.cxx`, using the same shared `PtHatKn45Projected` constants.
  Output (new `plots/projected_stats/` subdirectory, FULL-sample plot root): 4 PNGs
  (`{truth,reco}_pair_pt_kn_projected{,_150GeV}.png`). Verified the ORIGINAL
  `plots/reco_pair_pt_kn.png` etc. (mtime 2026-07-20) untouched.

- 2026-09-03 Step 5 DONE. `/review-plot` on D1 (`plot_impl_projected` + the 4 PNGs): PASS at
  iteration 1 (2 INFO-only notes, no fix required) --
  `.claude/logs/review-plot-20260903-121729-kn45-projected-stats-plot.md`. Reviewer independently
  reproduced kn4/kn5 bin content/error via a standalone RDataFrame check: central values
  untouched, error = nominal/sqrt(kSf) with kSf=3.75001 exactly as coded; confirmed original
  nominal PNGs untouched (mtime unchanged).
  `/review-analysis-code` on D2-prep (`FillMCTrigEffHists.cxx` book_kn_stats) +
  `PtHatKn45ProjectedStats.h` + `write_mc_pair_statistics_tables_projected.cxx`: FAIL at
  iteration 1 (1 WARNING -- `CheckSameAxes` never validated the opposite-sign all-slice
  histogram against the reference axis; 2 INFO, not fixed: no table-of-zeros guard, CSV
  delivers counts only not the full count/sumw/sumw2 triad -- both deliberate, documented). Fixed
  the WARNING (added `CheckSameAxes(ref, s.all)` to the loop), recompiled, reran -- byte-identical
  CSV output. PASS at iteration 3 -- `.claude/logs/review-analysis-code-20260903-121611-pthat-stats-projection.md`.
  Reviewer(s) independently re-derived every reported histogram integral and projected total
  (SS 318649->511515.6 x1.605, OS 2143213->2863688.6 x1.336) and confirmed the
  `MakeDrEtaGroups` partition is exact (no double-count/leak) via its internal cover-count
  assertion plus an independent CSV re-sum.

## Remaining Work

None. Deliverables complete and reviewed.

## Latest Stage

(cleared -- see Completion summary below)

## Completion

DONE 2026-09-03. Both deliverables built, reviewed (PASS) and on disk:
* D1 -- `plots/projected_stats/{truth,reco}_pair_pt_kn_projected{,_150GeV}.png` under
  `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/plots/`; original nominal plots
  confirmed untouched.
* D2 -- `same_sign_pt_vs_abs_eta_counts_projected.csv`,
  `opposite_sign_pt_vs_abs_eta_counts_projected.csv` under
  `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/`.
Headline numbers: current 319999 events/slice (confirmed against the production log, matches the
user's "320k"); target 1.2M/slice, sf=3.75001. Same-sign projected pair count 318649 -> 511515.6
(x1.605); opposite-sign 2143213 -> 2863688.6 (x1.336). Central-value cross sections are UNCHANGED
by the projection throughout (Design Decision D1) -- this is a statistics/decision-aid exercise,
not a new measurement.
