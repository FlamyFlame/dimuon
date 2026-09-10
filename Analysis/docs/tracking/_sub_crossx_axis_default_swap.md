# SUB-AGENT scratch: crossx pair-pT axis default swap (150 becomes nominal)

Owner: subagent. Canonical doc = `mu_pt45_gap125_pairpt9_adoption.md` (orchestrator merges).

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run to DONE; do NOT pause to confirm progress.
- Done = DEFAULT crossx pair-pT view is the 9->150 axis; the 9->120 view survives as opt-in
  `_pt_120`; old `_pt_150` output dirs cleaned; every consumer/pipeline/plot label consistent;
  plus the coverage-gap report (which `_pt_150` twins are NOT booked in the hist-filling).
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → record prominently
  here and report rather than guessing).

## Hard constraints
- Edit ONLY `Analysis/plotting_codes/**` and `Analysis/pipelines/**`. Read everything else.
- No git, no compile, no pipeline/macro runs that write the real output tree.
- Do NOT rename any histogram. Unsuffixed names = pT_bins_120 (now ALTERNATIVE);
  `_pt_150` names = pT_bins_150 (now DEFAULT).

## Step 0 — verify nesting arithmetic (DONE)
`ParamsSet.h:457` `fillLogBinningArray` is log10-uniform:
edge_i = 10^(log10(low) + i*(log10(high)-log10(low))/nBins).
- `pair_pt_coarse_bins` = fillLogBinningArray(8, 9.0, 150.0)   [line 773], N_COARSE_PAIR_PT_BINS=8 (line 106)
- `pT_bins_150`        = fillLogBinningArray(16, 9.0, 150.0)   [line 753]
- `pT_bins_120`        = fillLogBinningArray(16, 9.0, 120.0)   [line 752]
2 * 8 = 16 OK. coarse_k = 10^(L + k*D/8); fine_{2k} = 10^(L + 2k*D/16) = 10^(L + k*D/8).
D/16 = (D/8)/2 is exact in binary FP and 2k*(D/16) = k*(D/8) exactly → coarse edges are
bit-identical to every second fine edge of pT_bins_150. pT_bins_120 does NOT nest (different
high edge) → must never bin a correction.

## Step 1 — plan
Invert the `use_pt_bins_150` switch into `use_pt_bins_120` (default false), so the default
path draws the `_pt_150` histogram family into the NOMINAL output dir, and the opt-in draws
the unsuffixed (120) family into a `_pt_120`-suffixed dir. Files to touch (owned dirs only) —
enumerated below as I go.

## Step 1 — DESIGN (decided, before edits)
Two names, both truthful:
- class member `use_pt_bins_120` (default **false** = the 9->150 DEFAULT axis).
- macro argument `also_pt_120` (default false) = "ALSO refresh the opt-in `_pt_120` variant in
  this same invocation". Preserves the old one-call-refreshes-both property (the pipelines'
  comments say the two dirs drifted for months when only one was refreshed) WITHOUT keeping a
  flag whose name means the opposite of its effect.
- New helper in `SingleBCrossxPlotterBase`: `PtAxisHist(name_120)` maps the canonical UNSUFFIXED
  (120) histogram name onto the selected family by replacing the substring `pair_pt` -> `pt_150`.
  Verified against every booked twin: h2d_crossx_**pair_pt**_pair_eta_binned_... ->
  h2d_crossx_**pt_150**_pair_eta_binned_...; h3d_crossx_dr_vs_pair_eta_vs_**pair_pt**_w_signal_cuts
  -> ..._vs_**pt_150**_w_signal_cuts; h2d_counts_**pair_pt**_pair_eta_... ->
  h2d_counts_**pt_150**_pair_eta_...; PbPb h2d_op_crossx_..._vs_**pair_pt**_ctr ->
  ..._vs_**pt_150**_ctr; h2d_op_crossx_dsigma_..._vs_**pair_pt**_ -> ..._vs_**pt_150**_.
  So ONE rule reproduces every existing twin name; call sites type each name once.
- `PtAxisDirSuffix()` = "" (default) or "_pt_120" (opt-in).

Consequence, deliberate: the DEFAULT path now asks for `_pt_150` twins that are NOT booked
(minv/dR 2D). Those calls will THROW with the exact missing name rather than silently drawing a
9->120 figure next to a 9->150 one in the same directory (that mixing is precisely what the
BLOCKING binnings rule forbids). Reported to the orchestrator as the coverage gap.

## Step 2 — EDITS MADE (all inside plotting_codes/ and pipelines/ only)

`plotting_codes/single_b_analysis/SingleBCrossxPlotterBase.cxx`
- NEW public `bool use_pt_bins_120 = false;` (hoisted from the two derived plotters, which each
  carried their own `use_pt_bins_150`).
- NEW `std::string PtAxisHist(const std::string& name_120) const` — returns `name_120` unchanged
  when `use_pt_bins_120`, else swaps the first `pair_pt` for `pt_150`. Throws (naming the string)
  if there is no `pair_pt` to swap.
- NEW `std::string PtAxisDirSuffix() const` — "" or "_pt_120".

`plot_single_b_crossx_pp.cxx`
- Removed the member `use_pt_bins_150` and the whole duplicated `if (use_pt_bins_150){...}` block.
- `Run()` now does `output_dir += PtAxisDirSuffix();` then draws all FIVE figures through
  `PtAxisHist(...)`, so the directory holds one axis only.
- Driver signature `(int run_year=24, const std::string& input_file="", bool also_pt_120=false)`;
  runs the nominal (150) plotter, then optionally a second instance with `use_pt_bins_120=true`.

`plot_single_b_crossx_pbpb.cxx`
- `OutDirName()` appends `PtAxisDirSuffix()`; `DetectAvailableCtrBins()` and `HasCountsHists()`
  now probe the SELECTED family; the per-centrality loop (counts + TAA_weighted) goes through
  `PtAxisHist(...)`; the duplicated `_pt_150` block DELETED.
- Driver `plot_single_b_crossx_pbpb(bool also_pt_120=false)`, same two-instance pattern.

`plot_pp_counts_pair_pt_in_eta.cxx`
- `Run()` appends `PtAxisDirSuffix()` and reads `PtAxisHist("h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts")`.
- CSV provenance line: `pair pT > 8 GeV` -> `ParamsSet::signal_pair_pt_min`; the hard-coded
  `(ParamsSet::pT_bins_120)` -> the axis actually drawn, plus the bin count and the low/high edge
  READ OFF the histogram axis. Header comment updated.
- Driver gains `bool also_pt_120 = false` (4th arg, after `wp`).

`SignalAcceptancePlotter.cxx` + `plot_signal_acceptance_{pythia,powheg}.cxx`
- `use_pt_bins_150` -> `use_pt_bins_120` with the sense inverted: DEFAULT now reads
  `h2d_sig_accept_pt_150_eta` / `_num_pt_150_eta` / `_denom_pt_150_eta` into `pythia/`, `powheg/`;
  opt-in reads the unsuffixed trio into `*_pt_120/`. Drivers take `also_pt_120`.

`plot_dr_vs_pair_pt_diagnostic.cxx`
- Deleted the now-dead private `MakeLogBins()` (it re-derived 15 log bins 8->120); header comment
  corrected to ParamsSet::pT_bins_150, 16 log bins 9 -> 150 GeV (the body already read ParamsSet).

`run_all_crossx.sh`, `pipelines/pipeline_pp_crossx.sh`, `pipelines/pipeline_pbpb_crossx.sh`
- Stage-6 comments rewritten (nominal = 9->150, opt-in = `_pt_120`); the `true` argument now
  means `also_pt_120`. pp Stage 6 also passes `also_pt_120=true` to the counts macro.
- pp Stage 5 liveness probe switched to `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts` (the
  histogram the NOMINAL figure is built from).

## Step 3 — COVERAGE GAP REPORT (orchestrator action required; these files are NOT mine)
The `_pt_150` twins below do NOT exist. With the default inverted, the first two PP entries and
the first four PbPb entries are hit by the NOMINAL plot stage and will THROW
("Missing histogram: ...") until booked.

### `RDFBasedHistFillingPP.cxx` — 120-axis booking line -> required twin name
| 120 booking | line | required `_pt_150` twin | consumer | blocks? |
|---|---|---|---|---|
| `h2d_crossx_pair_pt_minv_w_signal_cuts` | 887 | `h2d_crossx_pt_150_minv_w_signal_cuts` | plot_single_b_crossx_pp (pp24_crossx_pair_pt_minv.png) | **YES** |
| `h2d_crossx_pair_pt_dr_w_signal_cuts` | 890 | `h2d_crossx_pt_150_dr_w_signal_cuts` | plot_single_b_crossx_pp (pp24_crossx_pair_pt_dr.png) | **YES** |
| `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr` | 823 | `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts_no_trig_corr` | plot_crossx_trig_corr_sanity.C | no (macro left on 120) |
| `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts` + every `CrossxCorrectionStages()` suffix (`_corr_raw`, `_corr_unfolded`, `_corr_unfolded_reco`, `_corr_unfolded_reco_trig`) | 855 | same names with `pair_pt` -> `pt_150` | plot_crossx_reco_eff_stages.C | no (macro left on 120) |
| `h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts` | 877 | `h2d_ss_crossx_pt_150_pair_eta_binned_w_signal_cuts` | OS-SS subtraction / R_AA | no consumer yet |
| `h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts` | 894 | `h3d_crossx_minv_vs_pair_eta_vs_pt_150_w_signal_cuts` | none found | no |

### `RDFBasedHistFillingPbPb.cxx` (all per-centrality, `_<ctr>` suffix)
| 120 booking | line | required twin | consumer | blocks? |
|---|---|---|---|---|
| `h2d_crossx_pair_pt_minv_w_signal_cuts_<ctr>` | 1230 | `h2d_crossx_pt_150_minv_w_signal_cuts_<ctr>` | pbpb TAA_weighted `*_pair_pt_minv.png` | **YES** |
| `h2d_crossx_pair_pt_dr_w_signal_cuts_<ctr>` | 1235 | `h2d_crossx_pt_150_dr_w_signal_cuts_<ctr>` | pbpb TAA_weighted `*_pair_pt_dr.png` | **YES** |
| `h2d_crossx_pair_pt_minv_w_signal_cuts_<ctr>_counts` | 1293 | `h2d_crossx_pt_150_minv_w_signal_cuts_<ctr>_counts` | pbpb counts `*_pair_pt_minv.png` | **YES** |
| `h2d_crossx_pair_pt_dr_w_signal_cuts_<ctr>_counts` | 1298 | `h2d_crossx_pt_150_dr_w_signal_cuts_<ctr>_counts` | pbpb counts `*_pair_pt_dr.png` | **YES** |
| `h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>_no_trig_corr` | 1273 | `..._vs_pt_150_<ctr>_no_trig_corr` | plot_crossx_trig_corr_sanity.C | no |
| `h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>` + `CrossxCorrectionStages()` suffixes | 1281 | `..._vs_pt_150_<ctr><suffix>` | plot_crossx_reco_eff_stages.C | no |
| `h2d_crossx_pair_pt_minv_dsigma_<ctr>` / `h2d_crossx_pair_pt_dr_dsigma_<ctr>` | 1247 / 1251 | `h2d_crossx_pt_150_minv_dsigma_<ctr>` / `h2d_crossx_pt_150_dr_dsigma_<ctr>` | none found (the eta/dR dsigma 150 twins DO exist) | no |
| `h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts_<ctr>` (+ `_counts`) | 1256 / 1303 | `..._vs_pt_150_...` | none found | no |
| `h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt` | 1170 | `h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pt_150` | R_AA input (global 3D) | no consumer yet, but R_AA would sit on the ALTERNATIVE axis |
| `h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt` | 1204 | `..._vs_pt_150` | R_AA OS-SS subtraction | same |

Twins that DO exist (no action): PP 907 / 915 / 918; PbPb 1320 / 1324 / 1330 / 1334 / 1339 / 1343;
signal acceptance `h2d_sig_accept_{,num_,denom_}pt_150_eta` (PythiaTruth 495-500 + 578-583,
PowhegTruth 313-318 + 340-...); the 1D `h1d_crossx_pair_pt_150_w_signal_cuts_{op,ss}_dsigma`
(PP 798, already on pT_bins_150).

### Sanity macros left on the 120 axis DELIBERATELY
`plot_crossx_reco_eff_stages.C` and `plot_crossx_trig_corr_sanity.C` have NO `_pt_150` inputs at
all, so they were NOT switched — a switch pointing at unbooked histograms would just crash, and
faking it by keeping them on 120 while the crossx moved to 150 is a real inconsistency, not a
cosmetic one. Once the twins above are booked, both macros should get the same
`use_pt_bins_120` / `PtAxisHist` treatment. Flagged, not guessed at.

## Step 4 — nesting arithmetic (recomputed, not assumed) — see Step 0
2 * N_COARSE_PAIR_PT_BINS (8) = 16 = pT_bins_150 bin count. Both are log10-uniform over the
IDENTICAL range 9 -> 150, so coarse edge k == fine edge 2k bit-for-bit
(D/16 = (D/8)/2 exactly in binary FP; 2k*(D/16) == k*(D/8) exactly).

## Step 5 — old `_pt_150` output directories
`rm -rf` was DENIED by the permission system, so nothing was deleted. All ten dirs contain ONLY
PNGs (0 non-png files each). Under
`/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/`:
SUPERSEDED, safe to delete (5): `pp24_pt_150` (3 PNGs), `pbpb_23_24_25_combined_pt_150` (24),
`pbpb_23_24_combined_pt_150` (24, obsolete 2-year combination), `pythia_pt_150` (2),
`powheg_pt_150` (2).
DELIBERATE dated snapshots — FLAGGED, not deleted (5): `pp24_pt_150_backup_20260615`,
`pp24_pt_150_backup_20260616_pre_reco_nominal`,
`pbpb_23_24_25_combined_pt_150_backup_20260615`,
`pbpb_23_24_25_combined_pt_150_backup_20260616_pre_reco_nominal`,
`pbpb_23_24_combined_pt_150_backup_20260505`.

## Latest Stage — DONE (subagent finished 2026-09-08)
All edits are in place inside the two owned directories; nothing was compiled, run or committed
(per the constraints). `grep -rn "use_pt_bins_150" plotting_codes/ pipelines/` returns nothing.
Two items are handed BACK to the orchestrator: (1) book the missing `_pt_150` twins listed in
Step 3 (the four PbPb + two PP "BLOCKS" rows are hard blockers for the nominal plot stage), and
(2) delete the five superseded `_pt_150` plot directories in Step 5 (the `rm` was permission-denied
here) and decide what to do with the five dated `_backup_` snapshots.

# =============================================================================================
# ROUND 2 (2026-09-09) — MAPPING INVERTED at the coordinator's instruction
# =============================================================================================

## Why the design changed
My Round-1 coverage-gap report showed the partial `_pt_150` family had no member for the two
Pb+Pb R_AA global 3Ds (`h3d_{op,ss}_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt`,
RDFBasedHistFillingPbPb.cxx:1170/1204), so promoting it to nominal would have left the R_AA input
and the cross-section on DIFFERENT pair-pT binnings. Instead the producers reboot the COMPLETE
UNSUFFIXED family onto `ParamsSet::pT_bins_150` (9->150) and demote the partial family onto
`pT_bins_120` (9->120), renaming its token `pt_150` -> `pt_120`. The coverage gap disappears.
NOTE: the ParamsSet VECTORS are unchanged — `pT_bins_150` is still 16 log bins 9->150 and is the
NOMINAL axis. Only which histogram family is booked on which vector changed.

## Inverted mapping (Round-2 edits)
`SingleBCrossxPlotterBase.cxx`
- `PtAxisHist(nominal_name)`: DEFAULT (`use_pt_bins_120 == false`) is now the IDENTITY — returns
  the unsuffixed canonical name. Opt-in maps `pair_pt` -> `pt_120` (was `pt_150`). Comment block
  rewritten with the reason above.
- NEW `bool SkipInAltView(hname, png_name)` called at the top of `Save2DColz`,
  `DrawPairPtByEtaWithDrLines` and `DrawPairPtByEta`: in the OPT-IN view only, a figure whose
  histogram was never booked is SKIPPED with a log line; in the NOMINAL view it still THROWS.
  Needed because the `pt_120` family is deliberately partial (no minv, no dR-2D, no correction
  stages, no SS) and `also_pt_120=true` would otherwise abort Stage 6 after drawing the nominal
  figures. The asymmetry is deliberate: a missing NOMINAL histogram is a producer bug.

`plot_single_b_crossx_pp.cxx`, `plot_single_b_crossx_pbpb.cxx` — no code change needed; they
already call `PtAxisHist()` with canonical unsuffixed names, which are now the nominal ones.

`plot_pp_counts_pair_pt_in_eta.cxx` — header input list inverted (unsuffixed = nominal,
`h2d_counts_pt_120_...` = opt-in); `WriteCsv` guarded by `SkipInAltView`. The CSV's binning label
`use_pt_bins_120 ? pT_bins_120 : pT_bins_150` was already correct and is unchanged.

`SignalAcceptancePlotter.cxx` — DEFAULT now reads the UNSUFFIXED trio `h2d_sig_accept_pt_eta`,
`_num_pt_eta`, `_denom_pt_eta`; opt-in reads `..._pt_120_eta`.

`plot_mc_data_pair_pt_in_eta.cxx` — `kDataHist` `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`
-> `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`; header provenance updated (also the stale
"pair pT > 8" -> `ParamsSet::signal_pair_pt_min`).

`pipelines/pipeline_pp_crossx.sh` — Stage-5 liveness probe reverted to the unsuffixed
`h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts` (the `pt_150` name will not exist); Stage-6 and
counts comments rewritten. `pipeline_pbpb_crossx.sh` Stage-6 comment rewritten (now also states
that the R_AA globals belong to the nominal family — the reason for the inversion).

Stale bin-count text corrected (comments only; the vector named was already right):
`plot_mc_trig_eff_closure.cxx:30`, `plot_mc_trig_eff_closure_compare.cxx:25`,
`pipelines/run_mc_trigeff_closure.sh:28`: "15 log bins 8-150" -> "16 log bins 9-150".

## Item 2 — the two sanity macros
`plot_crossx_reco_eff_stages.C` and `plot_crossx_trig_corr_sanity.C` reference ONLY unsuffixed
names (`h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts{,_no_trig_corr,+stage suffixes}` and the
Pb+Pb `h2d_op_crossx_..._vs_pair_eta_vs_pair_pt_ctr*`). Under the new design those ARE the nominal
9->150 objects, so both macros are correct with NO change and need no alternative-axis handling.
Verified by grep: neither file contains the token `pt_150`.

## Item 5 — coordinator's pre-existing edits: ALL INTACT, no conflict
- `plot_sig_accept_cutflow_above_60GeV.cxx:44-48` — `kCuts()` reads `SignalPairPtCutExpr`,
  two `FiducialGapCutExpr`, `PairFiducialEtaCutExpr` from ParamsSet. Untouched by me.
- `pp24_secondary_vertex_stats.cxx:156-159` — `signal_cuts` uses `SignalPairPtCutExpr` +
  fiducial/pair-eta windows. Untouched by me.
- `plot_dr_vs_pair_pt_diagnostic.cxx:133-144` — `truth_cuts` / `data_cuts` on the fiducial +
  pair-level windows and `SignalPairPtCutExpr`. INTACT. My only edits to that file were a header
  comment and deleting the dead private `MakeLogBins()` (which re-derived 15 log bins 8->120 and
  was referenced nowhere after the body moved to `pms.pT_bins_150`). The cut block was not touched.

## OPEN ITEMS FOR THE COORDINATOR
### (a) var1D / producer entries — I did NOT rename these; renaming them to `_120` would be WRONG
These are single-member variables ALREADY bound to `"binning": "pT_bins_150"`, i.e. already on the
NOMINAL 9->150 axis. There is no second, alternative-axis member. Renaming them `_150` -> `_120`
would rebind the ONLY 1D/truth pair-pT spectrum to the 9->120 alternative and silently move the
whole MC-vs-data comparison off the nominal axis — the failure this inversion exists to avoid.
RECOMMENDATION: drop the token (`_150` -> nothing), keep `"binning": "pT_bins_150"`.
  - `RDFBasedHistFilling/var1D_pythia_truth.json:28`  `pair_pt_log_150`      (binning pT_bins_150)
  - `RDFBasedHistFilling/var1D_pythia_fullsim.json:250-253` `truth_pair_pt_log_150`
  - `RDFBasedHistFilling/var1D_powheg_fullsim.json:208-211` `truth_pair_pt_log_150`
  - `RDFBasedHistFilling/var1D_powheg_truth.json:173`  (binning pT_bins_150)
  - `RDFBasedHistFilling/var1D_pp.json:93`             (binning pT_bins_150)
  - `RDFBasedHistFilling/var1D_pbpb.json:123`          (binning pT_bins_150)
  - `RDFBasedHistFillingPP.cxx:798` signal_region_1d_vars entry `{"pair_pt_150","pair_pt","pT_bins_150",...}`
  - `RDFBasedHistFillingPythiaFullsim.cxx:40,63` and `RDFBasedHistFillingPowhegFullsim.cxx:54,62`
    reference `truth_pair_pt_log_150`.
My consumers (`PlotMCDataComprBaseClass.c:38` data_var `pair_pt_150`, mc_var
`truth_pair_pt_log_150`; `plot_mc_data_pair_pt_in_eta.cxx` `kMcHist`/`kPowhegHist`) still use the
CURRENT names — correct on the nominal axis today either way. One-line follow-up on my side if the
token is dropped.

### (b) Signal-acceptance producers are NOT in the coordinator's stated scope
The `h2d_sig_accept_*` families live in `RDFBasedHistFillingPythiaTruth.cxx` (495-500, 578-583) and
`RDFBasedHistFillingPowhegTruth.cxx` (313-318, 340+), not in PP/PbPb. For my inverted plotter to be
right they need the SAME treatment: unsuffixed `h2d_sig_accept_{,num_,denom_}pt_eta` rebooked on
`pT_bins_150`, and `*_pt_150_eta` -> `*_pt_120_eta` on `pT_bins_120`.

### (c) Plot directories — NOT deleted (permission denied; coordinator will handle)
Superseded `*_pt_150` dirs under plots/single_b_analysis/: `pp24_pt_150`,
`pbpb_23_24_25_combined_pt_150`, `pbpb_23_24_combined_pt_150`, `pythia_pt_150`, `powheg_pt_150`.
Dated snapshots to keep or decide separately: `pp24_pt_150_backup_20260615`,
`pp24_pt_150_backup_20260616_pre_reco_nominal`, `pbpb_23_24_25_combined_pt_150_backup_20260615`,
`pbpb_23_24_25_combined_pt_150_backup_20260616_pre_reco_nominal`,
`pbpb_23_24_combined_pt_150_backup_20260505`. All contain PNGs only.

## Latest Stage — ROUND 2 DONE
Nothing compiled, run, git-ed or deleted. No live code in plotting_codes/ or pipelines/ names a
`pt_150` histogram any more except the var1D variable names in (a), which are pending the
coordinator's decision.
