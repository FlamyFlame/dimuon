# _sub_stats_tables_1 — MC pair-statistics CSV tables (subagent scratch, append-only)

## Objective
Write ONE new macro `plotting_codes/trig_effcy/mc_based/write_mc_pair_statistics_tables.cxx`
that consumes the already-filled Step-3 statistics TH2Ds and writes 6 CSVs per working point
(counts + crossx, for: pT-integrated-over-eta 2-column table, SS pT-vs-eta, OS pT-vs-eta).
Run for pp_full x {Tight, Medium}. NEVER run git. Own only that one source file.

## Context gathered (step 0)
- Input keys verified present in
  `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/mc_trig_eff_hists_pp24_full_step3.root`:
  h_mc_paircount_vs_pt_eta_{ss,os}, h_mc_pairsumw_vs_pt_eta_{ss,os},
  h_mc_pairsumw2_vs_pt_eta_{ss,os} (all TH2D).
- Axes read off the file (NOT retyped): X = 8 bins
  8 / 11.5402 / 16.6472 / 24.0141 / 34.641 / 49.9707 / 72.0843 / 103.984 / 150 GeV;
  Y = 9 bins -2.4 -2 -1.5 -1 -0.5 0.5 1 1.5 2 2.4.
- ANCHORS confirmed by direct ROOT read: SS in-range count 323349, sumw 2.23379 nb;
  OS 2169349, 18.1117 nb. SS all-inclusive (with under/overflow) = 793043 -> the
  out-of-range fraction is large, must be stated in the header.
- Step-3 selection source: `RDFBasedHistFilling/FillMCTrigEffHists.cxx:1097-1104`
  (`step3 selection`): both legs `pass_tight|pass_medium`, pt>4, |eta|<2.4,
  kTruthFidPair (`truth_pt>4 && |truth_eta|<2.4` both legs, l.606),
  kFwdVetoPair (`(pt>7 || q*eta>-2)` both legs, l.637), kGapPair
  (`ParamsSet::FiducialGapCutExpr` both legs, l.662; windows
  ParamsSet.h:436 = (-1.20,-1.05), (-0.06,0.06), (2.30,2.40)); overlay adds
  avg_centrality in [0,5). NO trigger requirement (denominator node, l.1150),
  NO signal-selection cut (doc §3.0(c)).
- Pair tree = `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root`
  (FillMCTrigEffHists.cxx:164) -> resonance vetoes NOT applied.
- Weight units: sumw = sigma_slice*eps_filt*r_isospin/N_slice in **nb**
  (FillMCTrigEffHists.cxx:1159-1162; AMI crossSection is nb, x1000 for pb).
- Output dir convention: WP lives in the TOP-LEVEL dir name (dr_correction_sample_cfg.h
  "VARIANT LAYOUT" block) -> `.../pp_trigger_efficiency/mc_statistics_pt8bins[_medium]/`.
  Derived at runtime from `cfg.out_base`'s parent + the ACTUAL nbinsx of the histogram,
  so the directory name can never lie about the binning.

## Step 1 — macro written (DONE)
`plotting_codes/trig_effcy/mc_based/write_mc_pair_statistics_tables.cxx`
(the ONLY repo file created/edited by this subagent; ACLiC build artifacts
`write_mc_pair_statistics_tables_cxx.{so,d}` + `_ACLiC_dict_rdict.pcm` sit beside it).
Signature: `void write_mc_pair_statistics_tables(const std::string& sample = "pp_full", bool use_tight_wp = true)`.
Design points:
- Input path built exactly like plot_mc_trig_eff.cxx:
  `cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + DrCorrWpSuffix + MCTrigEffPairPt::FileSuffix() + "_step3.root"`.
- Output dir derived at runtime: parent of `cfg.out_base` (asserted to end in a component
  starting with `mc_based`) + `"mc_statistics_pt" + <nbinsX read off the histogram> + "bins"`
  + (`""`|`"_medium"`). WP in the TOP-LEVEL name per the VARIANT LAYOUT block.
- Guards: missing key -> throw with a pointer to the STATISTICS BOOKKEEPING block; axis
  equality across all 6 hists; npt vs MCTrigEffPairPt::NBins(); empty in-range integral -> throw.
- One `Emit(path, writer)` calls the same lambda with std::cout and with the ofstream.
- Gap-cut windows printed via `ParamsSet::FiducialGapCutExpr("q*eta")` (never retyped);
  bin edges printed from the axes; labels formatted from the axes.
- Cross-section error format chosen: **paired columns**. 1D file ->
  `pair_pT_GeV,same_sign_nb,same_sign_stat_err_nb,opposite_sign_nb,opposite_sign_stat_err_nb`.
  2D files -> each pT column is followed by a `<range>_err` column. Documented in the header.
- Extra header facts computed, not assumed: per-sign in-range vs all-inclusive pair counts
  (Tight: SS 323349/793043 = 40.77 %, OS 2169349/3256770 = 66.61 %) explaining why these
  totals are smaller than the dR-histogram totals.

## Step 2 — ran pp_full x {Tight, Medium} (DONE)
Tight  -> /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_statistics_pt8bins/
Medium -> /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_statistics_pt8bins_medium/
6 CSVs each (12 total), all present:
  pair_pt_counts.csv  pair_pt_crossx.csv
  same_sign_pt_vs_eta_counts.csv  same_sign_pt_vs_eta_crossx.csv
  opposite_sign_pt_vs_eta_counts.csv  opposite_sign_pt_vs_eta_crossx.csv

ANCHOR CHECK (Tight) — PASS, exact:
  SS 323349 pairs / 2.23379 +- 0.004427 nb
  OS 2169349 pairs / 18.1117 (18.11167) +- 0.01362 nb
Independent awk row-sum cross-check of the written CSVs reproduces the same three numbers.
Medium (for the record, no anchor supplied): SS 352665 / 2.46256 +- 0.004667 nb;
OS 2385330 / 20.1386 +- 0.01443 nb.

Error path verified: `write_mc_pair_statistics_tables("overlay", true)` throws
"histogram 'h_mc_paircount_vs_pt_eta_ss' is MISSING from .../mc_trig_eff_hists_hijing_overlay_pbpb23_step3.root"
and writes NOTHING (no pbpb .../mc_statistics_* directory created). noovl not run (out of scope).

## Not done / notes
- No git run (per contract). No other repo file touched.
- The header's muon-cut wording is prose citing FillMCTrigEffHists.cxx + §3.0(c); only the gap
  windows and bin edges are machine-generated. If the generic muon cuts change, that prose is
  the one thing in the file that would need a manual update.
