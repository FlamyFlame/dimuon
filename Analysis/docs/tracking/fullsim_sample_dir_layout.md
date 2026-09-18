# Pythia fullsim sample directories: per-sample layout for the MC-only derived products

## Objective

The `~/usatlasdata/pythia_fullsim_*` sample directories (pp24 FULL, pp24 TEST, HIJING overlay
TEST pbpb23, r17663 no-overlay, and the new HIJING overlay TEST pbpb24 + the future overlay FULL
sample) hold ~200 flat ROOT files each: raw NTUP, ntuple-processing output, RDF hist-filling
output, and the several dozen MC-trigger-efficiency / reco-efficiency products. Put the MC-ONLY
derived products into subtrees, in **both** the code that writes them and every code that reads
them (MC trig-eff Steps 1–4, ΔR fits, closure, pair-eff, reco-eff, the pp24 crossx RDF stage),
move the existing files on disk, and retire the redundant `mc_trig_eff_fit_plots/` PNGs.
Also: make the HIJING-overlay sample year an explicit knob of the MC trig-eff code (default
pbpb24 = `pythia_fullsim_hijing_overlay_test_sample/`, 23 = `..._pbpb23/`).

## Autonomy Contract (DONE 2026-09-16 — all 6 Done items met; review PASS at iteration 3)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. One layout helper (single source of truth) in `MuonObjectsParamsAndHelpers/FullSimSampleType.h`;
     `DrCorrSample` and every producer/consumer compose paths from it — zero retyped
     `mc_trig_eff_hists_`/`dr_correction_`/`single_mu_effcy_pT_fit_mc`/`pair_trig_eff`/
     `pair_reco_eff`/`mc_trig_eff_closure` paths outside the helper (grep-verified, C++ AND shell).
  2. `single_mu_effcy_pT_fit_mc` carries the sample label (mc_trigger_efficiency.md RW 3c closed).
  3. `mc_trig_eff_fit_plots/` PNG production removed from `FitMCSinglesEffcy.cxx`; existing dirs
     archived under `backup/`, not deleted.
  4. Overlay year knob (`int overlay_year = 24`, settable 23) in every MC trig-eff macro +
     pipeline; `FullSimSampleType.h` default overlay TEST dir = pbpb24; sanity-check codes that
     explicitly study r17618/r17662/r17663 point at `_pbpb23/` explicitly.
  5. Disk: all 5 sample dirs reorganized (products moved, `.bak_*` → `backup/`, `*.log` → `logs/`,
     CSVs → `plots/projected_stats/`), every moved file md5-identical, nothing deleted.
  6. Verification: every C++ consumer recompiles (ACLiC); a read-only smoke test of each reader
     family resolves and opens its inputs at the new paths (pp_full Tight+Medium, overlay 23,
     noovl); the pp24 crossx RDF stage `OpenPairEfficiencyInputs` loads; `/review-analysis-code`
     on the code change; docs (`Analysis/README.md`, root `README.md`, pipeline docs) updated;
     committed.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

No physics changes. The invariant this task must preserve:

1. **Byte-identity.** Every consumer opens the SAME file it opened before (moved, md5-identical),
   so every number and plot downstream is unchanged. Verified by md5 before/after the move and by
   resolving each consumer's path at the new location.
2. **Sample identity is explicit.** The HIJING-overlay TEST sample now exists in two productions
   (pbpb23 conditions r17618 in `_pbpb23/`; pbpb24 conditions r17864 in the unsuffixed dir).
   Which one a step reads is a declared knob (`overlay_year`), never an implicit default that
   changed underneath the code. The AMI rule (`Analysis/docs/ami_weights.md`) is enforced by
   directory: `PythiaAlgCoreT` reads a sample's OWN `ami_info/` whenever the sample directory
   has one (FULL pp24; the pbpb24 overlay, whose e8613 evgen is registered as table B2), and
   falls back to the truth production's `ami_info/` only for the e8599 test samples that ship
   none (pp24 test, pbpb23 overlay). Fixed 2026-09-16 (review iteration 1): before, EVERY test
   sample read the truth production's file, so the pbpb24 sample's own AMI was never read.
   **SUPERSEDED 2026-09-17 (user ruling):** the "sample's own `ami_info/`" preference was wrong
   in principle — the weight of every fullsim/overlay event is its PYTHIA EVGEN's σ·ε_filt, never
   the AOD chain's (e8613 is the HIJING evgen; the pbpb24 overlay is 802776 e8599 = table A). AMI
   files are now keyed by evgen (`FullSimSampleType.h` `PythiaEvgen`: `ami_info_nPDF/` /
   `ami_info_PDF/` under the truth sample dir); no sample dir holds an `ami_info/`. See
   `ami_weights.md` and `hijing_overlay_pbpb24_test_sample_skim.md` §"AMI-weight correction".
3. **Negative constraints.** No file is deleted. No cut, binning, weight, fit range or sample
   content changes. `.bak_*` copies are archived, never re-read by any code. The NTUP files, the
   LGD symlink farm, `ami_info/` (2026-09-16 state; since renamed `ami_info_AOD_record_*`, see above), `merging-record.txt` and the grid-monitor state stay at the
   sample root (SkimCode/grid_monitor and the LGD farm scripts are NOT touched).

## Context

- Data directories (`dimuon_data/pp_2024/`, `pbpb_2024/`) keep `muon_pairs_*` and
  `histograms_*` flat and put only the turn-on fits in a subdir (`trg_effcy_pT_fitting_to_*`).
  User decision 2026-09-16: MC mirrors that — NTP output + RDF hist-filling output stay flat;
  MC-only products (trigger & reco efficiencies) go into subtrees.
- `mu_pt45_gap125_pairpt9_adoption.md`: the 4.5 GeV rerun is COMPLETE, so this is a pure move.
- `hijing_overlay_pbpb24_test_sample_skim.md`: pbpb24 test sample (single slice pTH125_300,
  r17864, calorimeter broken) is in the unsuffixed dir since 2026-09-15; the pbpb23 sample moved
  to `_pbpb23/`. Until this task, every overlay code path still pointed at the unsuffixed dir.
- `mc_trigger_efficiency.md` RW 3c: append the label to `single_mu_effcy_pT_fit_mc`.
- `FitMCSinglesEffcy.cxx` PNGs are a strict subset of
  `step1_singles_data_mc/step1_eff_pt_in_q_eta_bins_<charge>.png` (same MC points + same fit
  TF1, plus data + ratio there) → retired (user decision 2026-09-16).

## Scope

Sample dirs: `pythia_fullsim_full_sample`, `pythia_fullsim_test_sample`,
`pythia_fullsim_hijing_overlay_test_sample_pbpb23`, `pythia_fullsim_no_overlay_test_sample`,
`pythia_fullsim_hijing_overlay_test_sample` (pbpb24; NTUP only today — layout applies when
products appear). Code: everything under `Analysis/` that composes a path into these dirs
(`old/`, `OldRDF*`, `pythia_backup/`, `reference_codes/`, `docs/tracking/scripts_*` excluded —
dead code, git-backed).

## Design Decisions

- **D1 Layout** (user, 2026-09-16):
  ```
  <sample>/                       raw NTUP, ami_info/, merging-record.txt, r*.txt, grid_monitor_*  (unchanged)
    muon_pairs_*.root, hists_pythia_ntuple_processing_*.root, histograms_pythia_fullsim_*.root  (flat, as data)
    mc_trig_eff/hists/            mc_trig_eff_hists_<label>*.root      (FillMCTrigEffHists: steps 1-4, sanity, corrected, nogapcut)
    mc_trig_eff/singles_fits/     single_mu_effcy_pT_fit_mc_<label>*.root (FitMCSinglesEffcy)
    mc_trig_eff/dr_correction/    dr_correction_plateaus_*.root, dr_correction_fits_*.root (plot_mc_trig_eff / fit_dr_corrections / combine)
    mc_trig_eff/pair_eff/         pair_trig_eff_<label>*.root          (FillMCTrigEffPairEff)
    mc_trig_eff/closure/          mc_trig_eff_closure_<label>*.root    (FillMCTrigEffClosure)
    reco_eff/                     pair_reco_eff_<label>.root           (build_pp24_fullsim_pair_reco_eff.C)
    plots/                        unchanged (+ projected_stats/ CSVs)
    backup/                       *.bak_* ROOT copies, retired mc_trig_eff_fit_plots/
    logs/                         *.log
  ```
- **D2 Single source of truth**: layout functions live in `FullSimSampleType.h` (already the one
  place naming the sample dirs). `DrCorrSample.mc_dir` is RENAMED to `sample_dir` so the compiler
  finds every C++ consumer; each composes `FullSimMCTrigEff*Dir(sample_dir)`. Shell pipelines
  cannot be compiler-checked → grep sweep + two independent enumeration subagents, cross-checked.
- **D3 Overlay year knob** (user, 2026-09-16): `int overlay_year = 24` parameter (WP-config-var
  pattern) on `GetDrCorrSample`, `FillMCTrigEffHists`, `FitMCSinglesEffcy`, the plotters and the
  pair-eff/closure fillers; `OVERLAY_YEAR=${OVERLAY_YEAR:-24}` in pipelines. Selects dir
  (`..._test_sample/` vs `..._test_sample_pbpb23/`), label (`hijing_overlay_pbpb24` /
  `hijing_overlay_pbpb23`), headline text ("2024 conditions"/"2023 conditions") and plot root
  (`pbpb_trigger_efficiency/mc_based<tag>/pbpb<yy>/`, the year as the leaf exactly like the
  data-side `mu4/no_corr/pbpb<yy>`). `FullSimSampleType.h` gets the same `pbpb_year` parameter
  (default 24) on `FullSimSampleInputDir/Label/PlotDir` for the hijing type.
- **D4 Fit PNGs retired** (user): PNG block deleted from `FitMCSinglesEffcy.cxx`; the Step-1
  panel is the one figure of the MC turn-on fit.
- **D5 Fit-file label** (user; RW 3c): `single_mu_effcy_pT_fit_mc_<label>[_corrected[_sfclosure]][_medium_wp].root`,
  composed by ONE helper `FullSimMCSinglesFitFile(...)`.

- **D6 Overlay data reference follows the conditions year** (agent decision 2026-09-16, flagged
  to the user in the report): the Step-1 / corrected-MC DATA tag-and-probe reference of the overlay
  is the Pb+Pb data of the overlay's conditions year, 0-5 % (`DrCorrDataRefDir/Tag/Text`):
  pbpb23 → `dimuon_data/pbpb_2023/` (= the D2 choice of `mc_trigger_efficiency.md`, unchanged),
  pbpb24 → `dimuon_data/pbpb_2024/` (like-for-like; both files exist). Physics rationale: the
  comparison exists to validate the simulation of the trigger under the conditions it simulates.
- **D7 NTUP `.bak_*` and grid-monitor / LGD files stay at the root.** `grid_monitor.sh:423` scans
  `<NTUP>.bak_YYYYMMDD.root` by name pattern; `lgd_farm.log`, `lgd_rules.txt`,
  `grid_monitor_*.{log,txt}` are written by SkimCode scripts that are out of scope. Only
  PRODUCT `.bak_*` copies and pipeline logs move.
- **D8 Overlay plot root** = `pbpb_trigger_efficiency/mc_based<tag>/pbpb<yy>/` (year as the
  leaf, the data-side `mu4/no_corr/pbpb<yy>` convention). No live overlay MC plot tree existed
  (all archived under `old_4GeV_cut/` on 2026-09-15), so nothing on disk had to move.

## Implementation Plan

1. [x] Write plan (this doc), register in INDEX.md.
2. [x] Enumerate every path-composing site: my grep + 2 independent Explore subagents; cross-check
       (both agents agree on the C++/shell set; agent A added the doc/`.claude`/SkimCode README
       mentions and `IntNotes/figures/figure_provenance.json` (plots/ paths — unchanged)).
3. [x] Layout helper + year knob in `FullSimSampleType.h`; `dr_correction_sample_cfg.h` rebuilt on it.
4. [x] Producers: FillMCTrigEffHists, FitMCSinglesEffcy (+PNG removal), plot_mc_trig_eff (plateau file),
       fit_dr_corrections, combine_dr_correction_fits, FillMCTrigEffPairEff, FillMCTrigEffClosure,
       build_pp24_fullsim_pair_reco_eff.C, write_mc_pair_statistics_tables_projected (CSVs →
       plots/projected_stats/).
5. [x] Consumers: RDFBasedHistFillingPP (crossx), dr_correction_apply.h, PairTrigEffEvaluator,
       plot_* (mc_based, reco_effcy), write_* tables; pipelines via the new shell twin
       `pipelines/fullsim_sample_layout.sh` (run_dr_correction_fits, run_mc_trigeff_{round7,
       corrected,closure,pair_eff}, run_pp_fullsim_trigeff_fullsample, run_ntp_singles_round7_nvtx,
       pipeline_pythia_fullsim_pp (backups → backup/), pipeline_pythia_fullsim_overlay).
6. [x] Overlay repoint: `FullSimSampleType.h` default = pbpb24; `overlay_year`/`overlay_pbpb_year`/
       `OVERLAY_YEAR` knobs; explicit `_pbpb23` in plot_{pair,single_muon}_reco_effcy_r17618_vs_r17662,
       run_pythia_fullsim_overlay_r17662_nodr.sh, plot_mc_trig_eff noovl comparison series,
       run_ntp_singles_round7_nvtx (default 23).
7. [x] Disk move (md5 before/after), rename singles-fit files, archive PNG dirs, backup/, logs/.
8. [x] Recompile everything touched (ACLiC); smoke-test each reader family; `/review-analysis-code`.
9. [x] Docs (README, docs/, ami_weights.md), INDEX, commit(s).

## Progress Log

### 2026-09-16 steps 1-6 (code)
- Enumeration: two independent Explore subagents (bottom-up from writers; top-down from directory
  names) + own grep. Cross-checked: identical C++/shell site lists; extras from A = docs,
  `.claude/{agents/plot-reviewer.md,commands/review-plot.md}`, `SkimCode/README.md` (already
  correct), `IntNotes/figures/figure_provenance.json` (plots/ paths, unchanged by D1).
- `FullSimSampleType.h`: `FullSimCheckPbPbYear`, `pbpb_year = 24` on InputDir/Label/PlotDir,
  layout functions, `FullSimMCSinglesFitFile`. `dr_correction_sample_cfg.h`: `mc_dir` →
  `sample_dir`, `overlay_year`, `DrCorrHistFile/SinglesFitFile/ClosureFile/PairRecoEffFile`,
  `DrCorrDataRefDir/Tag/Text` (D6), plateau/fit files under `mc_trig_eff/dr_correction/`.
- `FitMCSinglesEffcy.cxx`: PNG block (TCanvas/TLegend/gStyle, `mc_trig_eff_fit_plots/`) deleted;
  output `single_mu_effcy_pT_fit_mc_<label>...` (RW 3c closed); `overlay_year` last arg.
- `FillMCTrigEffHists.cxx`: `GetSampleConfig(sample, wp, overlay_year)` built on `GetDrCorrSample`;
  data refs via `DataFitTmpl/DataHistTmpl(products)`; output via `DrCorrHistFile` + mkdir.
- Year knob added (last positional argument, default 24) to: FillMCTrigEffHists, FitMCSinglesEffcy,
  plot_mc_trig_eff, plot_mc_trig_eff_corrected, fit_dr_corrections(+_all), plot_dr_correction_fits(+_all),
  combine_dr_correction_fits, plot_mc_singles_2d_effcy. pp-only macros (closure, pair-eff, tables,
  pthat stats) unchanged in signature.
- `PythiaAlgCoreT.h/.c` (`overlay_pbpb_year`), `RDFBasedHistFillingPythiaFullsimOverlay.cxx`,
  `PythiaFullsimRecoEffPlotter.cxx` (overlay child), `plot_single_muon_reco_effcy.cxx` (label now
  from the helper — it had a pbpb23 label hard-coded against a helper-resolved dir),
  `plot_muon_truth_q_eta_spectrum.cxx`, `plot_pythia_fullsim_overlay_kn_pt_crossx.cxx`.
- Shell: `pipelines/fullsim_sample_layout.sh` (new) sourced by 9 scripts; `OVERLAY_YEAR` plumbed
  through `run_pythia_fullsim_overlay*.sh`, `run_pythia_fullsim_overlay_condor.sh` (2nd arg),
  `run_pythia_fullsim_overlay.sub` comment, `pipeline_pythia_fullsim_overlay.sh`.
- Docs: root `README.md` (layout), `Analysis/README.md`, `docs/pythia_fullsim_overlay.md`,
  `docs/systematic_uncertainties.md`, `docs/placeholder.md`, `.claude/agents/plot-reviewer.md`,
  `.claude/commands/review-plot.md`.

### 2026-09-16 steps 7-9 (disk, verification, review)
- Disk move (`<scratchpad>/reorg_sample_dirs.sh`, manifest committed as
  `fullsim_sample_dir_layout_move_manifest_20260916.tsv`, 377 lines): 373 files moved by
  same-filesystem rename, md5 identical before/after, 0 sources left, 0 destinations missing;
  4 `mc_trig_eff_fit_plots/` dirs → `backup/mc_trig_eff_fit_plots/`. Per sample: full_sample
  hists 29 / dr_correction 128 / singles_fits 6 / pair_eff 3 / closure 9 / reco_eff 1 / backup 14
  / logs 4 / CSV 4; pbpb23 hists 29 / dr_correction 88 / singles_fits 6 / backup 1; noovl hists 16 /
  dr_correction 28 / singles_fits 2; test_sample hists 4 / singles_fits 2 / logs 1; pbpb24: empty
  skeleton only. Left at the roots by D7: NTUP `.bak_*` (grid_monitor-managed), `grid_monitor_*`,
  `lgd_farm.log`, `lgd_rules.txt`, `merging-record.txt(.bak_*)`, `local_dev/`.
- Compile: 29 (iteration 1) + 25 (iteration 2) + 3 (iteration 3) ACLiC builds, all OK; `bash -n`
  on 16 shell scripts OK; the shell twin resolves every product to an existing file.
- Smoke (`<scratchpad>/smoke_readers.C`): 63/63 inputs open at the new paths (pp_full T+M,
  overlay-23 T+M, noovl T+M, pp); `DrCorrectionCrossxEvaluator::Load` + `PairRecoEffEvaluator::Load`
  (the pp24 crossx loaders) succeed for pp_full Tight (eps_dR(0.3,20,0.2)=0.955689,
  eps_reco(20,0.2,0.3)=0.739951). pp_full **Medium** throws "pair-pT edge 0 is 8 in the fit file
  but 9 canonically" — PRE-EXISTING (`mc_trigger_efficiency.md` R34(d): Medium ΔR fits still on
  the 8 GeV axis), not caused here.
- End-to-end producers: `FitMCSinglesEffcy("pp", true)` and
  `fit_dr_corrections("pp_full", true, 3, "expo", false, "", "nocorr")` read from and wrote into
  the subtrees; the refit is bin-identical to the previous file (1118 bins, 0 differences);
  both original files restored byte-identically afterwards. `write_mc_pair_statistics_tables
  ("pp_full", true)` regenerated `plots/pp_trigger_efficiency/mc_statistics_pt8bins/` (6 CSVs;
  the directory had been archived under `old_4GeV_cut/` on 2026-09-15, so this is a fresh
  product on the current 9 GeV Step-3 file, not an overwrite).
- `/review-analysis-code` (log `.claude/logs/review-analysis-code-20260916-173143-fullsim-sample-dir-layout.md`):
  iteration 1 FAIL — 3 CRITICAL (fit-stage `_step3.root.root` from `MakeStepCfg`; unbound
  `${DATA_ROOT}` in `pipeline_pythia_fullsim_pp.sh` Stage 10; `write_mc_pair_statistics_tables`
  out_base parser vs the overlay year leaf) + 2 WARNING (AMI dir of test samples always the
  truth production's → own `ami_info/` preferred, r17864 registered as `ami_weights.md` B2;
  three overlay consumers without a year knob) + 3 INFO — all fixed; iteration 2 FAIL — 1
  WARNING (kn-crossx headline literal "Pb+Pb 2023") + 3 INFO — fixed; iteration 3 PASS.
- Reviewer's D6 assessment: for `overlay_year = 23` every data-reference path is byte-identical
  to before (D2 unchanged); for 24 it is a new like-for-like choice (pbpb24 data 0-5 %) that the
  user should ratify before the first pbpb24 Step-1 comparison is read.

## Results & Observations

- R1. The layout is now ONE table (`FullSimSampleType.h` + the shell twin); `DrCorrSample.mc_dir`
  no longer exists, so any stale consumer fails to compile instead of reading a wrong path.
- R2. Pre-existing, untouched: pp_full Medium ΔR fits are stale (8 GeV axis) → the crossx
  evaluator throws for Medium; the pp24 TEST sample has no plateau file (never re-measured since
  the plateau file was introduced), so `fit_dr_corrections("pp", ...)` cannot run there.
- R3. The pbpb24 overlay sample dir has only the NTUP + `ami_info/` + empty subtree skeleton;
  the NTP on it needs `allow_missing_slices` (single slice) — loud, by design.

## Remaining Work

- User ratification of D6 (pbpb24 overlay ↔ Pb+Pb 2024 data 0-5 % reference).
- Not this task: Medium ΔR refit (R2), pbpb24 overlay NTP/analysis (blocked on the r-tag question
  in `hijing_overlay_pbpb24_test_sample_skim.md`).

## Latest Stage

DONE 2026-09-16. Code + disk + docs complete, review PASS (3 iterations), committed. Open for the
user: ratify D6. Doc closed.
