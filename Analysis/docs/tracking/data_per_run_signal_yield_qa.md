# Per-run signal-region yield / luminosity (abnormal-run QA)

## Objective
For pp24 and Pb+Pb 2023–26 data, plot (# signal-region events) / (integrated luminosity) per run to expose
abnormal runs. One PNG per dataset (5 total). Evidence only — **no run-exclusion recommendation**
(memory `feedback_run_exclusion_needs_experts`).

## Procedure (as run, 2026-10-08)
- Input: nominal `muon_pairs_*_mindR_0_02.root`, tree `muon_pair_tree_sign2` (= opposite sign); Pb+Pb NTuples
  already carry the 5-cut event selection (applied at NTuple stage, `PassEventSelHook`).
- Cuts (verbatim from `ParamsSet`, same as the crossx fill): `pair_pass_tight` (Tight WP),
  `SignalMinvCutExpr` (1.08–2.9), `SignalPairPtCutExpr` (>9), `FiducialGapCutExpr` on both muons,
  `PairFiducialEtaCutExpr` (|η^pair|<2.2); Pb+Pb also centrality acceptance 0 ≤ avg_centrality < 80.
- Yield: RAW, unweighted, counted as distinct events (run, LB, BCID, entry); pairs differ by <0.3 %.
- Luminosity: `Prescale Corrected` column of the analysis-trigger table in `IntNotes/data/luminosity/`
  (pp24 `HLT_2mu4_L12MU3V`, pb⁻¹; Pb+Pb `HLT_mu4`, µb⁻¹ → nb⁻¹ ÷1000).
- Not drawn: runs with zero lumi (461633) and `PbPbBadRuns()` (2023: 461674, 462964). Data and lumi run lists
  otherwise match exactly for all five datasets (pp24 8 runs; PbPb23 58; 24 14; 25 42; 26 35).

## Code / outputs
- `plotting_codes/data_qa/per_run_signal_yield_table.cxx` (RDF per-run counts), `plot_signal_yield_per_lumi_per_run.cxx` (plots); commit `c296228`.
- `~/usatlasdata/dimuon_data/plots/data_qa/signal_yield_per_lumi_per_run/` — PNGs + `*_yield_per_lumi.csv` (N, L, N/L, err, pull vs dataset average).
- Rerun: `root -l -b -q` a macro that `#include`s the plot .cxx and calls `plot_signal_yield_per_lumi_per_run()`.

## Results (|pull| > 4 vs dataset average N/L)
| Dataset | Run | N/L vs avg | Pull |
|---|---|---|---|
| pp24 | 488239 (N=699, L=0.26 pb⁻¹) | ≈2.1× | +20σ |
| PbPb23 | 463315 | ≈0.43× | −23σ |
| PbPb24, PbPb25 | none | | |
| PbPb26 | 522384 | ≈0.48× | −24σ |
| PbPb26 | 522408 | ≈0.50× | −30σ |
| PbPb26 | 522200 | ≈0.70× | −12σ |
| PbPb26 | 522298, 522628, 522925, 523138 | ≈1.07–1.09× | +4–5σ |

522384/522408 are the runs already covered by `pbpb2026_runs_522384_522408_ntrk_lowband.md`.
Causes of the others are NOT investigated.

## Remaining work
Hand the plots/CSVs to the experts; investigate 488239, 463315, 522200 only if requested.
`/review-plot` explicitly waived by the user (2026-10-08).

---
# Round 2 (REOPENED 2026-10-12): 2026 exclusion + outlier-run investigation

## Requests (user, 2026-10-12)
R1. Remake the PbPb26 per-run plot with 522384/522408 dropped and only the good LBs of 522336 (drop LB 344–521)
    and 522355 (drop LB ≥ 448) (user ruling 2026-10-12, single source `PbPbBadRuns`/`PbPbBadLBRanges` in
    `PbPbEventSelConfig.h`). **User: WAIT for their new HLT_mu4 lumitable (made from the updated GRL), and go ahead
    with the other tasks meanwhile.** Also needs the 2026 NTuple rerun with the LB exclusion, run by the
    pbpb2026_analysis_support Step-17 session (`pipelines/run_pbpb26_badlb_evsel_ntuple_20261012.sh`).
R2. Outlier runs: PbPb23 463315, PbPb26 522200, PbPb26 522041 (added by the user); pp24 488239 stays on the list.
  a. Per run: standalone event-selection PNGs (one per cut, cut lines drawn) under
     `plots/data_qa/<dataset>_run<run>/evsel/`, plus the % of the run's events that pass all event-selection cuts.
  b. Code check: does an event that fails any event-selection cut contribute zero muon pairs?
  c. Which bad PbPb 2023 runs with ZDC problems are skipped?
  d. Muon-pair distributions (pair pT, η, mass) plus other basic checks, outlier run vs the rest of its dataset.

## Plan
- R2a: new `plotting_codes/data_qa/plot_evsel_single_run.cxx`. It reads raw NTUP `HeavyIonD3PD` and requires
  `b_HLT_mu4_L1MU3V` (same as the cut derivation). It applies the year's DERIVED cuts from
  `event_sel_cuts_pbpb_20YY.root`, mirroring `PbPbExtras::PassEventSel` exactly: bad run/LB first, then cuts 1–5
  in sequence, with the per-run cut-3 thresholds when the TTree is present. Each standalone plot shows the stage
  BEFORE that cut, as in `plot_pbpb_event_sel_cuts.cxx::DrawStandaloneCut`. Per-run pass fractions are computed
  for ALL runs of the year so the outlier can be compared with the rest (also 0–80 %). 2026 waits for the
  Step-17 re-derived cuts file.
- R2d: from the nominal NTuple output (`muon_pairs_*_mindR_0_02.root`, sign2 = OS, Tight), comparing the outlier
  run with the rest of its dataset, shape-normalised: pair pT, η^pair, m_μμ, single-muon η/φ/pT, centrality,
  N_signal/N_HLT-event per run (independent of lumi), and signal yield per LB. No re-derivation from raw.
  2026 waits for the Step-17 NTuple rerun (the file is overwritten in place).
- Evidence only, no run-exclusion recommendation.

## Results R2 (2026-10-12)
**R2b — a failed event contributes NO pairs (verified in code).** `DimuonDataAlgCoreT.c::ProcessDataHook` (≈l.682):
`if (!this->PassEventSelHook()) continue;` runs BEFORE the pair loop. For Pb+Pb, `PassEventSelHook` → `PbPbExtras::PassEventSel`
(bad run/LB, then cuts 1–5, any failure → `false`). So none of that event's pairs reach the trees (OS or SS) or the
cut-acceptance histograms. This is the same code path for the nominal (`_mindR_0_02`) and trig-eff (`_res_cut_v2`)
NTuples. Downstream, the RDF only adds the centrality acceptance [0,80). Note: unlike the cut derivation, the NTuple stage
does NOT require `b_HLT_mu4_L1MU3V` at event level; the trigger enters through the per-pair trigger match.

**R2c — 2023 runs with ZDC problems, skipped:** 461674 and 462964. They are in `PbPbBadRuns(23)`
(`PbPbEventSelConfig.h`, commit `1e2b1d6`, 2026-05-05), excluded from the cut derivation, from the NTuple stage, and from the
lumi (1.18365 → 1.17576 nb⁻¹, `IntNotes/data/luminosity/README.md`). Diagnostic plots:
`plotting_codes/event_selection/plot_zdc_diag_groups_pbpb23.cxx` (ZDC energy, time A vs C, preamp, ZDC vs FCal: these runs vs
the rest). ⚠ Several docs call them "b-hadron runs" (README, analysis_metadata.md, roadmap, status summary). That label is
wrong; they are the ZDC-problem bad runs.

**Infrastructure note:** `data_pbpb23_part2.root` is a symlink into `/pnfs` (LGD), and NFS reads of it hang in D state
(another session's `GetEntries` process has been stuck for 5 days). Read it over xrootd instead:
`root://dcgftp.usatlas.bnl.gov:1096//pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/9a/b4/data_pbpb23_part2.root`.
All 46 216 614 entries read OK this way.

**R1 — PbPb26 replotted (2026-10-12)** with `plot_signal_yield_per_lumi_per_run("pbpb26")`. The NTuple is still the
2026-09-17 one, and the new bad-run/LB guard in `per_run_signal_yield_table.cxx` (`PbPbIsBadRunLB`) removes 522384/522408
and the bad LBs (2.85 % of OS pairs). 522336 and 522355 are NOT drawn: their N covers good LBs only, but the only 2026
lumitable (`lumitable_pbpb_26_HLT_mu4.csv`, original GRL) has full-run lumi. The updated GRL
(`SkimCode/xmls/physics_HI2026_50ns_noIBL_updated.xml`, in place) lists LB ranges but has no luminosity. To draw them, set
`lumi_has_good_lb_only=true` with the HLT_mu4 lumicalc table made from that GRL. Result: the dataset average is
3.497×10⁴ nb (up from the low-band runs dropping out); remaining outliers are 522200 (pull −13.9) and 522041 (≈0.57×, ≈−5σ).
Rerun once the Step-17 NTuple (pbpb2026_analysis_support) is hadd'ed: the re-derived cuts change the counts slightly.

**R2a/R2d — the per-run N/L deficit splits into three factors** (`run_norm_diagnostics.cxx` →
`plots/data_qa/signal_yield_factored_per_run/pbpb{yy}_norm_factors_per_run.png` + the two stand-alone factor plots `pbpb{yy}_signal_per_evsel_event_per_run.png`, `pbpb{yy}_evsel_events_per_lumi_per_run.png`, CSV `pbpb{yy}_run_norm_factors.csv`; user 2026-10-12: factored plots in their OWN subdir, the original N/L plots stay in `signal_yield_per_lumi_per_run/`):
N_sig/L = (N_HLT/L) × (N_evsel/N_HLT) × (N_sig/N_evsel).
| Run | N_HLT/L vs neighbours | evsel pass (HLT_mu4; 0–80 %) | N_sig/N_evsel pull | ⇒ deficit sits in |
|---|---|---|---|---|
| 463315 (2023) | **0.43×** (4.1e7 vs ≈9.5e7 nb) | 99.65 % (year 99.63 %) | −0.7 | events/lumi BEFORE any analysis cut: skim incomplete for this run, or wrong lumi row |
| 522200 (2026) | 1.00× | **65.2 %** (year 98.5 %) | +0.5 | event selection, cut 2 (ZDC time) |
| 522041 (2026) | 1.01× | **72.3 %** | −2.0 | event selection, cut 2 (ZDC time) |
- 2026 cut 2: the ZDC collision peak is shifted off-centre, (t^A, t^C) ≈ (+1.0, +0.7) ns in 522041; 522200 has a second core at
  t^A ≈ +1.3 ns (the Step-16 "t_A ≈ +1.3 ns cluster"). The fixed |t| < 1.5 ns box cuts into collisions. Only these two runs have
  cut-2 survival < 97 % (2023 none). 522041 is flat vs LB (72 %). 522200 changes in steps: LB 105–164 ≈ 75→61 %,
  165–244 ≈ 78 %, ≥ 245 ≈ 45 %. In 522200 the loss DEPENDS ON CENTRALITY: signal pairs 0–10 % ≈ 0.75×, 30–50 % ≈ 1.35× the
  other runs (χ² p = 2×10⁻⁹). Other pair shapes are consistent.
- 463315: every cut and every pair shape is normal (all p > 0.09); only the event count per lumi is low. Next check: the skim's
  input-file / event bookkeeping for run 463315 (AMI/rucio event count of the run's HardProbes input vs the events in
  `data_pbpb23_part*.root`).
- 462759 (2023): N_HLT/L 0.74× (pull −237, stat only), small lumi, so the N/L plot didn't flag it. Same signature as 463315.
- pp24 488239: low-lumi run (139 LBs, 0.26 pb⁻¹ ≈ 1/40 of the normal lumi/LB); pair shapes show no systematic difference and
  signal/all-OS-Tight = 0.207 vs 0.211 ⇒ normalization (lumi or trigger) issue. No raw-event normalization was done for pp.
- Event-selection per run (cut lines = the year's derived cuts, applied as in PassEventSel): `plots/data_qa/pbpb{23,26}_run<run>/evsel/`;
  per-run pass table for all runs: `plots/data_qa/signal_yield_factored_per_run/pbpb{yy}_evsel_pass_per_run.csv`. Pair plots: `.../pairs/`
  (χ² "UU"; single-muon panels count 2 entries per pair, so their p-values are approximate).
- 2026 per-run evsel used the Step-17 re-derived `event_sel_cuts_pbpb_2026.root` (19:58); the 2026 pair plots used the
  2026-09-17 NTuple (these runs are not affected by the LB exclusion).
- Evidence only — NO run-exclusion or cut-change recommendation (memory feedback_run_exclusion_needs_experts).

**Factored plots, all years (2026-10-12)**, in `plots/data_qa/signal_yield_factored_per_run/`:
- **Definitions.** N_HLT_mu4 = raw-skim events with event-level `b_HLT_mu4_L1MU3V`, good runs/LBs, BEFORE the 5 event-selection
  cuts. The skim has already applied: ≥ 1 configured trigger, original GRL, no detector error flags, ≥ 1 vertex, required ZDC data
  present. N_evsel,0–80% = of those, passing cuts 1–5 with centrality [0,80): no muon-level or pair-level cut. N_sig/N_evsel
  therefore holds every muon/pair effect, and N_evsel/L holds trigger, skim/DAQ completeness, lumi and event-selection efficiency.
- **PbPb24, PbPb25:** no deviating run (all within 8 % of the median N_HLT/L, evsel > 95 %, |sig/ev pull| < 3). Year evsel pass:
  99.30 % (2024), 97.20 % (2025). (2023: 99.63 %, 2026: 98.49 %.)
- **pp24** (user choice: N_HLT from the unprescaled analysis trigger HLT_2mu4_L12MU3V; HLT_mu4 is prescaled ×400–870, which would
  make N_sig/N_HLT_mu4 scale with the prescale). New code: `hlt_events_per_run.cxx` (raw-NTUP count; `/pnfs` read over xrootd)
  and `pp_norm_factors.cxx` → `pp24_hlt_events_per_lumi_per_run.png`, `pp24_signal_per_hlt_event_per_run.png`,
  `pp24_run_norm_factors.csv`. **488239: N_HLT_2mu4/L = 5.80×10⁴ pb vs 2.82×10⁴ pb for the full runs (2.06×); N_sig/N_HLT =
  0.0463 vs ≈0.0457, normal.** So the 2× N/L excess sits entirely in triggered events per lumi (lumi or trigger configuration of this
  low-lumi run), not in the analysis. 488210 (0.04 pb⁻¹): 1.05×, consistent.
- Shared plot style: `per_run_factor_plot.h` (`DrawPerRunFactor`).

**463315 single-muon η–φ** (`muon_eta_phi_run_vs_rest.cxx` → `plots/data_qa/pbpb23_run463315/muons/`). Muons taken from the pair NTuple
(no single-muon tree), each muon counted once per event; panels: run, other runs, ratio. Grid 0.25 × π/8 (QA display axis).
(a) muons of all Tight pairs OS+SS: 21 509 muons; (b) muons of OS Tight signal pairs: 1 348. χ²/ndf = 373/319 (p = 0.019) for both. The
near-identical χ² was checked: the histograms differ and a manual χ² matches `Chi2TestX`, so it is a coincidence. The ratio map is flat,
with NO localized hole or deficit, only a mild excess spread over the map. Consistent with 463315 being normal in everything except events/lumi.

**Rule recorded (user 2026-10-12):** diagnosis plots are exempt from `/review-plot` (`.claude/CLAUDE.md` §Auto-Dispatch Rules).

---
# Round 3 (2026-10-12): skim bookkeeping of PbPb23 463315 / 462759

## Requests (user, 2026-10-12)
R3a. Find the skimming logs / grid tasks for 463315 and 462759 and check whether something went wrong in the skim.
R3b. Submit a small grid re-skim of just these two runs (`SkimCode/run_23hi/`).

## Plan
Task IDs from `~/usatlasdata/dimuon_data/may2026_skim.txt`; BigPanDA task/job JSON per input dataset; per-LB skim counts
vs GRL; download/hadd record in `grid_monitor_status.log`; per-file test = read EventInfo (LB, eventNumber, detector
error flags) of sampled AOD files and count what fraction of each file's events is in the skim.

## Findings R3
- Tasks: 462759 → 50267236 (part2), 463315 → 50267250 (part3). JEDI: every input file of both runs `finished`
  (462759 393/393, 1 818 081 ev; 463315 2302/2302, 10 456 893 ev); no failed/abandoned files (the 78 abandoned files of
  part2 are all run 462969 → recovery part5). Grid output datasets are now empty (expired).
- Download + hadd lossless: 210/210 merged part3 files downloaded, pre-merge entries = 44 525 454 = `data_pbpb23_part3.root`.
- Lumi and skim GRL agree (463315: 402 LBs, 462759: 76 LBs, same in GRL v120-pro33-03 and the HLT_mu4 lumitable).
- AOD events / lumi is NORMAL for both runs (0.248 and 0.237 ×10⁹/nb vs 0.23–0.26 for the neighbouring runs); the deficit
  is skimmed HLT_mu4 events / AOD event: 463315 0.165, 462759 0.296 vs 0.35–0.41 for neighbours.
- Per-LB skim counts: 101 of 402 GRL LBs of 463315 and 14 of 76 of 462759 have ZERO skimmed events; present LBs have
  quantized-looking counts (~1 000 / ~3 500 / ~7 000 / ~19 000). Reference runs 463308 / 462763 have no such gaps.
- Skim code (`fe3b18b` `TrigRates::execute`): events are silently dropped (no job error) by `CleaningCuts`
  (errorState==Error in Pixel/SCT/TRT/Muon/LAr/Tile/Core, or Core flag bit 18) or by failing all listed triggers.
- R3b submitted: jediTaskID **53029101**, `user.yuhang.TrigRates.dimuon.PbPb2023data.Oct2026.v1.rerun463315_462759.`
  (`run_23hi/grid_sub_rerun_463315_462759.sh`, `InDstxt_PbPb2023data_rerun_463315_462759.txt`; HEAD source, same GRL/mode;
  local 200-event test OK). Name has NO `partN` on purpose: never auto-merge into `data_pbpb23_partN` (double counting).

- **ROOT CAUSE (R3a): JEDI `--mergeOutput` left finished job outputs unmerged.** Per-file test (EventInfo of sampled
  AODs vs the skim's eventNumbers): each AOD file is all-or-nothing, either 0.33–0.44 of its events in the skim (the
  normal HLT_mu4 fraction) or exactly 0. 463315: 37/61 files (61 %) absent, 0 detector-error-flagged events in any file.
  462759: 10/40 (25 %) absent. Expected N_HLT/AOD = 0.39×0.39 ≈ 0.15 and 0.40×0.75 ≈ 0.30, observed 0.165 and 0.296, so
  the missing files account for the whole deficit. Example: AOD `_000550/_000551` → job 7139059832 (DESY-HH, finished,
  exit 0, 2.95 MB output `panda.um...EXT0._002502`). That output has status `ready` and NO pmerge job consumed it, whereas
  the good job's output `_002354` was merged by pmerge 7140576483 into `EXT0._000053`. Counting the inputs of all merge jobs
  (210 + 228): part3 merged 9 666 of 10 378 finished job outputs (**712 orphaned, 6.9 %**), part2 9 506 of 9 568 (62).
  Confirmed `ready`-never-merged outputs so far: 109 (part3) + 49 (part2). 608 more part3 lookups were rate-limited by
  BigPanDA; not yet attributed to runs. The skim CODE is not at fault (the jobs ran and wrote output); this is a merge-step
  loss on the grid, invisible to task status (`done` 100 %) and to the download/hadd check.
- Caveat for 53029101: it also uses `--mergeOutput`, so it can lose files the same way. Check its per-LB completeness
  against the AOD before trusting it.
- Possible other victims: 463380 (N_HLT/AOD 0.350, N_HLT/L 0.91× the neighbours) and 463364 (0.367), both part3. The
  other May-2026 tasks (2023 parts 1/4, 2024, 2025, pp24) have NOT been audited for orphaned outputs.
- Scratch evidence (session scratchpad, not kept): task/job JSONs, `match_out.txt`, `perlb_*.csv`.

# Round 4 (2026-10-13): recover the lost 463315 / 462759 events

## Requests (user, 2026-10-13)
R4a. List the other May-2026 skim tasks.
R4b. Replace these runs' events in part2/3, or (easier, since part1-3/5 are rucio files on BNL LOCALGROUPDISK) add the
missing events. Be very careful to avoid duplication.

## Plan R4b
Rucio files on LGD cannot be rewritten, so the lost events go into a NEW supplementary part `data_pbpb23_part6.root`:
1. Wait for 53029101 to finish (background poller, BigPanDA every 20 min). Download its merged outputs to a separate dir.
2. Completeness of the re-skim: BigPanDA merge-input count = finished trf outputs (no orphans); every May event of the two
   runs is in the re-skim (May set is a subset); per-LB: no GRL LB empty; N_HLT/AOD back to the 0.35-0.41 band.
3. Supplement = re-skim events of runs 463315/462759 whose (RunNumber, eventNumber) is NOT in part3/part2 (event-level
   dedup, independent of which job outputs were orphaned). Drop the new branch `muon_match_L1MU3V` so the schema equals
   part1-5. Abort if the re-skim has internal duplicates.
4. Verify: union (part2/3 + part6) of each run has no duplicate eventNumber and equals the re-skim set exactly.
5. Wire part6 into the PbPb23 part count (`file_batch_max{23}`, `queue N` in the run_pbpb_23*.sub, plotting part lists).

## Findings R4
- Builder: `SkimCode/scripts/make_pbpb23_part6.C` (tested on fakes: dedup 58→30 May+28 new, union unique; aborts on
  re-skim duplicates and on May events missing from the re-skim, writing nothing).
- R4a: May-2026 = 15 tasks (`may2026_skim.txt`): PbPb23 50270743 (part1, v2), 50267236, 50267250, 50267265; PbPb24
  50267371, 50267396; PbPb25 50267423, 50267435, 50267444, 50267472, 50267482, 50267490; pp24 50267518, 50267529,
  50267539. Plus the Sep-2026 recovery of the abandoned files (`sep2026_may26_recovery.txt`): 52569154 (PbPb23 part5),
  52569156 (PbPb25 part7), 52569161 (pp24 part13). All used `--mergeOutput`; none except part2/3 audited for orphans.
- Re-skim (HEAD TrigRates) vs May skim, 58 overlapping events of 463315 LB 669 (AOD `_000003`): all 168 common branches
  identical element by element (23 752 values, TTreeFormula); the re-skim has one extra branch `muon_match_L1MU3V`.
  Skim files hold only the `HeavyIonD3PD` tree. May event lists have no duplicate eventNumbers (462759: 537 740,
  463315: 1 728 372). The two runs occur only in part2/part3 input lists (not part1/4, not the part5 file list).
- 53029101 `done` 05:38 (2695/2695 files, 0 failed). Merge audit: 1 247 trf outputs (`finished`) and all 1 247 are
  inputs of the 34 pmerge jobs (no orphan); the other 3 "finished trf" jobs are build jobs (3 lib inputs); the 17 gaps in
  the EXT0 serials are outputs of failed (retried) attempts. Output 34 files, 5.165 GB, 4 784 435 events, downloading to
  `~/usatlasdata/dimuon_data/pbpb_2023/rerun_463315_462759/` (subdir: not matched by `data_pbpb23_part*.root`).
- Download complete (34/34 files, 5 164 654 463 B). Re-skim completeness per LB (`perlb.C`): 463315 402/402 GRL LBs
  filled, 0 outside GRL, 4 060 027 events, N_HLT/AOD 0.388; 462759 76/76, 724 408 events, 0.398 (neighbour band
  0.35-0.41). No LB has fewer events than in May; 91 (463315) / 47 (462759) LBs are identical to May.
- part6 built (`make_pbpb23_part6.C`): 462759 re-skim 724 408 = May 537 740 (all found) + NEW 186 668; 463315 4 060 027
  = May 1 728 372 (all found) + NEW 2 331 655. `data_pbpb23_part6.root` = 2 518 323 events, 168 branches.
- Independent verification (`verify6.C`, May lists from the earlier `dumpev.C` dumps, not from the builder): part6 only
  runs 463315/462759; 0 duplicates in part6, the re-skim and May; part6 ∩ May = 0; every part6 event is in the re-skim
  with the same LB; May + part6 = re-skim exactly (724 408, 4 060 027); 168/168 part3 branches with the same class.
- Installed as `~/usatlasdata/dimuon_data/pbpb_2023/data_pbpb23_part6.root` (local file, not LGD). The raw re-skim
  stays in `pbpb_2023/rerun_463315_462759/` (no glob in the code reaches into it; it must never become a part).
  Bookkeeping: `data-merging-record.txt`, `oct2026_pbpb23_rerun_diag.txt`.
- Wiring 5 → 6 parts (sites found by two independent read-only agents, identical lists): `PbPbExtras.c` `{23, 6}`;
  `queue 6` in the 6 live `run_pbpb_23*.sub` (not `_part4_only`); `QUEUE_COUNTS[23]=6` in `pipeline_pbpb_crossx.sh` /
  `pipeline_pbpb_trig_eff.sh`; `ScrambGen.h` `case 23: return 6`; part6 added to the file lists of
  `plot_pbpb_event_sel_cuts`, `_event_level`, `plot_pbpb_fcal_comparison`, the 3 `plot_zdc_preamp_*`,
  `plot_overlay_rtag_compr_pbpb_data`; nparts 6 in `plot_npairs_vs_centrality`, `plot_muon_q_eta_spectrum`,
  `plot_muon_q_eta_pt_dependence`; loops `p <= 6` in the 3 one-off `plot_zdc_*_pbpb23`. `preflight_pbpb_year.sh 23`:
  ok=25, count-mismatch=0, 6 parts on disk; all 13 edited macros pass `g++ -fsyntax-only`. Docs: `pbpb_pipelines.md`
  job table, `SkimCode/README.md` skim table.
- Lumi unchanged: the runs' full GRL lumi was already in the lumitable / `PbPbSampledLumi`; part6 only restores events.
- NOT done (needs user go-ahead): PbPb23 NTuple processing for part6 and everything downstream (event-selection cuts
  derived from the 2023 part list, PbPb23 crossx/trig-eff hadd, combined PbPb crossx, per-run QA replot).

# Round 5 (2026-10-13): PbPb23 rerun with part6

## Request (user, 2026-10-13 ~16:30)
Rerun the PbPb23 NTuple processing and the other PbPb23-only steps, and update the PbPb23 data-QA plots. Two other
agents work on the pbpb26 NTuple process and the MC reco efficiencies; combined crossx + R_AA wait for all three.

## Done = (all PbPb23-only; nothing combined)
1. Event-selection cuts re-derived from parts 1-6: `event_level(23)` → `cuts(23)` → `cuts_alt(23)`
   (`event_sel_cuts_pbpb_2023{,_alt}.root`, cut-2 keys present); old files backed up.
2. Condor NTuples (6 parts) `_mindR_0_02_res_cut_v2` (trig-eff, `run_pbpb_23.sub`) + `_mindR_0_02` (crossx,
   `run_pbpb_23_nominal.sub`); validated (ExitCode 0, clean .err, fresh non-empty parts, run/LB exclusion check) + hadd;
   463315/462759 OS counts up; old combined files in `pbpb_2023/pre_rerun_20261013/`.
3. Trig-eff 23 (pipeline_pbpb_trig_eff stages 5-8 for 23 only): Tight RDF fill, Fermi+log fit, TF1/TH2D validation,
   `trig_effcy_plot_PbPb(23)`; Medium WP RDF fill + fit (run_data_trigeff_medium_wp.sh commands for 23).
4. Crossx RDF fill 23 (run_crossx_hist_filling_pbpb23.sh commands + freshness/filled-histogram probe) and the OS−SS
   twin (`run_crossx_hist_filling_pbpb_ss_subtr.sh 23`).
5. PbPb23 data-QA: per-run N/L and factored N_HLT/L, N_evsel/N_HLT, N_sig/N_evsel (`qa_pbpb.sh` recipe, 6 parts).
Concurrency: every ACLiC compile in a private build dir (other sessions compile the same classes); no combined-year
stage (Stage 0 FCal all-year comparison, crossx Stages 6-8, R_AA) is run. Not in scope (user decision 2026-10-13 for the
rerun: nominal + res_cut_v2 only): `_no_res_cut`, single-muon trees/ScrambGen, template fits.

## Latest Stage
2026-10-13 R5: PbPb23 rerun with part6 started (plan above); step 1 evsel.
2026-10-13 R4: DONE through step 5 — part6 (2 518 323 lost events) built, verified, installed, wired (preflight OK).
Open: user go-ahead for the PbPb23 NTuple/downstream rerun; audit of the other 16 May/Sep-2026 tasks for merge orphans.
2026-10-12 R3: root cause = JEDI merge orphaned finished outputs (part3 712 jobs, part2 62). Grid re-skim 53029101
submitted. Open: audit every May-2026 task for orphaned outputs; recover the missing files of 463315/462759 (e.g. from
53029101 once checked complete) → user decision on how to integrate (replace runs in part2/3, no double counting).
Previous: R1 replotted (522336/522355 waiting for good-LB lumi); R2a–d delivered. Open: (1) good-LB lumitable from the updated GRL → redraw pbpb26; (2) rerun pbpb26 per-run + pair plots after the Step-17 NTuple hadd; (3) experts/user: 463315 skim bookkeeping, 2026 ZDC-time shift in 522041/522200.
