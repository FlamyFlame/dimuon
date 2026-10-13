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

**Rule recorded (user 2026-10-12):** diagnosis plots are exempt from `/review-plot` (`.claude/CLAUDE.md` §Auto-Dispatch Rules).

## Latest Stage
2026-10-12: R1 replotted (522336/522355 waiting for good-LB lumi); R2a–d delivered. Open: (1) good-LB lumitable from the updated GRL → redraw pbpb26; (2) rerun pbpb26 per-run + pair plots after the Step-17 NTuple hadd; (3) experts/user: 463315 skim bookkeeping, 2026 ZDC-time shift in 522041/522200.
