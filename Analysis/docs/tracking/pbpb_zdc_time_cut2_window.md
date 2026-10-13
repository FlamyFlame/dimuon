# Pb+Pb event-selection cut 2: per-year, per-side ZDC time window

## Objective
Replace the symmetric ZDC time cut |t^A|, |t^C| < 1.5 ns (cut 2 of the 5-cut Pb+Pb event selection) with a
per-year, per-side window, and change the A-side edges as decided by the user (2026-10-12, **PROVISIONAL**: the user
checks the event-selection plots to finalize cut 2 before it is propagated to the NTuples):
- 2025: t^A lower edge −1.5 → **−0.6 ns**
- 2026: t^A upper edge +1.5 → **+0.8 ns**
- 2023, 2024 and every t^C edge stay at ±1.5 ns.

## Physics Procedure
1. **Motivation.** The ZDC time distributions have clusters beside the in-time collision peak at (0, 0) ns: 2025 has a
   t^A ≈ −1.4 ns shoulder (12.4 M events) and a −3 ns blob; 2026 has a t^A ∈ [1, 2] ns cluster (2.6 M events = 1.1 % of the
   main box). All sit at nominal vertex z, so they are not main–satellite collisions (those would sit at c·Δt/2 ≈ 200 mm)
   and no vertex-z cut can remove them (pbpb2026_analysis_support.md Step 16 item 3). The symmetric ±1.5 ns box keeps about
   half of each. The new A-side edges remove them.
2. **Equation.** Cut 2 passes an event iff `a_lo < t^A < a_hi  AND  c_lo < t^C < c_hi` (strict inequalities, as the
   old `|t| < 1.5`), with ZDC index [1] = A and [0] = C (memory reference_zdc_sides).
3. **Method.** The window is a fixed per-year constant (NOT fitted): single source `PbPbZdcTimeCut(yr)` in
   `NTupleProcessingCode/PbPbEventSelConfig.h`. The derivation (`plot_pbpb_event_sel_cuts.cxx`) applies it in the
   sequential chain (cut 1 → cut 2 → cut 3 …), so cuts 3–5 (per-run preamp μ+7σ for 25/26, N_trk Gaussian-slice fits) are
   RE-DERIVED on the new post-cut-2 population (user). The derivation writes the four edges into
   `event_sel_cuts_pbpb_20YY{,_alt}.root`; the NTuple stage (`PbPbExtras::PassEventSel`) reads them from there.
4. **Negative constraints.**
   - Do NOT change t^C, cut 1 or the cut 3–5 procedures; do NOT change 2023/2024 (their files are not re-derived).
   - Do NOT rerun NTuples / RDF / crossx in this step (user: derivation + plots only until cut 2 is final). Until then the
     2025/2026 NTuples were made with the OLD cut 2 and do not match the new cuts files.
   - Do NOT keep writing the old symmetric key `ZDC_time_cut_ns`: a reader that still uses it would silently apply ±1.5.
     Readers must use the four edges (and throw if missing).

## Context
- Previous Pb+Pb event-selection definition: memory project_pbpb_event_selection; `docs/data_analysis.md`.
- 2026 bad runs/LBs (D7) are applied before any step; unchanged here.
- The Step-17 driver (`pipelines/run_pbpb26_badlb_evsel_ntuple_20261012.sh`) hung in Stage 0 on a `/pnfs` NFS read of
  `data_pbpb23_part2.root` and was stopped (user: this session takes over 2026). Its 2026 NTuple stage is on hold until
  cut 2 is final.

## Design Decisions
- D1. Window constants live in the config header, not in the cuts file. The cuts file records what was applied
  (provenance), and `PassEventSel` throws if the file and the header disagree (catches a stale file).
- D2. The old key `ZDC_time_cut_ns` is removed from newly written files. 2023/2024 files (not re-derived) still carry
  only the old key, so readers fall back to it ONLY if the four edge keys are absent AND the header window is symmetric
  and equal to it. Otherwise they throw.
- D3. `/pnfs` inputs of the chain are read over xrootd: NFS reads of `data_pbpb23_part2.root` hang in D state.

## Implementation Plan
1. Config header: window struct, `PbPbZdcTimeCut(yr)`, `PbPbPassZdcTime`, 4 keys, a reader helper. — per §2/§3
2. Derivation macro: apply, draw, write. — §3
3. `PassEventSel` + all other readers of the cut-2 key. — §3, D2
4. xrootd path resolution for `/pnfs` inputs. — D3
5. Back up the old 25/26 cuts files; rerun event_level / cuts / cuts_alt (25, 26) + the 4 FCal comparisons. — §3
6. Report before/after cut 2 survival and cut 3–5 values; `/review-analysis-code` on steps 1–4.

## Progress Log
### 2026-10-13 — Steps 1–4 DONE (code), step 5 running
- `PbPbEventSelConfig.h`: `PbPbZdcTimeWindow`, `PbPbZdcTimeCut(yr)` (23/24 ±1.5; 25 a_lo −0.6; 26 a_hi +0.8; other years
  throw), `PbPbPassZdcTime`, `PbPbWriteZdcTimeWindow` (4 keys `ZDC_time_{A,C}_{lo,hi}_ns`, deletes legacy key),
  `PbPbReadZdcTimeWindow` (checks file == header; legacy symmetric key accepted only for a symmetric header window),
  `PbPbRawInputPath` (/pnfs symlink → `root://dcgftp.usatlas.bnl.gov:1096/` + target, door override `PNFS_XROOTD_DOOR`).
- Cut-2 sites switched to the window: `plot_pbpb_event_sel_cuts.cxx` (application in the preamp-per-run pre-pass and in the
  sequential chain, writer, the two box drawings + label), `PbPbExtras.{h,c}` (`PassEventSel`), `plot_pbpb_fcal_comparison.cxx`,
  `fill_pbpb_evsel_run_lb_diag.cxx`, `plot_zdc_preamp_{cut_over_mean,cut_vs_zdcamp,gauss_fit,per_run_pbpb25}.cxx`,
  `overlay_rtag_checks/plot_overlay_rtag_compr_pbpb_data.cxx`, `data_qa/evsel_single_run.cxx`. No other reader of ZDC time
  applies a cut (event_level, run_lb_diag plotter, zdc_diag/run_dependence 2023 macros only histogram it).
- Raw inputs via `PbPbRawInputPath` in `plot_pbpb_event_sel_{cuts,event_level}.cxx`, `plot_pbpb_fcal_comparison.cxx`.
- Tests: helper unit test ALL OK (edges strict; 2023 window == old |t|<1.5 over a 0.001 ns scan; round trip; legacy key
  removed; stale-file and mismatch throw; real 2023 file accepted). ACLiC compile OK for all 10 macros.
  `PbPbAnalysis` `InitEventSel`: 2023 loads; 2025 (old file) THROWS "only the legacy symmetric ZDC time cut … re-derive".
- Old cuts files backed up: `~/usatlasdata/dimuon_data/pbpb_20{25,26}/pre_cut2_window_20261012/`.
- `event_level` NOT rerun: it derives cut 1 only (no cut-2 dependence) and the nominal cuts file is RECREATEd by `cuts(yr)`
  anyway (the standard chain order gives the same file).
- Docs: `docs/data_analysis.md` cut table; memory `project_pbpb_event_selection`.

## Results & Observations
**Cut-2 preview (2026-10-13 00:00)** — scratch RDF macro `cut2_preview.cxx` (HLT_mu4, good run/LB, cut 1 from the pre-change
file, same t^A–t^C binning as the derivation) → `plots/single_b_analysis/event_selection/pbpb_20{25,26}/cut2_window_preview/`
(`event_sel_cut2_ZDC_time_standalone_*.png`, `event_sel_cut2_ZDC_time_A_1D_*.png`).
| Year | events after cut 1 | new cut 2 keeps | t^A window new/old (t^C in window) |
|---|---|---|---|
| 2025 | 260 593 971 | 88.76 % | 0.9076 |
| 2026 | 247 162 145 | 97.83 % | 0.9865 |
- 2026: the +0.8 ns edge sits in the dip between the main peak and the +1.3 ns cluster: a clean cut.
- 2025: the −0.6 ns edge sits on the rising flank of the main peak (≈2×10⁶ per 0.05 ns vs ≈1.4×10⁷ at the peak). It removes the
  −1.5…−0.6 ns plateau (≈1.3×10⁶/bin) AND part of the main peak's low tail (≈9 % of the old-window events). Open question
  for the user: if the plateau comes from specific runs with shifted ZDC timing, the cut removes good collisions there while
  their lumi stays in.
**Per-run effect** (`cut2_per_run.cxx`; events passing cut 1, t^C inside its window; new t^A window / old |t^A|<1.5):
- 2025: ALL 42 runs 0.887–0.963 (median 0.909). The −1.4 ns shoulder is in every run → NOT a bad-run effect. Since those
  events sit at nominal vertex z (Step 16), the new edge plausibly removes ≈9 % of REAL collisions in every run, i.e. a ≈9 %
  inefficiency that the lumi normalization does not correct. **Flagged to the user.**
- 2026: 522200 0.145, 522041 0.193 (whole collision peak shifted to ≈+1 ns, see data_per_run_signal_yield_qa.md R2);
  522949 0.985; all others ≥ 0.991 (median 0.9995). The +1.3 ns "cluster" is ≈ these two runs: the +0.8 edge effectively
  removes them while their lumi (50.72 + 2.48 µb⁻¹ = 2.2 % of 2026) stays in. The user has meanwhile decided (analysis-3d
  session) to DROP 522200 + 522041 "like D7, before evsel" → after that edit to `PbPbBadRuns(26)`, the 2026 derivation must be
  rerun (the running one includes them). analysis-3d informed 2026-10-13 00:05; it reruns NTuples for pp24/25/26 after the
  event selection is done → ping it on completion.

**New cuts files (verified 02:45)**: 2025 nominal 00:42 / alt 01:40 (this session), 2026 nominal 01:55 / alt 02:42
(analysis-3d, "Excluding 4 bad run(s) and 2 bad-LB range(s)"). `PbPbReadZdcTimeWindow` on all four: 2025 tA (−0.60, 1.50),
2026 tA (−1.50, 0.80), tC (−1.50, 1.50); legacy key gone (also from the UPDATE-mode alt files).
**Acceptance old vs new** (derivation counts, HLT_mu4, centrality 0–79; overlay `plot_pbpb_event_sel_survival_old_vs_new.cxx` →
`plots/single_b_analysis/event_selection/pbpb_20YY/event_sel_survival{,_accum}_old_vs_new_pbpb_20YY.png`; counts saved as
`event_sel_survival_counts_pbpb_20YY.txt` (new) and `pre_cut2_window_20261012/…` (old)):
| Year | C2 per-cut old → new | cumulative after C5 old → new |
|---|---|---|
| 2025 | 97.883 % → **88.820 %** | 97.328 % → **88.296 %** |
| 2026 | 99.190 % → 99.836 % | 98.546 % → 99.192 % (old INCLUDES 522200/522041, new excludes them) |
2025: cuts 3–5 per-cut survival unchanged to ≤ 0.02 %. 2026: the cut-2 improvement is the dropped runs, whose events failed cut 2.
**Per-run kept fraction plots**: `plots/data_qa/cut2_window_per_run/pbpb{25,26}_cut2_window_kept_fraction_per_run*.png` (+ CSVs).

**2025 per-run ZDC t^C vs t^A (2026-10-13 13:07)** — `plotting_codes/event_selection/plot_pbpb_zdc_time_per_run.cxx` →
`plots/single_b_analysis/event_selection/pbpb_2025/zdc_time_per_run/zdc_time_A_vs_C_per_run_pbpb_2025.png` (42 runs, 7×6, unit area,
common log z; red box ±1.5 ns, blue line −0.6 ns) + `zdc_time_regions_per_run_pbpb_2025.csv` (HLT_mu4, good run/LB, after cut 1,
t^C in ±1.5 ns). Shoulder (−1.5 < t^A < −0.6): median 9.4 %, range 3.7–11.6 %; blob (−4 < t^A < −2): median 0.47 %, range
0.00–0.71 %; the largest single-run share of either is ≈5 %. Only 510510 (small run) differs (3.7 %, no blob). → **UNLIKE 2026:
the 2025 off-centre clusters are in EVERY run at the same level, so dropping a few runs cannot remove them.** Options for a
common window across years: keep ±1.5 ns in 2025 (clusters stay in) or keep −0.6 ns and treat the ≈9 % as an efficiency. User
decision.

## Remaining Work
- User: check the plots, finalize cut 2 → then NTuple rerun (25, 26) and downstream.

## Latest Stage
2026-10-13 17:01: **2026 revert DONE end to end.** Event-selection chain: cuts_alt(26) 14:04, cuts(26) 15:10 (±1.5 ns, 4 bad
runs), FCal comparisons 16:20/16:33. pbpb26 NTuples rerun (user): Condor clusters 22 (trig-eff) + 23 (nominal), 7/7 procs exit 0
each. The driver `run_ntuple_rerun_20261013.sh` died at its validation step with a bash syntax error because another session
committed an edit to the script (1d5e044, +pbpb23 case; pbpb26 unchanged) at 16:35 while it ran (bash reads scripts incrementally).
Recovered with VERIFY_ONLY=22,23 (T0 = 16:33:00) from a private copy of the script: all parts validated, hadd OK, run/LB
exclusion check clean (nlb0 0, nbad 0, 31 runs). OS pairs: nominal 829 331 (was 828 345 with +0.8 → +0.12 %, matching cut-2
survival 99.94 % vs 99.84 %), trig-eff 1 537 226. Log `pipelines/logs/ntuple_rerun_20261013b_pbpb26{,_verify}.log`.
STALE downstream for 2026: RDF hist filling, trig-eff fits, crossx/R_AA, data-QA plots (not requested).
Still open (user): 2025 cut 2 — keep −0.6 ns (≈9 % per run, efficiency) or return to ±1.5 ns (clusters in every run).

(previous) 2026-10-13 16:35: 2026 cuts re-derived with ±1.5 ns (nominal 15:10, alt 14:04; "Excluding 4 bad run(s)"; cut 2 keeps
99.94 % of post-cut-1 events vs 99.84 % with +0.8). FCal comparisons finishing. User: rerun the pbpb26 NTuples after the FCal
plots → queued (waits for the FCal "ALL DONE", then `pipelines/run_ntuple_rerun_20261013.sh pbpb26` with CUTS_AFTER = revert
start; log `pipelines/logs/ntuple_rerun_20261013b_pbpb26.log`; nominal + trig-eff, validate, hadd, run/LB exclusion check).

(previous) 2026-10-13 13:50: 2026 revert chain — running cuts(26) and cuts_alt(26) in PARALLEL got cuts(26) OOM-killed (9.3 GB RSS each;
node 23 GB; dmesg "Out of memory: Killed process 837582 (root.exe)"). Lesson: never run two 2026-size cut derivations on this node
at once (25 ∥ 26 overnight worked only marginally). cuts_alt(26) survived; continuation `run_26_revert2.sh`: wait for it →
cuts(26) ALONE → 4 FCal comparisons in parallel.

(previous) 2026-10-13 13:00 — USER DECISION (supersedes the 2026 part of the Objective): with 522200 + 522041 dropped (D8, analysis-3d,
commit 41cd9a1) the +0.8 ns 2026 t^A upper edge is unnecessary → REVERT 2026 to ±1.5 ns and rerun the 2026 event-selection
chain (cuts + cuts_alt; event_level not needed — cut 1 does not depend on cut 2 and was re-derived with D8 at 01:04; then the
4 FCal comparisons). The 2026 NTuples of 03:59 (analysis-3d, +0.8 window) become STALE → flag, not rerun (not requested).
Backup of the +0.8 cuts files/plots: `cut2_a_hi_0p8_20261013/`. Also: per-run ZDC t^C vs t^A for 2025 (one panel per run,
after cut 1) to see whether the off-centre t^A clusters (−1.4 shoulder, −3 blob) come from a few runs (then drop runs and keep
a common ±1.5 window for all years) or from all runs.

(previous) 2026-10-13 02:50: Step 5 cuts DONE for both years (above); overlays + per-run plots delivered. analysis-3d continues with the 4 FCal
comparisons and the pbpb25/26 NTuple reruns. The user is finalizing cut 2 from the plots. Open flag: 2025 loses ≈9 % of
HLT_mu4 0–80 % events at cut 2 (every run 4–11 %); the shoulder events look like real collisions → uncorrected yield loss
unless treated as an efficiency.

(previous) 2026-10-13 00:45: HAND-OVER. The analysis-3d session (tracking doc `ntuple_rerun_pp24_pbpb2526_20261013.md`) killed my 2026
redo (SIGTERM 00:39) and runs `pipelines/run_evsel26_fcal_ntuple_20261013.sh`: 2026 event_level → cuts → cuts_alt (own
build dir, checks "Excluding 4 bad run(s)"), then waits on my 2025 chain PID and checks the 2025 cuts files are fresh, then the
4 FCal comparisons (parallel), then the pbpb25 + pbpb26 NTuple reruns. This session owns ONLY the running 2025 cuts + cuts_alt
(PID 3824921); 2026 / FCal / NTuples are analysis-3d's. Old event-selection plots backed up (`pre_cut2_window_20261012/` in
each 2025/2026 plot dir + FCal comparison dirs). Next here: old-vs-new acceptance overlays (per-cut + cumulative, 25 and 26)
once both years' new cuts exist (old counts from `restore_alt_cuts_20260922.log` / `run_pbpb_all_23_25_26_20260916b.log` /
`logs/pbpb26_badlb_evsel_ntuple_20261012.log`; new from this session's cuts logs (25) and analysis-3d's (26)).

(previous) 2026-10-13 00:20: 2025 cuts + cuts_alt still running (started 23:45). 2026 RESTARTED ~00:15 on the reduced run list: the
analysis-3d session added 522200 + 522041 to `PbPbBadRuns(26)` (user "D8", temporary; UNCOMMITTED edit 00:03 on top of
e17c97e). The first restart still printed "Excluding 2 bad run(s)" because ACLiC reused the in-source `.so` compiled at 23:45
despite the header change. Now force-rebuilt (`.cxx++`, separate build dir): "Excluding 4 bad run(s) and 2 bad-LB range(s)";
the driver aborts otherwise. Driver `run_cut2_26_redo.sh` then waits for 2025 and checks all 4 cuts files are fresh → ping
analysis-3d (my 00:05 heads-up EXPIRED undelivered, held by that session's permission mode) → FCal comparisons (4 parallel).
The 2026 cut-2 preview / per-run numbers were made BEFORE the drop (they include 522200/522041).

(previous) 2026-10-13 ~00:00: step 5 running. cuts + cuts_alt for 25 ∥ 26 started 23:45 (reference durations: cuts(26) 47 min,
cuts_alt(26) 45 min). The top-level driver was stopped so the 4 FCal comparisons could run IN PARALLEL (distinct output
files; serially they took 3 h 44 min on 2026-09-16). A waiter (exits on both year chains' OK/FAIL or on their PIDs
vanishing) starts them. Separately: a fast RDF preview of the cut-2 standalone (+ 1D t^A) for 25/26 →
`plots/single_b_analysis/event_selection/pbpb_20YY/cut2_window_preview/` (user asked for the cut-2 plot first).
Commits: e17c97e (code), 3a0a1bf (INDEX restore + row; e17c97e had accidentally reverted other sessions' INDEX rows
through a truncate-before-read bug in my staging snippet; restored exactly, net diff = 1 added row).
