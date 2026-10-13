# Data NTuple rerun pp24 / Pb+Pb 2025 / Pb+Pb 2026 (new GRLs, ZDC cut-2 window, 2026 bad runs) + data-QA refresh

**Created:** 2026-10-13 · **Mode:** Implementation · **Orchestrator session:** analysis-3d

## Objective
Bring the data NTuples of pp24, Pb+Pb 2025 and Pb+Pb 2026 in line with the current run/LB selection, luminosity and
event selection, then refresh the per-run data-QA plots. Inputs that changed:
- pp24: new pp-ref GRL v144-pro42-02, run 488456 LB 84–220 dropped at NTuple stage (`pp24_grl_v144_lumi_lb_exclusion.md`).
- Pb+Pb 2025/2026: per-year A-side ZDC time window (cut 2; `pbpb_zdc_time_cut2_window.md`; 25: t^A > −0.6 ns, 26: t^A < +0.8 ns).
- Pb+Pb 2026: D7 bad runs/LBs (`pbpb2026_analysis_support.md`) + **D8 (this request, TEMPORARY)**: runs 522200 and 522041
  also dropped, before any step incl. event selection.
Downstream (RDF, fits, crossx, R_AA) is NOT in scope: another session propagates it.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) D8 recorded (pbpb2026 doc D8, code `PbPbBadRuns(26)`, L(2026) = 2.35549 nb⁻¹ in `PbPbSampledLumi.h` +
  `PbPbBaseClass.h`, docs/IntNotes); (2) `SkimCode/run_24pp` and `run_26hi` CA.py + xml point at / hold the latest GRLs;
  (3) proofread of code, analysis docs and IntNotes for GRL / lumi / bad runs / cut-2, every stale site fixed;
  (4) after the other agent's 25+26 event-selection chains are final (incl. 2026 re-derivation with D8): Condor NTuple
  reruns DONE + hadd + verified (fresh mtimes, non-empty trees, `lb>0`, bad runs/LBs absent, clean .err) for
  pp24 {2mu4 `_mindR_0_02`, single_mu4 `_res_cut_v2`}, pbpb25 + pbpb26 {`_mindR_0_02`, `_res_cut_v2`};
  (5) data-QA plots for pp24, pbpb25, pbpb26 regenerated with the new GRL/lumi/evsel; (6) reviews closed
  (`/review-analysis-code` for code edits; data-QA = diagnosis plots, exempt), docs + INDEX updated, committed by explicit path.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure
1. **Motivation.** σ = N/(ε·L): N and L must cover the same runs/LBs; event-selection cuts must be the final ones.
2. **Equation.** crossx factor ∝ 1/L_yr; L_yr = Prescale-Corrected lumitable total over exactly the runs/LBs kept at the NTuple stage.
3. **Method.**
   a. Run/LB exclusions (pp: `PPBadLBRanges`; PbPb: `PbPbBadRuns`/`PbPbBadLBRanges`) are applied first, before any cut,
      in the NTuple stage AND in the PbPb event-selection derivation macros.
   b. D8: 2026 runs 522200 (50.7233 µb⁻¹) + 522041 (2.47881 µb⁻¹) join `PbPbBadRuns(26)`;
      L(2026) = 2408.69 − 53.2021 = 2355.49 µb⁻¹ = **2.35549 nb⁻¹** (31 runs).
   c. Event-selection cuts files for 25/26 are re-derived (other agent) before the NTuple rerun.
   d. NTuple variants: nominal (`_mindR_0_02`, crossx) and trigger-efficiency (`_res_cut_v2`) only (user 2026-10-13).
4. **Negative constraints.** No change to any cut value, binning, weight. D8 is not a GRL change: the GRL xml stays the
   ATLAS one; the exclusion is analysis-level code + lumi. Never set a lumi without the matching exclusion. Do not run
   RDF/crossx/R_AA. Do not touch pbpb23/24 NTuples.

## Design Decisions
- DD1. D8 applied exactly like D7 (user 2026-10-13: "Like D7, before evsel"), so the 2026 cuts are re-derived without the runs.
- DD2. Lumi from Total row minus excluded runs (same convention as Step 18; per-run sum differs by 0.009 µb⁻¹, rounding).

## Implementation Plan
1. D8 code + lumi + docs (§3b). /review-analysis-code. — code DONE; docs via step 3; review pending
2. Skim GRL/CA.py check (`run_24pp`, `run_26hi`). — DONE
3. Proofread sweep (2 independent read-only subagents, scratch `_sub_proofread_{A,B}.md` in session scratchpad) → fixes. — DONE except IntNotes .tex (/review-note pending)
4. Wait for the other agent's evsel chains (25, 26 incl. D8 re-derivation). — TODO
5. Condor NTuple reruns + hadd + verification (§3a, §3d). — TODO
6. Data-QA regeneration (pp24, 25, 26). — TODO
7. Close reviews, docs, INDEX, commits. — TODO

## Progress Log
### 2026-10-13 — Step 1 (code) + Step 2 DONE; part13 finding
- D8 code: `PbPbEventSelConfig.h` `PbPbBadRuns(26)` += 522200, 522041; `PbPbSampledLumi.h` 26 → 2.35549;
  `PbPbBaseClass.h` `kLumiPbPb2026Nb` = 2.35549 (+ comments). Unit test (scratch `t_d8.C`): 522200/522041/522384/522408
  bad, 522056/523437 kept, 2023 list unchanged, L26 = 2.35549. D8 recorded in `pbpb2026_analysis_support.md`.
  Other agent (cut-2 chain owner) told to re-derive cuts(26)/cuts_alt(26) after the edit.
- Step 2: `SkimCode/run_26hi` already had the updated xml + make_grid_sub.sh; its `TrigRates_CA.py` synced to
  `scripts/TrigRates_CA.py` (differed only in the pp24 GRL line). `SkimCode/run_24pp`: `TrigRates_CA.py` synced (pp24 GRL
  `physics_ppref2024_all_good.xml`), xml copied in, old copy renamed `physics_2024ppRef_25ns_legend_old_lumi_calib.xml`,
  `grid_sub.sh` + `grid_sub_recovery_part13.sh` `--extFile` → new xml. (run dirs are not git-tracked.)
- FINDING: pp24 part13 (Sep-2026 recovery of run 488534, 652 340 entries, 0 duplicates; documented as "NTuple rerun
  pending") was NEVER NTuple-processed: `PPExtras.c file_batch_max = 12`, pp24 `.sub` `queue 12`, pipelines
  `QUEUE_COUNT=12`. L = 395.836 pb⁻¹ includes run 488534 in full → numerator was short. Fixed: file_batch_max 13,
  `run_pp_24{,_nominal,_no_res_cut,_no_res_cut_trkqcut,_output_single_muon_tree}.sub` queue 13,
  `pipeline_pp_{crossx,trig_eff}.sh` QUEUE_COUNT=13.

### 2026-10-13 00:08 — Step 5 (pp24) launched
Driver `pipelines/run_ntuple_rerun_20261013.sh <pp24|pbpb25|pbpb26>` (backup of combined files to
`<data>/pre_rerun_20261013/`; Condor with `periodic_hold` after 2 h running = watchdog against /pnfs hangs (typical job
10–20 min); held → FATAL; .err must contain only TStreamerInfo warnings; parts fresh+non-empty; hadd; OS pairs > 0; last
line `ALL DONE <ds>` / `FATAL <ds>`). pp24 does not depend on the PbPb evsel chain → started now: clusters 16
(`run_pp_24.sub`, single_mu4 res_cut_v2) + 17 (`run_pp_24_nominal.sub`, 2mu4 nominal), 13 parts each.
Log `pipelines/logs/ntuple_rerun_20261013_pp24.log`. Interpreted build checked under LCG_107a (PPAnalysis(24,13) constructs).

### 2026-10-13 ~00:45 — Step 3 proofread: sweeps A + B cross-checked, fixes applied (except IntNotes .tex)
- Both sweeps: every year-keyed CODE constant correct (lumi 23/24/25/26 = 1.17576/0.85112/2.59933/2.35549, Σ 6.98170;
  pp 395.836; D7+D8; cut-2 window; skim CA.py/xml copies identical to `SkimCode/xmls`). No CRITICAL in code.
- Found + fixed by me: pp bad-LB filter missing in `data_qa/per_run_signal_yield_table.cxx` (pp branch) and
  `hlt_events_per_run.cxx` (default `pp_year=24`, raw `RunNumber`/`lbn`) → run 488456 would have been ~8 % high in the
  per-run N/L QA. Comments: "b-hadron" → ZDC-problem runs (PbPbEventSelConfig.h, PbPbSampledLumi.h, PbPbBaseClass.h);
  pp "no event selection" comments mention the LB exclusion (pipeline_pp_*.sh, pp_norm_factors.cxx). Tracking docs:
  pbpb2026 §1/§3d/§4/P3/Remaining/Latest + D8, INDEX row + scope, pp24_crossx_rerun_2026_08.md:80, roadmap Q2.1,
  ntrk_lowband banner, muon_wp_registry, slides_instruction (cut 2 was ±1.8 ns → per-year window).
- Doc subagent (log `_sub_docfix.md`, merged here): analysis_overview.md:35 (only that line), analysis_metadata.md,
  luminosity README (also its "GRL total" column: 7.03490 was the R_AA sum; GRL sum 7.04279), data_analysis.md
  (bad-run line incl. D7/D8/LB ranges/pp24), placeholder.md, analysis_status_summary.md, plots_normalization_summary.md
  (PbPb lumis were stale 1.3896/1.5411), signal_selection_change_impact.md (+2026), pbpb_pipelines.md, pp_pipelines.md
  (queue 13), Analysis/README.md, KB analysis/overview.md (cut 2 window; cut 3 wording corrected vs PbPbExtras.c:240),
  KB data/variables.md, KB data/samples.md, SkimCode/README.md GRL table, SkimCode/run_26hi/grid_sub.sh regenerated with
  the unmodified make_grid_sub.sh (diff = date, GRL header, --extFile → _updated.xml ×5).
- HAZARD (user, not fixed): `SkimCode/run_26hi/make_grid_sub.sh` run in place now picks up `part6/7_recovery` lists as
  ordinary parts (no --inputFileList → double counting) and overwrites historical grid_sub_part1..5.sh; tag v2 used.
- OPEN: IntNotes tex (datasets.tex has no 2026; introduction.tex:35 "2023--2025"; trigger_efficiency.tex:183 PHtext) →
  /review-note; data_per_run_signal_yield_qa.md Latest Stage (other session's doc) → at data-QA step.

### 2026-10-13 00:10–00:45 — Code review APPROVED (it.3); 2026 evsel re-derivation restarted; pp24 NTuple DONE
- `/review-analysis-code` log `.claude/logs/review-analysis-code-20261013-001426-d8-part13-ppqa-drivers.md`: it.1 FAIL
  (CRITICAL: my first 2026 driver skipped `event_level(26)` → cut 1 (ZDC–FCal band) would keep 522200/522041; + driver
  validation gaps, ScrambGenPP kNParts 12), it.2 FAIL (new driver not executable; no output-level exclusion check),
  it.3 PASS. Commits 7e46aff (D8), 7b2aa5b (part13), 973f537 (data-QA pp filter), a807216 (drivers +
  `pipelines/check_runlb_exclusion.C`).
- Cross-session: my two messages to the cut-2 chain owner ("Abnormal runs detection plots") expired undelivered
  (permission-mode hold). Its first cuts(26) died silently (log stops at "scanning…"). It then started
  `run_cut2_26_redo.sh` (cuts + cuts_alt only, no event_level) concurrently with mine on the same files → USER
  (2026-10-13 ~00:40): keep mine, stop theirs → killed pid 3908957 + its root before any write (files still 10-12).
  Its 2025 chain (pid 3824921, cuts+cuts_alt 25) continues (user); my driver waits for it.
- 2026 driver `run_evsel26_fcal_ntuple_20261013.sh` relaunched 00:21:55 (event_level → cuts → cuts_alt, all logs must
  show "Excluding 4 bad run(s)"; then wait for 2025; 4 FCal comparisons; pbpb25+26 NTuples). MANUAL CHECK owed:
  `event_level_26.log` must contain "Cut TGraphs saved to" and "Alt cut TGraphs saved to" (review it.2 INFO 3).
- pp24 NTuple DONE: clusters 16/17 13/13 ExitCode 0; re-validated with the new driver (VERIFY_ONLY, log
  `pipelines/logs/ntuple_rerun_20261013_pp24_verify.log`): single_mu4 res_cut_v2 OS 780,280; 2mu4 nominal OS 2,652,530
  (was 2,687,070); exclusion check nlb0 = 0, nbad = 0, max LB 1613, 8 runs, both variants. Run 488456 OS 473,508 → 435,818
  (LB 84–220 dropped); run 488534 OS 440,430 → 443,580 (+3,150 from part13). Old combined files in
  `~/usatlasdata/dimuon_data/pp_2024/pre_rerun_20261013/`.
- /review-note (log `.claude/logs/review-note-*-pbpb26-datasets-grl.md`) running: datasets.tex 2026 (events 270,087,106;
  L 2.35549; GRL + D7/D8 with \PHtext for D8), introduction 2023–2026, trigger_efficiency PHtext removed; it.1 FAIL
  (placeholder.md row 6 stale) → fixed; it.2 pending.
- Data QA pp24 started (`scratchpad/qa_pp24.sh`, log `pipelines/logs/qa_pp24_20261013.log`); pre-change plots backed up
  to `~/usatlasdata/dimuon_data/plots/data_qa/pre_rerun_20261013/`.

## Remaining Work
All steps.

## Latest Stage
2026-10-13 ~00:50: pp24 NTuple DONE+verified; pp24 QA running; 2026 evsel driver in stage A (then waits for the other session's
2025 chain, FCal comparisons, pbpb25/26 NTuples); /review-note it.2 running. Owed: manual "saved to" check of event_level_26.log;
pbpb25/26 output checks; QA 25/26; review-analysis-code-20261012-182614 it.3 (Step-17 verification); commits.
(previous) 2026-10-13 ~00:30: steps 1–2 done, step 3 sweeps running, step 4 waiting on session "Abnormal runs detection plots" (cut-2 chain owner).
