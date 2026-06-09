# New Skim: Event Selection Pipeline Integration & Full Rerun

## Objective

Integrate the PbPb event selection step into both pipeline scripts
(crossx and trig_eff), fix the PbPb23 file map (add part4), then rerun
both pipelines for all 3 years to update all outputs with the May 2026 skim.

## Context

- New skim recorded in `~/usatlasdata/dimuon_data/data-merging-record.txt`
- All existing event selection outputs (ROOT files with TGraphs, plots) date from May 5 — before the new skim (May 19–27)
- PbPb23 part4 (8.8M entries) was missing from event selection file maps
- Neither pipeline currently includes the event selection step

## Scope

1. Fix PbPb23 file maps in 3 event selection source files (add part4)
2. Add event selection stage to `pipeline_pbpb_crossx.sh` and `pipeline_pbpb_trig_eff.sh`
   - Before NTuple processing
   - Skippable via `SKIP_CONDOR=1` (which also skips event selection)
3. Update `docs/pbpb_pipelines.md`
4. Smoke test the event selection stage
5. Run both pipelines for all 3 years
6. Scan repo for non-pipeline data code producing final results
7. Create analysis status tracking doc

## Implementation Plan

1. [x] Fix PbPb23 file map in `plot_pbpb_event_sel_event_level.cxx` — add part4
2. [x] Fix PbPb23 file map in `plot_pbpb_event_sel_cuts.cxx` — add part4
3. [x] Fix PbPb23 file map in `plot_pbpb_fcal_comparison.cxx` — add part4
4. [x] Add event selection stage to `pipeline_pbpb_crossx.sh` before condor submit
5. [x] Add event selection stage to `pipeline_pbpb_trig_eff.sh` before condor submit
6. [x] Update `docs/pbpb_pipelines.md` with new stage
7. [x] Smoke test: run event selection for one year, verify ROOT output + plots
8. [x] Run both pipelines for all 3 years via /steer-pipeline
9. [x] Scan repo for non-pipeline data code
10. [x] Create analysis status tracking doc

## Design Decisions

- Event selection is controlled by `SKIP_EVSEL` (defaults to `SKIP_CONDOR`).
  Allows skipping event selection independently of Condor submission.
- Created `run_pbpb_all.sh` master wrapper: runs shared event selection once,
  then launches both pipelines in parallel with `SKIP_EVSEL=1`. Prevents the
  OOM issue from running two concurrent 260M-event ROOT processes.
- **OOM incident (Jun 8):** Both pipelines running `plot_pbpb_event_sel_cuts(25)`
  simultaneously caused OOM kill (each ~9 GB RSS). Crossx ROOT process was
  killed by kernel. Trig_eff survived and completed event selection solo.
  Crossx restarted with `SKIP_EVSEL=1`.

## Progress Log

### Step 1-3: Fix PbPb23 file maps
- Added `data_pbpb23_part4.root` to all 3 files, removed stale "needs to be rerun" comments

### Step 4-5: Add event selection stage to both pipelines
- Added "Stage 0: Event selection" before condor submit in both pipelines
- Runs `plot_pbpb_event_sel_event_level(YY)` → produces cuts ROOT files + plots
- Runs `plot_pbpb_event_sel_cuts(YY)` + `plot_pbpb_event_sel_cuts_alt(YY)` → produces all cut plots
- Runs 4 FCal comparison functions (all-year) after per-year loop
- Skipped when `SKIP_CONDOR=1`
- Validates `event_sel_cuts_pbpb_20YY.root` and `_alt.root` exist after each year

### Step 6: Update pipeline docs

### Step 7: Smoke test
- Ran `plot_pbpb_event_sel_event_level(24)` for PbPb 2024 (smallest dataset, 2 parts ~93M entries)
- Completed in ~12 min CPU time
- Both output ROOT files updated: `event_sel_cuts_pbpb_2024.root` (16.7 KB) and `_alt.root` (9.6 KB) — Jun 8 16:16
- 10 event-level diagnostic plots updated to Jun 8 in `event_selection/pbpb_2024/`
- 2 alt banana plots updated to Jun 8 in `event_selection_alternative_banana/pbpb_2024/`
- No errors

## Results & Observations

### Repo Scan: Non-Pipeline Data Code (Step 9)

**Core analysis code NOT in pipelines (produces final/publication results):**

| Code | Location | Reads | Produces |
|------|----------|-------|----------|
| `plot_npairs_vs_centrality.cxx` | `plotting_codes/single_b_analysis/` | PbPb 23/24/25 hadded muon_pairs | Dimuon counts vs centrality (1% bins) |
| `RAA_plotting.cxx` | Analysis root | PbPb 23/24 + pp24 crossx hists | R_AA vs pT, eta, centrality |
| `plot_sig_accept_cutflow_above_60GeV.cxx` | `plotting_codes/single_b_analysis/` | PbPb 23/24/25 hists_cut_acceptance | Signal acceptance cutflow |
| `draw_single_b_statistics_pp_PbPb_23_24_combined.cxx` | `plotting_codes/` | PbPb 23/24 + pp24 | Combined pair pT distributions |
| `trig_effcy_plot_pp.cxx` | `plotting_codes/trig_effcy/` | pp24 RDF output | pp24 trigger efficiency plots |
| `plot_single_b_crossx_pp.cxx` | `plotting_codes/single_b_analysis/` | pp24 RDF crossx output | pp24 crossx plots |

**pp24 pipeline gap:** No production pipeline exists for pp24. Crossx hist
filling and plotting are done via ad-hoc test scripts
(`test_crossx_pp24.sh`, `run_all_crossx.sh`).

**Event selection diagnostics (ad-hoc, not final results):**
8 ZDC/FCal diagnostic codes in `plotting_codes/event_selection/`
(per-run monitoring, preamp fits, etc.) — run manually when investigating
event selection behavior. Not candidates for pipeline inclusion.

**Other diagnostic/investigation macros:** `investigate_H1_fit_vs_data.C`,
`investigate_pass2mu4_vs_both_passmu4.C`, `trig_leg_check.C` — debugging
tools, not final results.

**MC-data comparison:** `plot_mc_data_compr.cxx` and
`plot_mc_data_2D_hists_and_1D_proj.cxx` — infrastructure classes that
read both MC and data; called from other plotting codes.

### Step 8: Pipeline runs (Jun 8–9)

**Event selection (Stage 0) — all 3 years:**
- PbPb23: 124,473,950 events (4 parts incl. part4). Cuts ROOT files + all plots updated.
- PbPb24: 92,650,031 events (2 parts). Cuts ROOT files + all plots updated.
- PbPb25: 260,392,022 events (6 parts). Cuts ROOT files + all plots updated.
- FCal comparison plots (4 functions): all updated.
- OOM incident: crossx ROOT process killed during yr25 cut plots (both pipelines
  running concurrently, ~9 GB each). Trig_eff completed all event selection solo.
  Crossx restarted with SKIP_EVSEL=1.

**Crossx pipeline (P1) — COMPLETE (Jun 8 23:34):**
- Condor clusters 800/801/802 for yr 23/24/25 — all completed.
- hadd: yr23 (4 parts), yr24 (2 parts), yr25 (6 parts) — all validated.
- RDF crossx hist filling: yr23/24/25 — all completed.
- Crossx plots: generated via `plot_single_b_crossx_pbpb.cxx()`.

**Trig eff pipeline (P2+P3) — COMPLETE (Jun 9 00:33):**
- Condor + hadd: all 3 years complete.
- P2 RDF hist filling: yr23/24/25 complete.
- P2 plotting (no-corr trig eff): yr23/24/25 complete.
- P3 RDF hist filling: yr23/24/25 complete.
- P3 plotting (dR corr + plateau norm): all 3 years complete.
- dR correction ROOT files: 6 files (3 years × pt_int + pt_binned).

## Remaining Work

None — all steps complete.

## Latest Stage

Complete. All 10 implementation plan steps done. Both pipelines finished
for all 3 PbPb years. Docs updated. Committed Jun 9.
