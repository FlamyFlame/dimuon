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

## Remaining Work
- User: check the plots, finalize cut 2 → then NTuple rerun (25, 26) and downstream.

## Latest Stage
2026-10-13: steps 1–4 done; step 5 chain running (scratch driver `run_cut2_chain.sh`: cuts + cuts_alt 25 ∥ 26, then 4 FCal comparisons). Next: validate outputs (freshness + keys), before/after table, /review-analysis-code (user: no /review-plot for these single-value cut changes).
