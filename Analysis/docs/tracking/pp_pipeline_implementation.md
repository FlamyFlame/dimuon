# PP Pipeline Implementation

## Objective

Build two end-to-end pp24 pipeline scripts (`pipeline_pp_crossx.sh` and
`pipeline_pp_trig_eff.sh`) mirroring the PbPb pipeline structure, with all
necessary source code modifications, pp17 legacy cleanup, and consistency
fixes.

## Physics Procedure

### Motivation

The pp 2024 data analysis needs the same automation, validation, and
reproducibility as PbPb. Currently pp runs are entirely manual: individual
condor submits, ad-hoc hadd, manual RDF execution — no pipelines exist.

### PP vs PbPb: Key Physics Differences

1. **No event selection stage.** PP has no ZDC (not installed), no
   nTrk_HITight (HI-specific branch), no FCal-based centrality. PbPb's
   Stage 0 (5-cut sequential event selection) is entirely inapplicable.

2. **No centrality.** PP has no centrality bins; PbPb trig eff is binned
   in centrality.

3. **Different crossx trigger.** PbPb crossx uses `trigger_mode=1`
   (single_mu4) with `mu4_nominal_pbpb_NO_trig_calc=true`. PP crossx
   uses `trigger_mode=3` (2mu4), with `pp_crossx_lumi_factor` set from
   `PPBaseClass::CrossxFactorMap`.

4. **Different trig eff procedure.** PbPb measures single-muon mu4
   efficiency (tag-probe-like, one muon fires mu4 → measure efficiency of
   second). PP measures pair-level 2mu4 efficiency. Specifically:
   - PbPb P2: per-muon ε(pT, q·η, centrality) → P3: dR corrections
   - PP P2: same per-muon structure, `FillHistogramsDimuTrigGivenMu4()`

5. **Different pT fitting function.** PbPb uses `fermi_plus_log`; PP uses
   `erf_plus_log`.

6. **No P3 (dR corrections) for either.** PbPb P3 is commented out due
   to bias in the procedure. The same bias applies to PP. PP has
   `FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` implemented
   (§3c: PP 2mu4 variant) but this should NOT be included in the pipeline.

### Top-level equations

**PP Crossx (trigger_mode=3, 2mu4):**
```
dσ/dX = (1/L_int) × N_observed(X)
```
where `L_int` comes from `PPBaseClass::CrossxFactorMap()`: pp24 2mu4 →
410.815 pb⁻¹.

**PP Trig Eff (trigger_mode=1, single_mu4):**
Same inclusion-exclusion as PbPb (see mu4_trig_effcy_implementation.md),
but without centrality binning.

### Step-by-step method

**Pipeline 1 — PP Crossx:**
1. NTuple processing: `trigger_mode=3` (2mu4), all 12 parts
2. hadd into combined muon_pairs + hists files
3. RDF crossx hist filling: `trigger_mode=3`, `FillHistogramsCrossx()`
4. Crossx plotting: `plot_single_b_crossx_pp.cxx()`

**Pipeline 2 — PP Trig Eff:**
1. NTuple processing: `trigger_mode=1` (single_mu4), all 12 parts, with
   res_cut_v2
2. hadd into combined files
3. RDF P2 hist filling: `trigger_mode=1`, `useCoarseQEtaBin=false` (fine),
   `FillHistogramsSingleMuonEffcy()` → `FillHistogramsDimuTrigGivenMu4()`
4. pT turn-on fitting: `erf_plus_log` mode
5. Validate fitting output (TF1s + TH2D fallback)
6. Trig eff plotting: `trig_effcy_plot_pp.cxx()`

### Negative constraints

- **No event selection.** Do not add ZDC/FCal/preamp/nTrk cuts. PP has
  none of these.
- **No P3 inverse-weighted dR corrections.** Biased procedure; do not
  include in pipeline (same as PbPb).
- **No centrality bins.** PP efficiency is not centrality-dependent.
- **PP crossx uses trigger_mode=3, NOT trigger_mode=1.** The
  `mu4_nominal_pbpb_NO_trig_calc` mechanism is PbPb-specific.

## Context

- 12 pp24 data files from May 2026 skim: `data_pp24_part{1..12}.root`
- `.sub` queue=12 and `PPExtras.c` file_batch_max=12 are consistent
- Only parts 1–4 have NTuple output for single_mu4 mode; none for 2mu4/crossx
- Existing RDF outputs on disk are from old skim (pre-May 2026)
- pp17 is legacy; run scripts to be deleted, source code cleaned up

## Scope

### Phase 1: Source Code Changes & Cleanup

1. Delete pp17 legacy run scripts
2. Clean pp17 from source code (if bulky/confusing)
3. Create separate nominal/trig-eff NTuple run scripts for pp24
4. Create RDF runner script for pp24 crossx
5. Fix any pp/PbPb inconsistencies found in audit

### Phase 2: Pipeline Scripts

6. Write `pipeline_pp_crossx.sh`
7. Write `pipeline_pp_trig_eff.sh`
8. Documentation

## Design Decisions

(None yet — to be populated during implementation.)

## Implementation Plan

### Prerequisite: pp17 Legacy Cleanup

1. [ ] Delete `run_pp_17.sh`, `run_pp_17.sub`, `run_pp_17_MB.sh` (if
   exists) — per §Cleanup
2. [ ] Assess pp17 in source code (PPExtras.c, RDFBasedHistFillingPP.cxx,
   TrigEffPlotterPP.cxx, SingleMuEffcyPtTurnOnFitter.cxx): remove if
   confusing/bulky, leave if inert — per §Cleanup
3. [ ] Delete `single_muon_trig_effcy_pT_fitting_pp.c` (standalone
   legacy; now superseded by unified `SingleMuEffcyPtTurnOnFitter.cxx`
   entry point `single_muon_trig_effcy_pT_fitting()`)

### NTuple Processing Setup

4. [ ] Create `run_pp_24_nominal.sh` + `run_pp_24_nominal.sub`:
   `trigger_mode=3` (2mu4) — per §Pipeline 1 step 1.
   Use `/review-analysis-code`.
5. [ ] Verify existing `run_pp_24.sh` / `run_pp_24.sub` are correct for
   trig eff mode (trigger_mode=1, default) — per §Pipeline 2 step 1.

### RDF Runner Scripts

6. [ ] Create `run_crossx_hist_filling_pp24.sh`: `trigger_mode=3`,
   `useCoarseQEtaBin=false` — per §Pipeline 1 step 3.
   Use `/review-analysis-code`.

### Source Code Fixes (from Audit)

7. [ ] Fix `RDFBasedHistFillingPP::SetIOPathsHook()`: add file candidate
   for `_2mu4_mindR_0_02` suffix (crossx NTuple files will have 2mu4
   suffix) — per §Pipeline 1 step 3.
8. [ ] Update `TrigEffPlotterPP::configureDataFiles()`: pp24 path should
   use `_fine_q_eta_bin` suffix (consistent with useCoarseQEtaBin=false)
   — per §Pipeline 2 step 6.
9. [ ] Update `SingleMuEffcyPtTurnOnFitterPP::configureIO()`: already
   expects `_fine_q_eta_bin` input — verify correct.
10. [ ] Fix `DatasetTriggerMap`: no pp24 entry for trig eff mode
    (single_mu4) — the crossx plotter uses `GetTrigger(24, "pp")` →
    "2mu4" which is correct for crossx but not for trig eff context.
    Assess whether an entry is needed.
11. [ ] Fix `plot_single_b_crossx_pp.cxx`: input path looks for
    `_coarse_q_eta_bin` first, but pipeline will produce
    `_no_trg_plots_fine_q_eta_bin` (or `_coarse_q_eta_bin` depending on
    config). Ensure consistency.

### Pipeline Scripts

12. [ ] Write `pipeline_pp_crossx.sh` — per §Pipeline 1.
    Use `/review-pipeline`.
13. [ ] Write `pipeline_pp_trig_eff.sh` — per §Pipeline 2.
    Use `/review-pipeline`.

### Documentation

14. [ ] Write `docs/pp_pipelines.md` — per §Documentation
15. [ ] Update `README.md` with pp pipeline info

## Progress Log

(Append-only.)

### 2026-06-10: Research & Planning

- Completed full audit of pp vs PbPb analysis across all stages
- Identified all gaps, inconsistencies, and needed changes
- Created this tracking doc with implementation plan

## Results & Observations

### pp17 Legacy Assessment

**PPExtras.c (NTuple):** pp17 appears in:
- Line 9: `if (self().run_year != 17 && self().run_year != 24)` — year validation
- Line 14: `std::map<int, int> run_year_to_file_batch_max = {{17, 3}, {24, 12}}`

**Impact:** Minimal. The `{17, 3}` entry in the map is inert since no
run script will call `PPAnalysis(17, ...)`. The year validation check
is a guard — removing 17 would make it cleaner. **Recommendation: remove
pp17 from PPExtras.c** — the map entry and guard clause add unnecessary
clutter and could confuse a reader into thinking pp17 is still active.

**RDFBasedHistFillingPP.cxx (RDF):** pp17 appears in:
- Lines 30–51: `SetIOPathsHook()` — entire `else if (run_year == 17)`
  block with hardcoded paths
- Lines 55–58: `InitializePPExtra()` — `if (run_year == 17)` overrides
  trigs and trigs_pair

**Impact:** Moderate. The IO paths block is ~20 lines of dead code.
The trigger override (`trigs = {"_mu4", "_2mu4"}`) could cause confusion
about which trigger list is the default. **Recommendation: remove pp17
blocks from RDFBasedHistFillingPP.cxx.**

**TrigEffPlotterPP.cxx (Plotting):** pp17 appears in:
- Lines 5, 13–17: `configureDataFiles()` — `case 17:` with hardcoded
  pp17 paths and plot dirs
- Lines 134, 452–453, 528–530, 712–713, 1074–1076: `isRun2pp` flag
  changes trigger list from `{mu4_mu4noL1, 2mu4}` to `{2mu4}` only

**Impact:** Moderate. The `isRun2pp` flag threads through multiple
drawing functions, changing trigger pair lists. Since pp24 also only
uses 2mu4 (not mu4_mu4noL1), the Run2-specific paths are conceptually
similar but still dead code. The `case 17:` in configureDataFiles has
wrong paths anyway. **Recommendation: remove pp17 from
TrigEffPlotterPP.cxx** — it simplifies trigger list logic.

**SingleMuEffcyPtTurnOnFitter.cxx (Fitter):** pp17 does NOT appear.
`SingleMuEffcyPtTurnOnFitterPP` hardcodes pp_2024 paths. No cleanup
needed.

**single_muon_trig_effcy_pT_fitting_pp.c (standalone):** Legacy wrapper
that calls the same `SingleMuEffcyPtTurnOnFitterPP` class. The unified
`SingleMuEffcyPtTurnOnFitter.cxx` has `single_muon_trig_effcy_pT_fitting()`
as the equivalent entry point. **Recommendation: delete this file.**

### .sub File Consistency

`run_pp_24.sub`: `queue 12` — **matches** 12 data files on disk.
`run_pp_24_MB.sub`: `queue 4` — **inconsistent with 12 data files.**
    This is a MB (min bias) script, likely using a different dataset.
    Needs verification whether MB scripts are needed in the pipeline.
    (MB is for trigger_mode=0 which is not part of the pipeline.)

`PPExtras.c`: `file_batch_max = 12` for year 24 — **consistent.**

### PP Trigger Efficiency: RDF & Downstream Status

**P2 (no-correlation single-muon efficiency):**
- `FillHistogramsSingleMuonEffcy()` → `FillHistogramsDimuTrigGivenMu4()`
  in `RDFBasedHistFillingPP.cxx` — **fully implemented.**
- Fills histograms for `{_mu4, _mu4_mu4noL1, _2mu4,
  _2mu4_AND_mu4_mu4noL1}` triggers with `{_sepr}` bias selection.
- Post-processing: `SumSingleMuonTrigEffHistsPP()`,
  `CalculateSingleMuonTrigEffcyRatios()`,
  `MakeAndWriteSingleMuonTrigEffPtGraphs()` — all implemented.
- pT fitting: `SingleMuEffcyPtTurnOnFitterPP` with `erf_plus_log` —
  **implemented,** reads `_fine_q_eta_bin` input.
- Plotting: `trig_effcy_plot_pp.cxx` → `TrigEffPlotterPP` — **implemented.**

**However, no fine_q_eta_bin RDF output exists on disk for pp24.** Only
`_coarse_q_eta_bin` files exist. The pipeline must produce fine binning.

**P3 (inverse-weighted dR corrections):**
- `FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` — **implemented**
  in `RDFBasedHistFillingPP.cxx` (lines 190–237).
- This is a PP-specific variant (§3c: "PP 2mu4 — one term only, no
  tag/probe, no cross-term"), which differs from PbPb's two-term
  tag/probe structure.
- `MakeAndWriteDRTrigEffGraphs()` → `MakeAndWriteDRTrigEffGraphsHelper({})`
  — implemented (empty categories = no centrality).
- **This code should NOT be included in the pipeline** (biased procedure),
  but the code itself is kept for reference, same as PbPb.

**Consistency with PbPb P3:** The PP P3 implementation is NOT the same
procedure as PbPb. PbPb uses two terms (tag + cross-term), PP uses one
term (pair-level 2mu4). This is correct for the respective physics but
confirms both are kept as reference only.

### PP/PbPb Inconsistencies Found (Beyond Physics Differences)

1. **PP crossx NTuple mode is different from PbPb.** PbPb crossx uses
   `trigger_mode=1` + `mu4_nominal_pbpb_NO_trig_calc=true` (shared
   single_mu4 NTuples with a "nominal" flag). PP crossx uses
   `trigger_mode=3` (2mu4 NTuples are a different dataset). This means
   PP needs **separate** NTuple processing for crossx vs trig eff, while
   PbPb uses the same trigger but different flags.
   **→ Need `run_pp_24_nominal.sh` with `trigger_mode=3`.**

2. **PP crossx plotter expects `_coarse_q_eta_bin` or bare suffix.**
   `plot_single_b_crossx_pp.cxx:66-83` looks for
   `histograms_real_pairs_pp_2024_2mu4_coarse_q_eta_bin.root` first,
   then falls back to `_2mu4.root`. But the RDF output suffix is built
   from `out_file_suffix = trig_suffix + qEtaBin_suffix`. With
   `trigger_mode=3` + `useCoarseQEtaBin=true`: suffix = `_2mu4_coarse_q_eta_bin`. With
   `useCoarseQEtaBin=false`: suffix = `_2mu4_fine_q_eta_bin`.
   The plotter doesn't look for `_fine_q_eta_bin`.
   **→ Decision needed: use coarse or fine for crossx? PbPb crossx uses fine.**

3. **PP RDF `SetIOPathsHook` doesn't handle 2mu4 NTuple files with
   mindR suffix.** The input file search only looks for
   `base_trig_suffix + input_mindR_suffix` patterns. With
   `trigger_mode=3`, `base_trig_suffix = "_2mu4"`. But NTuple output
   filenames follow `muon_pairs_pp_2024_2mu4_mindR_0_02.root` pattern
   only if `mindR_trig > 0` is used in NTuple processing.
   **→ Check: does the NTuple code apply mindR for trigger_mode=3?
   If not, the input_mindR_suffix search is fine as empty.**

4. **PP trig eff plotter `TrigEffPlotterPP` hardcodes
   `_fine_q_eta_bin` in configureDataFiles for pp24,** which is correct
   for the pipeline (fine binning), but the file doesn't exist yet.

5. **No `run_crossx_hist_filling_pp24.sh`.** PbPb has
   `run_crossx_hist_filling_pbpbYY.sh` per year; PP has none.

6. **PP `trig_effcy_plot_pp.cxx` hardcodes `run_year=24`.** PbPb's
   `trig_effcy_plot_PbPb.cxx(YY)` takes year as argument. PP doesn't
   need multi-year but the function signature should match for pipeline
   consistency.

## Remaining Work

All items in Implementation Plan are pending.

## Latest Stage

**Research & planning complete.** Report delivered to user. Awaiting
approval before proceeding to implementation.
