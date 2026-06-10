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

1. **Crossx uses coarse q·η binning; trig eff uses fine.** Both PP and
   PbPb crossx pipelines use `useCoarseQEtaBin=true`. Both trig eff
   pipelines use `useCoarseQEtaBin=false`. Reason: crossx doesn't need
   fine q·η for pT fitting; trig eff does.
   **Old approach:** PbPb crossx was using fine (default). **New:** set
   `useCoarseQEtaBin=true` in PbPb crossx runners + pipeline.

2. **DatasetTriggerMap gets `{24, "pp_single_mu4"}` entry.** Needed for
   trig eff context where the plotter needs to resolve pp24 → single_mu4.

3. **pp17 legacy code removed from all 3 source stages.** PPExtras.c,
   RDFBasedHistFillingPP.cxx, TrigEffPlotterPP.cxx — all pp17 branches
   removed as dead confusing code.

4. **No mindR suffix search for pp 2mu4 NTuple paths.** Documented as
   future non-urgent item. Current NTuple code for trigger_mode=3 does
   not produce mindR suffix files.

## Implementation Plan

### Step 1: pp17 Legacy Cleanup
- Delete `run_pp_17.sh`, `run_pp_17.sub` from NTupleProcessingCode
- Delete `single_muon_trig_effcy_pT_fitting_pp.c` (legacy wrapper)
- Remove pp17 from PPExtras.c: simplify year validation to pp24-only,
  remove `{17, 3}` from file_batch_max map
- Remove pp17 from RDFBasedHistFillingPP.cxx: delete `else if (run_year == 17)`
  block in SetIOPathsHook(), delete `if (run_year == 17)` trigger override
  in InitializePPExtra()
- Remove pp17 from TrigEffPlotterPP.cxx: delete `case 17:` in
  configureDataFiles(), remove `isRun2pp` flag and all references
  (simplify trigger lists to pp24-only)
- Use `/review-analysis-code` with Physics Procedure §3-§5.

### Step 2: NTuple Processing — Create Nominal Scripts
- Create `run_pp_24_nominal.sh`: sets `pp_24.trigger_mode = 3` (2mu4)
- Create `run_pp_24_nominal.sub`: queue 12, points to run_pp_24_nominal.sh
- Verify `run_pp_24.sh` uses default trigger_mode=1 (trig eff mode) — OK as-is
- Use `/review-analysis-code` with Physics Procedure §Pipeline 1 step 1.

### Step 3: PbPb Crossx Binning Fix
- Add `pbpb.useCoarseQEtaBin = true` to all 3 PbPb crossx runner scripts
  (`run_crossx_hist_filling_pbpb{23,24,25}.sh`)
- Add `useCoarseQEtaBin = true` to PbPb crossx pipeline inline RDF calls
  (Stage 5 of `pipeline_pbpb_crossx.sh`)
- Update PbPb crossx plotter search order
  (`plot_single_b_crossx_pbpb.cxx`): prefer `_coarse_q_eta_bin` first
- Use `/review-analysis-code`.

### Step 4: DatasetTriggerMap + PP Crossx Plotter
- Add `{24, "pp_single_mu4"} → "single_mu4"` to both FileSuffixMap and
  LabelMap in DatasetTriggerMap.h
- Verify `plot_single_b_crossx_pp.cxx` search order: already expects
  `_coarse_q_eta_bin` first → correct, no change needed.
- Use `/review-analysis-code`.

### Step 5: RDF Runner Script for PP Crossx
- Create `run_crossx_hist_filling_pp24.sh`: trigger_mode=3,
  useCoarseQEtaBin=true. Pattern after PbPb crossx runners.
- Use `/review-analysis-code`.

### Step 6: Write `pipeline_pp_crossx.sh`
- No event selection stage (PP has none)
- Stage 1: submit NTuple nominal condor jobs (run_pp_24_nominal.sub)
- Stage 2: wait for condor clusters
- Stage 3: validate per-batch outputs
- Stage 4: hadd per year (combined muon_pairs + hists)
- Stage 5: RDF crossx hist filling (trigger_mode=3, coarse q·η)
- Stage 6: crossx plotting
- Use `/review-pipeline`.

### Step 7: Write `pipeline_pp_trig_eff.sh`
- No event selection stage
- Stage 1: submit NTuple trig eff condor jobs (run_pp_24.sub)
- Stage 2: wait for condor clusters
- Stage 3: validate per-batch outputs
- Stage 4: hadd
- Stage 5: RDF P2 hist filling (trigger_mode=1, fine q·η)
- Stage 6: pT turn-on fitting (erf_plus_log)
- Stage 7: validate fitting output (TF1s + TH2D fallback)
- Stage 8: trig eff plotting
- No P3 (dR corrections excluded)
- Use `/review-pipeline`.

### Step 8: Documentation
- Write `docs/pp_pipelines.md`
- Update `README.md` with pp pipeline info

## Progress Log

(Append-only.)

### 2026-06-10: Research & Planning

- Completed full audit of pp vs PbPb analysis across all stages
- Identified all gaps, inconsistencies, and needed changes
- Created this tracking doc with implementation plan

### 2026-06-10: Implementation (Steps 1–8)

**Step 1 (pp17 cleanup):**
- Deleted `run_pp_17.sh`, `run_pp_17.sub`, `run_pp_17_MB.sh`, `run_pp_17_MB.sub`
- Deleted `single_muon_trig_effcy_pT_fitting_pp.c`
- PPExtras.c: simplified year validation to `!= 24`, replaced map with `const int file_batch_max = 12`
- RDFBasedHistFillingPP.cxx: removed `else if (run_year == 17)` block in SetIOPathsHook,
  removed `if (run_year == 17)` trigger override in InitializePPExtra
- TrigEffPlotterPP.cxx: removed `case 17:` in configureDataFiles, removed
  `isRun2pp` flag + all 6 references
- TrigEffPlotterBaseClass.h: removed `isRun2pp` member + its use in buildTriggerList

**Step 2 (nominal scripts):**
- Created `run_pp_24_nominal.sh` (trigger_mode=3) and `run_pp_24_nominal.sub` (queue 12)

**Step 3 (PbPb crossx binning fix):**
- Added `useCoarseQEtaBin = true` to `run_crossx_hist_filling_pbpb{23,24,25}.sh`
- Pipeline Stage 5 calls those scripts directly — no inline RDF to fix
- Fixed `plot_single_b_crossx_pbpb.cxx`: reordered search candidates to prefer
  `_coarse_q_eta_bin` first (was `_fine_q_eta_bin` first)

**Step 4 (DatasetTriggerMap):**
- Added `{24, "pp_single_mu4"} → "single_mu4"` to FileSuffixMap
- Added `{24, "pp_single_mu4"} → "mu4"` to LabelMap
- PP crossx plotter already expects coarse — no change needed

**Step 5 (RDF runner):**
- Created `run_crossx_hist_filling_pp24.sh` (trigger_mode=3, useCoarseQEtaBin=true)

**Step 6 (pipeline_pp_crossx.sh):**
- Created with 6 stages: condor submit → wait → validate → hadd → RDF crossx → plotting
- NTuple suffix: `_2mu4_mindR_0_02` (mindR detected at runtime from skim branches)
- RDF output: `_2mu4_coarse_q_eta_bin`

**Step 7 (pipeline_pp_trig_eff.sh):**
- Created with 8 stages: condor → wait → validate → hadd → RDF P2 → pT fitting → validate → plotting
- NTuple suffix: `_single_mu4_mindR_0_02_res_cut_v2` (resonance_cut_mode forced to 2)
- RDF output: `_single_mu4_fine_q_eta_bin`
- Fitter: `single_muon_trig_effcy_pT_fitting()` (erf+log, no year argument)
- No P3 (dR corrections excluded)

**Step 8 (documentation):**
- Created `docs/pp_pipelines.md`
- Updated `README.md` with pp pipeline entries

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

- Commit all changes (split by logical unit per feedback)

## Latest Stage

**Steps 1–8 implemented.** All source code changes, NTuple scripts, RDF
runners, pipeline scripts, and documentation complete. Ready for smoke
test and review.
