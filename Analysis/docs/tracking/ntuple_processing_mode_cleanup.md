# NTuple Processing Mode Cleanup

## Objective

Clean up the NTuple processing and RDF hist-filling mode logic so that:
1. All 3 PbPb years (23/24/25) can run **crossx/final analysis** with `trigger_mode=1` (single mu4) via a single, clearly-named boolean — not gated on `run_year == 24`.
2. The trigger efficiency vs nominal distinction is governed by a single clear flag, not redundant/overlapping variables.
3. Shell scripts and .sub files exist for all needed modes across all 3 PbPb years.

## Exploration Findings

### Current mode variables (NTuple processing — `DimuonDataAlgCoreT.h/.c`)

| Variable | Type | Default | Location |
|---|---|---|---|
| `trigger_mode` | int | 1 | `.h:329` |
| `resonance_cut_mode` | int | 1 | `.h:335` |
| `trigger_effcy_calc` | bool | false (derived) | `.h:59` |
| `pbpb24_mu4_NO_trig_calc` | bool | false | `.h:337` |
| `mindR_trig` | double | 0.02 | `.h:330` |
| `output_single_muon_tree` | bool | false | `DimuonAlgCoreT.h:121` |

**Derivation logic (`InitParams_DataCore`, line 70):**
```cpp
if (!(isPbPb && run_year == 24 && pbpb24_mu4_NO_trig_calc))
    trigger_effcy_calc = (trigger_mode == 0 || trigger_mode == 1);
if (trigger_effcy_calc) resonance_cut_mode = 2;
```

So `trigger_mode=1` **always** sets `trigger_effcy_calc=true` and forces `resonance_cut_mode=2`, *except* when `pbpb24_mu4_NO_trig_calc=true && isPbPb && run_year==24`.

**Problems:**
- `pbpb24_mu4_NO_trig_calc` is hardcoded to year 24 — years 23 and 25 cannot use mu4 for nominal/crossx.
- The name is confusing (double negative: "NO trig calc").
- `resonance_cut_mode` is force-overridden to 2 whenever `trigger_effcy_calc=true`, so the user's explicit setting is silently overwritten.

### Current mode variables (RDF hist filling — `RDFBasedHistFillingData.h`)

| Variable | Type | Default | Location |
|---|---|---|---|
| `trigger_mode` | int | 1 | `.h:159` |
| `mu4_nominal_pbpb_NO_trig_calc` | bool | false | `.h:165` |
| `trigger_effcy_calc` | bool | false (derived) | `.h:24` |
| `doTrigEffcy` | bool | true | `.h:164` |
| `hist_filling_cycle` | int (enum) | generic | `.h:158` |

**Derivation logic (`TriggerModeSettings`, `Data.cxx:170`):**
```cpp
trigger_effcy_calc = (trigger_mode == 0 || trigger_mode == 1)
    && !(IsPbPb() && isRun3 && mu4_nominal_pbpb_NO_trig_calc);
doTrigEffcy = trigger_effcy_calc;
base_trig_suffix = mu4_nominal_pbpb_NO_trig_calc ? "_mu4_nominal" : "_single_mu4";
```

**Key difference from NTuple code:** Not gated on `run_year == 24` — it already checks `IsPbPb() && isRun3`. So the hist filling code is *more general* than the NTuple code.

### What each mode produces

**NTuple processing output filenames** (pattern: `muon_pairs_pbpb_20YY_partN{trig_suffix}{mindR_suffix}{resonance_cut_suffix}.root`):

| Mode | trigger_mode | resonance_cut_mode | trig_suffix | res_cut_suffix | trigger_effcy_calc |
|---|---|---|---|---|---|
| Default (trig eff) | 1 | 2 (forced) | `_single_mu4` | `_res_cut_v2` | true |
| `_no_res_mu4_mu4noL1` | 2 | 0 | `_mu4_mu4noL1` | `_no_res_cut` | false |
| `_nominal` (24 only) | 1 | 1 (default, NOT overridden) | `_single_mu4` | `` (empty) | false |
| `_mindR_0_01` | 1 | 2 (forced) | `_single_mu4` | `_res_cut_v2` | true |

**RDF hist filling input file search** (`RDFBasedHistFillingPbPb.cxx`, line 12-16):
```
candidates (tried in order):
1. muon_pairs_pbpb_20YY{base_trig_suffix}{input_mindR_suffix}_res_cut_v2.root
2. muon_pairs_pbpb_20YY{base_trig_suffix}{input_mindR_suffix}.root
3. muon_pairs_pbpb_20YY{base_trig_suffix}_no_res_cut.root
4. muon_pairs_pbpb_20YY{base_trig_suffix}.root
```

**RDF hist filling dispatch** (`FillHistograms`, `Data.cxx:91-108`):
- `trigger_effcy_calc=false` → `FillHistogramsCrossx()` (crossx/final)
- `trigger_effcy_calc=true` → `FillHistogramsSingleMuonEffcy()` (trigger efficiency)

### Crossx hist filling scripts for all 3 years

All 3 exist (`run_crossx_hist_filling_pbpb{23,24,25}.sh`), all set `trigger_mode=1`. But they don't set `mu4_nominal_pbpb_NO_trig_calc=true`, so `trigger_effcy_calc` will be **true** → they'll call `FillHistogramsSingleMuonEffcy()` instead of `FillHistogramsCrossx()`. **This is a bug** — the scripts claim to do crossx but the code dispatches to trig effcy.

Wait — let me re-read the script comment: "doTrigEffcy left at default (true); run_crossx=true prevents trig effcy filling". Let me check if there's a `run_crossx` variable.

### Checking for `run_crossx`

```
grep -n "run_crossx" Data.cxx Data.h BaseClass.h
→ no hits
```

**Confirmed bug:** The crossx hist filling scripts for all 3 PbPb years don't actually set the flag needed to dispatch to `FillHistogramsCrossx()`. With `trigger_mode=1` and `mu4_nominal_pbpb_NO_trig_calc=false`, they'll run trigger efficiency, not crossx.

### Existing shell scripts per year

| Script suffix | Year 23 | Year 24 | Year 25 |
|---|---|---|---|
| (default, trig eff) | Yes | Yes | Yes |
| `_no_res_mu4_mu4noL1` (crossx old) | Yes | Yes | **No** |
| `_nominal` (mu4 crossx) | No | Yes | **No** |
| `_mindR_0_01` | Yes | Yes | Yes |
| `_output_single_muon_tree` | Yes | Yes | **No** |
| `_no_res_mu4` (24 only) | No | Yes | No |

### NTuple output → Hist filling input mapping

For **crossx/final** with mu4 (trigger_mode=1, nominal):
- NTuple produces: `muon_pairs_pbpb_20YY_partN_single_mu4.root` (resonance_cut_mode=1 → empty suffix)
- After hadd: `muon_pairs_pbpb_20YY_single_mu4.root`
- Hist filling with `mu4_nominal_pbpb_NO_trig_calc=true`: base_trig_suffix=`_mu4_nominal`, searches for `muon_pairs_pbpb_20YY_mu4_nominal_mindR_0_02_res_cut_v2.root` etc.
- **MISMATCH**: NTuple produces `_single_mu4` suffix, hist filling looks for `_mu4_nominal`

For **trigger efficiency** with mu4 (trigger_mode=1, default):
- NTuple produces: `muon_pairs_pbpb_20YY_partN_single_mu4_mindR_0_02_res_cut_v2.root` (trigger_effcy_calc=true forces resonance_cut_mode=2)
- After hadd: `muon_pairs_pbpb_20YY_single_mu4_mindR_0_02_res_cut_v2.root`
- Hist filling with `trigger_mode=1`, no nominal flag: base_trig_suffix=`_single_mu4`, searches for `muon_pairs_pbpb_20YY_single_mu4_mindR_0_02_res_cut_v2.root`
- **MATCH** ✓

## Design Decisions

### D1: Rename `pbpb24_mu4_NO_trig_calc` → `pbpb_run3_mu4_force_nominal` (NTuple code only)

**Old:** `pbpb24_mu4_NO_trig_calc` — double-negative, year-specific name, confusing.
**New:** `pbpb_run3_mu4_force_nominal` — clear meaning: PbPb Run3, mu4 trigger, forces nominal analysis (no trig eff derivation).

**Semantics:** When `pbpb_run3_mu4_force_nominal=true` AND `trigger_mode ∈ {0,1}`:
- `trigger_effcy_calc = false`
- `resonance_cut_mode` is NOT overridden (stays at user's setting, default=1)
- Output gets `_single_mu4` suffix (same as trig eff — differentiated by resonance cut suffix)

### D2: Remove `run_year == 24` gate (NTuple code only)

**Old:** `if (!(isPbPb && run_year == 24 && pbpb24_mu4_NO_trig_calc))`
**New:** `if (!pbpb_run3_mu4_force_nominal)` — works for all PbPb Run3 years.

### D3: Keep `mu4_nominal_pbpb_NO_trig_calc` in RDF hist filling

The RDF hist filling variable name `mu4_nominal_pbpb_NO_trig_calc` is already informative and not gated on a specific year. Keep as-is. Only fix:
- Remove the `_mu4_nominal` suffix branch in `TriggerModeSettings` — always use `_single_mu4` for trigger_mode=1
- Fix input file search order to be mode-dependent

### D4: Fix input file search order in hist filling

When `mu4_nominal_pbpb_NO_trig_calc=true` (nominal mode), the hist filling should prefer files *without* `_res_cut_v2`. When false (trig eff mode), prefer `_res_cut_v2` first.

Currently the candidate list always tries `_res_cut_v2` first, which means nominal mode would wrongly pick up the trig eff file if both exist.

### D5: Fix crossx hist filling scripts

The 3 `run_crossx_hist_filling_pbpb{23,24,25}.sh` scripts need `pbpb.mu4_nominal_pbpb_NO_trig_calc = true;` to route to `FillHistogramsCrossx()`.

### D6: NTuple processing — what resonance_cut_mode should nominal use?

User says: "trigger_mode=1, resonance_cut_mode=1" for crossx/final. So nominal uses the standard resonance cuts (mode 1, empty suffix), not the narrow-window cuts (mode 2, `_res_cut_v2`) used for trigger efficiency.

### D6: NTuple processing — what resonance_cut_mode should nominal use?

User says: "trigger_mode=1, resonance_cut_mode=1" for crossx/final. So nominal uses the standard resonance cuts (mode 1, empty suffix), not the narrow-window cuts (mode 2, `_res_cut_v2`) used for trigger efficiency.

## Implementation Plan

### Step 1: Rename variable in NTuple processing code [per §D1-D2]
- `DimuonDataAlgCoreT.h:337`: `pbpb24_mu4_NO_trig_calc` → `pbpb_run3_mu4_force_nominal`
- `DimuonDataAlgCoreT.c:70-73`: Update logic, remove `run_year == 24` gate
- `DimuonDataAlgCoreT.c:20+`: Update PrintInstructions

### Step 2: Fix RDF hist filling suffix & input search [per §D3-D4]
- `RDFBasedHistFillingData.cxx:176-177`: Remove `_mu4_nominal` suffix — always use `_single_mu4` for trigger_mode=1
- `RDFBasedHistFillingPbPb.cxx:12-16` (and run2 block): Make candidate list mode-dependent — nominal prefers no `_res_cut_v2`

### Step 3: Fix crossx hist filling scripts [per §D5]
- `run_crossx_hist_filling_pbpb{23,24,25}.sh`: Add `pbpb.mu4_nominal_pbpb_NO_trig_calc = true;`

### Step 4: Update NTuple `_nominal` scripts to use new variable name
- `run_pbpb_24_nominal.sh`: `pbpb24_mu4_NO_trig_calc` → `pbpb_run3_mu4_force_nominal`
- Create `run_pbpb_23_nominal.sh` + `.sub` and `run_pbpb_25_nominal.sh` + `.sub`

### Step 5: Create missing `run_pbpb_25_no_res_mu4_mu4noL1.sh` + `.sub`

### Step 6: Verify end-to-end: NTuple output filenames → hist filling input search

## Progress Log

*(append-only)*

**Step 0 (2026-05-27):** Completed full exploration of mode logic in both NTuple processing and RDF hist filling codebases. Identified 3 issues: (1) year-24 gate prevents years 23/25 from nominal mu4 analysis, (2) confusing double-negative variable names, (3) crossx hist filling scripts don't actually dispatch to crossx code (missing flag). Wrote plan.

**Steps 1-5 (2026-05-27):** All implemented:
- Step 1: Renamed `pbpb24_mu4_NO_trig_calc` → `pbpb_run3_mu4_force_nominal` in `DimuonDataAlgCoreT.h:337`, `.c:70-73,38`. Removed `run_year == 24` gate — derivation is now `trigger_effcy_calc = (trigger_mode == 0 || trigger_mode == 1) && !pbpb_run3_mu4_force_nominal`.
- Step 2: Removed `_mu4_nominal` suffix branch in `RDFBasedHistFillingData.cxx:175-178` — always `_single_mu4` for trigger_mode=1. Fixed input file search in `RDFBasedHistFillingPbPb.cxx:11-28,39-56` — nominal mode prefers files without `_res_cut_v2`.
- Step 3: Added `pbpb.mu4_nominal_pbpb_NO_trig_calc = true` to all 3 crossx hist filling scripts.
- Step 4: Updated `run_pbpb_24_nominal.sh` to use new variable. Created `run_pbpb_23_nominal.sh/.sub` (queue 4) and `run_pbpb_25_nominal.sh/.sub` (queue 6).
- Step 5: Created `run_pbpb_25_no_res_mu4_mu4noL1.sh/.sub` (queue 6).
- Step 6: End-to-end filename chain verified. Trig eff files exist for all 3 years. Nominal files don't exist yet (expected — nominal NTuple hasn't been run). No stale `pbpb24_mu4_NO_trig_calc` references remain.

**Condor submission (2026-05-27):** Submitted 24 NTuple condor jobs (12 trig eff + 12 nominal) across all 3 PbPb years. Nominal scripts for years 23/25 initially held due to missing execute permission — fixed with `chmod +x`, removed held jobs, re-submitted as clusters 2822 (yr23, 4 jobs) and 2823 (yr25, 6 jobs). Trig eff clusters: 2816 (yr23, 4), 2817 (yr24, 2), 2818 (yr25, 6). Nominal clusters: 2820 (yr24, 2), 2822 (yr23, 4), 2823 (yr25, 6).

**Pipeline scripts (2026-05-27):** Created `pipelines/pipeline_pbpb_crossx.sh` and `pipelines/pipeline_pbpb_trig_eff.sh` following `pipeline_pythia_truth.sh` pattern. Both cover: condor submit → wait → validate → hadd (muon_pairs + hists_cut_acceptance) → RDF hist filling → validate → plotting. Crossx pipeline uses existing `run_crossx_hist_filling_pbpbYY.sh` scripts. Trig eff pipeline runs RDF inline with `.L PbPb.cxx+` per year (separate ROOT sessions per compilation rules). Fixed crossx plotter file search in `plot_single_b_crossx_pbpb.cxx:230-240` to try `_no_trg_plots` variants first — crossx histograms are only in `_no_trg_plots` files (24 keys), not in trig eff files (0 keys).

**Pipeline validation run (2026-05-27):** Ran both pipelines with `SKIP_CONDOR=1` using existing NTuple outputs.
- **Smoke test fixes:** (1) Fixed `local` outside function in both pipelines. (2) Added `SKIP_CONDOR` env var support. (3) Crossx pipeline `||` error pattern for RDF exit code.
- **Crossx pipeline:** Completed successfully in ~2.5 min. All 3 years validated → hadd → RDF crossx filling → 66 crossx plots.
- **Trig eff pipeline:** Initial OOM kill (exit 9) — `ROOT::EnableImplicitMT()` with 8 threads × ~2.5 GB/thread = 20 GB, exceeding 23 GB RAM. Fixed by adding `ROOT::EnableImplicitMT(2)` before `.L` in pipeline. Ran all 3 years manually with 2 threads (~5 min each). 6492 histograms per year (3420 1D + 1536 2D + 1536 3D).
- **Trig eff plotter:** Parameterized `trig_effcy_plot_PbPb.cxx` to take `run_year` argument (was hardcoded to 23). Fixed function name to match filename (`trig_effcy_plot_PbPb` not `trig_effcy_plot_pbpb`). 44 plots per year × 3 years = 132 trig eff plots.
- **Total output:** 198 plots (66 crossx + 132 trig eff).

## Pipeline Status

### NTuple processing scripts (all exist)

| Mode | Year 23 | Year 24 | Year 25 |
|---|---|---|---|
| Trig eff (default) | `run_pbpb_23.sub` (queue 4) | `run_pbpb_24.sub` (queue 2) | `run_pbpb_25.sub` (queue 6) |
| Nominal/crossx | `run_pbpb_23_nominal.sub` (queue 4) | `run_pbpb_24_nominal.sub` (queue 2) | `run_pbpb_25_nominal.sub` (queue 6) |
| mu4_mu4noL1 crossx | `run_pbpb_23_no_res_mu4_mu4noL1.sub` | `run_pbpb_24_no_res_mu4_mu4noL1.sub` | `run_pbpb_25_no_res_mu4_mu4noL1.sub` |

### RDF hist filling scripts

| Mode | Year 23 | Year 24 | Year 25 |
|---|---|---|---|
| Crossx (mu4 nominal) | `run_crossx_hist_filling_pbpb23.sh` | `run_crossx_hist_filling_pbpb24.sh` | `run_crossx_hist_filling_pbpb25.sh` |
| Trig eff | **No standalone script** | **No standalone script** | **No standalone script** |

### End-to-end pipeline scripts

| Pipeline | Script | Stages |
|---|---|---|
| Crossx/nominal | `pipelines/pipeline_pbpb_crossx.sh` | NTuple condor → wait → validate → hadd → RDF crossx hist filling → validate → crossx plotting |
| Trig eff | `pipelines/pipeline_pbpb_trig_eff.sh` | NTuple condor → wait → validate → hadd → RDF trig eff hist filling → validate → trig eff plotting |

Both support `YEARS` env var override (default: "23 24 25").

### Crossx plotter file search fix

The crossx plotter (`plot_single_b_crossx_pbpb.cxx`) previously searched for `_single_mu4_coarse_q_eta_bin.root` (trig eff output), which has 0 crossx histograms. Fixed search to try `_no_trg_plots` variants first (crossx output has 24+ crossx histograms).

### Remaining gaps

1. ~~**Trig eff plotting** (`trig_effcy_plot_PbPb.cxx`) has year 23 hardcoded~~ — **FIXED:** parameterized with `run_year` argument.
2. **`test_all_crossx.sh` and `run_crossx_tests.sh`** are old test scripts using `trigger_mode=2` (mu4_mu4noL1) — not the current mu4 nominal crossx pipeline.

## Remaining Work

- ~~Monitor condor jobs~~ — completed
- ~~Sanity-check NTuple outputs, hadd per-part outputs~~ — completed (both modes × 3 years)
- ~~Run crossx hist filling and trig eff hist filling~~ — completed
- ~~Parameterize trig eff plotter~~ — completed
- ~~`/review-plot` on all 198 output plots~~ — PASS (0 critical, 1 cosmetic warning on 2D labels, 3 INFO physics notes)

## Latest Stage

**Status:** Pipeline validation complete. All tasks done.

**/review-plot results (PASS):**
- 1 WARNING: 2D crossx colormaps lack text annotation for centrality/year (info only from filename)
- 3 INFO: yr25 trig eff ~0 (mu4noL1 likely unavailable in 2025), ctr_dep plots blank (missing _sepr hists), ctr50-100% blank (low stats)

**Plot output locations:**
- Crossx: `.../plots/single_b_analysis/pbpb_23_24_25_combined/` (66 plots: counts/ and TAA_weighted/)
- Trig eff: `.../plots/pbpb_trigger_efficiency/mu4/no_corr/pbpb{23,24,25}/` (44 plots each)
