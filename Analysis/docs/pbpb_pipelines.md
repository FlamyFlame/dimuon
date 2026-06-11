# PbPb Analysis Pipelines

## Overview

Two end-to-end pipeline scripts drive the PbPb dimuon analysis. They share
the same NTuple processing stage (condor batch jobs via `DimuonDataAlgCoreT`)
but differ in how that stage is configured and what follows it:

| Script | Alias | Purpose |
|--------|-------|---------|
| `run_pbpb_all.sh` | Master | Runs shared event selection once, then launches P1 and P2+3 in parallel |
| `pipeline_pbpb_crossx.sh` | Pipeline 1 | Crossx/nominal measurements: mass spectra, differential distributions |
| `pipeline_pbpb_trig_eff.sh` | Pipeline 2 + 3 | Single-muon trigger efficiency (no-correlation) and dR corrections |

Both pipelines process three PbPb years (2023, 2024, 2025) by default.
They use the same `trigger_mode=1` (single mu4) for event selection, but
Pipeline 1 disables trigger efficiency derivation via
`mu4_nominal_pbpb_NO_trig_calc=true` at both the NTuple and RDF stages.

The NTuple processing code (`DimuonDataAlgCoreT`) produces per-batch ROOT
files containing muon pair TTrees and cut-acceptance histograms. The RDF
hist-filling code (`RDFBasedHistFillingPbPb`) reads the hadded muon pair
trees and fills analysis histograms. The fitter
(`SingleMuEffcyPtTurnOnFitter`) fits pT turn-on curves to extract trigger
efficiency parametrizations.


## Environment Variables

All variables are optional. Defaults are shown.

| Variable | Default | Pipelines | Description |
|----------|---------|-----------|-------------|
| `YEARS` | `"23 24 25"` | Both | Space-separated 2-digit years to process |
| `POLL_SECONDS` | `45` | Both | Condor polling interval (seconds) |
| `CONDOR_TIMEOUT_SECONDS` | `0` | Both | Max seconds to wait for condor (0 = no limit) |
| `SKIP_CONDOR` | `0` | Both | Set to `1` to skip event selection and condor submit+wait; reuses existing cuts and NTuple outputs but still validates and hadds |
| `SKIP_EVSEL` | `${SKIP_CONDOR}` | Both | Set to `1` to skip event selection only (still runs condor+RDF+plots). Useful when event selection was already run (e.g. by `run_pbpb_all.sh`) |
| `RDF_NTHREADS` | `2` | Trig eff only | ROOT implicit MT thread count. Value of 8 needs ~48 GB RAM due to per-thread histogram duplication; 2 is safe on standard condor nodes |

`DATA_BASE` is hardcoded to `/usatlas/u/yuhanguo/usatlasdata/dimuon_data`.
Per-year data lives in `${DATA_BASE}/pbpb_20YY/`.


## Master Wrapper: `run_pbpb_all.sh`

Runs the shared event selection once, then launches both pipelines in
parallel with `SKIP_EVSEL=1`. This avoids the OOM risk of two concurrent
260M-event ROOT processes writing to the same output files.

```bash
./run_pbpb_all.sh                     # full run: event sel + both pipelines
SKIP_EVSEL=1 ./run_pbpb_all.sh        # skip event sel, run both pipelines
SKIP_CONDOR=1 ./run_pbpb_all.sh       # skip event sel + condor, reuse existing outputs
```

Sub-pipeline logs are written to `pipelines/crossx_$$.log` and
`pipelines/trigeff_$$.log`. The wrapper waits for both to complete and
reports success/failure.


## Pipeline 1: Crossx / Nominal (`pipeline_pbpb_crossx.sh`)

### Purpose

Produce invariant mass spectra and differential cross-section distributions
for the PbPb dimuon analysis. Events are selected using the mu4 trigger
(`trigger_mode=1`). The crossx histograms are corrected for the per-pair
no-correlation trigger efficiency: `w_trig = 1/(ε₁ + ε₂ − ε₁·ε₂)`,
where εᵢ = ε^{nc}(pTᵢ, q·ηᵢ) from Pipeline 2 TF1 fits. The P2 fit
file must exist before running Pipeline 1's RDF step. Raw (uncorrected)
histograms are also saved with `_no_trig_corr` suffix for comparison.

### Stages

#### Stage 0 -- Event selection (derive cuts + plots)
- Runs `plot_pbpb_event_sel_event_level(YY)` per year to derive
  ZDC/FCal/preamp/nTrk cuts from the raw data and write them to
  `event_sel_cuts_pbpb_20YY.root` and `_alt.root`
- Runs `plot_pbpb_event_sel_cuts(YY)` and `plot_pbpb_event_sel_cuts_alt(YY)`
  per year to produce all cut visualization plots (nominal + alternative banana)
- Runs 4 FCal comparison functions (nominal/alt x all-years/2524) to produce
  FCal scaling comparison plots
- Validates that both `event_sel_cuts_pbpb_20YY.root` and `_alt.root` exist
  after each year
- **Skipped when `SKIP_EVSEL=1`** (defaults to `SKIP_CONDOR`) — the event
  selection only needs to be rerun when the input data changes
- The NTuple processing code (`PbPbEventSelConfig.h`) reads TGraphs from
  these ROOT files to apply event-level cuts

#### Stage 1 -- Submit NTuple condor jobs
- Submits `run_pbpb_YY_nominal.sub` for each year
- Each condor job runs `run_pbpb_YY_nominal.sh`, which sets
  `pbpb_run3_mu4_force_nominal = true` on the `PbPbAnalysis` class
- This flag skips resonance-veto and trigger-efficiency branches in the
  NTuple processing, producing output files **without** the `_res_cut_v2`
  suffix

#### Stage 2 -- Wait for condor clusters
- Polls `condor_q` every `POLL_SECONDS` seconds
- Fails immediately if any job enters held state

#### Stage 3 -- Validate per-batch outputs
- Checks every per-batch ROOT file (muon_pairs + hists_cut_acceptance)
  is non-empty and openable

#### Stage 4 -- hadd per year
- Merges per-batch outputs into combined files per year
- Validates the combined muon_pairs file has non-empty
  `muon_pair_tree_sign1` and `muon_pair_tree_sign2` TTrees

#### Stage 5 -- RDF crossx hist filling + validate
- Runs `run_crossx_hist_filling_pbpbYY.sh` per year
- Each script calls `RDFBasedHistFillingPbPb` with:
  - `trigger_mode = 1`
  - `mu4_nominal_pbpb_NO_trig_calc = true`
  - `mindR_trig = 0.02`
- Validates the RDF output file

#### Stage 6 -- Crossx plotting (combined)
- Runs `plot_single_b_crossx_pbpb.cxx()` which auto-discovers all
  available year files, sums histograms, and produces combined PbPb plots
- PbPb crossx plots are always combined across all years; per-year plots
  are not produced

### RDF Configuration (Pipeline 1)

| Setting | Value | Effect |
|---------|-------|--------|
| `trigger_mode` | `1` | Uses single mu4 trigger; sets `base_trig_suffix = "_single_mu4"` |
| `mu4_nominal_pbpb_NO_trig_calc` | `true` | Sets `trigger_effcy_calc = false`, appends `_no_trg_plots` to `trig_suffix` |
| `hist_filling_cycle` | `1` (generic, default) | Routes to `FillHistogramsCrossx()` |
| `mindR_trig` | `0.02` | Selects input files with `_mindR_0_02` suffix |
| `useCoarseQEtaBin` | `false` (default) | Uses fine q*eta binning; appends `_fine_q_eta_bin` to output |

Derived `trig_suffix`: `_single_mu4_no_trg_plots`
Derived `out_file_suffix`: `_single_mu4_no_trg_plots_fine_q_eta_bin`

### Output Files

Per year (YY = 23, 24, 25):

| File | Path |
|------|------|
| NTuple per-batch muon pairs | `pbpb_20YY/muon_pairs_pbpb_20YY_partN_single_mu4_mindR_0_02.root` |
| NTuple per-batch cut histograms | `pbpb_20YY/hists_cut_acceptance_pbpb_20YY_partN_single_mu4_mindR_0_02.root` |
| Combined muon pairs (hadd) | `pbpb_20YY/muon_pairs_pbpb_20YY_single_mu4_mindR_0_02.root` |
| Combined cut histograms (hadd) | `pbpb_20YY/hists_cut_acceptance_pbpb_20YY_single_mu4_mindR_0_02.root` |
| RDF crossx histograms | `pbpb_20YY/histograms_real_pairs_pbpb_20YY_single_mu4_no_trg_plots_fine_q_eta_bin.root` |

All paths relative to `DATA_BASE`.


## Pipeline 2+3: Trigger Efficiency (`pipeline_pbpb_trig_eff.sh`)

### Purpose

Derive single-muon trigger efficiency as a function of pT, charge * eta,
and centrality, then measure dR-dependent corrections via inverse
weighting. Pipeline 2 produces the no-correlation efficiency; Pipeline 3
applies it as inverse weights to measure the cross-term dR correction
(P(both muons fire) / (ε₁ε₂)). Note: on a mu4-triggered sample, the
cross-term ratio is biased by the event selection (see D9 in tracking
doc); it is kept as a reference for evaluation on unbiased MC samples.

### Pipeline 2: No-Correlation Single-Muon Efficiency

#### Stage 0 -- Event selection (derive cuts + plots)
- Identical to Pipeline 1 Stage 0 (same code, same outputs)
- Both pipelines share the same event selection cuts; running either
  pipeline updates the cuts for all years
- Skipped when `SKIP_CONDOR=1`

#### Stage 1 -- Submit NTuple condor jobs (trig eff mode)
- Submits `run_pbpb_YY.sub` for each year (note: no `_nominal` suffix)
- Each condor job runs `run_pbpb_YY.sh`, which does **not** set
  `pbpb_run3_mu4_force_nominal`, so trig eff derivation is active
- Output files include the `_res_cut_v2` suffix (photoproduction resonance
  veto applied)

#### Stage 2 -- Wait for condor clusters
- Same polling logic as Pipeline 1

#### Stage 3 -- Validate per-batch outputs
- Validates every per-batch ROOT file (muon_pairs + hists_cut_acceptance)

#### Stage 4 -- hadd per year
- Merges per-batch outputs; validates combined TTree entry counts

#### Stage 5 -- RDF Pipeline 2 hist filling (fine q*eta)
- Runs inline ROOT macro per year with:
  - `ROOT::EnableImplicitMT(RDF_NTHREADS)`
  - `trigger_mode = 1`
  - `mindR_trig = 0.02`
  - `hist_filling_cycle` left at default (`1` = generic)
- With `trigger_mode=1` and no `mu4_nominal_pbpb_NO_trig_calc`,
  `trigger_effcy_calc = true`, so `FillHistograms()` routes to
  `FillHistogramsSingleMuonEffcy()`
- Post-processing (`HistPostProcessDataCommon`) projects and divides
  trigger efficiency histograms: `SumSingleMuonTrigEffHists()`,
  `MakeAndWriteSingleMuonTrigEffPtGraphs()`,
  `CalculateSingleMuonTrigEffcyRatios()`
- Fine q*eta binning (`useCoarseQEtaBin=false`) is required for the pT
  fitting stage and for TH2D fallback efficiency lookup

#### Stage 6 -- pT turn-on fitting (Fermi+log)
- Runs `SingleMuEffcyPtTurnOnFitter.cxx` **interpreted** (`.L` without `+`)
  because the CRTP class hierarchy causes rootcling compilation failures
- Calls `single_muon_trig_effcy_pT_fitting_PbPb(YY)` per year
- Uses `fermi_plus_log` fitting mode (PbPb-specific; PP uses `erf_plus_log`)
- Input: the RDF output from Stage 5
- Output: fitted TF1 objects per (q*eta bin, centrality, trigger) stored in
  `trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root`

#### Stage 7 -- Validate fitting output
- Checks the fit file has at least 10 TF1 objects
- Checks the RDF histogram file has at least 4 `TH2D` `_divided` histograms
  (`pt2nd_vs_q_eta2nd`) for TH2D fallback efficiency lookup

#### Stage 8 -- Pipeline 2 plotting
- Runs `trig_effcy_plot_PbPb.cxx(YY)` per year
- Produces 1D and 2D trigger efficiency ratio plots per centrality bin

### Pipeline 3: dR Corrections via Inverse Weighting

#### Stage 9 -- RDF Pipeline 3 hist filling (inv_weight_by_single_mu_effcy)
- Runs inline ROOT macro per year with:
  - `trigger_mode = 1`
  - `mindR_trig = 0.02`
  - `hist_filling_cycle = 2` (inv_weight_by_single_mu_effcy)
- `FillHistograms()` routes to
  `FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()`
- Post-processing produces dR correction graphs via
  `MakeAndWriteDRTrigEffGraphs()`
- The output file is opened in **UPDATE** mode (not RECREATE), appending
  Pipeline 3 histograms into the same file that Pipeline 2 created in
  Stage 5. This means Pipeline 2 must complete before Pipeline 3 runs
  for each year.

#### Stage 10 -- Pipeline 3 plotting (dR corrections)
- Runs `plot_dR_trig_corr.C` (cross-term dR correction overlays, all
  years combined)

**Note (D9):** Only the cross-term dR correction is produced. On a
mu4-triggered sample, all inverse-weighted dR corrections are biased by
the event selection (denominator depleted by P(at least one muon fires)).
The cross-term is kept as a reference for unbiased MC samples. The
single-muon and pair-level terms, and the plateau normalization script,
have been removed.

### RDF Configuration (Pipeline 2)

| Setting | Value | Effect |
|---------|-------|--------|
| `trigger_mode` | `1` | Single mu4; defines `trigs` list and `trigs_pair` combinations |
| `mu4_nominal_pbpb_NO_trig_calc` | `false` (default) | `trigger_effcy_calc = true`; full trig eff derivation active |
| `hist_filling_cycle` | `1` (generic) | Routes to `FillHistogramsSingleMuonEffcy()` |
| `useCoarseQEtaBin` | `false` (default) | Fine q*eta binning for fitting + TH2D fallback |
| `mindR_trig` | `0.02` | Selects `_mindR_0_02_res_cut_v2` input files |

Derived `trig_suffix`: `_single_mu4`
Derived `out_file_suffix`: `_single_mu4_fine_q_eta_bin`

### RDF Configuration (Pipeline 3)

| Setting | Value | Effect |
|---------|-------|--------|
| `trigger_mode` | `1` | Same as Pipeline 2 |
| `hist_filling_cycle` | `2` (inv_weight_by_single_mu_effcy) | Routes to inverse-weighted dR histograms |
| Output file mode | UPDATE | Appends to Pipeline 2's output file |

### Output Files

Per year (YY = 23, 24, 25):

| File | Path |
|------|------|
| NTuple per-batch muon pairs | `pbpb_20YY/muon_pairs_pbpb_20YY_partN_single_mu4_mindR_0_02_res_cut_v2.root` |
| NTuple per-batch cut histograms | `pbpb_20YY/hists_cut_acceptance_pbpb_20YY_partN_single_mu4_mindR_0_02_res_cut_v2.root` |
| Combined muon pairs (hadd) | `pbpb_20YY/muon_pairs_pbpb_20YY_single_mu4_mindR_0_02_res_cut_v2.root` |
| Combined cut histograms (hadd) | `pbpb_20YY/hists_cut_acceptance_pbpb_20YY_single_mu4_mindR_0_02_res_cut_v2.root` |
| RDF histograms (P2 + P3 combined) | `pbpb_20YY/histograms_real_pairs_pbpb_20YY_single_mu4_fine_q_eta_bin.root` |
| pT fit TF1s | `pbpb_20YY/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root` |

All paths relative to `DATA_BASE`.


## File Naming Conventions

Output filenames are constructed from components:

```
histograms_real_pairs_pbpb_20YY_{trig_suffix}_{qEtaBin_suffix}.root
```

Where:
- `trig_suffix` is built by `TriggerModeSettings()`:
  - Base: `_single_mu4` (trigger_mode=1)
  - If `trigger_effcy_calc=false`: append `_no_trg_plots`
  - If `filter_out_photo_resn_for_trig_effcy=false`: append `_no_photo_resn_cuts`
- `qEtaBin_suffix`:
  - `_fine_q_eta_bin` when `useCoarseQEtaBin=false` (default)
  - `_coarse_q_eta_bin` when `useCoarseQEtaBin=true`

NTuple input file search order depends on the pipeline:
- **Nominal** (`mu4_nominal_pbpb_NO_trig_calc=true`): prefers files
  **without** `_res_cut_v2`, falls back to `_res_cut_v2`
- **Trig eff** (default): prefers files **with** `_res_cut_v2`, falls back
  to without


## Data Flow

```
Pipeline 1 (crossx):
  Event sel (cuts + plots) --> NTuple (nominal) --> hadd --> RDF (crossx hists)
                               --> combined crossx plots

Pipeline 2+3 (trig eff):
  Event sel (cuts + plots) --> NTuple (trig eff) --> hadd
                               --> RDF (P2: single-mu eff, RECREATE)
                               --> pT fitting --> validate --> P2 plots
                               --> RDF (P3: inv-weight dR, UPDATE same file)
                               --> P3 plots (cross-term dR corr)
```

Pipeline 2 must complete fully (including fitting and validation) before
Pipeline 3 begins for each year, because:
1. Pipeline 3 opens the RDF output file in UPDATE mode (appends to it)
2. Pipeline 3 reads the pT fit TF1s produced by the fitter in Stage 6


## Queue Counts

Condor queue counts per year (must match `.sub` files). The count reflects
the number of input file batches for that year's dataset.

| Year | Queue count | Condor jobs per pipeline |
|------|-------------|------------------------|
| 2023 | 4 | 4 per pipeline (8 total across both) |
| 2024 | 2 | 2 per pipeline (4 total across both) |
| 2025 | 6 | 6 per pipeline (12 total across both) |

Total per pipeline run: 12 condor jobs (4 + 2 + 6).

The queue counts are identical between Pipeline 1 and Pipeline 2 (same
data, different processing mode), and are set in both the pipeline scripts
(`QUEUE_COUNTS` associative array) and the `.sub` files (`queue N`).
