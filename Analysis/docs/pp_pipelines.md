# PP Analysis Pipelines

## Overview

Two pipeline scripts drive the pp 2024 dimuon analysis, mirroring the PbPb
pipeline structure:

| Script | Purpose |
|--------|---------|
| `pipeline_pp_crossx.sh` | Crossx/nominal: mass spectra, differential distributions (trigger_mode=3, 2mu4) |
| `pipeline_pp_trig_eff.sh` | Single-muon trigger efficiency (trigger_mode=1, single_mu4, erf+log fitting) |

Unlike PbPb, PP has **no event selection** (no ZDC, no FCal centrality, no
nTrk_HITight, no preamp cuts) and **no centrality binning**.

PP crossx uses `trigger_mode=3` (2mu4 trigger), fundamentally different from
PbPb's `trigger_mode=1` + `mu4_nominal_pbpb_NO_trig_calc=true`. This means
PP requires separate NTuple processing for crossx vs trig eff.


## Environment Variables

| Variable | Default | Pipelines | Description |
|----------|---------|-----------|-------------|
| `POLL_SECONDS` | `45` | Both | Condor polling interval (seconds) |
| `CONDOR_TIMEOUT_SECONDS` | `0` | Both | Max wait for condor (0 = no limit) |
| `SKIP_CONDOR` | `0` | Both | Skip condor submit+wait; reuse existing NTuple outputs |
| `RDF_NTHREADS` | `2` | Trig eff only | ROOT implicit MT thread count |

`DATA_BASE` is hardcoded to `/usatlas/u/yuhanguo/usatlasdata/dimuon_data`.
PP data lives in `${DATA_BASE}/pp_2024/`.


## Pipeline 1: Crossx / Nominal (`pipeline_pp_crossx.sh`)

### Purpose

Produce invariant mass spectra and differential cross-section distributions
for pp 2024. Events are selected using the 2mu4 trigger (`trigger_mode=3`).
The luminosity factor comes from `PPBaseClass::CrossxFactorMap()`:
pp24 2mu4 → 1/410.815 pb⁻¹.

### Stages

| Stage | Description |
|-------|-------------|
| 1 | Submit NTuple nominal condor jobs (`run_pp_24_nominal.sub`, queue 12) |
| 2 | Wait for condor cluster |
| 3 | Validate per-batch ROOT outputs (12 muon_pairs + 12 hists) |
| 4 | hadd into combined files; validate TTree entry counts |
| 5 | RDF crossx hist filling via `run_crossx_hist_filling_pp24.sh` (trigger_mode=3, coarse q·η) |
| 6 | Crossx plotting via `plot_single_b_crossx_pp.cxx()` |

### NTuple Configuration

- `trigger_mode = 3` (2mu4)
- `resonance_cut_mode = 1` (default, no suffix)
- `mindR_trig = 0.02` (runtime detection adds `_mindR_0_02` suffix)

### RDF Configuration

| Setting | Value |
|---------|-------|
| `trigger_mode` | `3` → `base_trig_suffix = "_2mu4"` |
| `useCoarseQEtaBin` | `true` → `_coarse_q_eta_bin` suffix |
| `pp_crossx_lumi_factor` | Set from `PPBaseClass::CrossxFactorMap()` |

### Output Files

| File | Path (relative to DATA_BASE) |
|------|------|
| NTuple per-batch muon pairs | `pp_2024/muon_pairs_pp_2024_partN_2mu4_mindR_0_02.root` |
| Combined muon pairs | `pp_2024/muon_pairs_pp_2024_2mu4_mindR_0_02.root` |
| RDF crossx histograms | `pp_2024/histograms_real_pairs_pp_2024_2mu4_coarse_q_eta_bin.root` |


## Pipeline 2: Trigger Efficiency (`pipeline_pp_trig_eff.sh`)

### Purpose

Derive single-muon trigger efficiency as a function of pT and charge × η
(no centrality binning for PP). The erf+log fitting function is used (PP-
specific; PbPb uses fermi+log).

No Pipeline 3 (dR corrections) — the inverse-weighted procedure is biased
on triggered samples, same as PbPb. The code exists as reference only.

### Stages

| Stage | Description |
|-------|-------------|
| 1 | Submit NTuple trig eff condor jobs (`run_pp_24.sub`, queue 12) |
| 2 | Wait for condor cluster |
| 3 | Validate per-batch ROOT outputs |
| 4 | hadd into combined files |
| 5 | RDF P2 hist filling (trigger_mode=1, fine q·η, inline ROOT macro) |
| 6 | pT turn-on fitting via `single_muon_trig_effcy_pT_fitting()` (erf+log) |
| 7 | Validate fitting output: ≥10 TF1s + ≥4 TH2D fallback histograms |
| 8 | Trig eff plotting via `trig_effcy_plot_pp.cxx()` |

### NTuple Configuration

- `trigger_mode = 1` (single_mu4, default)
- `resonance_cut_mode` forced to `2` when `trigger_effcy_calc=true` → `_res_cut_v2` suffix
- `mindR_trig = 0.02` → `_mindR_0_02` suffix

### RDF Configuration

| Setting | Value |
|---------|-------|
| `trigger_mode` | `1` → `base_trig_suffix = "_single_mu4"` |
| `trigger_effcy_calc` | `true` (auto from trigger_mode=1) |
| `useCoarseQEtaBin` | `false` (default) → fine q·η for fitting |

### Output Files

| File | Path (relative to DATA_BASE) |
|------|------|
| NTuple per-batch muon pairs | `pp_2024/muon_pairs_pp_2024_partN_single_mu4_mindR_0_02_res_cut_v2.root` |
| Combined muon pairs | `pp_2024/muon_pairs_pp_2024_single_mu4_mindR_0_02_res_cut_v2.root` |
| RDF histograms | `pp_2024/histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin.root` |
| pT fit TF1s | `pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root` |


## Key Differences from PbPb

| Aspect | PP | PbPb |
|--------|-----|------|
| Event selection | None | 5-cut sequential (ZDC, FCal, preamp, nTrk) |
| Centrality | None | Centrality-dependent efficiency |
| Crossx trigger | trigger_mode=3 (2mu4) | trigger_mode=1 + NO_trig_calc |
| NTuple split | Separate scripts for crossx vs trig eff | Same trigger, different flags |
| pT fitting | erf+log | fermi+log |
| Years | 24 only | 23, 24, 25 |
| Queue count | 12 (all from one year) | 4+2+6 = 12 (across 3 years) |
| Pipeline 3 | Not included | Commented out (biased on triggered sample) |


## Queue Count

12 condor jobs per pipeline (matching 12 pp24 data files from the May 2026
skim). Queue count is set in both pipeline scripts and `.sub` files.
