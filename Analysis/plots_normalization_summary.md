# Plots Normalization Summary

Summary of which years are plotted, what normalization factors are applied, and whether
distributions are differential, across the three sets of plots.

Trigger assignments are fixed by `MuonObjectsParamsAndHelpers/DatasetTriggerMap.h`:

| (run_year, data_type) | trigger |
|---|---|
| (23, PbPb) | mu4_mu4noL1 |
| (24, PbPb) | mu4 |
| (25, PbPb) | mu4 |
| (24, pp) | mu4_mu4noL1 |
| (24, pp_2mu4) | 2mu4 |

---

## 1. MC/data comparison (`plotting_codes/mc_data_compr/`)

### Years / datasets plotted
- MC: POWHEG (bb + cc combined), Pythia
- Data: **pp 2017** (Run 2), **pp 2024 2mu4**, **pp 2024 mu4mu4noL1**

### Normalization factors
Applied in `helper_functions.c:56` as `h->Scale(norm, "width")`, where `norm` comes from
`PlotMCDataComprBaseClass.h:36` (`norm_factor` array). pp luminosity values are the
single source of truth in `MuonObjectsParamsAndHelpers/PPBaseClass.h`.

| Dataset enum | Dataset | `norm_factor` |
|---|---|---|
| `powheg_bb/cc` | POWHEG bb+cc | `1.` (already in cross-section units) |
| `pythia` | Pythia | `1.` + explicit `Scale(1e6)` at `plot_mc_data_compr.cxx:177` (unit conversion) |
| `pp_2017` | pp 2017 (Run 2) | `1/256.8` pb⁻¹ |
| `pp_2024_2mu4` | pp 2024, 2mu4 | `1/395.836` pb⁻¹ |
| `pp_2024_mu4mu4noL1` | pp 2024, mu4mu4noL1 | `1/113.999` pb⁻¹ |

### Differential?
**Yes** — `hist_helper` always calls `h->Scale(norm, "width")`, dividing by bin width.
When `norm_unity=true`, additionally normalizes to unit integral (also with "width").

---

## 2. RAA plots (`RAA_plotting.cxx` + `SingleBAnalysis/`)

### Years / datasets plotted
Active modes:
- **PbPb 2023 vs pp 2024** (mode 2: pp mu4mu4noL1, mode 3: pp 2mu4)
- **PbPb 2024 vs pp 2024** (mode 4: pp mu4mu4noL1, mode 5: pp 2mu4)

Legend labels are generated dynamically from `DatasetTriggerMap` in `InputOutputPrepare()`.

### Normalization factors
**Hist-filling stage** (`simplified_single_b_analysis_pp.cxx`, `simplified_single_b_analysis_PbPb.cxx`):

*pp* — `h2d->Scale(crossx_factor)` where factor comes from `PPBaseClass::GetCrossxFactor()`:

| (run_year, trigger) | crossx_factor | L_int |
|---|---|---|
| (17, 2mu4) | `1/256.8` pb⁻¹ | 256.8 pb⁻¹ |
| (24, mu4_mu4noL1) | `1/113.999` pb⁻¹ | 113.999 pb⁻¹ |
| (24, 2mu4) | `1/395.836` pb⁻¹ | 395.836 pb⁻¹ |

*PbPb* — per-event `weight × crossx_factor(ctr)` via `CalculateWeightForRAA`
(`PbPbBaseClass.h`), formula per centrality bin:
`crossx_factor = 1 / (Δctr × L_int×1000 × T_AA × L_run/1000)`

| Run | 0–5% formula | Source |
|---|---|---|
| Run 2 (2015+2018) | `1/(0.05 × 7670 × 26.2339 × 1.8172/1000)` | `PbPbBaseClass.h` |
| 2023 | `1/(0.05 × 7800 × 26.1428 × 1.17576/1000)` | `PbPbBaseClass.h` |
| 2024 | `1/(0.05 × 7800 × 26.1428 × 0.85112/1000)` (2023 T_AA placeholder) | `PbPbBaseClass.h` |
| 2025 | `1/(0.05 × 7800 × 26.1428 × 2.59933/1000)` (2023 T_AA placeholder) | `PbPbBaseClass.h` |
| 2026 | `1/(0.05 × 7800 × 26.1428 × 2.35549/1000)` (2023 T_AA placeholder; L after D8, 2026-10-13) | `PbPbBaseClass.h` (`kLumiPbPb2026Nb`) |

**Plot stage** (`RAA_plotting.cxx`): `Scale(1., "width")` applied to both pp and PbPb
projections for modes 1 & 2.

### Differential?
- Modes 1 & 2 (RAA vs pair pT / pair η): **yes** — `Scale(1.,"width")` at plot stage
- Mode 3 (RAA vs centrality): **mixed** — pp integral is width-weighted; PbPb histogram is not width-scaled

---

## 3. Cross-section plots (`plotting_codes/single_b_analysis/` + `RDFBasedHistFilling/`)

### Years / datasets plotted
- pp: run_year argument (default 2024); trigger auto-determined from `DatasetTriggerMap`
- PbPb: run_year argument (default 2024; supports 2023, 2024, 2025, 2026); trigger auto-determined from `DatasetTriggerMap`

### Normalization factors
**Hist-filling stage:**

*pp* (`RDFBasedHistFillingPP::FillHistogramsCrossx()`):
- `pp_crossx_lumi_factor` set in `InitializePPExtra()` via `PPBaseClass::CrossxFactorMap()`
- `crossx_weight = weight × pp_crossx_lumi_factor`; factor is `1/L_int` (same table as RAA above)
- Throws at runtime if called with an unconfigured `(run_year, trigger_mode)` combination

*PbPb* (`RDFBasedHistFillingPbPb::FillHistogramsCrossx()`):
Two histogram sets are filled per centrality bin:
1. **TAA-weighted** (existing names): `weight_for_RAA = weight × crossx_factor(ctr)` — same normalization as RAA hist-filling
2. **Event counts** (`_counts` suffix): raw `weight` (= 1 for data), no luminosity normalization

**Plot stage** (`SingleBCrossxPlotterBase::Save2DColz`, `DrawPairPtByEtaWithDrLines`):
- `Scale(1., "width")` applied to all histograms (1D and 2D) at plot time
- PbPb produces two PNG sets: `<tag>_<var>.png` (event counts) and `<tag>_TAA_weighted_<var>.png` (crossx-normalized)
- Z-axis / Y-axis titles set by caller to reflect differential quantity and units

### Differential?
- All 1D and 2D plots: **yes** — `Scale(1., "width")` applied at plot time in the plotter
- Z-axis title for 2D (pp, TAA-weighted): `d²σ/dp_T d[var] [pb GeV⁻¹]`
- Z-axis title for 2D (PbPb, event counts): `d²N_events/dp_T d[var] [GeV⁻¹]`
- Z-axis title for 2D (PbPb, TAA-weighted): `d²σ/dp_T d[var] [pb GeV⁻¹]`
