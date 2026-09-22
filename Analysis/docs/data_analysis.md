# Data Analysis (pp and Pb+Pb)

<!-- TODO: Full documentation to be written. -->

## Overview

Processes real collision data from ATLAS Run 2 (pp17, PbPb15, PbPb18) and Run 3 (pp24, PbPb23, PbPb24, PbPb25, PbPb26).

**Active datasets (Run 3):** pp24, PbPb23, PbPb24, PbPb25, PbPb26. Run 2 datasets exist for cross-checks but are no longer maintained for new analysis decisions.

## Key files

| File | Role |
|------|------|
| `NTupleProcessingCode/DataAnalysisClasses.h` | `PPAnalysis`, `PbPbAnalysis` |
| `NTupleProcessingCode/DimuonDataAlgCoreT.{h,c}` | Data algorithm core |
| `NTupleProcessingCode/PPExtras.{h,c}` | pp-specific processing |
| `NTupleProcessingCode/PbPbExtras.{h,c}` | PbPb-specific (centrality, event selection) |
| `NTupleProcessingCode/PbPbEventSelConfig.h` | PbPb event selection cut configuration |
| `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` | pp RDF filler |
| `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` | PbPb RDF filler |
| `plotting_codes/trig_effcy/TrigEffPlotterPP.cxx` | pp trigger efficiency |
| `plotting_codes/trig_effcy/TrigEffPlotterPbPb.cxx` | PbPb trigger efficiency |

## Trigger modes

| Dataset | Trigger | `trigger_mode` | Suffix |
|---------|---------|----------------|--------|
| PbPb23/24/25 | single mu4 | 1 | `_single_mu4` |
| pp24 | 2mu4 | 3 | `_2mu4` |

See `DatasetTriggerMap.h` for the full mapping.

## Resonance-cut modes (`resonance_cut_mode`)

The OS resonance veto is applied **OS-only at the ntuple stage** (`DimuonAlgCoreT::ResonanceTagging`). `resonance_cut_mode` selects the veto list and the output-file suffix:

| Mode | Suffix | Veto list (`ParamsSet`) | Used for |
|------|--------|--------------------------|----------|
| 1 (default) | `` (file ends `_mindR_0_02`) | `minv_cuts` (V1): `{0,1.06},{2.9,3.3},{3.55,3.8},{9.08,10.5}` — removes the **entire** light-resonance continuum below 1.06 | **Nominal / cross-section / R_AA** |
| 2 | `_res_cut_v2` | `minv_cuts_v2`: narrow per-peak `{0,0.6},{0.72,0.85},{0.94,1.06},{2.9,3.3}` | **Trigger efficiency** (`trigger_effcy_calc=true` auto-sets mode 2) |
| 0 | `_no_res_cut` | none | **Resonance templates** (φ/J/ψ shapes), MC studies needing resonances present |

**Why nominal uses V1, trigger-eff uses V2:** for the **generic** (no-signal-cut) same-/opposite-sign histograms, V1 removes the light-resonance continuum *cleanly*, whereas V2 cuts only the narrow peak windows and **leaves continuum tails between the windows** that would skew the generic SS/OS plots. Trigger efficiency *wants* V2 precisely so it can probe efficiency in the low-mass region *between* the peaks.

**Key fact:** after the **signal selection** (minv ∈ [1.08, 2.9]) V1 / V2 / none are **identical** — the signal window already excludes every resonance window. So the choice only matters for non-signal-cut (generic / 0–4 GeV template-fit) histograms. (Convention also in memory `project_resonance_cut_modes`; OS-only veto detail in `project_os_resonance_veto`.)

## Condor submit files

- `run_pp_17.sub`, `run_pp_24.sub` -- pp data
- `run_pbpb_23.sub`, `run_pbpb_24_nominal.sub` -- PbPb data
- Various specialized `.sub` files for trigger studies, single-muon trees, etc.

## PbPb event selection

**Five cuts, applied in this order.** The order and the names are fixed by the
`PbPbEvSelCut` enum in `NTupleProcessingCode/PbPbEventSelConfig.h`, which is included by
BOTH the derivation and the application so the two cannot drift; every threshold is
*derived per year from that year's own data* and stored in
`~/usatlasdata/dimuon_data/pbpb_20YY/event_sel_cuts_pbpb_20YY.root`.

| # | `kPbPbEvSelCutLabel` | What is required | Form of the threshold |
|---|---|---|---|
| 1 | `ZDC_FCal_banana` | ZDC E_total below the banana curve at that FCal E_T^{A+C} | TGraph `g_ZDC_FCal_cut`, per-0.1-TeV-FCal-slice two-band Gaussian fit: `max(main μ+5σ, pile-up μ_bg−3σ_bg)` |
| 2 | `ZDC_time` | \|t_A\| < 1.5 ns AND \|t_C\| < 1.5 ns | scalar `ZDC_time_cut_ns` |
| 3 | `ZDC_preamp` | NOT (both sides' preamp sums above threshold) — fails only if A **and** C exceed it | scalars `ZDC_preamp_{A,C}_cut_ADC`; per-run μ+7σ from the optional `t_preamp_per_run` tree where the year's cuts file provides it (2025, 2026) |
| 4 | `nTrk_frac` | N_trk^HItight / N_trk^total above a lower bound, vs N_trk^total | TGraph `g_ntrk_frac_cut_lo`, per-slice Gaussian μ−5σ |
| 5 | `nTrk_FCal_band` | N_trk^HItight inside a band vs FCal E_T^{A+C} | TGraphs `g_ntrk_fcal_cut_{lo,hi}`, per-slice Gaussian [μ−5σ, μ+5σ] |

Cut 1 targets the out-of-time / pile-up band that sits above the main hadronic band in the
ZDC–FCal plane; cuts 4 and 5 target events whose track multiplicity does not match the
calorimeter activity. The **nominal** cut-1 procedure is the two-band Gaussian fit above; an
**alternative quadratic** ("alt banana") procedure exists as a cross-check only — it writes
`event_sel_cuts_pbpb_20YY_alt.root` and its own plot directory, and **no analysis stage reads
it**. Both are documented slide-by-slide, with the physics of the two bands, in
`docs/tracking/event_selection_banana_cut_comparison.md`.

**Where things live:**
- Definition / key names / cut order: `NTupleProcessingCode/PbPbEventSelConfig.h`
- Application: `PbPbExtras::PassEventSel()` (and `InitEventSel()`, which throws if the year's
  cuts file is missing) in `NTupleProcessingCode/PbPbExtras.c`
- Derivation + per-cut figures: `plotting_codes/event_selection/plot_pbpb_event_sel_event_level.cxx`
  (event-level distributions, cut 1) then `plot_pbpb_event_sel_cuts.cxx` (cuts 2–5)
- Run-quality exclusions: `PbPbBadRuns()` in `PbPbEventSelConfig.h` (2023: 461674, 462964;
  no other year) — any change must move the luminosity in `Utilities/PbPbSampledLumi.h` and
  `PbPbBaseClass.h::make_crossx_factors_pbpb_<yr>()` in the same commit.

<!-- TODO: Document pipeline steps, validation, plotting -->
