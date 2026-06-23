# Data Analysis (pp and Pb+Pb)

<!-- TODO: Full documentation to be written. -->

## Overview

Processes real collision data from ATLAS Run 2 (pp17, PbPb15, PbPb18) and Run 3 (pp24, PbPb23, PbPb24, PbPb25).

**Active datasets (Run 3):** pp24, PbPb23, PbPb24, PbPb25. Run 2 datasets exist for cross-checks but are no longer maintained for new analysis decisions.

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

5-cut sequential selection (see `PbPbEventSelConfig.h`):
1. ZDC timing cut
2. FCal/ZDC correlation
3. ZDC preamp saturation (per-year threshold)
4. Centrality range
5. Vertex / track multiplicity

<!-- TODO: Document pipeline steps, validation, plotting -->
