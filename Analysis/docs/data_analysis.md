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
