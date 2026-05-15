# Samples and Data Paths

## Data file locations

Base directory: `/usatlas/u/yuhanguo/usatlasdata/`

| Category | Directory | Contents |
|----------|-----------|----------|
| Data (all periods) | `dimuon_data/` | Subdirs: pp_run2, pp_2024, pbpb_run2, pbpb_2023, pbpb_2024, pbpb_2025 |
| Data plots | `dimuon_data/plots/` | Data and data-MC comparison plots |
| Pythia | `pythia/` | Sample files, muon-pair trees, hist root files |
| Pythia plots | `pythia/plots/` | Pythia plots |
| Powheg bb | `powheg_full_sample/bb_full_sample/` | Powheg bb samples, trees, hists |
| Powheg cc | `powheg_full_sample/cc_full_sample/` | Powheg cc samples, trees, hists |
| Powheg plots | `dimuon_data/plots/powheg/` | Powheg plots |

## Trigger modes and file naming

`DatasetTriggerMap.h` maps (year, data_type) to file-name suffix. Must stay in sync with `base_trig_suffix` from `trigger_mode` in `RDFBasedHistFillingData.cxx`:

| `trigger_mode` | `base_trig_suffix` | File pattern |
|----------------|-------------------|--------------|
| 1 | `_single_mu4` | `*_single_mu4_*.root` |
| 2 | `_mu4_mu4noL1` | `*_mu4_mu4noL1_*.root` |
| 3 | `_2mu4` | `*_2mu4_*.root` |

Current assignments:
- `{23,"PbPb"}` -> `"single_mu4"` (trigger_mode=1)
- `{24,"PbPb"}` -> `"single_mu4"` (trigger_mode=1)
- `{25,"PbPb"}` -> `"single_mu4"` (trigger_mode=1)
- `{24,"pp"}` -> `"2mu4"` (trigger_mode=3)

If the map entry doesn't match the actual output filename, the plotter fails to find the file and skips the year silently.

## Skimming run modes

See `SkimCode/README.md` for full details. Key modes:

| Mode | Year | Sample | ZDC |
|------|------|--------|-----|
| `hi2023` | 2023 | data23_hi Pb+Pb | yes |
| `hi2024` | 2024 | data24_hi Pb+Pb | yes |
| `hi2025` | 2025 | data25_hi Pb+Pb | yes |
| `pp2024` | 2024 | data24_5p36TeV pp | no |
| `ppmcfullsim2024` | -- | Pythia8 pp fullsim | no |
| `ppmcfullsim_hioverlay24` | 2024 | overlay (no ZDC in AOD) | no |

## Powheg weighting

Per-event weight: `event_weight = EventWeights[0] * filter_effcy`
- `filter_effcy_bb = 0.003`
- `filter_effcy_cc = 0.001108`
