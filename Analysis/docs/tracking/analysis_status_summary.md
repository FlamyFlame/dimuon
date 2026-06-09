# Analysis Status Summary

## Objective

Track the current status of all analysis steps across all PbPb years and
data types, documenting what has been rerun with the latest skim and what
remains.

## Data Skim Reference

- **Skim record:** `~/usatlasdata/dimuon_data/data-merging-record.txt`
- **Skim tag:** May 2026 skim (v1/v2)
- **Datasets:**
  - PbPb 2023: 4 parts, 124.5M entries total
  - PbPb 2024: 2 parts, 92.6M entries total
  - PbPb 2025: 6 parts, 260.4M entries total
  - pp 2024: 12 parts (chunked from 2 grid tasks)

## PbPb Analysis Status (All 3 Years: 2023, 2024, 2025)

| Step | Pipeline | Status | Last Run | Notes |
|------|----------|--------|----------|-------|
| Event selection (cuts ROOT files + plots) | Both | **DONE** | 2026-06-08 | All 3 years + FCal comparisons; PbPb23 part4 included |
| NTuple processing (nominal) | P1 crossx | **DONE** | 2026-06-08 | Condor clusters 800/801/802 |
| NTuple processing (trig eff) | P2+3 trig eff | **DONE** | 2026-06-08 | — |
| hadd (nominal) | P1 crossx | **DONE** | 2026-06-08 | yr23: 4 parts, yr24: 2 parts, yr25: 6 parts |
| hadd (trig eff) | P2+3 trig eff | **DONE** | 2026-06-08 | — |
| RDF crossx hist filling | P1 crossx | **DONE** | 2026-06-08 | — |
| Crossx plotting (combined) | P1 crossx | **DONE** | 2026-06-08 | — |
| RDF P2 hist filling (fine q·η) | P2 trig eff | **DONE** | 2026-06-09 | — |
| pT turn-on fitting | P2 trig eff | **DONE** | 2026-06-09 | (run as part of P3 cycle) |
| P2 plotting (no-corr trig eff) | P2 trig eff | **DONE** | 2026-06-09 | All 3 years |
| RDF P3 hist filling (inv weight) | P3 dR corr | **DONE** | 2026-06-09 | — |
| P3 plotting (cross-term dR corr) | P3 dR corr | **DONE** | 2026-06-09 | Single-muon/pair-level/plateau removed (D9: biased on mu4 sample) |

## Non-Pipeline Data Code (Producing Final Results)

These codes read PbPb/pp data and produce final results but are **not**
run as part of either PbPb pipeline. They must be run manually after the
pipelines complete.

| Code | What it produces | Depends on |
|------|-----------------|------------|
| `plotting_codes/single_b_analysis/plot_npairs_vs_centrality.cxx` | Dimuon counts vs centrality | hadded muon_pairs (P1 or P2) |
| `RAA_plotting.cxx` | R_AA vs pT, eta, centrality | PbPb crossx + pp24 crossx |
| `plotting_codes/single_b_analysis/plot_sig_accept_cutflow_above_60GeV.cxx` | Signal acceptance cutflow | hists_cut_acceptance (P1) |
| `plotting_codes/draw_single_b_statistics_pp_PbPb_23_24_combined.cxx` | Combined pair pT distributions | PbPb + pp24 outputs |
| `plotting_codes/trig_effcy/trig_effcy_plot_pp.cxx` | pp24 trigger efficiency plots | pp24 RDF trig eff output |
| `plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx` | pp24 crossx plots | pp24 RDF crossx output |

**Note:** pp24 has no production pipeline. Crossx hist filling uses ad-hoc
test scripts (`RDFBasedHistFilling/test_crossx_pp24.sh`).

## Latest Update

2026-06-09: **All PbPb analysis steps updated with May 2026 skim for all 3 years.**

Both pipelines completed successfully:
- **Crossx pipeline (P1):** completed 2026-06-08 23:34
- **Trig eff pipeline (P2+P3):** completed 2026-06-09 00:33

Event counts (new skim):
- PbPb 2023: 124,473,950 events (4 parts incl. new part4)
- PbPb 2024: 92,650,031 events (2 parts)
- PbPb 2025: 260,392,022 events (6 parts)

Data skim record: [`~/usatlasdata/dimuon_data/data-merging-record.txt`](file:///usatlas/u/yuhanguo/usatlasdata/dimuon_data/data-merging-record.txt)
