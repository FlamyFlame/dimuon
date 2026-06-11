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
| RDF crossx hist filling | P1 crossx | **DONE** | 2026-06-10 | w_trig = 1/ε_pair applied; output: `_nominal.root` (was `_coarse_q_eta_bin`) |
| Crossx plotting (combined) | P1 crossx | **DONE** | 2026-06-10 | 66 plots, all 3 years combined |
| RDF P2 hist filling (fine q·η) | P2 trig eff | **DONE** | 2026-06-09 | — |
| pT turn-on fitting | P2 trig eff | **DONE** | 2026-06-09 | (run as part of P3 cycle) |
| P2 plotting (no-corr trig eff) | P2 trig eff | **DONE** | 2026-06-09 | All 3 years |
| RDF P3 hist filling (inv weight) | P3 dR corr | **DONE** | 2026-06-09 | — |
| P3 plotting (cross-term dR corr) | P3 dR corr | **DONE** | 2026-06-09 | Single-muon/pair-level/plateau removed (D9: biased on mu4 sample) |

## PP Analysis Status (pp 2024)

| Step | Pipeline | Status | Last Run | Notes |
|------|----------|--------|----------|-------|
| NTuple processing (nominal, trigger_mode=3) | PP crossx | **DONE** | 2026-06-10 | Condor cluster 12 (12 jobs); `_2mu4_mindR_0_02` suffix |
| NTuple processing (trig eff, trigger_mode=1) | PP trig eff | **DONE** | 2026-06-10 | Condor cluster 10 (12 jobs); `_single_mu4_mindR_0_02_res_cut_v2` suffix |
| hadd (nominal) | PP crossx | **DONE** | 2026-06-10 | 589M combined |
| hadd (trig eff) | PP trig eff | **DONE** | 2026-06-10 | 133M combined |
| RDF crossx hist filling | PP crossx | **DONE** | 2026-06-10 | trigger_mode=3, nominal binning, w_trig = 1/ε_pair correction; output: `_2mu4_nominal.root` |
| RDF P2 hist filling (fine q·η) | PP trig eff | **DONE** | 2026-06-10 | trigger_mode=1, fine q·η for fitting |
| pT turn-on fitting (erf+log) | PP trig eff | **DONE** | 2026-06-10 | 61K fitted TF1s + TH2D fallback |
| P2 plotting (no-corr trig eff) | PP trig eff | **DONE** | 2026-06-10 | 39 plots; fixed `trig_pair`→`trg_pair` typo in TrigEffPlotterPP.cxx |
| Crossx plotting | PP crossx | **DONE** | 2026-06-10 | 5 plots; review PASSED |

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

## Latest Update

2026-06-10: **PP crossx and trigger efficiency pipelines completed with May 2026 skim.**

PP pipelines (both new — first production runs):
- **PP trig eff pipeline:** completed 2026-06-10. 12 Condor jobs (cluster 10), erf+log fitting (61K TF1s), 39 plots. Fixed `trig_pair`→`trg_pair` typo in TrigEffPlotterPP.cxx.
- **PP crossx pipeline:** completed 2026-06-10. 12 Condor jobs (cluster 12), RDF with per-pair 2mu4 trig eff correction (`w_trig = 1/(ε₁·ε₂)`), nominal binning, 5 plots. Recovered from chmod +x fix on `run_pp_24_nominal.sh`.

Key code changes in this session:
- `RDFBasedHistFillingPP.cxx`: added per-pair trigger efficiency correction to crossx weight
- `run_crossx_hist_filling_pp24.sh`: switched from coarse q·η to nominal binning
- `pipeline_pp_crossx.sh`: updated RDF output name to `_2mu4_nominal.root`
- `plot_single_b_crossx_pp.cxx`: input file search uses candidate list (`_nominal.root` → `_coarse_q_eta_bin.root` → bare)

---

2026-06-09: **All PbPb analysis steps updated with May 2026 skim for all 3 years.**

PbPb pipelines completed:
- **Crossx pipeline (P1):** completed 2026-06-08 23:34
- **Trig eff pipeline (P2+P3):** completed 2026-06-09 00:33

Event counts (new skim):
- PbPb 2023: 124,473,950 events (4 parts incl. new part4)
- PbPb 2024: 92,650,031 events (2 parts)
- PbPb 2025: 260,392,022 events (6 parts)

Data skim record: [`~/usatlasdata/dimuon_data/data-merging-record.txt`](file:///usatlas/u/yuhanguo/usatlasdata/dimuon_data/data-merging-record.txt)
