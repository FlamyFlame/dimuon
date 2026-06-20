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

2026-06-19: **PbPb 2023 luminosity corrected (wrong/old GRL) + two b-hadron runs
subtracted from the R_AA luminosity.** The 2023 `mu4` lumi table was re-exported
with the correct GRL (v120-pro33-03): `Prescale Corrected` total 1029.52 →
**1183.650457 µb⁻¹ = 1.18365 nb⁻¹**. The R_AA luminosity subtracts the two
b-hadron runs **461674** (2.623332 µb⁻¹) + **462964** (5.262407 µb⁻¹) — both
already excluded at event level (`PbPbBadRuns(23)` in `PbPbExtras::PassEventSel`),
so numerator and denominator stay consistent — giving **1.17576 nb⁻¹** (was
1.02426). Previously only 462964 was subtracted (461674 was below the old table's
run range). Updated `PbPbBaseClass.h::make_crossx_factors_pbpb_2023` (6 factors) +
`Utilities/PbPbSampledLumi.h` case 23 (coupled; `/review-analysis-code` PASS); 2023
crossx refilled; combined PbPb crossx + R_AA replotted. Net: combined R_AA /
TAA-weighted crossx **× 0.967** (ΣL 4.475 → 4.626 nb⁻¹); per-2023-year crossx
× 0.871. Obsolete 2023 `mu4_mu4noL1` lumi CSV deleted. Docs updated
(analysis_metadata.md, luminosity/README.md, datasets.tex, roadmap). See
`docs/tracking/raa_from_rdf_crossx.md` "REOPENED 2026-06-19 (2)".

---

2026-06-18→19: **PbPb reco-eff placeholder upgraded to the colleague's EXACT Run 2
Medium-μ fits, evaluated directly as TF1.** Replaced the eyeball-digitized F.2
arrays with the colleague's logistic fits
`MuonRecoEffcyRun2MC_medium.root::tf1_eff_fit_cent{C}_eta{E}` (centrality map
{12,13,4,5,6,7,8}; q·η slice i↔eta{i}). The fits are stored in
`run2_reco_eff_placeholder.root` as **TF1** (`tf1_reco_eff_medium_pbpb_...`, 63 of
them) and **evaluated at the exact muon pT** in the lookup — no resampling (the
first iteration sampled into TGraphs; switched to TF1-direct on user feedback that
resampling an analytic fit is needless). `RDFBasedHistFillingData` loads TF1 (PbPb)
+ TGraph (pp). Builder + lookup `/review-analysis-code` PASS (×2). Reran PbPb
crossx RDF (23/24/25), crossx plots, R_AA, stage plots (`/review-plot` PASS on the
resampled version; TF1-direct numerics within ~0.2%). dσ shifts vs pre-reco backup:
ctr0_5 +4.04%, ctr30_50 −5.76%, ctr50_80 −4.44%; reco/raw inflation 1.43–1.81
(central>peripheral, matching the new fits). pp UNCHANGED (colleague file
PbPb-only). Backup `dimuon_data/crossx_hist_backup_20260618_pre_run2_real_fits/`.
Still a placeholder (single-μ ε₁·ε₂ proxy) pending Run 3 3D pair ε_reco. See
`docs/tracking/reco_eff_placeholder_run2.md` Steps 11–12, `docs/placeholder.md` item 3.

---

2026-06-16: **Reco-efficiency PLACEHOLDER applied to nominal crossx (pp + all PbPb years).**
Run 2 single-muon ε_reco proxy (ε₁·ε₂; PbPb dimuon-note F.2, pp HF R_AA Fig.31)
folded into the nominal corrected weight (`*_trig_corr` = base·w_reco·w_trig), so
crossx histograms and the R_AA 3D input are now reco+trig corrected. Crossx RDF +
nominal crossx plots reran. A correction-stage framework (`CorrectionStages.h`)
also saves raw/unfolded/reco/reco+trig variants; before/after 3-line plots in
`plots/sanity_check_crossx/`. Pre-reco backup
`crossx_hist_backup_20260616_pre_reco_nominal/`. **Still preliminary/placeholder**
(replace with proper 3D pair ε_reco when Run 3 MC lands; R_AA *plotting* still
needs task_06 modernization). See `docs/tracking/reco_eff_placeholder_run2.md`,
`docs/placeholder.md` item 3.

---

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
