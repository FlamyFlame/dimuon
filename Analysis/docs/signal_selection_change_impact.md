# Signal-Selection Change — Impact & Rerun Map

**Type:** Repo reference doc (durable). **Scope:** what must be recompiled,
rerun, and re-plotted whenever the **single-b dimuon signal-region selection**
changes — a cut added, removed, or its value changed. Applies to *any* such
change, not one specific edit; it is also the checklist for **systematic
uncertainty** variations of the selection (e.g. ΔR, minv, pair-pT, q·η bounds).

> **First executed run of this map:** the ΔR>0.05 removal (2026-06-22) — see
> `docs/tracking/remove_dr_cut_signal_selection.md` for a worked example of the
> full chain below (recompile 7 classes → rerun pp+pbpb crossx + pythia truth →
> replot crossx/R_AA/acceptance/cutflow). Key lesson logged there: a cut can live
> in MORE than one signal-region definition per file (PythiaTruth had both the
> acceptance `signal_cuts` AND the template-fit `kin_cuts`) — always grep, don't
> trust an enumerated line list.
>
> **Ground truth for the selection itself:** `docs/analysis_overview.md` §2
> (signal region) and the relevant tracking-doc Physics Procedure. This doc
> owns only the *engineering dependency graph*. Pipeline/stage details:
> `README.md` + `pipelines/`.

---

## 0. The current signal region (reference)

Data (reco), identical in pp and PbPb crossx:
```
minv > 1.08 && minv < 2.9 && pair_pt > 8
  && m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2 && dr > 0.05
```
Truth analog (Pythia/Powheg): same with `truth_*` variables + `from_same_b`.

> **Note on q·η:** the cut is one-sided per muon (`q*eta < 2.2`, no explicit
> lower bound; the floor is the muon |η|≈2.5 detector edge, mapped into the
> −2.4…2.2 efficiency binning). A change here has the same blast radius below.

> **Note on the ΔR cut (motivation, for systematics):** `dr > 0.05` was added
> because the **data-based** dR-dependent trigger-efficiency inverse-weighting
> lacked statistics at small ΔR (the `minv > 1.08` cut drives the ΔR
> distribution toward zero at small ΔR). It induces a strong pair-pT–dependent
> signal acceptance (collimation: ΔR shrinks as pair pT rises — see the cutflow
> `…/single_b_analysis/pythia/pythia_sig_accept_above_60GeV*.png`). The dR
> trigger correction is being **moved to MC**; the correct (possibly absent)
> pair-pT–dependent ΔR cut can only be fixed once that MC exists.

---

## 1. What does NOT change (boundaries of the blast radius)

- **NTuple processing / `muon_pairs_*` trees** — UNCHANGED. The selection is
  applied downstream in RDF hist-filling, not at tree creation. Do **not**
  reprocess ntuples or resubmit Condor.
- **Trigger-efficiency derivation (P2 ε^nc fits, P3)** — UNCHANGED. The
  single-muon tag-and-probe ε^nc(pT, q·η) fits and their `_fine_q_eta_bin`
  inputs are independent of the *pair* signal-region cuts. `pipeline_pp_trig_eff.sh`,
  `pipeline_pbpb_trig_eff.sh`, `SingleMuEffcyPtTurnOnFitter` and all trig-eff
  plots stay valid. (The *applied* per-pair trigger weight inside crossx is
  recomputed automatically when crossx re-runs — no separate action.)
- **Reco-eff PLACEHOLDER file** (`run2_reco_eff_placeholder.root`, Run 2 TF1,
  function of pT & q·η only — **no ΔR dependence**) — UNCHANGED. It is what the
  nominal crossx currently folds in as `w_reco`. See §4 for why the fullsim
  reco-eff re-derivation is code-consistency-only today.

---

## 2. [RECOMPILE] — code carrying the signal cut (ACLiC `.L file.cxx+`)

Edit the cut in **all** of these (keep them in sync):

| File | Where | Feeds |
|---|---|---|
| `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` | :382 (crossx); :552 (ΔR-binned "no-minv" variant) | pp data crossx |
| `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` | :920 (crossx); :1195 (ΔR-binned variant) | PbPb data crossx |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaTruth.cxx` | :412 | truth signal acceptance |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegTruth.cxx` | :158 | Powheg truth acceptance |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsim.cxx` | :129/:131 (`pass_signal_truth`/`_reco`) | pp reco-eff / det-resp |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsimOverlay.cxx` | :69/:71 | PbPb (overlay) reco-eff |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsim.cxx` | :185/:186 | Powheg reco-eff (obsolete demo) |
| `plotting_codes/single_b_analysis/plot_sig_accept_cutflow_above_60GeV.cxx` | :40 (`kCuts`) | cutflow diagnostic — **must mirror the cut list/order** |

**Legacy / retired (leave, or fix for hygiene only):**
`SingleBAnalysis/SingleBAnalysisBase.cxx:25`,
`RDFBasedHistFilling/FillHistogramsCrossx_PP_clean.cxx:8` — both still use the
old `pair_eta < 2.2` form; not in the active chain.

> **⚠ ΔR-axis floors (only relevant when changing/removing the ΔR cut).** The
> ΔR-binned 2D histograms are booked with a **uniform axis from 0.05 to 1.0**:
> PP `:457/:519/:534`, PbPb `:1000/:1073/:1089/:1136/:1154`
> (`make_unif_edges(50, 0.05, 1.0)`). If `dr > 0.05` is removed, extend these
> axes to **0.0** or the newly-admitted small-ΔR pairs fall in underflow and the
> 2D ΔR distributions stay clipped. (Not needed for non-ΔR cut changes.)

---

## 3. [RERUN HIST FILLING] — order, scripts, outputs

**A. Data crossx (always required):**
1. `RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh`
   → `dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root`
2. `run_crossx_hist_filling_pbpb23.sh`, `…pbpb24.sh`, `…pbpb25.sh`
   → `dimuon_data/pbpb_20YY/histograms_real_pairs_pbpb_20YY_single_mu4_no_trg_plots_nominal.root`
   (all three years — crossx & R_AA are always year-combined)

These regenerate the OS and SS signal-region histograms
(`h2d_crossx_…_w_signal_cuts` + `h2d_ss_…` for pp; `h3d_op/ss_…_vs_centr…` for
PbPb) that feed crossx plots **and** R_AA.

**B. Truth signal acceptance (required if acceptance/α is reported):**
3. `run_rdf_pythia_truth.sh` (modes: private, nonprivate 5.36, nonprivate 5.02)
   → `histograms_pythia_*_no_data_resonance_cuts.root` (+ `_near_away_divided`)
   — regenerates `h2d_sig_accept_{num,denom}_pt_eta` (and `_pt_150_eta`).

**C. MC fullsim reco-eff / det-response — code-consistency now, numerically
inert until real MC (see §4):**
4. `pipeline_pythia_fullsim_overlay.sh` (hijing / zmumu / data) — PbPb reco-eff.
5. `pipeline_powheg_fullsim_single_muon.sh` — **single-muon, NO rerun** (does
   not use the pair signal cut).

---

## 4. Reco-eff / det-response nuance (read before rerunning §3.C)

The nominal crossx currently applies a **placeholder** reco efficiency (Run 2
TF1, pT·q·η, ΔR-independent) — **not** the fullsim pair ε_reco. So:

- Removing/altering the signal cut changes the **set of pairs filled** in data
  crossx and the **truth acceptance** → §3.A and §3.B genuinely change results.
- The fullsim/overlay `pass_signal_*` definitions (§2) must be edited for
  consistency, but **re-deriving them does not change today's nominal crossx/
  R_AA** (placeholder is used). Treat §3.C as *bookkeeping until the real Run 3
  MC pair ε_reco(pair pT, pair η, ΔR) and det-response replace the placeholder*,
  at which point §3.C becomes result-affecting. Flagged in `docs/placeholder.md`.

---

## 5. [RERUN PLOTTING / FITTING] — consumers of the above

**Data crossx (pipeline stages, `pipelines/pipeline_{pp,pbpb}_crossx.sh`):**
- `plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx`
  → `plots/single_b_analysis/pp24/pp24_crossx_*.png`
- `plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx` (combined yrs)
  → `plots/single_b_analysis/pbpb_23_24_25_combined*/{TAA_weighted,counts}/*.png`
- `plotting_codes/single_b_analysis/plot_crossx_trig_corr_sanity.C` → sanity plots
- `plotting_codes/single_b_analysis/plot_crossx_reco_eff_stages.C`
  → `plots/sanity_check_crossx/` (raw/unfolded/reco/reco+trig stage overlays)

**R_AA (run manually, not in pipeline):**
- `RAA_plotting.cxx` (mode 6: PbPb 23+24+25 vs pp24, OS−SS) → `plots/single_b_analysis/RAA/*.png`

**MC–data comparison (PP crossx pipeline stage 8, optional):**
- `plotting_codes/mc_data_compr/plot_mc_data_compr.cxx`
  (+ `plot_mc_data_2D_hists_and_1D_proj.cxx`) → `plots/mc_data_compr/*.png`

**Signal acceptance / cutflow:**
- `plot_sig_accept_cutflow_above_60GeV.cxx` → `plots/single_b_analysis/{pythia,powheg}/*_sig_accept_above_60GeV*.png`
- `plot_signal_acceptance_pythia.cxx`, `plot_signal_acceptance_powheg.cxx`
  (`SignalAcceptancePlotter`) → 2D α(pair pT, pair η) + α-vs-pT-by-η subplots
  in `plots/single_b_analysis/{pythia,powheg}/`
  *(the 2D α(pair pT, pair η) plot itself is on the near-term to-add list.)*

**Reco-eff / det-resp plots (only meaningful once §3.C is result-affecting):**
- `plotting_codes/reco_effcy/PythiaFullsimRecoEffPlotter.cxx`,
  `plot_reco_effcy_pythia_fullsim_pp24.cxx`, det-response plotters.

**Internal note:** any synced figure embedding the above (crossx, R_AA,
acceptance) goes stale → re-run `/sync-note-figures` then `/check-note-sync`
(Gate G4). Figures live in `IntNotes/figures/` via the provenance manifest.

---

## 6. Ordered execution checklist (copy per change)

```
[ ] Edit cut in all §2 files (+ cutflow kCuts; + ΔR-axis floors if ΔR changes)
[ ] /review-analysis-code on the edits (quote analysis_overview §2 signal region)
[ ] Recompile (ACLiC) PP, PbPb, PythiaTruth, fullsim/overlay classes
[ ] Rerun data crossx: pp24 + pbpb23/24/25            (§3.A)
[ ] Rerun pythia truth acceptance (3 modes)           (§3.B)
[ ] (Real MC only) rerun fullsim/overlay reco-eff      (§3.C/§4)
[ ] Replot: pp crossx, pbpb crossx, sanity, stages    (§5)
[ ] Rerun RAA_plotting (mode 6)                        (§5)
[ ] (optional) MC–data comparison                      (§5)
[ ] Replot signal acceptance + cutflow                 (§5)
[ ] /review-plot on regenerated plots
[ ] /sync-note-figures + /check-note-sync              (§5)
[ ] Update docs (analysis_overview §2 if region changed; placeholder.md; roadmap)
```

## 7. Use for systematics

For a selection **systematic**, run §2–§6 with the varied cut into a **separate
output tag / backup dir** (do not overwrite nominal), then compare the varied
crossx/R_AA against nominal. The same blast radius applies; only data crossx
(§3.A) + truth acceptance (§3.B) + their plots are needed for a pure
acceptance/selection systematic (reco-eff placeholder is ΔR-independent, so it
cancels in the variation today).
