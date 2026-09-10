# Placeholder registry

**Purpose.** Single list of every value/correction currently standing in as a
**placeholder** in the analysis, so none silently leaks into final results or the
internal note. Each entry says what is used now, what is actually needed, where it
lives in the code, and what must be disclosed in the note. The roadmap
(`docs/tracking/analysis_roadmap_2026_06.md`, §Q2 and §Q4) is the authoritative
status source; this file is the consolidated quick index.

A symlink to this file exists in the IntNote repo (`IntNotes/placeholder.md`) so
the note writer always has the disclosure list at hand.

Last updated: 2026-08-24. (Items 3 and 7 rescoped: the pp side of BOTH is no longer a
placeholder. See `docs/tracking/intnote_trig_reco_eff_sections.md` for the code evidence.)

---

## Summary table

| # | Placeholder | Scope | Used now | Needs | Code / data location | Roadmap |
|---|-------------|-------|----------|-------|----------------------|---------|
| 1 | **Centrality classification** | PbPb 2024 (provisional own bins), 2025 (2023 thresholds) | 2024: `GetCentralityPbPb2024` with `FCal_ET_Bins_PbPb2024`; 2025: the 2023 Glauber FCal-ET thresholds. **The cross-year FCal scaling is documented but implemented in NO source file** (`grep fcal_scale` over `*.c *.C *.h *.cxx` returns nothing). | Official 2024/2025 Glauber centrality calibration | `MuonPairPbPb.h:147` (2024), `:150` (2025) | Q2.3, ledger step 6 |
| 2 | **⟨T_AA⟩ values** | PbPb 2024, 2025 | 2023 ⟨T_AA⟩ = {26.1428, 20.3241, 14.0502, 8.5074, 3.7733, 0.6716} mb⁻¹ | Official 2024/2025 ⟨T_AA⟩ + uncertainties + citable ref | `PbPbBaseClass.h` `make_crossx_factors_pbpb_2024/2025`; source `IntNotes/data/centrality/TaaValues2023.txt` | Q2.3 |
| 3 | **Reconstruction efficiency** — Run 2 single-muon ε_reco proxy (ε₁·ε₂), standing in for the proper 3D pair ε_reco(pair pT, pair η, dR). **pp24 IS NO LONGER A PLACEHOLDER (2026-08-18):** it uses the measured pp24-fullsim 3D pair efficiency (`pair_reco_eff_pp24_full.root`, inclusive 0.6620 Tight / 0.7301 Medium) via `Utilities/PairRecoEffEvaluator.h`. PbPb only. | PbPb (all yrs) | **Applied** in crossx as a correction STAGE (`_corr_unfolded_reco[_trig]` hists); **PbPb = colleague's EXACT Run 2 Medium-μ TF1 fits** (`MuonRecoEffcyRun2MC_medium.root`, evaluated at exact pT; 2026-06-19); **pp is NOT a placeholder** (measured pp24-fullsim 3D pair eps_reco since 2026-08-18) | Full Pythia fullsim HIJING-overlay (r17662) + full pp24 fullsim → proper 3D pair ε_reco | `RDFBasedHistFilling{Data,PP,PbPb}` (`EvaluateSingleMuonRecoEffPlaceholder`, `CorrectionStages.h`); `EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root` | Q4, ledger 12; tracking `reco_eff_placeholder_run2.md` |
| 4 | **σ_PbPb (total hadronic)** | PbPb (all yrs) | 7.8 b at 5.36 TeV — **unvalidated guess** | Citable 5.36 TeV reference | `PbPbBaseClass.h` (guess comment on 2023 helper) | Q2.2 |
| 5 | **Luminosity uncertainty** | all years | not set | Official per-year Run 3 lumi uncertainty | — (note systematics) | Q2.1 |
| 6 | **PbPb 2026 lumi + GRL** | PbPb 2026 | placeholder (data not yet in skim) | 2026 data in skim, then lumi/GRL | `IntNotes/analysis_metadata.md` | Q2.1/Q2.5 |
| 7 | **dR trigger-correlation correction** — **pp24 IS NO LONGER A DUMMY (2026-08-18):** a measured MC ε_dR is applied (`DrCorrectionCrossxEvaluator`, expo / opposite sign / last-two-pair-pT-bins-merged; inclusive 0.818 at dR→0). **PbPb still applies the BARE union with NO ε_dR at all** (`RDFBasedHistFillingPbPb.cxx:1008`) — the dR-corrected union in `mc_trigger_efficiency.md` §2 is documented but not implemented. | PbPb | dummy ε_dR ≡ 1; **and the former ΔR>0.05 SIGNAL cut was removed 2026-06-22 (interim nominal)** — both pending the MC-based dR trigger correction that decides whether any ΔR cut is needed | Fullsim overlay with trigger sim | crossx pipeline (dR corr = 1; no ΔR signal cut); `docs/tracking/remove_dr_cut_signal_selection.md` | Q4 |
| 8 | **Detector response / unfolding** | PbPb + pp | test-sample shapes | Full Pythia fullsim (pp24 + overlay) | unfolding inputs | Q4 |
| 9 | **Δp/p significance + template-fit purity** | PbPb + pp | **DEFERRED** (2026-06-22) — not built; fake-muon purity treated as ~flat (Run 2 dimuon note >98%) for the preliminary chain | **Future TODO needing MC:** signal/bkg Δp/p templates from fullsim + a π/K-enriched MC (Run 3 equivalent not yet identified). Build the fake fit (at RECO level — Δp/p is intrinsically reco) once that MC exists | (to be built) | Q4; tracking `low_mass_dimuon_template_fit.md` §3a |
| 10 | **MCP scale factors** | MC | omitted (pure MC-driven) | Run 3 HI MCP recommendation | — | Q2.7 |

---

## Detail on the two primary placeholders (2024 / 2025)

### 1–2. Centrality & ⟨T_AA⟩ for PbPb 2024 and 2025

Official Glauber centrality calibrations and ⟨T_AA⟩ for 2024 and 2025 **do not
exist yet**. Until they do:

- **Centrality:** **2024** events are classified with their own provisional
  `FCal_ET_Bins_PbPb2024` via `MuonPairPbPb::GetCentralityPbPb2024`
  (`MuonPairPbPb.h:147`); **2025** events fall back to the **2023** Glauber FCal-ET
  thresholds via `GetCentralityPbPb2023` (`:150`, *"use pbpb2023 thresholds until
  pbpb2025 are derived"*). **The per-year FCal cross-year scale factors
  (`fcal_scale_pbpb_20YY.root`) are documented in several docs but appear in NO
  source file** — `grep fcal_scale` over `*.c *.C *.h *.cxx` returns nothing, so no
  cross-year FCal rescaling is actually applied. Corrected 2026-08-24 against the
  code; the previous text claimed 2023 thresholds for both years and a scaling that
  is not implemented.
- **⟨T_AA⟩:** the crossx normalization for 2024, 2025 and 2026
  (`make_crossx_factors_pbpb_2024/2025/2026` in `PbPbBaseClass.h`) uses the **2023**
  ⟨T_AA⟩ array as a placeholder; only the per-year luminosity is year-specific
  (2024 = 0.85112 [GRL ≥489703, corrected 2026-06-19], 2025 = 2.59933,
  2026 = 2.62316 nb⁻¹ [GRL `physics_HI2026_50ns_noIBL.xml`, added 2026-09-10]).
  Source values: `IntNotes/data/centrality/TaaValues2023.txt`.

**Affects:** centrality binning, cross-section normalization, and R_AA for 2024,
2025 and 2026. **Note disclosure (required):** state explicitly that 2024/2025/2026
centrality and ⟨T_AA⟩ are 2023 placeholders pending official calibrations.

### 3. Reconstruction efficiency

**pp24 IS NO LONGER A PLACEHOLDER (2026-08-18).** The pp cross-section applies the measured
**3D pair** efficiency ε_reco(pair pT, pair η, dR) from the pp24-condition Pythia8 fullsim FULL
production: 8 × 9 × 4 = 288 cells, built by
`plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C` into
`~/usatlasdata/pythia_fullsim_full_sample/pair_reco_eff_pp24_full.root` and read at fill time by
`Utilities/PairRecoEffEvaluator.h` (`RDFBasedHistFillingPP.cxx:396-399`). Inclusive value
**0.6620 Tight / 0.7301 Medium**; it is a FIDUCIAL efficiency (the detector-gap cut sits on both
the truth denominator and the reco numerator), so ε_acc = 0.9133 is a separate, unapplied factor.
Figures: `<sample>/plots/pp24_reco_effcy_plots/{tight,medium}/applied/`.

**Pb+Pb IS STILL A PLACEHOLDER.** The proper Run 3 correction is the same 3D pair efficiency from
the Pythia+HIJING overlay, but only a small TEST overlay production exists, so Pb+Pb instead uses
a Run 2 **single-muon** product proxy ε_reco(p_a)·ε_reco(p_b), with **no dR dependence at all**:

- the colleague's **exact Run 2 fits** from the Run 2 dimuon analysis (ATL-COM-PHYS-2021-1094) —
  `EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_{tight,medium}.root`, logistic TF1
  `tf1_eff_fit_cent{C}_eta{E}`, single-muon ε_reco(pT, q·η) per centrality (0–10…60–80%), HIJING
  overlay 5.02 TeV. **BOTH working points are built and the consumer picks the WP-matched key**
  (`RDFBasedHistFillingData.cxx:783`); the muon pT is clamped to [4, 19] GeV at the point of use.
- Written into `EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root` by
  `plotting_codes/reco_effcy/build_run2_reco_eff_placeholder.C` (63 TF1s per WP).
  `EvaluateSingleMuonRecoEffPlaceholder` loads them and dispatches by centrality; applied in
  `RDFBasedHistFillingPbPb.cxx:1010-1023` as ε₁·ε₂, floored at 0.05 before inversion.
- The pp TGraphs in the same file (HF R_AA note HION-2019-58 Fig. 31, Medium) are reached only
  through the centrality < 0 sentinel, which the Pb+Pb path never takes; they are effectively
  dead now that pp has its own measured map.

**Why this placeholder is poor (must be replaced):** (1) Run 3 muon reco is
expected considerably better than Run 2 (New Small Wheel + other Run 3 muon
upgrades). (2) ε₁·ε₂ does **not** factorize for our signal — the two muons are
nearby, so pair reconstruction is correlated; the proper correction is the 3D
pair ε_reco(pair pT, pair η, dR). For 2024/2025 the eventual overlay ε_reco also
inherits the 2023 centrality placeholder (item 1).

**Nominal-result note (updated 2026-06-16):** the reco placeholder is **now in
the nominal** — `w_reco` is folded into the nominal corrected weight
(`weight_for_RAA_trig_corr` = `weight_for_RAA·w_reco·w_trig`;
`crossx_weight_trig_corr` = `crossx_weight·w_reco·w_trig`), so all standard
crossx histograms and the R_AA 3D input (`h3d_op_crossx_..._vs_centr...`) are
reco+trig corrected (== the validated `_corr_unfolded_reco_trig` stage). Nominal
crossx plots reran (pp24, pbpb_23_24_25_combined). **R_AA (task_06 DONE
2026-06-16):** `RAA_plotting.cxx` modernized to read the RDF crossx outputs
(combined PbPb 23+24+25+26 vs pp24), with SS signal-region histos added to the RDF
(`h3d_ss_...`, `h2d_ss_...`) so R_AA does the OS−SS combinatorial subtraction;
reco-corrected R_AA plots (vs pair pT/η/centrality) in
`dimuon_data/plots/single_b_analysis/RAA/`. See `docs/tracking/raa_from_rdf_crossx.md`.
**Still placeholder/preliminary** — do not quote as final; pre-reco trig-only
nominal preserved in `dimuon_data/crossx_hist_backup_20260616_pre_reco_nominal/`.

> **Year-combination normalization (RESOLVED 2026-06-16):** combined-year Pb+Pb
> results now use the **luminosity-weighted average** `Σ(L_y·h_y)/ΣL_y` (HF R_AA
> note HION-2019-58 §4.1 Eq.3), via single-source `Utilities/PbPbSampledLumi.h`,
> in the crossx combined plotter, R_AA, and the stage plotter. R_AA absolute
> scale is now physical (~0.1–0.9). Still preliminary due to the 2023 ⟨T_AA⟩
> placeholder (items 1–2) and σ_PbPb guess (item 4).

**Needs:** full Pythia fullsim HIJING-overlay (r17662) + full pp24 fullsim →
proper 3D pair ε_reco. **Note disclosure:** quote reco-eff and any
efficiency-derived result as placeholder/preliminary.

---

For items 4–10 see roadmap §Q2 / §Q4 for the dummy strategy and current status.
