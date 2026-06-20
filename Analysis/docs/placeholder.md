# Placeholder registry

**Purpose.** Single list of every value/correction currently standing in as a
**placeholder** in the analysis, so none silently leaks into final results or the
internal note. Each entry says what is used now, what is actually needed, where it
lives in the code, and what must be disclosed in the note. The roadmap
(`docs/tracking/analysis_roadmap_2026_06.md`, §Q2 and §Q4) is the authoritative
status source; this file is the consolidated quick index.

A symlink to this file exists in the IntNote repo (`IntNotes/placeholder.md`) so
the note writer always has the disclosure list at hand.

Last updated: 2026-06-15.

---

## Summary table

| # | Placeholder | Scope | Used now | Needs | Code / data location | Roadmap |
|---|-------------|-------|----------|-------|----------------------|---------|
| 1 | **Centrality classification** | PbPb 2024, 2025 | 2023 Glauber FCal-ET thresholds + cross-year FCal scaling | Official 2024/2025 Glauber centrality calibration | `MuonPairPbPb.h:144,150` (`GetCentralityPbPb2023`); `fcal_scale_pbpb_20YY.root` | Q2.3, ledger step 6 |
| 2 | **⟨T_AA⟩ values** | PbPb 2024, 2025 | 2023 ⟨T_AA⟩ = {26.1428, 20.3241, 14.0502, 8.5074, 3.7733, 0.6716} mb⁻¹ | Official 2024/2025 ⟨T_AA⟩ + uncertainties + citable ref | `PbPbBaseClass.h` `make_crossx_factors_pbpb_2024/2025`; source `IntNotes/data/centrality/TaaValues2023.txt` | Q2.3 |
| 3 | **Reconstruction efficiency** — Run 2 single-muon ε_reco proxy (ε₁·ε₂), standing in for the proper 3D pair ε_reco(pair pT, pair η, dR) | PbPb (all yrs) + pp | **Applied** in crossx as a correction STAGE (`_corr_unfolded_reco[_trig]` hists); **PbPb = colleague's EXACT Run 2 Medium-μ TF1 fits** (`MuonRecoEffcyRun2MC_medium.root`, evaluated at exact pT; 2026-06-19); pp = Medium-μ digitized from HF R_AA Fig.31 | Full Pythia fullsim HIJING-overlay (r17662) + full pp24 fullsim → proper 3D pair ε_reco | `RDFBasedHistFilling{Data,PP,PbPb}` (`EvaluateSingleMuonRecoEffPlaceholder`, `CorrectionStages.h`); `EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root` | Q4, ledger 12; tracking `reco_eff_placeholder_run2.md` |
| 4 | **σ_PbPb (total hadronic)** | PbPb (all yrs) | 7.8 b at 5.36 TeV — **unvalidated guess** | Citable 5.36 TeV reference | `PbPbBaseClass.h` (guess comment on 2023 helper) | Q2.2 |
| 5 | **Luminosity uncertainty** | all years | not set | Official per-year Run 3 lumi uncertainty | — (note systematics) | Q2.1 |
| 6 | **PbPb 2026 lumi + GRL** | PbPb 2026 | placeholder (data not yet in skim) | 2026 data in skim, then lumi/GRL | `IntNotes/analysis_metadata.md` | Q2.1/Q2.5 |
| 7 | **dR trigger-correlation correction** | PbPb | dummy ε_dR ≡ 1 | Fullsim overlay with trigger sim | crossx pipeline (dR corr = 1) | Q4 |
| 8 | **Detector response / unfolding** | PbPb + pp | test-sample shapes | Full Pythia fullsim (pp24 + overlay) | unfolding inputs | Q4 |
| 9 | **Δp/p significance + template-fit purity** | PbPb + pp | placeholder; framework only | Signal/bkg templates from fullsim + π/K MC | (to be built) | Q4 |
| 10 | **MCP scale factors** | MC | omitted (pure MC-driven) | Run 3 HI MCP recommendation | — | Q2.7 |

---

## Detail on the two primary placeholders (2024 / 2025)

### 1–2. Centrality & ⟨T_AA⟩ for PbPb 2024 and 2025

Official Glauber centrality calibrations and ⟨T_AA⟩ for 2024 and 2025 **do not
exist yet**. Until they do:

- **Centrality:** 2024 and 2025 events are classified with the **2023** Glauber
  FCal-ET thresholds (`FCal_ET_Bins_PbPb2023`) via
  `MuonPairPbPb::GetCentralityPbPb2023`, after applying the per-year **FCal
  cross-year scale factors** (`fcal_scale_pbpb_20YY.root`) so the 2024/2025 FCal
  scale is mapped onto the 2023 reference. Code: `MuonPairPbPb.h:144` (2024),
  `:150` (2025, *"use pbpb2023 thresholds until pbpb2025 are derived"*).
- **⟨T_AA⟩:** the crossx normalization for 2024 and 2025
  (`make_crossx_factors_pbpb_2024/2025` in `PbPbBaseClass.h`) uses the **2023**
  ⟨T_AA⟩ array as a placeholder; only the per-year luminosity is year-specific
  (2024 = 0.85112 [GRL ≥489703, corrected 2026-06-19], 2025 = 2.59933 nb⁻¹). Source values:
  `IntNotes/data/centrality/TaaValues2023.txt`.

**Affects:** centrality binning, cross-section normalization, and R_AA for 2024
and 2025. **Note disclosure (required):** state explicitly that 2024/2025
centrality and ⟨T_AA⟩ are 2023 placeholders pending official calibrations.

### 3. Reconstruction efficiency

The proper Run 3 correction is the **3D pair** efficiency ε_reco(pair pT, pair η,
dR) (method settled; index-mapping/HIJING bugs fixed; residual low-pT overlay
deficit shown physical — investigation Steps 23–24), but only **test** MC exists
(bug-fixed r17618 60k-evt + r17662 10k-evt overlay; pp24 fullsim test only), so
it is not yet derivable for real.

**Interim placeholder NOW APPLIED (2026-06-15):** a Run 2 **single-muon** ε_reco
**product proxy** ε_reco(p_a)·ε_reco(p_b), per the Run 2 dimuon-note treatment
(`w⁻¹ = ε_trig·ε_reco(p_a)·ε_reco(p_b)`). Placeholder source values (Medium muons):

- **PbPb (updated 2026-06-18):** the colleague's **EXACT Run 2 Medium-μ fits**
  actually used in the Run 2 note — `EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_medium.root`
  (logistic TF1 `tf1_eff_fit_cent{C}_eta{E}`), single-muon ε_reco(pT, q·η) per
  centrality (0–10…60–80%), HIJING overlay 5.02 TeV. This **replaces** the earlier
  eyeball digitization of App. F.2 (same source physics, real fitted numbers).
- **pp:** Run 2 HF-muon R_AA note **HION-2019-58 / arXiv:2109.00411 Fig. 31** —
  data-driven Medium-μ ε_reco^pp(pT), barrel (|η|<1.05) + endcap (1.3<|η|<2.1).
  This is the source the dimuon note itself cites for its pp efficiencies; chosen
  over the peripheral-PbPb fallback. (See `reading`/memory `project_pp_reco_eff_placeholder`.)

Written into `EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root`
(`plotting_codes/reco_effcy/build_run2_reco_eff_placeholder.C`): PbPb stored as the
colleague's Medium **TF1 fits** (`tf1_reco_eff_medium_pbpb_ctr{lo}_{hi}_q_eta_{suffix}`,
63 of them) and **evaluated at the exact muon pT** in the lookup — no resampling
(2026-06-19; pp still eyeball HF R_AA Fig.31 as 2 TGraphs). `EvaluateSingleMuonRecoEffPlaceholder`
loads TF1 (PbPb) + TGraph (pp) and dispatches by side. Applied in the crossx RDF as a
correction **stage** (`CorrectionStages.h`): histograms saved at each stage —
`_corr_raw` → `_corr_unfolded` (identity) → `_corr_unfolded_reco` →
`_corr_unfolded_reco_trig`; before/after 3-line plots in
`dimuon_data/plots/sanity_check_crossx/*_reco_eff_stages_*`.

The reco correction is also folded into the **pp generic analysis weight**
(`FillHistogramsGeneric` → `generic_weight_col = w_reco_trig`), so the **MC-data
comparison** (POWHEG/Pythia vs pp24) reflects it. **INVARIANT:** regenerate the
MC-data comparison (`plot_mc_data_compr.cxx`; pp crossx pipeline stage 8) after
ANY pp efficiency / detector-response / unfolding change. The PbPb crossx also
has a genuine **differential cross-section** dσ/dp_T = (1/L)·dN [nb/GeV]
(`differential_crossx/` plots), distinct from the T_AA-weighted R_AA input.

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
(combined PbPb 23+24+25 vs pp24), with SS signal-region histos added to the RDF
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
