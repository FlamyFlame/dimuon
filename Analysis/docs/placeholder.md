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
| 3 | **Pair reconstruction efficiency** ε_reco(pair pT, pair η, dR) | PbPb (all yrs) + pp | **Not yet applied** to crossx; method settled, only test MC exists | Full Pythia fullsim HIJING-overlay (r17662) + full pp24 fullsim | `RDFBasedHistFilling*` (no reco-eff weight yet) | Q4, ledger step 12 |
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
  (2024 = 1.59663, 2025 = 2.59933 nb⁻¹). Source values:
  `IntNotes/data/centrality/TaaValues2023.txt`.

**Affects:** centrality binning, cross-section normalization, and R_AA for 2024
and 2025. **Note disclosure (required):** state explicitly that 2024/2025
centrality and ⟨T_AA⟩ are 2023 placeholders pending official calibrations.

### 3. Reconstruction efficiency

The pair reco efficiency ε_reco(pair pT, pair η, dR) is methodologically settled
(3D pair efficiency; index-mapping/HIJING bugs fixed; residual low-pT overlay
deficit shown physical — investigation Steps 23–24) but **is not yet applied** to
the cross-section/R_AA filling, because only **test** MC samples exist (bug-fixed
r17618 60k-evt and r17662 10k-evt overlay; pp24 fullsim test only). Current
crossx/R_AA figures are therefore **preliminary** — reco-eff-uncorrected.

For 2024/2025 there is an additional layer: ε_reco is derived from the HIJING
overlay whose **own** centrality is classified with the 2023 placeholder above,
so the 2024/2025 efficiency-vs-centrality inherits the centrality placeholder
until official calibrations land.

**Needs:** full Pythia fullsim HIJING-overlay production (r17662 recommended) and
full pp24 fullsim. **Note disclosure:** quote reco-eff and any efficiency-derived
result as placeholder/preliminary.

---

For items 4–10 see roadmap §Q2 / §Q4 for the dummy strategy and current status.
