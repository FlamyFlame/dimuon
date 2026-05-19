# Trigger Efficiency Plot Structure & Physics Procedures Review

## Objective

1. Establish a clean, hierarchical directory structure under `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/trig_effcy_plots/` for all trigger efficiency outputs.
2. Add `isBNL` toggle to the TrigEffPlotter classes (currently Mac-only paths).
3. Verify that the physics procedures in Pipelines 2 and 3 are correct — in particular the role-swap summing, cross-term definition, and 2mu4 vs mu4 distinction.

## Context

User requirements:
- mu1tag/mu2tag histograms should always be **combined** (role-swapped sum), never plotted separately
- For dR corrections to 2mu4 (PP), there should be no tag/probe structure: the pair-level 2mu4 trigger decision is used directly
- 2mu4 should have no cross term: the final correction is ε₁^{no-corr} × ε₂^{no-corr} × ε_dR(dR)
- Both no-correlation and dR-correction plots should be available in centrality bins

## Proposed Directory Structure (INITIAL — WRONG, needs revision)

This was the first proposed structure before the user's corrections. Kept for reference.

```
plots/trig_effcy_plots/
├── trg_effcy_no_corr/
│   ├── pbpb_2023/
│   │   ├── pair_obs/
│   │   │   ├── op_and_sig/
│   │   │   └── sepr/
│   │   ├── single_muon_effcy/
│   │   │   ├── ctr_inclusive/
│   │   │   └── ctr_dep/
│   │   ├── single_muon_effcy_mu4noL1/
│   │   │   ├── ctr_inclusive/
│   │   │   └── ctr_dep/
│   │   └── pT_fits/
│   │       ├── ctr_inclusive/
│   │       └── ctr_dep/
│   ├── pbpb_2024/
│   ├── pbpb_2025/
│   └── pp_2024/
├── dR_corr_to_mu4/                      # PbPb
│   ├── pbpb_2023/
│   │   ├── ctr_inclusive/
│   │   │   ├── mu1tag/       ← WRONG: should be combined
│   │   │   ├── mu2tag/       ← WRONG: should be combined
│   │   │   └── cross_term/
│   │   └── ctr_dep/
│   │       ├── ctr0_5/
│   │       │   ├── mu1tag/   ← WRONG
│   │       │   ├── mu2tag/   ← WRONG
│   │       │   └── cross_term/
│   │       └── ...
│   ├── pbpb_2024/
│   └── pbpb_2025/
├── dR_corr_to_2mu4/                     # PP
│   └── pp_2024/
│       ├── mu1tag/            ← WRONG: no tag for 2mu4
│       ├── mu2tag/            ← WRONG: no tag for 2mu4
│       └── cross_term/        ← WRONG: no cross term for 2mu4
└── dR_corr_to_mu4_cross_terms/          # PbPb cross-term (separate)
    └── ...                    ← WRONG: should be under dR_corr_to_mu4
```

Problems with this structure:
1. mu1tag/mu2tag are shown as separate directories — they should always be combined (summed)
2. 2mu4 has tag/probe structure — incorrect, 2mu4 is a pair-level trigger with no tag
3. 2mu4 has a cross-term — incorrect, the 2mu4 dR correction has only one term
4. mu4 cross-term is a separate top-level — should be under dR_corr_to_mu4

## isBNL Requirement

The TrigEffPlotter classes (`TrigEffPlotterPbPb.cxx`, `TrigEffPlotterPP.cxx`) currently use **Mac-only** paths:
- `data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pbpb_2023/"`
- Plots: `data_dir + "trig_effcy_plots/..."` → only on local Mac

The `SingleMuEffcyPtTurnOnFitter.cxx` uses **BNL paths**:
- `data_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20YY/"`
- Fitting plots: `data_dir + "trg_effcy_pT_fitting_to_fermi_plus_log/..."`

Need: an `isBNL` bool config in `TrigEffPlotterBaseClass` (default false). When true, rewrite `data_dir` to BNL path and redirect plot output to the centralized `plots/trig_effcy_plots/` structure.

## Physics Procedures — Current Implementation

### Pipeline 2: No-correlation single-muon trigger efficiency ε^{no-corr}(pT, q·η)

**Goal:** Measure the probability that a single muon fires a given trigger (mu4, mu4_mu4noL1, 2mu4), as a function of the muon's pT and q·η, with no requirement on proximity to another muon.

**Method: Tag-and-probe on muon pairs**

For each muon pair:
- One muon ("1st") is required to fire single-mu4 (the "tag")
- The other muon ("2nd") is the probe — its pT, q·η, φ are the kinematic variables
- Test whether the probe passes the trigger of interest

The code loops over **both role assignments** (`_mu1passmu4`, `_mu2passmu4`): in each, the "1st" muon fires mu4, and the "2nd" is the probe. This effectively doubles statistics and symmetrizes.

**Histograms filled (per role × per pair_sign × per musign × per trigger):**
- 2D: `h_pt2nd_vs_q_eta2nd{pair_sign}{mu4sel}{musign}{trg}{bias}` (pT vs q·η of the 2nd muon)
- Also 1D pair observables, 3D (with φ), etc.

**Post-processing:**
1. **Role-swap summation** (levels {0,1}): The code sums `_mu1passmu4` and `_mu2passmu4` histograms for each (pair_sign, musign, trg, bias). After this, the "probe" is a generic muon from the pair, not specifically muon 1 or 2. **This is correct — the role swap is always combined.**
2. **Efficiency ratio:** For each q·η projection bin:
   - ε(pT | q·η) = N(probe passes trigger && probe in q·η bin) / N(probe in q·η bin)
   - Specifically: `h_pt2nd_vs_q_eta2nd_{trg}_sepr` / `h_pt2nd_vs_q_eta2nd_mu4_sepr`
   - The denominator is "probe passes mu4" (same as the tag condition), with `_sepr` = dR > 0.8

**Triggers measured:**
- `_mu4`: denominator condition = tag passes mu4. Numerator condition = (none; mu4 in denominator is the tag). This gives ε(mu4) = N(all pairs where tag passes mu4) / N(all pairs where tag passes mu4) = 1. → Actually `_mu4` is the denominator for other triggers.
- `_mu4_mu4noL1`: numerator = pair passes mu4_mu4noL1 AND probe fires the noL1 leg
- `_2mu4`: numerator = pair passes 2mu4
- `_2mu4_AND_mu4_mu4noL1`: numerator = pair passes both

So what's measured:
- **ε(mu4_mu4noL1 | mu4, dR > 0.8)** = P(probe fires mu4noL1 | tag fires mu4, dR > 0.8) as function of (pT_probe, q·η_probe)
- **ε(2mu4 | mu4, dR > 0.8)** = P(pair passes 2mu4 | tag fires mu4, dR > 0.8) — this is P(probe fires mu4 | tag fires mu4, dR > 0.8), i.e., the single-muon mu4 efficiency of the probe

The `_sepr` (dR > 0.8) requirement ensures no dR-dependent correlation contaminates the "no-correlation" efficiency.

**pT turn-on fitting:**
- The 2D efficiency ratio is projected onto pT in q·η bins
- Fitted with Fermi+log function per (centrality, musign, q·η bin)
- Output: TF1s stored in `single_mu_effcy_pT_fit.root`
- Centrality-binned fits are produced (PbPb)

### Pipeline 3: Inverse-weighted dR correlation correction ε_dR(dR)

**Goal:** Measure how much the trigger efficiency changes as a function of dR (and other pair observables), relative to the no-correlation baseline.

**Physics model:**
```
P(pair passes trigger | dR) = ε₁^{no-corr}(pT₁, q·η₁) × ε₂^{no-corr}(pT₂, q·η₂) × ε_dR(dR)
```

To extract ε_dR, divide out the no-correlation part via inverse weighting.

**Current PbPb implementation (three terms):**

**Term 1: Single-muon, mu1 = tag (mu2 = probe)**
- Denominator: all pairs where mu1 fires mu4 AND ε₂ is defined (valid2). Weight = 1. Fill dR histogram.
- Numerator: pairs where mu1 fires mu4 AND mu2 fires mu4. Weight = 1/ε₂^{no-corr}. Fill dR histogram.
- Ratio = Σ(1/ε₂) for pairs where both fire / Σ(1) for pairs where mu1 fires

**Term 2: Single-muon, mu2 = tag (mu1 = probe)**
- Same as Term 1 but with mu1/mu2 roles swapped.
- Ratio = Σ(1/ε₁) for pairs where both fire / Σ(1) for pairs where mu2 fires

**Term 3: Cross-term**
- Denominator: all pairs where both ε₁ and ε₂ are defined (valid_both). Weight = 1.
- Numerator: pairs where BOTH mu1 and mu2 fire mu4. Weight = 1/(ε₁ × ε₂).
- Ratio = Σ(1/(ε₁ε₂)) for pairs where both fire / Σ(1) for all valid pairs

**What these terms physically represent:**

Term 1 (mu1 tag): At a given dR, among pairs where mu1 fires mu4, what fraction of the time does mu2 also fire mu4, after dividing out ε₂^{no-corr}? If there's no dR correlation, this ratio = 1. If proximity suppresses triggering, ratio < 1 at small dR.

Term 2 (mu2 tag): Same but with swapped roles.

Terms 1+2 should be **combined** (summed numerators and summed denominators separately, then divided) to get the combined single-muon dR correction. The current code produces them separately in `MakeAndWriteDRTrigEffGraphsHelper` and the plotting code `plot_dR_trig_corr.C` plots them individually.

**Issue: mu1tag and mu2tag are never combined in the current code.** The ratio graphs `g_DR_*_mu1tag_*_divided` and `g_DR_*_mu2tag_*_divided` exist separately. They should be summed before taking the ratio to get the symmetric single-muon correction.

Term 3 (cross): Measures the pair-level trigger correlation. If ε(both fire) = ε₁ε₂ (no correlation), the inverse-weighted ratio = 1. Deviations from 1 indicate genuine two-body trigger correlation (L1 sector sharing, etc.).

**Current PP implementation:**
- `FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` is an **empty stub** — not implemented yet.
- Design decision D2 says: "For PP24 Pipeline 3, the inverse weighting fills numerator for pairs passing 2mu4 trigger (not 'both muons individually pass mu4'), weighted by 1/(ε₁ × ε₂)."

## Issues

### Issue 1: Single-muon dR correction has WRONG denominator filter (mu4 requirement on "other" muon)

**Root cause identified:** The code applies a mu4 trigger requirement on the "other" muon in the denominator, which is physically wrong for the mu4 single-muon dR correction.

**What the current code does** (RDFBasedHistFillingPbPb.cxx lines 456–476):
```
"mu1tag" denom: Filter("m1.passmu4 && valid2")  ← WRONG: requires m1 fires mu4
"mu1tag" num:   additionally Filter("m2.passmu4"), weight = 1/ε₂
"mu2tag" denom: Filter("m2.passmu4 && valid1")  ← WRONG: requires m2 fires mu4
"mu2tag" num:   additionally Filter("m1.passmu4"), weight = 1/ε₁
```

**What's physically wrong:**

The mu4 single-muon dR correction measures: does the **mere existence** of another offline muon at angular distance dR affect P(this muon fires mu4)?

Correct: P(2nd muon passes mu4 | 1st offline muon exists at dR) = ε^{no-corr}(pT₂, q·η₂) × ε_dR(dR)

The current code measures something different: P(2nd muon passes mu4 | **1st muon passes mu4** at dR). This conditions on the other muon's trigger decision, which is the approach for mu4_mu4noL1/2mu4 corrections (where the 1st muon passing mu4 is a trigger requirement), NOT for the mu4 single-muon correction.

**Correct implementation:**

For each pair, with "2nd muon" as the probe:
- **Denominator:** `Filter("valid2")` — all pairs where ε(2nd muon) is defined. **No mu4 requirement on either muon.**
- **Numerator:** `Filter("m2.passmu4")` from the denominator population, weight = 1/ε^{no-corr}(2nd muon).
- Role-swap (each muon is "2nd" in turn): **always combined** — sum numerators, sum denominators, then divide.

**Code locations to fix:**
- PbPb.cxx line 457: `node.Filter("m1.passmu4 && valid2")` → `node.Filter("valid2")`
- PbPb.cxx line 468: `node.Filter("m2.passmu4 && valid1")` → `node.Filter("valid1")`
- Remove "tag" naming throughout — rename to role-swap notation (e.g., `_mu2probe`, `_mu1probe`)
- `MakeAndWriteDRTrigEffGraphsHelper` (Data.cxx lines 564–575): must combine the two role-swap terms before dividing (sum numerator hists, sum denominator hists, then single divide)

### Issue 1b: mu1tag/mu2tag are never combined

`MakeAndWriteDRTrigEffGraphsHelper` (Data.cxx lines 564–575) produces separate ratio graphs for `_mu1tag` and `_mu2tag`. The plotting code (`plot_dR_trig_corr.C`) also plots them individually. Per the physics procedure, the two role-swap terms must **always** be combined (sum num, sum denom, then divide) to get the symmetric dR correction.

### Issue 2: PP 2mu4 Pipeline 3 not yet implemented

`FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` is an **empty stub** in RDFBasedHistFillingPP.cxx (line 190).

**User requirement:** For PP 2mu4, the correct procedure is:
- **No tag/probe**: 2mu4 is a pair-level trigger. There is no muon that serves as a "tag."
- **No cross-term**: There is only ONE dR correction term for 2mu4.
- **Factorization:** P(pair passes 2mu4 | dR) = ε₁^{no-corr}(pT₁, q·η₁) × ε₂^{no-corr}(pT₂, q·η₂) × ε_dR^{2mu4}(dR)
- **Implementation:** Denominator = all pairs (both ε defined, wt=1). Numerator = pairs passing `pass2mu4` (wt = 1/(ε₁ × ε₂)). Ratio = ε_dR^{2mu4}(dR).

### Issue 3: Misleading `_2mu4_` in Pipeline 3 histogram suffix names

The histogram suffixes use `_2mu4_` (e.g., `_mu1tag_2mu4_denom`, `_mu1tag_2mu4_invw_num`). This "2mu4" refers to design decision D5 (the TF1 fits were derived from 2mu4-conditional efficiency), but in context it reads as if these histograms are for the 2mu4 trigger correction. They are actually for the **mu4 single-muon** dR correction. Should be renamed to `_mu4_` to match the physics meaning.

Code locations: PbPb.cxx lines 458, 463, 469, 474; Data.cxx line 566 (`"_2mu4_invw_num"`, `"_2mu4_denom"`).

### Issue 4: No unfitted 2D histogram fallback for gap regions

Per the physics procedure, for muons with q·η in gap regions (where the pT-fitted TF1 doesn't exist), the no-correlation efficiency should be evaluated from the **unfitted 2D histogram** bin content. The current code (`EvaluateSingleMuonEffcyPtFitted`, Data.cxx line 541) returns -1.0 for gap regions, which causes the muon's weight to be set to 0 (excluded). This biases the dR correction by systematically removing all muons in gap regions from the sample.

Code location: Data.cxx lines 537–551. The function has no path to fall back to the 2D histogram.

### Issue 5: Pipeline 3 has no centrality binning (was Issue 3)

`FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` only fills centrality-inclusive histograms. `MakeAndWriteDRTrigEffGraphs()` calls `MakeAndWriteDRTrigEffGraphsHelper({})` — empty categories = centrality-inclusive only.

**User requirement:** Both no-correlation efficiencies AND dR corrections should be available in centrality bins.

## Output Directory Structure Change Request

**Current problems:**
- dR correction plots in `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/dR_trig_corr_plots/` — flat, no hierarchy, wrong location
- TrigEffPlotter uses Mac-only paths (`/Users/yuhanguo/...`), no BNL equivalent
- No centralized plot repository

**User requirement:** Move all trigger efficiency plots under:
```
/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/trig_effcy_plots/
```
with hierarchical subdirectories:
1. `trg_effcy_no_corr/` — no-correlation mu4 efficiencies ε(pT, q·η), including fitting, including mu4noL1 from P(mu4_mu4noL1 | mu4 && dR > 0.8)
2. `dR_corr_to_mu4/` — PbPb dR correlation corrections to mu4
3. `dR_corr_to_2mu4/` — PP dR correlation corrections to 2mu4
4. `dR_corr_to_mu4_cross_terms/` — PbPb mu4 cross-term (requested as separate top-level)

All should have per-year and centrality-bin subdirectories. An `isBNL` toggle is needed in TrigEffPlotter. Exact revised structure pending resolution of Issues 1–3.

## Accumulated Findings

(append-only)

### F1: TrigEffPlotter uses Mac-only paths (2026-05-18)
TrigEffPlotterPbPb.cxx line 13: `data_dir = "/Users/yuhanguo/Documents/..."`. No BNL equivalent. SingleMuEffcyPtTurnOnFitter already uses BNL paths.

### F2: Pipeline 2 role-swap is correctly combined (2026-05-18)
`levels_trg_effcy_to_be_summed = {0,1}` sums over pair_sign AND mu4_selection (mu1passmu4 + mu2passmu4). Post-summed histograms are what's used for fitting. This is correct.

### F3: Pipeline 3 single-muon denominator has WRONG mu4 filter (2026-05-18)
Root cause identified: the code conditions the denominator on the "other" muon passing mu4 (`m1.passmu4` for mu1tag, `m2.passmu4` for mu2tag). This is the legacy mu4_mu4noL1/2mu4 approach where the 1st muon must fire mu4. For the mu4 single-muon dR correction, there should be **no trigger requirement** on the other muon — we are measuring whether the mere presence of another offline muon at dR affects triggering. Fix: remove `m1.passmu4`/`m2.passmu4` from the denominator filters. Additionally, the two role-swap terms are never combined (separate ratio graphs produced).

### F6: Misleading `_2mu4_` suffix in Pipeline 3 histogram names (2026-05-18)
Histogram suffixes like `_mu1tag_2mu4_denom` contain "2mu4" which refers to D5 (TF1 fits from 2mu4-conditional efficiency), but reads as if these are for 2mu4 trigger. Should be `_mu4_` since this is the mu4 dR correction.

### F7: No unfitted 2D fallback for gap regions in EvaluateSingleMuonEffcyPtFitted (2026-05-18)
Data.cxx line 541: returns -1.0 for gap-region muons (no TF1 fit). Per physics procedure, should fall back to unfitted 2D histogram bin content. Currently all gap-region muons are excluded (weight=0), biasing the dR correction.

### F4: PP Pipeline 3 not implemented yet (2026-05-18)
`FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` is empty in RDFBasedHistFillingPP.cxx (line 190). User requirement: for 2mu4, no tag/probe, no cross-term — one dR correction term only using pair-level `pass2mu4`.

### F5: Pipeline 3 centrality binning not implemented (2026-05-18)
`MakeAndWriteDRTrigEffGraphs()` calls `MakeAndWriteDRTrigEffGraphsHelper({})` — empty categories = centrality-inclusive only. `FillTrigEffcyHistsInvWeightedbySingleMuonEffcies()` also only fills centrality-inclusive histograms. User requires centrality-binned output.

### F8: Inverse-weighted dR correction ratio systematically ~1.15-1.28, not ~1.0 (2026-05-19)

**Observation:** After Step 8 fixes (new code, rerun), the inverse-weighted ratio
`Σ(1/ε × pass_mu4) / Σ(all valid)` is ~1.15 (centrality-inclusive OS) across ALL dR bins,
including high-dR (0.6-0.8) where no proximity effect is expected. The cross-term shows
similar behavior (~1.14).

**Asymmetry between roles:** Role1 (mu2 probe) ratio = 1.088; Role2 (mu1 probe) ratio = 1.283.
mu1 (leading pT) has raw mu4 pass rate ~84%, mu2 (subleading) ~54%.

**Data (PbPb 2023 OS, DR_zoomin):**
- Combined: Num=24838, Den=21459, Ratio=1.157
- Bin 1 (dR 0.00-0.04): ratio=6.15 (very few events, expected noise)
- Bin 17 (dR 0.64-0.68): ratio=1.16 (well above 1 at high dR)
- Full-range DR (0-5.75): overall ratio=1.186

**Mathematical expectation:** If ε^{nc}(pT,q·η) is the true P(fire mu4|pT,q·η), then
E[pass/ε] = E[1] = 1 for any dR. Ratio > 1 means fitted ε systematically underestimates
the true single-muon mu4 efficiency.

**TF1 plateau values:** ε ≈ 0.6-0.8 at high pT. These are P(2mu4 | mu4, dR > 0.8) from
Pipeline 2, which should equal P(probe fires mu4) for well-separated pairs. If the actual
mu4 efficiency is higher than the fitted value, 1/ε inflates the numerator.

**Hypotheses for ratio > 1:**
1. **Fitted ε too low:** The Fermi+log fit may systematically undershoot the data at high pT
   (where most muons are), especially if the fit is pulled down by noisy low-pT bins.
2. **Centrality-averaged efficiency mismatch:** Pipeline 2 fits are per-centrality, but the
   centrality-inclusive Pipeline 3 averages over all centralities. If events in different
   centrality bins have different weights, the effective average ε might not match. (RULED OUT:
   Pipeline 3 uses per-event centrality to pick the right ctr-bin TF1.)
3. **dR > 0.8 selection bias in Pipeline 2:** Pipeline 2 only uses pairs with dR > 0.8 for
   tag-and-probe. If high-dR pairs have systematically different kinematics (e.g., lower
   average pT for the probe), the fitted ε(pT) may not represent the true efficiency at
   each pT. (UNLIKELY: the fit is a function of pT, so it should capture pT dependence.)
4. **Gap-region 2D histogram fallback:** Muons in q·η gap regions fall through to the 2D
   histogram, which has noisy bin-by-bin efficiency. If the 2D values are systematically
   lower than the true efficiency, the inverse weight inflates the ratio.

**Key diagnostic from pair_pt_log plot (2026-05-19):**
The offset is **pT-dependent**: ~1.3 at pair pT ≈ 6-8 GeV, decreasing to ~1.0 at pair pT > 30 GeV.
This strongly suggests the fitted ε systematically underestimates the true efficiency at low pT
(near threshold), where the Fermi+log fit shape matters most. At high pT (well above threshold),
the fit plateau is closer to the true efficiency and the ratio approaches 1.

The full-range DR plot confirms the offset is **dR-independent** (flat from dR~0.5 to dR~5.75).
Deta_zoomin and Dphi_zoomin are also flat. This rules out any dR/angular structure — the issue
is purely a normalization offset driven by the pT-dependent mismatch.

**Sub-steps and results:**
- [x] H1: Compare fitted TF1 values against raw data efficiency in pT bins
  - INITIAL RESULT (flawed): macro projected the 2D *divided* hist onto pT, which sums
    efficiencies across q_eta bins → values >1, meaningless comparison.
  - CORRECTED: project numerator and denominator separately, divide → proper ε(pT) per q_eta bin.
  - **Result: 2205 bins compared, 48.9% undershoots, 51.1% overshoots, avg mismatch = +3.9%.**
    Per centrality: 0-5% = +4.9%, 5-10% = +4.5%, 10-20% = +1.8%, 20-30% = +3.0%,
    30-50% = +4.0%, 50-80% = +6.8%. All small, no systematic undershoot. **H1 RULED OUT.**
  - Macro: `plotting_codes/investigate_H1_fit_vs_data.C` (has the flawed approach; the corrected
    analysis was done interactively)
- [x] H2: Check normalization offset across centrality bins (pre-compaction)
  - Offset varies significantly: 0-5% → 1.265, 30-50% → 1.072 (single-muon);
    0-5% → 1.420, 30-50% → 1.014 (cross-term). NOT ruled out — centrality matters.
- [x] H3: dR > 0.8 selection bias in Pipeline 2
  - The full-range DR graph shows the ratio is flat at ~1.15-1.20 from dR=0.2 to dR=5.0.
    If dR>0.8 selection bias caused a pT-spectrum mismatch, we would expect a dR-dependent
    ratio (higher at small dR, flat at large dR). The flatness rules out a dR-specific bias.
    **H3 RULED OUT.**
- [x] H4: Gap-region 2D histogram fallback bias (pre-compaction)
  - ~7% of lookups fall in gap regions. Gap-region avg ε = 0.906 vs non-gap avg = 0.871.
    If anything, gap regions have *higher* ε, which would *reduce* the ratio. **H4 RULED OUT.**

**H5: mu6 support trigger mismatch — DISPROVEN (2026-05-19)**

Initially identified as root cause: PbPb 2023 uses `use_mu6_for_trg_eff=true` (line 65,
DimuonDataAlgCoreT.c), adding HLT_mu6 to `passmu4`. Pipeline 2's ε^{nc} is measured from
2mu4 (no mu6), creating a mismatch. However, PbPb 2024/2025 have `use_mu6_for_trg_eff=false`
AND unprescaled mu4, yet show the SAME ~1.2 offset (2024: 1.23, 2025: 1.24). Therefore mu6
is NOT the primary cause. The mu6 mismatch may contribute a small additional offset for 2023,
but is not the dominant effect.

**ROOT CAUSE IDENTIFIED — H7: Event-level trigger selection bias (2026-05-19)**

Events enter the ntuple only if `b_HLT_mu4 && (m1 fires mu4 || m2 fires mu4)` (lines 440-446,
DimuonDataAlgCoreT.c). Pipeline 3 role1 (probe = m2) has NO trigger requirement on m1 in the
denominator, but the event was selected because at least one muon fired. This creates a
selection bias:

```
P(m2 fires | event in sample) = p2 / (p1 + p2 - p1*p2)  >  p2
```

The inverse-weighted ratio becomes:
```
E[1(m2 fires)/ε2] / E[1] = 1/(p1 + p2 - p1*p2)
```

For typical values p1 ≈ 0.7, p2 ≈ 0.5 → ratio ≈ 1/(0.85) = 1.18, matching observations.

**Why this explains all observations:**
- **Year-independent:** The event-level OR selection exists for all years.
- **Centrality-dependent:** p1, p2 vary with centrality (lower efficiency in central →
  larger bias). 0-5%: 1.265 vs 30-50%: 1.072.
- **pT-dependent:** Low-pT muons have lower ε → P(OR) is lower → 1/P(OR) is larger.
  At high pT (>30 GeV), both muons fire efficiently → P(OR) ≈ 1 → offset ≈ 1.
- **dR-independent:** Trigger-level event selection, not spatial — flat from dR=0.2 to 5.0.

**H6: pass2mu4 vs (m1.passmu4 && m2.passmu4) mismatch — RULED OUT (2026-05-19)**

Tested on PbPb 2025 ntuples: `pass2mu4` and `m1.passmu4 && m2.passmu4` agree to within 0.07%
(0 pairs have pass2mu4 but not both passmu4; only 41 out of 56k pairs have both passmu4 but
not pass2mu4). Trigger matching is not the issue.

**Implication for the dR correction:**

Pipeline 2's ε^{nc} = P(pass_2mu4 | tag fires mu4, dR > 0.8) is measured under the SAME
event-level selection (events must have b_HLT_mu4 && at least one muon fires). The selection
bias affects both the efficiency measurement (Pipeline 2) and the inverse weighting
(Pipeline 3). Whether the offset cancels in the final dR correction ratio depends on whether
Pipeline 2's tag-and-probe conditioning already accounts for the OR selection. This requires
further analysis.

### F9: plot_dR_trig_corr.C uses wrong axis ranges for Deta_zoomin and Dphi_zoomin (2026-05-19)

The plotting code has `{0.0, 0.4}` for both Deta_zoomin and Dphi_zoomin axis ranges.
The histograms actually have range [-0.8, 0.8] (matching var1D_pbpb.json). The plot only
shows 1/4 of the data. Must fix to [-0.8, 0.8].

### F10: No full-range DR plot in dR correction plotting code (2026-05-19)

`plot_dR_trig_corr.C` only plots DR_zoomin (0-0.8) and DR_0_2 (0-0.2). The full-range DR
histogram (0-5.75, 40 bins) exists in the ROOT file but is never plotted. Need to add it to
see the dR correction across the full kinematic range.

## Ruled Out

(append-only)

### R1: Coarse vs fine q·η mismatch in TF1 lookup (2026-05-19)
Initially suspected: `EvaluateSingleMuonEffcyPtFitted` uses `q_eta_proj_ranges_fine_excl_gap`
while TF1s were produced with coarse binning. Checked: the TF1 names contain fine bin edges
(`0_10_TO_0_50`, etc.), confirming the fitter used fine binning regardless of
`useCoarseQEtaBin`. Ruled out as cause of normalization offset.

### R2: Missing centrality-inclusive TF1 (2026-05-19)
`FindCtrSuffix` returns per-event centrality bin suffix even for ctr-inclusive fills, so
per-ctr-bin TF1s are always used. This is correct — no centrality-inclusive TF1 needed.

### R3: H1 — Fitted ε too low (Fermi+log undershoot) (2026-05-19)
Corrected comparison (projecting num/den separately then dividing) shows 2205 bins with
48.9% undershoots, avg mismatch +3.9%. The Fermi+log fit describes the data well. Ruled out
as cause of the ~15-28% normalization offset.

### R4: H3 — dR > 0.8 selection bias in Pipeline 2 (2026-05-19)
Full-range DR graph shows ratio is flat (~1.15-1.20) from dR=0.2 to dR=5.0. A dR-dependent
selection bias would produce a dR-dependent ratio. Flatness rules this out.

### R5: H4 — Gap-region 2D histogram fallback bias (2026-05-19)
~7% of lookups in gap regions. Gap-region avg ε (0.906) is higher than non-gap (0.871), so
gap fallback would reduce the ratio, not increase it. Ruled out.

### R6: H5 — mu6 support trigger mismatch (2026-05-19)
Initially identified as root cause for PbPb 2023: `use_mu6_for_trg_eff=true` adds mu6 to
passmu4, but ε^{nc} is from 2mu4 (no mu6). Disproven: PbPb 2024/2025 have
`use_mu6_for_trg_eff=false` AND unprescaled mu4, yet show the same ~1.2 offset
(2024: 1.23, 2025: 1.24). mu6 may contribute marginally for 2023 but is not the primary cause.

### R7: H6 — pass2mu4 vs (m1.passmu4 && m2.passmu4) mismatch (2026-05-19)
Tested on PbPb 2025: 0.07% disagreement (41/56k pairs). Essentially identical. Ruled out.

## Remaining Work

- Issues 1-5: FIXED in Step 8 (code changes done, compilation verified)
- Directory restructure + isBNL toggle: DONE (Step 8.5)
- Step 9: Rerun Pipeline 3: DONE
- Step 10: Fix dR correction plotting issues: DONE
- **Issue 6: Normalization offset ~1.15-1.28 — ROOT CAUSE: event-level OR selection bias (F8, H7)**
  - Events require `b_HLT_mu4 && (m1 || m2 fires mu4)`. Pipeline 3 denominator has no trigger
    req on "other" muon, but event selection biases P(probe fires) upward by factor 1/(p1+p2-p1*p2).
  - Need to determine: does this offset cancel when Pipeline 2 ε^{nc} is also measured under the
    same event selection? If so, the offset is a bookkeeping artifact, not a physics bias.
- Issues 7-10: DONE (plot fixes committed)
- **Investigation macros to clean up:** `plotting_codes/investigate_H1_fit_vs_data.C` (has
  flawed approach; consider deleting or annotating)

## Latest Stage

**Investigation complete (2026-05-19).** Root cause: **event-level OR trigger selection bias (H7).**
Events require (m1 OR m2 fires mu4), but old separate single-muon term tests only one muon →
P(m2 fires | event) = p2/(p1+p2-p1*p2) > p2 → ~1.2 offset.

**Update (2026-05-19, post pair-level implementation):** The pair-level approach (D6)
reduces the offset (from ~1.24 to ~1.15-1.20) but does NOT eliminate it. The remaining
offset is a kinematic selection bias: E[1/ε_pair | pair passes] = 1/P(pair passes) > 1
because high-ε pairs are over-represented in the sample. This is inherent to measuring
a trigger's efficiency on a sample selected by that same trigger, regardless of whether
the measurement is per-muon or per-pair.

**Key insight from plots:** The pair-level dR ratio shows clear dR-dependent structure —
~1.08 at dR < 0.1, rising to ~1.2 at dR > 0.3. This means the pair trigger IS suppressed
at small dR (~10% relative to large dR). The pT-sliced plot confirms: the 8-12 GeV slice
shows the largest offset (~1.2) and strongest dR structure, while the 40-120 GeV slice is
nearly flat at ~1.07.

The full-range DR plot shows the ratio peaks at dR ~ 1.5 then decreases at very large dR
(>3), with a consistent shape across all years.

**Path forward:** The flat offset is a normalization artifact from the selection bias.
The physical dR correction can be extracted by normalizing to the large-dR plateau
(ε_dR = ratio / plateau). Alternatively, a supporting trigger (e.g., mu4_mu4noL1 or 2mu4)
could be used to select events for unbiased efficiency measurement — but that requires
the supporting trigger to work (currently broken due to passmu4noL1 skimming bug).

Implementation tracked in mu4_trig_effcy_implementation.md, Steps 12-14.
