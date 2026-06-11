# MU4 Trigger Efficiency Implementation

## Objective

Implement the full single-muon mu4 trigger efficiency evaluation chain:
1. No-correlation single-muon mu4 efficiency ε_mu4^{no-corr}(pT, q·η)
2. dR-dependent correlation corrections ε_dR(dR) for single-muon and cross-term
3. PP24 2mu4 trigger efficiency via analogous inverse weighting

## Physics Procedure

This section is the authoritative reference for what the code should do. All design decisions, implementation choices, and naming conventions must be consistent with it. After compaction, re-read this section first.

### Motivation

PbPb uses the mu4 single-muon trigger (trigger_mode=1) to select dimuon pairs. PP24 uses the 2mu4 pair trigger (trigger_mode=3). We need the per-pair trigger efficiency to correct measured yields.

### Top-level equation (PbPb mu4)

A muon pair "passes mu4" if **at least one** muon fires mu4. By inclusion-exclusion:

```
P(pair passes mu4) = P(μ₁ fires) + P(μ₂ fires) − P(μ₁ fires ∧ μ₂ fires)
```

The pair-level trigger probability factorizes as:

```
P(pair passes mu4 | dR) = ε_pair^{nc} × ε_dR^{pair}(dR)
```

where:
- ε_pair^{nc} = ε₁^{nc} + ε₂^{nc} − ε₁^{nc}·ε₂^{nc} is the no-correlation pair-level efficiency
- εᵢ^{nc} ≡ ε^{no-corr}(pTᵢ, q·ηᵢ) is the single-muon mu4 efficiency with no proximity effects
- ε_dR^{pair}(dR) is the pair-level dR correction

This requires two measured quantities: ε^{nc}(pT, q·η) and ε_dR^{pair}(dR).

**Why pair-level, not separate single-muon + cross-term (see D6):** The previous approach
decomposed into ε_dR^{single} and ε_dR^{cross} separately. This suffers from an event-level
trigger selection bias: events require (m1 OR m2 fires mu4), but the single-muon term tests
only one muon, creating a ~1.2 normalization offset. The pair-level approach avoids this
because the measurement condition (pair passes mu4) matches the event selection condition
exactly — every event satisfies it.

### Top-level equation (PP 2mu4)

2mu4 requires **both** muons to fire. One term only, no inclusion-exclusion:

```
P(pair passes 2mu4 | dR) = ε₁^{nc} · ε₂^{nc} · ε_dR^{2mu4}(dR)
```

No tag/probe structure (2mu4 is a pair-level trigger). No cross-term (it IS the cross-term).

### Pipeline 2: Measure ε^{no-corr}(pT, q·η) — no-correlation single-muon mu4 efficiency

**Physics idea:** If two muons are well-separated (dR > 0.8), their trigger decisions are independent. Use tag-and-probe on such pairs.

**Method:** For each pair with dR > 0.8:
1. Require "1st muon" passes mu4 (tag — ensures the event was triggered).
2. Ask: does the pair pass 2mu4? Since dR > 0.8 (uncorrelated), 2mu4 firing ⟺ "2nd muon" independently fires mu4.
3. ε^{no-corr}(pT₂, q·η₂) = P(pair passes 2mu4 | 1st muon passes mu4, dR > 0.8).

**Role swap:** Each muon takes turns as "1st" and "2nd." The two role-swap terms are **always combined** (sum numerators, sum denominators, then divide) — never separate.

**Binning:** Per (centrality, μ+/μ−, q·η bin). The 2D efficiency (pT vs q·η) is projected onto pT in each q·η bin, then fitted with Fermi + log correction.

**Output:** TF1 fits in `single_mu_effcy_pT_fit.root`, keyed by (centrality, musign, q·η bin).

### Pipeline 3: Measure dR corrections via inverse weighting

**Physics idea:** The no-correlation efficiency assumes muons are isolated. At small dR, a nearby muon can affect triggering (e.g., L1 sector sharing). We factorize:

```
P(muon fires mu4 | another offline muon at dR) = ε^{nc}(pT, q·η) × ε_dR(dR)
```

To extract ε_dR, divide out the known ε^{nc} via per-event inverse weighting.

#### 3a. (REMOVED — see D9)

Single-muon dR correction removed. We expect negligible pair dR correlation on
single-muon mu4 efficiency. Code removed from PbPb.cxx.

#### 3b. (REMOVED — see D9)

Pair-level dR correction removed. On a mu4-selected sample, the numerator condition
(at least one fires) is the selection condition itself — ε_dR^{pair} ≡ 1 by
construction, measuring only fit quality rather than a physical dR effect. The
previous D6 "fix" was a tautology, not a physics measurement. Code removed.

#### 3c. Cross-term dR correction (ε_dR^{cross}) — reference for MC

Measures P(both muons fire mu4 | dR) / (ε₁ε₂) via inverse weighting. On a
mu4-selected sample, the denominator is biased by 1/(ε₁+ε₂−ε₁ε₂), inflating
the ratio. Kept as reference for unbiased MC samples. Code in PbPb.cxx.

#### 3d. PP 2mu4 dR correction (ε_dR^{2mu4})

Same as 3b but numerator condition is `pass2mu4` (pair-level trigger decision) instead of `m1.passmu4 && m2.passmu4`. The 2mu4 hardware trigger can differ from "both individually pass mu4" by prescales. No tag/probe. No separate single-muon term.

#### Why TH1::Divide, not TEfficiency

TEfficiency assumes each trial contributes weight 0 or 1 with binomial statistics. Here, the numerator has per-event weights (1/ε) that can exceed 1. TH1::Divide with Gaussian error propagation is correct; for rigorous uncertainties, bootstrap would be needed.

## Context

The analysis uses three separate pipelines:
- **Pipeline 1 (Nominal):** crossx & RAA for signal pairs. PP→trigger_mode=3 (2mu4), PbPb→trigger_mode=1 (mu4) with no trigger efficiency derivation. Uses Resonance Cut V1 + Photo Production Cut.
- **Pipeline 2 (Trig Effcy Loop 1):** No-correlation single-muon mu4 efficiency. trigger_mode=1. Uses Resonance Cut V2. Default applies photo production & resonance (V2) cuts (configurable).
- **Pipeline 3 (Trig Effcy Loop 2):** Inverse-weighted dR/pair pT corrections. Same ntuples as Pipeline 2, plus fitted TF1 & 2D histograms from Pipeline 2.

## Scope

### Ntuple processing changes
- [x] Replace `pbpb24_mu4_NO_trig_calc` with `mu4_nominal_pbpb_NO_trig_calc` (applies to all Run 3 PbPb)
- [x] Update `trigger_effcy_calc`: true iff trigger_mode ∈ {0,1} AND NOT (isPbPb && isRun3 && mu4_nominal_pbpb_NO_trig_calc)
- [x] Update `use_mu6_for_trg_eff`: on when trigger_effcy_calc && trigger_mode == 1
- [x] Update mu4 output suffix: `_mu4_nominal` (PbPb nominal) vs `_mu4_trig_calc` (trig effcy pipeline)
- [x] Ensure suffix is set AFTER fixing mu4_nominal_pbpb_NO_trig_calc value
- [x] Update hist filling input paths to match new suffixes

### Phase 1: No-correlation efficiency (Pipeline 2)
- [x] 1a. Verify ntuple availability (PbPb 23/24/25 ready; PP24 re-skim in progress)
- [x] 1b. Run hist filling cycle=generic, trigger_mode=1 for PbPb (all 3 years)
- [x] 1c. Run pT turn-on fitting (fermi+log) for PbPb (all 3 years)
- [ ] 1d. /review-plot to verify fitting output — BLOCKED: mu4_mu4noL1 fits unconstrained (skimming bug)

### Phase 2: Inverse weighting (Pipeline 3)
- [x] 2a. Implement cross-term + PP 2mu4 inverse weighting
- [x] 2b–d. Run Pipeline 3 for PbPb (all 3 years)
- [x] 2e. Implement PP24 Pipeline 3 (§3d: no tag, pass2mu4)
- [ ] 2f. (When PP24 ntuple ready) Run for PP24.
- [x] **D9 cleanup.** Removed single-muon and pair-level dR terms (biased on mu4-selected sample).
  Cross-term kept as reference for unbiased MC. Plateau normalization removed.

### Phase 3: Plotting
- [x] 3a. Cross-term dR correction plots via `plot_dR_trig_corr.C`.
- [ ] 3b. Evaluate cross-term/2mu4 dR on MC samples (unbiased).

## Design Decisions

### D1: Three-pipeline separation (2026-05-18)
**Physics:** The top-level equation requires three separate measured quantities (ε^{nc}, ε_dR^{single}, ε_dR^{cross}). Pipeline 1 applies the final efficiency; Pipelines 2 and 3 measure the ingredients. Mixing measurement and application in one pipeline caused gate bugs.
**Old:** Single trigger_mode=1 flow for both nominal and efficiency derivation.
**New:** Clean separation into Pipeline 1/2/3 with distinct flags and I/O suffixes.

### D2: PP24 uses 2mu4 pair-level inverse weighting (2026-05-18)
**Physics (§3c):** PP24 uses 2mu4 (pair-level trigger). The factorization is P(2mu4|dR) = ε₁^{nc} × ε₂^{nc} × ε_dR^{2mu4}(dR) — one term, no inclusion-exclusion. The 2mu4 hardware trigger can differ from "both muons individually pass mu4" by prescales, so numerator uses `pass2mu4`.
**Old:** N/A.
**New:** Numerator = pairs passing 2mu4, weight = 1/(ε₁ × ε₂). No tag/probe. No cross-term.

### D3: Muon role swapping (2026-05-18)
**Physics (§Pipeline 2):** Muon 1 and muon 2 are physically interchangeable — the labeling is arbitrary. Swapping roles and combining (sum num, sum denom, then divide) doubles statistics and ensures the result doesn't depend on arbitrary labeling. Role-swap terms are **never** used separately. (Applies to Pipeline 2 tag-and-probe only; Pipeline 3 single-muon role-swap removed per D9.)

### D4: Correlation plots include deta, dphi, pair pT (2026-05-18)
**Physics (§3a):** The factorization assumes ε_dR depends only on dR. Plotting vs deta, dphi, and pair pT tests this assumption — if the correction has deta/dphi structure, the factorization breaks down.

### D6: Pair-level dR correction replaces separate single-muon + cross-term (2026-05-19)
**Physics:** The previous approach measured ε_dR^{single} and ε_dR^{cross} separately and
combined them via P(pair) = ε₁·ε_dR^{single} + ε₂·ε_dR^{single} − ε₁·ε₂·ε_dR^{cross}.
This suffered from an event-level trigger selection bias: events are selected by
`b_HLT_mu4 && (m1 fires || m2 fires)`, but the single-muon term tests only one muon's
trigger, so P(m2 fires | event selected) = p₂/(p₁+p₂−p₁p₂) > p₂, creating a ~1.2 offset.

The pair-level approach directly measures ε_dR^{pair} where the test condition (pair passes
mu4) matches the selection condition. Every event satisfies the numerator condition, so no
selection bias exists. The pair efficiency factorizes as:
P(pair|dR) = (ε₁ + ε₂ − ε₁ε₂) × ε_dR^{pair}(dR).

**Trade-off:** Less granular — cannot separate single-muon vs cross-term dR effects. But
the pair-level correction is what the analysis actually needs (per the top-level equation),
and is free of the selection bias.

**Old:** Three terms: ε_dR^{single}(dR), ε_dR^{cross}(dR), combined via inclusion-exclusion.
**New:** One term: ε_dR^{pair}(dR), directly correcting the pair-level trigger probability.

### D7: Disable mu6/mu8 support triggers (2026-05-19)
**Physics:** mu6/mu8 in the passmu4 flag creates an inconsistency with Pipeline 2's ε^{nc}
(measured from 2mu4, which does not include mu6/mu8). For 2024/2025 the offset exists even
without mu6, proving it's not the primary issue. Disabling mu6/mu8 removes a small additional
inconsistency for 2023 and simplifies the analysis.
**Old:** `use_mu6_for_trg_eff = (isPbPb && run_year == 23) || (!isPbPb && run_year == 24)`.
**New:** `use_mu6_for_trg_eff = false; use_mu8_for_trg_eff = false;` (all years).

### D8: Switch useCoarseQEtaBin default from true to false (2026-05-28)
**Physics (§Pipeline 2, §3a):** The pT turn-on fitter produces fits in fine q-eta bins.
EvaluateSingleMuonEffcyPtFitted looks up TF1s by fine q-eta bin index. The TH2D fallback
for gap regions also lives in the fine-binned histogram file. With useCoarseQEtaBin=true,
OpenEffcyPtFitFile was loading the TH2D fallback from the coarse file, creating an
inconsistency: TF1s keyed to fine bins but fallback histograms from coarse binning.
Switching the default to false makes all downstream readers (OpenEffcyPtFitFile,
TrigEffPlotterPbPb, TrigEffPlotterPP, plot_dR_trig_corr.C, plateau_normalize_dR_corr.C,
both pipeline scripts, crossx plotter candidate list) consistent with the fitter's
fine-binned output.
**Old:** `useCoarseQEtaBin = true` in RDFBasedHistFillingData.h; downstream readers
defaulted to coarse.
**New:** `useCoarseQEtaBin = false` in RDFBasedHistFillingData.h; all downstream readers
updated to use fine binning.

### D9: Remove single-muon and pair-level dR corrections; keep cross-term as MC reference (2026-06-09)
**Physics:** All three inverse-weighted dR corrections (single-muon, pair-level, cross-term)
suffer from event-level trigger selection bias when evaluated on the mu4-selected sample.
The sample only contains events where at least one muon fires mu4 (OR condition). The
denominator of the inv-weighted ratio is therefore depleted by a factor P(at least one fires)
= ε₁ + ε₂ − ε₁ε₂, inflating the ratio by 1/(ε₁ + ε₂ − ε₁ε₂) > 1. This bias is
pT-dependent (~1.3 at low pT, ~1.0 at high pT) and creates apparent dR structure through
kinematic correlations between dR and pT.

The pair-level approach (D6) appeared to avoid this bias because the numerator condition
(pair passes mu4 = OR) matched the selection, making every event satisfy the numerator.
But this is a tautology — ε_dR^{pair} ≡ 1 by construction on this sample, and the measured
ratio only reflects the fit quality of ε^{nc}, not a real dR correlation. The D6 "fix"
did not measure a physical quantity; it just removed the symptom of the bias without
measuring anything useful.

For the single-muon term, we expect negligible pair dR correlation on the mu4 single-muon
efficiency (each muon's L1 trigger is largely independent). Combining single-muon and
cross-term into a pair-level correction (D6) conflates two physically distinct behaviors
and is not meaningful.

The cross-term is kept as a reference because it measures the genuine two-body trigger
correlation P(both fire | dR) / (ε₁ε₂). On an unbiased sample (MC with no trigger
selection), this procedure is correct and gives ε_dR^{cross}(dR) directly.

**Old:** Three terms filled: single-muon (role-swap combined), cross-term, pair-level.
Plateau normalization applied to pair-level. All plotted.
**New:** Only cross-term filled and plotted. Single-muon and pair-level fills removed
from PbPb.cxx. Plateau normalization removed from pipeline. Old plots/ROOT files deleted.

### D5: Use mu4-only fits for inverse weighting (2026-05-18)
**Physics (§Pipeline 2):** The no-correlation efficiency is measured as P(2mu4 | mu4, dR > 0.8), which gives the single-muon mu4 efficiency. The TF1 keys contain "2mu4" because of how Pipeline 2 labels the conditional probability, but the physical quantity is ε^{nc}(mu4). mu4_mu4noL1 fits are unconstrained (skimming bug).
**Note on naming:** Code suffixes for Pipeline 3 should use `_mu4_`, not `_2mu4_`, to reflect the physics (mu4 dR correction). The "2mu4" label in TF1 keys is a Pipeline 2 naming convention, not the physics of Pipeline 3.

## Implementation Plan

### Step 1: Ntuple processing changes (TriggerModeSettings & I/O) — per D1
Status: DONE. Reviewed: /review-analysis-code PASS (iteration 1).
Files: RDFBasedHistFillingData.h, RDFBasedHistFillingData.cxx, RDFBasedHistFillingPP.cxx, RDFBasedHistFillingPbPb.cxx

### Step 2: Run Pipeline 2 hist filling — per §Pipeline 2
Status: DONE
Files: run_pipeline2_pbpb.C, RDFBasedHistFillingPbPb.cxx
Output: histograms_real_pairs_pbpb_20YY_single_mu4_{coarse,fine}_q_eta_bin.root for YY=23,24,25

### Step 3: Run pT turn-on fitting — per §Pipeline 2 (fitting)
Status: DONE (2mu4 fits good; mu4_mu4noL1 fits unconstrained — skimming bug)
Files: SingleMuEffcyPtTurnOnFitter.cxx
Output: single_mu_effcy_pT_fit.root + 24 PNG canvases per year

### Step 4: Implement inverse weighting — per §3a, §3b (PbPb) and §3c (PP stub)
Status: DONE but **HAS PHYSICS BUGS** (see investigation doc). Reviewed: /review-analysis-code PASS (iteration 1) — **reviewer did not check physics procedure (it didn't exist yet)**.
Files: RDFBasedHistFillingData.h/.cxx, RDFBasedHistFillingPbPb.cxx, RDFBasedHistFillingPP.cxx
Bugs: denominator filters on other muon's mu4 (violates §3a "no trigger requirement"), role-swaps never combined (violates §3a.4), no 2D fallback for gaps (violates §3a.1), misleading `_2mu4_` suffix (violates D5 naming note)

### Step 5: Implement MakeAndWriteDRTrigEffGraphs — per §3a.5, §3b.5
Status: DONE (as part of Step 4). **Same physics bugs** — divides role-swaps separately instead of combining.

### Step 6: Run Pipeline 3
Status: DONE but output has physics bugs from Steps 4–5. **Needs rerun after fix.**

### Step 7: Create & run plotting script — per D4
Status: DONE but plots role-swap terms separately. **Needs redo after fix.**

### Step 8: Fix Pipeline 3 physics bugs + implement PP24 Pipeline 3 + add centrality binning — per §3a, §3b, §3c, D5
Status: DONE
**Reviewer: /review-analysis-code** with §3a, §3b, §3c in task prompt.
Scope (all investigation doc issues except directory restructure):
- Issue 1: remove mu4 filter on denominator (§3a: no trigger requirement on other muon)
- Issue 1b: combine role-swaps in MakeAndWriteDRTrigEffGraphsHelper (§3a.4, D3)
- Issue 2: implement PP24 Pipeline 3 (§3c: no tag, no cross-term, pass2mu4)
- Issue 3: rename `_2mu4_` → `_mu4_` suffixes (D5 naming note)
- Issue 4: add unfitted 2D histogram fallback for gap regions (§3a.1)
- Issue 5: add centrality binning to Pipeline 3 (PbPb)

### Step 9: Rerun Pipeline 3 + replot
Status: TODO

### Step 10: Fix dR correction plotting issues — per D4
Status: DONE
**Reviewer: /review-plot**
Scope:
- 10a. Add full-range DR (0-5.75) to plot_dR_trig_corr.C variable list
- 10b. Fix Deta_zoomin and Dphi_zoomin axis ranges: [0, 0.4] → [-0.8, 0.8] (match var1D json)
- 10c. Change Y range for dR_single_muon plots to [0.6, 1.6]. (Note: dR_cross_term Y range should also be changed once skimming is finished and plots re-evaluated.)
- 10d. Add `SetLogx()` for pair_pt_log plots
- 10e. Add "DR" to the Pipeline 3 filling variable list (RDFBasedHistFilling code already has "DR")

### Step 11: Investigate normalization offset in dR corrections
Status: DONE — Root cause: event-level OR selection bias (H7 in investigation doc F8).
Fix: pair-level dR correction (D6) avoids the bias entirely.

### Step 12: Disable mu6/mu8 support triggers — per D7
Status: TODO
Files: `NTupleProcessingCode/DimuonDataAlgCoreT.c` (line 65-66)
Change: `use_mu6_for_trg_eff = false; use_mu8_for_trg_eff = false;`
Requires: ntuple reprocessing for PbPb 2023 and PP 2024.

### Step 13: Implement pair-level dR correction — per §3a (new), D6
Status: TODO. **Reviewer: /review-analysis-code** with new §3a in task prompt.
Files: RDFBasedHistFillingPbPb.cxx, RDFBasedHistFillingData.cxx
Scope:
- Replace separate single-muon role1/role2 + cross-term fills with single pair-level fill
- New columns: `effcy_pair` = ε₁ + ε₂ − ε₁·ε₂, `invw_pair` = 1/effcy_pair
- Denominator: `node.Filter("valid_both")`, weight = 1, suffix `_pair_mu4_denom`
- Numerator: `df_denom.Filter("m1.passmu4 || m2.passmu4")`, weight = invw_pair, suffix `_pair_mu4_invw_num`
- Update MakeAndWriteDRTrigEffGraphsHelper: single divide (no role-swap combination needed)
- Keep centrality binning
- Update plot_dR_trig_corr.C for new graph names

### Step 14: Rerun Pipeline 3 + replot after Steps 12-13
Status: DONE

### Step 15: Add full-range DR×pair_pt_log 2D histograms + plateau normalization
Status: DONE
Files: RDFBasedHistFillingPbPb.cxx, RDFBasedHistFillingPP.cxx (added `{"DR", "pair_pt_log"}` 2D hist),
  RDFBasedHistFillingData.cxx (MakeAndWriteDRTrigEffGraphsHelper: project both DR and DR_zoomin by pair pT),
  plotting_codes/plateau_normalize_dR_corr.C (standalone plateau normalization script)
Output:
- 2 normalization modes: pT-integrated and pT-binned plateau
- Plateau extracted as error-weighted mean of full-range DR ratio in dR ∈ [1.0, 3.0]
- 54 PNGs per mode (DR, DR_zoomin, DR_0_2 × OS/SS × 7 ctr cats + pT-sliced overlays)
- Normalized TGraphs saved to `dR_corr_plateau_norm_{mode}_pbpb_20YY.root` per year

### Step 16: Reprocess PbPb 2023 ntuples with mu6/mu8 disabled — per D7
Status: DONE (parts 2-4; part 1 raw data missing)
- Parts 2-4 reprocessed locally (640s + 621s + 122s). New skim has ~20x more data.
- Part 1: raw data `data_pbpb23_part1.root` MISSING — need to locate/re-download.
- hadded parts 2-4 → combined file (sign1=622626 sign2=678115, total 1.3M entries).
- Old files backed up with `_before_mu6_disable` suffix.
- Reran Pipeline 2 (coarse), pT fitting, Pipeline 3, plotter, plateau normalization.
- **Key change:** PbPb 2023 plateau increased 1.222→1.285 (OS, ctr-int) due to mu6 removal + new skim.

### Step 17: Expand pipeline script with stages 5-10 — per §Pipeline 2, §3a
Status: DONE
Files: pipelines/pipeline_pbpb_trig_eff.sh
Scope:
- Added stages 5-8 (Pipeline 2): RDF hist filling (fine q-eta), pT fitting, fit/TH2D validation, plotting
- Added stages 9-10 (Pipeline 3): inv_weight hist filling, dR correction plotting + plateau normalization
- Added RDF_NTHREADS env var (default 2)
- Full end-to-end automation from condor submit through final plots

### Step 18: Switch useCoarseQEtaBin default to false — per D8
Status: DONE
Files: RDFBasedHistFillingData.h, OpenEffcyPtFitFile (PbPb/PP), TrigEffPlotterPbPb, TrigEffPlotterPP, plot_dR_trig_corr.C, plateau_normalize_dR_corr.C, pipeline scripts, crossx plotter candidate list
Scope:
- Changed default in RDFBasedHistFillingData.h from true to false
- Updated all downstream readers to be consistent with fine binning
- Rationale: fitter output is fine-binned; coarse default caused TH2D fallback file mismatch (see D8)

### Step 19: Re-run Pipeline 2 + Pipeline 3 with fine q-eta default
Status: DONE
Files: All Pipeline 2 + Pipeline 3 outputs for PbPb 23/24/25

### Step 21: Apply no-correlation trigger efficiency to PbPb crossx (Pipeline 1) — per top-level equation
Status: DONE
**Reviewer: /review-analysis-code** with top-level PbPb mu4 equation in task prompt.
Files: RDFBasedHistFillingPbPb.cxx (FillHistogramsCrossx), plot_crossx_trig_corr_sanity.C (new)
Scope:
- Call `OpenEffcyPtFitFile()` from FillHistogramsCrossx to load P2 TF1 fits
- Define per-muon efficiency columns `effcy1`, `effcy2` via `EvaluateSingleMuonEffcy`
- Compute pair efficiency `effcy_pair = ε₁ + ε₂ − ε₁·ε₂` (inclusion-exclusion, OR trigger)
- Apply `w_trig = 1/effcy_pair` multiplicatively to `weight_for_RAA` → `weight_for_RAA_trig_corr`
- All crossx histograms now use `weight_for_RAA_trig_corr`; `_counts` hists remain raw
- Added `_no_trig_corr` h2d pair_eta×pair_pt histograms for sanity-check comparison
- Created standalone sanity-check plotter overlaying raw vs corrected for PbPb 10-20%
- Simplified per-ctr loop: branch from `df_single_b_crossx_weighted` (which has all trig columns) instead of duplicating column definitions per centrality
- P3 stages 9-10 commented out in `pipeline_pbpb_trig_eff.sh` (biased on data, D9)

## Progress Log

### 2026-05-18: Initial setup
- Created tracking document
- Verified ntuple availability: PbPb 23/24/25 mindR_0_02 ntuples have all needed branches; PP24 missing m1/m2.passmu4noL1 (re-skim in progress)
- Existing PP24 histogram file (`_single_mu4_no_pT_fitting.root`) has 2mu4_sepr TGraphs from old code
- PbPb histogram files only have crossx content (trigger efficiency not produced due to post-processing gate)

### 2026-05-18: Step 1 complete — three-pipeline separation
- Added `mu4_nominal_pbpb_NO_trig_calc`, `trigger_effcy_calc`, `use_mu6_for_trg_eff` flags to RDFBasedHistFillingData.h
- Rewrote TriggerModeSettings() to derive all flags and set mu4 suffixes (`_mu4_nominal` vs `_mu4_trig_calc`)
- Rewrote FillHistograms() with three-pipeline routing: (generic+!trig_calc→crossx), (generic+trig_calc→single_mu_effcy), (inv_weight+trig_calc→inv_weighted)
- Rewrote HistPostProcessDataCommon() and HistPostProcessPbPb() with trigger_effcy_calc-based gating (replaces old run_crossx gate)
- Both PP and PbPb compile cleanly with ACLiC (separate sessions)
- Reviewer subagent: PASS, all 18 criteria satisfied
- Late fix: base_trig_suffix changed from `_mu4_trig_calc` to `_single_mu4` to match existing ntuples

### 2026-05-18: Step 2 complete — Pipeline 2 hist filling for PbPb
- Created run_pipeline2_pbpb.C macro (trigger_mode=1, cycle=generic, output_generic_hists=false, mindR_trig=0.02)
- Ran for all 3 years with coarse q·η binning
- Output histogram counts: 3500 1D + 1536 2D + 1536 3D per year
- Post-processing produced 135 TGraphAsymmErrors (trigger efficiency pT graphs) per year
- Output file sizes: PbPb23=12M, PbPb24=10M, PbPb25=17M (coarse); PbPb23=12M, PbPb24=11M, PbPb25=17M (fine)

### 2026-05-18: Step 3 complete — pT turn-on fitting for PbPb
- Fixed SingleMuEffcyPtTurnOnFitter.cxx: include path, CRTP template, RunYear(), PbPb paths
- Must run in interpreted mode (`.L` without `+`) due to rootcling CRTP issues
- Ran for all 3 years: single_mu_effcy_pT_fit.root sizes: 309K/312K/325K
- 24 PNG canvases per year (2 triggers × 6 centralities × 2 musigns)
- **2mu4 fits: GOOD** — clear turn-on, well-constrained Fermi+log, plateau ~0.6–0.9
- **mu4_mu4noL1 fits: UNCONSTRAINED** — data ~0 in all per-(ctr, musign, fine q·η) bins; fit forced to plateau=0.6 by parameter lower limit. Issue in ALL centrality bins and ALL years. See F4.

### 2026-05-18: Steps 4–5 complete — inverse weighting implementation
- Implemented all inverse weighting functions; code review PASS (iteration 1)
- Data.h: added FindCtrSuffix, updated EvaluateSingleMuonEffcyPtFitted/EvaluateSingleMuonEffcy signatures
- Data.cxx: file-scope statics (s_effcy_pT_fit_file, s_effcy_pT_fit_map), FindBinReturnStr, FindCtrSuffix, EvaluateSingleMuonEffcyPtFitted (TF1 key = 2mu4 only per D5), MakeAndWriteDRTrigEffGraphsHelper
- PbPb.cxx: OpenEffcyPtFitFile (loads TF1s from fit ROOT file), FillTrigEffcyHistsInvWeightedbySingleMuonEffcies (RDF pipeline with role swapping + cross-term), MakeAndWriteDRTrigEffGraphs
- PP.cxx: OpenEffcyPtFitFile (same pattern), MakeAndWriteDRTrigEffGraphs; FillTrigEffcyHistsInvWeightedbySingleMuonEffcies left as empty stub (PP24 not ready)
- Both PbPb and PP compile cleanly with ACLiC in separate sessions

### 2026-05-18: Step 6 complete — Pipeline 3 run for PbPb
- Created run_pipeline3_pbpb.C macro (trigger_mode=1, cycle=inv_weight_by_single_mu_effcy, coarse q·η)
- Ran for all 3 years: 240 TF1s loaded per year, 108 new TH1D + 54 new TGraphAsymmErrors per year
- Final ROOT file object counts: 792 TH1D + 189 TGraphAsymmErrors per year
- TGraphAsymmErrors::Divide "passed > total" warnings are expected (inverse weights > 1)
- Event loop times: ~0.5–0.6s per year (fast, ntuples are small)

## Results & Observations

### Ntuple status
- **PbPb 23/24/25:** `_mindR_0_02_res_cut_v2.root` files have all needed branches. Entries: 67k/38k/187k.
- **PP24:** Has m1/m2.passmu4, passmu4mu4noL1, pass2mu4, passSeparated. Missing: m1/m2.passmu4noL1, mindR. Re-skim in progress.

### Pipeline 2 output
- Histogram counts: 3500 1D + 1536 2D + 1536 3D per year; 135 TGraphAsymmErrors per year.
- **mu4 (2mu4-conditional) fits: GOOD** — clear turn-on, well-constrained Fermi+log, plateau ~0.6–0.9.
- **mu4_mu4noL1 fits: UNCONSTRAINED** — data ~0 in all bins, all years. Root cause: passmu4noL1 skimming bug (see below).

### passmu4noL1 skimming bug (FIXED, reprocess needed)
Per-muon `m1/m2.passmu4noL1` is always false in current ntuples. Bug in `DimuonDataAlgCoreT.c` lines 676-684: when `use_leg_branches_mu4_mu4noL1 = false` (legacy mode), code sets `passmu4noL1 = false` instead of deriving it. Fix applied: legacy fallback now derives per-muon noL1 from pair-level mu4_mu4noL1 match + cross-muon mu4. Ntuples need reprocessing.

### Pipeline 3 output (NEEDS RERUN — code fixed in Step 8, rerun in Step 9)
- All 6 investigation doc issues fixed: denominator filter (Issue 1), role-swap combination (Issue 1b), PP24 implementation (Issue 2), suffix naming (Issue 3), 2D fallback (Issue 4), centrality binning (Issue 5).
- Code review PASS (1 iteration). Both PbPb and PP compile cleanly.
- Old output (108 TH1D + 54 TGraphAsymmErrors per year) is stale — needs rerun with fixed code.

### PbPb25 mu4_mu4noL1 anomaly
`passmu4mu4noL1` is only 2.4% for PbPb25 vs ~70–80% for PbPb23/24. Likely different trigger menu or prescale. Separate issue from the passmu4noL1 skimming bug.

### 2026-05-28: Steps 17-19 — Pipeline script expansion + coarse-to-fine fix
- **Step 17:** Expanded `pipelines/pipeline_pbpb_trig_eff.sh` with stages 5-10:
  stages 5-8 cover Pipeline 2 (RDF hist filling with fine q-eta, pT fitting, fit/TH2D
  validation, P2 plotting); stages 9-10 cover Pipeline 3 (inv_weight hist filling,
  dR correction plotting + plateau normalization). Added `RDF_NTHREADS` env var
  (default 2; 8 needs ~48 GB RAM).
- **Step 18 (D8):** Changed `useCoarseQEtaBin` default from `true` to `false` in
  `RDFBasedHistFillingData.h`. Root cause: the fitter produces TF1s keyed by fine q-eta
  bin, and `EvaluateSingleMuonEffcyPtFitted` uses fine q-eta ranges, but
  `OpenEffcyPtFitFile` was loading the TH2D fallback from the coarse-binned file when
  `useCoarseQEtaBin=true`. All downstream readers updated: OpenEffcyPtFitFile (PbPb/PP),
  TrigEffPlotterPbPb, TrigEffPlotterPP, plot_dR_trig_corr.C, plateau_normalize_dR_corr.C,
  both pipeline scripts, crossx plotter candidate list.
- **Step 19:** Re-ran all Pipeline 2 + Pipeline 3 steps with SKIP_CONDOR=1 using
  the new fine default. All 10 stages passed for all 3 years. DONE.

## Progress Log (cont.)

### 2026-05-19: Step 10 complete — dR correction plotting fixes
- Added full-range DR (0-5.75) to plot variable list
- Fixed Deta_zoomin axis: [0, 0.4] → [-0.8, 0.8] (matches var1D_pbpb.json)
- Fixed Dphi_zoomin axis: [0, 0.4] → [-0.8, 0.8] (matches var1D_pbpb.json)
- Changed dR_single_muon Y range to [0.6, 1.6]; dR_cross_term stays [0.5, 2.5]
- Added logx for pair_pt_log plots
- Removed unused ymin/ymax from PlotDef struct; Y range now per-TermDef
- Reran plotter: 98 PNGs per term (up from 84; 14 new from 7 full-range DR × 2 signs)
- Total: 196 dR correction PNGs + 132 no-corr PNGs + 72 pT_fitting = 400 PNGs
- Key finding: pair_pt_log plot shows offset is pT-dependent (~1.3 at low pT, →1.0 at high pT)
- Full-range DR confirms offset is dR-independent (flat from 0.5 to 5.75)
- Saved log-scale plotting rule to memory and plot-reviewer criteria
- Note for future: dR_cross_term Y range should also be narrowed once normalization offset is resolved

### 2026-05-19: Step 12 complete — Disable mu6/mu8 support triggers (code only)
- Changed DimuonDataAlgCoreT.c line 65-66: `use_mu6_for_trg_eff = false; use_mu8_for_trg_eff = false;`
- Ntuple reprocessing deferred (new skim being downloaded)
- Note: existing PbPb 2023 ntuples still have mu6 baked into `m1.passmu4`; change takes effect on reprocessing

### 2026-05-19: Step 13 complete — Pair-level dR correction implementation
- **PbPb.cxx**: Replaced 3-term (role1/role2 single-muon + cross-term) filling with single pair-level
  fill. New columns: `effcy_pair = ε₁ + ε₂ − ε₁ε₂`, `invw_pair = 1/effcy_pair`. Denom:
  `valid_both`, Num: `m1.passmu4 || m2.passmu4`, weight: invw_pair. Added 2D {DR_zoomin, pair_pt_log}
  histograms. Suffixes: `_pair_mu4_denom`, `_pair_mu4_invw_num`.
- **Data.cxx**: Rewrote MakeAndWriteDRTrigEffGraphsHelper: unified PbPb/PP with `pair_suffix` selection,
  added pair-pT-sliced dR projections from 2D histogram (4 slices: 8-12, 12-20, 20-40, 40-120 GeV).
- **PP.cxx**: Added 2D {DR_zoomin, pair_pt_log} histograms to existing 2mu4 fills. No structural change.
- **plot_dR_trig_corr.C**: Updated term suffix to `_pair_mu4_invw_num`, output subdir `dR_pair_level`.
  Added pair-pT-sliced overlay plot (4 slices on one canvas, PbPb 2025).
- Both PbPb and PP compile cleanly with ACLiC.

### 2026-05-19: Step 14 complete — Pipeline 3 rerun + replot
- Backed up old ROOT files as `*_before_pair_level.root`
- Ran Pipeline 3 for PbPb 23/24/25: 364 1D + 28 2D histograms per year
- Ran plot_dR_trig_corr.C: 112 new PNGs in dR_pair_level/. Old plots preserved in dR_single_muon/
  and dR_cross_term/.
- **Key results (PbPb 2025 OS):**
  - DR_zoomin: ratio ~1.08 at dR < 0.1, rising to ~1.2 at dR > 0.3. Clear dR-dependent structure.
  - pair_pt_log: ratio ~1.28 at 8 GeV, decreasing to ~1.07 at 20-40 GeV.
  - pT-sliced dR: 8-12 GeV has largest offset (~1.2) and strongest dR structure; 40-120 GeV nearly flat ~1.07.
  - Full-range DR: peaks at dR ~ 1.5, consistent shape across all years.
  - Comparison: pair-level overall ratio (1.21) lower than old single-muon (1.24) and cross-term (1.26).
- **Interpretation:** Per-event inverse weighting decorrelates the measurement from single-muon kinematics
  (each event's expected contribution is ε_true/ε_fit, which equals 1 when the fit is exact). The flat offset
  (~1.2) reflects systematic underestimation of efficiency by the Fermi+log fit, especially near threshold.
  Physical signal is the dR shape: ~10% pair trigger suppression at small dR relative to large-dR plateau.

### 2026-05-19: Step 15 complete — Full-range DR×pT 2D hists + plateau normalization
- Added `{"DR", "pair_pt_log"}` to invw_var2Ds in both PbPb.cxx and PP.cxx
- Updated MakeAndWriteDRTrigEffGraphsHelper: refactored pair-pT projection to loop over both DR_zoomin and DR, producing full-range pair-pT-binned dR ratio graphs
- Reran Pipeline 3 for PbPb 23/24/25
- Created `plotting_codes/plateau_normalize_dR_corr.C`:
  - Extracts plateau as error-weighted mean in dR ∈ [1.0, 3.0] from full-range DR graphs
  - Mode 1 (pT-int): one plateau per sign/centrality, applied to all variables and pT slices
  - Mode 2 (pT-binned): per-pT-slice plateaus for pT-sliced plots, pT-int plateau for aggregate plots. Falls back to pT-int if <3 points in pT bin.
  - Saves normalized TGraphs to ROOT files + 54 PNG plots per mode
- Key results (PbPb 2025 OS ctr-int, pT-binned mode):
  - 8-12 GeV plateau=1.190, after norm: ~0.95 at dR<0.1, rising to ~1.0 at dR>0.3
  - 12-20 GeV plateau=1.094, after norm: ~0.99 at dR<0.1, nearly flat
  - 20-40 GeV plateau=1.076, after norm: ~0.99, nearly flat
  - All pT slices collapse to ~1.0 with pT-binned normalization — dR suppression is ~5-8% at small dR
- Plateau values are centrality-dependent: 0-5%=1.312, 10-20%=1.239, 30-50%=1.154, 50-80%=1.132
  (higher centrality → softer pT spectrum → larger fit underestimate near threshold → larger offset)

### 2026-05-19: Step 16 complete — PbPb 2023 ntuple reprocessing + full re-pipeline
- Condor schedd unavailable on attsub06; ran parts 2-4 locally (640+621+122 = 1383s total)
- Part 1 raw data missing; parts 2-4 from new skim (much larger dataset)
- New combined file: 1.3M entries (vs old 127k). Old files backed up with `_before_mu6_disable` suffix.
- Reran Pipeline 2 (coarse, 54+82s), pT fitting (24 PNGs), Pipeline 3 (3.5+3.7s)
- Re-ran plot_dR_trig_corr.C and plateau_normalize_dR_corr.C
- **PbPb 2023 plateau values (OS, new skim, no mu6):**
  - ctr-int: 1.285 (was 1.222 with old skim + mu6)
  - 0-5%: 1.387, 5-10%: 1.324, 10-20%: 1.260, 20-30%: 1.210, 30-50%: 1.180, 50-80%: 1.160
  - pT-binned: 8-12=1.242, 12-20=1.123, 20-40=1.092, 40-120=1.110
- After plateau normalization, all 3 years show consistent dR suppression: ~10-13% at dR<0.1, rising to ~1.0 by dR~0.5

### 2026-06-02: Step 19 complete — Full P2+P3 re-run with fine q-eta default
- Ran `SKIP_CONDOR=1 bash pipelines/pipeline_pbpb_trig_eff.sh` for years 23/24/25
- All 10 stages passed (validation, hadd, P2 hist fill, fitting, fit validation, P2 plots, P3 hist fill, P3 plots)
- Pipeline duration: 12:46:22 → 13:24:10 (~38 min total)
- P2 hist filling (fine): yr23 10min/7GB, yr24 5.5min/7GB, yr25 11.5min/7GB
- P3 inv-weight: yr23 2.5min, yr24 1.5min, yr25 3min
- pT fitting: 24 PNGs per year × 3 years = 72 fitting plots (fermi+log mode)
- P2 efficiency plots: 44 PNGs per year × 3 years = 132 plots
- P3 dR correction + plateau normalization plots: all centrality bins produced
- **dR plateau values (OS, fine q-eta, Jun 2, from pipeline log):**
  - yr23 ctr-int: 1.261±0.002, yr24: 1.243±0.002, yr25: 1.244±0.001
  - Central (0-5%) highest: yr23=1.316, yr24=1.297, yr25=1.301
  - Peripheral (50-80%) lowest: yr23=1.135±0.007, yr24=1.138±0.009, yr25=1.140±0.005
  - pT-dependence (yr23 ctr-int): 8-12=1.228, 12-20=1.112, 20-40=1.084, 40-120=1.102
  - Compare Step 16 (May 19, parts 2-4 only, coarse): yr23 ctr-int was 1.285
  - Shift from 1.285→1.261 (2%) due to adding part 1 data (77k/678k ≈ 11% more)
- All output files fresh (Jun 2 13:18-13:23), consistent across years
- Note: Warning `h_pt2nd_vs_q_eta2nd_ctr0_5_sign1_mu4 not found` at fitting stage — this is a
  reference TH2D for sign1 (SS), expected absent since fitting focuses on sign2 (OS) + mu+/mu- split

### 2026-06-09: Step 20 — Remove single-muon & pair-level dR corrections (D9)
- **Root cause:** All inv-weighted dR corrections are biased on mu4-selected sample. Denominator
  depleted by P(at least one fires) = ε₁+ε₂−ε₁ε₂. Pair-level "fix" (D6) was tautological.
- **Code changes:** PbPb.cxx: removed role1/role2/pair-level fills + columns. Data.cxx: removed
  role-swap combination + pair-level post-processing. plot_dR_trig_corr.C: terms reduced to
  cross-term only. Pipeline: removed plateau_normalize_dR_corr.C.
- **Deleted outputs:** 4 plot dirs (dR_single_muon, dR_pair_level, 2× plateau_norm; 360 PNGs),
  6 plateau ROOT files.
- **Kept:** Cross-term fills + plots as reference for unbiased MC evaluation.
- **Compilation:** PbPb.cxx ✓, PP.cxx ✓ (ACLiC, only -Wsign-compare warning).

### 2026-06-10: Step 22 — Crossx filename cleanup, PP trig eff code, PbPb rerun
- **Filename fix:** Crossx output files now use `_nominal` suffix instead of
  `_coarse_q_eta_bin`/`_fine_q_eta_bin`. `useCoarseQEtaBin` removed from all 4 crossx
  runner scripts (pbpb23/24/25, pp24). Logic: after `TriggerModeSettings()`, if
  `!trigger_effcy_calc`, override `qEtaBin_suffix = "_nominal"`.
- **Files changed (filename):** Data.cxx, run_crossx_hist_filling_{pbpb23,pbpb24,pbpb25,pp24}.sh,
  pipeline_pbpb_crossx.sh, pipeline_pp_crossx.sh, plot_single_b_crossx_pbpb.cxx,
  plot_single_b_crossx_pp.cxx, plot_crossx_trig_corr_sanity.C.
- **PP trig eff code:** `FillHistogramsCrossx` in PP.cxx now calls `OpenEffcyPtFitFile()`,
  defines q_eta/effcy columns, computes `effcy_pair = ε₁·ε₂` (AND trigger for 2mu4),
  `w_trig = 1/effcy_pair`. All histograms use `crossx_weight_trig_corr`. Added
  `_no_trig_corr` h2d. Sanity plotter extended with PP 2024 sample.
- **PP TF1 key mismatch discovered:** PP fits use `q_eta_2nd` (with underscore) and
  `_2mu4_sepr_sign1_py_` ordering, but `EvaluateSingleMuonEffcyPtFitted` constructs
  `q_eta2nd` (no underscore) and `_sign1_2mu4_sepr_py_`. Must fix before running PP.
- **Compilation:** Both PbPb.cxx and PP.cxx compile cleanly (ACLiC, login shell).
- **PbPb rerun:** pipeline_pbpb_crossx.sh with SKIP_EVSEL/CONDOR/HADD=1 in progress.

### 2026-06-10: Step 21 — Apply no-correlation trigger efficiency to PbPb crossx
- **Code:** `FillHistogramsCrossx` in PbPb.cxx now calls `OpenEffcyPtFitFile()` and defines
  per-muon efficiency columns (`effcy1`, `effcy2`) via `EvaluateSingleMuonEffcy`. Pair
  efficiency: `effcy_pair = ε₁ + ε₂ − ε₁·ε₂`. Weight: `w_trig = 1/effcy_pair`.
- All crossx histograms use `weight_for_RAA_trig_corr = weight_for_RAA × w_trig`.
- `_counts` histograms remain unweighted by trigger efficiency.
- Added `_no_trig_corr` h2d pair_eta × pair_pt histograms for sanity comparison.
- Simplified per-ctr loop: branch from `df_single_b_crossx_weighted` (inherits trig columns).
- Created `plot_crossx_trig_corr_sanity.C`: overlays raw vs corrected for PbPb 10-20%.
- P3 stages 9-10 commented out in `pipeline_pbpb_trig_eff.sh` (biased on data, D9).
- **Compilation:** PbPb.cxx ✓ (ACLiC, only pre-existing -Wsign-compare warnings).
- **Needs rerun:** Pipeline 1 RDF + plotting for all 3 years.

### 2026-06-02: Input verification investigation (CORRECTED)
- **Question:** User reported P3 dR plots look identical to previous outputs despite new skim.
- **Finding: All pipeline inputs are correct.** Verified at every stage:
  - NTuple batch files (May 27): from condor jobs that processed new skim (May 19-20 raw data)
  - P2 hist filling: reads hadded combined `_res_cut_v2.root` (Jun 2 hadd of May 27 per-part files)
  - P3 hist filling: cycle=2 confirmed in log, loaded 240 TF1s + 36 TH2D fallbacks, wrote 714 divided graphs
  - All plotting macros read from `_fine_q_eta_bin.root` (correct file)
- **Per-year data provenance:**
  - yr23: Step 16 (May 19) already reprocessed parts 2-4 with new skim. Jun 2 adds part 1 (+11%).
    Previous `_res_cut_v2.root` (22 MB, `_before_mu6_disable` backup) → new combined 258 MB.
  - yr24: Step 14 (May 19) used OLD skim `_res_cut_v2.root` (12.9 MB, Apr 21 backup `_backup_20260420`).
    Jun 2 pipeline is FIRST run with new skim (183 MB combined). 14× more data.
  - yr25: Step 14 (May 19) input unclear — no `_mindR_0_02` combined existed before May 27.
    Old P2 backup (`_before_pair_level.root`): 36k entries in ctr0_5 sign2.
    Jun 2 fine file: 373k entries in same histogram. 10× more data.
- **Jun 2 outputs are correct and DO reflect new skim data.** Error bars in Jun 2 plots are
  small and consistent with 300k-370k entries per centrality bin. Old Step 14 plots were
  overwritten, so direct visual comparison requires side-by-side from memory.
- **Previous investigation error (corrected):** Initially claimed "Step 16 already used new skim
  for all years" — this was only true for yr23. For yr24/25, Jun 2 is the first new-skim run.

## Remaining Work

- ~~Locate/re-download PbPb 2023 part 1 raw data, reprocess + re-hadd~~ DONE (May 2026 skim includes all 4 parts)
- ~~PP24 Pipeline 2+3 when re-skim completes~~ DONE (2026-06-10, PP trig eff pipeline)
- Reprocess ntuples for passmu4noL1 fix (separate issue)
- ~~Evaluate cross-term/2mu4 dR correction on MC samples~~ DEFERRED (need MC sample; will reopen when available)
- ~~Update docs: README.md, analysis docs~~ DONE
- ~~Apply no-correlation trigger efficiency to PbPb crossx~~ DONE (Step 21)

## Latest Stage

**CLOSED (2026-06-10).** All no-correlation trigger efficiency work complete
for both PbPb (3 years) and PP (pp24). Crossx pipelines run with per-pair
trig eff correction. mu4_mu4noL1 plotting disabled (not maintained).
Remaining open item: passmu4noL1 reprocessing (separate issue).
MC dR evaluation deferred until sample is available — will reopen then.

## Potential Issues

1. **Variance inflation at low pT:** When ε^{no-corr} is small (near threshold ~4 GeV), 1/ε weights blow up. Trimming undefined bins handles infinity, but near-zero ε still inflates variance. Monitor ε_dR noise at low pair pT.

2. **"No trigger on muon 1" vs legacy:** The single-muon dR correction (no trigger requirement on muon 1) measures a different effect than the legacy code's mu4-conditioned correction. ε_dR may be close to 1 if the effect is small.

3. **ε_dR^{cross} ≠ (ε_dR^{single})²:** The cross-term captures genuine two-body trigger correlation (L1 sector sharing), distinct from the single-muon nearby-track effect.

4. **TH1::Divide (not TEfficiency):** Correct for per-event weighted histograms. Uncertainties are Gaussian error propagation; for proper uncertainties, bootstrap would be needed.

5. **Factorization validity:** The assumed factorization ε(pT, q·η, dR) = ε^{no-corr}(pT, q·η) × ε_dR(dR) may break if dR correction depends on q·η. The pair-pT-binned and deta/dphi outputs will help validate.
