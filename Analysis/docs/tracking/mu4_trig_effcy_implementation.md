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

Each term factorizes into a no-correlation baseline × dR-dependent correction:

```
P(pair passes mu4 | dR) = ε₁^{nc} · ε_dR^{single}(dR)
                        + ε₂^{nc} · ε_dR^{single}(dR)
                        − ε₁^{nc} · ε₂^{nc} · ε_dR^{cross}(dR)
```

where εᵢ^{nc} ≡ ε^{no-corr}(pTᵢ, q·ηᵢ) is the single-muon mu4 efficiency with no proximity effects.

This is the final result applied to the analysis. It requires three measured quantities: ε^{nc}(pT, q·η), ε_dR^{single}(dR), ε_dR^{cross}(dR).

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

#### 3a. Single-muon dR correction (ε_dR^{single})

**What we measure:** Does the mere **existence** of another offline muon at angular distance dR affect P(this muon fires mu4)?

**CRITICAL:** No trigger requirement on the other muon. This is different from the legacy mu4_mu4noL1/2mu4 approach where the other muon was required to pass mu4. For the mu4 single-muon correction, we correct for: P(2nd muon passes mu4 | 1st offline muon exists at dR) = ε^{nc}(pT₂, q·η₂) × ε_dR^{single}(dR).

**Per-event procedure** (for each pair, with "2nd muon" as probe):
1. Evaluate ε^{nc}(2nd muon) from Pipeline 2: use pT-fitted TF1 if q·η is in a fitted bin; use unfitted 2D histogram bin content if q·η is in a gap region.
2. **Denominator:** all pairs where ε(2nd muon) is defined. Weight = 1. **No trigger requirement on either muon.**
3. **Numerator:** from the denominator population, pairs where 2nd muon passes mu4. Weight = 1/ε^{nc}(2nd muon).
4. Role-swap: each muon takes the role of "2nd muon" in turn. **Always combine** (sum numerators, sum denominators) before dividing.
5. TH1::Divide(numerator, denominator) → ε_dR^{single}(dR). Also fill as function of pair pT, deta, dphi, and dR in pair-pT bins (since dR and pair pT are correlated).

If ε_dR^{single} = 1 at all dR, there is no proximity effect.

#### 3b. Cross-term dR correction (ε_dR^{cross})

**What we measure:** P(μ₁ fires mu4 ∧ μ₂ fires mu4) = ε₁^{nc} × ε₂^{nc} × ε_dR^{cross}(dR).

**Per-event procedure** (double inverse weighting):
1. Evaluate ε^{nc} for both muons from Pipeline 2 (same fallback logic as 3a).
2. **Denominator:** all pairs where both ε₁ and ε₂ are defined. Weight = 1.
3. **Numerator:** pairs where BOTH muons pass mu4. Weight = 1/(ε₁^{nc} × ε₂^{nc}).
4. No role-swap needed (symmetric by construction).
5. TH1::Divide → ε_dR^{cross}(dR).

#### 3c. PP 2mu4 dR correction (ε_dR^{2mu4})

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
- [x] 2a. Implement (done, has bugs) → **2a-FIX below**
- [x] 2b–d. Run (done, output has bugs) → **2d-RERUN below**
- [ ] **2a-FIX.** Fix physics bugs per investigation doc (§3a, §3b, D5 naming). **Reviewer: /review-analysis-code** with §3a, §3b, §3c in task prompt.
- [ ] **2d-RERUN.** Rerun Pipeline 3 for PbPb after fixes.
- [ ] 2e. Implement PP24 Pipeline 3 (§3c: no tag, no cross-term). **Reviewer: /review-analysis-code** with §3c in task prompt.
- [ ] 2f. (When PP24 ntuple ready) Run for PP24.

### Phase 3: Plotting
- [ ] 3a. Create SingleMuTrigDRCorrPlotter. **Reviewer: /review-plot.**
- [ ] 3b. /review-plot on all output.

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
**Physics (§Pipeline 2, §3a):** Muon 1 and muon 2 are physically interchangeable — the labeling is arbitrary. Swapping roles and combining (sum num, sum denom, then divide) doubles statistics and ensures the result doesn't depend on arbitrary labeling. Role-swap terms are **never** used separately.

### D4: Correlation plots include deta, dphi, pair pT (2026-05-18)
**Physics (§3a):** The factorization assumes ε_dR depends only on dR. Plotting vs deta, dphi, and pair pT tests this assumption — if the correction has deta/dphi structure, the factorization breaks down.

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
Status: TODO (tracked in investigation doc F8)
See `mu4_trig_effcy_investigation.md` F8 for hypotheses and sub-steps.

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

## Remaining Work

- Rerun Pipeline 3 after Step 8 fixes (Step 9) — may not be needed if current output already uses new code
- Investigate normalization offset (Step 11, tracked in investigation doc F8)
- Update docs: README.md, analysis docs
- PP24 Pipeline 2+3 when re-skim completes
- Reprocess ntuples for passmu4noL1 fix

## Latest Stage

**Step 10 COMPLETE.** dR correction plotting fixes done (full-range DR, axis ranges, Y range, logx). Key diagnostic: normalization offset is pT-dependent (~1.3 at low pair pT → 1.0 at high pT), dR-independent. Next: investigate normalization offset root cause (Step 11).

## Potential Issues

1. **Variance inflation at low pT:** When ε^{no-corr} is small (near threshold ~4 GeV), 1/ε weights blow up. Trimming undefined bins handles infinity, but near-zero ε still inflates variance. Monitor ε_dR noise at low pair pT.

2. **"No trigger on muon 1" vs legacy:** The single-muon dR correction (no trigger requirement on muon 1) measures a different effect than the legacy code's mu4-conditioned correction. ε_dR may be close to 1 if the effect is small.

3. **ε_dR^{cross} ≠ (ε_dR^{single})²:** The cross-term captures genuine two-body trigger correlation (L1 sector sharing), distinct from the single-muon nearby-track effect.

4. **TH1::Divide (not TEfficiency):** Correct for per-event weighted histograms. Uncertainties are Gaussian error propagation; for proper uncertainties, bootstrap would be needed.

5. **Factorization validity:** The assumed factorization ε(pT, q·η, dR) = ε^{no-corr}(pT, q·η) × ε_dR(dR) may break if dR correction depends on q·η. The pair-pT-binned and deta/dphi outputs will help validate.
