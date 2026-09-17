# pp24 trigger efficiency: HYBRID application (ΔR procedure below 74 GeV, single-value pair efficiency × single-muon SF above)

**Mode:** Implementation. **Created:** 2026-09-17. **Session:** "MC trigger efficiency apply".
**Parents:** `mc_trigger_efficiency.md` (ACTIVE — Step 3, the ΔR fits and their screens),
`mc_trigeff_single_value_pair_eff.md` (ACTIVE — the single-value pair efficiency this doc wires
in; its Remaining Work 2 is settled here), `mc_trigeff_dr_binning_approaches.md` (ACTIVE — the
|η^pair| fold this doc adopts for the ΔR cells), `mu_pt45_gap125_pairpt9_adoption.md` (ACTIVE —
the 9 GeV / 4.5 GeV selection every input here is on), `pp24_crossx_rerun_2026_08.md` (CLOSED —
the previous application and its Physics Procedure §2).

---

## Objective

Replace the pp24 per-pair 2mu4 trigger weight — until now the product of the DATA single-muon mu4
turn-ons and an MC ΔR correction in every pair-pT bin, with an expo → poly → raw-bin fallback —
by a **hybrid**: the ΔR procedure on **sign-separated, |η^pair|-folded, no-a-priori-plateau** fit
cells in the first six of the 8 coarse pair-pT bins, and the **single-value MC pair efficiency
times the product of the two single-muon data/MC scale factors** in the last two bins. Propagate
to the pp24 cross-section and regenerate the correction-stage comparison plots (no correction /
after trigger correction / after trigger + reconstruction correction).

## Autonomy Contract (DONE 2026-09-17 — every Done item met; both reviews PASS; kept for the record)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `ParamsSet` is the single source of the signal mass window (`signal_minv_min/max` +
     `SignalMinvCutExpr`), read by every live signal-selection site (data crossx, MC analogs,
     `MCTrigEffPairSel`, `PairTrigEff::Windows()`), bit-identical selection strings.
  2. A crossx-side evaluator delivering the §2 hybrid weight, with the four user decisions
     (D1–D4) implemented and MARKED TEMPORARY where they are; the old raw-bin tier unreachable;
     load-time census + throws as specified in §3.
  3. `RDFBasedHistFillingPP` applies it in the ONE weight-column helper (crossx OS, SS, template
     pass, generic) and books the trigger-only correction stage `_corr_unfolded_trig`.
  4. pp24 crossx pipeline (Stages 5–8, no Condor) rerun: `plots/single_b_analysis/pp24/`,
     `plots/sanity_check_crossx/PP_2024_*` (the stage plot now: uncorrected → +trigger →
     +trigger+reco), `plots/mc_data_compr/`.
  5. `/review-analysis-code` on the C++, `/review-plot` on the regenerated figures; docs
     (this doc, INDEX, `analysis_overview.md`, `signal_selection_change_impact.md`, parent docs'
     hand-off lines) updated; committed by explicit path.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The factorized weight `ε(pT1,q·η1) ε(pT2,q·η2) ε_ΔR(ΔR; cell)` is the right object where the
ΔR-shape fit has statistics: it is evaluated per pair and carries the within-cell kinematic
dependence. In the top two coarse pair-pT bins the pp24 fullsim sample cannot constrain that
shape (`mc_trigger_efficiency.md` R26/R32/R35: rejected fits, a runaway baseline, a raw-bin
placeholder at ε_ΔR = 2.16), while the same events DO fix one binomial number per cell — the
single-value pair efficiency `ε^pair_2mu4` (`mc_trigeff_single_value_pair_eff.md` §PP-1). That
number is pure MC, and the analysis takes its single-muon efficiencies from DATA tag-and-probe;
the MC/data pair over-efficiency is measured at `<r1 r2> = 1.27` (`mc_trig_eff_closure.md` R1).
The user's decision (parent D4) is to correct it by the **product of the two single-muon data/MC
scale factors**, per pair — not by the cell-averaged K, which would re-introduce the
factorization the single-value procedure exists to avoid.

Below 74 GeV the ΔR fit cells move to the **3-group |η^pair| fold** (statistics ×3 per fit; the
correlation barely depends on the SIGN of pair η — `mc_trigeff_dr_binning_approaches.md` D11,
closure ranking R4-retracted/R5: the fold is the BEST of the four approaches after the lookup-bug
fix) and to the **no-a-priori-plateau** fit family (free baseline C determined by the ΔR < 1 data
alone; the far plateau window never enters). The fit-form choice per cell is the user's, made on
the regenerated fit figures: **expo everywhere except three forward cells where the polynomial
describes the points and expo rails** (§3b).

### 2. Top-level equation

For a pair in the pp24 signal region (`analysis_overview.md` §2) with pair pT p, pair η, opening
ΔR, and legs i = 1,2 at (pT_i, q_i·η_i):

```
                    ⎧ ε^nc_data(1) · ε^nc_data(2) · ε_ΔR(ΔR ; cell_A(p, |η|))          p_lo ≤ p < p_split   [region A]
ε_trig^pair(pair) = ⎨
                    ⎩ ε^pair_MC(cell_B(p, |η|); OS, signal window) · SF(1) · SF(2)      p_split ≤ p < p_hi   [region B]

SF(i) = ε^nc_data(pT_i, q_i·η_i) / ε_MC(pT_i, q_i·η_i)          w_trig = 1 / ε_trig^pair
```

- `p_lo, p_split, p_hi` = `ParamsSet::pair_pt_coarse_bins` edges 0, N−2 and N (N = 8): 9,
  **74.24**, 150 GeV. Region B = the LAST TWO coarse bins, [74.24, 105.53) and [105.53, 150).
  Never retyped; read from the vector.
- `ε^nc_data` = the data tag-and-probe single-muon mu4 turn-on (unchanged;
  `RDFBasedHistFillingData::EvaluateSingleMuonEffcy`, TF1 at the exact pT, cap 1, floor 0.01).
- `ε_MC` = the MC direct single-muon mu4 efficiency, `FitMCSinglesEffcy` fit file
  (`FullSimMCSinglesFitFile`), the SAME ε_MC the Step-3 inverse weighting and K divide by
  (`SingleMuEffEvaluator::Src::kMCDirect`: TF1 at exact pT clamped into the fit range, cap 1,
  floor 0.02).
- `ε_ΔR(ΔR; cell) = f_cell(ΔR)/C_cell` for ΔR < 1, exactly 1 for ΔR ≥ 1 (unchanged form,
  `dr_correction_apply.h`).
- `ε^pair_MC(cell)` = S1/S0 of `mc_trigeff_single_value_pair_eff.md` §2, the PURE form, in the
  SIGNAL mass window, opposite sign, on the un-merged 8 × 3 cells (`nomerge`).
- `ε_trig^pair` is capped at 1 (an efficiency), with the cap COUNTED and printed.

### 3. Step-by-step method

**(a) Region A cells** = `ParamsSet::pair_pt_coarse_bins` bins 1–6 × the 3 sign-independent
|η^pair| groups `MakeDrEtaGroups(..., true)` builds from
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` — |η| < 1 / [1, 2) / [2, 2.2) — i.e.
the `nocorr_etamerge` plateau mode, **opposite-sign series**, Step-3 fit files
`dr_correction_fits_pp24_full_step3_<method>_os_nocorr_etamerge.root`. 18 cells.

**(b) Fit form per region-A cell.** Primary = `expo` (`f = C + A exp[−(ΔR/λ)^p]`, A ≤ 0,
p ≥ 1) in every cell EXCEPT the three **(pair-pT bins 2, 3, 4 = [12.8,18.2), [18.2,25.8),
[25.8,36.7)) × (the last |η| group, [2.0, 2.2))**, whose primary is the constrained polynomial
`polyu_fixedRp`. The cells are named by their INDEX on the canonical axes (1-based pT bin 2..4,
the LAST η group) and their physical ranges are printed from the axes at load time.
**Fallback** for every region-A cell = the linear **interpolation** (`interp`) if the primary is
rejected by the producer's `fit_ok` or by the two baseline screens (`DrCorrPlateauUsable`,
`DrCorrBaselineConsistent`). **No polynomial fallback for expo cells, no raw-bin tier**: a
region-A cell with neither its primary nor `interp` accepted THROWS at load time (user, 2026-09-17;
never triggers on the 2026-09-10 fits — 18/18 primaries accepted).

**(c) Region B cells** = pair-pT bins 7, 8 × the same 3 |η| groups, from
`pair_trig_eff_pp24_full.root` via `PairTrigEffEvaluator` (sign `os`, window `sig`, `kPure`,
mode `nomerge`), gated by its delivery rules (≥ 50 raw pairs, value in (0, 1]). 6 cells, 5
delivered. **D1 (TEMPORARY):** the one refused cell, [105.53, 150) × [2.0, 2.2) (3 raw pairs), is
served from the **pT-merged** cell [74.24, 150) × [2.0, 2.2) of the SAME file (mode `ptmerge`, 62
raw pairs). Any other refusal THROWS.

**(d) The SF** is evaluated per leg at the exact (pT, q·η) from the two turn-on sets of §2, so
it carries the same kinematic resolution as the region-A singles product. It is NOT clipped
(a SF > 1 is a legitimate data-over-MC region); only the final ε_trig^pair is capped at 1.

**(e) Application.** ONE helper (`RDFBasedHistFillingPP::AddPairEfficiencyWeightColumns`) builds
`effcy_pair` and `w_trig` for every pull — the OS signal-region crossx, the SS crossx histogram,
the template-fit pass and the generic / MC-data-comparison histograms — so there is exactly one
definition of the weight. Consequences, both user decisions (2026-09-17):
- **D2 (TEMPORARY):** SAME-SIGN pairs are weighted with the OPPOSITE-sign-derived numbers
  (the SS single-value cells above 74 GeV hold 33/6/0 and 6/0/0 raw pairs — unmeasurable).
- **D3:** the region-B number is the SIGNAL-WINDOW one for EVERY pair, including pairs outside
  1.08–2.9 GeV in the template-fit and generic passes (the 1–4 GeV window differs by 10–21 %,
  parent R2 item 3; stated limitation).

**(f) Correction stages** for the comparison plot: `cw_raw` (1/L) → `cw_unfolded_trig`
(× w_trig) → `cw_unfolded_reco_trig` (× w_reco), i.e. the trigger correction is shown FIRST,
as the user asked; the existing reco-first stage `cw_unfolded_reco` stays booked.

### 4. Negative constraints

- Do NOT apply ε_ΔR^single (Step 4) to pp — the 2mu4 product absorbs it (`mc_trigger_efficiency.md` §2).
- Do NOT multiply ε^pair_MC by the single-muon efficiencies: it is a PAIR-level efficiency, not a
  ΔR correction factor. Only the data/MC SCALE FACTORS multiply it.
- Do NOT use K (the calibrated form) anywhere in the application (parent D4).
- Do NOT route ε^pair through `DrCorrectionEvaluator` (parent §4).
- Do NOT change any binning: pair-pT cells from `ParamsSet::pair_pt_coarse_bins`, |η| groups from
  the live `MakeDrEtaGroups` fold, region split = the last two canonical bins.
- Do NOT retype the mass window: `ParamsSet::signal_minv_min/max` only.
- Pb+Pb is out of scope (mu4 UNION weight, a different formula).

---

## Context

- The previous application (`pp24_crossx_rerun_2026_08.md` §3c; `dr_correction_sample_cfg.h`
  `DrCorrCrossxMode() = nocorr_ptmerge`, expo → polyu → raw) was flagged TEMPORARY at the time.
  `mc_trigger_efficiency.md` R35 (2026-09-10) found its raw-bin tier delivering ε_ΔR down to 0.044
  and the polynomial tier never reached.
- Inputs verified on disk 2026-09-17, all on the 9 GeV / 4.5 GeV selection, Tight WP:
  `mc_trig_eff_hists_pp24_full_step3.root` 2026-09-10 06:01; Step-3 `nocorr_etamerge` OS fits
  (expo 14:51, polyu 14:52, interp 14:53); `pair_trig_eff_pp24_full.root` 16:12;
  `single_mu_effcy_pT_fit_mc_pp24_full.root` 06:00; data T&P `single_mu_effcy_pT_fit.root`.
  **Medium WP inputs are STALE** (pair-eff file 2026-09-08 on the 8 GeV axis; Medium ΔR fits not
  refit after the 2026-09-10 Medium refill) — the code keeps the WP switch, a Medium run throws on
  the binning guards until those are regenerated (parent R35 scope decision).
- A concurrent session's `run_pbpb_all.sh YEARS="23 25 26"` (PID 2579325) is running in this
  checkout (event selection at 2026-09-17 00:00). It will compile `RDFBasedHistFillingPbPb.cxx`
  and run `plot_crossx_trig_corr_sanity.C()` later. Shared-file edits here are value-neutral for
  Pb+Pb (mass window constant, one extra stage column) and are compile-checked immediately.

## Scope

**In:** ParamsSet mass window + all live signal-selection sites; the crossx-side hybrid evaluator;
`RDFBasedHistFillingPP` weight columns + trigger-first stage; `CorrectionStages.h` + the PbPb
one-line Define; the two sanity-plot macros; pp24 pipeline Stages 5–8; docs; reviews; commits.

**Out:** template-fit changes (hand-over note below); Pb+Pb; R_AA (stale after this — see
Remaining Work); Medium-WP regeneration; the open R12 bound on the delivered f/C; any binning.

## Design Decisions

- **D1 (user, 2026-09-17) — refused region-B cell → pT-merged value, TEMPORARY.** The
  [105.53,150) × [2.0,2.2) OS cell has 3 raw MC pairs. Served from the [74.24,150) × [2.0,2.2)
  merged cell of the same file. **Marked in code and doc: a request for additional high-pT MC
  statistics is being prepared; once that sample is in, revert this |η| group to the nominal
  un-merged cells** (delete the fallback, not widen it).
- **D2 (user, 2026-09-17) — SS pairs keep the OS-derived weights, TEMPORARY.** A same-sign-filtered
  high-pT MC sample is being requested specifically for this region; once available, SS pairs get
  the SS series (ΔR fits below 74 GeV, SS single-value cells above).
- **D3 (user, 2026-09-17) — signal-window number everywhere.** The template fit is to be
  performed in the SAME 1.08–2.9 GeV window in future (no φ / J/ψ peaks in the fit); if resonance
  leakage into the window is found to be significant the window may shrink and the new window
  applies everywhere — determining that is NOT this task. Template-fit changes are OUT of scope →
  hand-over note in Remaining Work.
- **D4 (user, 2026-09-17) — fallbacks.** poly → interp in the three poly cells, expo → interp
  elsewhere, no raw-bin tier, throw when nothing usable.
- **D5 — the mass window becomes a `ParamsSet` constant** (`signal_minv_min = 1.08`,
  `signal_minv_max = 2.9`, doubles, formatted with `%g` so every generated selection string is
  BYTE-IDENTICAL to the retyped `minv > 1.08 && minv < 2.9` → bit-identical outputs, no rerun
  blast radius from this step alone). `PairTrigEff::Windows()` reads it; `CheckSignalWindowMirror`
  stays as a second, independent guard.
- **D6 — trigger-first stage as an ADDITIONAL stage**, not a re-ordering of the existing table:
  `cw_unfolded_trig` is added to `CorrectionStages.h` (shared) with the one-line Define in both
  fillers; the stage plot draws raw → +trig → +trig+reco for pp and falls back to the old
  reco-first triple for a file that lacks the new stage (Pb+Pb on disk until its next refill).

## Implementation Plan

*(all six steps DONE 2026-09-17; `/review-analysis-code` PASS at iteration 2, `/review-plot` PASS at iteration 3)*

1. [x] `ParamsSet.h` mass window + `SignalMinvCutExpr`; replace the live sites (§Done 1). Per §4.
   Reviewer: `/review-analysis-code`.
2. [x] `dr_correction_sample_cfg.h`: the crossx configuration (mode `nocorr_etamerge`, primary
   `expo`, poly cells, fallback `interp`, region-B sign/window/mode, the D1 merged fallback flag)
   — named in ONE place. `Utilities/DrCorrectionCrossxEvaluator.h`: three-method per-cell routing
   + throw (§3a/b). NEW `Utilities/PairTrigEffCrossxEvaluator.h`: the §2 hybrid (region split,
   region-B lookup + D1, SF, cap, census). Reviewer: `/review-analysis-code` with §2, §3, §4.
3. [x] `RDFBasedHistFillingPP.cxx`: `AddPairEfficiencyWeightColumns` → the hybrid;
   `cw_unfolded_trig`; `CorrectionStages.h` + PbPb Define. Reviewer: `/review-analysis-code`.
4. [x] `plot_crossx_reco_eff_stages.C` (trigger-first order + fallback) and Stage 7 of
   `pipeline_pp_crossx.sh` runs it. Reviewer: `/review-plot`.
5. [x] Run `SKIP_CONDOR=1 pipeline_pp_crossx.sh`; validate; record the census + numbers here.
6. [x] Reviews, docs, commits.

## Progress Log

*(append-only; newest at the END)*

- 2026-09-17 — Doc created after triage + input verification (Context) and the four user
  decisions D1–D4 (AskUserQuestion). Plan written; no code touched yet.

- 2026-09-17 01:30–01:50 — **Steps 1–4 written and compiled.**
  * Step 1: `ParamsSet::signal_minv_min/max` + `SignalMinvCutExpr()` (+ `PassSignalMinv`);
    verified `SignalMinvCutExpr()` == `"minv > 1.08 && minv < 2.9"` byte-for-byte. Sites moved:
    `RDFBasedHistFillingPP.cxx:592`, `PbPb.cxx:991`, `PythiaTruth.cxx:471`,
    `PythiaFullsim.cxx:298/300`, `PowhegFullsim.cxx:285/287/328`, `PowhegTruth.cxx:291`,
    `PythiaFullsimOverlay.cxx:84/89`, `Utilities/MCTrigEffPairSelection.h` (SingleBSignalCutsReco),
    `Utilities/PairTrigEffEvaluator.h` (`Windows()` sig), `pp24_secondary_vertex_stats.cxx:157`,
    `plot_dr_vs_pair_pt_diagnostic.cxx:153/159`, `plot_sig_accept_cutflow_above_60GeV.cxx:42-43`.
    All 9 affected ACLiC targets compile. Left alone (dead / legacy, still retyped):
    `SingleBAnalysis/*`, `OldMuonPairHistFillingCodePreRDF/*`, `FillHistogramsCrossx_PP_clean.cxx`,
    `plot_reco_distr_singleb_vs_op_pp24.C` (gated behind a STALE banner).
  * Step 2: `dr_correction_sample_cfg.h` crossx block rewritten (`DrCorrCrossxMode()` =
    `nocorr_etamerge`, `DrCorrCrossxMethod()` = expo, `DrCorrCrossxPolyMethod()` = polyu_fixedRp,
    `DrCorrCrossxFallbackMethod()` = interp, `DrCorrCrossxPolyPtBins()` = {2,3,4},
    `PairTrigEffCrossxSign/Window/Mode()` = os/sig/nomerge,
    `PairTrigEffCrossxRefusedCellFallbackMode()` = ptmerge (D1, TEMP)).
    `Utilities/DrCorrectionCrossxEvaluator.h` rewritten (three loaders, per-cell primary, interp
    fallback, throw, region-A-only extrema scan, canonical-binning guard on all three files,
    refuses a pT-merged mode). NEW `Utilities/PairTrigEffCrossxEvaluator.h` (the §2 hybrid).
    `Utilities/proj_range_to_suffix.cxx` got an include guard (it is included twice in the PP TU
    now, via the base class and via `SingleMuEffEvaluator.h`).
  * Step 3: `RDFBasedHistFillingPP.cxx` — `s_pair_trig_eff`, `effcy_pair` from the hybrid (11
    input columns), `cw_unfolded_trig`; `CorrectionStages.h` 5th stage `_corr_unfolded_trig`;
    `RDFBasedHistFillingPbPb.cxx` one-line Define. PP + PbPb compile.
  * Step 4: `plot_crossx_reco_eff_stages.C(include_pbpb)` — trigger-first triple with the
    reco-first fallback for files lacking the stage; `pipeline_pp_crossx.sh` Stage 7 now runs it.
  * **Standalone load test** (`scratchpad/test_hybrid.C`, pp_full Tight): region A **15 expo +
    3 poly (the three named cells) + 0 interp**, delivered ε_ΔR over the 18 cells spans
    [0.284 (bin 6, |η| group 2), 1]; region B **5 cells delivered + the D1 fallback for
    [105.5,150)×[2,2.2)** (ε^pair = 0.441 ± 0.096 from the merged 62-pair cell). Probe pairs
    behave (±η symmetric; dR ≥ 1 → singles product; pT ≥ 150 → singles product, counted;
    region B: ε^pair × SF1·SF2 with SF1·SF2 = 0.829 for 0.8/0.8576 × 0.8/0.9002).
  * Pre-run input state: `pp_2024/muon_pairs_pp_2024_2mu4*.root` 2026-09-10 16:05, data T&P fits
    2026-09-09 23:38, crossx output to be replaced `histograms_real_pairs_pp_2024_2mu4_nominal.root`
    2026-09-10 16:07 (696 kB).

- 2026-09-17 01:42–01:45 — **Step 5 DONE: `SKIP_CONDOR=1 INCLUDE_PBPB_SANITY=false
  pipelines/pipeline_pp_crossx.sh` exit 0** (log `pipelines/logs/pp_crossx_hybrid_20260917.log`).
  RDF stage 01:42:47–01:43:22; output `histograms_real_pairs_pp_2024_2mu4_nominal.root` fresh and
  filled. Regenerated: `plots/single_b_analysis/pp24/` (crossx 5 PNG + counts PNG/CSV),
  `pp24_pt_120/`, `plots/sanity_check_crossx/PP_2024_pair_pt_in_eta_subplots.png` +
  `PP_2024_reco_eff_stages_pair_pt_in_eta.png` (now uncorrected → +trigger → +trigger+reco),
  `plots/mc_data_compr/{signal,generic}/` (15 PNG). Census in R1.

## Results & Observations

### R1 — Load-time census and per-pair routing (pp24 nominal RDF stage, Tight)

- Region A: **18 cells routed — 15 expo, 3 polyu_fixedRp (exactly the three named cells:
  [12.8,18.2), [18.2,25.8), [25.8,36.7) × |η| [2.0,2.2)), 0 interp**. Delivered ε_ΔR over the 18
  cells spans [0.284 (bin 6 × |η| group 2, at ΔR → 0), 1]. Per pair: 1 295 666 region-A
  evaluations, 96.5 % expo / 3.5 % poly, 0 outside the grid; 40 % of the expo evaluations are at
  ΔR ≥ 1 (no correction), 0 floored, 0 capped.
- Region B: **5 cells delivered + the D1 fallback** for [105.53,150) × [2.0,2.2) (ε^pair =
  0.441 ± 0.096 from the 62-pair merged cell); 934 region-B evaluations (0.04 % of all pairs
  seen by the helper, incl. the generic sample), **14 on the D1 fallback**; MC ε for the SF: 1868
  evaluations, 26 capped at 1 (1.4 %), 0 floored; **0 pairs capped at ε_trig = 1**.
- 849 662 evaluations (39.6 %) outside [9,150) GeV — the generic sample's pair pT < 9 GeV, singles
  product only, exactly as before this change.

### R2 — The delivered weights and the spectrum across the 74.24 GeV split (signal region, OS)

Mean per-pair weights `<w_trig>` (= Σ trig-corrected / Σ raw) and `<w_all>` (× reco) per fine
pair-pT bin (`ParamsSet::pT_bins_150`), and N = raw pair counts:

| pT bin [GeV] | \|η\|<1: N / <w_trig> / <w_all> | 1≤\|η\|<2 (+): N / <w_trig> | 2≤\|η\|<2.2 (+): N / <w_trig> | all η: N / <w_trig> / <w_all> |
|---|---|---|---|---|
| [43.8,52.2) | 711 / 3.21 / 4.04 | 474 / 1.48 | 49 / 1.81 | 1745 / 2.20 / 2.73 |
| [52.2,62.3) | 334 / 4.06 / 5.17 | 203 / 1.63 | 20 / 1.86 | 768 / 2.69 / 3.43 |
| [62.3,74.2) | 114 / 4.46 / 5.68 | 64 / 1.73 | 12 / 1.64 | 273 / 2.84 / 3.62 |
| **[74.2,88.5)** | 44 / **6.60** / 8.39 | 30 / 1.89 | 5 / 2.32 | 112 / 3.77 / 4.77 |
| [88.5,105.5) | 15 / 6.49 / 8.15 | 11 / 1.90 | 2 / 2.26 | 32 / 4.07 / 5.11 |
| [105.5,125.8) | 7 / 6.15 / 7.73 | 2 / 2.46 | 0 | 12 / 4.59 / 5.86 |
| [125.8,150) | 2 / 5.79 / 7.15 | 1 / 2.42 | 0 | 4 / 4.11 / 5.21 |

Reading (corrected after the two 2026-09-17 reviews). The barrel trigger weight steps ×1.48
(4.46 → 6.60) between the fine bins on either side of the split, ×1.58 between the coarse
cells (cell-mean `<w_trig>` 3.10 / 4.16 / 6.57 / 6.07 for the four barrel cells from 36.7 GeV).
Most of that is the MC's own statement: the ΔR-integrated correlation factor K of the
single-value file is 0.575 / 0.400 / 0.273 / 0.282 in the same cells (cell ratios 1.44 / **1.47** /
0.97 — the close-by loss grows with pair pT), so of the ×1.58 at the split, ×1.47 is K and
**~×1.07 is the offset between the two procedures**. That offset is real and measured
(parent `mc_trigeff_single_value_pair_eff.md` R9, on the superseded 8 GeV selection): the
**fold-ΔR** procedure the hybrid uses in region A corrects *less* than the single value by
C_D/C_sv = 0.949 (56–68 GeV) and 0.908 (68–83 GeV); the "~1 %" agreement quoted there is the
pT-MERGED ΔR mode (C_B/C_sv), not the fold. Inside region B the weight is flat across the two
fine bins of each coarse cell (one number per cell, modulated only by SF₁SF₂). The η-integrated
spectrum falls monotonically through the split (0.655, 0.206, 0.0935, 0.0240, 0.0087, 0.0022
pb/GeV from 52 to 150 GeV). See R3 for the seam measured on the data pairs themselves and R4 for
the fine-bin imprint; both are user-facing.

### R3 — The seam at 74.24 GeV measured on the DATA pairs (code review, 2026-09-17)

On the selected OS signal-region pairs of the LAST region-A cell (coarse bin 6, [52.2,74.2)),
the applied region-A efficiency `<ε₁ε₂ε_ΔR>` versus what the region-B recipe would give the SAME
pairs, `ε^pair_MC(bin 6) × <SF₁SF₂>`:

| \|η^pair\| | N | <ΔR> | <ε_A> | ε^pair × <SF₁SF₂> | ratio |
|---|---|---|---|---|---|
| < 1 | 448 | 0.087 | 0.2533 | 0.3181 × 0.6950 = 0.2211 | **1.146** |
| 1–2 | — | — | 0.6174 | 0.6097 × 0.9520 = 0.5805 | 1.064 |
| 2–2.2 | — | — | 0.6122 | 0.6066 × 0.9644 = 0.5850 | 1.046 |

(bin 5, [36.7,52.2): 1.049 / 1.066 / 0.941; region B reproduces itself exactly, as it must.)
The two recipes disagree by **5–15 % on the same data cell**, most in the barrel where the
weight is largest (the MC ε^pair stat error there is 3.4 %). Not a code bug — it is what §2
implies: the parent R9 closure identity is MC-internal and cell-integral-preserving, whereas on
data the two recipes differ through the pairs' ΔR distribution and are not pinned to each other.
**Open, user decision:** whether the seam disagreement becomes a systematic on the trigger
correction, or whether the split edge / a recipe overlap is revisited when the requested high-pT
MC arrives. RUN2-CROSSCHECK UNVERIFIED (no Run 2 nearby-pair dσ/dpT to compare to).

### R4 — Fine-bin imprint of the coarse-cell corrections, and the [52.2,62.3) outlier (plot review, 2026-09-17)

The corrections are one number per coarse cell (region B) or one curve per coarse cell (region
A) while the presentation axis is 2:1 finer. Barrel fine-bin `<w_trig>` = 3.05, 3.21 | 4.06, 4.46
| 6.60, 6.49 | 6.15, 5.79: within-cell rise 5–10 %, across-cell jumps 26 % (52.2 GeV) and 48 %
(74.2 GeV). Local-power-law pulls of the η-integrated corrected spectrum, [52.2,62.3) …
[125.8,150): +3.8, −2.5, +2.2, −1.0, +0.6, −0.8σ, against raw counts +2.7, −0.8, +1.6, −1.0,
+0.4, −0.3σ — the corrected−raw difference alternates in phase with the cell edges (+1.1, −1.7,
+0.6, 0, +0.2, −0.5σ), i.e. a ±5–12 % sawtooth that is a systematic of the granularity, not a
fluctuation (no single bin exceeds ~1.7σ of it). At coarse granularity the spectrum is smooth
through the split: rebinned to the 8 canonical cells, [74.2,105.5) sits at 0.87 (−1.6σ) and
[105.5,150) at 0.68 (−1.8σ) of the local power law, and Pythia/data per coarse cell is
continuous, 1.41 ± 0.03 / 1.26 ± 0.05 / 1.20 ± 0.12 / 1.19 ± 0.33 (η-integrated).
**User decision:** present the corrected spectrum above 52 GeV at coarse granularity, or carry
the imprint as a systematic. No binning is changed here.

The **[52.2,62.3) η-integrated bin remains a +3.8σ (×1.18) local-power-law outlier** — R13 item 1
of `mu_pt45_gap125_pairpt9_adoption.md`, still open, NOT introduced here: the raw counts are
already +2.7σ (×1.11) there, the rest is the ×1.26 weight step at the 52.2 GeV coarse edge; in
the [0.5,1.0) panel [43.8,52.2) is ×0.76 (−4.1σ) and [52.2,62.3) ×1.46 (+3.0σ) with raw −2.4σ /
+2.8σ, and Pythia/data drops 1.69 ± 0.15 → 0.98 ± 0.11 across it.

### R5 — Amendments after the two reviews (both FAIL with WARNINGs only, no CRITICAL)

Code: inner census counters made atomic (`SingleMuEffEvaluator`, `DrCorrectionEvaluator`; the
first run's inner counts were ~0.4 % low under ImplicitMT — delivered values never affected);
D1 fallback now THROWS if more than ONE region-B cell is refused (§3c as written); the SF
evaluation asymmetry is documented in code (data turn-on = formula at exact pT, uncapped in pT,
cap 1, floor 0.01; MC twin clamped into [4.5,60], cap 1, floor 0.02 — ≤ 1 % per leg above
60 GeV, floors never active); stale header comments refreshed; comment-quoted edges marked
illustrative. Plots: `plot_crossx_reco_eff_stages.C` decides the triple per FILE (a mixed-year
Pb+Pb set falls back as a whole), gained `use_tight_wp` (+ `_medium_wp` file/PNG suffix) and a
"pp 2024, 2mu4, tight WP" label per panel, ratio pads keep the line at 1 in frame;
`plot_crossx_trig_corr_sanity.C` gained the same WP argument and label;
`SingleBCrossxPlotterBase.cxx` 2D maps draw log-z down to the smallest positive bin (the default
10⁻³·z_max floor had hidden every bin above ~62 GeV).

## Remaining Work

- **R_AA is STALE**: it uses this pp24 file as denominator; the concurrent Pb+Pb driver
  (`run_pbpb_all.sh YEARS="23 25 26"`, PID 2579325) refills the Pb+Pb crossx inputs and does NOT
  run R_AA itself — rerun R_AA (`run_year_trigger_mode 6`) after it finishes, never mid-refill.
- **User decisions from the reviews (R3, R4):** the 5–15 % seam disagreement at 74.24 GeV as a
  systematic vs. a recipe change; coarse vs fine presentation above 52 GeV (or a granularity
  systematic); the still-open [52.2,62.3) outlier (`mu_pt45_gap125_pairpt9_adoption.md` R13 item 1).
- **HAND-OVER (D3): the low-mass template fit must move to the signal mass window** (1.08–2.9 GeV,
  `ParamsSet::signal_minv_min/max`) — no φ / J/ψ / ψ(2S) peaks in the fit. Not done here.
- R_AA becomes STALE the moment the pp24 crossx file is refilled (pp is the denominator). The
  concurrent Pb+Pb run regenerates R_AA at its end; if it ran BEFORE this refill, rerun R_AA.
- Medium WP: pair-eff file + ΔR fits to be regenerated before a Medium crossx run.
- Revert D1 and D2 when the requested high-pT (and SS-filtered) MC samples arrive.

## Latest Stage

**2026-09-17 02:30 — DONE.** Code, run, both reviews (code PASS iter 2; plots PASS iter 3), docs.
Committed by explicit path (see Progress Log). Open for the user: Remaining Work (seam systematic
vs recipe, coarse vs fine presentation above 52 GeV, the [52.2,62.3) outlier, R_AA rerun, D1/D2
reverts when the requested MC arrives, the template-fit window hand-over).
