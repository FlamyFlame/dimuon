# pp24 cross-section rerun with gap cuts, MC dR-corrected trigger efficiency, and fullsim reco efficiency

Mode: **IMPLEMENTATION**. Opened 2026-08-17.

## Objective

Update the code and rerun the **pp24 cross-section pipeline** end-to-end so that the pp24
crossx plots and the pp-data/MC comparison plots are produced with:

1. A **detector-gap fiducial cut applied to BOTH muons** of every pair
   (`ParamsSet::single_mu_fiducial_gap_cuts`, latest values), replacing the standalone
   per-muon `q*eta < 2.2` signal cut.
2. The **one-sided Delta-p/p cut** (`dP_overP < 0.12`, negative tail kept) everywhere.
3. Trigger efficiency = **data-derived single-muon mu4 eps^nc(pT, q*eta)** per leg, times the
   **MC-derived 2mu4 dR-correlation correction eps_dR(dR)** in the configuration:
   no-plateau-correction, **opposite-sign**, **7 pair-pT bins** (canonical 8-bin binning with
   the last two bins combined), **exponential** fitted function, falling back to the
   **polynomial** fit where the exponential is rejected and to the **raw (absolute) bin
   values** where both fits are rejected.
4. Reconstruction efficiency derived from the **pp24-condition Pythia fullsim sample**
   (replacing the Run-2 placeholder).

Explicitly OUT of scope: unfolding, signal-selection acceptance, template fitting, PbPb
T_AA-weighted event yields, R_AA.

## Autonomy Contract (DONE 2026-08-18 — all six Done items met)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. Gap cut (both muons, from `ParamsSet::single_mu_fiducial_gap_cuts`) wired into the pp24
     crossx signal selection and the pp24 fullsim reco-efficiency selection, replacing the
     `q*eta < 2.2` cut; all affected classes recompiled.
  2. One-sided Delta-p/p verified in force in every pp24 stage; pp24 NTuple outputs
     regenerated if they predate the fix.
  3. pp24 crossx run with the trigger weight = data eps^nc per leg x MC eps_dR in the
     configuration of Objective item 3 (7 pair-pT bins; the 7-bin input comes from the
     sibling agent working in `mc_trigger_efficiency.md` — use the 8-bin file only as an
     interim and re-run once the 7-bin file exists).
  4. pp24 crossx run with the pp24-fullsim-derived reco efficiency in place of the Run-2
     placeholder.
  5. Regenerated at their canonical paths: pp24 crossx plot set
     (`plots/single_b_analysis/pp24/` + the `pp24_pt_150/` opt-in variant), the trigger-
     correction sanity plots, and the MC-data comparison plot set
     (`plots/mc_data_compr/`).
  6. Docs updated (this doc, `signal_selection_change_impact.md` if the map changed,
     `analysis_overview.md` §2 signal region, INDEX.md) and the work committed.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

The pp24 single-b dimuon cross-section is the denominator of R_AA. Three of its inputs have
moved and one is a placeholder:

* **Detector gaps.** In the barrel crack (q*eta ~ 0), the barrel/endcap transition
  (q*eta ~ -1.1) and the forward edge (q*eta > 2.3) the muon reconstruction and trigger
  efficiencies are small and vary fast, so an efficiency *correction* there carries a large
  systematic. The standard alternative is a **fiducial cut** plus an acceptance factor.
  `muon_gap_cuts_acceptance.md` settled the cut set and showed (F12, from TRUTH) that all
  three features are 100 % detector, i.e. the cut is a pure acceptance loss.
* **Delta-p/p** is a ONE-SIDED cut (`dp/p < 0.12`; the negative tail is physical and is kept).
  The pp24 nominal trees on disk were produced before that fix.
* **Trigger.** The 2mu4 pair efficiency is NOT the product of two single-muon efficiencies:
  the two legs are correlated at small opening angle. The MC-derived correction eps_dR
  supplies that correlation; until now the cross-section applied eps_1*eps_2 with no eps_dR.
* **Reconstruction.** The cross-section applies a Run-2 single-muon PLACEHOLDER
  (`run2_reco_eff_placeholder.root`, eps_1*eps_2, no dR dependence). The pp24-condition
  Pythia fullsim sample now exists and supplies the genuine PAIR efficiency.

### 2. Top-level equation

For opposite-sign pairs in the single-b signal region, per (pair p_T, pair eta) cell:

```
dsigma/dX = (1/L_int) * SUM_pairs  1 / [ eps_trig^pair * eps_reco^pair ]
```

with `L_int` = 400.412 pb^-1 (pp24, 2mu4) and

```
eps_trig^pair = eps^nc(p_T1, q*eta1) * eps^nc(p_T2, q*eta2) * eps_dR^2mu4(dR ; pair p_T, pair eta)
eps_reco^pair = eps_reco(pair p_T, pair eta, dR)
```

* `eps^nc` — the DATA tag-and-probe single-muon mu4 turn-on, fitted vs p_T in contiguous
  coarse q*eta bins (`CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`, -2.4 <= q*eta < 2.3).
* `eps_dR^2mu4` — the MC-derived pair-trigger correlation correction (see §3c).
* `eps_reco` — a genuine PAIR reconstruction efficiency, NOT eps_1*eps_2 (see §3d).

**No unfolding** (`w_unfold == 1`), **no signal-selection acceptance** alpha, **no template-fit
background subtraction**, and **no acceptance efficiency eps_acc** are applied. The result is
therefore a **fiducial, reco-level cross-section** inside the gap-cut fiducial region.

### 3. Step-by-step method

**(a) Signal region (both muons).** Opposite-sign pair with

```
m_uu in (1.08, 2.9) GeV,  pair p_T > 8 GeV,
BOTH muons pass the detector-gap fiducial cut: q*eta = charge*eta is NOT inside any window of
ParamsSet::single_mu_fiducial_gap_cuts = {(-1.20,-1.05), (-0.10,+0.06), (2.30,2.40)}.
```

No dR cut (removed 2026-06-22). The gap cut **REPLACES** the former per-muon one-sided
`q*eta < 2.2`. Since the ntuple stage already requires |eta| < 2.4, the forward window makes
the effective forward edge **2.30**, which is exactly the top edge of the coarse q*eta turn-on
binning — so every surviving muon has a fitted eps^nc and no pair can be silently dropped.
The window values are READ from `ParamsSet` at every call site and never retyped.

**(b) Delta-p/p.** `dp/p < 0.12` per muon, ONE-SIDED, applied at the NTuple stage
(`DimuonDataAlgCoreT.c:601`). The negative tail is KEPT. Requires regenerating the pp24
nominal pair trees, which predate the fix.

**(c) Trigger efficiency.** Per pair, `eps_trig^pair = eps^nc_1 * eps^nc_2 * eps_dR`.
`eps_dR` is measured in the pp24 Pythia fullsim by inverse weighting
(`mc_trigger_efficiency.md` §3.3): the ratio of pairs firing 2mu4, each weighted
`1/(eps_MC,1 * eps_MC,2)`, to all pairs with no trigger requirement, as a function of dR, in
cells of (pair p_T, pair eta). The configuration this analysis applies is, by user instruction:

| knob | value | token |
|---|---|---|
| plateau correction | **none** — the raw eps_dR is fitted with a free additive baseline C and the delivered correction is `f(dR)/C` | `nocorr` |
| pair sign | **opposite-sign** (truth charges) | `_os` |
| pair-p_T cells | the canonical 8-bin `ParamsSet::pair_pt_coarse_bins` with the **last two bins merged** -> **7 cells**, top cell [72.08, 150) GeV | `_ptmerge` |
| functional form | **exponential** `f = C + A exp[-(dR/lambda)^p]`; where that fit is REJECTED fall back to the **polynomial** `f = C + u^2 (a2 + a3 u + a4 u^2)`, u = max(0, 1 - dR/R_p); where BOTH are rejected use the **raw (absolute) binned values** normalised by their own mean over dR in [0.5, 1] | `expo` -> `polyu_fixedRp` -> raw |

"Rejected" is the producer's own verdict, `h_step3_fit_ok(cell) == 1` (fit converged, >= nfree+2
informative points, and `f(dR) > 0` over [0, R_p]); consumers must also require the fitted
baseline C to be usable (`DrCorrPlateauUsable(C, err_C)`), because a near-zero C inflates `f/C`.
`eps_dR` is applied for dR < 1 only and is exactly 1 above (the fit domain is dR in [0,1]).
Pair-eta cells are `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` (9 bins).

**NEGATIVE CONSTRAINT.** The single-leg correction `eps_dR^single` (Step 4) must **NOT** be
applied to pp: the 2mu4 product already absorbs the single-leg dR dependence, and applying both
double-counts. Step 4 is a Pb+Pb (mu4 union) object.

**NEGATIVE CONSTRAINT.** MC supplies the dR RATIO only. The single-muon efficiencies in the
cross-section are the DATA eps^nc, never eps_MC (the MC trigger simulation is ~1.2x
over-efficient per leg; that cancels inside the eps_dR ratio and must not be imported).

**(d) Reconstruction efficiency.** A single PAIR efficiency measured in the pp24-condition
Pythia fullsim FULL sample (`_pdf` production, DSIDs 803015-803020), binned in
**(pair p_T, pair eta, dR)**:

```
eps_reco(pair p_T, pair eta, dR) =
      N[ single-b OS pair, both muons truth-matched to reco, pair passes the Tight WP,
         and the RECO pair passes the §3a signal region (gap cut on RECO q*eta) ]
    / N[ single-b OS pair whose TRUTH pair passes the §3a signal region
         (gap cut on TRUTH q*eta) ]
```

both binned in TRUTH kinematics and weighted by the MC event weight.

* It is **not** eps_1 * eps_2: the two close-by muons share ID hits and compete in ambiguity
  resolution, and the dR axis is what captures that correlation
  (`analysis_overview.md` §4b).
* **The gap cut sits in BOTH numerator (reco q*eta) and denominator (truth q*eta).** Therefore
  `eps_reco` is a *fiducial* efficiency and does **not** contain the truth-level gap acceptance
  `eps_acc = 0.9133` (`muon_gap_cuts_acceptance.md` F12). `eps_acc` remains a separate factor
  and is deliberately NOT applied here (the user deferred acceptance).
* It is evaluated at the **reco** kinematics of the data pair. This is the standard
  pre-unfolding approximation and is exactly what the placeholder it replaces does; it becomes
  exact once detector-response unfolding lands, at which point the correct ordering is
  unfold-then-apply-eps_reco-in-truth-bins (`CorrectionStages.h`).

### 4. Negative constraints (things the code must NOT do)

1. Do **not** apply `eps_acc` — it is not built, and folding it into `eps_reco` (by putting the
   gap cut only on the reco leg) would silently turn the fiducial result into a full-acceptance
   one and would double count the day `eps_acc` is added.
2. Do **not** apply the Step-4 single-leg `eps_dR^single` to pp (see §3c).
3. Do **not** use `eps_MC` single-muon turn-ons in the cross-section weight (see §3c).
4. Do **not** re-invent any binning. Pair p_T = `ParamsSet::pair_pt_coarse_bins`; pair eta =
   `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`; single-muon q*eta =
   `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`; gap windows =
   `ParamsSet::single_mu_fiducial_gap_cuts`. Read them, never retype them.
5. Do **not** touch the Pb+Pb crossx, R_AA, the truth signal acceptance, or the template fits in
   this task. Consequence to keep visible: after this change the pp24 and Pb+Pb signal regions
   DIFFER (pp = gap cut, Pb+Pb = still `q*eta < 2.2`), so R_AA must not be recomputed until
   Pb+Pb is brought over.

## Latest Stage

**DONE 2026-08-18.** All six Autonomy-Contract Done items are met; the pipeline ran clean end to
end and both mandated reviews returned PASS-WITH-COMMENTS with every finding closed or explicitly
carried forward (see Remaining Work). Nothing is in flight.

## Progress Log

- 2026-08-17 Step 0: doc created; Autonomy Contract pinned.

- 2026-08-17 Step 1: **pp24 nominal NTuple Condor rerun SUBMITTED (cluster 17, 12 jobs)**.
  Reason: `muon_pairs_pp_2024_part*_2mu4_mindR_0_02.root` on disk are from 2026-06-10, i.e.
  they PREDATE the one-sided Delta-p/p fix `71fcf1c` (2026-07-16 17:00). The cut lives at the
  NTuple stage (`NTupleProcessingCode/DimuonDataAlgCoreT.c:601`,
  `m1.dP_overP > thrsh || m2.dP_overP > thrsh` — one-sided, no fabs), so the trees, not just
  the histograms, are stale. `run_pp_24_nominal.sub` (queue 12) -> `run_pp_24_nominal.sh`
  (`PPAnalysis(24, batch); trigger_mode = 3`). No NTuple-code change is needed for the gap
  cut: the gap cut replaces the per-muon `q*eta < 2.2` SIGNAL cut, which lives in the RDF
  stage, not in the trees.

- 2026-08-17 Step 2 (gap cut wired in, per Physics Procedure §3a). The per-muon one-sided
  `q*eta < 2.2` was REPLACED by `ParamsSet::FiducialGapCutExpr(...)` on BOTH muons at every site
  in scope; the window values are read from `ParamsSet::single_mu_fiducial_gap_cuts` and are
  never retyped:
  * `RDFBasedHistFillingPP.cxx` `signal_cuts` (gates the OS crossx AND the SS crossx) and BOTH
    `signal_cuts_no_minv` (template-fit pass + nominal low-mass pass).
  * `RDFBasedHistFillingPythiaFullsim.cxx` `pass_signal_truth` (TRUTH q*eta) and
    `pass_signal_reco` (RECO q*eta).
  * `Utilities/MCTrigEffPairSelection.h` `SingleBSignalCutsReco()` — the data-like mirror; it now
    calls the same `FiducialGapCut()` the Step-3 selection uses, so mirror and original cannot
    drift.
  * `ParamsSet.h` STATUS block rewritten (it claimed the vector was trigger-efficiency-only).
  NOT changed, and deliberately so (out of scope): PbPb crossx, PythiaFullsimOverlay,
  Pythia/Powheg truth acceptance, PowhegFullsim. **The pp24 and PbPb signal regions therefore
  DIFFER until Pb+Pb is brought over — R_AA must not be recomputed in between.**
  Both classes recompiled clean (ACLiC).

- 2026-08-17 Step 3 (eps_reco from the pp24 fullsim, Physics Procedure §3d).
  * The fullsim reco-efficiency chain had NO persisted ROOT product at all — only PNGs — which is
    why the crossx was still on the Run-2 placeholder. Added:
    `plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C` ->
    `pythia_fullsim_full_sample/pair_reco_eff_pp24_full.root`, and the consumer
    `Utilities/PairRecoEffEvaluator.h`.
  * BINNING: the pre-existing 3D reco-eff view uses `pT_bins_80` (12 log bins, stops at 80 GeV)
    x 48 uniform eta x 20 uniform dR — none of them canonical, and the pT axis does not even
    reach the crossx range. A SECOND, APPLIED 3D view was added on the canonical axes
    (`ParamsSet::pair_pt_coarse_bins` x
    `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` x
    `RDFBasedHistFillingPythia::dr_bins_edges_for_reco_effcy`), registered in
    `hist_binning_map` and `var1D_pythia_fullsim.json` as
    `truth_pair_pt_coarse` / `truth_pair_eta_coarse` / `truth_dr_effcy`. The fine views are
    untouched — they are diagnostics, not the applied correction.
  * **The dR axis needs NO extension above 1.0.** Measured on the pp24 fullsim single-b OS truth
    signal region: max truth dR = 0.701, and the fraction above 1.0 is exactly 0. This is
    KINEMATIC, not a cut: m_uu < 2.9 GeV with pair pT > 8 GeV forces dR ~< 2 m/pT = 0.725.
    (Without the mass cut the same pairs reach dR = 1.86.)
  * Fullsim RDF hist filling rerun (8 slots, 20 s event loop, 414 keys). eps_reco built:
    8 x 9 x 4 = 288 cells; INCLUSIVE eps_reco = **0.6620 (Tight)**, 0.7301 (Medium).
  * 144 of 288 3D cells (50 %) carry no measurement — but that is the SAME kinematic band
    (dR ~ 1/pT), not a statistics failure: only **9.8e-5 %** of the MC denominator weight sits in
    an empty cell. A pair landing there would have been left UNCORRECTED (w_reco = 1), a silent
    one-sided bias, so the evaluator falls back 3D cell -> dR-integrated (pair pT, pair eta)
    cell -> inclusive, each level counted and printed. The 2D fallback level is fully populated
    (0 of 72 empty).
  * Cross-check against an independent measurement on the pair tree: inclusive eps_reco with the
    gap cut on both legs = 0.6624 there vs 0.6620 here (0.06 %). Truth-gap and reco-gap
    definitions of the numerator agree to 4e-4.
  * Loss at the top edges is negligible: truth pair pT > 150 GeV is 3e-6 of the weight,
    |truth pair eta| > 2.4 is 1.1e-4.

- 2026-08-17 Step 4 (eps_dR wired into the crossx, Physics Procedure §3c).
  * `Utilities/DrCorrectionCrossxEvaluator.h` — the exponential -> polynomial -> raw-values
    cascade the user asked for, built ON TOP of the existing `DrCorrectionEvaluator` (which
    already owns the fit_ok gate, the baseline-C screen, the raw-bin fallback, the dR clamp and
    the floor/cap counters). Routing is decided once per (pair pT, pair eta) cell at load time.
    Configuration is read from the named crossx defaults `DrCorrCrossxMethod/Sign/Mode()` =
    expo / os / nocorr_ptmerge; only the BACKUP method name is defined locally.
    A hard binning guard compares the fit file's own axes against
    `ParamsSet::pair_pt_coarse_bins` with the last two bins merged (7 cells) and
    `pair_eta_proj_ranges_coarse_incl_gap` (9 cells), and THROWS on any mismatch.
  * `RDFBasedHistFillingPP.cxx`: the per-pair weight block existed FIVE times verbatim (crossx
    OS, crossx SS, the two no-minv passes, and FillHistogramsGeneric). All five now call ONE new
    method, `AddPairEfficiencyWeightColumns`, which defines
    `eps_dr`, `effcy_pair = eps1*eps2*eps_dr`, `w_trig = 1/effcy_pair`,
    `effcy_reco_pair` (the fullsim pair efficiency) and `w_reco = 1/effcy_reco_pair`.
    `OpenRecoEffPlaceholderFile` is replaced by `OpenPairEfficiencyInputs`.
  * The sibling session's 7-pair-pT-bin fits landed at 23:15-23:22 on 2026-08-17:
    `dr_correction_fits_pp24_full_step3_{expo,polyu_fixedRp}_os_nocorr_ptmerge.root` (+ Medium WP).

- 2026-08-18 Step 5 (USER DECISIONS, asked and answered 2026-08-17/18):
  1. **eps_reco = the FIDUCIAL efficiency, gap cut on BOTH legs** (truth in the denominator, reco
     in the numerator). Inclusive 0.6620 Tight / 0.7301 Medium. The pp24 result is therefore a
     FIDUCIAL cross-section and `eps_acc = 0.9133` stays a separate, not-yet-applied factor.
     The rejected alternative (gap cut on the reco leg only) would have folded eps_acc in,
     giving eps_reco = 0.5569 and a ~19 % higher cross-section, at the price of double counting
     once eps_acc is built.
  2. **The MC-data comparison reads the PLAIN data histograms, not `_wgapcut`.** The legacy
     `_wgapcut` suffix applies `PassSingleMuonGapCut` (|eta| < 0.135 at all pT + the pT<6
     `charge_eta_gap_cuts`), a different and wider object than the analysis fiducial cut. With the
     fiducial cut now in the BASE pp24 selection, the plain histograms already describe exactly
     the region the cross-section measures. Changed in `plot_mc_data_compr.cxx:63`; the sibling
     macro `plot_mc_data_2D_hists_and_1D_proj.cxx` was reading `_gapcut1`, a suffix the RDF hist
     filling has not produced in a long time (it could only have thrown), and now reads plain too.

- 2026-08-18 Step 6: the fiducial gap cut is applied to the GENERIC dataframes as well
  (`RDFBasedHistFillingPP::FillHistogramsGeneric`), which is MANDATORY rather than cosmetic: the
  generic histograms are efficiency-corrected, and a muon at q*eta in [2.30, 2.40) has no fitted
  turn-on on the contiguous coarse binning, so `EvaluateSingleMuonEffcyPtFitted` would THROW.
  (The crossx has not been rerun since 2026-07-08, i.e. since before the coarse top edge moved to
  2.30 on 2026-08-04 -- so this was a latent failure independent of the present change.) It also
  puts the generic / MC-data-comparison histograms in the SAME fiducial region as the crossx.

- 2026-08-18 Step 7: evaluator smoke test (scratchpad `test_evaluators.C`) on the real inputs.
  `DrCorrectionCrossxEvaluator` on the 7x9 = 63 merged cells routes **60 to the exponential fit,
  2 to the polynomial (exponential rejected), 1 to the raw measured bins (both rejected)** --
  exactly the requested cascade. Delivered eps_dR spans 0.061-1.32; values are below 1 almost
  everywhere, as a 2mu4 close-by correction must be (it is a LOSS). One forward-eta cell exceeds
  1 and is flagged by the producer as a probable fit artefact (parent doc R23/R26, OPEN).
  `PairRecoEffEvaluator` returns 0.35-0.80 over the populated cells.

- 2026-08-18 Step 8 — **BUG FOUND AND FIXED IN THE FIRST FULL RUN: a boundary hole between the
  fiducial cut and the trigger-efficiency binning.** The first pp24 crossx run threw:
  `EvaluateSingleMuonEffcyPtFitted: no fitted turn-on for q_eta=2.300 pt=8.10`. Root cause: the
  gap windows were rejected on OPEN intervals (`q*eta > lo && q*eta < hi`) while the coarse q*eta
  turn-on bins are HALF-OPEN (`[lo, hi)`), and the forward window's lower edge 2.30 IS the top
  bin's upper edge. A muon at exactly q*eta = 2.300 therefore survived the cut and had no fitted
  turn-on. The same hole exists at q*eta = 2.40 (the ntuple keeps |eta| <= 2.4, so 2.40 is
  reachable). FIX at the source, `ParamsSet::PassSingleMuFiducialGap` and
  `ParamsSet::FiducialGapCutExpr`: the windows are now rejected CLOSED, `[lo, hi]`, which makes
  the surviving region exactly `[-2.4, 2.30)` minus the two interior windows -- precisely the
  region the turn-on fits cover. The change is measure-zero everywhere else (it moves only exact
  boundary values), so the published gap-cut cost numbers (3.79 % pp muons etc.) are unaffected.
  It DOES touch the shared helper the MC trigger-efficiency chain uses, at the same measure-zero
  level. All dependent classes recompiled and the fullsim reco-eff rebuilt for consistency.

- 2026-08-18 Step 9 — pipeline hardening. `pipelines/pipeline_pp_crossx.sh` reported
  `OK RDF crossx pp24` for the throwing run above: ROOT's TRint CATCHES an exception raised inside
  the RDF event loop, prints it, and the process still exits 0, leaving a freshly-RECREATEd
  759-byte histogram file that opens fine. Stage 5 now additionally requires (a) the output to be
  NEWER than a stamp taken before the run and (b) the output to actually contain
  `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`. Without this the failure surfaced only two
  stages later, in the plotter, and a rerun of a DIFFERENT stage could have gone unnoticed.

- 2026-08-18 Step 10 — **SECOND, DEEPER CAUSE OF THE SAME BOUNDARY BUG: a float/double mismatch.**
  Closing the intervals (Step 8) did NOT fix the throw; the second run failed identically at
  `q_eta=2.300`. Reason: `FiducialGapCutExpr` built a JIT expression comparing `charge*eta`
  (a float, PROMOTED TO DOUBLE by the comparison) against a DOUBLE decimal literal `2.300000`.
  Since `2.30f = 2.2999999523...` is strictly less than the double `2.3 = 2.2999999999...`, the
  edge muon was NOT rejected — while `FindBinReturnStr` compares in FLOAT against the float bin
  edge `2.30f` and found no bin, so the evaluator threw. FIX: the expression is now built entirely
  in float — `(float)(q*eta) >= 2.300000f && (float)(q*eta) <= 2.400000f` — which also makes the
  RDF Filter agree bit-for-bit with `PassSingleMuFiducialGap`, which always compared in float.
  Verified directly: `PassSingleMuFiducialGap(2.3f, +1) = false`, `(2.2999f, +1) = true`,
  `(2.4f, +1) = false`. Everything recompiled; the fullsim RDF and eps_reco were rebuilt and are
  numerically IDENTICAL (inclusive eps_reco 0.661979 Tight / 0.730083 Medium before and after),
  confirming the fix is measure-zero as argued.

- 2026-08-18 Step 11 — **PIPELINE RAN CLEAN, END TO END** (`pipeline_pp_crossx.sh`, SKIP_CONDOR=1,
  exit 0). Delivered:
  * `histograms_real_pairs_pp_2024_2mu4_nominal.root` (689 kB), 705 404 OS signal-region pairs.
  * `plots/single_b_analysis/pp24/` — 5 PNGs; `plots/single_b_analysis/pp24_pt_150/` — 3 PNGs.
  * `plots/sanity_check_crossx/PP_2024_pair_pt_in_eta_subplots.png` (before/after trigger corr).
  * `plots/mc_data_compr/` — 7 PNGs (Dphi, Dphi_unity, DR, DR_unity, DR_jacobian_corrected,
    DR_zoomin, DR_zoomin_jacobian_corrected).
  * eps_dR routing in the production run: **60 cells on the exponential fit, 2 on the polynomial,
    1 on the raw bins** (of 7 x 9 = 63).

### Results & Observations

**Correction-stage decomposition** (integrals of
`h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_corr_*`, pb):

| stage | integral | ratio to previous stage | meaning |
|---|---|---|---|
| `corr_raw` | 1761.68 | — | 1/L only |
| `corr_unfolded` | 1761.68 | 1.000 | w_unfold == 1 (unfolding is a deliberate identity) |
| `corr_unfolded_reco` | 2601.61 | **1.477** | mean 1/eps_reco -> mean eps_reco = 0.677, consistent with the inclusive 0.662 measured in the fullsim |
| `corr_unfolded_reco_trig` | 5784.00 | **2.223** | mean 1/(eps_1 eps_2 eps_dR) -> mean pair trigger efficiency ~0.45, as expected for a 2mu4 AND of two legs |
| nominal (`..._w_signal_cuts`) | 5784.00 | — | equals the last stage, as the invariant in CorrectionStages.h requires |
| `..._no_trig_corr` | 1761.68 | — | equals `corr_raw`, as it must |

**Sanity plot:** the trigger-corrected spectrum sits ABOVE the raw one in every pair-eta bin, by a
smooth factor that shrinks with pair pT and converges at the turn-on plateau. Neither known
pathology is present: no high-pT blow-up (`pp_trig_eff_highpt_jump.md`) and no corrected/raw < 1
(the silent-pair-drop signature).

- 2026-08-18 Step 12 — **/review-analysis-code: PASS-WITH-COMMENTS**, and every finding is now
  closed or explicitly carried forward.
  * **H1 (fixed + ESCALATED).** The reviewer scanned all 63 merged cells and found the cascade
    DELIVERS eps_dR = **2.156** at dR = 0.36 in (pair pT bin 7 [72.1,150), pair eta bin 1
    [-2.4,-2.0)) — the raw-bin-placeholder cell. eps_dR is the probability that a CLOSE-BY pair
    fires 2mu4 relative to two independent legs, so it can only be a LOSS: any delivered value
    above 1 is an artefact, and since w_trig = 1/(e1 e2 eps_dR) it SUPPRESSES those pairs. It was
    invisible because `DrCorrectionEvaluator` prints its span over the cells IT FITTED only — the
    raw branch bypasses that scan, and the Step-7 record here (0.061-1.32) was that fitted range,
    not what the cascade delivers. FIX: `DrCorrectionCrossxEvaluator::ReportDeliveredExtrema`
    now scans EVERY cell over every route at load time and prints the extremum plus an explicit
    `> 1` artefact warning. **Nothing is clipped** — that is the user's call, coupled to the OPEN
    R26 in `mc_trigger_efficiency.md`, whose own hand-off says do not paper over it.
  * **H2 (fixed).** Neither evaluator's `PrintStats()` was ever called, so the whole "counted and
    printed, never silent" contract produced nothing in a production run. Added
    `RDFBasedHistFillingPP::PrintPairEfficiencyStats()`, hooked to `HistPostProcessDataExtra()`
    (after the event loop).
  * **M3 (fixed).** `docs/signal_selection_change_impact.md` was stale in three places and
    actively misleading. Rewritten: §0 now carries the pp24 fiducial region AND the Pb+Pb old
    form side by side with a top banner that the two DIFFER; §1 and §4 record that pp24 no longer
    uses the Run-2 placeholder, so §3.C is result-affecting for pp; §2 gained the three missing
    call sites (`MCTrigEffPairSelection.h`, `plot_dr_vs_pair_pt_diagnostic.cxx`,
    `plot_reco_distr_singleb_vs_op_pp24.C`).
  * **M4 (partly fixed, partly declared).** The low-mass template pass shares the helper, so its
    `crossx_weight_trig_only` now contains eps_dR and its `signal_cuts_no_minv` carries the gap
    cut. It was NOT rerun (template fitting is out of scope), so **the on-disk
    `*_template_fit*.root` outputs are STALE** — recorded here rather than silently changed. Its
    `_nosel` dataframes WOULD have thrown (they were the one un-gap-cut node still feeding the
    trigger weight); they now carry the fiducial cut, with a comment stating that "_nosel" means
    no SIGNAL-REGION cuts, not "no cuts" — the fiducial cut is a prerequisite of the weight.
  * **M5 (fixed where it was a bug, flagged where it is physics).** BUG: the MC-data comparison
    normalised the data by 1/410.815 pb^-1, a value corrected to 400.412 on 2026-06-15 (`009f3b8`)
    and never propagated — every data curve was 2.6 % low against the cross-section it is compared
    with. It now READS `PPBaseClass::GetCrossxFactor(24, "2mu4")`. PHYSICS, flagged not fixed: the
    data is now gap-cut and the truth MC is not, an eps_acc-sized (~16 % of pairs) offset that is
    exactly the factor this task deliberately does not apply.
  * **M6 (documented).** `AddPairEfficiencyWeightColumns` runs on the whole generic df_op/df_ss,
    which has no pair-pT threshold, so eps_reco (measured only for signal-region pairs) is CLAMPED
    into the [8, 11.5) cell below 8 GeV while eps_dR returns 1 outside its grid. The two
    extrapolate differently. Both are now counted and printed by H2's fix.
  * **L7/L8/L9 (fixed).** Counters in both new evaluators are `std::atomic` (Eval runs under
    ImplicitMT); the static evaluator cache is keyed on the working point and throws on a WP
    switch; the crossx binning guard reads `pms.pair_pt_coarse_bins` directly instead of
    `MCTrigEffPairPt::Edges`, so a trigger-efficiency study's `MCTRIGEFF_PAIRPT_4BIN` environment
    variable can no longer move it.
  * Reviewer confirmed clean: no lost or duplicated RDF column across the five refactored sites
    (and the refactor removed a latent dangling `RNode&` into an erased `df_map` entry); the
    cascade routing; the merged pair-pT guard; `Clamp` returning a bin centre; `Project3D("yx")`
    axis order; `Divide(...,"B")` validity.

- 2026-08-18 Step 13 — **/review-plot: PASS-WITH-COMMENTS. No CRITICAL physics failure.**
  * The 2026-07 high-pT blow-up did **not** recur. Per-cell 1/eps_reco spans [1.14, 2.09] and
    1/eps_trig [0.94, 10.6]; the TOTAL correction spans [1.40, 13.7] and is **> 1 in every
    populated cell**. Pair-eta mirror-bin ratios all lie in [0.964, 1.045]. The gap dips are
    present at exactly the cut windows and the right depth (-59 % at |eta| < 0.11, -31 % at
    1.09-1.20, -85 % at 2.29-2.40) — correct for a fiducial result with eps_acc withheld.
  * dsigma/dDeltaR is exactly zero above **0.70**, independently reproducing the kinematic bound.
  * Corrected/raw is a smooth U in pair pT (6.4 -> 2.36 at 17-34 GeV -> 6.2 at 120 GeV) with no
    step or spike. The RISE at high pT is physical, not a blow-up: m < 2.9 GeV forces
    dR ~ 2m/pT, so at 100 GeV the legs sit at dR ~ 0.04 and the 2mu4 two-RoI loss is severe. The
    correction FLATTENS the tail (raw local power-law index 8.06 -> corrected 7.4); the 2026-07
    bug flattened it by x13-88 in a single bin, here the largest single-bin change is x1.25.
  * **P1, the one physics WARNING:** exactly two fine cells have 1/eps_trig < 1 (0.942 and 1.097),
    both at pair eta ~ -2.24 — the same H1 eps_dR > 1 artefact, worth **3 pairs**, 4e-6 of the
    sample. Nothing on any plot moves; carried forward with H1.
  * sigma_fid = **5784 pb** for OS pairs in the fiducial signal region, no background subtraction.
    **RUN2-CROSSCHECK UNVERIFIED** — the KB's Run 2 back-to-back dimuon analysis measures the
    opposite topology and the HF-muon note measures single muons, so neither bounds this. What
    does compare and passes: pair eps_reco 0.662 vs a Run 2 Tight single-muon 0.75-0.80 squared
    (0.56-0.64), and inclusive eps_trig 0.45 for a 2mu4 AND of two near-threshold legs.
  * PRESENTATION, all six fixed and replotted: the sanity plot's legend said the gap between the
    curves was the TRIGGER correction when it is trigger x reco (now "Uncorrected" /
    "Corrected (eps_trig x eps_reco)"); its 5x2 grid for 9 panels is now 3x3 via the shared
    `DetermineSubplotGrid`; its y-title was clipped so the plot displayed the WRONG unit
    ("[pb GeV ']") and now renders "[pb GeV^-1]"; the MC-data 1D distributions had NO log y at
    all, which on the jacobian-corrected dR panels put data and POWHEG flat on zero — all of them
    are log y now, with a positive lower limit taken from the smallest non-zero content; the sign
    label moved from a legend ENTRY (which pushed the box onto the markers) into the legend
    HEADER, and an unset legend position now builds a default-constructed TLegend that ROOT
    AUTO-PLACES into a gap it finds itself (no fixed corner works across a top-peaked dR and a
    U-shaped dphi); and `SingleBCrossxPlotterBase` no longer retypes the reco-efficiency dR edges
    — they moved to `CommonEffcyConfig::dr_bins_edges_for_reco_effcy`, now the single source read
    by BOTH the producer (`RDFBasedHistFillingPythia`) and the plotter.
  * Reviewer's completeness note, accepted: no pp24 plot indicates that the result is FIDUCIAL
    (gap windows removed, eps_acc not applied, no background subtraction), so the visible dips at
    |eta| ~ 0, 1.1, 2.35 are uninterpretable from the figure alone. Given the no-prose-on-plots
    rule this belongs in the note caption; recorded as Remaining Work.

## Remaining Work

1. **USER DECISION — eps_dR > 1 in the forward top-pair-pT cell** (delivered max 2.156). Clip to
   1, route the cell to no correction, or accept. Coupled to `mc_trigger_efficiency.md` OPEN R26.
   Affects 3 data pairs today, but it is in the nominal weight.
2. **eps_acc is not applied** — the pp24 result is a FIDUCIAL cross-section. Building eps_acc
   (0.9133 per muon from truth, `muon_gap_cuts_acceptance.md` F12) is the next step if a
   full-acceptance number is wanted; it also removes the MC-data comparison's data/MC gap-cut
   mismatch (M5).
3. **Pb+Pb has NOT been brought over** — its signal region is still `q*eta < 2.2` and it still
   uses the Run-2 reco-eff placeholder. **R_AA must not be recomputed until it is.**
4. **Stale on disk, deliberately not rerun** (out of scope): the pp24 template-fit outputs
   (`*_template_fit*.root`, `*_no_res_cut`, `*_scrambled` trees), the Pythia/Powheg truth signal
   acceptance, and the MC-closure "data-like" variant (its selection mirror moved with the crossx).
5. **Note/caption** must state that the pp24 spectra are fiducial, eps_acc-free and
   background-unsubtracted (plot-review completeness note).
6. Validate eps_dR below dR = 0.1 before the top three pair-pT points are quoted — they hold
   175 / 60 / 17 raw pairs and rely on the sparsest end of the fit domain.

## Follow-up 2026-08-18 — the MC-data comparison plot set

**Trigger:** the user reported the `plots/mc_data_compr/` PNGs as stale. Diagnosis: only 7 of the
34 were pipeline-produced; the other 27 were ORPHANS from the pre-`0a9eaba` macro (a pp17 series
plus a `_ratio_to_pp17` mode, reading the Run-2 `histograms_real_pairs_pp.root`). Nothing in the
repo could regenerate them, so they had been frozen since Feb 2025 while sitting next to current
output. The root cause of the narrow produced set was UPSTREAM, not in the plotter:
`RDFBasedHistFillingData::BuildFilterToVarListMapDataCommon` filled the generic 1D list for the
`_ss`/`_op` categories with exactly `{Dphi, DR, DR_zoomin}`, so no other 1D data histogram existed.

**Done (user-directed):**
1. **pp17 removed from the code.** The active macro never had it (`s_nDtTypes = 4` = POWHEG bb/cc,
   Pythia, pp24). The four quarantined Run-2 macros under `mc_data_compr/others/` that still drew
   it (normalized by 1/256.8 pb^-1, reading `athena/runMCV2`) are DELETED; nothing referenced them.
2. **Generic 1D list extended** with `Dphi_zoomin`, `Deta_zoomin`, `minv_zoomin` — each already
   defined in `var1D_pp.json` with a fixed binning AND already present under the same name on the
   Pythia side, so no binning is invented. Crossx hist filling rerun; the three new comparison
   plots are produced.
3. **Pair-pT cross-section comparison, 8-150 GeV** (user: "use pair pT 8-150 as default"),
   pair-eta integrated AND in the 9 canonical pair-eta bins (user: "same as
   pp24_crossx_pair_pt_in_eta_subplots.png but with MC"). New macro
   `plotting_codes/mc_data_compr/plot_mc_data_pair_pt_in_eta.cxx` ->
   `pair_pt_mc_data_compr.png` + `pair_pt_in_eta_subplots_mc_data_compr.png`.

**Why the pair-pT plot uses a DIFFERENT Pythia sample from the 1D plots.** The private sample the
1D comparison reads cannot serve it, for four independent reasons: its pair-pT axis is 30 uniform
bins 0-30 GeV (2.8 % of its weight already in the pT overflow, and the whole 30-150 region
absent); 4 of the 9 canonical pair-eta boundaries (+-0.5, +-1.5) fall at BIN CENTRES on its 24-bin
eta axis, so the panels could not be cut without straddling and double-counting; its selection is
`from_same_b` with no signal region; and its weight never touches AMI
(`PythiaAlgCoreT.c:983`, `eventWeight / njobs`). The partner used instead is
`h2d_sig_accept_num_pt_150_eta` from `pythia_truth_full_sample/pythia_5p36TeV/`, whose axes are
BIT-IDENTICAL to the data's `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts` (all 16 pT edges of
`ParamsSet::pT_bins_150`, all 45 pair-eta edges — the macro re-verifies this at load time and
THROWS on divergence). One factor puts it on the data footing: **nb -> pb = 1000**, the AMI
`crossSection` being in nb.

**Result: MC/data = 1.63 integrated** (data 5784.06 pb, MC 9407.74 pb); per pair-eta panel
1.19, 1.33, 1.87, 1.60, 1.92, 1.59, 1.85, 1.35, 1.18. The deliberate mismatches behind that number
are listed in the macro header and must be quoted with it: MC is pure truth single-b while the
data OS carries unsubtracted gluon-splitting and combinatorial background; the MC truth acceptance
is still on `q*eta < 2.2` while the data moved to the gap windows; the MC is the 4-beam
4:6:6:9 Pb-isospin-averaged NN cross-section compared against pp data; and the data is
reco-level-corrected but NOT unfolded.

### OPEN — two normalization/domain issues found on the way, NOT yet decided

**O1. The 1D comparison plots' Pythia normalization is an undocumented empirical scale.**
`plot_mc_data_compr.cxx` applies `1e3 (NB_TO_PB) x 1e3 (TRUTH_COMBINE_RENORM) = 1e6` to the
PRIVATE sample, justified in-code as "the stored MC weight is in nb (AMI crossSection is nb)".
That rationale does not hold for that file — the private weight never touches AMI. Measured, after
each side is treated as `hist_helper` treats it: private Pythia / data = **1.342** on dR, whereas
the AMI-normalized FULL sample with the clean x1000 gives **1.597**, consistent with the 1.626 the
2D cross-section gives independently. So the Pythia curve in the 1D plots sits ~16 % low on a
tuned number with no derivation. **Recommendation: move the 1D comparison onto the same full
AMI-normalized sample the pair-pT plot uses.** Not done — it changes every 1D plot's MC curve.

**O2. eps_reco is extrapolated over most of the GENERIC sample.** The generic dataframes have no
signal-region cut, but eps_reco was measured only inside it. Census from the production run:
**73.6 % of evaluations clamped in dR**, 52.5 % in pair pT, 21.8 % in pair eta; 12.2 % fall back to
the dR-integrated map. So in the 1D comparison plots the away-side pairs (dR ~ 3, which dominate
dR and dphi) are corrected by the dR in [0.6, 1.0) efficiency. The correction is applied as the
user asked, and every clamp is now counted and printed, but the 1D plots' data points are NOT on
the same footing as the signal-region cross-section. Options: measure eps_reco over the generic
domain, restrict the comparison to the signal region, or accept and document.
