# Low-mass dimuon template fit (OS/SS minv 0–4 GeV + Pythia truth templates)

**Mode:** Implementation. **Created:** 2026-06-18. **Session:** "template fit prep".
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.
**Siblings:** `analysis_overview.md` (§2,§4), `reco_eff_placeholder_run2.md` (Q1 dσ
weights), `raa_from_rdf_crossx.md` (OS−SS methodology), `project_overlay_pair_structure`,
`reference_sign_convention` (sign1=SS, sign2=OS).

---

## Objective
Build the ingredients for a template fit of the **low-mass (0–4 GeV) opposite-sign
dimuon mass spectrum** to separate the **single-b signal** from the **gluon-splitting
(open-HF) background** and the **uncorrelated combinatoric** background, using the
same-sign spectrum to constrain the combinatoric and **Pythia evgen truth** minv
templates (split by flavor/origin category) for the correlated components. POWHEG
is intentionally excluded for now (NLO template biased against g→QQ̄ — see
`project_mc_sample_roles`).

First concrete deliverable (this cycle): add **efficiency-corrected,
crossx-normalized OS and SS minv histograms over 0–4 GeV with the dimuon
mass-window cut REMOVED** (other single-b kinematic cuts kept), for pp and PbPb.
These are the data inputs `D_OS(m)`, `D_SS(m)` of the fit.

**Umbrella scope (added 2026-06-21).** This doc is the authoritative plan for the
**full background-subtraction program that REPLACES the provisional `OS − SS`**
currently used in crossx and R_AA (`raa_from_rdf_crossx.md` task_06). The minv
template fit here extracts the signal *fraction/yield* per bin; that yield is then
**acceptance-corrected** and fed into the existing crossx/R_AA normalization. The
program has three physics stages that must be ordered correctly: (1) Δp/p
fake-muon removal; (2) origin-blind detector corrections (efficiency + unfolding)
applied to the *whole mixture* BEFORE the fit; (3) the minv template fit, then
(4) signal acceptance applied to the signal yield AFTER the fit. See Physics
Procedure §3e–§3g and §4.

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation
The single-b signal (both muons from one b-hadron: direct b→μ + cascade b→c→μ) is
opposite-sign, collinear, and low-mass. Under it sit two backgrounds: (a)
uncorrelated **combinatoric** pairs (charge-symmetric), and (b) **gluon-splitting /
open-HF** pairs (g→bb̄, g→cc̄: two muons from two different HF hadrons), also
opposite-sign and low-mass. Same-sign subtraction removes (a) but NOT (b), so the
residual correlated HF background must be removed with an MC-shape template fit.

### 2. Top-level model (per minv bin, fixed single-b kinematic region, 0–4 GeV)
```
OS(m) = S(m) + C(m) + G(m)
SS(m) = C(m) + k · G(m)
```
- `S(m)` — single-b signal (truth flavor category `_single_b` / `from_same_b`).
- `C(m)` — uncorrelated combinatoric, **charge-symmetric**: identical in OS and SS.
  In the fit, `C(m) = N_C · T_mix(m)` where `T_mix` is the unit-normalized
  **mixed-event** shape and **N_C is a free normalization parameter** (the
  mixed-event template has an arbitrary generation normalization, so it MUST be
  scaled — answering the user's "factor on the mixed-event template": yes, that
  factor IS N_C, and it lives inside C). Nominal: the SAME `N_C·T_mix` in OS and SS
  (charge symmetry); if the 5b closure shows a charge asymmetry, free a separate SS
  combinatoric scale `N_C^SS = (1+ε)·N_C` as a systematic.
- `G(m)` — opposite-sign gluon-splitting / open-HF background (truth `_bb`+`_cc`
  +`_one_b_one_c`; origin distinguishes flavor-creation `_FC` from gluon-splitting `_gs_*`/`_GS_*`).
- `k` (0<k<1) — fraction of OS gluon splitting that appears as SS. Only **g→bb̄**
  contributes via either **B⁰–B̄⁰ oscillation** (one neutral-B mixes, flipping a
  muon sign) or **direct+cascade pairing** (one b→μ direct, the partner b̄→c̄→μ →
  same sign); g→cc̄ has no charge-flip mechanism. (NB: cascade+cascade is OS, not
  SS.) So `k` ≈ (osc + direct·cascade prob)·(bb̄ fraction of G).
  **Why k is worth using (robustness):** k is a *ratio* of SS to OS gluon-splitting
  drawn from the *same* g→QQ̄ population, so the absolute (poorly-known, large-
  theory-uncertainty — KB `gluon_splitting_flavour_excitation`) g→QQ̄ cross-section
  largely cancels. k is driven by *measured* quantities (semileptonic BRs, χ_d) plus
  one theory-sensitive piece — the **cc̄:bb̄ ratio in OS** (cc̄ enters OS but not SS).
  This is why anchoring the physics-background normalization via k + SS (data) beats
  taking the G normalization from Pythia directly. Mitigate the cc̄:bb̄ sensitivity:
  cc̄ and bb̄ have different minv shapes (the OS fit partly constrains the mix), and
  bound the rest as a systematic. **k may depend on (pair pT, η, minv)** — measure
  the dependence in MC, do not assume constant.

  **k is EMPIRICAL — what it absorbs vs what can break the model (refined 2026-06-22).**
  k is *measured* as the total correlated-background ratio G_SS/G_OS in MC, so it
  **automatically includes** (i) **flavour excitation (FE)** — ISR g→QQ̄ has the SAME
  QQ̄ charge structure as FSR gluon-splitting (two HF hadrons → identical
  direct/cascade/mixing combinatorics), so FE lives inside G and inside the measured
  k; and (ii) **ALL bb̄ charge-flip mechanisms** — direct+cascade, single AND double
  B⁰ mixing, b-baryons. None of these are "neglected" or separately modeled; they are
  in k by construction. Therefore `SS = C + k·G_OS` can fail ONLY from:
  - (a) **SS sources outside {correlated-HF, combinatoric}** — residual jet/fake muons
    (π/K, punch-through surviving the Δp/p fit) and **>2 HF muons in a jet** (extra
    splittings → wrong-charge pairings not following the simple k). The >2-HF-μ term is
    expected **rare**; if the data closure (5b) shows the OS/SS factor roughly holds, it
    need not be modeled — **carry it as a systematic** (per user 2026-06-22). Residual
    fakes are the part the Δp/p fit must control upstream.
  - (b) **MC mismodeling of the ratio** (cc̄:bb̄ fraction, FE:GS fraction, mixing) — k_data ≠ k_MC. This is the systematic, bounded by the data closure test (Step 5).
  - (c) **combinatoric charge asymmetry** (C_OS ≠ C_SS) — must be checked.
  - (d) **k kinematic dependence** not captured if a single k is used.
  These four — not cascade/mixing/FE — are what Step 5 validation must rule out.

Subtraction identity (consequence): `OS(m) − SS(m) = S(m) + (1−k)·G(m)` — combinatoric
cancels exactly; a fraction (1−k) of gluon splitting survives and is removed by the template.

### 3. Step-by-step method

**The two fits live at two DIFFERENT kinematic levels (resolved 2026-06-22, per user):**
- **Δp/p fake-muon fit → RECO level.** Δp/p = (p_ID − p_MS)/p_ID is the ID-vs-MS
  momentum imbalance — an **intrinsically reconstructed** quantity (no truth analogue),
  so the fake/real separation is done at reco level, BEFORE unfolding.
- **Combinatoric + physics (minv) fit → TRUTH level.** The Pythia templates (signal,
  G, resonances) are **truth-level** objects, so the minv template fit is performed on
  the **unfolded** data (after efficiency + unfolding bring the data to truth level).
  This settles the earlier "fit at reco vs truth" question: minv fit at truth.

a. **Data spectra (this cycle):** fill OS and SS `minv` over **0–4 GeV** in the
   single-b kinematic region **without the mass window** — selection
   `pair_pt>8 && m1.charge*m1.eta<2.2 && m2.charge*m2.eta<2.2`
   (= `signal_cuts` MINUS `minv>1.08 && minv<2.9`; **the ΔR>0.05 cut is REMOVED** —
   per user 2026-06-22, see Design Decisions). Weight = the existing
   crossx differential-cross-section weight, reco+trig corrected:
   - pp: `crossx_weight_trig_corr` (= weight·(1/L_pp)·w_reco·w_trig).
   - PbPb: `weight_for_dsigma_trig_corr` (= weight·(1/L_year)·w_reco·w_trig), per
     centrality, lumi-combined at plot time (ΣN/ΣL). Binning **50 bins, 0–4 GeV**
     to match the Pythia truth `minv_zoomin` template exactly.
b. **Templates:** Pythia evgen truth `minv_zoomin` (0–4 GeV) per flavor category
   (`_single_b` = S; `_bb`,`_cc`,`_one_b_one_c` = G; `_resonance(_contaminated)` =
   J/ψ etc.) and per origin category. Truth = no reco, so compared to
   efficiency-corrected data (data unfolded back to truth within acceptance).
c. **Baby-step fit (later cycle):** assume k≈0 → form `D_OS−D_SS` (removes
   combinatoric), fit `= N_S·T_S + N_G·T_G` (TFractionFitter or binned χ²) over
   0–4 GeV, vetoing the J/ψ window [2.95,3.25]. Extract N_S.
d. **Refinement (later):** restore k; coupled OS+SS fit with shared C,G and a
   combinatoric template (event-mixing); float/constrain k.
e. **Correction ordering (origin-blind BEFORE the fit).** The reco+trigger
   efficiencies and the detector response (bin migration / unfolding) are
   properties of the ATLAS muon detector/trigger/offline reconstruction; given the
   pair kinematics they are **independent of pair origin** (the detector cannot
   know signal vs background). They are applied to the **whole mixture before the
   fit**. The PRECISE internal order matters because each correction lives in a
   different kinematic frame (refined 2026-06-22, per user):
   1. **Trigger efficiency** — evaluated as a function of **RECONSTRUCTED**
      kinematics (the trigger fires on reco objects). Apply as a per-pair weight at
      reco level, fill the reco spectrum.
   2. **Unfolding** (reco→truth, bin migration) — a **spectrum-level** operation
      (response matrix / iterative Bayes), NOT a per-pair multiplicative weight.
      Applied to the trigger-corrected reco spectrum.
   3. **Reconstruction efficiency** — evaluated as a function of **TRUTH**
      kinematics (the proper 3D pair ε_reco(truth pair pT, η, ΔR) is defined per
      truth pair), so it is applied **after** unfolding, on the truth spectrum.
   Then the template fit runs at truth level (sharpest signal mass peak; templates
   need no detector folding). **Run 2 precedent (HF-muon R_AA, arXiv:2109.00411 /
   ANA-HION-2019-58 §4.1–4.2):** the ρ and d0 template fits are run on
   **efficiency-weighted** distributions — efficiency FIRST, then fit → N_corr
   directly. (The dimuon note is no counter-precedent: its Δp/p was a purity
   demonstration with no yield subtraction.) **Current-code status (RDF, verified
   2026-06-22):** trigger weight ✓ at reco; reco-eff applied at **reco** kinematics
   (placeholder ε₁·ε₂, NOT truth); unfolding ABSENT (`w_unfold ≡ 1.0` identity).
   The order is therefore not yet implemented — but it is **degenerate now** (both
   weights are reco-level placeholders multiplied together, no unfolding). The
   prescribed trig(reco)→unfold→reco(truth) order is a **structural change** to make
   when (i) real unfolding and (ii) the truth-binned 3D pair ε_reco land (roadmap
   Q4 / task_05). These corrections are NOT applied after the fit.
f. **Signal acceptance (origin-specific, AFTER the fit).** A_sig(pair pT, pair η)
   = P(a truth single-b dimuon passes the signal selection), from signal MC truth
   (Pythia/Powheg). It is signal-only, so it is applied to the **extracted signal
   yield only**, after the fit: `N_sig(X) = N_S^fit(X) / A_sig(X)`. Build it on a
   **fine 2D (pair pT, pair η) grid** (or per-event MC weight) — MC stats are not
   the bottleneck — and project to the 1D observable. Default: **fiducial**
   measurement (A_sig corrects only cuts not in the binning — the mass window, ΔR,
   single-μ thresholds — and resolution at the pair-pT>8 boundary), no full
   phase-space extrapolation unless decided otherwise.
   - **Binning (refined 2026-06-22, per user):** use **log pair-pT binning** — even
     with Pythia pT-hat slices the pair-pT count spectrum falls as a power law, so
     log bins give roughly uniform stats per bin. Pair-η **coarse** (pair-η
     resolution for this study is low). Before trusting the map, **plot
     A_sig(pair pT) projected in a few pair-η bins** to confirm per-bin statistics
     and smoothness. Expect to **fit/interpolate the pair-pT dependence** with a
     smooth function and **Eval at the exact pair pT** rather than using raw bin
     values (continuous-function convention — see memory `feedback_use_continuous_fit_function`),
     so statistical bin-to-bin fluctuations do not propagate into the correction.
   - **What A_sig covers vs what reco-eff covers (clarified 2026-06-22, per user
     Q11):** the single-muon **kinematic/fiducial** cuts — **pT>4, |η|<2.4** — are
     ACCEPTANCE (truth-level: did the truth muons land in the fiducial region). The
     **reco-eff** is conditional on that: P(reconstruct+ID a muon ALREADY in the
     fiducial region), and it also carries the **resolution smearing across the pT>4 /
     η edges** (in/out migration). The **Δp/p / quality / ID** cuts are reco-level
     selections → they belong to **reco-eff / fake-removal, NOT acceptance** (no truth
     analogue). **CRITICAL caveat — verify the pair reco-eff denominator** (fullsim):
     if its denominator is "truth pairs already inside the single-μ fiducial
     (pT>4,|η|<2.4)", then A_sig and ε_reco factorize cleanly and A_sig owns the pT>4/
     |η|<2.4 cuts; if the denominator is "ALL truth single-b pairs", then ε_reco
     ALREADY contains the single-μ kinematic acceptance → A_sig must NOT re-apply it
     (double-count). Settle this before building A_sig.
g. **Feed crossx & R_AA (replace OS−SS).** The acceptance-corrected N_sig(X)
   replaces the `OS − SS` count that `raa_from_rdf_crossx.md` §3 currently feeds
   into the normalization (pp: 1/L; PbPb: crossx_factor carrying 1/N_evt·1/⟨T_AA⟩).
   No OS−SS anywhere in the final chain.
h. **Resonance leakage into the signal window (OS-only — added 2026-06-22, per
   user).** The OS muon-pair trees are resonance-VETOED at **reco** minv in windows
   [0,1.06], [2.9,3.3] (J/ψ), [3.55,3.8] (ψ′) (see Results & Observations; veto is
   OS-only — SS has no resonances). Reco smearing widens the narrow truth peaks, so
   a φ (1.019) or J/ψ (3.097) can be reconstructed *outside* its veto window and
   **leak into the signal mass window** — an OS-only contamination.
   - **Unfolding does NOT remove it.** The response matrix is built from *signal*
     MC, which contains no resonances; unfolding migrates a leaked-resonance reco
     count as if it were signal → contaminates the truth signal spectrum. Leakage
     must be handled in the fit, not by unfolding.
   - **Effect on the combined OS+SS fit:** the leaked resonance tail is OS-only,
     exactly like signal (SS has none) → without modeling it, the fit attributes it
     to signal and **biases N_S up**.
   - **Treatment (DECIDED 2026-06-22, per user): build OS-only resonance mass
     templates that ACCOUNT FOR THE SMEARING and include them in the minv template
     fit.** Even though the fit is performed after unfolding, the leaked-resonance
     contamination is a reco-smearing effect that the signal-only unfolding does not
     remove, so a dedicated resonance component is added to the OS fit model:
     `OS(m) = S + C + G_OS + Σ_r R_r(m)`, with `R_r` the φ/J/ψ/ψ′ templates (SS has
     none). Build `R_r` from a **no-veto OS pass** (`trigger_effcy_calc=true` skips
     the resonance cut — R&O): take the **reconstructed (smeared)** resonance peak
     shape and carry it through the SAME signal unfolding as the data, so `R_r`
     represents the residual the resonance leaves in the fit space. Each `R_r`
     normalization can be constrained from the peak region outside the signal window.
   - **Signal-window upper edge = 2.9 (CONFIRMED 2026-06-22):** the window is
     **1.08–2.9** (the earlier "2.6" was a user typo). With the J/ψ veto starting at
     2.9 there is no buffer, so the smeared-J/ψ template above is REQUIRED, not
     optional — its tail leaks directly below 2.9 into the signal window.

**Granularity (applies to the whole program):** the template fits (§3c/§3d and the
Δp/p fake fit) are **data-statistics-limited → coarse, 1D in the plotted
observable** (per pT bin for R_AA vs pT; per η bin in coarse pT slices for R_AA vs
η) — NOT a fine 2D (pT,η) fit grid. The MC corrections (efficiency, acceptance)
are **fine** (MC-limited, not data-limited). Do not couple their granularities.

**Δp/p fake-muon removal (precedes everything).** Fakes (π/K decay-in-flight,
mis-ID) are removed first via the Δp/p significance template fit (task_07 / roadmap
step 16), because the muon efficiency is defined for real muons. The Run 2 dimuon
note found >98% purity ⇒ likely a flat/coarse purity factor or demonstration-only.

### 4. Negative constraints
- The signal selection is `pair_pt>8 && per-muon q·η<2.2` ONLY (mass window applied
  separately; **NO ΔR cut** — ΔR>0.05 removed 2026-06-22). The minv histograms drop
  the minv-window part; OS and SS use this identical selection.
- **Do NOT apply any ΔR cut in signal selections** (data spectra, truth templates,
  acceptance, mixed-event). ΔR is still a reco-eff binning axis, NOT a selection cut.
- OS and SS MUST use the SAME selection, binning, and efficiency+lumi weight, so
  OS−SS is a subtraction of identically-normalized spectra.
- MUST NOT reuse the T_AA-weighted `weight_for_RAA_trig_corr` for the spectra — the
  cross-section normalization is the 1/L `dsigma` weight (pp `crossx_weight_trig_corr`).
- Template (truth) selection MUST match the data kinematic region
  (truth pair_pt>8, per-muon truth q·η<2.2; **NO ΔR cut**); the ORIGINAL truth
  per-category minv fills are INCLUSIVE (no kinematic cut) — see Results & Observations.
- k is bb̄-only; do not model g→cc̄ as contributing to SS beyond the lumped approximation.
- **Do NOT apply efficiency/unfolding after the fit** — they are origin-blind, so
  they go first, on the whole mixture (§3e). **Do NOT apply signal acceptance to the
  mixture** — it is signal-only, applied to the fitted signal yield AFTER the fit
  (§3f); applying it to a contaminated yield corrects background with a signal factor.
- **Do NOT leave OS−SS as the final background subtraction** in crossx/R_AA. OS−SS
  cancels combinatoric but leaves `(1−k)·G` (correlated physics bkg); it is retained
  only as the provisional first-look (task_06) until this program replaces it (§3g).
- **Combinatoric C comes from the mixed-event method**, not from SS. SS's role is the
  G-normalization anchor via k (§2), not the combinatoric estimate.
- **Do NOT take the G (gluon-splitting/open-HF) normalization from Pythia as nominal**
  once k is validated — use the SS-anchored normalization; keep pure-MC normalization
  as a systematic cross-check.
- k and the combinatoric charge symmetry (C_OS ≈ C_SS) MUST be validated in MC before
  the coupled OS+SS fit is trusted; if k fails, fall back to pure-MC templates and
  carry the full g→QQ̄ theory uncertainty.

## Context
- pp crossx fill: `RDFBasedHistFillingPP.cxx` — `signal_cuts` (L382), OS weighted
  node `df_single_b_crossx_weighted` (L421+, `crossx_weight_trig_corr` L433), SS block
  (L476–513). 1D minv currently only as 2D `h2d_crossx_pair_pt_minv_w_signal_cuts` (L515).
- PbPb crossx fill: `RDFBasedHistFillingPbPb.cxx` — `signal_cuts` (L920), OS weighted
  node `df_single_b_crossx_weighted` (L962+, `weight_for_dsigma_trig_corr` L983), SS
  block (L1011–1044), per-centrality loop (L1048+), dσ minv as 2D (L1083).
- Pythia truth templates: `RDFBasedHistFillingPythiaTruth.cxx` — `FillHistogramsFlavorBinned`
  (L254), `FillHistogramsOriginBinned` (L281); truth `minv_zoomin` 0–4 GeV (50 bins,
  `var1D_pythia_truth.json`) filled as 2D `{pair_pt, minv_zoomin}` per category, OS
  (`_sign2`/df_op_weighted) and SS (`_sign1`/df_ss_weighted). Category maps:
  `MuonObjectsParamsAndHelpers/muon_pair_enums_MC_utils.h` (flavor `_single_b/_bb/_cc/...`,
  origin `_FC/_gs_FSR/_GS_ISR_*/...`).

## Scope
In (this cycle): OS+SS 0–4 GeV no-mass-cut minv data histos (pp + PbPb per centrality),
eff+lumi (dσ) weighted. Out (later cycles): truth templates with matching kinematic
selection; the fitter; combinatoric (event-mixing) template; k determination; plots.

## Design Decisions
- **Binning 50×[0,4] GeV** to match Pythia truth `minv_zoomin` for apples-to-apples fit.
- **Self-contained new branches** off `df_op`/`df_ss` filtered with `signal_cuts_no_minv`,
  replicating the existing trig+reco+dσ weight chain (mirrors the existing self-contained
  SS block), rather than restructuring the shared weighted node — lower risk, no change
  to existing histograms. Reviewer may suggest factoring the weight chain into a helper.
- **Names:** `h1d_crossx_minv_0_4_{op,ss}_dsigma` (pp); `..._{op,ss}_dsigma_<ctr>` (PbPb).
  "0_4" = full 0–4 GeV spectrum (no mass window); "dsigma" = 1/L cross-section weight.
- **ΔR>0.05 REMOVED from the signal selection (2026-06-22, user).** Old: signal
  required ΔR(μ,μ)>0.05. New: no ΔR cut anywhere in the signal selection (data, truth
  templates, acceptance, mixed-event). Reason (user decision): keep the full collinear
  single-b regime. ΔR remains a **reco-eff binning axis**, not a selection. Consequence:
  data histos (Steps 2/3) and truth templates (4a), both filled with ΔR>0.05, must be
  regenerated; `analysis_overview` §2 and the RDF `signal_cuts` must drop ΔR too.
  Implication to watch: very-low-ΔR pairs (near-merged muons) now enter — the reco-eff
  ΔR binning must cover ΔR→0 reliably.
- **Reco-eff / acceptance boundary VERIFIED — no double-count, no reco-eff change
  (2026-06-22, investigation).** pT>4/|η|<2.4 is NOT folded into either reco-eff value:
  (a) fullsim 3D pair ε_reco — the truth fiducial (`PythiaAlgCoreT::PassCuts_PythiaCore`,
  truth pT>4 & |truth η|<2.4) is in BOTH numerator and denominator → cancels → ε_reco is
  conditional; (b) placeholder `EvaluateSingleMuonRecoEffPlaceholder` clamps pT to [4,19]
  as a range floor (not a cut), η only for barrel/endcap. ⇒ **acceptance owns
  pT>4/|η|<2.4; do NOT change reco-eff.** CRITICAL for the acceptance build (Step 7): the
  A_sig DENOMINATOR must come from the `h_cutAcceptance` cutflow (all truth pairs in the
  `nocut` bin, filled BEFORE `PassCuts`), NOT the pair tree (already fiducial-restricted
  → A≈1). Edge migration across pT=4/|η|=2.4 is owned by **unfolding** (detec_resp), not
  reco-eff or acceptance.
- **ScrambGen revival approach (2026-06-22, plan):** rewrite ScrambGen to the modern
  object model (read `MuonObj`, WRITE `muon_pair_tree_sign1/sign2` with `MuonPairObj`) so
  the **existing RDF `signal_cuts` fills the mixed-event minv template** — identical
  selection to the data D_OS/D_SS. ScrambGen becomes a pure pair generator: DROP its
  resonance/photoproduction vetoes + dR bucketing (cuts move to RDF). Defaults adopted:
  **20 centrality intervals** (0–100%, `ParamsSet::nCtrIntvls`), **per-year** PbPb mixing
  (combine downstream), **×5 oversampling** for `nScramb` (template is shape-only; N_C
  floats). Same-centrality (≤5%) mixing kept. Requires a single-muon-tree production pass
  first (stale/missing for the May 2026 skim). ScrambGen currently does NOT compile (uses
  the retired `class Muon`); the rewrite fixes that and ensures consistency.

## Implementation Plan
1. Tracking doc + Physics Procedure (this file). DONE.
2. Add OS+SS 0–4 GeV no-mass-cut minv histos to PP + PbPb crossx RDF (per §3a,§4).
   → `/review-analysis-code` (quote §2,§3a,§4). ACLiC compile. **DONE** (2026-06-21,
   reviewer PASS iter 1).
3. Recompile + rerun crossx RDF (pp + 3 PbPb, skip ntuple); sanity-check the new
   histos (non-empty, OS≥SS in continuum). **DONE** (2026-06-21, all 4 years).
4a. Truth templates: 1D `minv_zoomin` (0–4 GeV) per flavor + origin category, OS+SS,
   with the MATCHING truth single-b kinematic selection (`truth_pair_pt>8 &&
   m1.truth_charge*m1.truth_eta<2.2 && m2.truth_charge*m2.truth_eta<2.2`,
   **NO ΔR cut**, no minv cut). New method in `RDFBasedHistFillingPythiaTruth`
   mirroring `FillHistogramsFlavorBinned/OriginBinned`, suffix `_sigsel`. (per §3b,§4)
   → `/review-analysis-code` (quote §3b,§4). **DONE & REGENERATED WITHOUT ΔR (2026-06-22).**
   Method `FillHistogramsTemplateMinvSignalRegion` added; the ΔR-cut-removal sweep
   (`remove_dr_cut_signal_selection.md`, `/review-analysis-code`+`/review-plot` PASS)
   dropped `truth_dr>0.05` from its `kin_cuts` and refilled pythia truth (output
   `pythia_5p36TeV/histograms_pythia_5p36TeV_no_data_resonance_cuts.root`, 01:28; backup
   `hist_backup_20260622_pre_dr_cut_removal/`). Verified: 32 `_sigsel` histos
   (2 signs × {9 flavor + 7 origin}); OS `single_b` integral 12.54 / 2.23M ent, OS `bb`
   1.32 / 1.40M ent, 50 bins 0–4 GeV. Step 2/3 data histos likewise refilled without ΔR
   by the same sweep → templates and data spectra are selection-consistent.
4b. (later) Fitter; plots.
5. (later) **k-ratio validation study — the GATE for the combined fit** (per §2
   "k is EMPIRICAL"). TWO levels, both REQUIRED; save all plots as the
   justification (or, if it fails, the explanation to colleagues for falling back
   to pure-MC templates):
   - **(5a) MC truth-level k.** Using Pythia truth ORIGIN/FLAVOR labels, compute
     k(bin) = G_SS(bin)/G_OS(bin) from the truth-G categories (`_bb`,`_cc`,
     `_one_b_one_c`; FSR gs + ISR FE together) over the fit minv window, **as a
     function of pair pT and pair η** (and minv). Deliverable: k maps + projections
     → is k stable, or does it need (pT,η) dependence? Decompose k into the robust
     part (BRs, χ_d) vs the theory-sensitive cc̄:bb̄ part for the systematic.
     **Sufficiency note (answering the user):** the origin-categorized 0–4 GeV minv
     distributions ARE sufficient to measure k itself — but NOT to validate the full
     `SS=C+kG` model, because they cannot see the jet/fake/>2-HF terms (§2 a) or the
     real combinatoric. So 5a is necessary but not sufficient; 5b is required.
   - **(5b) Data-level closure.** Predict SS_data ?= C_mixed-event + k·G_OS,MC and
     overlay vs the ACTUAL SS data, **per (pair pT, pair η) bin**. Agreement across
     all bins ⇒ the assumption holds (no significant jet/fake/>2-HF SS excess, k_MC
     ≈ k_data). Disagreement localizes the failure. This is the test of the user's
     real concern (the non-HF SS terms).
   - Also verify **C_OS ≈ C_SS** (combinatoric charge symmetry) from the mixed-event
     template.
   → `/review-investigation` then `/review-plot` (quote §2, §3c/§3d, §4).
   **Raw material:** 5a needs the truth origin/flavor minv ALSO binned in (pair pT,
   pair η) — Step 4a produced 1D minv per category with the matching selection but
   may be pT/η-integrated; CHECK and, if so, add a (pT,η)-binned truth fill first.
   5b needs the mixed-event template (Step 6) → so Step 6 partly precedes 5b.
5'. (later) **Resonance-leakage quantification (§3h).** No-veto OS pass
   (`trigger_effcy_calc=true`); fit φ/J/ψ/ψ′ peaks; extrapolate tails into the
   signal window; report tail/signal ratio. Decide: add OS-only resonance template
   vs negligible-with-systematic. Settle the 2.6-vs-2.9 signal-window edge first.
   → `/review-investigation`, `/review-plot` (quote §3h). **Raw material:** a
   no-resonance-veto ntuple/RDF pass (exists via `trigger_effcy_calc=true`).
6. (later) **Coupled OS+SS fit + event-mixing combinatoric**, per coarse RAA bin →
   N_S^fit(X). **Mixed-event template is NOT ready** (Step 6 prereq, verified
   2026-06-22): `ScrambGen/` exists (PbPb + PP) and does same-centrality (≤5%)
   mixing, but (i) its resonance veto removes exactly our low-mass fit region —
   must switch to our signal selection (pair_pt>8, per-muon q·η<2.2, dR>0.05) with
   NO/looser mass window; (ii) it reads FLAT single-muon branches while the current
   producer writes a `MuonObj` object branch — reconcile before re-running; (iii)
   existing outputs are the 2023 skim (`scrambled_muon_pairs_*.root`, local GPFS) —
   must re-run on the May 2026 skim; (iv) regenerate the hard-coded `nScramb` target
   counts. ⇒ revive+modify+rerun ScrambGen before 5b/6.
   → `/review-analysis-code`, `/review-plot` (quote §3c/§3d, §4).
7. (later) **Signal acceptance** 2D MC map; `N_sig = N_S^fit / A_sig` after the fit.
   → `/review-analysis-code`, `/review-plot` (quote §3f).
8. (later) **Wire into crossx & R_AA**, replacing OS−SS (`raa_from_rdf_crossx.md`).
   → `/review-analysis-code`, `/review-plot` (quote §3g, §4).
9. (later) **Systematics** — k uncertainty, cc̄:bb̄, template shapes, mixed-event
   normalization, fit model, fiducial vs extrapolated. → `/review-analysis-code`.

## Progress Log
- 2026-06-18 — Step 1: doc created. Established OS/SS background model (§2), grounded
  truth categories (`_single_b/_bb/_cc`, origin `_FC/_gs_*`). Confirmed existing dσ
  weights (`crossx_weight_trig_corr` pp, `weight_for_dsigma_trig_corr` PbPb) are the
  correct crossx normalization. Found existing truth per-category minv fills are
  KINEMATICALLY INCLUSIVE (no pair_pt/q·η/dr cut) and stored 2D — see R&O.

- 2026-06-21 — **Steps 2+3 DONE** (`/review-analysis-code` PASS iter 1; log
  `.claude/logs/review-analysis-code-20260621-235115-low-mass-minv-os-ss.md`). Added the
  no-mass-cut OS+SS minv histos to PP (`RDFBasedHistFillingPP.cxx`, end of
  `FillHistogramsCrossx`: `h1d_crossx_minv_0_4_{op,ss}_dsigma`, weight
  `crossx_weight_trig_corr`) and PbPb (`RDFBasedHistFillingPbPb.cxx`, per centrality:
  `h1d_crossx_minv_0_4_{op,ss}_dsigma_<ctr>`, weight `weight_for_dsigma_trig_corr`),
  50 bins 0–4 GeV, selection `signal_cuts` minus the minv window, OS=df_op/SS=df_ss
  identical. ACLiC-clean (separate sessions). Reran crossx pp24 + pbpb23/24/25 (read
  existing ntuples): pp 2 new 1D histos, each PbPb year 12 (6 ctr × OS+SS), all
  non-empty. Sanity: PP OS integral 6209 ≥ SS 990.6 (OS≥SS all continuum bins, 0
  violations); PbPb ctr0_5 SS/OS≈0.72 (large central combinatoric, expected).

- 2026-06-21 — **Scope extended to the full background-subtraction program**
  (planning only, no code). Added: correction-ordering as authoritative procedure
  (origin-blind eff+unfolding BEFORE the fit; signal acceptance AFTER — §3e/§3f);
  R_AA integration replacing OS−SS (§3g); k-ratio robustness via the cc̄:bb̄
  decomposition (§2); granularity split (coarse data fits, fine MC corrections);
  Δp/p fake removal as the first stage; extended negative constraints (§4) and plan
  steps 5–9. Grounded in KB `concepts/muon_source_template_fits`,
  `physics/background/gluon_splitting_flavour_excitation`. A duplicate doc
  (`template_fit_background_subtraction.md`) drafted this session was deleted in
  favor of extending this one (no-duplication rule). **Awaiting explicit user
  approval before any implementation of steps 5–9; steps 2–4 remain as previously
  planned.** Roadmap step 16 updated to point here (OS−SS marked provisional).

- 2026-06-22 — **Physics Procedure refined from 9-point user review (planning only).**
  (1) §3e efficiency ORDER made precise: trigger(reco) → unfold(spectrum-level) →
  reco-eff(truth), then fit at truth; confirmed by Run-2 HF-muon precedent
  (efficiency-weight FIRST, then ρ/d0 fit → N_corr; ANA-HION-2019-58 §4.1–4.2) and
  by a code check (RDF currently: trig ✓ reco-level, reco-eff at RECO not truth,
  unfolding `w_unfold≡1.0` absent → order degenerate now, structural change due when
  real unfolding + 3D truth-binned pair ε_reco land, task_05/Q4). (2) §2 reframed:
  **k is EMPIRICAL** → FE (ISR g→QQ̄, same charge structure as FSR GS) and ALL bb̄
  charge-flip mechanisms (direct+cascade, single+double mixing, b-baryons) are IN k,
  not neglected; the model breaks only from non-HF SS (jet/fake, >2 HF μ), k_MC≠k_data
  mismodeling, C_OS≠C_SS, or uncaptured k(pT,η). (3) New §3h: OS-only **resonance
  leakage** (φ/J/ψ smeared past the reco veto into the signal window) — unfolding
  does NOT remove it (signal-only response matrix), biases N_S up in the combined
  fit; treat via OS-only resonance template or prove negligible; 2.6-vs-2.9 signal
  edge flagged. (4) §3f acceptance: log pair-pT binning + projection check +
  fit/interpolate. (5) Step 5 split into 5a (MC truth k) + 5b (data closure) — origin-
  categorized minv alone validates k but NOT the full SS=C+kG (can't see jet/fake);
  new Step 5' resonance quantification. (6) Step 6: mixed-event ScrambGen NOT ready
  (stale, wrong mass veto, flat-vs-MuonObj input mismatch, 2023 outputs) → revive+
  modify+rerun. KB `gluon_splitting_flavour_excitation` synced (empirical-k framing).
  **No code; awaiting user approval.**

- 2026-06-22 (2) — **Second 9-point user round applied (planning + 1 code comment).**
  (1) RDF `w_unfold` placeholder comment in PP+PbPb now states unfolding is a
  spectrum-level structural change, not a weight (code comment only, no behavior
  change). (2) §3 lead note: Δp/p fake fit at RECO (intrinsic ID−MS), minv physics+comb
  fit at TRUTH (Pythia templates) — settles the fit-level question. (3) >2-HF-μ reframed
  as rare → systematic if k holds (§2a). (4) §3h resonance: DECIDED to build smeared
  OS-only resonance templates and include them in the fit (`OS=S+C+G+ΣR_r`); window
  CONFIRMED 1.08–2.9 (2.6 was a typo). (5) §2 C: clarified `C=N_C·T_mix`, N_C a free
  fit parameter (the mixed-event normalization the user asked about); C_OS=C_SS nominal.
  (9) **ΔR>0.05 REMOVED** from the signal selection everywhere (Design Decisions); data
  histos + truth templates must be regenerated; analysis_overview §2 updated; RDF
  `signal_cuts` flagged (other agent in that code). (11) §3f: pT>4/|η|<2.4 = acceptance,
  Δp/p/quality = reco-eff/fake; clean split requires verifying the pair reco-eff
  denominator. KB `gluon_splitting_flavour_excitation` re-synced (>2-HF-μ→systematic).
  **No physics code/behavior change; awaiting user approval to start steps 5–9.**

- 2026-06-22 (3) — **AUTONOMOUS IMPLEMENTATION APPROVED — orchestration started.** User:
  defer Δp/p (mark future-TODO needing MC); run preliminary with placeholder reco-eff +
  identity unfolding; start with reco-eff & ScrambGen; proceed autonomously as orchestrator,
  parallelize independent tasks, tracking-doc + commit each step. **T0 reco-eff/acceptance
  boundary VERIFIED** (parallel read-only investigation): no double-count — acceptance owns
  pT>4/|η|<2.4, denominator from `h_cutAcceptance` cutflow; no reco-eff change (Design
  Decisions). **T0 Δp/p deferred** in `placeholder.md` #9 + roadmap step 16. **ScrambGen
  plan** captured (Design Decisions; object-model rewrite, cuts in RDF, 20 ctr intervals,
  per-year, ×5 oversampling; needs single-muon-tree production first). Confirmed
  `signal_cuts` is already ΔR-free (RDFBasedHistFillingPP.cxx:382). Next: T1a single-muon-
  tree production (→ Condor) ∥ T1b ScrambGen rewrite (→ /review-analysis-code). Orchestration
  graph in Latest Stage.

- 2026-06-22 (4) — **INCIDENT + RECOVERY: T1a single-muon jobs clobbered PbPb nominal
  muon_pairs.** Root cause: `output_single_muon_tree` is a **protected** member
  (`DimuonDataAlgCoreT.h:32` re-exposes the public base member under protected — a
  regression vs the Oct-2025 era when the existing single-muon scripts worked), so the
  macro assignment `pbpb_2X.output_single_muon_tree = true;` failed to compile in cling;
  `Run()` then executed in DEFAULT mode (output_single_muon_tree=false) AND without
  `pbpb_run3_mu4_force_nominal=true` → `trigger_effcy_calc=TRUE` (no resonance veto) →
  each job opened the NOMINAL `muon_pairs_pbpb_2X_part*_single_mu4_mindR_0_02_res_cut_v2.root`
  with `recreate` and wrote wrong-mode (or truncated, for killed jobs) content. **Scope:**
  PbPb 23/24/25 nominal muon_pairs (+ hists_cut_acceptance) clobbered (Jun 22 02:35–02:40).
  **PP nominal SAFE** (pp single-muon script lacked `trigger_mode=3` → wrote a different
  filename, not the nominal `_2mu4_*`; pp nominal files remain Jun 10). **Response:**
  (1) `condor_rm yuhanguo` to stop the 17 still-running jobs; (2) resubmitted nominal PbPb
  production run_pbpb_2{3,4,5}_nominal.sub (clusters 86/87/88, 12 jobs) to deterministically
  restore from intact raw data (`pbpb_run3_mu4_force_nominal` is PUBLIC so nominal scripts
  work). **Pending fixes before re-running T1a:** make `output_single_muon_tree` publicly
  settable (header), and add `pbpb_run3_mu4_force_nominal=true` (pbpb) / `trigger_mode=3`
  (pp) to the single-muon scripts so muon selection matches nominal. Single-muon jobs are
  SAFE once output_single_muon_tree works (they write `single_muon_trees_*`, not muon_pairs).
  **Lesson:** verify a Condor job's first output is the intended filename before submitting
  a full fleet; the protected-member assignment failed SILENTLY in cling (non-fatal error,
  Run() proceeded).

- 2026-06-22 (5) — **INCIDENT recovery is INCOMPLETE — restore wrote the WRONG suffix;
  exact nominal recipe needs confirmation.** Findings:
  - The RDF crossx reads `muon_pairs_pbpb_2X_part*_single_mu4_mindR_0_02_**res_cut_v2**.root`
    (RDFBasedHistFillingPbPb.cxx:18/24/49/55). `_res_cut_v2` ⟺ `resonance_cut_mode=2`
    (DimuonDataAlgCoreT.c:398-401). These `_res_cut_v2` files are the ones clobbered.
  - My restore used `run_pbpb_2X_nominal.sub`, which sets only `pbpb_run3_mu4_force_nominal=true`
    → `resonance_cut_mode` stays default 1 → suffix "" → it produced PARALLEL
    `muon_pairs_pbpb_2X_part*_single_mu4_mindR_0_02.root` (02:55-56) and did NOT restore the
    `_res_cut_v2` files. The `_res_cut_v2` files remain CLOBBERED (zombie/missing-trees from
    the original incident).
  - **Ambiguity:** no run script sets `resonance_cut_mode=2`; `_res_cut_v2` is auto-set only
    when `trigger_effcy_calc=true` (plain `run_pbpb_2X.sh`, no force_nominal), but that mode
    is documented to SKIP the resonance cut — conflicting with the earlier R&O that the
    nominal OS tree IS resonance-vetoed. So the exact original recipe (force_nominal+res2 vs
    plain-trigger_effcy) is not certain from the scripts.
  - **No published results affected:** the downstream `histograms_real_pairs_*` (crossx/R_AA
    histogram outputs) were NOT regenerated and remain intact; only the upstream muon_pairs
    ntuples are damaged, which matters only for FUTURE RDF re-runs.
  - **April `_backup_20260420` backups** exist only for pbpb24 (incomplete, older skim) → not
    a reliable restore source.
  - **Action in progress:** running a no-clobber validation test (pbpb23 part4,
    `force_nominal=true + resonance_cut_mode=2 + is_test_run=true` → `_res_cut_v2_test.root`)
    to check it reproduces the expected OS resonance-veto structure. If it matches → full
    restore with that recipe. **Open question for user: confirm the exact recipe used to
    produce the nominal `_res_cut_v2` crossx ntuples on the May skim.**
  - Also: `run_pbpb_2X_nominal.sh` appear INCOMPLETE (don't set `resonance_cut_mode=2`), so
    they do NOT reproduce the crossx inputs — flag for fixing once recipe confirmed.

## Results & Observations

### TEMPLATE INVENTORY (for the fitting agent — 1D minv templates, Step 4a) ###
**Files** (Pythia evgen truth, `_no_data_resonance_cuts` = the truth sample has NO
resonance veto, unlike OS data):
- 5.36 TeV (nominal): `/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/histograms_pythia_5p36TeV_no_data_resonance_cuts.root`
- 5.02 TeV (variant): `.../pythia_5TeV/histograms_pythia_5p02TeV_no_data_resonance_cuts.root` (regenerate if needed — Step 4a only reran 5.36)

**Histogram name pattern:** `h_minv_zoomin_sign{1|2}_<category>_sigsel`
- `sign1` = SS, `sign2` = OS (`reference_sign_convention`).
- 1D TH1D, var `truth_minv`, **50 bins, 0–4 GeV**. Weight = generator `weight` (SHAPE
  templates — area-normalize in the fit). Truth level (no reco/eff/unfolding).
- Selection applied: `truth_pair_pt>8 && m1.truth_charge*m1.truth_eta<2.2 &&
  m2.truth_charge*m2.truth_eta<2.2` (NO ΔR cut, NO minv window) — matches the
  ΔR-removed data spectra D_OS/D_SS.
- 32 histos = 2 signs × 16 categories (9 flavor + 7 origin).

**Flavor categories** (`<category>` ∈): `single_b`, `bb`, `cc`, `one_b_one_c`,
`resonance`, `resonance_contaminated`, `photon_splitting`, `drell_yan`, `other_flavors`.
**Origin categories**: `FC`, `gs_FSR`, `phs_FSR`, `GS_ISR_no_HS`, `gs_ISR_one_HS`,
`diff_GS_same_HS`, `others`. (origin = production mechanism of the both-open-HF subset.)

**Fit-component mapping (§2,§3b):** with the OS spectrum (`sign2`):
- `T_S` (signal) = `h_minv_zoomin_sign2_single_b_sigsel`.
- `T_G` (gluon-splitting / open-HF) = sum of `sign2_{bb,cc,one_b_one_c}_sigsel`
  (use the origin set `{gs_FSR,phs_FSR,GS_ISR_no_HS,gs_ISR_one_HS,diff_GS_same_HS,FC}`
  if FC vs GS must be separated; `FC` = flavor creation, the rest = gluon splitting).
- `R_r` (smeared resonances, §3 resonance treatment) relate to `sign2_{resonance,
  resonance_contaminated}_sigsel` — but those need the no-veto reco-smeared build, NOT
  these truth templates directly.
- SS analogues (`sign1_..._sigsel`) drive the k-ratio study (Step 5: k=G_SS/G_OS).
**Caveat:** these truth templates have NO resonance veto; the OS *data* does (R&O below).
The fit must mask the OS data veto windows ([0,1.06],[2.9,3.3],[3.55,3.8]) consistently.

- **OS muon-pair tree is RESONANCE-VETOED (OS-only); SS is NOT — affects the fit
  (found 2026-06-21).** The new PP OS spectrum is EXACTLY 0 in minv ∈ [0,1.06],
  [2.9,3.3], [3.55,3.8] GeV while SS is smooth across all 0–4. Root cause:
  `DimuonAlgCoreT::ResonanceTaggingImpl` (NTupleProcessingCode/DimuonAlgCoreT.c:128)
  removes resonance-tagged muons **only for `op_sign`** (`if(op_sign)`); SS pairs are
  never tagged. Veto windows = `pms.minv_cuts` (ParamsSet.h:404-407): **[0,1.06]**
  (ρ/ω/φ + sub-threshold), **[2.9,3.3]** (J/ψ 3.097), **[3.55,3.8]** (ψ(2S) 3.686),
  [9.08,10.5] (Υ, outside 0–4). The nominal crossx ntuples are the `_no_res` variant.
  **Implications for the fit (§3c/§3d):** (1) the OS continuum the histos cover IS
  S+C+G as modeled — the J/ψ/ψ(2S)/low-mass resonances the §3c veto would remove are
  already gone, so no extra J/ψ veto is needed; (2) OS−SS is only valid in
  OS-populated bins — the veto windows MUST be masked in BOTH OS and SS before any
  subtraction or fit (else OS−SS = −C−kG < 0 there); (3) to ever SEE the J/ψ in OS
  (e.g. a validation plot of the raw 0–4 spectrum) needs a no-veto ntuple pass
  (`trigger_effcy_calc=true` skips the resonance cut). The data histos themselves are
  correct as specified (no mass cut in the RDF; identical OS/SS filter — the asymmetry
  is purely the upstream tree).
- **Pythia truth minv templates — usability finding:** truth `minv_zoomin` (0–4 GeV,
  50 bins) IS filled & plotted per flavor AND origin category, OS+SS — but (1) stored
  as 2D `{pair_pt, minv_zoomin}` (needs a projection to 1D), and (2) filled on the raw
  `df_op`/`df_ss` with NO single-b kinematic selection (`df_op` has no kinematic Filter,
  BaseClass.cxx:59; fill loops apply only the category filter). So they are NOT
  drop-in usable as fit templates that match the data selection. Fix (later cycle):
  add a truth fill with the matching selection (pair_pt>8, per-muon q·η<2.2, **no ΔR
  cut**, no minv cut) producing 1D `minv_zoomin` per category, OS+SS.

## Remaining Work
- Truth templates ALSO binned in (pair pT, pair η) for the k(pT,η) study (check if
  Step 4a output is pT/η-integrated; add a binned fill if so) — §3c, Step 5a.
- k validation: 5a MC truth k(pT,η,minv) + 5b data closure (SS ?= C+k·G_OS) — the GATE.
- Resonance-leakage quantification (§3h, Step 5'); settle 2.6-vs-2.9 signal edge.
- Mixed-event combinatoric: revive+modify+rerun `ScrambGen/` (Step 6 prereq —
  swap mass veto for signal selection, reconcile flat-vs-MuonObj input, rerun May skim).
- Coupled OS+SS fitter; closure test.
- Signal acceptance 2D MC map (log pT bins, fit/interpolate); apply after fit (§3f).
- Efficiency-order structural change (trig→unfold→reco-truth) when real unfolding +
  3D truth-binned pair ε_reco land (§3e; roadmap Q4/task_05).
- Wire acceptance-corrected signal yield into crossx & R_AA, replacing OS−SS (§3g).
- Δp/p fake-muon removal coordination (task_07 / roadmap step 16).
- Systematics (k incl. cc̄:bb̄ & FE:GS, template shapes, mixed-event norm, resonance
  leakage, fit model, fiducial-vs-extrapolated).

## Latest Stage
**AUTONOMOUS IMPLEMENTATION APPROVED & UNDERWAY (2026-06-22).** User approved running
the full chain to a PRELIMINARY result with placeholder reco-eff + identity unfolding;
Δp/p deferred (needs π/K MC); start with reco-eff (DONE: verified, no change) & ScrambGen.

**Orchestration / dependency graph** (✅=done, ▶=in progress, ⏳=queued):
- ✅ **T0 reco-eff/acceptance boundary** — verified no double-count; acceptance owns
  pT>4/|η|<2.4 via the `h_cutAcceptance` cutflow denominator (Design Decisions).
- ✅ **T0 Δp/p deferred** — placeholder #9 + roadmap step 16 marked "needs π/K MC".
- ▶ **T1 ScrambGen mixed-event (Step 6 prereq)** — biggest item, two sub-tracks:
  - T1a produce centrality-binned single-muon trees (May skim) → Condor (long pole, run first/in background).
  - T1b rewrite ScrambGen to object model (Design Decisions) → `/review-analysis-code` (parallel with T1a Condor).
  - T1c regenerate `nScramb`; T1d run ScrambGen (Condor); T1e fill T_mix via RDF.
- ⏳ **T2 data + truth template REGEN without ΔR** — confirm other agent's 1D Pythia
  histos dropped ΔR + are (pT,η)-binned; refill data D_OS/D_SS if still ΔR>0.05.
- ⏳ **T3 signal-acceptance map** (log pT, cutflow denominator) — INDEPENDENT of ScrambGen,
  parallelizable; needs the cutflow-based A_sig fill → `/review-analysis-code`.
- ⏳ **T4 resonance templates** (no-veto OS pass) — INDEPENDENT.
- ⏳ **T5 k-validation 5a/5b** — needs T1e + T2.
- ⏳ **T6 combined OS+SS fitter** — needs T2 + T1e + T4.
- ⏳ **T7 wire into crossx/R_AA** — needs T6 + T3.

**Current action:** T1a single-muon-tree production SUBMITTED to Condor 2026-06-22
(clusters 82=pbpb23×4, 83=pbpb24×2, 84=pbpb25×6, 85=pp24×12; 24 jobs). While they run:
T1b ScrambGen rewrite (code, no data dependency). Then monitor T1a → sanity-check
non-empty `muon_tree_ctr*`/`muon_tree` → hadd parts per year → T1c/T1d/T1e.
Parallelizable independent tracks: T1 (ScrambGen), T3 (acceptance), T4 (resonance) do not
conflict file-wise and can interleave; commits sequential (orchestrator).
