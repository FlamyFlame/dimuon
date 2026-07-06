# Low-mass dimuon template fit (OS/SS minv 0–4 GeV + Pythia truth templates)

**Mode:** Implementation. **Created:** 2026-06-18. **Session:** "template fit prep".
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.
**Siblings:** `analysis_overview.md` (§2,§4), `reco_eff_placeholder_run2.md` (Q1 dσ
weights), `raa_from_rdf_crossx.md` (OS−SS methodology), `project_overlay_pair_structure`,
`reference_sign_convention` (sign1=SS, sign2=OS).

---

## Sub-document map & current status (READ THIS FIRST)

This doc is the **umbrella + authoritative Physics Procedure + index**. Specific work lives in
five sub-docs (all in `Analysis/docs/tracking/`). New agents (new chat / post-compaction /
subagent): read this doc's Physics Procedure below, then go to the sub-doc that owns your thread.

| Sub-doc | Scope | Owns live status of |
|---|---|---|
| **A** `tf_nominal_fit_build.md` | **CURRENT NOMINAL method**: OS−SS + MC(S+G) 2-template fit + §3h smeared resonance templates; reco-level correction ordering; signal acceptance; wire N_sig into crossx/R_AA (replaces OS−SS). | The production fitter build → R_AA. |
| **B** `tf_k_factor_mixed_event.md` | OS→SS factor **k** study (5a/5b), mixed-event/ScrambGen T_mix, the **ABANDONED** combined OS+SS fit (kept as systematic/history). | k-validation + why the combined fit was dropped. |
| **C** `tf_bkg_composition_normalization.md` | Background **composition**, data/Pythia **normalization** (nb→pb unit fix), provenance classifier, tight-vs-medium WP study, extended-mass control region, V1 charge / V2 near-side charge-symmetry checks. | Data/MC comparison + composition. |
| **D** `tf_dpop_fake_muon_program.md` | Δp/p fake-muon **fit-and-subtract** yield program (D0–D5). | The muon-level Δp/p yield fit. |
| **E** `tf_upfront_bkg_reduction.md` | **Cut-upfront** counterpart to D: muon working points (χ² recommendation), |d0| discrimination, Δp/p distribution — reduce hadronic/fake by tighter selection instead of fitting. | Working-point + |d0| + Δp/p distribution study. |

**Current nominal background-subtraction method (authoritative):** the combined OS+SS/k-anchored
fit was **abandoned** (2026-06-24; the SS combinatoric has a soft near-side component no template
models — details in B). The nominal is **OS−SS (removes the charge-symmetric combinatoric
data-drivenly) + a 2-template MC fit `OS−SS = N_S·S + N_G·G_OS` + §3h smeared-resonance
templates**, per coarse R_AA bin, at RECO level (fit → then reco-eff/unfold → acceptance). This
**supersedes Physics Procedure §2's combined-fit model** (§2 is kept below for the k physics,
which now feeds only a systematic). Build + status: sub-doc **A**.

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
program has stages that must be ordered correctly (ordering CORRECTED 2026-07-01):
(1) trigger-efficiency correction (reco) + Δp/p fake-muon removal (reco); (2) the minv
template fit **at RECO level** on the trigger-corrected reconstructed spectrum, extracting
N_S; (3) origin-blind detector corrections (reco-eff + unfolding) applied to the
**extracted signal yield AFTER the fit**; (4) signal acceptance on that corrected yield.
The fit and all its templates use RECONSTRUCTED quantities (mixed-event = real-data muons;
fakes have no truth). See Physics Procedure §3 lead, §3e–§3g and §4.

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

**BOTH fits live at RECO level (CORRECTED 2026-07-01, per user — REVERSES the 2026-06-22
"minv fit at truth" decision).** The whole template-fit program is performed on
**RECONSTRUCTED quantities**, AFTER trigger-efficiency correction but BEFORE
reconstruction-efficiency correction and unfolding. Reco-eff + unfolding (and then signal
acceptance) are applied to the **extracted signal yield AFTER the fit**, never to the fit
input (see §3e for the ordering, §3f for acceptance).
- **Δp/p fake-muon fit → RECO level.** Δp/p = (p_ID − p_MS)/p_ID is the ID-vs-MS
  momentum imbalance — an **intrinsically reconstructed** quantity (no truth analogue),
  so the fake/real separation is done at reco level.
- **Combinatoric + physics (minv) fit → RECO level.** Two decisive reasons the minv fit
  MUST be at reco, not truth: (i) the **mixed-event combinatoric** template is built from
  **real-data muons** — a reconstructed object with no MC/truth analogue, so it is
  intrinsically a reco-minv shape; (ii) the **fake/hadronic background has NO truth match**
  (fakes are unmatched; punch-through/decay-in-flight hadrons have no signal-muon truth), so
  no truth minv exists for it. A truth-level fit could not represent either background.
  **Consequence — templates move to RECO:** the signal (S) and gluon-splitting/open-HF (G)
  MC templates must be **reconstructed minv from Pythia FULLSIM** (reco muon pairs,
  truth-labelled by category), NOT the truth `minv_zoomin` templates (Step 4a). The φ/J/ψ
  resonance leakage (§3h) is likewise a reco-smearing effect, naturally represented at reco.
  This settles the earlier "fit at reco vs truth" question: **minv fit at RECO.**

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
   - **⚠ SUPERSEDED (2026-06-22):** Steps 2/3 filled these histos in the NOMINAL crossx
     RDF, whose OS tree is **V1 resonance-vetoed** (holes at [0,1.06],[2.9,3.3],…) — wrong
     for the fit, which needs the resonances PRESENT (φ/J/ψ leakage templates, §3h). They
     MUST be re-filled from the **`_no_res_cut`** ntuples in a SEPARATE template-fit pass
     (see Design Decisions: "Template-fit input = `_no_res_cut`"). The `_no_res_cut`
     May-skim production is running (2026-06-22, Condor clusters 42–45).
b. **Templates (RECO level — CORRECTED 2026-07-01).** The fit is at reco (§3 lead, §3e), so
   the S and G templates must be **reconstructed `minv` from Pythia FULLSIM** (reco muon
   pairs, truth-labelled by flavor category: `_single_b` = S; `_bb`,`_cc`,`_one_b_one_c` = G),
   compared directly to the trigger-corrected reco data — NO unfolding of the data.
   Rationale: the mixed-event and fake backgrounds are reco-only, so ALL fit components must
   share the reco frame. **The truth `minv_zoomin` `_sigsel` templates (Step 4a) are for a
   truth-level fit and are SUPERSEDED for the nominal reco fit** (kept for the k=G_SS/G_OS
   MC-truth study, which is a ratio and frame-insensitive). Rebuild S/G as fullsim reco
   templates — see Remaining Work.
c. **Baby-step fit (later cycle):** assume k≈0 → form `D_OS−D_SS` (removes
   combinatoric), fit `= N_S·T_S + N_G·T_G` (TFractionFitter or binned χ²) over
   0–4 GeV, vetoing the J/ψ window [2.95,3.25]. Extract N_S.
d. **Refinement (later):** restore k; coupled OS+SS fit with shared C,G and a
   combinatoric template (event-mixing); float/constrain k.
e. **Correction ordering (CORRECTED 2026-07-01, per user — REVERSES the 2026-06-22
   "eff+unfold BEFORE the fit at truth" ordering).** The template fit runs at **RECO
   level** on the **trigger-corrected reconstructed spectrum**; reco-eff + unfolding are
   applied to the **extracted signal yield AFTER the fit**. Ordering:
   1. **Trigger efficiency** — evaluated as a function of **RECONSTRUCTED** kinematics
      (the trigger fires on reco objects). Applied as a per-pair weight at reco level →
      the **trigger-corrected reco spectrum** = the fit input. (Code: template-fit
      histograms use the TRIGGER-ONLY weight `crossx_weight_trig_only` (pp) /
      `weight_for_dsigma_trig_only` (PbPb) = 1/L·w_trig, **NO w_reco** — RDF
      `low_mass_template_calc` block.)
   2. **Template fit** — performed HERE, at reco, per coarse R_AA bin, with RECO-level MC
      templates (S, G from Pythia fullsim reco minv; §3b/§3 lead) + mixed-event
      combinatoric (real-data muons, reco) + §3h resonance leakage templates → extract the
      signal yield **N_S^fit** (still at reco level).
   3. **Reconstruction efficiency + unfolding** — applied to **N_S^fit AFTER the fit**, per
      bin. They are **origin-blind** (properties of the detector, the same function of pair
      kinematics for signal and background), and the fit is linear, so correcting the
      extracted signal yield per bin is equivalent to correcting the whole spectrum first —
      AND it is the only correct order given (i) mixed-event and (ii) fakes have no truth
      minv (§3 lead). reco-eff = the 3D pair ε_reco (truth-binned; §3f boundary); unfolding
      = spectrum-level bin migration (response matrix / iterative Bayes), a **structural**
      operation, not a per-pair weight.
   4. **Signal acceptance** — applied last, to the corrected signal yield (§3f).
   **Run 2 precedent (HF-muon R_AA, arXiv:2109.00411 / ANA-HION-2019-58):** the ρ (=Δp/p)
   template fit is run at **reco** to extract the HF-muon yield, and efficiency corrections
   are applied to that yield — i.e. **fit at reco → correct the yield after**, exactly this
   ordering. **Current-code status (RDF, 2026-07-01):** template-fit histograms are now
   trigger-only at reco (the fit input is correct); reco-eff/unfolding after-fit correction
   of N_S is future work (needs the real 3D pair ε_reco + unfolding, roadmap Q4/task_05).
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
- **Do the minv template fit at RECO level, on the trigger-corrected reconstructed
  spectrum (§3e, §3 lead).** Do NOT apply reco-eff or unfolding to the fit INPUT — they are
  applied to the extracted signal yield N_S AFTER the fit. Do NOT fit at truth (the
  mixed-event and fake backgrounds have no truth minv). The MC templates (S, G) must be
  RECO minv from Pythia fullsim, NOT truth `minv_zoomin`. **Do NOT apply signal acceptance
  to the mixture** — it is signal-only, applied to the corrected signal yield AFTER the fit
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
## Design Decisions (foundational, program-wide)

Method / k / ScrambGen / correction-ordering decisions were moved to the sub-docs:
OS−SS+MC(S+G) method + reco-level correction ordering → **A**; combined-fit-abandoned, gate
autonomy, k(m,pT,η) criterion, two code sets, ScrambGen rewrite/revival → **B**. The
foundational, program-wide decisions remain here:
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
- **Template-fit input = `_no_res_cut`, in a SEPARATE pass — NOT the nominal V1 crossx
  (DECIDED 2026-06-22 from the pipeline investigation).** Background: the investigation
  enumerated the nominal crossx histograms and found it DOES produce generic non-signal-cut
  histos that span the resonance region (case (b)): the 0–4 GeV `h1d_crossx_minv_0_4_{op,ss}`
  template inputs, plus (pp only) `minv_zoomin`/`minv_log`. So:
  - **Keep nominal/crossx on V1** (`_mindR_0_02`, `resonance_cut_mode=1`) for the signal-region
    crossx/R_AA program AND the diagnostic generic plots that expect clean V1 resonance removal.
    (Confirmed: the 06-08 PbPb nominal ran in NOMINAL mode — photoproduction cut APPLIED, V1
    OS veto — and the crossx reads the intact hadded V1 `_mindR_0_02.root`.)
  - **The low-mass template fit needs the resonances PRESENT** (to build the φ/J/ψ leakage
    templates, §3h) ⇒ it MUST read **`_no_res_cut`** (`resonance_cut_mode=0`), which CANNOT be
    folded into nominal (that would leave resonances in the nominal OS generic histos and break
    the OS−SS combinatoric logic that assumes the V1-vetoed OS tree). ⇒ build the template fit
    as a **SEPARATE pipeline mode/pass** that reads the `_no_res_cut` ntuples.
  - **Consequence:** the Step 2/3 `h1d_crossx_minv_0_4_{op,ss}_dsigma[_<ctr>]` histos (filled in
    the nominal V1 crossx, so their OS has V1 resonance HOLES) are **superseded** — they must be
    re-filled from `_no_res_cut` in the separate template-fit pass. (§3a updated.)
  - **Logic flow (answering the user):** separate no-res template-fit pass → per coarse RAA bin,
    fit {signal, GS, FE, mixed-event comb, φ/J/ψ resonance templates} on the 0–4 GeV OS+SS →
    signal fraction/yield per bin → acceptance-correct → feed crossx/R_AA (§3f/§3g). The nominal
    V1 crossx still provides the signal-region normalization; the no-res fit provides the
    background-subtracted signal yield that replaces OS−SS.
- **`_res_cut_v2` NOMINAL fallback REMOVED (2026-06-23, user directive 1).** Background:
  the data RDF `SetIOPathsHook` nominal/crossx branch had a multi-candidate fallback list
  (`_mindR_0_02.root` → `_res_cut_v2.root` → `_no_res_cut.root` → bare `.root`) built when
  the analysis was still on the OLD skim that lacked `mindR` branches, so it could limp
  along on a V2/no-res file if V1 was missing. The May-2026 skim HAS the mindR branches and
  V1 `_mindR_0_02.root` is always produced for nominal ⇒ the fallback is **obsolete and
  unsafe** (it could silently feed the wrong resonance-cut variant into crossx/R_AA). New
  behavior: the **nominal/crossx** branch (`!trigger_effcy_calc` for PP;
  `mu4_nominal_pbpb_NO_trig_calc` for PbPb) requires ONLY V1 `_mindR_0_02.root` and **throws
  a descriptive error if absent** — no silent fallback. The **trigger-efficiency** branch
  (`_res_cut_v2` first) is UNCHANGED (V2 is the correct trig-eff input — see
  `data_analysis.md` resonance-cut convention, memory `project_resonance_cut_modes`).
- **Subagent scratch-doc procedure (2026-06-23, CLAUDE.md update, user directive 3).** For
  any delegated semi-complex investigation/implementation, each subagent checkpoints to its
  OWN scratch doc (`Analysis/docs/tracking/_sub_<task>_<n>.md`), append-only, never the
  canonical doc; subagents NEVER run git or edit shared files. The orchestrator (this agent)
  merges each scratch doc into this canonical doc's Progress Log on return, THEN deletes the
  scratch file. Reviewer subagents (`/review-*`) are stateless and exempt.

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

## Progress Log (foundational; later per-area detail lives in the sub-docs)

Detailed progress moved with its topic: OS−SS+S+G fitter + correction-ordering reversal → **A**;
T_mix/ScrambGen/5a k-validation → **B**; nb→pb fix, plot batch, V1/V2 → **C**; Δp/p program → **D**;
upfront-reduction (WP/|d0|/Δp/p) → **E**. Foundational early entries kept here:
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

- 2026-06-23 — **T2 DONE: low-mass template-fit data mode (D_OS/D_SS from `_no_res_cut`)** —
  `/review-analysis-code` PASS iter 1 (log `review-analysis-code-20260623-181845-low-mass-template-fit-mode.md`;
  all numbers MATCH). New PUBLIC flag `low_mass_template_calc` (RDFBasedHistFillingData.h) → SetIOPathsHook
  (PP + both PbPb blocks) reads `_no_res_cut`, TriggerModeSettings appends `_template_fit` to the output;
  FillHistogramsCrossx early-return block fills `h1d_crossx_minv_0_4_{op,ss}_dsigma` + 2D
  `h2d_crossx_minv_0_4_vs_{pair_pt_log_150,pair_eta}_{op,ss}_dsigma` (PbPb per-ctr `_<ctr>`), dσ weight,
  selection signal_cuts MINUS minv (no dR), then returns (no signal-region crossx from `_no_res_cut`). New
  run scripts `run_template_fit_{pp24,pbpb23,pbpb24,pbpb25}.sh`. Ran all 4 → distinct
  `histograms_real_pairs_*_template_fit.root`. **OS resonances PRESENT** (pp24 [0,1.06]=3330, J/ψ=7131;
  pbpb23 ctr0_5 [0,1.06]=15950, J/ψ=21130 — all 0 under V1); SS smooth; nominal outputs UNTOUCHED.
  **INCIDENT+RECOVERY:** flag first placed PROTECTED → cling assignment silently failed → pp ran nominal mode
  (output_generic_hists=false) and overwrote the nominal pp output → made flag PUBLIC + restored nominal pp via
  the canonical nominal crossx (678KB, verified). Same root cause as the 2026-06-22 incident; LESSON: macro-set
  members MUST be public. **D_OS/D_SS now ready for the 5b closure + the fit.**
- 2026-06-23 — **Step A DONE (Task #1): `_res_cut_v2` nominal fallback removed** (user
  directive 1). `/review-analysis-code` PASS iter 1 (log
  `review-analysis-code-20260623-160554-remove-res-cut-v2-nominal-fallback.md`). Edited
  `RDFBasedHistFillingPP.cxx` SetIOPathsHook (nominal `!trigger_effcy_calc` branch) and
  `RDFBasedHistFillingPbPb.cxx` SetIOPathsHook (BOTH nominal `mu4_nominal_pbpb_NO_trig_calc`
  blocks: run_year 23/24/25 and 15/18): nominal now resolves ONLY V1
  `*_mindR_0_02.root` and throws a descriptive runtime_error if absent (states V1 required,
  fallback removed). Trigger-efficiency `else` branches (V2 `_res_cut_v2` first) UNCHANGED in all
  three places. ACLiC-clean both classes (separate sessions; PP `.so` 16:09, PbPb `.so` 16:13).
  **No rerun needed:** for the May-2026 skim the nominal hook now selects the SAME V1 file it
  already preferred (V1 was first in the old candidate list), so crossx/R_AA output is bit-identical;
  the change only removes the dangerous silent-fallback path. Verified all 4 V1 nominal files exist
  (pp24 589M; pbpb 23/24/25 245/174/492M). Committed to master.

## Remaining Work
- **Rebuild S/G MC templates at RECO level (fullsim reco minv, truth-labelled by category)**
  to match the reco-level fit (§3b/§3e, 2026-07-01 reversal). The truth `_sigsel` templates are
  superseded for the nominal fit (kept for the k=G_SS/G_OS ratio only). The bkg_mc_provenance
  reco-seeded provenance machinery already produces reco-minv category histos — reuse it.
- **After-fit correction of N_S:** apply reco-eff (3D pair ε_reco) + unfolding to the extracted
  signal yield per bin, THEN acceptance (§3e steps 3–4). Needs real ε_reco + unfolding (roadmap Q4).
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


---

## Appendix — ntuple production, pipeline & incidents (cross-cutting)

These blocks are cross-cutting infrastructure (not owned by one sub-area). Kept verbatim.

**✅ SESSION 2026-07-01 — orchestrated multi-task batch COMPLETE.** All 8 tasks done (details in the
2026-07-01 Progress Log entries + roadmap 2026-07-01 row). Git commits: `a257787` (Task A code+docs),
`3bfc9cc` (Task B ParamsSet coarse binning), `aee3ada` (Task H review-plot label rules). Data-area plot work
(Tasks C/F, D, G) is not git-tracked. Every plot batch passed `/review-plot` (updated with the new label rules)
after ≤1 cosmetic amend. Plan (below) kept for the record.
User directives this session (orchestrator = main agent; subagents scratch-doc + never git/shared-file per directive 3):

- **Task A (HIGHEST PRIORITY) — CORRECTION-ORDERING REVERSAL.** The minv template fit MUST be
  done **after trigger-efficiency correction but BEFORE reconstruction-efficiency correction and
  unfolding, using RECONSTRUCTED quantities** — NOT at truth level after unfold+reco-eff. This
  REVERSES the earlier §3/§3e decision ("minv fit at truth"). **Physics reasons (user):** (i) the
  mixed-event combinatoric template uses REAL-DATA muons (reco quantities, no MC/truth); (ii) the
  fake/hadronic background has NO truth match, so no truth quantities exist for it. New ordering:
  trigger-eff (reco) → **template fit (reco) → extract N_S** → reco-eff + unfolding → acceptance.
  Scope: correct pp+pbpb crossx pipeline (template-fit data-input weights = trigger-only, drop
  w_reco), this doc's Physics Procedure (§3 lead, §3e, negative constraints), and sweep the WHOLE
  repo for other affected code/docs. **Single git commit.** (Δp/p fake fit was already reco-level;
  this now makes the minv fit reco-level too — the two fits are now BOTH at reco.)
- **Task B** — nominal COARSE pair-pT binning [8-15, 15-27, 27-50, 50+] (4 bins, may shrink to 3)
  into ParamsSet.h as the SINGLE source of truth + `N_COARSE_PAIR_PT_BINS` (or auto-len); comment +
  memory that coarse binning is read from this file, never re-made-up per plot set. Fine log binning
  also lives there.
- **Task C** — bkg_mc_provenance_20260624: delete all pythia-TRUTH-quantity plots; Data/Pythia
  comparison → correct data by TRIGGER ONLY (not reco), because reco-quantity Pythia fullsim muons
  carry the same reco-eff (correcting data for reco-eff invalidates the comparison). (Consistent with Task A.)
- **Task D** — OS_to_SS_factor plots: fix (1) k_of_m_in_pt_slices legend showing BIN NUMBERS → use
  physical pair-pT ranges; (2) shape_overlay "sign1/2" → "same sign"/"opposite sign". New: OS/SS
  shape overlays in coarse pT bins (2 layouts); k(m,pair pT) template-fit 4-subplot ratio plots.
- **Task E** — explain hadronic/fake background-extraction procedure (bkg_mc_provenance); PP vs PbPb;
  R17618 vs R17662 overlay tag (see sibling `hijing_overlay_truth_barcode_duplicate_investigation.md`).
- **Task F** — inclusive 0-4 GeV minv plots must ALSO require pair-pT>8 & q·η<2.2 unless stated exception.
- **Task G** — extended-mass (0-20 GeV reco minv) Pythia dσ by origin (fullsim + overlay), hadronic/fake
  only, mixed-event minv, alt ranges → select a combinatoric-dominated control region to fix mixed-event norm.
- **Task H** — /review-plot skill rules: never sign1/sign2, never OS/SS as main text (subscripts OK),
  ranges never bin numbers, flag arbitrary-vs-physical wordings, analysis vars want equation defs;
  verify auto-trigger for template-fit plots; test with old OS_to_SS code via executor-review loop.

ORCHESTRATION: 3 read-only Explore agents dispatched (weights/ordering map; plot-macro map;
ParamsSet binning map). Implementation subagents partitioned by file to avoid concurrency conflicts;
orchestrator owns git + this doc. `.claude/logs` confirmed tracked-by-convention (56 files; prior
commit c503ec4) → committing them (fc87300) was correct.

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

**✅ INCIDENT RESOLVED — zero result-impact (2026-06-22).** The pipeline investigation showed the
crossx reads the **hadded** `_mindR_0_02.root` (V1) and trig-eff reads the **hadded** `_res_cut_v2.root`
(V2). Both hadded files are INTACT (Jun 8–10, untouched by the incident; V1 OS=745k/526k/1.49M valid).
The single-muon-job clobber only hit the intermediate **V2 part** files (`_res_cut_v2_part*.root`),
which the RDF does NOT read (regenerate only if trig-eff is re-hadded from parts; trig-eff is done).
Root cause (protected `output_single_muon_tree`) FIXED + committed. The nominal "restore" was
unnecessary (hadded V1 was always fine) but harmless — it produced V1 part files. **Cleanup later
(non-urgent):** zombie `_res_cut_v2_part*` (regenerate on demand), spurious `_mindR_0_02_part*` (02:55),
`_res_cut_v2_test`, the parallel `_mindR_0_02.root` re-hadd if any.

**File-usage clarified:** PbPb crossx reads intact hadded V1 (correct nominal, photoproduction+V1
veto applied); pp on the new skim since 2026-06-22. The d4c1197 `minv_cuts_v2` narrowing is real but
only affects V2/trigger-eff and the (separate) template-fit path — NOT nominal V1.

**Design DECIDED (investigation):** keep nominal/crossx on **V1**; build the low-mass template fit as
a **SEPARATE pass reading `_no_res_cut`** (resonances present for φ/J/ψ templates). The Step 2/3
`h1d_crossx_minv_0_4_*` histos (filled from V1 nominal → OS resonance holes) are superseded; re-fill
from `_no_res_cut`. See Design Decisions + §3a.

**DONE this round (2026-06-22):** (a) **PP `SetIOPathsHook` fix** committed — `!trigger_effcy_calc`
→ V1 first (mirrors PbPb); ACLiC-clean; `/review-analysis-code` PASS. (b) **`_no_res_cut` May-skim
production COMPLETE** (clusters 42–45, all batches wrote output) + **hadded** to (sizes):
`pbpb_2023/muon_pairs_pbpb_2023_single_mu4_no_res_cut.root` (290M),
`pbpb_2024/...` (206M), `pbpb_2025/...` (581M), `pp_2024/muon_pairs_pp_2024_2mu4_no_res_cut.root` (863M).
Verified resonances PRESENT in OS (pbpb23 p1: J/ψ=7493, φ=582, low[0,1.06]=5197 — vs 0 under V1).
These are the inputs for the separate template-fit pass.

**Non-urgent cleanup (unchanged):** zombie `_res_cut_v2_part*` (regenerate only if trig-eff
re-hadded), spurious `_mindR_0_02_part*` + `_res_cut_v2_test` from the incident.

---

## Latest Stage

**2026-07-06 — doc FACTORIZED into sub-docs A–E** (this umbrella + `tf_nominal_fit_build` A /
`tf_k_factor_mixed_event` B / `tf_bkg_composition_normalization` C / `tf_dpop_fake_muon_program` D
/ `tf_upfront_bkg_reduction` E). This doc now holds the authoritative Physics Procedure + index +
program-wide decisions; live per-area status is in the sub-docs (table at top).

- **Nominal fit build → R_AA:** sub-doc **A** (Latest Stage there).
- **Active new work (2026-07-06):** upfront hadronic/fake reduction (WP done; |d0| + Δp/p
  distribution plots in progress) — sub-doc **E**.
