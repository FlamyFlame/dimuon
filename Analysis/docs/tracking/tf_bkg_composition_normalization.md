# Low-mass template fit — background composition, data/Pythia normalization, provenance, WP, control region

**Parent:** low_mass_dimuon_template_fit.md (sub-doc C)
**Scope:** Background composition studies (provenance classifier, near-side charge symmetry V1/V2, tight-vs-medium WP), the data/Pythia-fullsim differential-cross-section normalization (the nb→pb unit bug + systematic cleanup + data efficiency-correction), and the extended-mass control-region selection. Supporting/history for the nominal OS−SS + MC(S+G) method (sub-doc A) and the Δp/p fake-muon program (sub-doc D).
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.

---

## Summary / current status

- **Data/Pythia normalization**: the 300–460× data/MC discrepancy was a **unit bug** — AMI `crossSection` is in **nb**, treated as pb → Pythia dσ 1000× too small. After nb→pb and the missing efficiency-correction of data, data and Pythia fullsim agree to ~O(1) (OS data/MC≈0.9, SS≈1.5). Systematic cleanup done (comments/labels/named-consts, behavior-preserving).
- **Composition**: with pair_pt>8, the MC low-mass background is dominated by (real μ + within-jet π/K hadronic fake) — OS-enhanced (does NOT fully cancel in OS−SS) and NOT the g→QQ̄ G template. Tight WP suppresses fake ~25–30% / hadronic ~20% but is insufficient; the Δp/p fake-muon program (sub-doc D) is required.
- **Charge symmetry (V1/V2)**: charge misreconstruction RULED OUT (agreement 99.994% of selected muons). The near-side low-mass background is OS-correlated (does not cleanly mirror); the OS−SS charge-symmetry assumption for the near-side peak is NOT settled by MC — flagged for the user.
- **Control region / combinatoric normalization (updated 2026-07-06)**: proposed m∈[6,9] GeV (signal-free, no resonances, at the mixed-event peak) but the correlated bb̄/cc̄ continuum persists there — a pure combinatoric-dominated MASS window does not exist. RESOLUTION (KB `run2_dimuon_backtoback_paper`): the combinatoric is charge-symmetric and **geometrically enhanced ∝ N_coll²**; it CANNOT come from the HIJING overlay (single embedded hard event → the HF×HF-from-different-scatters combinatoric is absent). So fix N_C from the **CHARGE structure of the data** (OS−SS needs none; T_mix by pair-counting / matching same-sign at wide angle), NOT from MC or a mass window. Data-driven evidence (extended-mass THStack + centrality study): central Pb+Pb data sit ~1 order of magnitude above the anchored correlated MC across 4–20 GeV (= the combinatoric the overlay can't make; T_mix reproduces the excess shape), and the same-sign [4,9] GeV yield grows **×19.3 peripheral→central vs ×2.7 for the correlated OS−SS** — the N_coll² vs N_coll signature.

---

## Progress Log

- 2026-07-06 — **Combinatoric control-region follow-up: extended-mass 0-20 GeV THStack + centrality-differential
  cross-check (user physics discussion → build).** KB-grounded (`run2_dimuon_backtoback_paper`: combinatoric =
  same-sign pedestal, geometrically enhanced in Pb+Pb, event-mixing only a systematic). PHYSICS settled with the
  user: combinatoric = charge-symmetric UNCORRELATED pairs ∝ N_coll² → NOT obtainable from the HIJING overlay (one
  embedded hard event; the dominant HF×HF-from-different-hard-scatters combinatoric is absent — only UE-hadronic-
  involving combinatoric partly present). "Combinatoric" (pair-correlation) ≠ "hadronic" (single-muon origin): a
  combinatoric pair can be real+real. ⇒ combinatoric is DATA-DRIVEN (same-sign/OS−SS/mixed-event), shows as a data
  EXCESS over the correlated MC growing toward central; normalize N_C via the charge structure, not MC/a mass window
  (no clean combinatoric-only mass window — bb̄/cc̄ persists to 20 GeV).
  - **CODE (committed `3f481ce`, `/review-analysis-code` PASS):** added extended-mass 0-20 GeV (40 bins) DATA +
    mixed-event histos to the template-fit RDF (`h1d_crossx_minv_0_20_{op,ss}_dsigma[_ctr]` + 2D vs coarse pair-pT,
    from `pms.pair_pt_coarse_bins`), trigger-only, PP + PbPb per-ctr. Reran all 8 passes; verified (pp OS int 1.69e4;
    PbPb central ctr0_5 1.07e5 ≫ peripheral ctr50_80 1.73e4).
  - **PLOTS (`extended_mass_control_region_20260701/thstack/`, `/review-plot` PASS after 1 cosmetic amend; all
    numbers MATCH):** `thstack_extmass.C` → (1) THStack pp baseline (correlated MC stack + data + T_mix); (2) THStack
    PbPb-central (overlay stack ANCHORED to central OS data in [1.1,2.85] GeV, factor **×6.07**), data excess =
    combinatoric; (3) centrality-differential data-driven cross-check. **RESULT:** central Pb+Pb data ~1 order of
    magnitude above the anchored correlated stack across 4–20 GeV; T_mix reproduces the excess; **same-sign [4,9] GeV
    yield ×19.34 peripheral→central vs ×2.71 for correlated OS−SS ([1.1,2.85])** → combinatoric ∝ N_coll² vs signal
    ∝ N_coll. CAVEATS (SUMMARY.md): HF-filtered MC → pp absolute not 1:1; centrality trend is a qualitative relative
    cross-check (dσ-normalized), not an absolute N_coll² fit; 10-20% bin high (per-ctr norm artifact). IMPLICATION:
    fix mixed-event N_C from the data charge structure + validate its N_coll² scaling via the data-vs-overlay
    high-mass excess — not from a mass control region.

- 2026-07-01 (2) — **Plot batch (Tasks C/F, D, G) — 3 executor subagents + `/review-plot` loops, all PASS after
  ≤1 cosmetic amend; scratch docs merged here + deleted.** Data-area plots (not git). Reviewer used the NEW
  label-semantics rules (committed `aee3ada`).
  - **bkg_mc_provenance_20260624 (Tasks C+F):** deleted the ONLY Pythia-truth-quantity plot
    (`data_vs_fullsim_truth_xsec.png`, both plots/ + plots_tight/) and removed its code block. Data-vs-fullsim
    comparison now uses TRIGGER-ONLY data (the regenerated `h1d_crossx_minv_0_4_{op,ss}_dsigma` — trigger-only
    after Task A) and applies the signal cuts (pair_pt>8, q·η<2.2) to BOTH data and MC identically (switched
    data `_nosel`→selected, MC `h_nosel_*`→`h_minv_*_all`); the one no-selection composition (category_stack)
    is kept but LABELLED "no pair selection" with data+MC both no-sel. Physics reason recorded: fullsim MC is
    reco → already carries reco-eff, so correcting data for reco-eff would invalidate the comparison. Data/MC
    (0–4 GeV, selected): medium OS 0.71 / SS 1.12, tight OS 0.79 / SS 1.26 (order-1; reviewer re-derived, MATCH).
    `/review-plot` PASS after 1 amend (fixed a bare `(OS+SS)` y-title → `(opposite + same sign)` + legend/data
    overlaps).
  - **OS_to_SS_factor_..._20260623 (Task D):** fixed all label bugs the Task-H test flagged (bin-number legends
    → physical pT ranges; `(sign1)/(sign2)` → `(same sign)/(opposite sign)`; `k(m)=SS/OS` → `k(m)=G_{SS}/G_{OS}`;
    `one_b_one_c` → "one b, one c"). NEW plots via `code/os_ss_shapes_and_kfit_coarse.C` (reads
    `pms.pair_pt_coarse_bins`/`N_COARSE_PAIR_PT_BINS` from ParamsSet.h, coarse edges snapped to fine
    `pair_pt_log_150` boundaries, labelled by physical ranges 8.00–14.38 / 14.38–25.84 / 25.84–46.44 /
    46.44–150 GeV): (A) `shape_overlay_by_coarse_pt_OSandSS.png` (2 subplots OS/SS × 4 coarse-pT lines),
    (B) `os_vs_ss_shape_per_coarse_pt.png` (4 subplots, OS vs SS per range), (C) `k_of_m_fit_ratio_per_coarse_pt.png`
    (4 subplots, ratio: top G_SS vs G_OS·k(m) with pol3 k(m) fit, bottom k(m)=G_SS/G_OS). k_int=0.3078; per-coarse-bin
    k=0.305/0.314/0.331/0.363 (rises with pT, bb̄ hardening). `/review-plot` PASS after 1 amend (2 cosmetic). Old
    macros backed up to `code_backup_pre_label_fix_20260701/`.
  - **extended_mass_control_region_20260701 (Task G, NEW dir):** reco m_{μμ} 0–20 GeV, signal cuts, coarse-pT from
    ParamsSet. 24 PNGs: dσ/dm by origin (flavor {single_b,bb,cc,b+c,resonance,other} + hadronic + fake), pT-int +
    4 coarse bins, fullsim + overlay; backgrounds-only; mixed-event; a 0–10 zoomed subdir. **Findings (reviewer
    C1–C4 PASS, numbers MATCH):** single-b SIGNAL dies at m≈4 GeV (pT-independent; frac(m>5)≈0); the CORRELATED
    bb̄ continuum does NOT die out — dominates 4–20 GeV (bb̄ OS pT8–15 ∫≈155/144/227/301 pb over [4,6]/[6,8]/[8,12]/
    [12,20]); cc̄ ~4–6× smaller; `one_b_one_c` empty; hadronic flat ~3–10 pb/GeV no peaks; fake low/spiky; mixed-event
    OS≈SS (charge-symmetric), smooth, peak-free, broad max m≈7–9 GeV. **Proposed control region m∈[6,9] GeV
    (pair_pt>8):** signal-free, no resonances (above ψ(2S) veto 3.8, below Υ 9.08), at the mixed-event peak —
    **CAVEAT (correctly stated): the correlated bb̄/cc̄ continuum persists there, so the mixed-event normalization
    must be fit to DATA in [6,9] with bb̄/cc̄ subtracted or fit simultaneously** (a pure combinatoric CR does not
    exist in signal MC). `/review-plot` PASS after 1 amend (legend-off-data). SUMMARY.md in the dir.

## Results & Observations — data/Pythia normalization (2026-06-27 → 2026-06-30)

**✅ BLOCKER RESOLVED (2026-06-30) — it WAS a unit bug: AMI `crossSection` is in nb, treated as pb → Pythia dσ 1000× too small.**
FIX (macro-only, per user "fix macro & regenerate first, then suggest systematic cleanup"): `bkg_mc_provenance_20260624/code/fill_weighted_fullsim.C` — `ami_weight()` now returns `xs·eff·NB_TO_PB` (NB_TO_PB=1000.0, nb→pb), so weighted_fullsim_{signal,overlay}.root are in pb; `plot_diff_xsec.C` — panel (b) raw production-weight path ×1000 (nb→pb), data side untouched, on-plot `data/MC` annotation `%.0f`→`%.2f`. Re-ran fill + plots: signal no-sel OS 19.22→1.922e4 (exactly ×1000); **data/MC OS=0.31, SS=0.46 (was 314/459) → SAME ORDER OF MAGNITUDE.** 5 PNGs regenerated (`data_vs_fullsim_{reco,truth}_xsec`, `provenance_stack_xsec`, `survives_quality_cuts_xsec`, `category_stack_vs_data`). MC ends ~2-3× ABOVE data (acceptable order-of-mag; fullsim test-sample normalization caveat + resonances present in both). `/review-analysis-code` loop opened; reviewer subagent stopped by user → loop closed at iter 1 with executor verification (pure ×1000 scaling: no shape/discontinuity change, data untouched). Log `review-analysis-code-20260630-152104-fullsim-nb-pb-unit-fix.md`. NOT git-tracked (data area).

**SYSTEMATIC CLEANUP — DONE (2026-06-30, single commit; user-requested).** Decision: keep stored MC `weight` in **nb** (AMI native) — do NOT change the production weight (would force a ntuple rerun) — and document the convention + convert nb→pb at the consumer. Items:
- (1) ✅ `PythiaAlgCoreT.c` AMI-read comment/print "pb"→"nb" + convention note at the weight definitions (fullsim `fullsim_weight_factor`, truth `w_norm`). Comment-only → no recompile/rerun.
- (4) ✅ `RDFBasedHistFillingPP.cxx` PP-data crossx dσ y-axis labels `nb GeV⁻¹`→`pb GeV⁻¹` (4 occ: lines 519/520/674/677; values ARE pb since L is pb⁻¹). Label-only; existing .root titles stale until next crossx re-fill (cosmetic, values unaffected).
- (2) ⚠ DOCUMENTED, not removed: `plot_mc_data_compr.cxx` `Scale(pow(10,6))` → `NB_TO_PB(1e3)·TRUTH_COMBINE_RENORM(1e3)` named consts + comment (identical value, behavior unchanged). Cannot reduce to the clean 1e3 until (3) is fixed; removing it would break the known-good data-vs-Pythia overlay.
- (3) ⏸ DEFERRED: the Pythia TRUTH-combined per-slice-combine under-normalization (~1e3; weight branch 0→3.46e6) needs fixing at the truth-combine step (likely a ntuple/combine rerun). Immaterial to the analysis (area-normalized templates, k=G_SS/G_OS ratio). Flagged in the mc_data_compr comment + here.
BLAST RADIUS: none realized — all edits behavior-preserving (comments/labels/named-const). Memory `project_ami_crosssection_nb_units` records the nb fact.
NEXT: method 6 — signal-region (pair_pt>8, q·η<2.2) data/MC comparison with `_no_res_cut` data; then resume BUILD (OS−SS + S+G fitter → acceptance → R_AA).

**✅ DATA EFFICIENCY-CORRECTION (2026-06-30, user: "data SHOULD have both trigger & reconstruction efficiency applied").** Added no-pair-selection eff-corrected 0–4 GeV data dσ to the `low_mass_template_calc` RDF block (`RDFBasedHistFillingPP.cxx`): `h1d_crossx_minv_0_4_{op,ss}_dsigma_nosel`, weighted by `crossx_weight_trig_corr = 1/L·w_reco·w_trig` (reuses the validated efficiency lambda; muon-level selection only, NO pair_pt/q·η). Reran `run_template_fit_pp24.sh` (reads `_no_res_cut`). `plot_diff_xsec.C` now reads these as the data points (legend "ε-corr."). **RESULT: data & Pythia fullsim now agree to ~1** — OS data/MC=0.94 (shapes track across m<4, J/ψ peak matches), SS=1.48 (data > MC: extra uncorrelated SS combinatoric, expected). ε-correction raised data ~3× (ε_total≈0.33) — confirming the prior 0.31/0.46 residual WAS the missing trigger/reco efficiency, not a bug. CAVEAT (in code): soft low-mass bins carry turn-on/clamp artifacts (1/ε well-defined only on the pair-pT≳8 plateau). Tight plots also regenerated (`plots_tight/`, data/MC OS=1.05, SS=1.66; data WP stays medium — only MC WP varies).

**✅ TIGHT-WP fake/hadronic study DONE (2026-06-30, user; delegated subagent, scratch merged + deleted).** Added `useTight` (quality bitmask MEDIUM `1|8|32|256`→TIGHT `1|16|32|256`) to `fill_weighted_fullsim.C` + `bkg_mc_provenance.C` (data-area, NOT git). NO ntuple rerun — `muon_quality` already encodes Tight (bit16). Tight outputs `*_tight.root` + `plots_tight/` (medium preserved). `SUMMARY_tight_vs_medium.md` written. **RESULT (HIJING overlay, sig sel, near-side m<1.5):** tight/medium yield ratio OS fake 0.72 / had 0.77 / real 0.91; SS fake 0.40* / had 0.80 / real 0.94 (*low stat). Tight preferentially removes fake (~25–30%) & hadronic (~20%) over real (~10%), BUT absolute bkg fraction stays substantial (hadronic ~2% OS / ~13–18% SS). **⇒ tight helps but is NOT sufficient; the Δp/p fake-muon program is still required.** Dominant low-mass bkg = hadronic (π/K decay-in-flight, punch-through), only ~20% suppressed by tightening. data/MC absolute dσ rises OS 0.31→0.35, SS 0.46→0.51 (tight removes MC bkg).

**🔴 (resolved above) OPEN BLOCKER — DATA/PYTHIA NORMALIZATION DISCREPANCY (treat as a BUG, not physics) — 2026-06-27 (user).**
**ISSUE:** the differential-crossx data-vs-Pythia-fullsim plots (`bkg_mc_provenance_20260624/plots/data_vs_fullsim_{reco,truth}_xsec.png`, `category_stack_vs_data.png`) show **data ≈ 300–460× the Pythia fullsim** at m<4, no selection (data/MC ≈ 459 SS, 314 OS). Data dσ/dm = N/L/Δm (L=400.412 pb⁻¹); Pythia dσ/dm = Σw/Δm (w = crossSection·genFiltEff·beam_ratio/N_slice). **My earlier "expected physics" reasoning was WRONG and rejected by the user:** Pythia fullsim DOES contain resonances and the fake/hadronic backgrounds (we have been using Pythia fullsim for the fake & hadronic templates all along), so it MUST agree with data at least to the same order of magnitude / crossx factor. A 0.2–0.3% ratio is **almost certainly a normalization/unit bug**, not a physical statement about the Pythia model.
**GOAL:** do ALL investigation needed to get data and Pythia on the **same order of magnitude** (rule out a code/normalization bug) BEFORE any conclusion about whether the Pythia model itself agrees with data. Only after the units/normalization are proven correct is the data-vs-Pythia comparison valid for judging the generator.
**PROPOSED METHODS (investigate, in order — NORMALIZATION FIRST):**
1. **Luminosity / unit bookkeeping (FIRST GUESS).** Cross-check L=400.412 pb⁻¹ (PPBaseClass::CrossxFactorMap {24,"2mu4"}=1/400.412) and the Pythia `ami_weight` units (AMI `crossSection` — pb? the file gives 4.816e6 for pp pTH8_14; PythiaAlgCoreT.c:450 prints "pb") against the **EXISTING pp data–MC comparison plots run in the PP crossx pipeline** (Stage 8; POWHEG+Pythia vs pp24 data — memory `project_mc_data_compr_plots`; they presumably AGREE). Match that pipeline's exact normalization recipe; find where my standalone macro differs (a pb↔nb↔µb factor, a missing 1/N or prescale, etc.).
2. **Data normalization.** Does the data need a per-event **prescale / trigger weight** beyond 1/L? Is my all-reco-pair count (TChain over muon_pair trees, m<4, no sel) consistent with how the pipeline counts di-muon candidates? (the existing crossx D_OS integral was 1.71e4 with the SIGNAL selection — compare scales.)
3. **Pythia normalization.** Re-derive `fullsim_weight_factor = ami_w·nom_ratio/N_beam` (PythiaAlgCoreT.c:665); confirm no missing factor; confirm AMI crossSection units. Confirm the reco-seeded fake/hadronic fills (`fill_weighted_fullsim.C`) use the SAME weight.
4. **Resonances present in Pythia.** Ensure the Pythia processing is NOT applying a resonance cut that removes J/ψ etc. (rerun pythia/fullsim ntuple processing with `resonance_cut_mode=0`/no-res-cut if needed) so the OS J/ψ peak is present in the MC — else the OS shape/normalization is wrong.
5. **Fake/hadronic the production way (optional, cleaner).** Add a mode to the Pythia fullsim (+overlay) ntuple-processing that **pairs and saves ALL reconstructed muons** (incl. fakes/hadronic passing the single-muon cuts), weighted to differential crossx — instead of (or to validate) the standalone reco-seeded `fill_weighted_fullsim.C`.
6. Once data & Pythia are within ~O(1) → THEN judge generator agreement; THEN also make the **signal-region** (pair_pt>8, q·η<2.2) data/MC comparison where the single-b HF signal dominates and Pythia is the intended model.
**Plots already made (differential dσ/dm, weighted, L=400.412):** data_vs_fullsim_{reco,truth}_xsec, provenance_stack_xsec, survives_quality_cuts_xsec, category_stack_vs_data (all in `bkg_mc_provenance_20260624/plots/`). G audit answered: G=bb+cc+one_b_one_c (all open-HF non-signal = GS+FE+FC), weighted. Labels "real, truth-matched μ" (NOT prompt — memory `feedback_hf_muons_not_prompt`). Legacy unweighted/unity plots deleted. **DO NOT trust the magnitude of these data/MC plots until the normalization bug is found.**

**▶ INVESTIGATION CHECKPOINT (2026-06-29) — ROOT CAUSE IDENTIFIED: AMI `crossSection` is in nb, not pb (1000× unit bug). Verification + reconciliation of the truth path in progress.**
Plan: (method 1) compare the standalone macro's normalization recipe against the known-good crossx pipeline + mc_data_compr; (then) verify the AMI unit; reconcile fullsim-vs-truth; document; route fix via /review-analysis-code + /review-plot.
FINDINGS so far (ROOT, this session):
- **DATA is NOT the bug.** Standalone data dσ/dm OS (m<4, no sel, /L=400.412) ≈ 6035 pb ≈ the known-good crossx pipeline `h1d_crossx_minv_0_4_op_dsigma` = **6136.7** (signal-sel, w_reco·w_trig, V1 veto), SS pipeline = 915.95. The standalone uses *only* 1/L (a SUBSET of the pipeline's `crossx_weight·w_reco·w_trig`), so if anything data is *under*-corrected, not inflated. Data dσ OS ≈ 6×10³ pb is real and pipeline-consistent.
- **Discrepancy PERSISTS in the signal region** (NOT a no-selection/phase-space artifact): fullsim signal-sel (pair_pt>8) OS = 18.32, SS = 0.622 → **data/MC OS=335, SS=1474**. So applying pair_pt>8 does not fix it → genuine normalization, as the user said.
- **AMI unit = nb, not pb (the 1000× bug).** `ami_info_..._pp_..._pTH8_14.txt`: `crossSection=4816000`, `genFiltEff=4.938681e-6` (= sumOfWeights 2.221e5 / sumOfWeightsNoFilter 4.497e10 — genuine xAODMultiMuonFilter eff). The macro parses correctly (ami_w=xs·eff=23.78). **PHYSICS:** for HardQCD:All pTHat 8–14 GeV at 5.36 TeV, σ≈few mb; 4.816e6 **nb** = 4.816 mb ✓; 4.816e6 pb = 4.816 µb is absurdly small ✗. ATLAS AMI convention = nb. The standalone macro AND `PythiaAlgCoreT.c:449` comment treat it as **pb** → MC dσ is 1000× too small. After nb→pb: fullsim OS = 18.32 nb = 1.832e4 pb vs data 6137 pb → **data/MC = 0.34** (same order ✓; 335 = 1000×0.335).
- **Corroboration:** `plotting_codes/mc_data_compr/plot_mc_data_compr.cxx:94` already applies a hand-tuned `h_pythia->Scale(pow(10,6))` to reach data scale (y-axis `[pb]`) — a symptom of the SAME unit confusion (truth weight too small). Also `RDFBasedHistFillingPP.cxx:519` labels the crossx y-axis `nb GeV⁻¹` while the factor is 1/400.412 with L in **pb⁻¹** → that label is inconsistent too (separate cosmetic issue; magnitude is pb).
- **TRUTH-path normalization is a SEPARATE quirk — does NOT affect the blocker or the analysis.** The Pythia TRUTH combined (`pythia_private_sample/...combined...root`) OS `h_minv_sign2` integral=0.0415 (weight-sum 0.0426), vs fullsim macro OS=18.3 — a ~440× gap despite the same `weight=ami_w·nom_ratio/N_beam` formula. The truth COMBINED tree's `weight` branch spans 0→3.46e6 (vs fullsim 0→1198): the per-pTHat-slice combine left an inconsistent absolute normalization, which `mc_data_compr` papers over with the hand-tuned `×10⁶`. **This is immaterial because every place the truth weight is consumed it is either area-/shape-normalized (the fit templates — doc §3b "SHAPE templates, area-normalize"; the k=G_SS/G_OS ratios) or a RATIO (reco-eff), so the absolute unit cancels.** It only bites where an ABSOLUTE truth dσ is shown vs data without compensation. ⇒ NOT on the blocker's critical path; flag for a later standalone cleanup, do not block.
- **ROOT CAUSE (confirmed): standalone fullsim macro `fill_weighted_fullsim.C` treats AMI `crossSection` (nb) as pb → MC dσ 1000× too small.** After ×1000 (nb→pb), signal-region OS data/MC=0.34, SS≈1.5 — same order of magnitude. Residual (MC ≳ data) expected: data crossx is V1 resonance-VETOED (J/ψ/φ/ρω removed — large in OS m<4) while fullsim INCLUDES resonances; clean comparison wants `_no_res_cut` data (method 4/6).
- **FIX SCOPE (decision pending user):** (a) MINIMAL/blocker — `fill_weighted_fullsim.C` ami_w nb→pb (×1000) + honest `[pb]` labels; regenerate the 5 plots → demonstrate O(1). (b) SYSTEMIC — also fix the `PythiaAlgCoreT.c:449` "pb"→"nb" comment, decide whether the production weight should carry the unit, and remove the `mc_data_compr` `×10⁶` hack once units are coherent (BLAST RADIUS: reco-eff plots, all MC dσ plots — most cancel, but must be audited). Route: /review-investigation (verdict) → /review-analysis-code (fix) → /review-plot (regen). Awaiting user on (a)-only vs (a)+(b).

## Results & Observations — near-side charge symmetry (V1/V2, 2026-06-24)

The user requested two sanity checks on the foundational claim (low-mass SS peak <1 GeV = soft
near-side combinatoric from within-jet fakes/hadrons → cancels in OS−SS):

- **V1 — charge-misreconstruction ruled out (user directive: "change turn_on_track_charge to true in
  ntuple process").** Concern: at low minv (small opening angle) a muon's reco (combined) charge could
  be misreconstructed, so a true OS pair leaks into SS, distorting the SS low-mass spectrum. Independent
  handle: the ID-**track** charge `m.trk_charge` (from `muon_pair_muon{1,2}_trk_pt` sign) vs the combined
  `m.charge` (from `muon_pair_muon{1,2}_pt` sign). The cut already exists —
  `DimuonDataAlgCoreT.c:605` requires `m1.charge==m1.trk_charge && m2.charge==m2.trk_charge` — but is
  gated by `turn_on_track_charge` (default **false**, `DimuonDataAlgCoreT.h:343`; never set true anywhere
  → trk_charge stored as 0 in all existing pairs). **Plan:** (a) ntuple-processing change → set/enable
  `turn_on_track_charge=true` so trk_charge is populated and the consistency cut applies; **WRITE TO A
  DISTINCT OUTPUT SUFFIX** (e.g. `_trkqcut`) so nominal/`_no_res_cut` pairs are NOT clobbered (incident
  history: nominal clobbered twice). → `/review-analysis-code`. (b) Re-process pp24 `_no_res_cut` (Condor)
  with the flag → distinct pairs. (c) RDF `low_mass_template_calc` fill → D_OS/D_SS with the cut.
  (d) Compare low-mass (<1 GeV) SS minv WITH vs WITHOUT the trk-charge cut + the per-SS/OS
  `pass_muon_trk_charge` cutflow bin → fraction of SS pairs removed = charge-flip contamination. Expect
  SMALL (user: SS peak shows no resonance structure → unlikely charge misreco). → `/review-plot`.
  NOTE: keeping `turn_on_track_charge` default false (nominal unchanged); the flag is set true only in the
  dedicated check run. User may adopt it as nominal afterward (then blast radius = re-run crossx/R_AA).
- **V2 — OS/SS mirror of the low-mass peak, with background MC.** Claim to support: the low-mass peak
  mirrors in SS and OS (charge-symmetric combinatoric) → OS−SS removes it. **Plan:** (a) DATA overlay
  D_OS vs D_SS (0–4 GeV, from existing `*_template_fit.root`: `h1d_crossx_minv_0_4_{op,ss}_dsigma`) — show
  the soft peak present in both. (b) Background MC: HIJING overlay (`pythia_fullsim_hijing_overlay_test_sample/
  histograms_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts.root`) and/or pythia fullsim — show
  the low-mass SS≈OS shape in the background/combinatoric-dominated region (charge symmetry), and that the
  near-side soft component is NOT signal. → `/review-investigation` + `/review-plot`. (Connects to R&O
  "ROOT CAUSE RESOLVED": 40% of low-mass SS pairs ΔR<0.3, peak minv 1.25.)

These gate confidence in the BUILD-phase OS−SS+S+G method; BUILD steps 2–5 (PbPb fitter, acceptance, R_AA
wiring, systematics) resume after V1/V2.

**V1 LATENT BUG FOUND + FIX (2026-06-24).** The first trkqcut re-process produced **OS=0**, SS halved
(part1: SS 145852→73911, OS 576516→0). Root cause: the PRE-EXISTING `turn_on_track_charge` code derived
`m.trk_charge` from `sign(muon_pair_muon{1,2}_trk_pt)`, but `trk_pt` is stored **UNSIGNED** in the skim
(verified: 0 negative) → trk_charge≡+1 → the consistency cut `trk_charge==charge` rejected EVERY negative
muon → all OS pairs (and −− SS) gone. This path was NEVER exercised before (flag default false), so no past
result is affected. **The real fix:** the raw skim DOES have signed `muon_pair_muon{1,2}_trk_charge`
(`vector<int>`, ±1; EXISTS in `data_pp24_part1.root`) — the combined charge is `sign(muon_pair_muon{1,2}_pt)`
(no `_charge` branch), the track charge is the `_trk_charge` branch. So `m.trk_charge` must be read from the
`muon_pair_muon{1,2}_trk_charge` branch, NOT `sign(trk_pt)`. The 12 garbage trkqcut outputs will be
overwritten after the fix. (Fullsim has the SAME latent bug — `muon_trk_pt` unsigned — but is out of V1 scope;
note for future.) → fix via `/review-analysis-code`, re-run Condor.

**V1 RESULT — CHARGE MISRECONSTRUCTION RULED OUT (2026-06-24, DONE).** After the trk_charge-branch fix,
re-processed pp24 with the charge-consistency cut (`turn_on_track_charge=true`, cluster 58, `_trkqcut` output).
The cut removes ≈0 pairs: trkqcut SS/OS entry counts = nominal `_no_res_cut` to 4 decimals across 5 completed
parts (OS 831515, SS 210383; OS-kept=SS-kept=1.0000). Direct disagreement rate `sign(muon_pair_muon1_pt) ≠
muon_pair_muon1_trk_charge` (part1, 5M events): **0.376% over ALL muons** (confirms the two charges ARE
independent measurements — cut not a no-op) but **0.0058% after quality+|η|<2.4+pt>4 selection** (32/556407).
⇒ combined-muon charge and ID-track charge agree for **99.994%** of selected muons → charge misreconstruction is
NEGLIGIBLE → it does NOT cause the low-mass SS peak; OS→SS charge-flip leakage ≈0. (User's expectation confirmed.)
Note: the per-muon disagreers preferentially fail the other cuts (dP/p etc.), so the pair-level removal is even
smaller (~exactly 0). The remaining trkqcut Condor parts finish but only reconfirm. The signed-track-charge fix
is now in the codebase for any future use.

**V1 STATUS (2026-06-24):** suffix code DONE — `/review-analysis-code` PASS iter 1 (log
`review-analysis-code-20260624-021327-trk-charge-check-pass.md`): `DimuonDataAlgCoreT.c` adds
`trkqcut_suffix="_trkqcut"` (gated on `turn_on_track_charge`) to BOTH output_file_path + output_hist_file_path
→ nominal byte-identical when flag false; new `run_pp_24_no_res_cut_trkqcut.{sh,sub}`. **Condor SUBMITTED**
cluster 57 (12 jobs) re-processing pp24 raw skim → `muon_pairs_pp_2024_part*_2mu4_no_res_cut_trkqcut.root`
(cut ON). Monitor → hadd → compare low-mass SS minv with vs without trk cut + `hists_cut_acceptance` SS/OS
`pass_muon_trk_charge` bin. Nominal `_no_res_cut` untouched (baseline = existing `*_template_fit.root` D_SS).

**V2 BACKGROUND MC — CORRECTED (2026-06-24, user).** My earlier "no background MC" was WRONG: I conflated "no
pre-FILLED background histograms" with "no background sample." The fullsim/overlay reco muon collection IS
skimmed exactly like data → it contains fake + hadronic-background muons; only the EXISTING histos are
signal-truth-matched. Starting from RECO muons (not truth-seeded), the truth-match branches classify provenance.
The current `PythiaFullSimExtras.c` fill is TRUTH-SEEDED (loops `truth_muon_list` pairs, lines 210-224) → never
pairs fakes/hadronic; the test needs a NEW **reco-seeded** fill.
- **Provenance classifier (VERIFIED on the actual NTUP, 2026-06-24).** Per reco muon, branches `muon_truth_prob`,
  `muon_truth_id` (pdgId), `muon_truth_IsPrimary`. NOTE: the `barcode>200k`=Geant4 convention is OVERLAY/HIJING
  only — in pythia fullsim 0% of matched muons have bc>200k, so use pdgId+IsPrimary instead:
  * **fake** = `muon_truth_prob ≤ 0.5` (no truth match)
  * **hadronic** = prob>0.5 AND ( `|muon_truth_id|≠13` [punch-through hadron, e.g. K±=321] OR (`|id|==13` AND
    `IsPrimary==0`) [π/K decay-in-flight muon] )
  * **prompt** = prob>0.5 AND `|id|==13` AND `IsPrimary==1`
- **Breakdown (pTH8–14 file, 20175 reco muons):** prompt 97.8%, fake 1.6%, punch-through 0.6%, decay-in-flight
  ~0%. ⇒ background exists + is labelable, but ~2% in this DiMu-filtered signal sample → LOW stats for low-mass
  background PAIRS. Sum all pTH slices first; if still marginal → **POWHEG Run2 pp17 fullsim** (dcache backup; valid
  for a physics test though NOT for detector-condition procedures like reco-eff/unfolding — pp condition suffices,
  background seen in both pp & PbPb). Data-magnitude won't be reproduced by signal MC (data SS 40% near-side vs MC
  2.7%) → use overlay for magnitude; signal-MC tests the PRINCIPLE (charge symmetry) + PROVENANCE.
- **The test (user, valid):** reco-seeded OS/SS low-mass minv split by provenance → (1) recover the SS low-mass
  shape [qualitative in signal MC], (2) confirm a genuine background with ≥1 fake/hadronic muon despite passing
  ALL muon-quality cuts, (3) SS↔OS mirror (charge symmetry → OS−SS valid).
- **RESULT (2026-06-24, `bkg_mc_provenance.C` on pythia signal fullsim + HIJING overlay; `/review-analysis-code`
  PASS).** macro: `dimuon_data/plots/template_fitting/bkg_mc_provenance_20260624/`. Check (2) CONFIRMED: genuine
  fake+hadronic background exists and SURVIVES all muon-quality cuts (signal MC: prompt 99.3%/fake 0.02%/had
  0.66%; overlay post-selection prompt 96.6%/fake 0.38%/had 3.05% — note the UE's 34.8% raw fakes are mostly
  removed by medium-quality+dP/p+d0/z0, leaving hadronic-dominated bkg). **Check (3) — the mirror is OPENING-ANGLE
  DEPENDENT:** overlay background SS/OS ≈ **0.98 at wide-angle/high-mass (minv>4)** [charge-symmetric uncorrelated
  combinatoric — mirrors ✓] but ≈ **0.20–0.29 at near-side/low-mass (minv<1.5)** [OS-ENHANCED]. Same in signal MC
  (low-mass bkg SS/OS≈0.29≈k). PHYSICS: the near-side low-mass background = (prompt signal muon + near-side
  hadronic from the same HF jet) → CORRELATED, OS-enhanced, ≈ the g→QQ̄ G template (k≈0.3); the symmetric
  combinatoric lives at wide angle. **IMPLICATION:** in MC the near-side low-mass background (analog of the data
  SS peak) is OS-correlated, NOT charge-symmetric → it does NOT cleanly mirror. This is CONSISTENT with the
  adopted OS−SS + MC(S+G) method (OS−SS cancels the symmetric part; the OS-enhanced (1−k)·G residual is removed
  by the MC G template) but COMPLICATES the simple "SS peak mirrors OS" claim — the near-side part does not.
  **CAVEAT:** these are hard-scatter SIGNAL samples (every event has a prompt muon → near-side bkg biased to
  prompt+hadronic); data near-side may differ (UE-only could be more symmetric). Robust data charge-symmetry
  evidence remains T_mix (SS/OS=0.999) — but that is wide-angle/uncorrelated by construction, NOT near-side.
  **OPEN for user:** whether the data near-side SS peak is truly charge-symmetric is NOT settled by this MC; flags
  a real consideration for the OS−SS charge-symmetry assumption (possible next looks: fake-fake/UE-only component;
  min-bias sample; or accept that OS−SS + G template absorbs the residual). Plots pending /review-plot.
- **USER DECISION (2026-06-24): PROBE the near-side charge symmetry further**, with rigorous per-stage documentation
  (Hypothesis -> Methodology -> Observation -> Verdict[support/challenge/rule out]) tracked here AND saved to a results
  dir with a clarifying doc. Structured investigation doc = `dimuon_data/plots/template_fitting/bkg_mc_provenance_20260624/SUMMARY.md`.
  **Stage 1 (DONE)** = the result above (near-side OS-enhanced, confounded by prompt-muon bias of signal samples).
  **Stage 2 (DONE, /review-analysis-code PASS)** = prompt-content decomposition. RESULT: overlay background = 98%
  prompt_plus_bkg (real muon + near-side jet hadronic), which CARRIES the OS-enhancement (near-side SS/OS=0.198);
  the genuinely-combinatoric bkg_bkg is too rare under pair_pt>8 (33 overlay / ~15 signal) → its charge symmetry is
  STATISTICALLY INCONCLUSIVE. **KEY:** with pair_pt>8 the MC low-mass background is dominated by (real μ + π/K
  hadronic FAKE) — OS-enhanced (does NOT fully cancel in OS−SS) AND not the g→QQ̄ G template (two real HF muons) →
  possibly a background the S+G fit doesn't model; this IS the Δp/p fake-muon population (connects low-mass ↔ Δp/p,
  per user). **Stage 3 options** (need user steer): (3a) relax pair_pt>8 for the bkg_bkg charge-symmetry test
  (cheap stats); (3b) higher-stat UE/min-bias sample; (3c) characterize the (prompt+fake) background vs G-template /
  Δp/p; (3d) data-driven near-side composition. **Plots (4 requested)** ready to make from Stage 1+2 → /review-plot.
- **dp/p coupling (user):** the hadronic decay-in-flight background IS the Δp/p-fit target (ID-MS kink); this same
  background MC enables the Δp/p template fit (reco quantities, before unfolding, in the nominal pipeline). Priority:
  finish the current combinatoric+truth-bkg template fit FIRST; defer the dp/p-vs-RAA ordering decision until after
  this investigation. Data evidence still ready (D_OS 1.71e4 / D_SS 935; ρ/ω+φ peaks; T_mix SS/OS=0.999).
