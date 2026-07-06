# Low-mass template fit — NOMINAL background-subtraction method + fit build

**Parent:** low_mass_dimuon_template_fit.md (sub-doc A)
**Scope:** The CURRENT nominal background-subtraction method — OS−SS + MC(S+G) 2-template fit + §3h smeared resonance (J/ψ/ψ′) templates; per-coarse-R_AA-bin signal yield N_S; correction ordering (RECO-level fit, 2026-07-01 reversal); signal acceptance A_sig (§3f); wiring N_sig into crossx/R_AA replacing OS−SS (§3g); the `OSminusSS_SG_fit` build. This is the LIVE method doc.
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.

---

## Summary / current method

The nominal background subtraction is **OS−SS + a 2-template MC fit** `OS−SS = N_S·S + N_G·G_OS` (+ §3h smeared resonance templates for J/ψ/ψ′ leakage), per coarse R_AA bin, extracting the signal yield N_S. This REPLACES the earlier combined OS+SS / k-anchored fit (which was abandoned — see sub-doc B `tf_k_factor_mixed_event.md`) and REPLACES the provisional OS−SS-only used in crossx/R_AA (`raa_from_rdf_crossx.md`). The whole program runs at **RECO level** (2026-07-01 correction-ordering reversal): trigger-eff → template fit (reco) → extract N_S → reco-eff + unfolding on N_S after the fit → signal acceptance A_sig last → wire N_sig into crossx/R_AA.

---

## Design Decisions

- **BACKGROUND METHOD CHANGED: combined OS+SS/k-anchored fit ABANDONED → OS−SS + MC (S+G) fit (2026-06-24,
  user-approved after the root-cause investigation).** OLD approach (this doc's original §2): coupled OS+SS
  template fit with the combinatoric from mixed-event T_mix and the correlated-G normalization SS-anchored via
  k=G_SS/G_OS. NEW approach: **OS−SS** (data-driven, removes the charge-symmetric combinatoric INCLUDING the
  soft near-side part) **+ a 2-template MC fit `OS−SS = N_S·S + N_G·G_OS`** (+ §3h smeared resonance templates
  for J/ψ/ψ′ leakage), per coarse R_AA bin; extract N_S. REASON (root-cause investigation, /review-investigation
  PASS iter 3 + 2 independent T_mix verifications, R&O "ROOT CAUSE RESOLVED"): the SS combinatoric has a large
  SOFT NEAR-SIDE component (40% of low-mass SS pairs ΔR<0.3, peak minv 1.25, ~15× the MC correlated near-side
  rate) from within-jet muons (π/K/fakes) that NO template models — event-mixing is uncorrelated→wide, the MC
  g→QQ̄ is correctly wide-angle→hard. So the combined fit (needing a combinatoric template) cannot close
  (N_C→0/χ²=368), but OS−SS removes that combinatoric data-drivenly and the residual S+(1−k)·G fits with the two
  well-modeled, well-separated MC templates (validated: N_S=5710,N_G=1639 positive, G/(S+G)=0.22, χ²/ndf=72).
  CONSEQUENCE: the ScrambGen/T_mix mixed-event machinery is SUPERSEDED for the nominal background (keep the code,
  it's correct + may serve a systematic / the C_OS≈C_SS check); §2's k-anchoring is no longer the nominal. The
  `low_mass_template_calc` _no_res_cut D_OS/D_SS fills + the (pT,η) MC S/G templates REMAIN the inputs (now for
  OS−SS+S+G). **This supersedes Physics Procedure §2's combined-fit model; the authoritative method is now
  OS−SS + MC (S+G) + resonance templates.**

- **CORRECTION-ORDERING REVERSAL — minv fit at RECO, not truth (2026-07-01, user).** OLD (2026-06-22):
  the minv template fit is done at TRUTH level, after reco-eff + unfolding bring the data to truth
  (§3e prescribed trig(reco)→unfold→reco(truth)→fit). NEW: the fit is done at **RECO** level, after
  trigger-eff but BEFORE reco-eff + unfolding; reco-eff/unfolding correct the extracted signal yield
  N_S AFTER the fit; then acceptance. PHYSICS REASON (user): the mixed-event combinatoric template uses
  **real-data muons** (a reconstructed object, no truth analogue) and the fake/hadronic background has
  **no truth match** — neither can be represented at truth, so the fit MUST be at reco. Consistent with
  the Run 2 HF-muon Δp/p reco-level yield fit. **Code impact (DONE):** template-fit RDF fills
  (`low_mass_template_calc` block, PP+PbPb) now use a TRIGGER-ONLY weight (`crossx_weight_trig_only` /
  `weight_for_dsigma_trig_only` = 1/L·w_trig, dropped w_reco); nominal signal-region crossx/R_AA
  UNCHANGED (still reco+trig). Bkg_mc_provenance data-vs-Pythia comparison likewise uses trigger-only
  data (fullsim MC carries the same reco-eff → correcting data for reco-eff would invalidate the
  comparison). **Template impact (REMAINING):** the S/G MC templates move truth→fullsim reco minv;
  the truth `_sigsel` templates are kept only for the frame-insensitive k=G_SS/G_OS ratio study.
  Verified against Physics Procedure §3 lead + §3e (both updated to reco).

## Progress Log

- 2026-07-01 — **CORRECTION-ORDERING REVERSAL: minv template fit moved to RECO level (Task A, user
  HIGHEST priority).** Reversed the 2026-06-22 "fit at truth after eff+unfold" decision → the fit is at
  RECO, after trigger-eff but BEFORE reco-eff/unfolding (§3 lead, §3e, §4, Objective, Design Decisions all
  updated). Physics: mixed-event = real-data muons (reco, no truth); fakes have no truth match. **Code
  (`/review-analysis-code`-style review PASS, 0 issues):** `RDFBasedHistFillingPP.cxx` +
  `RDFBasedHistFillingPbPb.cxx` `low_mass_template_calc` blocks now weight the template-fit histograms
  with a TRIGGER-ONLY dsigma weight — PP `crossx_weight_trig_only = crossx_weight·w_trig`, PbPb
  `weight_for_dsigma_trig_only = weight·(1/L_year)·w_trig` — dropping w_reco (reco-eff Defines removed from
  the template lambdas). Nominal signal-region crossx/R_AA UNCHANGED (still `*_trig_corr` = reco+trig).
  ACLiC-clean both (exit 0). **Reran all 8 template-fit passes** (pp24 + pbpb23/24/25 data, + 4 mixed-event
  T_mix), all rc=0. Sanity (pp24): entries IDENTICAL to the reco+trig baseline (same pairs) — op 3382719,
  ss 516004 — while integrals dropped ~0.757× (op 1.713e4→1.296e4, ss 935.4→694.1; op_nosel 1.807e4→1.356e4,
  ss_nosel 973.8→716.6), i.e. exactly the removed w_reco≥1 inflation (⟨w_reco⟩≈1.32 pp). Docs swept:
  roadmap (2026-07-01 row + old-ordering row superseded), analysis_overview §6 (reco-level note). S/G MC
  templates truth→fullsim-reco flagged as Remaining Work. Committed as a single commit. NOT git-tracked:
  the refilled `_template_fit.root` data histos (data area).

## Results & Observations

- 2026-06-24 — **BUILD step 1 DONE: production OS−SS + MC(S+G) fitter (pp), per R_AA pT bin.**
  `OSminusSS_SG_fit_20260624/code/osss_sg_fit_pp.C`: per coarse pT bin, OS−SS = N_S·S + N_G·G_OS (NNLS,
  non-negative, [1.08,4], J/ψ[2.82,3.32]+ψ′[3.55,3.82] PEAKS masked — modeling them with fixed Gaussians
  failed, χ²→5919, since the sharp data peaks aren't Gaussian; the small smeared-tail leakage below 2.85 →
  §3h systematic). **Per-bin signal yields N_S (the R_AA inputs):** integrated 5716; pT[8,21]=5298;
  pT[21,56]=408; pT[56,150]=2.9 (falls steeply ~ cross-section). G/(S+G)=0.18–0.32 (sensible). χ²/ndf high at
  low pT (60, 58) — MC G/S shape vs data continuum at the %-level amplified by tiny dσ errors → a fit-model/
  template-shape SYSTEMATIC; the yields are physical (vs the combined fit's N_C=0). Plots
  `osss_sg_fit_pp_*.png`, numbers `numbers_pp.txt`. NOTE: plot title still says "+ resonances" (stale — now
  masked) → fix at /review-plot. **Remaining BUILD: /review-plot the fits → PbPb fitter (per centrality) →
  signal acceptance (N_sig=N_S/A_sig) → wire R_AA → systematics.**

- 2026-06-24 — **ROOT CAUSE RESOLVED — soft near-side combinatoric; robust path = OS−SS + MC S+G
  (/review-investigation PASS iter 3).**
  WHY the data is softer than every template: a **soft NEAR-SIDE same-sign combinatoric**. In the low-mass region
  (minv<4, pp24), **40% of SS data pairs are near-side (ΔR<0.3)**, peaking at minv 1.25 — ~15× the MC correlated
  near-side rate (2.7%). The MC g→QQ̄ template is correctly **WIDE-angle/hard** (⟨ΔR⟩=2.50; two separate HF
  hadrons → wide muons → high minv) — NOT mismodeled; the MC signal is collimated (⟨ΔR⟩=0.30, 64% near-side). The
  near-side soft SS combinatoric (within-jet π/K/fakes) is captured by NEITHER event-mixing (uncorrelated→wide)
  NOR the HF template, which is why the combined OS+SS/T_mix fit fails (no template models it). **ROBUST FIX
  (validated):** OS−SS removes the charge-symmetric combinatoric data-drivenly (incl. near-side); then
  **OS−SS = N_S·S + N_G·G_OS** fits with the two well-modeled, well-separated MC templates → N_S=5710, N_G=1639
  (both POSITIVE), G/(S+G)=0.22, χ²/ndf=72 (residual = smeared-J/ψ/ψ′ leakage → §3h resonance templates) — vs the
  combined fit's unphysical N_C=0. Deliverables: `rootcause_data_soft_vs_template_hard_20260624/`
  (summary.md, data_SS_minv_by_dR.png, OSminusSS_S_plus_G_fit.png).
  **STRATEGY DECISION (for the user):** abandon the combined-fit/T_mix path; adopt **OS−SS + 2-template MC
  (S + G_OS) fit + §3h resonance templates** (= the pre-approved MC-only-OS path, sharpened; connects to the
  existing provisional OS−SS — the refinement is the MC S+G removing the residual (1−k)·G). Then signal
  acceptance → wire R_AA. Caveats: MC truth low effective stats; OS−SS charge-symmetry → a systematic.

## Latest Stage

**▶▶ BUILD PHASE — APPROVED 2026-06-24 (user): OS−SS + MC (S+G) fit → acceptance → R_AA.** Gate resolved
(root-cause: soft near-side combinatoric → combined fit is the wrong tool; OS−SS+S+G validated). Plan (each
code step → /review-analysis-code or /review-plot; fits → /review-plot):
1. **§3h resonance templates** — extract the smeared OS-only J/ψ(3.097)/ψ′(3.686) (+ φ leakage) shapes from the
   `_no_res_cut` OS data peak regions, to model the leakage into the signal window (the χ²=72 residual).
2. **OS−SS + (S + G_OS + resonances) FITTER, per coarse R_AA bin** (pp: 3 pT groups; PbPb: per centrality, then
   pT): non-negative fit, extract N_S per bin. = the validated path made production + per-bin. → /review-plot.
3. **Signal acceptance** A_sig(pair pT, pair η) from the cutflow denominator (log pT bins; §3f, Design Decisions
   reco-eff/acceptance boundary) → N_sig = N_S/A_sig.
4. **Wire into crossx/R_AA** replacing the provisional OS−SS-only (raa_from_rdf_crossx.md §3g): N_sig is the
   background-subtracted signal yield feeding 1/L (pp) and crossx_factor (PbPb). Rerun R_AA.
5. **Systematics:** OS−SS charge-symmetry, MC G shape (cc̄:bb̄, FE:GS), resonance-leakage, fit model, fiducial.
Inputs READY: D_OS/D_SS (`_template_fit`), MC S/G 1D + (pT,η) templates, the validated OS−SS=S+G fit. The
ScrambGen/T_mix path is SUPERSEDED (code kept for a systematic / the C_OS≈C_SS check).

**Current BUILD status (as of 2026-07-01):**
- ✅ BUILD step 1 (pp OS−SS + S+G fitter, per R_AA pT bin) DONE — N_S extracted per bin (integrated 5716).
- ✅ Correction-ordering reversal (fit at RECO) DONE — template-fit RDF fills now trigger-only; all 8 passes reran.
- Remaining: /review-plot the pp fits (fix stale "+ resonances" title) → §3h resonance templates → PbPb fitter (per centrality) → signal acceptance (N_sig=N_S/A_sig) → wire into crossx/R_AA → systematics. Rebuild S/G MC templates at RECO level (fullsim reco minv) to match the reco-level fit (truth `_sigsel` superseded for the nominal fit, kept for the k=G_SS/G_OS ratio only). After-fit reco-eff + unfolding of N_S needs real ε_reco + unfolding (roadmap Q4).
