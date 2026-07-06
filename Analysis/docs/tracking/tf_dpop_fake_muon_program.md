# Low-mass template fit — Δp/p fake-muon fit-and-subtract program (D0–D5)

**Parent:** low_mass_dimuon_template_fit.md (sub-doc D)
**Scope:** The Δp/p (= ρ = (p_ID − p_MS)/p_ID) fake-muon FIT-AND-SUBTRACT yield program: muon-level Δp/p yield extraction (Run 2 HF-muon method) propagated to the pair-level low-mass template fit and to final R_AA, with the double-subtraction and sufficiency logic. This is the DOWNSTREAM fit program — DISTINCT from the upfront fake-reduction cut study in sub-doc E `tf_upfront_bkg_reduction.md`, which it should be cross-referenced against. Provenance/composition/WP context lives in sub-doc C `tf_bkg_composition_normalization.md`.
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.

---

## Summary / current status

Program follows V1/V2 (sub-doc C): the low-mass near-side combinatoric background is fake/hadronic (within-jet π/K), OS-enhanced, distinct from gluon-splitting G. The plan (D0–D5) does the muon-level Δp/p yield fit (real vs fake/hadronic), propagates per-muon P_real → pair weight → purified real–real sample, tests sufficiency in MC, and integrates without double-subtraction. **D1 DONE (2026-06-24)**: Δp/p templates built; real peaks ρ≈0, hadronic shifted +0.08–0.09; within-jet hadronic is significantly softer/less-separable in the realistic overlay ⇒ Δp/p likely UNDER-removes ⇒ a scoped pair-level residual-F template is LIKELY needed (quantify in D4). NEXT: D2 muon Δp/p yield fit on data.

---

## Latest Stage — Δp/p fake-muon program directive + plan

**▶▶▶ Δp/p FAKE-MUON PROGRAM + PROPAGATION TO THE PAIR FIT & R_AA — NEW DIRECTIVE 2026-06-24 (user).**
Follows V1/V2: V2 established (MC) that the low-mass near-side combinatoric background is **fake/hadronic
(within-jet π/K)**, OS-enhanced, distinct from gluon-splitting G (wide-angle, real μ). User directive: (1) confirm
the low-mass background is fake/hadronic-dominated (vs g→QQ̄), (2) detailed study of the fake contribution's MASS
distribution, (3) perform the **muon-level Δp/p template fit in RECO quantities** (Run 2 HF-muon method), (4) work
out how to **propagate** muon Δp/p → the pair-level low-mass template fit → final R_AA, (5) determine whether
extracting signal at the muon Δp/p level is **SUFFICIENT** to resolve the within-jet fake/hadronic background, or
whether a **pair-level** fake/hadronic template is ALSO needed — and if so, **prevent double background subtraction**.

**GROUNDING (KB + analysis_overview §6).** Δp/p ≡ ρ = (p_ID − p_MS)/p_ID (KB `concepts/muon_source_template_fits`,
`run2_hf_muon_raa/vn`): REAL muons → symmetric peak ρ≈0; π/K FAKE (decay-in-flight, punch-through) → broad,
shifted to **positive ρ** (MS momentum mismodeled). Run 2 HF-muon = **yield-extraction** multi-template fit
(signal/hadronic/fake) → HF-muon yield (a per-muon yield correction applied FIRST). Run 2 dimuon note = pair
significance (quadrature of the two muons), 3 templates (Sig–Sig/Sig–Bkg/Bkg–Bkg) → >98% purity, NO cut. The
user now wants the **yield-extraction (HF-muon) version at the muon level**, propagated to our pair fit.
analysis_overview §6 names the THREE low-mass backgrounds explicitly: combinatorial (SS), gluon-splitting G
(correlated, real μ), fake μ (π/K, Δp/p) — this program builds the third and integrates all three.

**DESIGN ANALYSIS (the propagation + sufficiency + double-subtraction logic — authoritative; promote to
Physics Procedure §3i).**
Low-mass OS taxonomy: **S** signal (real μμ, one b) · **G** gluon-splitting (real μμ, two HF hadrons; OS-corr,
wide-angle) · **F** fake/hadronic (≥1 π/K fake; within-jet, near-side, partly OS-enhanced per V2) · **C** real-real
combinatoric (two unrelated real μ; charge-symmetric). The two fits act at DIFFERENT levels on DIFFERENT
backgrounds: **muon Δp/p fit removes F** (only handle separating a real μ from a π/K — kinematics cannot); **pair
OS−SS + S+G fit removes C** (OS−SS) **and separates S from G** (templates).
- **PROPAGATION (muon→pair→R_AA), clean ordering = Δp/p FIRST (matches Run 2 HF-muon).** Per-muon real fraction
  P_real(ρ,pT,η) from the Δp/p fit → pair weight w = P_real(μ1)·P_real(μ2) (if the two muons' fake-ness is
  independent) → apply at muon/pair level to get a **real–real** purified sample BEFORE the pair fit. Then OS−SS
  removes C; the S+G template fit on the purified sample → N_S (fake-, combinatoric-, G-free). N_S → A_sig → R_AA.
- **DOUBLE-SUBTRACTION RISKS (must be prevented).** (i) OS−SS already removes the **charge-symmetric part of F**;
  Δp/p removes **all F** → if both applied naively, symmetric-F is subtracted TWICE. (ii) Adding a pair-level F
  template on top of the Δp/p removal subtracts F twice. ⇒ **SCOPE each step:** apply Δp/p FIRST (removes ALL F at
  muon level) → then OS−SS operates ONLY on the Δp/p-purified real–real sample (it now removes C only, no F left to
  double-count) → S+G fit, **NO pair-level F template** (it would double-subtract). This is the clean design IF
  Δp/p is sufficient.
- **SUFFICIENCY — NOT assumed; it is the TEST.** Per-muon Δp/p may UNDER-remove F when: (a) within-jet hadronic μ
  have a weaker/different Δp/p signature than inclusive fakes (softer, less MS mismodeling); (b) the per-muon→pair
  product weight assumes INDEPENDENT fake-ness, but a near-side pair (both in one jet) may have CORRELATED fake-ness
  → product mis-estimates the pair fake yield. If a residual F survives the muon Δp/p, a **pair-level residual-F
  template is needed — scoped to the RESIDUAL only** (F_pair − F_removed_by_dpp), never the full F, to avoid (ii).
- **SUFFICIENCY TEST (MC, decisive):** truth-label pairs real–real vs ≥1-fake/hadronic (provenance machinery, DONE);
  build muon Δp/p templates real vs fake/hadronic (truth-labeled); apply the per-muon Δp/p purity (product weight)
  to the pair sample; compare the purified pair yield to the truth real–real yield. MATCH → Δp/p SUFFICIENT (no pair
  F template). Residual fake/hadronic excess → quantify it → pair-level residual-F template (scoped) required.

**CONFIRM-COMPOSITION (premise, step 0):** V2 already established (MC) the low-mass near-side COMBINATORIC/SS
background = fake/hadronic (within-jet π/K), NOT g→QQ̄ (which is wide-angle real-μ G, the S+G template's job).
Step 0 quantifies F vs G in the low-mass OS to confirm "fake/hadronic-dominated" for the combinatoric part.

**PLAN (each stage: Hypothesis→Methodology→Observation→Verdict; results saved to a dated dir with a doc):**
- **D0 — Confirm composition + fake MASS distribution** (existing MC provenance machinery + truth-flavor G): is the
  low-mass combinatoric background F-dominated vs G? Characterize the F minv shape (OS/SS, near-side). → analysis + /review-plot.
- **D1 — Muon Δp/p TEMPLATES** (real vs fake/hadronic) from fullsim/overlay (truth-labeled `muon_deltaP_overP`),
  inclusive and in the low-mass region; check separation power (fake shifted to positive ρ?) esp. within-jet. → /review-analysis-code + /review-plot.
  **DONE 2026-06-24 (/review-analysis-code PASS; `dpop_fake_muon_20260624/{code/dpop_templates.C, dpop_templates_*.root, SUMMARY.md}`).**
  D0 feasibility + D1: real peaks ρ≈0 (mean +0.014), hadronic shifted +0.08-0.09 (43% beyond cut) → fit feasible.
  **KEY for D4:** within-jet (low-mass-pair) hadronic Δp/p MATCHES inclusive in signal MC (KS=1.00) but is
  SIGNIFICANTLY SOFTER/less-separable in the OVERLAY (realistic UE: KS=0.000, mean +0.027 vs +0.077 inclusive) ⇒
  (i) the D2 fit must use the WITHIN-JET template, (ii) weaker within-jet separation ⇒ Δp/p likely UNDER-removes
  ⇒ a scoped pair-level residual-F template is LIKELY needed (quantify in D4). Caveat: templates pT-hat-unweighted
  (revisit weighting at D2). **NEXT: D2 muon Δp/p yield fit on data** (needs a relaxed-Δp/p data muon sample).
- **D2 — Muon Δp/p YIELD FIT (data)**: fit data muon Δp/p = f_real·T_real + f_fake·T_fake per (pT,η) → P_real per muon. (Run 2 HF-muon yield method, reco.) → /review-analysis-code + /review-plot.
- **D3 — Propagate to pairs**: per-muon P_real → pair weight w=P_real,1·P_real,2 → purified real–real sample. → /review-analysis-code.
- **D4 — SUFFICIENCY TEST (MC)**: does the propagated Δp/p purity reproduce the truth real–real pair yield? → decides whether a pair-level residual-F template is needed (scoped, no double-subtract). → /review-investigation + /review-plot.
- **D5 — Integrate**: Δp/p-purified → OS−SS → S+G fit → N_S → A_sig → R_AA (wire in, no double subtraction). → /review-analysis-code + /review-plot.
**Tracking convention (user): each stage = Hypothesis/Methodology/Observation/Verdict, saved to a dated dir with a
clarifying doc** (as the V2 `bkg_mc_provenance_20260624/SUMMARY.md`). Δp/p program dir: `dimuon_data/plots/template_fitting/dpop_fake_muon_<date>/`.
