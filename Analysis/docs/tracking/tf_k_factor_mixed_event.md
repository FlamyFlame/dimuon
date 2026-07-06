# Low-mass template fit — OS→SS factor k, mixed-event/ScrambGen, and the abandoned combined OS+SS fit

**Parent:** low_mass_dimuon_template_fit.md (sub-doc B)
**Scope:** The OS→SS correlated-background factor k study (5a MC-truth + 5b data closure), the mixed-event / ScrambGen combinatoric machinery (T_mix), and the ABANDONED combined OS+SS / k-anchored fit. This method is SUPERSEDED for the nominal by sub-doc A `tf_nominal_fit_build.md` (OS−SS + MC S+G); the code and results are retained here as history and as a candidate systematic / the C_OS≈C_SS charge-symmetry check.
**Reviewer rules:** RDF/C++ code → `/review-analysis-code`; plots/fits → `/review-plot`.

---

## Summary / current status

The combined OS+SS coupled fit (combinatoric from mixed-event T_mix + correlated-G normalization SS-anchored via k=G_SS/G_OS) was the original nominal plan. The **5b data closure FAILED** (2026-06-23, /review-investigation PASS iter 2): the data SS combinatoric is SOFTER than both T_mix and the MC G_SS, so no non-negative combination reproduces D_SS (N_C→0, χ²/ndf=368). Two independent neutral reviewers (2026-06-24) confirmed T_mix is NOT buggy and IS a valid strictly-uncorrelated combinatoric estimator. The root cause (soft near-side combinatoric) was resolved in sub-doc A: the combined fit is the wrong tool; **OS−SS + MC(S+G) is now nominal**. The k-validation infrastructure (5a MC-truth k, favorable) and the ScrambGen/T_mix mixed-event machinery are complete and correct, and are retained for a possible systematic / charge-symmetry cross-check.

---

## Physics Procedure §2 — top-level model + k (COPIED for reference; §2 also remains in the main doc's Physics Procedure)

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

---

## Design Decisions

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

- **Gate-driven autonomy (2026-06-23, user directive 2).** The template-fit mode is built
  end-to-end to produce per-R_AA-bin signal yields wired straight into the R_AA inputs, BUT
  the **k-validation + a closure plot** (Step 5) are delivered as intermediate results and
  reviewed with `/review-plot` using physics. The GATE is the **major combined-fit
  assumption**: the SS correlated-physics background (GS + FE) equals the OS correlated-physics
  background times a factor, i.e. `G_SS = k·G_OS` (and the full `SS = C + k·G_OS`). **If
  validated → proceed autonomously and wire into R_AA without waiting for review.** **If it
  fails → STOP/BLOCK and wait for the user's decision** (the fallback the user named: switch
  to **MC-only templates with SEPARATE fits for SS and OS**, abandoning the coupled OS+SS fit).

- **Gate criterion CONFIRMED = k(m,pT,η) function, NOT constant scalar (2026-06-23, user).** The
  coupled OS+SS fit uses the SS and OS correlated backgrounds' OWN MC shapes (sign1/sign2 templates)
  tied by a per-bin normalization, so a mass-dependent k is EXPECTED — NOT a failure. The gate is
  therefore: (a) is the MC-predicted ratio ROBUST against the cc̄:bb̄ theory uncertainty, and (b)
  does it survive the **5b data closure**. **Focus now on k(m,pT,η): if it passes the closure test,
  proceed to Part 2c.** The first-look constant-scalar FAIL (k(m) 0.14→0.38) is retained only as the
  justification artifact (see "Save validation results" below), not as a gate verdict.

- **TWO candidate code sets for the FINAL background subtraction — consider, do NOT implement B yet
  (2026-06-23, user).** Because `G_SS(m)` does NOT share `G_OS(m)`'s shape:
  - **(A) SS+OP COMBINED fit** with k(m,pT,η) from MC — the current path (build this).
  - **(B) MC-ONLY fit** (NO SS+OP combined) — separate SS/OS fits with pure-MC templates. **NOT
    implemented yet; a possibility to keep in mind.** The deeper open question it addresses: *if the
    SS shape does not trace the OS shape, is SS a valid reference for OS at all, and does the combined
    fit actually buy us anything against the large g→QQ̄ theory uncertainty?* Revisit (A)-vs-(B) after
    the 5b closure result; if closure is marginal or the cc̄:bb̄ sensitivity is large, (B) becomes the
    serious alternative. Keep the two as SEPARATE code sets so either can be the nominal.

- **Save validation results to a dedicated, non-overwriting folder (2026-06-23, user).** All k-validation
  (5a/5b) deliverables — plots (NOT overwritten), backed-up analysis code, numerical values, and a
  written summary doc — live together under `dimuon_data/plots/template_fitting/` (dated subfolders per
  snapshot). This preserves the evidence justifying the constant-scalar-fails finding and the k(m,pT,η)
  characterization.

- **ScrambGen rewrite — input/output formats + 5 open-question resolutions (2026-06-23, from the read-only
  scoping subagent `_sub_scrambgen_rewrite_scope_1.md`, merged + verified).** INPUT: single-muon trees =
  one TTree `muon_tree`, one branch `MuonObj` (type `MuonPbPb` PbPb / `MuonPP` pp), **one entry = one muon**
  carrying its own `ev_num`, `pt/eta/phi/charge`, `ev_centrality` (percent 0–84, −1 if >85%/invalid),
  `ev_FCal_Et`. The per-ctr `muon_tree_ctr{N}` trees were NOT written → the rewrite **bins muons itself** by
  `ev_centrality`. The `_hadd.root` files are empty (hadd needs the dict) → read parts via `TChain("muon_tree")`.
  OUTPUT: the RDF reads `muon_pair_tree_sign1`(SS)/`muon_pair_tree_sign2`(OS), branch `MuonPairObj`
  (`MuonPairPbPb`/`MuonPairPP`). Construct a pair: `p.m1=muonA; p.m2=muonB; p.year=yr; p.weight=1.0; p.Update();`
  — `Update()` computes minv (m_μ=0.105658), pair_pt/eta, dr, same_sign, avg_centrality. Sign: same_sign→sign1,
  else sign2. ScrambGen must include `DataAnalysisClasses.h` (object-branch dict). **RESOLUTIONS:**
  (1) **nScramb** = `5 × N_muons(interval)` per 5% interval (shape-only, N_C floats — no stale tables);
  (2) **ev_centrality<0 muons EXCLUDED** (outside 0–80%); (3) **spill DROPPED** — strict same 5% interval
  (`nCtrIntvls=20`, `CtrStep=5`, idx=ev_centrality/5); (4) **avg_centrality set EXPLICITLY** post-`Update()`
  to `(m1.ev_centrality+m2.ev_centrality)/2` (year-independent; bypasses the ambiguous yr25 single-event FCal
  recalc; both muons same interval so the average is well-defined); (5) **RDF mixed-event IO hook** (read the
  scrambled file → fill `T_mix`) is a SEPARATE step after ScrambGen. Mixing core: pick i; pick j≠i with
  `ev_num[j]≠ev_num[i]` (no same-event pairs); NO resonance/photoprod/dR cuts (now in the RDF signal_cuts).
  Output `muon_pairs_pbpb_20YY_single_mu4_scrambled.root` / `muon_pairs_pp_2024_2mu4_scrambled.root`. Per-year
  (combine downstream); run locally (in-memory, fast). PP: single bin, `MuonPP`/`MuonPairPP`, no year/avg_centrality.

## Progress Log

- 2026-06-23 — **T_mix fill DONE: mixed-event combinatoric template** — `/review-analysis-code` PASS iter 1
  (log `review-analysis-code-20260623-223706-mixed-event-tmix-fill.md`; all numbers MATCH). Added PUBLIC
  `mixed_event_template` sub-flag of `low_mass_template_calc`: when both set, SetIOPathsHook reads
  `muon_pairs_*_scrambled.root` (ScrambGen output) and TriggerModeSettings writes a distinct
  `*_template_fit_mixed_event.root`; the fill block is UNCHANGED (same 0-4 GeV OS+SS minv 1D+2D, PbPb per-ctr,
  efficiency-weighted — scrambled weight=1 → lumi·w_reco·w_trig, the 1/L scale absorbed by the floating N_C).
  New run scripts `run_template_fit_mixed_*.sh`. Ran all 4 → T_mix produced. **Verified:** OS minv SMOOTH, NO
  resonance peaks (pp J/ψ-ratio 0.99, pbpb23 ctr0_5 1.00 — confirms it read the scrambled input, not _no_res_cut);
  **SS≈OS** (pp SS/OS=0.999, area-norm shape dev 0.002; pbpb SS/OS≈1.03 per ctr) → **C_OS≈C_SS** (the 5b
  charge-symmetry input). Nominal + _no_res_cut template outputs UNTOUCHED. Stale "_no_res_cut" cout/run-script
  labels fixed. **T_mix ready → 5b coupled-fit closure (the gate).**
- 2026-06-23 — **ScrambGen REWRITTEN to the object model + scrambled muon_pairs produced** —
  `/review-analysis-code` PASS iter 2 (log `review-analysis-code-20260623-202339-scrambgen-objectmodel-rewrite.md`;
  all numbers MATCH). Replaced the 4 old non-compiling `ScrambGen/{ScrambGen,ScrambGenPP}.{c,h}` (retired
  `class Muon`) with header-only object-model classes: read `MuonObj` single-muon trees via TChain, bin by
  ev_centrality into 20 5% intervals (exclude <0), mix two muons from DIFFERENT events (`ev_num` differs) same
  interval, build `MuonPairPbPb` (m1/m2/year/weight=1/`Update()`), route `same_sign`→sign1(SS)/else sign2(OS),
  write `muon_pair_tree_sign1/sign2` branch `MuonPairObj`. NO physics cuts (RDF applies signal_cuts).
  oversample=5, fixed seed. New `run_scrambgen.sh`. Produced `muon_pairs_pbpb_20YY_single_mu4_scrambled.root`
  (1.9/1.4/3.7 GB) + `muon_pairs_pp_2024_2mu4_scrambled.root` (1.7 GB). **Verified:** SS≈OS (charge-symmetric
  combinatoric: yr23 7.046M/7.049M etc.); **OS minv = smooth continuum, NO J/ψ/φ peaks** (ratio 0.997/0.968 —
  the defining mixed-event correctness check); per-ctr population central-weighted.
  **BUG found+fixed (iter 1 CRITICAL):** yr25 `Update()` recomputes centrality from the pair `FCal_Et` (a mixed
  pair has none → all yr25 avg_centrality=−1, would be dropped by the RDF filter). Fix: set
  `avg_centrality=(m1.ev_centrality+m2.ev_centrality)/2` explicitly after `Update()` (Design-Decision
  resolution #4; no-op for yr23/24). Re-ran: yr25 mean=15.22, 0% at −1. **Open INFO for user:** `MuonPairPbPb.h:51`
  yr25 "centrality all-zeros→recalc-from-FCal" premise looks STALE (yr25 single-muon ev_centrality now valid,
  mean 15.48); nominal real pairs unaffected (real FCal valid), but worth revisiting globally later.
  **Scrambled muon_pairs ready → next: RDF mixed-event mode to fill T_mix → 5b closure.**
- 2026-06-23 — **5a folder + code RENAMED (user, descriptive naming).** `k_validation_5a_20260623/` →
  `OS_to_SS_factor_validation_MC_truth_constant_k_20260623/`; macros `k_validation_5a_{minv,ptEta}.C` →
  `OS_to_SS_factor_MC_truth_{minv,ptEta}.C` (+ matching function names, OUT paths, header comments, cout tags);
  `summary.md` retitled + macro refs updated; this doc's references updated. Re-ran both macros → identical plots
  (k_int=0.307816), old dir removed. Name reflects: MC-truth validation of the OS→SS correlated-bkg factor k
  (headline: constant/m-independent k FAILS, smooth k(m,pT,η) holds). The data-closure counterpart stays "5b".
- 2026-06-23 — **5a (pT,η)-dependence DONE (k stable per R_AA bin)** — `/review-plot` PASS iter 1 (log
  `review-plot-20260623-...-k-validation-5a-ptEta.md`; 6 re-extractions MATCH). From the 2D templates:
  **k vs pair pT stable ~0.30–0.32 (8–40 GeV) = k_int**, mild modelable rise to ~0.35–0.41 (46–68 GeV,
  bb̄ hardening), large low-stat errors >68 GeV; **k vs pair η flat ~0.31, symmetric** (outermost |η|>2.2
  bins ~0.48–0.50 = forward-acceptance-edge, low-stat, NOT instability); k(m) shape consistent across pT
  slices. ⇒ the per-R_AA-bin normalization k is STABLE (the m-dependence is the template shape) — favorable
  for the coupled fit. Plots `k_vs_pair_pt/eta.png`, `k_of_m_in_pt_slices.png` + `numbers_ptEta.txt` +
  summary in the dedicated folder. **5a (MC truth) COMPLETE; gate now hinges on the 5b DATA closure.**
- 2026-06-23 — **(pT,η)-binned truth fill DONE (k(m,pT,η) raw material)** + **single-muon Condor
  production SUBMITTED** + **ScrambGen safety verified.** (1) `/review-analysis-code` PASS iter 1 (log
  `review-analysis-code-20260623-175612-ptEta-binned-truth-kfill.md`): added 2D
  `{pair_pt_log_150,minv_zoomin}`+`{pair_eta,minv_zoomin}` fills (signal-selected, `_sigsel`) to
  `FillHistogramsTemplateMinvSignalRegion`; 1D templates unchanged (G_OS=2.2138, G_SS=0.6814); new histos
  `h_minv_zoomin_vs_{pair_pt_log_150,pair_eta}_sign{1,2}_<cat>_sigsel` (15×50 / 24×50); ACLiC-clean; pythia
  truth refilled (backup `*_pre_ptEta_kfill_backup.root`); sign mapping verified (k_bb=0.512, k_cc=0.009).
  (2) **ScrambGen single-muon-tree production safety VERIFIED** (background read-only subagent, scratch doc
  merged here + deleted per directive 3): all 3 incident fixes CONFIRMED in committed source (git 47f81f0):
  `output_single_muon_tree` now PUBLIC (DimuonAlgCoreT.h:121 / DimuonDataAlgCoreT.h:27); scripts set
  `output_single_muon_tree=true`+`pbpb_run3_mu4_force_nominal=true` (pbpb) / `trigger_mode=3` (pp) →
  `trigger_effcy_calc=FALSE`; output base name `single_muon_trees` ≠ `muon_pairs`
  (DimuonDataAlgCoreT.c:410) ⇒ **clobber structurally impossible**. May-skim inputs ready (pbpb 23/24/25
  parts 4/2/6, pp24 12); existing single_muon_trees STALE (Oct-2025, old naming) → production needed.
  (3) **SUBMITTED** 24 jobs: `condor_submit run_{pbpb_23,pbpb_24,pbpb_25,pp_24}_output_single_muon_tree.sub`
  → clusters 49/50/51/52 (4/2/6/12 jobs). Monitor → sanity-check first job's `.out` (single_muon_trees path,
  no protected-member error) → hadd parts per dataset for ScrambGen.
- 2026-06-23 — **Step C / 5a m-dependence DONE (Task #3, partial): k(m) characterised, constant-scalar
  k disproven, robustness driver identified.** `/review-plot` PASS iter 1 (log
  `review-plot-20260623-173509-k-validation-5a-mdep.md`; all numbers independently re-extracted, MATCH).
  Deliverables saved to `dimuon_data/plots/template_fitting/OS_to_SS_factor_validation_MC_truth_constant_k_20260623/` (plots/, code/,
  numbers.txt, summary.md — the dedicated non-overwriting folder per user). **Findings:** k_int=0.308;
  **k_bb=0.512 and ~flat in m (the ROBUST piece: B⁰ mixing χ_d + direct/cascade combinatorics)**;
  **k_cc≈0.009 (no charge-flip, as expected)**; one_b_one_c genuinely empty under the signal selection
  (G≈bb+cc here); k_comb(m) rises 0.14→0.38 ENTIRELY because f_bb(m)=OS_bb/G_OS rises 0.31→0.78
  (cc̄ dilution at low mass). ⇒ constant-scalar k FAILS (SS/(k_int·G_OS) sweeps 0.45→1.25, ~3×); smooth
  k(m) valid & expected; robustness favorable (robust k_bb modulated by theory-sensitive f_bb, which the
  OS fit can itself constrain since bb̄/cc̄ have different minv shapes). **Remaining for the gate:**
  (pT,η)-binned k(m,pT,η) [needs a (pT,η)-binned truth fill, code → /review-analysis-code] + the 5b DATA
  closure [needs T2 data refill from `_no_res_cut` + ScrambGen mixed-event]. Not git-tracked (data area).

## Results & Observations

### k-VALIDATION (5a) — FIRST LOOK from existing integrated truth templates (2026-06-23) ###
Quick read of `G_SS/G_OS` from the pT/η-INTEGRATED `_sigsel` truth templates
(`histograms_pythia_5p36TeV_no_data_resonance_cuts.root`; G = flavor {bb,cc,one_b_one_c},
OS=sign2, SS=sign1). 50 bins, 0–4 GeV.
- **Integrated `k = G_SS/G_OS = 0.308`** — physically sensible (0<k<1; the bb̄ charge-flip
  fraction of the correlated background). For context: G_OS/S_OS = 0.176 (correlated bkg ≈ 18%
  of the OS single-b signal integral over 0–4 GeV).
- **`k(m)` is smooth but NOT constant:** rises monotonically ~0.14 (m≈0.2) → ~0.38 (m≈3.9),
  a ~2.7× variation. ⇒ `G_SS(m)` does NOT share a shape with `G_OS(m)`; a **single scalar k
  is NOT valid**, a smooth **k(m) function** trivially is. Physics: `k(m) ≈ k_bb·f_bb(m)` where
  f_bb = bb̄ fraction of G_OS(m); f_bb rises with minv (bb̄ heavier) while cc̄ (no charge-flip,
  k≈0) dominates low mass. So the m-dependence is DRIVEN by the cc̄:bb̄ mix — the theory-sensitive
  piece the doc flagged (§2). **Implication for the GATE:** the constant-scalar reading of "SS
  bkg = OS bkg × factor" FAILS; the gate hinges on whether a smooth, ROBUST k(m, pT, η) makes the
  combined fit usable (needs: bb̄-only vs cc̄-dilution decomposition for robustness; (pT,η)
  binning; and the 5b data closure). NOT yet a pass/fail — full 5a investigation pending.

### 5b CLOSURE — FIRST ATTEMPT (2026-06-23, pp24) — fit-setup issues, NOT a clean pass/fail ###
First coupled-fit closure macro (`closure_5b_OS_SS_kfactor_data_20260623/code/closure_5b_pp.C`): separate
linear-LS fits SS_data=N_C·T̂_mix+N_kG·Ĝ_SS over [1.08,4]; OS_data=N_S·Ŝ+N_C·T̂_mix+N_G·Ĝ_OS over [1.08,4]
masking J/ψ+ψ′. **Result (integrated):** k_data=N_kG/N_G=0.296±0.010 vs k_MC=0.308 (ratio 0.96 — PROMISING
on the physics) BUT **N_C came out NEGATIVE** (SS −562, OS −1353; unphysical) and χ²/ndf poor (OS=83). High-pT
bins erratic/degenerate. **Diagnosis (hypotheses):** (a) the [1.08,4] window removes the low-mass
combinatoric-dominated region that separates T_mix from G → degeneracy → large-cancellation (negative-N_C)
solution; (b) SEPARATE OS/SS fits give inconsistent N_C — should be a COUPLED fit with SHARED N_C anchored by
the clean, WIDE-range SS (SS has no resonances → fit [0.3,4]); (c) OS resonance leakage (smeared J/ψ tail
below the [2.85,3.3] mask) contaminates → needs wider mask or resonance templates (§3h). ⇒ NOT a closure
fail (k_data≈k_MC is encouraging); the fit MACHINERY needs development before the gate verdict. → /review-investigation.

### 5b CLOSURE — FINAL VERDICT: the GATE FAILS (2026-06-23, /review-investigation PASS iter 2) ###
After developing a proper COUPLED OS+SS fit (shared N_C, SS wide clean anchor [0.3,4], OS resonance-masked,
non-negative yields via NNLS), the closure does NOT close for pp24. **Decisive evidence (shape infeasibility):**
the data SS combinatoric D_SS peaks at m=**1.72 GeV then FALLS**, while BOTH background templates rise to high
mass (T_mix peaks 3.96, G_SS peaks 3.64) — **the data SS is SOFTER than both templates** (41/46 bins in [0.3,4]
lie outside the templates' [min,max] envelope). Therefore NO non-negative combination N_C·T_mix + N_kG·G_SS can
reproduce D_SS → the NNLS pins N_C=0 with χ²/ndf=368 (unphysical: SS is combinatoric-dominated yet the fit uses
zero combinatoric). Pearson: corr(D_SS,T_mix)=**0.40** (the mixed-event T_mix genuinely mismatches the data
combinatoric — it over-populates high mass, because the single-muon mixing pool is harder than the soft muons
in real low-mass b-decay pairs), corr(T_mix,G_SS)=0.89. Root causes: (1) mixed-event T_mix too hard; (2) data
softer than both bkg templates ⇒ the combined fit is infeasible as built. Deliverables in
`dimuon_data/plots/template_fitting/closure_5b_OS_SS_kfactor_data_20260623/` (numbers_pp{,_coupled}.txt;
DIAGNOSTIC_SS_shapes_degeneracy.png; closure_coupled overlays; macros). **GATE FAIL ⇒ STOPPED per the user's
gate condition.** Options for the user: **(A)** refine the mixed-event method (mix within pT/η classes so T_mix
matches the soft data combinatoric) then retry the combined fit; **(B)** switch to **code-set B = MC-only
separate SS/OS fits** (the pre-agreed fallback — avoids the data-driven combinatoric entirely); **(C)** keep
the combined fit with a large combinatoric-shape systematic. Recommendation: the T_mix mismatch is the
immediate blocker; (A) is the natural next try, but (B) is the robust fallback if (A) still struggles.
**[SUPERSEDED 2026-06-24 by the T_mix independent verification below — the "T_mix too hard → abandon" reading
was incomplete; T_mix is valid, the issue is upstream (MC correlated template / near-side).]**

### T_mix VALIDITY — INDEPENDENT VERIFICATION (2026-06-24, two neutrally-prompted reviewers) ###
Per the user directive (verify T_mix validity WITHOUT my conclusions before abandoning the combined fit):
- **Reviewer 1 (code correctness, neutral): VERDICT CORRECT — no bug.** Verified on the produced files: 0
  same-event pairs, no self-pairing, SS/OS sign routing correct (sign1=100% SS, sign2=100% OS, matches the RDF),
  pair kinematics correct, centrality binning correct (0 cross-interval pairs), yr25 avg_centrality override
  works, RNG fine, downstream fill correct + distinct output. Non-fatal caveats: ×5 oversampling-with-replacement
  UNDERSTATES template per-bin errors (matters only if they enter the fit); pp scrambled file unclean-close
  warning (data intact; re-confirm); cross-part ev_num collisions → harmless false rejections.
- **Reviewer 2 (physics validity, neutral): T_mix is a VALID strictly-UNCORRELATED estimator, near-side caveat.**
  Independently reproduced T_mix (1% agreement → genuine method property, not a code artifact). KEY: real
  same-event SS pairs in 0-4 GeV carry NEAR-SIDE (Δη,Δφ) COLLIMATION (⟨ΔR⟩ 0.375 real vs 0.454 mix) that mixing
  destroys → mixed pairs pushed HIGHER in mass (minv ∝ opening angle). This fully explains T_mix harder than
  D_SS (mean 2.63 vs 2.22; mix/data 0.46→1.89; SS shape corr 1-4 GeV = −0.46). Pool ~10% softer (sub-dominant).
  **R2's crucial interpretation:** if T_mix is the uncorrelated-ONLY component with the correlated near-side
  (g→QQ̄/FE) in SEPARATE templates, the T_mix–SS mismatch is EXPECTED, not fatal; if T_mix is meant to be all of
  SS, it's biased.

**SYNTHESIS (orchestrator) — NOT 100% sure the assumption is wrong ⇒ do NOT switch to MC-only yet.** T_mix is
NOT buggy (R1) and IS a valid uncorrelated combinatoric estimator (R2). My closure fit IS the "T_mix
uncorrelated + separate MC correlated G" design R2 describes, so T_mix≠D_SS is EXPECTED, not a failure. The fit
still gave N_C=0/χ²=368 because the data low-mass D_SS (peaks 1.72) is softer than BOTH T_mix (uncorrelated,
hard) AND the MC G_SS (correlated, peaks 3.64) — NEITHER template is soft enough. ⇒ the real issue is UPSTREAM
of T_mix: either (a) the MC correlated-HF template G is mismodeled (too hard vs the data's low-mass near-side
correlated bkg), or (b) a near-side combinatoric component missed by both mixing and the HF templates. This is
NOT cleanly "SS=C+k·G is wrong," and it would ALSO affect the MC-only fallback (same MC G template). ⇒ per the
user's "switch to MC-only only if 100% sure," I am NOT — escalated to the user; the next investigation is the
data-soft-vs-MC-hard correlated-background question, not T_mix.
[NOTE: the root cause was subsequently resolved — soft near-side combinatoric; the nominal path became OS−SS +
MC S+G. See sub-doc A `tf_nominal_fit_build.md` R&O "ROOT CAUSE RESOLVED".]

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

## Latest Stage (combined-fit / gate saga — historical)

**▶ RESUMED 2026-06-23 — autonomous build APPROVED with a GATE (user directives 1–3).**
Three directives this session (see Design Decisions for each):
1. **Remove the `_res_cut_v2` nominal fallback** (PP + PbPb data RDF `SetIOPathsHook`) — nominal
   requires V1 `_mindR_0_02.root`, else throw. → `/review-analysis-code`. (Task #1.)
2. **Build the template-fit mode end-to-end → per-R_AA-bin signal yields wired into R_AA**, but
   deliver **k-validation + closure plot** (Step 5) as intermediate `/review-plot` results. GATE =
   `G_SS = k·G_OS` (full `SS = C + k·G_OS`): PASS → proceed autonomously & wire into R_AA; FAIL →
   STOP for user (fallback = MC-only separate SS/OS fits). Template fit is part of the nominal/crossx
   pipeline (its changes always rerun crossx + R_AA). (Tasks #2–#4.)
3. **Subagent scratch-doc procedure** in force for all delegated work.

**Execution ordering this session (critical path to the GATE):**
- **Step A (Task #1):** `_res_cut_v2` nominal fallback removal → `/review-analysis-code` → recompile.
- **Step B (Task #2 = T2):** template-fit mode reads `_no_res_cut`, re-fills 0–4 GeV OS+SS data
  spectra `D_OS/D_SS` (superseded V1 versions had resonance holes) + per-R_AA-bin minv; distinct
  output file; pipeline-integrated. → `/review-analysis-code`, rerun (reads existing `_no_res_cut`
  hadded ntuples).
- **Step C (Task #3 = T5a, FAST path to the gate):** MC truth-level k = G_SS/G_OS from Pythia truth
  G categories (`_bb`,`_cc`,`_one_b_one_c`; FSR gs + ISR FE together) binned in (pair pT, pair η,
  minv). Needs the truth-G minv ALSO (pT,η)-binned — check Step 4a output; add a binned truth fill
  if pT/η-integrated. Deliver k maps + projections → `/review-plot`. This is the structural MC read
  on whether `G_SS = k·G_OS` holds (necessary, not sufficient — 5b needed for the full ansatz).
- **Step D (long pole, parallel/background):** T1 single-muon-tree production (Condor) + ScrambGen
  object-model rewrite → mixed-event combinatoric template `T_mix`, prerequisite for **5b data
  closure** (`SS_data ?= C_mixed + k·G_OS`) and the coupled fitter. VERIFY the incident fixes
  (`output_single_muon_tree` public; `force_nominal`/`trigger_mode=3`) are in place BEFORE
  resubmitting any single-muon Condor fleet.
- **Step E (gate decision at 5b):** PASS → T4 resonance templates → T6 coupled fitter → T3 signal
  acceptance → T7 wire into crossx/R_AA. FAIL → STOP for user.

**PROGRESS 2026-06-23:** Step A DONE+committed (28757cd). Gate criterion CONFIRMED by user = k(m,pT,η)
function (not constant scalar); save validation to a dedicated folder; keep code-set (B) MC-only as a
future possibility (Design Decisions). **Step C / 5a m-dependence DONE** (/review-plot PASS; Progress Log
2026-06-23): k_bb robust ~0.5, k_cc≈0, k(m) rise = cc̄ dilution f_bb(m); constant-scalar disproven; smooth
k(m) valid; saved to `dimuon_data/plots/template_fitting/OS_to_SS_factor_validation_MC_truth_constant_k_20260623/`.
**USER GO-AHEAD 2026-06-23:** run the FULL chain end-to-end to R_AA with reviews, autonomously,
stopping ONLY if the 5b closure fails. Binning decided: pair pT = `pair_pt_log_150`
(=`pT_bins_150`, 15 log bins 8–150 GeV, matches the R_AA pT axis); pair η = `pair_eta` (24 bins).
k-fill = add 2D `{pair_pt_log_150, minv_zoomin}` + `{pair_eta, minv_zoomin}` (signal-selected, `_sigsel`)
to `FillHistogramsTemplateMinvSignalRegion`. EXECUTION (parallel tracks): FG = (pT,η) truth-fill code
(/review-analysis-code) → recompile → rerun pythia truth → k(pT,η) plots (/review-plot). BG =
read-only verification of the ScrambGen single-muon-production incident fixes (scratch doc
`_sub_scrambgen_singlemuon_verify_1.md`; NO submit/git by the subagent — orchestrator owns those).

**STATUS 2026-06-23 (checkpoint — all 5b prerequisites DONE):**
- ✅ **Part 1** (`28757cd`): nominal RDF requires V1; `_res_cut_v2` fallback removed.
- ✅ **5a gate (MC truth)** (`671af99`, `53ebd1c`, `ff6c52d`): k STABLE per R_AA bin (k_int=0.308; k_bb≈0.5
  robust; k_cc≈0; k(pT)/k(η) flat ~0.31). Constant-scalar disproven; smooth k(m) valid. Gate-FAVORABLE. 2×/review-plot PASS.
- ✅ **(pT,η) truth k-fill** (`53ebd1c`): 2D templates for k(m,pT,η).
- ✅ **Single-muon production**: 24 Condor jobs done; ALL 24 parts valid (no zombies), 19.72M muon entries
  (`single_muon_trees_{pbpb_20YY,pp_2024}_part*_*_mindR_0_02.root`). Incident-safe (distinct filename).
- ✅ **T2 data refill** (`37c5de6`): `low_mass_template_calc` mode → D_OS/D_SS from `_no_res_cut`
  (`*_template_fit.root`, 1D + 2D pair pT/η, PbPb per-ctr). OS resonances present; nominal untouched. PASS.

**ACTIVE 2026-06-23 (user): (i) RENAME the 5a folder + code to a descriptive name, then (ii) ScrambGen rewrite.**
(i) `k_validation_5a_20260623/` → `OS_to_SS_factor_validation_MC_truth_constant_k_20260623/` (descriptive: it is the
MC-truth validation of the OS→SS correlated-background factor k, headline result = constant/m-independent k FAILS,
smooth k(m,pT,η) holds). Rename macros `k_validation_5a_{minv,ptEta}.C` → `OS_to_SS_factor_MC_truth_{minv,ptEta}.C`
(+ matching function names + OUT paths); update summary.md + this doc's references; re-run to verify identical plots.
(ii) Then the ScrambGen rewrite (below), continuing to R_AA, stopping only if the 5b closure fails.

**REMAINING to the 5b GATE (Task #4, in progress):**
1. ✅ **ScrambGen object-model REWRITE DONE** (`c31fb8b`, /review-analysis-code PASS iter 2). Scrambled
   muon_pairs produced for pbpb 23/24/25 + pp24 (`muon_pairs_*_scrambled.root`); OS minv smooth, no resonance
   peaks (mixing correct); yr25 avg_centrality bug fixed.
2. **Fill `T_mix` via the RDF** — a mixed-event mode that reads `muon_pairs_*_scrambled.root` and fills the
   SAME 0-4 GeV OS+SS minv histos (1D + 2D pair pT/η, PbPb per-ctr) as the T2 D_OS/D_SS, efficiency-weighted
   (w_reco·w_trig; the 1/L scale is absorbed in the floating N_C). Cleanest: extend `low_mass_template_calc`
   with a `mixed_event_template` sub-flag (scrambled input + distinct output). → /review-analysis-code.
3. **5b DATA closure = THE GATE.** Per coarse R_AA (pT,η) bin: (a) charge-symmetry C_OS≈C_SS from T_mix
   (already SS≈OS globally); (b) the coupled-fit closure — does SS_data ≈ N_C·T_mix + (k·G_OS) hold (fit SS
   with {T_mix, G_SS,MC}; OS with {S_MC, T_mix, G_OS,MC}; check N_C consistent OS↔SS and N_G,SS ≈ k·N_G,OS).
   This overlaps the coupled fitter (Step 6). Deliver closure overlays → /review-investigation + /review-plot.
   PASS → Part 2c (resonance templates §3h → finalize coupled fitter → signal acceptance → wire R_AA)
   autonomously; FAIL → STOP for user (consider code-set B, MC-only separate fits).

**STANDING DIRECTIVE (user, 2026-06-23): DO NOT STOP for checkpoints any more. Autonomously continue
through ALL remaining steps until everything is done and R_AA is FINALIZED. The ONLY stop condition is the
5b closure FAILING (then stop for the user's decision per the gate). Otherwise: T_mix fill → 5b closure →
[on PASS] resonance templates §3h → coupled fitter → signal acceptance → wire crossx/R_AA → finalize R_AA.
Keep the tracking-doc protocol (plan before, results after) and the /review-* loops throughout, but do not
pause to report between steps.**

**▶ USER DIRECTIVE 2026-06-24 — VERIFY T_mix VALIDITY INDEPENDENTLY BEFORE ABANDONING.** Before concluding the
combined fit must be abandoned, make 100% sure the mixed-event T_mix procedure is VALID (not buggy). Get
SEPARATE reviewers to verify INDEPENDENTLY — **without being given my assumptions/conclusions** (no "too hard",
no gate-fail framing) to avoid confirmation bias. Keep ALL records (procedure + evidence of why the assumption
fails → why we abandon). Only if 100% sure the assumption is wrong → do MC-only OS-only. Track status in this
doc throughout. **Plan:** (i) spawn ≥2 independent, neutrally-prompted reviewers on the ScrambGen mixing code +
the T_mix RDF fill + the outputs vs data — "is this mixed-event combinatoric correctly implemented & a valid
combinatoric estimate? report bugs/concerns + verdict"; (ii) if a reviewer finds a BUG → fix + re-test the
closure; (iii) if independent reviewers CONFIRM the procedure is valid AND the closure still fails → archive
the failure (separate documented dir, like the constant-k folder) → switch to MC-only OS-only → proceed to R_AA.

**(prior) RESUMED 2026-06-23 — user decision on the gate-fail: (A) REFINE mixed-event T_mix (pT-class mixing so it
matches the soft data combinatoric) → RETRY the coupled closure. If the refined combined fit WORKS → proceed
with the combined fit autonomously to R_AA. If it STILL FAILS → (1) archive the combined-fit-failure results in
a SEPARATE documented directory (like the constant-k folder) with a doc on validations done + how it fails;
(2) switch to MC-only SEPARATE fits (likely OS-ONLY) and proceed with that to R_AA.** Plan: (i) diagnostic —
compare the single-muon pT pool vs the muon pT in real data pairs (confirm pT-class mixing is the right fix);
(ii) refine ScrambGen mixing to match the data muon-pT structure → re-fill T_mix; (iii) re-run the coupled
closure; (iv) branch per the result.
  - **DIAGNOSTIC RESULT (2026-06-24) — pT-class mixing is CONTRADICTED.** Single-muon pool mean pT = **7.27**
    GeV (4.2% > 15 GeV); muons in real data SS pairs (signal selection) mean pT = **8.12** GeV (8.7% > 15 GeV).
    The mixing pool is SOFTER than the data-pair muons, NOT harder. So matching the mixing to the data muon-pT
    would make T_mix HARDER (worse). ⇒ the T_mix hardness is NOT from the muon-pT pool; it is from the
    pair-building KINEMATICS (uncorrelated random-Δφ event mixing + the pair_pt>8 selection produce a harder
    minv than the real same-event combinatoric). Literal pT-class mixing will not fix it. A proper fix needs
    the mixed pairs to reproduce the data combinatoric's ANGULAR/pair kinematics (more involved, uncertain),
    OR go to code-set B. Surfaced to the user (the chosen fix is invalidated by the data).

**⛔ (superseded by the RESUMED block above) STOPPED AT THE GATE — 5b CLOSURE FAILED (2026-06-23).** Per the standing directive (only stop condition =
closure fail) and the Gate-driven-autonomy design decision, the autonomous run halts here for the USER's
decision. The combined OS+SS template fit does NOT close in data: the SS combinatoric (data) is SOFTER than
both the mixed-event T_mix (Pearson 0.40 vs data; over-populates high mass) and the MC G_SS — no non-negative
N_C·T_mix+N_kG·G_SS fits (N_C→0, χ²/ndf=368). Full verdict: R&O "5b CLOSURE — FINAL VERDICT" (/review-investigation
PASS iter 2). **USER OPTIONS:** (A) refine mixed-event (pT/η-class mixing to soften T_mix) + retry; (B) MC-only
separate fits = code-set B (pre-agreed fallback); (C) combined + large combinatoric-shape systematic.
**Everything upstream of the gate is DONE & committed and remains valid** (Part 1; 5a MC-truth k-validation —
favorable; (pT,η) truth k-fill; single-muon production; T2 data D_OS/D_SS; ScrambGen+T_mix). Acceptance (T3)
and the R_AA wiring (T7) are NOT started (they were gated on this). Resume once the user picks A/B/C.

**CHECKPOINT 2026-06-23 (ScrambGen done):** All infrastructure for the gate is now in place — data D_OS/D_SS
(`_no_res_cut`, T2), MC k-validation (5a, favorable), and the mixed-event combinatoric (ScrambGen). The
remaining gate machinery is the `T_mix` RDF fill (small) + the coupled OS+SS fit that IS the 5b closure test
(the decision point). Commits this session: 28757cd, 671af99, 53ebd1c, ff6c52d, e3b1c3a, 37c5de6, 06e5a63,
21547e9, c31fb8b.

**Next (remaining for the gate):**
1. **(pT,η)-binned truth fill** of the G categories (+ single_b) → /review-analysis-code → then k(m,pT,η)
   maps/projections → /review-plot (extends 5a to per-RAA-bin).
2. **T2 data refill from `_no_res_cut`** (template-fit mode) → /review-analysis-code (D_OS/D_SS for closure).
3. **ScrambGen mixed-event** (long pole: single-muon Condor production + object-model rewrite) — VERIFY the
   incident fixes first.
4. **5b DATA closure** `SS_data ?= C_mixed + k·G_OS` → /review-investigation + /review-plot = the FINAL gate.
   PASS → Part 2c (fitter → acceptance → wire R_AA) autonomously; FAIL → STOP for user (consider code-set B).

**RESOLUTION (see sub-doc A):** the gate saga was resolved 2026-06-24 — the root cause was a soft near-side
combinatoric (not a T_mix bug), the combined fit is the wrong tool, and the nominal method became OS−SS + MC(S+G).
This combined-fit / k-anchored path is SUPERSEDED for the nominal; the code and the favorable 5a MC-truth k
result are kept for a possible systematic and the C_OS≈C_SS charge-symmetry cross-check.
