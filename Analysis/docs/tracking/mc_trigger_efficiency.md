# MC-Based Trigger Efficiency (mu4 / 2mu4)

**Mode:** Implementation. **Created:** 2026-07-10. **Session:** "MC-based trigger efficiency".
**Siblings:** `mu4_trig_effcy_implementation.md` (CLOSED — data-derived chain; top-level
equations, tag-and-probe ε^nc, D9 bias lesson, cross-term §3c inherited from there),
`mc_trigger_info_skim.md` (ACTIVE — produced the trigger-enabled MC NTUPs this doc consumes),
`analysis_roadmap_2026_06.md` Q4 (ΔR trigger-correlation correction).

---

## Objective

Use the trigger-enabled Run-3 MC skims (pp24 Pythia fullsim `_July2026`; PbPb23-conditions
HIJING overlay r17618 `_July2026`) to deliver the MC-based trigger-efficiency program:

1. **MC single-muon mu4 efficiency** ε_MC(pT, q·η) = P[muon fires mu4 | offline reco muon],
   validated against the data-derived tag-and-probe efficiency.
2. **Factorization cross-check:** does the presence of a second offline muon at distance ΔR
   change a muon's single-mu4 probability? (ΔR-binned singles efficiency.)
3. **The ΔR correlation correction** ε_ΔR(ΔR) for 2mu4 (= the mu4 cross-term), via inverse
   weighting on an unbiased (no-trigger-requirement) MC pair sample. This is the analysis
   deliverable that replaces the current dummy ε_ΔR ≡ 1 (roadmap Q4).

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The analysis corrects yields by a per-pair trigger efficiency. The per-muon no-correlation
efficiency ε^nc(pT, q·η) is measured **from data** (tag-and-probe, kept as-is). The missing
ingredient is the **two-muon correlation**: whether the pair-level trigger probability
factorizes into the product/union of the two single-muon probabilities. In data this cannot
be measured without bias — the sample only exists because a muon trigger fired, so every
inverse-weighted denominator is depleted by P(at least one fires) = ε₁+ε₂−ε₁ε₂
(mu4_trig_effcy_implementation.md D9). The new MC skims store **every** event with the
trigger decision merely recorded (`StoreAllEvents` fix, mc_trigger_info_skim.md R1b), so MC
provides the unbiased sample data cannot.

MC's role is the **correlation (ratio) only**. The known MC L1 over-efficiency (per-leg
ε_L1^MC/ε_L1^data ≈ 1.13, largest in the turn-on; skim doc §8d) largely cancels in a ratio
but not in an absolute efficiency — so the data-derived ε^nc stays in the analysis, and Step 1
below is a *validation*, not a replacement.

### 2. Top-level equations

Per-pair trigger probability (inherited from mu4_trig_effcy_implementation.md):

- **PbPb, mu4 (union — at least one muon fires):**
  `P(pair | ΔR) = ε₁ + ε₂ − ε₁·ε₂·ε_ΔR^cross(ΔR)`
  where `ε_ΔR^cross(ΔR) = P(both fire | ΔR) / (ε₁ ε₂)` is the cross-term correction.
- **pp24, 2mu4 (product — both legs fire):**
  `P(pair | ΔR) = ε₁ · ε₂ · ε_ΔR^2mu4(ΔR)`
  — the 2mu4 correction IS the cross-term; no inclusion-exclusion.

with εᵢ = ε^nc(pTᵢ, q·ηᵢ), the single-muon mu4 efficiency of an isolated muon. The two
assumptions this doc measures/tests: (i) the singles terms ε₁, ε₂ carry **no** ΔR dependence
(Step 2 — critical for the PbPb union formula, whose linear terms cannot be absorbed into the
cross term); (ii) the correlation is a function of ΔR alone (checked via kinematic binning).

### 3. Step-by-step method

#### §3.1 MC single-muon mu4 efficiency + data comparison

- **What it measures:** ε_MC(pT, q·η) = P[muon is mu4-trigger-matched | offline reconstructed
  muon in MC]. Direct conditional probability — **no tag-and-probe needed**: the MC sample has
  no trigger selection, so every reconstructed muon enters the denominator unbiased.
- **Numerator condition:** the muon itself is matched to an HLT mu4 object (per-muon trigger
  matching branch), NOT the event-level chain decision (which the *other* muon could have fired).
- **Selection:** mirror the data-side measurement exactly — Tight WP, same fiducial cuts
  (pT > 4 GeV, |η| < 2.4), per-charge (μ⁺/μ⁻ separately; the toroid bends the charges
  oppositely, hence the q·η variable). MC weighted per pTHat slice (σ·ε_filt/N; isospin
  4:6:6:9 for pp fullsim).
- **Comparison target:** the data tag-and-probe ε^nc — P[2mu4 | offline pair, 1st muon passes
  mu4, ΔR > 0.8] in the 2nd muon's (pT, q·η) — pp24 data vs pp fullsim; PbPb data vs HIJING
  overlay.
- **Plots:** overlaid data-derived vs MC-derived 1D efficiencies in pT, η, φ (μ⁺ left panel,
  μ⁻ right panel); pT turn-on overlays in q·η bins.
- **Expected:** MC above data by ~10–15% per leg in the turn-on, converging on the plateau
  (L1 muon simulation lacks real RPC/TGC chamber inefficiency). The comparison establishes
  whether the residual is a smooth per-leg effect (cancels in the Step-3 ratio) or has
  localized (η, φ) structure (would need care).

#### §3.2 Factorization cross-check: ΔR-binned singles efficiency

- **What it tests:** the assumption that P[muon fires mu4] is a property of that muon's own
  (pT, q·η) only, unaffected by another offline muon at distance ΔR — i.e. that the ΔR
  correction lives ONLY in the cross/pair term.
- **Method:** for muons in MC **pairs** (no trigger requirement anywhere), measure
  P[muon mu4-matched | offline reco muon ∧ ∃ other offline reco muon at ΔR ∈ bin] in three
  ΔR bins: **0–0.2** (inside L1 RoI / MS-sector-sharing scale), **0.2–1.0** (intermediate),
  **1.0+** (isolated reference — must reproduce §3.1).
- **Why kinematic binning is essential:** small-ΔR pairs are kinematically special (boosted,
  higher pair pT), so an inclusive comparison confounds kinematic correlation with genuine
  trigger correlation. The comparison must be at fixed (pT, q·η).
- **Plots:** 1D pT and 1D q·η dependence, and pT in q·η bins — the 3 ΔR lines overlaid.
- **Interpretation:** lines agree at fixed kinematics ⇒ assumption holds; the PbPb union
  linear terms need no ΔR correction, and the data tag-and-probe ε^nc (measured at ΔR > 0.8)
  applies to close pairs. A small-ΔR deviation quantifies the L1 RoI-merging / sector-sharing
  effect on the singles level.

#### §3.3 ΔR correction via inverse weighting (the deliverable)

- **What it measures:** ε_ΔR(ΔR) = P(pair fires | ΔR) / (ε₁ ε₂), the genuine two-body trigger
  correlation. pp: numerator condition = pair passes 2mu4. PbPb: numerator condition = both
  muons mu4-matched (cross term).
- **Method (inverse weighting, mu4 doc §3c, now on an unbiased sample):**
  - Denominator: **all** MC reco pairs, no trigger requirement, unit weight.
  - Numerator: pairs satisfying the trigger condition, each weighted 1/(ε₁·ε₂) with
    εᵢ = **MC-derived** ε_MC(pTᵢ, q·ηᵢ) from §3.1, evaluated continuously at the exact
    (pT, q·η) (TF1 Eval, never resample-to-nearest).
  - Ratio vs ΔR = ε_ΔR(ΔR).
- **Why MC ε in the weights:** self-consistency. If the factorization were exact and the
  §3.1 fits perfect, the ratio would be 1 at every ΔR. Using data ε^nc would instead park the
  plateau at (ε_MC/ε_data)² ≈ 1.28, conflating per-leg normalization with correlation.
- **Diagnostics:** (1) **plateau at large ΔR** — well-separated muons occupy different
  detector regions, decisions must decorrelate; failure to flatten means residual kinematic
  mis-parameterization of ε_MC leaking in. (2) **plateau = 1** — an offset measures fit
  quality, not physics. The physical content is the small-ΔR shape relative to the plateau.
- **Application:** ε_ΔR multiplies ε₁ε₂ in the pp 2mu4 weight, and dresses the ε₁ε₂ cross
  term in the PbPb union weight (both data-derived ε's unchanged).

### 4. Negative constraints

- **NO trigger requirement on any denominator** (the whole point of the MC sample; D9 lesson).
  Any per-muon numerator uses that muon's own trigger match, never the event-level decision.
- **Do NOT replace the data-derived ε^nc in the analysis with ε_MC.** MC contributes the ΔR
  correlation ratio only; §3.1 is validation.
- **Do NOT weight §3.3 numerators with data-derived ε^nc** (breaks the plateau=1 diagnostic).
- **Do NOT read raw NTUPs standalone** — consume/extend ntuple-processing output per the
  provenance rule (new trigger branches → propagate via a processing mode/flag, distinct
  output suffix).
- This differs from **mu4_mu4noL1**: no requirement on the other muon anywhere in the singles
  efficiency (§3.1, §3.2).
- MC is always weighted (per-slice σ·ε_filt/N; pp isospin 4:6:6:9).
- **Menu caveat (systematic, not a bug):** overlay MC uses HLT menu `Dev_HI_run3_v1`
  (L1 `MC_HI_run3_v1`, TriggerValidation prescale set, all relevant chains PS=1) vs PbPb23
  data HLT `PhysicsP1_HI_run3_v1` (L1 `Physics_HI_run3_v1`, physics prescale sets) — data
  side verified 2026-07-10 from the local test AOD run 462240 (mc_trigger_info_skim.md R2).
  pp fullsim uses `PhysicsP1_pp_lowMu_run3_v1`. Record as a possible systematic on the ΔR
  correction.

## Context (condensed from siblings)

- **Samples (mc_trigger_info_skim.md, all validated 10 000 entries/NTUP, non-empty trigger
  branches):** pp24 fullsim = 24 NTUPs (4 isospin beams × 6 pTHat slices, DSIDs 802758–802781,
  r16578) in `~/usatlasdata/pythia_fullsim_test_sample/`; HIJING overlay = 6 NTUPs (r17618)
  in `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/`. Full-statistics productions
  still in the pipeline — this doc's machinery reruns unchanged when they land.
- **Available branches:** event-level `b_HLT_*`; per-muon matching `muon_b_<chain>` (+ `_0_01`,
  `_V2`, `_V3`); dimuon `dimuon_b_<chain>_<dR>` + 8 per-leg branches
  `dimuon_b_HLT_2mu4_L12MU3V_mu{1,2}passLeg{1,2}_dR_*`. In MC, 2mu4 ⊆ mu4 and
  mu4_mu4noL1 ⊆ mu4 hold exactly (skim doc §8b).
- **Data-derived ε^nc (mu4 doc, Pipeline 2):** tag-and-probe P(2mu4 | mu4-tag, ΔR > 0.8),
  role-swap combined, fitted Fermi+log per (centrality, μ±, fine q·η bin), TF1s in
  `single_mu_effcy_pT_fit.root` (PbPb 23/24/25 + pp24).
- **D9 bias (why data failed):** on a mu4-selected sample every inv-weighted denominator is
  depleted by ε₁+ε₂−ε₁ε₂ (~1.3 low pT → 1.0 high pT), faking ΔR structure through kinematics;
  the pair-level "fix" (D6) was tautological. Cross-term code (§3c) was kept explicitly as the
  reference for unbiased MC — this doc is that use case.
- **Known MC/data trigger offset (skim doc §8d):** per-leg L1 ratio ≈ 1.13 (turn-on-dominated);
  MC/data ≈ 1.24–1.44 for 2mu4-needing-two-RoIs, ≈ 1.0–1.1 for one-RoI noL1 chain. Benign for
  a ratio-based correction; quantify as systematic.
- Current analysis status: crossx weights use w_trig = 1/(ε₁+ε₂−ε₁ε₂) (PbPb) and 1/(ε₁ε₂)
  (pp) with ε_ΔR ≡ 1 dummy (roadmap Q4).

## Scope

**In:** the three measurements above on both trigger-enabled MC families (pp fullsim → 2mu4;
HIJING overlay → mu4 cross-term vs PbPb data), the data-vs-MC comparison plots, and delivery
of ε_ΔR(ΔR) in a form the crossx weights can consume. Ntuple-processing extension (mode/flag)
to propagate the trigger branches.

**Out:** re-deriving the data tag-and-probe ε^nc; applying ε_ΔR to crossx/R_AA (follow-up once
the correction is validated); the full-statistics MC productions (rerun when they land).

## Design Decisions

### D1: pp is the primary tester; PbPb (overlay) is implemented+run but secondary (2026-07-10, user)
**Physics:** Two known limitations degrade the overlay-vs-data comparison quality:
(a) the overlay MC simulates the `Dev_HI_run3_v1` HLT menu with validation prescales while
PbPb23 data ran `PhysicsP1_HI_run3_v1` with physics prescales (verified from files, skim doc
R2); (b) the overlay TEST sample was produced with a single fixed vertex position instead of
the Gaussian sampling over 4 positions used in data-like production — vertex-z affects muon
trajectories through the trigger geometry. pp fullsim has neither issue at this severity.
**Decision:** validate the procedure and the MC trigger simulation primarily on pp
(pp24 fullsim vs pp24 data); run the full HI chain too, but interpret slightly larger MC-data
differences there as expected, not as failures of the method.

### D2: Overlay test sample compares ONLY to PbPb23 data, 0–5% centrality (2026-07-10, user)
**Physics:** The overlay test sample is generated with PbPb 2023 conditions and impact
parameter b ∈ 0–5 fm only ⇒ the sample is dominated by 0–5% centrality events. Comparing to
any other year or centrality class would mix detector conditions and occupancy regimes that
the sample does not represent.
**Test-sample-only caveat:** the upcoming full HIJING overlay production will cover all
centrality bins, use PbPb 2024 conditions, and the data-overlay will cover year-to-year
(2023–2026) variations — the centrality/year restriction lifts then.

### D3: Proceed with the pp fullsim slices already re-skimmed; skip in-flight ones (2026-07-10, user)
2 of the 24 pp fullsim `_July2026` grid jobs were still running at start of work. Use the
completed slices (weights are per-slice, so a missing slice biases nothing at fixed pTHat;
it only thins statistics in its pT range). Check `mc_trigger_info_skim.md` later for
completion and re-run the (cheap) downstream stages when all 24 are in. Not highest priority.

## Implementation Plan

*Work on dedicated git branch. pp first (D1), then PbPb overlay.*

1. [x] **Discovery** (results in R1/R2): trigger-branch inventory; NTP fullsim path; data
   tag-and-probe outputs (pp24 + pbpb23 0–5%); pp `_July2026` completion (23/24).
2. [x] **NTuple-processing extension** (per provenance rule): `store_mc_trigger` flag
   propagating per-muon mu4 matching + pair-level 2mu4 decisions into the fullsim
   muon-pair and single-muon trees, `_mc_trig` output suffix. `/review-analysis-code`
   **PASS iter 1** (0 CRITICAL / 0 WARNING, 5 INFO; provenance 359/359 bit-exact vs raw).
3. [x] **Run NTP** (full stats, 2026-07-10, both samples): pp pairs SS 41 096 / OS 152 219
   (kin-sum == global ⇒ double-fill fix verified at full stats), pp singles 455 440;
   overlay pairs SS 10 636 / OS 39 211, overlay singles 95 410. Old `pp_pTH8_14` skipped
   (D3). Outputs `*_mc_trig*.root` in the two sample dirs.
4. [x] **Step 1 (§3.1):** MC singles ε vs data tag-and-probe, both samples
   (/review-analysis-code PASS iter 2; /review-plot PASS iter 2).
5. [x] **Step 2 (§3.2):** ΔR-binned singles factorization check, both samples (same reviews).
6. [x] **Step 3 (§3.3):** ε_ΔR^2mu4 (pp, plateau 0.963≈1 ✓, small-ΔR suppression) and
   ε_ΔR^cross (overlay, plateau 0.864 flat systematic → plateau-normalize before use,
   small-ΔR enhancement ~2.2× rel. plateau).
7. [x] **Overlay chain** delivered together with pp in steps 4–6 (D2 restriction applied).
8. [x] **Bookkeeping:** roadmap Q4 row updated (measured; preconditions for application
   listed); INDEX scope updated; branch `mc-trigger-efficiency` left UNMERGED for user
   review (user requested a dedicated branch).

## Progress Log

*(append-only; newest entries at the END)*

- 2026-07-10 — Doc created. Physics Procedure written from user specification + inheritance
  from `mu4_trig_effcy_implementation.md` (equations, §3c cross-term, D9) and
  `mc_trigger_info_skim.md` (samples, branches, R1b/R2/§8d). Presented to user for approval;
  no code yet.
- 2026-07-10 — Physics Procedure approved by user with constraints recorded as D1–D3
  (pp primary tester; overlay ↔ PbPb23 0–5% only; skip the 2 in-flight pp slices).
  Sibling alert: `hijing_overlay_det_response_band.md` (branch
  `fix/hijing-overlay-det-response`) has pending fixes to `PythiaFullSimExtras.c`
  (exclusive truth→reco matching; global-pair-tree double-fill) and the
  `hijing_overlay_pp24`→`hijing_overlay_pbpb23` naming change. Trigger-efficiency
  measurements here are reco-level ratios — unaffected by the truth-matching bug; the
  double-fill affects absolute yields of the global fullsim pair tree, so ratio
  measurements are safe, but per-pair trees used here must be checked for the double-fill
  (dedupe or use kin trees if absolute denominators matter — they don't for efficiencies,
  numerator and denominator double alike). Naming: this doc's overlay outputs adopt
  whatever `FullSimSampleType.h` provides at run time; coordinate at merge.
- 2026-07-10 — **Step 2 (NTP trigger extension) DONE, /review-analysis-code PASS iter 1**
  (log `.claude/logs/review-analysis-code-20260710-010710-ntp-mc-trigger-propagation.md`).
  Code: `store_mc_trigger` public flag (`PythiaAlgCoreT.h`); `_mc_trig` output suffix
  (OutputTreePathHook/OutputHistPathHook); trigger binds + per-file skip-if-absent
  (`PythiaFullSimExtras.c::InitInputExtra`, one-file-per-chain asserted); per-muon
  `m.passmu4 = muon_b_HLT_mu4_L1MU3V[reco_ind]` (bare branch = mindR 0.02 data-mirror, NO
  mu6/mu8); new `MuonFullsimExtra::reco_ind`; new `PairMCTrigExtras{pass2mu4, ev_pass_mu4,
  ev_pass_2mu4}` on both fullsim pair structs (**ROOT gotcha: a single-member base class is
  not member-split into named leaves — keep ≥2 members**); pair lookup via the skim's (i<j)
  `muon_pair_muon{1,2}_index` block, fail-fast throw; single-muon tree gate in trigger mode
  = reco-based (reco_match && pt>3 && |η|<2.6; exact fiducial in RDF) so the Step-1
  denominator is an offline-reco condition (truth gate would sculpt the turn-on). 4 runners
  `run_pythia_fullsim{,_overlay}{,_single_muon}_mc_trig.sh`.
  **Design note (reviewer INFO-2):** `passSeparated` NOT added to the pair struct — `dr` is
  already stored and §3.2/§3.3 bin ΔR directly (0–0.2/0.2–1.0/1.0+), so RDF defines it when
  needed.
  Smoke tests (300 ev/file): pp — old trigger-off `pp_pTH8_14` skipped loudly; OS pairs
  both-reco 4379, both-legs-mu4 2732, pass2mu4 2406 (the 2406<2732 gap = the pair-level
  correlation we will measure); ALL subset relations exact (2mu4⊆both, match⊆ev_mu4,
  ev_2mu4⊆ev_mu4: 0 violations). Singles: 13674 entries all reco-matched; unweighted turn-on
  0.564/0.723/0.855/0.882/0.871 (plateau ~0.87–0.88, Run-2-consistent). Overlay: 733
  both-reco pairs, avg centrality 2.06% (D2 confirmed). Provenance: 359/359 pairs bit-exact
  vs raw NTUP; C6 weight = σ·ε_filt·isospin/N exact.

- 2026-07-10 — **Step 3 (full NTP runs) DONE.** 4 sequential runs, all rc=0: pp pairs
  (SS 41 096, OS 152 219; kin-tree sum == global tree ⇒ the imported double-fill fix
  verified at full stats), pp singles (455 440 reco-matched muons), overlay pairs
  (SS 10 636, OS 39 211 — exactly the previously-halved kin-sum count), overlay singles
  (95 410). Old `pp_pTH8_14` skipped per D3. The single-muon-mode "cut-acceptance ZERO
  bin" caught-error is pre-existing (pair histograms unused in that mode).

- 2026-07-10 — **Step 4 MC side DONE (delegated executor; merged from `_sub_mctrig_mc_side.md`).**
  New `RDFBasedHistFilling/FillMCTrigEffHists.cxx` (`FillMCTrigEffHists(sample, do_step3)`;
  Step-1 singles + Step-2 ΔR-binned pair-leg hists → `mc_trig_eff_hists_<label>.root`, Step-3
  → separate `..._step3.root`) and `RDFBasedHistFilling/FitMCSinglesEffcy.cxx` (MC turn-on
  fits mirroring SingleMuEffcyPtTurnOnFitter exactly: pp erf+log / overlay fermi+log, "QR",
  [4,60], BayesDivide graphs; TF1 keys `f_mc_pt_vs_q_eta_<muplus|muminus>_<lo>_TO_<hi>` +
  fallback ratio TH2D, in `<dir>/single_mu_effcy_pT_fit_mc.root`). Data binning conventions
  replicated exactly (pt2nd = pT_bins_8+pT_bins_60 incl. the zero-width duplicated-8.0 edge;
  q·η = eta_bins_trig_effcy 171 bins; φ 128; pair_pt = pT_bins_120). All MC-weighted; overlay
  restricted to 0–5% (D2). **Headline numbers (full stats):**
  - §3.1 P(mu4|Tight+fiducial): pp μ⁺/μ⁻ 0.786/0.776 integrated, turn-on 0.62–0.65 → plateau
    0.90–0.91; overlay 0.629/0.627 → plateau 0.81–0.84.
  - Fits: 40/40 converged (status 0), χ²/ndf ~20–53/36.
  - **§3.3 ε_ΔR plateau (dR∈[1,3] avg): pp 0.963 (≈1 within 4% ✓); overlay 0.864** (~14% low —
    low stats + D1-expected overlay degradation). Small ΔR: pp 2mu4 SUPPRESSION 0.77–0.84 below
    0.2 recovering by ~0.35 (two distinct L1 RoIs required → merging kills 2mu4); overlay mu4
    cross-term ENHANCEMENT 1.95/1.43/1.31/1.05 in [0,0.05)/…/[0.15,0.2) (one RoI can match both
    offline muons → both legs flagged), statistically weak. Nothing tuned.
  - ε-evaluation bookkeeping: TF1 clamp[4,60]+floor 0.02 fired 0.053% (pp) / 0.27% (overlay);
    ~18% of legs in q·η gap regions use the unfitted 2D-ratio fallback (mirrors data
    EvaluateSingleMuonEffcyPtFitted).
  Note: the global pair trees are no longer double-filled (kin-sum == global verified, step 3),
  so the scratch doc's residual double-fill caveat is obsolete — ratios were immune either way.

- 2026-07-10 — **Step 4 data side DONE (delegated executor; merged from `_sub_mctrig_data_side.md`).**
  Probe 1D histograms added to the P2 tag-and-probe filling: `pt2nd`/`eta2nd`/`phi2nd`/`q_eta2nd`
  APPENDED to `single_muon_trig_effcy_var1Ds` (`RDFBasedHistFillingData.h:38` — note the `={}`
  at Data.cxx:62 is the isForSoumya branch, not the nominal default); `eta2nd_bins` registered
  (48 uniform [-2.4,2.4]; per-ctr PbPb variants 48/48/48/48/24/12); `eta2nd` entries added to
  var1D_pp.json + var1D_pbpb.json. P2 regenerated for pp24 + pbpb23 (backups `*.bak_20260710`).
  **Validation:** graph name sets identical to backups (60 pp / 540 pbpb) with bitwise-identical
  spot-checked values; 2D integrals exactly equal (pp `h_pt2nd_vs_q_eta2nd_mu4_sepr` 1 021 432;
  pbpb `_ctr0_5_sign1_mu4_sepr` 184 180); added keys exclusively the probe 1Ds (48 pp / 240
  pbpb). Integrated P(2mu4 | mu4-tag, sepr): **pp24 sign1 0.6836, pbpb23 ctr0_5 sign1 0.5040**.

- 2026-07-10 — **Step 4 review (/review-analysis-code, 2 reviewer subagents).** Data side:
  **PASS iter 1** (0C/0W/2 INFO) — content preservation verified exhaustively (all 252 pp +
  1672 pbpb pre-existing objects bitwise identical to backups incl. Sumw2 and graph error
  arrays); all numbers reproduced; new probe-η/φ efficiencies show the expected L1 structure
  (η≈0 crack, |η|≈1.1–1.2 transition, feet dips, 16-sector φ modulation). MC side: **FAIL
  iter 1 (1 WARNING, 5 INFO) → amended → see next entry.** All MC numbers reproduced to 6
  digits; §4 constraints, D2, C5, C6, fit-mirror exactness all verified. Corrections from
  the review: (a) the q·η axis is **184** bins (not 171 as previously logged; the bit-exact
  match to the data histogram edges is what matters and holds); (b) pp fit χ²/ndf range is
  20.1–**56.9**/36 (worst bin `muplus_minus2_40_TO_minus2_00`), still acceptable.
- 2026-07-10 — **WARNING resolved: overlay ε_ΔR^cross plateau 0.864 is a flat SYSTEMATIC
  offset, not low statistics.** Reviewer's per-bin check: all 8 bins in dR∈[1,3] sit below 1
  (0.82–0.96), offset ≈3–4σ even with conservative errors. Per §3.3 diagnostic (2) this is
  per-leg fit/parameterization quality (fermi+log residuals in the turn-on where the pairs
  live; φ/occupancy structure absent from the (pT, q·η) parameterization) — pre-registered
  by the Physics Procedure as a plateau-offset metric, so no C4 investigation. **Consequence
  recorded: ε_ΔR^cross must be PLATEAU-NORMALIZED (shape relative to the large-ΔR plateau)
  — or the per-leg parameterization improved — before it dresses the PbPb union cross term.**
  pp needs the same treatment in principle (plateau 0.963, a 4% offset).
  **INFO-6 verified (orchestrator):** the overlay small-ΔR enhancement is NOT a floor
  artifact — recomputed with floor 0.02→0.10 and with floored pairs dropped entirely:
  [0,0.05) 1.952→1.900/1.934, [0.05,0.1) 1.434→1.337/1.360, plateau unchanged (0.864→0.861/
  0.863). Relative to plateau, ε_ΔR^cross(dR<0.05)/plateau ≈ 2.2. INFO code fixes applied
  (output-TFile zombie checks; MCEffEvaluator throws on missing efficiency source instead of
  silently flooring; SS+OS summing intent commented); both macros recompile clean.

- 2026-07-10 — **Steps 4–7 plots DONE, /review-plot APPROVED iter 2** (log
  `.claude/logs/review-plot-20260710-034514-mc-trig-eff-plots.md`; executor scratch
  `_sub_mctrig_plots.md` merged here, then deleted). Macro
  `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff.cxx` (`plot_mc_trig_eff(sample,
  use_tight_wp=true)`); 24 PNGs under
  `~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/mc_based/{step1_singles_data_mc,step2_dr_binned_singles,step3_dr_correction}/`.
  Data sign mapping verified in code: **sign1 = μ⁺, sign2 = μ⁻** (RDFBasedHistFillingPP.cxx:168,
  SingleMuEffcyPtTurnOnFitter.cxx:153-154). Reviewer independently re-derived both plateaus
  (0.9621±0.0123 / 0.8629±0.0317) and the MC/data ratio points (pp foot 1.207, plateau 1.110 —
  the expected L1 over-efficiency, flat ⇒ benign for the ΔR ratio); C3: pp data plateau ~0.81 /
  PbPb 0–5% ~0.75 consistent with Run 2 mu4; pp small-ΔR suppression recovery by ΔR≈0.3–0.35
  qualitatively matches the Run 2 ρ_ΔR close-by correction; the overlay union cross-term
  enhancement has no Run 2 analog (UNVERIFIED, neutral). Iter-1 WARNINGs were cosmetic
  (clipped legend/headline on overlay canvases) — fixed, all 24 regenerated, content unchanged.
  INFO for full-stat rerun: extreme-ΔR tail bins (ΔR>4) are 1.5–2.4σ fluctuations; recheck then.

- 2026-07-10 — **Step 9 (Step-1 SF variant) EXECUTED** (reviews in flight). Old step1 plots
  backed up (`*.bak_20260710_noSF`). NTP now propagates `muon_eff_SF_medium/tight` →
  `MuonFullsimExtra::eff_sf_{medium,tight}` (unfilled ≤0 → 1, counted); all 4 productions
  rerun. **SF fill fractions:** unfilled = EXACTLY the pT < 5 GeV muons (100% below, 0%
  above — official SF-map validity boundary; flat in η). Within the Step-1 Tight+fiducial
  selection: pp 22.68% (μ⁺) / 22.71% (μ⁻); overlay 23.18% / 23.25%. Filled SF_tight mean
  0.961 (pp) / 0.958 (overlay). `FillMCTrigEffHists`: SF-weighted Step-1 numerators
  (`h_mc_*_num_sf_<chg>`, w_sf = ev_weight × active-WP SF). New plots
  `step1_singles_data_mc_sf/` (5/sample: data vs MC vs MC×SF, SF-boundary note on canvas;
  zero-width-8-GeV bin skipped via DivideGraphClean).
  **Agreement check (MC×SF)/data at plateau:** pp μ⁺ 1.206 (4.3 GeV, uncorrected by
  construction) / 1.058 (5.5) / 1.075 (10) / 1.079 (20) / 1.081 (40) — the ~0.96 reco/ID SF
  closes about HALF of the ~10–11% no-SF offset; a flat ~6–8% excess REMAINS above 5 GeV.
  Overlay: 1.185 (4.3) / 1.001 (5.5) / 1.014 (7) / 1.048 (10) / 1.096 (20) / 0.805 (40,
  low-stat) — near-agreement at low pT, residual growing with pT. Physics reading: the
  reco/ID WP scale factor corrects the WP-selection efficiency difference, NOT the L1
  trigger simulation excess — the residual is the genuine trigger-simulation over-efficiency
  (skim doc §8d), which is exactly what stays data-driven in the analysis (ε^nc from data;
  MC contributes only the ΔR ratio, where a flat excess cancels).

- 2026-07-13 — **CORRECTION (user query): the Step-9 SF branches are RECO/ID WP scale
  factors, NOT trigger SFs — the MC×SF vs data-trigger-efficiency comparison is INVALID as
  a trigger-agreement test.** Verified in the skim: `TrigRates.h:447-449` =
  `ToolHandle<CP::IMuonEfficiencyScaleFactors>`; `TrigRates_CA.py:242-253` =
  `MuonEfficiencyCorrectionsCfg(WorkingPoint="Medium"/"Tight",
  CalibrationRelease="250418_Preliminary_r24run3")` — the MCP muon reco/ID efficiency SF
  tool. NO `MuonTriggerScaleFactors` tool exists anywhere in the skim. Both the data T&P
  and the MC Step-1 efficiencies CONDITION on an offline reconstructed Tight muon, so the
  reco/ID efficiency cancels out of both conditionals: multiplying the MC trigger
  efficiency by the reco SF (~0.96) injects an unrelated factor, and the apparent
  "half-gap closure" logged in the Step-9 entry is numerical coincidence, NOT a partial
  trigger correction — that interpretation is RETRACTED. Plots/hists kept as a record
  (they are labeled "reco/ID Tight-WP scale factor", so not mislabeled); do not use them
  as a trigger-agreement test. **Official Run-3 trigger SFs for our chains do not exist**
  (KB: `atlas_run3_muon_performance.md` — no Run-3 HI/low-mu muon performance
  recommendations; `atlas_run2_muon_trigger.md` — central trigger SFs are Z/J/ψ T&P for
  the standard pp menus): our menus are `PhysicsP1_pp_lowMu_run3_v1` /
  `PhysicsP1_HI_run3_v1` (verified from files, skim doc R2) and mu4/2mu4 there have no
  central SF product — which is exactly why ε_trig is measured from data T&P in this
  analysis; the Step-1 data/MC ratio is itself the analysis's own effective trigger SF.
  Awaiting user decision: keep the Step-9 outputs as a documented reco-SF systematics
  ingredient, or remove them.

## Results & Observations

### R1. NTP discovery (2026-07-10, Explore agent + orchestrator check)

**Chain:** `PythiaFullSimAnalysis` / `PythiaFullSimOverlayAnalysis`
(`NTupleProcessingCode/PythiaAnalysisClasses.h:37-85`) → `PythiaAlgCoreT` → `DimuonAlgCoreT`,
with `PythiaFullSimExtras` (reco-truth matching, `ProcessEventFullsim`) and, for overlay,
`PythiaFullSimOverlayExtras`. Runners `run_pythia_fullsim*.sh`, `run_pythia_fullsim_overlay*.sh`.
Output naming: `muon_pairs_pythia_fullsim_<label>` + `_no/with_data_resonance_cuts` +
`extra_output_suffix` (`PythiaAlgCoreT.c:488-504`); labels per `FullSimSampleType.h:54-62`
(pp=`pp24`, hijing=`hijing_overlay_pbpb23` — naming fix already in working tree).

**Trigger-field status:**
- Leg structs (`Muon.h:16-17`) ALREADY have `passmu4`/`passmu4noL1` — never filled in fullsim.
- Pair-level trigger fields (`pass2mu4`, `passmu4mu4noL1`, `passSeparated`, …) live in
  `PairDataExtras` (`MuonPairReco.h:5-14`) which is NOT a base of the fullsim pair structs
  (`MuonPairPythia.h:76-109`) — must be mixed in.
- `PythiaFullSimExtras::InitInputExtra` (`.c:5-30`) binds no trigger branches today.
- Reco index: `fill_reco_quantities` sets `cur_muon.ind = truth_ind`
  (`PythiaFullSimExtras.c:193`); the reco index `truth_to_reco[truth_ind]` is discarded —
  must be kept to index `muon_b_HLT_*` (data path indexes by reco ind,
  `DimuonDataAlgCoreT.c:740-746`).

**Data procedure to mirror (`DimuonDataAlgCoreT.c`):** Run-3 branch names — event
`b_HLT_mu4_L1MU3V`; per-muon `muon_b_HLT_mu4_L1MU3V` (mindR 0.02 default = bare name,
`_0_01` variant; `.c:178-200`); pair `dimuon_b_HLT_2mu4_L12MU3V_0_02` (order-insensitive,
indexed by NTUP pair index; `.c:214-228`); mu4_mu4noL1 3-tier incl. per-leg
`_mu{1,2}passLeg{1,2}_dR_0_02` (`.c:230-267`). mu6/mu8 force-disabled (D7, `.c:65-66`).
`passSeparated = dr>0.8` computed, not read.

**Pair-index bridge (orchestrator, `SkimCode .../TrigRates.cxx:1494-1506`):** the skim fills
dimuon branches over all (i<j) reco-muon combinations and stores
`muon_pair_muon1_index`/`muon_pair_muon2_index`. So fullsim can map its two matched reco
indices → skim pair index per event and read `dimuon_b_*` directly. Pairs with an
unreconstructed leg have no trigger info — irrelevant, all measurements condition on
offline reco muons.

**Single-muon tree:** fullsim already supports `output_single_muon_tree`
(`PythiaFullSimExtras.c:215-225`, truth-level fiducial gate pT>4, |η|<2.4; runners
`run_pythia_fullsim{,_overlay}_single_muon.sh`, suffix `_single_muon`).

**Weights:** pair `weight`==`crossx`== per-(slice×beam) `fullsim_weight_factor`
(σ·ε_filt·isospin/N, `PythiaAlgCoreT.c:683`); isospin 4:6:6:9 baked in, overlay pp-beam only.

**Double-fill hazard (sibling doc):** global fullsim pair trees `muon_pair_tree_sign{1,2}`
are filled 2× (`PythiaFullSimExtras.c:315`); kin trees correct. RDF fullsim readers use the
GLOBAL trees (`RDFBasedHistFillingBaseClass.h:82-83`). Efficiency ratios are immune
(numerator & denominator double alike), but do not read absolute yields from the global tree.

**RDF slot-in:** `RDFBasedHistFillingPythiaFullsim.cxx` (pp, input path `.cxx:15-17`) /
`...FullsimOverlay.cxx` (`SetIOPathsHook` `.cxx:19-30`, `extra_suffix` already threaded);
define+filter pattern to mirror at `RDFBasedHistFillingPythiaFullsim.cxx:128-139`.

### R2. Sample + data-reference discovery (2026-07-10, Explore agent)

**MC NTUPs (tree `HeavyIonD3PD`, 10 000 entries each):**
- pp `_July2026`: **23/24 done**; ONLY `pp_pTH8_14` (task 51360141) still running — its local
  file is still the Apr 13 trigger-off version (140 branches, no `.bak`). New files: 260
  branches, Jul 9. → trigger-mode NTP must SKIP files lacking trigger branches (loud log),
  never default-fill false (would bias ε). Rerun when the 24th lands (D3).
- Overlay: all present (224 branches; adds HI `_VTE50` chains; lacks pp's mu10-15 ladder).
- Per-muon matching: `muon_b_HLT_<chain>` `vector<bool>`, same length as `muon_pt`; variants
  `_0_01`, `_V2`, `_V3`; nominal mindR 0.02 = BARE name (mirrors data). Pair-level:
  `dimuon_b_HLT_2mu4_L12MU3V_0_02` + per-leg `_mu{1,2}passLeg{1,2}_dR_0_02`.
- `muon_charge` does NOT exist — `muon_trk_charge` (existing fullsim code already handles
  charge; only trigger branches are new reads).
- Sanity: `b_HLT_mu4_L1MU3V` overlay pTH8_14 78.3%, pp pTH14_24 86.4% ✓.

**Data references to overlay against:**
- TF1 turn-on fits: `~/usatlasdata/dimuon_data/pp_2024/trg_effcy_pT_fitting_to_*/single_mu_effcy_pT_fit.root`
  (also pbpb_2023 analog; several fit-mode subdirs — confirm nominal mode from pipeline
  before use). Key pattern: `f_pt2nd_vs_q_eta2nd[_ctr0_5]_sign{1,2}_2mu4_sepr_py_<lo>_TO_<hi>_divided`.
  pbpb23 has all 6 centrality tokens incl. required `ctr0_5`. NO TH2D fallback objects in
  either file (contrary to mu4-doc-era expectations).
- Tag-and-probe graphs: `histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin.root` (60
  TGraphAsymmErrors) / `..._pbpb_2023_...` (540). **Graphs exist ONLY vs pT in q·η bins** —
  no 1D eta/phi efficiency graphs anywhere ⇒ Step 1's eta/phi data-MC overlay requires
  adding probe-eta/phi num+denom hists to the data P2 filling and rerunning P2 for pp24 +
  pbpb23 (cheap; data ntuples on disk).
- Fine q·η binning: `RDFBasedHistFilling/CommonEffcyConfig.h:15-26`
  `q_eta_proj_ranges_fine_excl_gap`, 10 bins
  {-2.4,-2.0}{-2.0,-1.6}{-1.6,-1.3}{-0.9,-0.5}{-0.5,-0.1}{0.1,0.5}{0.5,1.0}{1.3,1.6}{1.6,2.0}{2.0,2.2}.
- Plot roots: `~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/` (mu4/,
  mu4_mu4noL1/ subdirs; per-year leaf dirs like `pT_fitting/pp24<wp_suffix>`). MC study dirs
  will follow this convention: one directory per step per sample family.

## Remaining Work

- Implementation Plan steps 1–8 (pp first, then overlay).
- Later: re-run downstream when the 2 in-flight pp slices land (D3) and when the full-stat
  productions arrive (all-centrality PbPb24-conditions overlay; lifts D2).

## Latest Stage

**2026-07-13 — Step 9 (Step-1 SF variant) COMPLETE & REVIEWED** (/review-analysis-code
PASS iter 1, /review-plot PASS iter 1; both with full independent numerical reproduction).
Branch `mc-trigger-efficiency` still unmerged, awaiting user. Standing preconditions for
applying ε_ΔR to crossx unchanged: plateau normalization + full-stat MC + `pp_pTH8_14`
slice. The SF study's conclusion (reco/ID SF closes ~half the pp trigger-eff offset; flat
~6-8% L1-simulation excess remains) reinforces keeping ε^nc data-driven.
