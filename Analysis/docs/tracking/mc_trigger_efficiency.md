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

## Autonomy Contract (DONE 2026-07-14 — items 0–4 complete; item 5 (merge) HELD by user decision)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (0) Step-9 SF variant fully removed: `step1_singles_data_mc_sf/` + the
  `*.bak_20260710_noSF` backups deleted, SF code reverted out of
  `plot_mc_trig_eff.cxx` and out of the `FillMCTrigEffHists` numerators (NTP tree
  propagation of `eff_sf_*` KEPT), SF entries struck from this doc.
  (1) ALL MC trig-eff outputs (pp24 + overlay), Steps 1–3, regenerated on top of
  `bda241a` (ΔR<0.05 truth→reco fallback deleted) + `8c917b4` (Pythia truth-index
  guard): NTP ×4 → FillMCTrigEffHists → FitMCSinglesEffcy → step3 → all plots;
  with a written statement of whether/how each fix moved pp and overlay.
  (2) Evidence-based answers to the two pp Step-2 anomalies (q·η∈(−2.4,−2) low-pT
  MC≫data shape mismatch; the ΔR<0.2 low-pT excess in q·η∈(2,2.2), μ⁺) — code bug
  vs genuine L1 over-efficiency, resolved either way.
  (3) Overlay Step-1 q·η∈(−2.4,−2) low-pT MC<data checked on the rerun outputs;
  overlay Step-2 plots rebinned coarser.
  (4) Step-3 plots: last pair-pT bin (40–120) recolored away from the dark red/blue;
  plateau estimated over ΔR ∈ [1,4].
  (5) Branch `mc-trigger-efficiency` reviewed and merged to master; follow-ups continue
  from master.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

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
- **Menu caveat — RESOLVED 2026-07-14 (trigger convener), NOT a systematic.** Overlay MC uses
  HLT menu `Dev_HI_run3_v1` (L1 `MC_HI_run3_v1`, TriggerValidation prescale set, all relevant
  chains PS=1) vs PbPb23 data HLT `PhysicsP1_HI_run3_v1` (L1 `Physics_HI_run3_v1`, physics
  prescale sets) — data side verified 2026-07-10 from the local test AOD run 462240
  (mc_trigger_info_skim.md R2); pp fullsim uses `PhysicsP1_pp_lowMu_run3_v1`.
  **Convener response:** `Dev_HI_run3_v1` is a SUPERSET of `PhysicsP1_HI_run3_v1` — it
  contains every physics chain of the PhysicsP1 menu plus additional development chains.
  The development chains neither affect our chains (mu4 / 2mu4 / mu4_mu4noL1, which are the
  same chains with the same definitions in both menus) nor are they skimmed. **Using the Dev
  menu in MC is therefore fine: no menu systematic on the ΔR correction, and no menu-driven
  penalty on the overlay-vs-data comparison.** (This does NOT touch the separate, real
  MC-vs-data L1 *simulation* over-efficiency of §8d — that is a detector-simulation effect,
  not a menu effect, and it remains the reason ε^nc stays data-driven.)

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
~~(a) the overlay MC simulates the `Dev_HI_run3_v1` HLT menu with validation prescales while
PbPb23 data ran `PhysicsP1_HI_run3_v1` with physics prescales~~ — **(a) WITHDRAWN 2026-07-14:
the trigger convener confirms `Dev_HI_run3_v1` ⊃ `PhysicsP1_HI_run3_v1` (same physics chains,
plus dev chains that we neither use nor skim), so the menu difference is a non-issue** (see
§4 Menu caveat); (b) the overlay TEST sample was produced with a single fixed vertex position
instead of the Gaussian sampling over 4 positions used in data-like production — vertex-z
affects muon trajectories through the trigger geometry. **D1 still stands on (b) alone**
(plus the overlay's ~4× smaller statistics and the D2 centrality restriction), but the
overlay-vs-data comparison is one systematic *less* degraded than originally recorded.
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

- 2026-07-14 — **Trigger-menu caveat RESOLVED by the trigger convener (user relayed): no
  systematic.** The MC/data HLT-menu difference flagged on 2026-07-10 (`Dev_HI_run3_v1` in
  overlay MC vs `PhysicsP1_HI_run3_v1` in PbPb23 data) is not a physics difference:
  `Dev_HI_run3_v1` **includes all physics chains of `PhysicsP1_HI_run3_v1`** and adds
  development chains on top. The extra chains do not affect our chains and are not skimmed.
  ⇒ §4 "Menu caveat" rewritten from "record as a possible systematic" to RESOLVED; D1's
  reason (a) WITHDRAWN (D1 still holds on the fixed-vertex-z reason (b) + overlay statistics
  + D2). **Unchanged:** the MC L1 *simulation* over-efficiency (skim doc §8d) — a
  detector-simulation effect, not a menu effect — so ε^nc stays data-driven and MC still
  contributes only the ΔR ratio.

- 2026-07-14 — **S10 (item 0) DONE: Step-9 SF variant REMOVED.** Plot dirs
  `step1_singles_data_mc_sf/` (both samples) and the `*.bak_20260710_noSF` backups deleted;
  `plot_mc_trig_eff.cxx` SF blocks reverted (364200c) and the SF-weighted Step-1 numerators +
  SF fill report reverted out of `FillMCTrigEffHists.cxx` (RDF half of 123a6ff) — commit
  `6d30f7e`. Step 1 is data vs MC, no SF factor. **KEPT:** the NTP propagation of
  `muon_eff_SF_{medium,tight}` → `MuonFullsimExtra::eff_sf_*` (genuine reco/ID WP SFs; an
  ingredient for the reco-efficiency systematic, simply not a trigger observable). The SF
  fill report still prints from NTP (22.7% of Tight muons unfilled = exactly pT<5 GeV).

- 2026-07-14 — **S11 (item 1) DONE: staleness audit + full rerun of BOTH samples.**
  **Audit:** the `_mc_trig` NTP trees (07-13 22:58–23:20) postdated `bda241a` but PREDATED
  `8c917b4` (07-14 00:10); the Step-3 hists (`*_step3.root`, 07-10 03:00) predated BOTH.
  ⇒ the PbPb (and pp) plots were NOT run with the fixes. Everything regenerated.
  **Isolation (IMPORTANT):** a sibling session is *actively* refactoring `PythiaAlgCoreT.{h,c}`
  + every `run_*.sh` in the shared tree (→ a single `isTestSample` switch). My first rerun
  compiled a half-edited tree and 3 of 4 NTP runs threw (`'isTestSample' is a protected
  member`, then a lookup in the FULL-sample dir). **ROOT exits 0 even when a macro throws —
  the wrapper's rc=0 was a lie; every stage is now log-grepped for errors.** Rerun redone from
  an isolated `git worktree` pinned to the committed branch tip → immune to the live edits.
  **LATENT BUG FOUND on the branch:** `8c917b4` changed the pp-fullsim isospin default to
  pp-beam-only, but its companion runner fix (`py.setIsospinBeams(true)`, restoring the
  4-beam test-sample behaviour) is still UNCOMMITTED in the sibling's tree ⇒ running the
  *committed* pp runners silently drops 3/4 of the pp statistics (119 209 vs 477 135 muons).
  Reran pp with 4 beams restored (user-confirmed: the isospin mix cannot bias a
  detector-response ratio at fixed (pT, q·η) — it only adds statistics). **Merge is held until
  the sibling commits** (user decision, 2026-07-14).
  **`pp_pTH8_14` (D3) IS IN:** the slice landed 07-10 trigger-enabled (260 branches) ⇒ all 24
  pp beam×slice files processed for the first time; **D3 is now LIFTED.**
  **Impact of the two fixes — measured, not asserted:**
  - `bda241a` (ΔR<0.05 fallback deleted): touches `PythiaFullSimExtras.c`, used by BOTH
    samples, but the production default already had the fallback OFF ⇒ **no change to either**.
  - `8c917b4` (truth-index guard): HIJING-specific by construction (the barcode-restart
    criterion cannot fire on pp's monotonic generator block) ⇒ **pp unchanged**; overlay
    singles 95 410 → 95 797 (+0.4%, the leaked HIJING truth muons removed from the
    reco-matched denominator).
  - pp singles 455 440 → **477 135** (+4.8%) — this is the `pp_pTH8_14` slice, NOT the fixes.
  **New headline numbers (full stats, both fixes in):**
  - pp ε_ΔR^2mu4 plateau (ΔR∈[1,4]) = **0.9588 ± 0.0094** (was 0.963 over [1,3]).
  - overlay ε_ΔR^cross plateau (ΔR∈[1,4]) = **0.8776 ± 0.0257** (was 0.864 over [1,3]).
    Still ~12% below 1 ⇒ **plateau normalization is still required** before it dresses the
    PbPb union cross term (unchanged conclusion).
  - pp Step-1 MC/data: 1.21 (4.3 GeV) → **flat ~1.11** above the turn-on, both charges — the
    expected, benign L1 over-efficiency (cancels in the ΔR ratio).
  - ε-evaluation bookkeeping: pp 17.6% gap-q·η 2D fallbacks, 0.048% floor firings; overlay
    18.6% / 0.27%.

- 2026-07-14 — **S12 (items 3b + 4) DONE, commit `d8235e1`.** Step-2 (overlay only): coarse
  rebinning to a strict SUBSET of the native edges (pT → 7 bins {4, 5.1, 6.5, 8, 10.8, 16.2,
  26.8, 60}; q·η → uniform 0.2) — the 0-5% pair sample on the native 41×184-bin axes was an
  unreadable error-bar forest. pp keeps the fine axes (~4× the pairs). **Physics now legible:
  the three ΔR series overlap within errors in essentially every q·η bin ⇒ the §3.2
  factorization assumption holds at the singles level for the overlay too.** Step-3: plateau
  window [1,3] → **[1,4]**; last pair-pT slice (40–120 GeV) `kMagenta+2` → `kMagenta` (the
  darkened shade read as a third dark red/blue against kRed+1 / kBlue+1).

- 2026-07-14 — **S13 + S14 (items 2 + 3a) DONE — q·η anomalies root-caused (delegated
  investigation; merged from `_sub_mctrig_qeta_anomaly.md`, then deleted). VERDICT: NOT a code
  bug.** Full evidence in R3/R4 below. Three headlines:
  1. **My "MC/data = 2.1" was wrong** — it compared MC's pT 4–6 *average* against the data
     *graph's first point* (pT ≈ 4.1 GeV). Correctly pT-matched: data(4–6) = 0.600,
     **MC/data = 1.50**; 1.04 at 6–8; **1.01 at 8–15**. Data stats are excellent
     (18 276 probes) ⇒ no data-side bias.
  2. **The MC/data offset is BARREL-only, not flat-global.** Plateau ratio: barrel
     (|q·η|<1.0) **1.18–1.26**; endcap (|q·η|>1.3) **1.01–1.03**. The doc's "known 10–15%
     per-leg L1 over-efficiency" is really a ~20% *barrel* effect — KB
     `atlas_run3_muon_performance.md` §B/§C: the barrel L1 is simulated with an optimistic
     lower-bound chamber efficiency; real barrel RPC L1 is degraded by gas-distribution leaks.
  3. The **only** genuine anomaly is a **turn-on-SHAPE failure confined to |η| > 2.0,
     pT < 6 GeV**, where the pp24 fullsim L1 is already *saturated* (0.889 at pT 4.0–4.5 vs
     0.450 in data). Not migration (ε vs truth pT identical, ⟨truth/reco⟩ = 1.0010); φ-flat;
     Medium WP and unweighted alike; present in all 24 slices; absent in the overlay (same
     code) ⇒ the discriminant is the **conditions / L1-muon (TGC) configuration** of r16578
     (pp24) vs r17618 (PbPb23), NOT the menu (convener) and NOT the code (5 independent kills,
     incl. the per-muon flag agreeing to <0.4% with the independently-read pair-level branch,
     `2mu4 && !(both legs mu4) = 0` exactly).
  **Item 3a answered:** the overlay's MC < data in that bin is the *same* forward-endcap
  turn-on story with the opposite sign — the two MC productions disagree with each other and
  each with its own data in that bin, while both match their data at the plateau to ~1–2%.
  **KB LIMIT (stated, not papered over):** the KB does not document the Run-3 endcap L1
  NSW/inner-station coincidence or its MC modelling, so the *microscopic* cause is NOT
  asserted. One open question for the trigger group (see Remaining Work).

- 2026-07-14 — **⚠ PHYSICS-PROCEDURE CONSEQUENCE (needs user decision): assumption (i) of §2
  is VIOLATED for the PbPb union weight.** The item-2 Step-2 "ΔR<0.2 excess" is NOT a
  q·η(2.0,2.2) peculiarity — it is a **ΔR-dependent SATURATION of the single-muon efficiency**
  (R4). §2 assumes "the singles terms ε₁, ε₂ carry **no** ΔR dependence — critical for the
  PbPb union formula, whose linear terms cannot be absorbed into the cross term." Measured
  violation at pT 4–6: **ε(ΔR<0.12)/ε(ΔR>1) = 1.20 (pp) / 1.34 (overlay)** inclusive, up to
  **1.9–2.1** in forward q·η bins.
  - **pp / 2mu4 (product form) — SAFE.** ε_ΔR^2mu4 is *defined* relative to ε₁ε₂ using the
    MC's own Step-1 fits, so it absorbs the singles enhancement by construction, and any
    per-leg normalization error cancels exactly in that ratio (§3.3 self-consistency). The pp
    deliverable is unaffected — including by the forward-saturation anomaly above.
  - **PbPb / mu4 (union form) — AT RISK.** In `P = ε₁ + ε₂ − ε₁ε₂·ε_ΔR^cross`, the **linear**
    ε₁, ε₂ are the data-derived ε^nc, measured at ΔR > 0.8 and applied ΔR-independently. The
    true single-leg efficiency is up to ~34% higher (inclusive) at ΔR < 0.12, and that error
    **cannot** be absorbed by ε_ΔR^cross. **Proposed (NOT yet implemented — awaiting user):**
    carry ε_single(ΔR)/ε_single(ΔR>1) — a *ratio*, in which the L1 normalization offset
    largely cancels — as a ΔR-dependent correction on the union's linear terms, or reformulate
    the union weight. Per the tracking-doc rules the Physics Procedure is NOT changed without
    user approval, so §2 is left as-is and this is flagged.
  - **Also refines the Step-3 note:** the overlay small-ΔR cross-term enhancement (1.95) was
    attributed to "one L1 RoI matching both offline muons". That is now shown to be wrong: the
    enhancement is flat out to ΔR ≈ 0.12, which is **4× the 0.02 matching cone**, so a shared
    HLT object could only double-flag ΔR < 0.04. It is genuine L1 hit-sharing / close-by
    recovery, not double-matching (85% both-legs-fire at ΔR<0.2 vs 37% expected for
    independent legs). KB-coherent: Run-2 ρ_ΔR is the *inefficiency of resolving two
    overlapping RoIs at L1* (and official J/ψ T&P rejects ΔR<0.2 for exactly this), while
    Run-3 added RPC close-by di-muon recovery + an inside-out HLT algorithm ⇒ single-leg match
    ENHANCED, 2-resolved-RoI 2mu4 SUPPRESSED — precisely the two signs measured.

- 2026-07-14 — **S12 plots re-reviewed: /review-plot APPROVED at iteration 3** (log
  `.claude/logs/review-plot-20260714-022830-mc-trig-eff-s12-rebin-ratio.md`). Iter-1 FAIL
  caught a real defect in my coarse rebinning: the pT edges were **rounded literals**
  ({4, 5.1, 6.5, …}) that do not coincide with the native 41-bin variable axis, so
  `TH1::Rebin` was printing "bin edge does not match … result can be inconsistent" and
  grouping by bin centre. Fixed with a `SnapToAxis()` that snaps each target to the nearest
  EXACT native edge and **throws** otherwise; reviewer verified zero ROOT warnings, exact
  content conservation, and num ≤ denom in every coarse bin. Iter-2 FAIL caught a regression
  from my own y-range cap (two overlay Step-3 points pushed off-frame **unmarked**) → every
  off-scale point now carries an up-arrow and is listed on the canvas. **Also added (criterion
  R3, which post-dates the last review of these plots): ratio panels on every Step-1
  (MC/data) and Step-2 (/ΔR≥1) canvas** — this is what made the anomaly quantitative rather
  than eyeballed. All 24 PNGs regenerated; both plateaus unchanged by all of it.

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

### R3. MC-vs-data Step-1 structure: barrel offset vs forward turn-on saturation (2026-07-14)

pp24, charge-averaged, MC (singles tree) / DATA (2D cell integrals over the SAME pT window):

| q·η bin | MC 4–6 | DATA 4–6 | ratio | **plateau ratio (8–15)** | turn-on excess = r(4–6)/r(8–15) |
|---|---|---|---|---|---|
| (−2.4,−2.0) | 0.899 | 0.600 | 1.50 | **1.01** | **1.49** |
| (−2.0,−1.6) | 0.916 | 0.842 | 1.09 | **1.03** | 1.06 |
| (−1.6,−1.3) | 0.734 | 0.754 | 0.97 | **1.01** | 0.96 |
| (−0.9,−0.5) | 0.810 | 0.643 | 1.26 | **1.23** | 1.02 |
| (−0.5,−0.1) | 0.850 | 0.684 | 1.24 | **1.22** | 1.02 |
| ( 0.1, 0.5) | 0.773 | 0.632 | 1.22 | **1.21** | 1.01 |
| ( 0.5, 1.0) | 0.599 | 0.529 | 1.13 | **1.18** | 0.96 |
| ( 1.3, 1.6) | 0.475 | 0.465 | 1.02 | **1.01** | 1.01 |
| ( 1.6, 2.0) | 0.659 | 0.609 | 1.08 | **1.02** | 1.06 |
| ( 2.0, 2.2) | 0.507 | 0.366 | 1.38 | **1.02** | **1.35** |

Two cleanly separated effects: (1) a **barrel-only** flat ~1.18–1.26 plateau offset
(KB-explained: optimistic barrel-RPC L1 simulation); (2) a **turn-on-width excess confined to
|q·η| > 2.0** (1.49 / 1.35) — everywhere else the turn-on excess is 0.96–1.06.

Forward turn-on, pp q·η ∈ (−2.4,−2.0): MC is already **at plateau in its first pT bin**
(0.8888 ± 0.0059 at 4.0–4.5) while data rises 0.450 → 0.924; plateaus agree to 1%. KB Run-2 §9
says the *endcap* turn-on should be BROAD (barrel is the steeper one) ⇒ the saturated MC curve
is the unphysical one. The overlay in the same bin turns on normally (0.355 → 0.856) and sits
*below* its data at low pT — opposite sign, same root cause (item 3a).

### R4. ΔR-dependent saturation of the SINGLE-muon efficiency (2026-07-14) — the §2 assumption-(i) violation

pp24, pT 4–6, Tight, per fine q·η bin — ε(ΔR<0.12) vs ε(ΔR>1.0):

| q·η | ε(ΔR<0.12) | ε(ΔR>1) | close/iso |
|---|---|---|---|
| (−2.4,−2.0) | 0.928 ± 0.012 | 0.905 | 1.03 |
| (−1.6,−1.3) | 0.824 ± 0.016 | 0.730 | 1.13 |
| ( 0.5, 1.0) | 0.697 ± 0.015 | 0.602 | 1.16 |
| ( 1.3, 1.6) | 0.769 ± 0.018 | 0.465 | 1.65 |
| ( 1.6, 2.0) | 0.931 ± 0.010 | 0.661 | 1.41 |
| ( 2.0, 2.2) | 0.928 ± 0.016 | 0.488 | **1.90** |
| **ALL** | **0.828 ± 0.004** | **0.689 ± 0.002** | **1.202** |

Overlay (0–5%): ALL ε(ΔR<0.12) = 0.684 ± 0.014 vs ε(ΔR>1) = 0.512 ± 0.006 → **1.337**;
per-bin up to **2.08** in (2.0,2.2).

**SATURATION, not a q·η peculiarity:** ε(ΔR<0.12) lands at ~0.85–0.94 in *every* q·η bin,
nearly independent of the isolated-muon efficiency (0.47–0.93). The apparent "excess" is just
(saturation − ε_isolated), largest where ε_isolated is smallest — which is why it looks huge in
(2.0,2.2) (ε_iso = 0.49, the smallest) and vanishes in (−2.4,−2.0) (ε_iso = 0.90, already
saturated). **R3 and R4 are two faces of the same MC behaviour.**

**Not double-matching.** Fine ΔR scan (pp, q·η 1.6–2.2, pT 4–6): 0.935 (ΔR<0.04) / 0.955
(0.04–0.08) / 0.907 (0.08–0.12) / 0.780 (0.12–0.20) / 0.611 (0.20–1.0) / 0.608 (>1.0). The
per-muon matching cone is ΔR < 0.02, so a shared HLT object could only double-flag ΔR < 0.04 —
yet the enhancement is flat and full out to ΔR ≈ 0.12. Correlation at ΔR<0.2 (forward legs,
N=320): both legs fire **85.0%**, exactly one 12.2% — vs 37% / 48% expected for independent
legs at ε = 0.61. ⇒ the HLT genuinely reconstructs both close muons.

### R5. Bookkeeping correction (2026-07-14)

The NTP "reco-matched muons" counter (477 135 pp / 95 797 overlay) is **not** the single-muon
TREE entry count (**475 108** pp / **95 389** overlay) — the tree applies a further
`pt>3 && |η|<2.6` gate. Earlier Progress-Log entries quoting 455 440 / 95 410 etc. are the
same class of counter; when quoting statistics, state which one. All physics numbers in this
doc reproduce from the trees (independently verified by the /review-plot reviewer).

**Step-1 comparison asymmetry worth quantifying before the MC/data ratio is quoted in the
note:** the MC single-muon tree denominator contains only **reco-matched (real)** muons
(`PythiaFullSimExtras.c:309-311`), while the data T&P denominator contains **all** offline
muons including fakes (~zero trigger efficiency), which dilutes ε_data at low pT. This inflates
MC/data at low pT *everywhere*; it cannot explain the |η|>2.0 localization (which survives a
plateau-normalized comparison), but it is a real asymmetry in the Step-1 validation.

## Remaining Work

**Blocking / needs user decision:**
1. **The §2 assumption-(i) violation (R4)** — decide how the PbPb **union** weight handles the
   ΔR-dependent singles efficiency (proposed: a ΔR-dependent ratio correction
   ε_single(ΔR)/ε_single(ΔR>1) on the linear terms; alternative: reformulate the union weight).
   The Physics Procedure §2 is deliberately NOT edited pending that decision. **pp/2mu4 is
   unaffected and can proceed.**
2. **Merge of branch `mc-trigger-efficiency` → master is HELD** (user, 2026-07-14) until the
   sibling session commits its `isTestSample`/isospin refactor: commit `8c917b4` changed the
   pp isospin default without its companion runner fix, so the committed pp runners silently
   drop 3/4 of the pp statistics (S11).

**Open questions / follow-ups:**
3. **One question for the trigger group** (the only thing standing between us and the
   microscopic cause of R3): which L1-muon **endcap** configuration (inner/NSW coincidence,
   TGC coincidence-window LUTs) do the **r16578** (pp24-conditions) and **r17618**
   (PbPb23-conditions) simulations use, and is the pp low-μ endcap coincidence applied in data
   but not in that MC?
4. **Double-matching cross-check for R4** (would close it outright): propagate the tighter-cone
   per-muon branch `muon_b_HLT_mu4_L1MU3V_0_01` (and/or `_V2`/`_V3`) via an NTP flag with a
   distinct output suffix (provenance rule). If ε(ΔR<0.12) is unchanged, double-matching is
   excluded.
5. **Quantify the Step-1 denominator asymmetry** (MC = real muons only; data = includes fakes)
   before quoting the Step-1 MC/data ratio as a number in the note (R5).
6. **Plateau-normalize ε_ΔR** before it is applied to crossx (unchanged precondition); pp
   plateau 0.9588 ± 0.0094, overlay 0.8776 ± 0.0257.
7. Re-run downstream when the full-stat productions arrive (all-centrality PbPb24-conditions
   overlay; lifts D2). **D3 is now LIFTED** — all 24 pp slices are in.

## Latest Stage

**2026-07-14 — Steps S10–S15 COMPLETE except the merge (S15), which is HELD by user decision.**
Items 0, 1, 2, 3a, 3b, 4 all done and reviewed (/review-plot APPROVED iter 3; the q·η
investigation returned NOT-a-code-bug with a KB-grounded evidence chain). Two things now sit
with the user: (a) **the §2 assumption-(i) violation** — the singles efficiency IS ΔR-dependent
(R4), which puts the **PbPb union** weight at risk while leaving **pp/2mu4 safe**; (b) the
**merge**, held until the sibling session commits its runner/isospin refactor (S11). See
Remaining Work items 1–2. Everything else regenerated and consistent.

---

**Original PLAN (written before work, retained for the record):**

- **S10 (item 0) — remove the Step-9 SF variant.** User: Step 1 is MC vs data, no SF
  factor. Delete `step1_singles_data_mc_sf/` (both samples) and the
  `step1_singles_data_mc.bak_20260710_noSF` backups; revert the SF blocks from
  `plot_mc_trig_eff.cxx` (364200c) and the SF-weighted numerators + SF fill-report from
  `FillMCTrigEffHists.cxx` (the RDF half of 123a6ff); strike the SF entries from this doc.
  **KEEP** the NTP half of 123a6ff (`muon_eff_SF_{medium,tight}` →
  `MuonFullsimExtra::eff_sf_*`): those are genuine reco/ID WP scale factors and are a
  reco-efficiency-systematics ingredient — they are simply not a trigger observable.
- **S11 (item 1) — staleness audit + full rerun.** Established: the `_mc_trig` NTP trees
  (07-13 22:58–23:20) postdate `bda241a` but PREDATE `8c917b4` (07-14 00:10); the Step-3
  hists (`*_step3.root`, 07-10 03:00) predate BOTH. So every stage is rerun for BOTH
  samples: NTP ×4 → `FillMCTrigEffHists` → `FitMCSinglesEffcy` → `FillMCTrigEffHists`
  (do_step3) → `plot_mc_trig_eff`. Expected impact — `bda241a` touches
  `PythiaFullSimExtras.c` truth→reco matching, which BOTH samples use, and the Step-1/2
  denominator is the reco-matched-muon gate ⇒ can move pp too (the production default
  already had the fallback OFF, so the expectation is *no change*; verify, don't assume).
  `8c917b4` is HIJING-specific by construction (the barcode-restart criterion cannot fire
  on pp's monotonic generator block) ⇒ overlay-only; it removes 643 HIJING truth muons
  that had leaked into the Pythia truth list, which contaminated the reco-matched
  denominator ⇒ overlay Step-1/2 efficiencies may rise slightly. Both to be measured,
  not asserted. Working-tree note: the sibling session's uncommitted isospin/AMI edits are
  behaviour-preserving for these two test samples (pp keeps 4 beams, overlay keeps pp-beam
  only) and add an AMI-DSID provenance guard — weights unchanged.
- **S12 (item 3b + 4) — plot changes.** Overlay Step-2: coarser x-binning (statistics too
  thin for the current pT/q·η bins). Step-3: recolor the last pair-pT bin (40–120) to
  gray/kMagenta (currently confusable with the dark red/blue); plateau estimate over
  ΔR ∈ [1,4] (was [1,3]). → /review-plot.
- **S13 (item 2) — pp Step-2 anomalies** (investigate on the RERUN plots): (a) q·η ∈
  (−2.4,−2), pT 4–6 GeV: MC ≫ data with a different shape, both charges; (b) q·η ∈ (2,2.2),
  μ⁺, pT 4–6 GeV: ΔR<0.2 clearly above the other two ΔR slices (present but much smaller
  elsewhere). Code bug vs genuine L1 over-efficiency. → /review-investigation.
- **S14 (item 3a) — overlay Step-1** q·η ∈ (−2.4,−2), pT 4–6 GeV: MC < data (the opposite
  sign to everywhere else). Real or artifact — check on the rerun.
- **S15 (item 5)** — review + merge `mc-trigger-efficiency` → master; continue from master.

Standing preconditions for applying ε_ΔR to crossx unchanged: plateau normalization +
full-stat MC + the `pp_pTH8_14` slice.
