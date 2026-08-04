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

## Autonomy Contract (round 7 — ACTIVE, opened 2026-08-03; re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a
  passing small test, or one pipeline stage is NOT a stopping point.
- Done = all of the following, in this order (item 1 lands BEFORE any other):
  1. **Truth fiducial as DEFAULT**: `truth_pt > 4 && |truth_eta| < 2.4` is part of the MC
     trig-eff sample selection for **every** sample (pp_full, overlay, noovl) and **every**
     step (1–4), and **all** MC trig-eff plots are regenerated on it (both WPs).
  2. **Sanity-check Step-1 plot set** in its OWN subdirectory, with two extra requirements —
     (i) the event has exactly ONE reconstructed (track-bearing) primary vertex, (ii) the muon
     passes `|truth_pt − reco_pt| / truth_pt < THRESHOLD` with THRESHOLD fixed by me and the
     choice justified — overlaid against the **original MC** (NOT data), plus the reported
     **percentage of MC muons passing each requirement** (and both).
  3. **Conditional**: if after (2) the MC still shows a saturated ≈1 plateau (no data-like
     turn-on) in q·η ∈ (−2.4,−2.0) — i.e. the anomaly is REAL and not a bad-muon artefact —
     then Steps **2, 3 and 4** are remade with the single-muon requirement
     **`pT > 7 GeV || q·η > −2`**, i.e. rejecting ONLY muons that are simultaneously
     `pT < 7 GeV` **and** `q·η < −2`. *(User correction 2026-08-03: the earlier "q·η > −2 for
     all pT" form was wrong — the anomaly is confined to the low-pT turn-on, so only the
     low-pT forward corner is removed; high-pT forward muons are kept.)* A pair is kept only
     if BOTH legs satisfy it, as for every other muon-level requirement.
  4. **ΔR-correction error bars explained**: exact formula (file:line), bug / not-a-bug verdict
     with quantitative evidence, and — if not a bug — why the errors exceed the bin-to-bin
     scatter (highest pair-pT bin 40.6–120 GeV).
  5. **Corrected-MC study**: per-muon weight `SF = ε_data/ε_MC(pT, q·η)`; (a) Step-1 replotted
     with corrected-MC vs data **and the corrected MC fitted**; (b) Step 3 and Step 4 replotted
     as **1 PNG per pair-pT bin, 1 subplot per pair-η bin**, original-MC vs corrected-MC overlaid.
  6. **ΔR-correction fits**: the large-ΔR plateau per (pair pT, pair η) cell is written to a
     **ROOT file** by the step that measures it and **read back** by the fit stage in one
     integrated flow — never hardcoded, never read from a .md/.txt. For any **FULL** sample
     (the HIJING overlay is a TEST sample ⇒ exempt) `|plateau − 1| > 0.1` in any cell **throws
     and exits**. Then fit/interpolate the plateau-normalized corrections (Step 4 flat for
     ΔR ≳ 0.3, Step 3 flat for ΔR ≳ 0.5); plots = absolute values (black) + fitted function
     (red), 1 PNG per pair-pT bin, 1 subplot per pair-η bin, one subdirectory per fit
     function / method.
  7. A **dedicated systematic-uncertainty document** exists and points to this doc / the
     plateau ROOT file for the large-ΔR plateau.
  8. **MC closure written as the next TODO** in this doc — NOT executed (user will check the
     above first).
  9. `/wrap-up` run (docs + git + `.claude` review, incl. whether the CLAUDE.md orchestrator
     instructions need updating), and a final summary listing every new plot path and
     subdirectory breakdown plus the items needing human judgement.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Autonomy Contract (round 6 — **DONE 2026-07-28**; Step 4 delivered + both reviews PASS)
**Every Done item met.** §2 reformulated + §3.4 written; `do_step4` in FillMCTrigEffHists +
Step-4 block in plot_mc_trig_eff (new `step4_dr_correction_singles/` dir, step1/2/3 untouched);
overlay + pp_full × {T,M} run + verified (pp plateau 0.993≈1 validates; overlay small-ΔR rise
reproduces R4 1.34); /review-analysis-code PASS iter 1 + /review-plot PASS iter 1; Remaining Work 1
RESOLVED. Committed. Retained below for the record.

## Autonomy Contract (round 6 — superseded header above; original text)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a
  passing small test, or one pipeline stage is NOT a stopping point.
- Done = **Step 4** delivered: the single-leg ΔR correction `ε_ΔR^single(ΔR) =
  ε_single(ΔR)/ε_single(ΔR>1)` measured by leg-level inverse weighting (§3.4), for
  **overlay (PbPb — the deliverable)** and **pp_full (validation only — NOT physically needed
  for pp; make this explicit in the doc)**, at **both WPs** (Tight + Medium). Concretely:
  (1) `do_step4` mode in `FillMCTrigEffHists.cxx` writing `..._step4.root`; (2) a Step-4
  rendering block in `plot_mc_trig_eff.cxx` (ε_single(ΔR) zoom+full + plateau-normalization +
  pair-pT / pair-η breakdown) → a NEW `step4_dr_correction_singles/` plot dir, nothing existing
  overwritten; (3) §2 union reformulation + §3.4 written (DONE); (4) /review-analysis-code PASS +
  /review-plot PASS; (5) chain run for overlay + pp_full × 2 WP, outputs verified (pp plateau ≈ 1
  validates the machinery; overlay ε_single(ΔR) rises at small ΔR consistent with R4's saturation
  1.20/1.34); (6) doc Progress Log + INDEX updated; committed.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).
- Resolved at start (2026-07-28, user): the RW-1 union-weight form is the `ε_ΔR^single·(ε₁+ε₂)`
  reformulation (§2 note); pp derived for validation and flagged not-needed-for-pp.

## Autonomy Contract (round 5 — **DONE 2026-07-22**; all 4 changes delivered + both reviews PASS)
**Every Done item met.** All 4 user changes applied + validated on all 3 samples (pp FULL,
overlay, r17663) × 2 WPs; grid re-skimmed all 3 (13 tasks, the L1 branch); /review-analysis-code
PASS (iter 1) and /review-plot PASS (iter 4); committed (`dea8c26`..`ec8ef7d`). Headlines:
pp ε_ΔR^2mu4 0.9855/0.9875, overlay ε_ΔR^cross 0.8795/0.8697; #1 overlay HIJING excluded
(99763→95389); #3 eff(L1)≥eff(chain), eff(HLT\|L1)≤1; #4 pp pair-η plateau ~0.97–1.00 (validates
the plateau→1 normalization). Retained below for the record.

## Autonomy Contract (round 5 — superseded header above; original text)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = the four user changes below, applied + validated, with ALL trigger-efficiency plots
  regenerated for pp (**FULL sample**), overlay, and r17663 (both WPs); reviews PASS.
  1. **Muon sample REVERTED to TRUTH-SEEDED Pythia-only** (undo round-2 D4/§3.0 reco-seeding).
     Physics (user): the reco-seeded overlay sample admits real HIJING muons, but Pythia and
     HIJING carry different event weights, so — exactly as for reco-eff / detector-response —
     only Pythia-truth muons (matched to reco + passing a WP) may enter the trigger-efficiency
     sample; HIJING is for template fits (hadronic/fake background) ONLY. Restore the
     pre-round-2 (`221f49a^`) unified truth-seeded `ProcessEventFullsim` + its `store_mc_trigger`
     conditionals; DELETE the reco-seeded `ProcessEventFullsimMCTrig` + `IsRealRecoMuon`; KEEP
     one-sided dp/p (`71fcf1c`) + the multi-file LGD-farm trigger check (`31c2724`). pp result
     identical (set-equality, T5); overlay/r17663 drop HIJING muons. Rewrite §3.0/D4. NTP re-run.
  2. **q·η fine bin (−2.4,−2.0) split → (−2.4,−2.2) + (−2.2,−2.0)** (Run-3
     `q_eta_proj_ranges_fine_excl_gap`). Edit `CommonEffcyConfig.h:16` + the 3 hardcoded lists
     (`SingleMuEffcyPtTurnOnFitter.cxx:45-48`, `plot_mc_trig_eff.cxx:395-401`; legend pad 11→12,
     10→11 bins). Regenerate MC fits AND data T&P (pp24+pbpb23, both WPs — RDF Stage-5 refill +
     Stage-6 refit, `SKIP_CONDOR`) + all comparison plots. The fine 2D q·η axis already has an
     edge at −2.2 (aligns exactly). Crossx would drift → OUT OF SCOPE, noted (Remaining Work).
  3. **Step-2 split into L1 / HLT / full-chain** — user decision (2026-07-21): **RE-SKIM** (no
     per-muon L1 branch exists). Add a per-muon `muon_match_L1MU3V` (offline↔L1 muon-RoI ≥ MU3V,
     prescale-free) branch to the skim, analogous to the existing `muon_match_mu4roi`; re-run
     grid skims for pp-full (6 DSIDs 803015–020) + overlay (802776–781) + r17663 (802781);
     re-farm pp-full to LGD; re-NTP. Step-2 plots: **P[muon fires L1 | offline μ ∧ other μ in
     ΔR bin]**, **P[muon fires HLT | fires L1 ∧ …]**, and the current **P[full mu4 chain | …]**,
     each in its OWN subdirectory. No L1 nor HLT prescale in any efficiency.
  4. **Step-3 pair-η dependence**: 2 plots (each subplot = a pair-η bin from the crossx
     `pair_eta_proj_ranges_coarse_incl_gap`; each line = a coarse pair-pT bin; one zoom-ΔR range,
     one full-ΔR range) + a **large-ΔR-plateau table** over all (pair pT, pair η) bins + a
     **plateau-fluctuation table** (std-dev / error-of-mean, whichever is the right measure).
     Physics: quantify the plateau's deviation from 1 and its pair-pT / pair-η dependence to
     validate the manual plateau→1 normalization and size its systematic. `pair_eta` is an
     existing reco column (`MuonPairReco.h:19`) → no NTP change for #4.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).
- Resolved at start (2026-07-21): #3 clean L1/HLT split needs a re-skim — user chose RE-SKIM.

## Autonomy Contract (round 4 — **DONE 2026-07-20**; retained for the record)
**Every Done item met:** skim submitted (jediTaskID 51600634) and completed; full chain run
(NTP ×2 → Fill/Fit/Step3 × 2 WPs → 24 PNGs) into a fully separate identity — the pp24 and
overlay outputs and their 48 PNGs are byte-untouched (verified); R10 verdict written;
/review-analysis-code and /review-plot both PASS. **Do NOT re-run the chain.**
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = the **r17663 no-overlay sample** (`mc23_5p36TeV.802781...e8599_s4614_r17663_r15970`,
  10 k events) is skimmed on the grid (trigger-enabled, HI menu — the R8 F4 plan), downloaded,
  and pushed through the **full chain to plots**: NTP `store_mc_trigger` → `FillMCTrigEffHists`
  → `FitMCSinglesEffcy` → Step-3 → plot set. **HARD CONSTRAINT (user): separate output file
  names and a separate plot directory — NOTHING existing may be overwritten** (r17663 gets its
  own sample identity: own input dir, own label/suffix, own plot dir; the pp24 and overlay
  outputs and their 48 PNGs stay untouched). Plus: the **R8 outcome-tree verdict answered** —
  ε in bin A (q·η∈(−2.4,−2), pT 4–6, per-muon mu4 match): ≈0.9 pp-like ⇒ NOT the r-tag
  conditions; ≈0.5 overlay-like ⇒ cause IS the r16578 config — written into the doc (R10) with
  the consequence for the trigger-group question; reviews passed; docs + git updated.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Autonomy Contract (round 2 — DONE 2026-07-14; A–D all met, both reviews PASS)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = **(A)** the MC trigger-efficiency muon sample is RECO-SEEDED and TRUTH-MATCHED-REAL
  (`prob>0.5 & |muon_truth_id|==13 & muon_truth_IsPrimary==1`; §3.0 below), replacing the
  truth-seeded Pythia-block loop, for BOTH the single-muon tree and the pair trees (pairs =
  all (i<j) combinations of selected reco muons), in `store_mc_trigger` mode only.
  **(B)** MC offline muons pass EXACTLY the data generic muon cuts (incl. the `fabs(dP/P)`
  change — *historical note: that fabs was itself reversed 2026-07-16, the cut is ONE-SIDED
  on both sides now; see D5 reversal*); every residual data/MC difference enumerated in D5.
  **(C)** BOTH working points produced end-to-end for BOTH samples — Medium-for-both and
  Tight-for-both — for Steps 1, 2 AND 3 (MC *and* the data reference; data Medium T&P via
  `isTight=false` → `_medium_wp`). 4 plot sets total.
  **(D)** All pp (all 24 beam×slice configs) + PbPb overlay test-sample results regenerated
  and verified correct; /review-analysis-code + /review-plot passed; docs + git updated.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Autonomy Contract — round 1 (DONE 2026-07-14; retained for the record)
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

Per-pair trigger probability (inherited from mu4_trig_effcy_implementation.md, with the
single-leg ΔR correction added 2026-07-28 — see the note below):

- **PbPb, mu4 (union — at least one muon fires):**
  `P(pair | ΔR) = ε_ΔR^single(ΔR)·(ε₁ + ε₂) − ε₁·ε₂·ε_ΔR^cross(ΔR)`
  where
  - `ε_ΔR^single(ΔR) = ε_single(ΔR) / ε_single(ΔR>1)` is the ΔR correction on the single-leg
    (marginal) trigger probability — the deliverable of **Step 4 (§3.4)**, plateau-normalized so
    it is 1 for well-separated muons; and
  - `ε_ΔR^cross(ΔR) = P(both fire | ΔR) / (ε₁ ε₂)` is the joint (cross-term) correction — **Step 3
    (§3.3)**.

  By inclusion–exclusion the union needs both the corrected **marginals** (each leg dressed by
  `ε_ΔR^single`) and the corrected **joint** (`ε_ΔR^cross`); the two are measured independently
  (single-leg vs both-legs), so there is no double counting.
- **pp24, 2mu4 (product — both legs fire):**
  `P(pair | ΔR) = ε₁ · ε₂ · ε_ΔR^2mu4(ΔR)`
  — the 2mu4 correction IS the cross-term; no inclusion-exclusion. **The single-leg ΔR
  dependence is fully absorbed into ε_ΔR^2mu4**, so `ε_ΔR^single` is NOT needed for pp (Step 4 is
  run on pp only to validate the procedure on the FULL sample's statistics; see §3.4).

with εᵢ = ε^nc(pTᵢ, q·ηᵢ), the **data-derived** single-muon mu4 efficiency of an isolated muon
(measured at ΔR > 0.8, i.e. already the ΔR>1 / plateau value, so the linear terms reduce to εᵢ
at large ΔR as required). MC supplies only the ΔR **ratios** ε_ΔR^single and ε_ΔR^cross — the
per-leg L1 over-efficiency cancels in each (§1, §4). The two assumptions this doc
measures/tests: (i) ~~the singles terms ε₁, ε₂ carry **no** ΔR dependence~~ **VIOLATED (R4) and
now CORRECTED**: the singles ARE ΔR-dependent, so the union's linear terms are dressed by
`ε_ΔR^single(ΔR)` (Step 4) rather than assumed flat; (ii) the correlation is a function of ΔR
alone (checked via kinematic binning).

**Decision note (2026-07-28, user — resolves Remaining Work 1).** The RW-1 union-weight
question is settled in favour of the "ΔR-dependent ratio correction on the linear terms" option:
`ε_ΔR^single(ΔR) = ε_single(ΔR)/ε_single(ΔR>1)`, measured by inverse weighting exactly as Step 3
measures ε_ΔR^cross but at the single-leg level (§3.4). The union weight is reformulated as
above; pp/2mu4 is unchanged.

### 3. Step-by-step method

#### §3.0 The MC muon sample (AUTHORITATIVE — added 2026-07-14 by user instruction)

Every MC muon entering Steps 1–3 (singles AND pair legs) is defined as follows. This
**replaces** the previous truth-seeded construction (which looped over the Pythia
truth-muon block and kept those with a reco match).

**(a) RECO-SEEDED.** Start from the **reconstructed (offline) muon collection**, exactly as
data does. Never from the truth-muon list — a truth-seeded loop makes the denominator a
truth object with reco decorations, and (in the overlay) silently discards every real muon
that is not in the Pythia signal block.

**(b) TRUTH-MATCHED REAL.** Keep a reco muon only if it is a **real** muon:
`muon_truth_prob > 0.5 && |muon_truth_id| == 13 && muon_truth_IsPrimary == 1`
— the reco-muon provenance classifier's real/hadronic/fake axis, identical to the
template-fit definition (`low_mass_dimuon_template_fit.md` §"HIJING muons ARE real",
`tf_upfront_bkg_reduction.md`). This removes **fake** (`prob ≤ 0.5`) and **hadronic**
(`|id| ≠ 13` punch-through, or `id = ±13 & IsPrimary = 0` decay-in-flight) muons, which
otherwise contaminate the sample **precisely in the regions of interest — low pT and low
ΔR** — and which have (near-)zero trigger probability, biasing ε and the ΔR correlation low.
- **The axis is `(prob, |id|, IsPrimary)` ONLY — NEVER Pythia-signal-block membership**
  (`muon_truth_index`/barcode). ⇒ in the HIJING overlay, **real HIJING muons ARE kept**:
  HIJING simulates hard scattering too, so a truth-matched primary muon from the HIJING
  underlying event is a real muon, and the trigger fires on it exactly as on a Pythia one.
  (The "Pythia-only truth" restriction used for reco-eff / detector-response exists *only*
  because the Pythia AMI slice weights do not apply to HIJING — irrelevant here, because a
  trigger efficiency is a per-event RATIO in which the event weight cancels.)

**(c) EXACTLY THE DATA GENERIC MUON CUTS** (`DimuonDataAlgCoreT::PassCuts_DataCore`) — the
*generic* muon selection, NOT the signal selection:
`quality&1` (combined) · WP bit (`&16` Tight / `&8` Medium) · `quality&32` (IDCuts) ·
`quality&256` (MuonCuts) · `|η| < 2.4` · `pT > 4 GeV` · `Δp/p < 0.12` (**ONE-SIDED** — the
negative tail is kept; user definition 2026-07-16, corrected from the historical `fabs`, see
D5 reversal) · `|d0| < 2 mm` ·
`|z0 sinθ| < 2 mm` · (track-charge agreement only if `turn_on_track_charge`, which is
`false` on both sides).
**No signal-selection cut is applied** (no pair-pT, no q·η, no ΔR, no m_μμ, no resonance
veto) — those define the measurement, not the muon.

**(c') TRUTH FIDUCIAL (added 2026-08-03 by user instruction — DEFAULT for every sample and
every step).** On top of (c), each MC muon must also satisfy
`truth pT > 4 GeV` and `|truth η| < 2.4`.
Because the sample is truth-SEEDED (round-5 change #1), every selected muon has a truth partner,
so this is a well-defined cut on the muon itself; for pairs it is required of **both legs**. Its
purpose is to remove muons that enter the reco fiducial region *only* through mismeasurement —
truth momentum below threshold, or truth direction outside acceptance — which is the same "bad
muon" population §3.5 probes. Implemented at the RDF stage (`FillMCTrigEffHists.cxx`), since
`truth_pt`/`truth_eta` are already stored on both trees.
- **It has NO data analogue** — the data tag-and-probe denominator cannot be truth-gated. The
  Step-1 MC/data comparison therefore carries one more MC-only selection than data does; this is
  a second, deliberate instance of the asymmetry noted at the end of this section. It does **not**
  affect the deliverables, which are MC-internal ΔR ratios.
- Measured impact (pp24 FULL, Tight): keeps 98.49% of the weighted reco-selected sample,
  ε(mu4) 0.7682 → 0.7722 (R13).

**(d) BOTH WORKING POINTS.** Steps 1–3 are produced for **Medium-for-both** and
**Tight-for-both** (MC *and* the data reference at the same WP). The WP must never be mixed
across the data/MC comparison, and the ΔR-correlation procedure (Steps 2–3) must be shown
valid at both.

**Known, accepted asymmetry (not a bug):** the data denominator cannot be truth-matched, so
it still contains fakes/hadronic muons (~zero trigger efficiency), which dilute ε_data at low
pT. The MC/data ratio of Step 1 is therefore NOT a pure trigger comparison; the contamination
fractions are reported so it can be interpreted. This does not affect the deliverable, which
is the MC-internal ΔR **ratio** (Steps 2–3).

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

#### §3.4 Single-leg ΔR correction via inverse weighting (Step 4 — the PbPb union linear-term deliverable)

**Added 2026-07-28 (user), resolving Remaining Work 1 / the §2 assumption-(i) violation (R4).**
Step 2 (§3.2) / R4 established that the single-muon efficiency itself depends on ΔR (saturation
at small ΔR: ε(ΔR<0.12)/ε(ΔR>1) ≈ 1.20 pp / 1.34 overlay, up to ~1.9–2.1 forward). The **pp
2mu4 product absorbs this into ε_ΔR^2mu4, but the PbPb union linear terms do not** — they need
their own ΔR correction. Step 4 measures it as a continuous, kinematics-divided-out curve.

- **What it measures:** `ε_single(ΔR)` = P(a muon fires the full mu4 chain | it is a leg of an
  offline reco pair at separation ΔR), with the leg's (pT, q·η) kinematic dependence divided
  out; plateau-normalized → `ε_ΔR^single(ΔR) = ε_single(ΔR)/ε_single(ΔR>1)`. This dresses the
  union's linear terms (§2).
- **Method (inverse weighting — the leg-level analog of §3.3):**
  - **Object:** each muon **leg** of every MC reco pair (role-swap: both legs of a pair are
    probes; SS + OS trees summed — the trigger response is a per-muon detector property, blind to
    the pair charge product), binned by the **pair ΔR** (fine bins, as Step 3).
  - **Denominator:** all selected legs, unit weight × MC event weight. **NO trigger requirement
    on the leg or on its partner** — the partner only defines ΔR (§4; identical to §3.2).
  - **Numerator:** legs whose **own** per-muon mu4 match fires (never the event-level chain,
    never the partner's decision), each weighted `1/ε_MC(pTleg, q·ηleg)` with ε_MC the **MC-derived**
    §3.1 fit, TF1-Eval'd continuously (clamp to fit range, floor, cap — as §3.3).
  - **Ratio vs ΔR = ε_single(ΔR)**; plateau-normalize over ΔR∈[1,4] → ε_ΔR^single(ΔR).
- **Why MC ε in the weights / plateau = 1 by construction:** same self-consistency as §3.3. If
  the single efficiency had no ΔR dependence (ε_single(ΔR)=ε_MC(pT,q·η)), then in each ΔR bin
  Σ_num 1/ε_MC ≈ Σ_den ε_MC/ε_MC = N → ratio = 1 at every ΔR. A large-ΔR plateau ≠ 1 is a fit-quality
  offset, normalized away; the physical content is the small-ΔR rise (R4's saturation, continuous).
- **Why inverse weighting and not R4's raw ΔR-binned ε:** small-ΔR pairs are kinematically
  special (boosted, higher pair pT → different (pT, η) mix than isolated muons); weighting by
  1/ε_MC(pT,q·η) divides out that kinematic shift, leaving pure ΔR — the continuous version of
  R4's coarse (pT, q·η) table, and the §3.2 "kinematic binning is essential" concern handled
  exactly.
- **Numerator = the full mu4 chain** (what the union weight uses). No inverse-weighted L1/HLT
  split here (it would need L1-only ε fits, which do not exist); Step 2's L1/HLT split (round-5
  #3) already shows the L1-RoI-merging mechanism qualitatively.
- **Diagnostics (mirror §3.3):** large-ΔR plateau (must flatten; offset from 1 sizes the
  normalization systematic), and the plateau's stability across (pair pT, pair η) cells.
- **Sample scope:** the deliverable is the **overlay** (PbPb, union). **pp is run only to
  VALIDATE the machinery on the FULL-sample statistics** (full HIJING-overlay stats not yet
  available); Step 4 results are physically NOT needed for pp (2mu4 product — see §2). Both WPs.

#### §3.5 Step-1 SANITY CHECK: is the forward MC≫data efficiency a "bad muon/event" artefact?

**Added 2026-08-03 (user).** R3/R8/R10 established that in the forward endcap
(`q·η ∈ (−2.4,−2.0)`, `pT ≲ 6 GeV`) the pp24 fullsim MC single-muon efficiency is far above data
and *saturated* — flat at ~0.9 with no data-like turn-on — and that this tracks the r16578
production **configuration**. Round-5 change #2 split that bin into (−2.4,−2.2) + (−2.2,−2.0) and
the behaviour is present in **both** halves. Before accepting it as a genuine simulation problem,
rule out the mundane explanation: that it is driven by badly reconstructed muons or by pile-up
muons that the d0/z0 cuts did not remove.

- **What it measures:** the Step-1 efficiency ε_MC(pT, q·η) recomputed on sub-samples defined by
  two extra requirements, imposed separately and together:
  1. **One-vertex events** — exactly ONE reconstructed, **track-bearing** primary vertex
     (`n_vtx == 1`, counting only vertices with `vtx_ntrk ≥ 2`). This removes any residual
     pile-up muon. The skim dumps `PrimaryVertices` unfiltered, so every event also carries one
     dummy beamspot vertex with `ntrk = 0`; counting it would make `n_vtx == 1` true for every
     event and the requirement vacuous.
  2. **Truth-reco pT match** — `|truth pT − reco pT| / truth pT < 0.10`. This removes badly
     mismeasured muons while keeping essentially all well-measured ones (threshold justification
     below).
- **Comparison target: the ORIGINAL MC, never data.** Neither requirement has a data analogue
  (data has no truth; and the vertex requirement changes the *event* sample rather than the muon
  selection). Overlaying against data would confound the question being asked.
- **Threshold choice (0.10), from the pp24 FULL sample after the §3.0(c') truth fiducial:**
  the `|ΔpT|/pT^truth` distribution has median 0.0144, 90% 0.0390, 95% 0.0489, 99% 0.0721,
  99.9% 0.1101; by region, |η| < 1.05 → 99% at 0.0537, 1.05–2.0 → 0.0779, |η| > 2.0 → 0.0904, so
  a single global threshold is not unfair to the forward region. Scan of candidate values:

  | thr | keeps (weighted) | ε(kept) | ε(rejected) | rejected |
  |---|---|---|---|---|
  | 0.05 | 95.36% | 0.7729 | 0.7581 | 4.64% |
  | **0.10** | **99.82%** | **0.7723** | **0.7311** | **0.18%** |
  | 0.15 | 99.98% | 0.7722 | 0.7209 | 0.02% |
  | 0.20 | 100.00% | 0.7722 | 0.7474 | 0.004% |

  0.05 cuts 4.6% of the sample while the rejected muons are barely worse than average (0.758 vs
  0.773) — it is eating the resolution core. 0.10 sits at ≈5σ of the core resolution: it keeps
  99.8% of the muons yet isolates a population with visibly degraded trigger efficiency
  (0.731 vs 0.772). Beyond 0.15 the rejected sample is too small to be informative. **0.10 is the
  smallest threshold that is clearly outside the resolution core** — good enough for a sanity
  check, which is all that is required.
- **Interpretation / decision rule.** If the saturated ≈1 plateau in `q·η ∈ (−2.4,−2.0)` SURVIVES
  both requirements, the anomaly is not a bad-muon or pile-up artefact and is treated as a real
  simulation problem ⇒ the affected muons are removed from Steps 2–4 by requiring
  **`pT > 7 GeV || q·η > −2`** (i.e. rejecting only muons that are *simultaneously* low-pT and
  forward-negative; both legs of a pair must satisfy it). If instead the plateau develops a
  data-like turn-on, the anomaly is an artefact of those muons and no selection change is needed.
- **Sample scope.** pp24 (the sample where MC ≫ data) is the decisive test — it is the only one
  of the three that has pile-up. Measured N(track-bearing vertices) per event:
  **pp24 fullsim** 1:5.0%, 2:14.7%, 3:21.8%, 4:22.7%, 5:17.1%, 6:10.3%, 7:5.2%, 8:2.1% (mean ≈ 4);
  **HIJING overlay r17618 and r17663** exactly ONE in 100% of events. So requirement 1 costs pp
  95% of its statistics (13.6 M → ~0.7 M muons, still ample) and is a **no-op** on the overlay —
  automatically satisfied, for the correct physical reason (no pile-up), not for lack of
  information. Only requirement 2 discriminates on the overlay.

### 4. Negative constraints

- **NO trigger requirement on any denominator** (the whole point of the MC sample; D9 lesson).
  Any per-muon numerator uses that muon's own trigger match, never the event-level decision.
- **Do NOT replace the data-derived ε^nc in the analysis with ε_MC.** MC contributes the ΔR
  correlation ratio only; §3.1 is validation.
- **Do NOT weight §3.3 numerators with data-derived ε^nc** (breaks the plateau=1 diagnostic).
- **Step 4 (§3.4):** the per-leg numerator uses that leg's OWN mu4 match, never the event-level
  or partner decision; NO trigger requirement on the partner (it only defines ΔR). Weight with
  MC-derived ε_MC (same reason as §3.3). **Do NOT apply ε_ΔR^single to the pp 2mu4 weight** —
  it would double-count the single-leg ΔR effect already inside ε_ΔR^2mu4.
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

### D4: Reco-seeded + truth-matched-real MC muons; real HIJING muons INCLUDED (2026-07-14, user)
**Physics:** see §3.0. The old truth-seeded loop (a) made the denominator a truth object, and
(b) in the overlay discarded every real muon outside the Pythia signal block. Fakes and
hadronic muons have ~zero trigger probability and cluster at low pT / low ΔR — exactly the
regions the ΔR correction is measured in — so leaving them in biases ε and the correlation low;
conversely, dropping real HIJING muons throws away real, triggerable muons that data contains.
**Decision:** reco-seeded, `prob>0.5 & |id|==13 & IsPrimary==1`, no index gate. Applies to the
single-muon tree AND the pair trees (pairs = all (i<j) combinations of the selected reco muons,
so overlay pairs may be Pythia×Pythia, Pythia×HIJING or HIJING×HIJING — as in data).
**Scope guard:** `store_mc_trigger` mode ONLY. The nominal fullsim path (reco-eff,
detector-response, template-fit MC) keeps its truth-seeded Pythia-only construction, which is
correct for *those* measurements (a reco-efficiency denominator MUST be a truth muon, and the
Pythia AMI slice weights do not apply to HIJING).

### D5: MC generic muon cuts must mirror data EXACTLY — the `fabs(Δp/p)` bug (2026-07-14, user)
**Found:** `PythiaFullSimExtras::PassMuonMediumCuts` had `if (muon.dP_overP > thrsh)` while data
has `if (fabs(...) > thrsh)`. `muon_deltaP_overP` is **signed** (41% of reco muons are negative;
**2.6% sit below −0.12**), so MC was ACCEPTING muons that data REJECTS. Fixed to `fabs`.
**Full cut-by-cut audit vs `PassCuts_DataCore` (the only differences found):**
| cut | data | MC (before) | action |
|---|---|---|---|
| `quality&1,&32,&256`, WP bit, `\|η\|<2.4`, `pT>4`, d0/z0 | ✓ | ✓ identical | — |
| `\|Δp/p\| < 0.12` | `fabs()` | **no `fabs`** | **FIXED** |
| track-charge agreement | `turn_on_track_charge=false` | `false` | consistent (both off) |
| Tight WP | `&1,&16,&32,&256` | `pass_medium && &16` (⊃ `&8`) | equivalent — quality bits are cumulative (`getQuality()` Tight=0 < Medium=1, so a Tight muon sets BOTH `&8` and `&16`) |
| WP application | on the PAIR (`m1.q & m2.q & bit`) | per muon | equivalent for pairs |
**⚠ D5 dP/P PART REVERSED (2026-07-16, user):** the dP/P cut is ONE-SIDED by definition
(`dp/p < thrsh`; the negative tail is KEPT) — for data AND MC. Data's `fabs` dated from the
repo's initial commit (2022-11-02) and was always wrong; MC was right all along, and the
round-2 change below copied data's bug into MC. Fixed both ways in `71fcf1c`:
`DimuonDataAlgCoreT.c:599` and `PythiaFullSimExtras.c:160` are now one-sided. The blast
radius below now applies to the DATA side (all data NTP outputs stale); per user decision
(2026-07-16) only the trig-eff data references (pp24 + pbpb23, both WPs) are rerun now, the
full data cascade is scheduled separately. The rest of D5 (cut-by-cut audit) stands.
**Blast radius of the round-2 `fabs` change (superseded by the above — the MC-side rerun
happened in round 3 as the REVERT):** `PassMuonMediumCuts` is shared by
the whole Pythia-fullsim path, so **reco-efficiency, detector-response and template-fit MC
results are now stale** by ~2.6% of reco muons. Out of scope for this task (trigger efficiency);
must be rerun before those are used. `PowhegFullSimExtras.c:25` carries the SAME bug — left
alone (POWHEG fullsim is obsolete, `project_mc_sample_roles`), noted here so it is not lost.
**NOT applied to MC:** the PbPb *event*-level selection (ZDC/FCal/pileup). That is an event
cut, not a generic muon cut; the overlay is by construction a hadronic PbPb sample, and MC
efficiency derivation does not apply data's event cleaning. Overlay stays restricted to 0–5%
centrality (D2).

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

### Round 6 (2026-07-28) — Step 4: single-leg ΔR correction ε_ΔR^single (§3.4)

9.  [x] **§2 reformulation + §3.4** written (Physics Procedure): union linear terms dressed by
    ε_ΔR^single; pp unchanged; Step-4 method = leg-level inverse weighting.
10. [x] **`do_step4` mode in `FillMCTrigEffHists.cxx`** — done; /review-analysis-code PASS iter 1.
11. [x] **Step-4 rendering block in `plot_mc_trig_eff.cxx`** → new `step4_dr_correction_singles/`
    dir; /review-plot PASS iter 1.
12. [x] **Run** step4 + plot for overlay + pp_full × {Tight, Medium} — done. pp plateau 0.993≈1
    (machinery valid); overlay small-ΔR rise reproduces R4 (1.20/1.34). Both reviews PASS. **Step 4
    NOT needed for pp (validation only).** Commit pending in this step.

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

- 2026-07-14 (round 2) — **T1 (NTP) DONE.** `PythiaFullSimExtras`: new
  `ProcessEventFullsimMCTrig()` (reco-seeded loop, `store_mc_trigger` ONLY — the nominal
  truth-seeded path is untouched, D4 scope guard) + `IsRealRecoMuon()` (the §3.0(b) axis) +
  `FillRecoQuantities()` (the lambda extracted to a member function so the two loops share ONE
  reco-muon definition and cannot drift). Bound `muon_truth_{id,IsPrimary,pt,eta,phi,charge}`.
  **`fabs(Δp/p)` FIXED** (D5). Pairs are now all (i<j) combinations of the selected real reco
  muons; **the truth fiducial gate on pairs (`PassCuts_PythiaCore`: truth pT>4, |η|<2.4) is GONE**
  — it was sculpting the reco-pT turn-on. `PerformTruthPairAnalysisHook` deliberately not called
  in trigger mode (it traces the PYTHIA parent map, undefined for HIJING muons; trigger
  efficiency needs no truth ancestry). New end-of-run provenance report.
  **R6 — Reco-muon provenance (full test samples, ALL reco muons, before fiducial/WP):**
  | sample | reco muons | REAL (kept) | fake | hadronic |
  |---|---|---|---|---|
  | pp24 fullsim | 522 522 | **477 135 (91.3%)** | 9 005 (1.7%) | 36 382 (7.0%) |
  | HIJING overlay (all centralities) | 11 321 (400-ev smoke) | **35.9%** | 33.8% | 30.3% |
  - **pp cross-check (important):** REAL = **477 135** is EXACTLY the old truth-seeded
    reco-matched count ⇒ for pp the old Pythia-barcode match and the new real-muon match select
    the **identical muon sample**. So pp's changes come only from the `fabs` fix and the removal
    of the truth fiducial gate on pairs — NOT from the muon definition.
  - **The overlay is the real story:** in central Pb+Pb the reco muon collection is only **~36%
    real** — **64% fake + hadronic**. Reco-seeding WITHOUT the truth match would have put that
    64% (near-zero trigger probability, concentrated at low pT / low ΔR) straight into the
    denominator and destroyed ε and the ΔR correlation. The old code escaped this only by being
    truth-seeded, at the price of discarding real HIJING muons.
- 2026-07-14 (round 2) — **T2 (Medium-WP DATA reference): pp24 DONE** (`isTight=false` →
  `histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin_medium_wp.root`, 20 MB). **pbpb23
  FAILED (OOM)** — the node has 23 GB / 8 cores and RDF ran `EnableImplicitMT(8)` on the large
  PbPb data alongside the NTP job → killed (rc=9, 448-byte stub). Retrying serially with fewer
  threads. No data re-skim or NTuple reprocessing was needed: the data pair trees are Medium
  supersets carrying a per-pair Tight flag.
- 2026-07-14 (round 2) — **T3:** `plot_mc_trig_eff` now selects the **WP-matched data file**
  (`MakeCfg(sample, use_tight_wp)` → `_medium_wp`). It previously hardcoded the Tight data path,
  so a Medium run would have silently compared Medium MC against **Tight data**.

- 2026-07-14 (round 2) — **T4 DONE: all four result sets regenerated** (2 samples × 2 WPs;
  pp = all 24 beam×slice configs). PbPb Medium data reference succeeded on retry (the OOM was
  RDF `EnableImplicitMT(8)` on a 23 GB / 8-core node next to the NTP job; rerun serially with
  3 threads → 53 MB, rc=0). 48 PNGs (12 per sample per WP).
  **R7 — headline numbers (round 2):**
  | quantity | Tight | Medium |
  |---|---|---|
  | pp ε_ΔR^2mu4 plateau (ΔR∈[1,4]) | **0.9583 ± 0.0095** | **0.9589 ± 0.0091** |
  | overlay ε_ΔR^cross plateau | **0.8426 ± 0.0234** | **0.8407 ± 0.0224** |
  | pp Step-1 MC/data (μ⁺, plateau) | ~1.11 | ~1.11 |
  | overlay Step-1 MC/data (μ⁺) | 1.22 (4.3 GeV) → 1.06–1.12 | 1.24 → 1.06–1.13 |
  - **The ΔR correction is WP-STABLE**: pp and overlay plateaus agree between Medium and Tight
    well within errors ⇒ §3.0(d) satisfied, the correlation procedure is valid at both WPs.
  - **ε_ΔR is essentially UNCHANGED by the whole round-2 rework** (pp 0.9588 → 0.9583). This is
    exactly what §3.3's self-consistency argument predicts: ε_ΔR is defined relative to ε₁ε₂
    using the MC's *own* Step-1 fits, so a per-leg change (muon definition, `fabs` cut) cancels
    in the ratio. The deliverable was robust; the Step-1 *validation* is what needed the fix.
  - Tree sizes: pp singles **475 108** (bit-identical to round 1 ⇒ pp muon sample unchanged, as
    predicted by the 477 135 = old-count cross-check); pp pairs SS 53 964 / OS 192 530 (+26–31%
    — the removed truth fiducial gate); overlay singles 95 389 → **99 763** (+4.6%, the real
    HIJING muons); overlay pairs SS 14 176 / OS 37 168.

- 2026-07-14 (round 2) — **⚠ KEY VALIDATION: the R4 small-ΔR singles enhancement SURVIVES the
  truth-match fix ⇒ it is NOT fake/hadronic contamination.** This was the central risk: fakes
  and hadronic muons have ~zero trigger probability and cluster at low pT AND low ΔR, so they
  could in principle have *manufactured* the ΔR-dependent singles effect that R4 reports (and
  which threatens the PbPb union weight). With the sample now **reco-seeded and truth-matched
  real** — i.e. with 64% of the overlay's reco muons (34% fake + 30% hadronic) REMOVED — the
  Step-2 ΔR<0.2 series still sits **~1.2–1.35× above the ΔR≥1 reference at low pT**, in both
  samples and at BOTH working points. The §2 assumption-(i) violation is therefore **genuine
  close-by L1 trigger correlation**, confirmed on a clean real-muon sample. (Contamination
  would also have pushed ε *down* at small ΔR, i.e. the opposite sign to what is observed.)
  ⇒ Remaining Work item 1 (the PbPb union-weight decision) STANDS and is now better founded.

- 2026-07-14 (round 2) — **T5 DONE: both reviews PASS.**
  **/review-analysis-code — PASS at iteration 2.** Iter-1 caught one genuine **silent trap**:
  the reco-seeded pair loop deliberately skips `PerformTruthPairAnalysisHook`, leaving the
  truth-ancestry branches at `Clear()` defaults — but `pair_origin_analysis_skipped` also
  defaulted to `false`, i.e. asserting *"the origin analysis ran"*, and `parent_group == 0` is
  a LEGITIMATE category (the hook's failure sentinel is −10). A future consumer would have read
  a valid-looking all-zeros ancestry block and believed it. Fixed: the flag is now set `true`
  explicitly (verified 100% of pairs in all four `_mc_trig` pair trees). Also removed three
  now-unreachable `store_mc_trigger` branches from the nominal loop, and `require()`d the six
  new truth branches. **The flag fix moved NO physics** — all four plateaus bit-identical.
  **Reviewer verification (stronger than counts — SET equality, event by event):**
  - **pp: `REAL \ OLD = 0` AND `OLD \ REAL = 0`** — the old Pythia-barcode match and the new
    real-muon match select the **identical 477 135 reco muons**. Provably coextensive in pp: with
    no HIJING block, every truth-matched *primary* muon necessarily carries a Pythia-block
    barcode, and every fake/hadronic muon does not. So the old code was *already* rejecting all
    9 005 fakes and 36 382 hadronic muons in pp.
  - **overlay: `OLD \ REAL = 0` (strict subset) but `REAL \ OLD = 6 822`** — the old truth-seeded
    loop never admitted a fake or hadronic muon, but it was **discarding 6 822 real HIJING
    muons** (+7.1% of the real sample). This is exactly the D4 failure mode, and it can only
    occur where a real primary muon lives outside the signal generator block — i.e. only in the
    overlay. **D4 confirmed by measurement.**
  - D4 scope guard verified airtight: `FillRecoQuantities` is a character-for-character
    extraction of the old lambda; the nominal path is byte-equivalent to HEAD except the
    intentional `fabs` fix; `IsRealRecoMuon` has no index/barcode gate; pairs carry no truth
    fiducial cut (tree counts equal Σ C(n_real,2) exactly); `dr` is reco ΔR; C5/C6 clean.
  **/review-plot — PASS at iteration 3.** Iter-1: 4 cosmetic WARNINGs (Step-3 off-scale note
  colliding with the legend; Step-2 legend fill washing out endcap data; ratio-pad points
  clipped without indication; a stale comment inviting a re-hardcoded Tight data path). Iter-2
  caught a **real defect I had introduced**: `DivideGraphClean` paired numerator/denominator
  **by index**, but `TGraphAsymmErrors::Divide` SKIPS empty-denominator bins, so indices shift
  and the x-guard then silently dropped every later point — **12 computable ratio points lost
  per WP in pp**; and the `yn<=0` skip discarded genuine ratio-0 points. Fixed by matching on
  the x value; ratio-0 points are now emitted and down-arrowed. Reviewer re-verified by
  exhaustive **set equality** against an independent derivation: **3703/3703 (pp) and
  1608/1608 (overlay) ratio points present, 0 missing**, no duplicate-x mis-pairing, every
  off-frame point arrowed in both directions.

- 2026-07-16 (round 3) — **U1 DONE (code): dP/P cut polarity corrected — ONE-SIDED everywhere**
  (commit `71fcf1c`). Archaeology: data's `fabs` present since the initial commit `eac35a3`
  (2022-11-02, `MuonNTupleFirstPass.C:28`) — never one-sided, no recent regression; MC was
  one-sided until round-2 D5 (`221f49a`) copied the data bug in. Fixed:
  `DimuonDataAlgCoreT.c:599` + `PythiaFullSimExtras.c:160` (Powheg was already one-sided; no
  other cut site in RDF/plotting layers — grepped). D5 annotated (dP/P part REVERSED).
  User decision on data blast radius: **rerun trig-eff data references only now** (pp24 +
  pbpb23, both WPs); full data cascade (all years nominal NTP, crossx, template fits)
  scheduled separately — tracked in Remaining Work.
- 2026-07-16 (round 3) — **U5 DONE: Tight/Medium WP wiring verified — NO BUG** (delegated;
  merged from `_sub_wp_verify.md`, then deleted). All code sites WP-keyed correctly
  (`FillMCTrigEffHists.cxx:282-283` wp_col pass_tight/pass_medium through singles, both pair
  legs, and Step-3; `FitMCSinglesEffcy.cxx:114-116`; `plot_mc_trig_eff.cxx:269,467` WP-matched
  data file; NTP `PythiaFullSimExtras.c:245-246` pass_medium=PassMuonMediumCuts,
  pass_tight=pass_medium&&(quality&16); data side `RDFBasedHistFillingData.cxx:108` +
  `.Filter("pair_pass_tight")` PP.cxx:149/PbPb.cxx:320). Files: all 8 Tight/Medium pairs differ
  at byte level; **all 1932 compared TH1s differ, 0 bit-identical**; Medium ⊇ Tight everywhere
  (e.g. pp `h_mc_pt_denom_muplus` 172 862 vs 180 829; pbpb data `h_pt2nd_ctr0_5_sign1_mu4_sepr`
  184 197 vs 208 285). **Why plots look identical: 94–95% (MC) / 88–90% (data) of Medium muons
  are also Tight, and the mu4 response is nearly WP-blind — per-pT-bin |ε_T − ε_M| = 0.001–0.011
  absolute (pp MC 0.7844 vs 0.7766; pbpb data 0.5040 vs 0.4937) — invisible at plot scale; the
  Step-3 ratio cancels even that.** Near-identical plots are the physically expected outcome.
- 2026-07-16 (round 3) — **U3(b) DONE: forward-endcap anomaly investigated (delegated; merged
  from `_sub_mctrig_fwd_l1_invest.md`, then deleted). VERDICT: NOT a code bug — a
  production-configuration issue in the pp24 r16578 L1 endcap trigger simulation (~95%
  confidence).** Full evidence in **R8** below: AMI tag chains fetched and compared (4 candidate
  discriminants flip together between r16578 and r17618); matching-variant, vertex-z and
  side/charge hypotheses all KILLED on the NTUPs; the saturation is already present at the
  **L1-item level** (L1_MU3V ≈ 0.95 where data's full chain = 0.60). r17663 (no-overlay clone of
  r17618) verified in AMI/rucio as the single-variable discriminator; one cheap grid skim —
  proposed to user, not submitted. Trigger-group question drafted (R8).
  *(→ PARTLY SUPERSEDED by R10, 2026-07-20: the "r16578 production configuration" conclusion
  STANDS and is now confirmed by r17663, but the localization to **L1** is WITHDRAWN. The
  L1_MU3V ≈ 0.95 measurement itself stands; the inference "⇒ the deficit is at L1" does not —
  R10 measures L1_MU3V ≈ 0.96 in r17663 alongside a chain match of only 0.560.)*

- 2026-07-16 (round 3) — **U1/U4 rerun DONE: everything regenerated on the one-sided dp/p cut**
  (results in R7b + R9). MC: NTP ×4 → RDF Fill/Fit/Step3 ×(2 samples × 2 WPs) → 4 plot sets,
  all rc=0, logs grepped. Data (user-approved scope): `pipeline_pp_trig_eff.sh` +
  `pipeline_pbpb_trig_eff.sh` (YEARS=23, 3 RDF threads) full reruns incl. condor NTP, T&P,
  turn-on fits, validation; Medium T&P refills (`isTight=false`) for pp24 + pbpb23. Plot
  layout: Medium → `step*/medium/` subdirs (U4), old `*_medium_wp.png` removed, 48 fresh PNGs.
  **/review-analysis-code: code + numerics + C1–C3 passed at iteration 1; documentation
  amendments reviewed in further iterations — final verdict recorded in the review log**
  (iter-1 WARNINGs: my overlay Tight/Medium plateau attribution was swapped in the report —
  corrected in R9, files/plots were always right; §3.0(c) + round-2 Done(B) updated to the
  one-sided definition; iter-2 WARNING: stale round-1 plateaus in Remaining Work item 6 —
  updated to R9). Fit status-1 note: R7b.
  Log: `.claude/logs/review-analysis-code-20260716-180312-dpop-one-sided-wp-subdir.md`.

- 2026-07-21 (round 5) — **START. Four user changes (see round-5 Autonomy Contract).**
  Doc triage + INDEX read; full doc re-read; sibling `pythia_fullsim_pp24_full_sample_skim.md`
  read (pp FULL-sample trig-eff infra: `pp_full` knob, LGD farm, plots in canonical
  `pp_trigger_efficiency/mc_based/`, TEST backup `mc_based_TESTSAMPLE_backup_20260720`).
  Verified NO concurrent session: sibling pp-full trig-eff DONE (Tight Jul 21 01:02 plateau
  0.9882±0.0012; Medium Jul 21 12:20 0.9875±0.0012), no active processes ⇒ I own all files.
  **Investigations (2 subagents):** (a) SkimCode has NO per-muon L1-RoI match branch — only
  `muon_b_HLT_mu4_L1MU3V` (full chain, per-muon) + `muon_match_mu4roi` (RoI-seeded HLT-CB,
  ΔR<0.1, per-muon) + event-level `b_HLT_mu4_L1MU3V_L1TBP`. Clean L1/HLT split ⇒ RE-SKIM
  (add `muon_match_L1MU3V`). **User chose RE-SKIM.** (b) #2 data side = cheap RDF re-fill
  (fine 2D q·η axis already has a −2.2 edge) via Stage-5/6 SKIP_CONDOR; MC fits auto-follow
  `CommonEffcyConfig`; 3 hardcoded q·η lists to edit; crossx would drift (out of scope). Run
  mechanics + all mc_trig NTUP input paths confirmed present.
  **Plan:** (Phase A, long pole) SkimCode `muon_match_L1MU3V` + local test + `/review-analysis-code`
  → grid re-skim all 3 samples → re-farm pp-full → re-NTP. (Phase B, parallel, no grid dep)
  #1 truth-seed revert (NTP), #2 q·η split (config+fitter+plot+data), #4 step-3 pair-η
  (Fill+plot+tables), #3 Fill/plot L1-HLT booking — coded + tested on existing NTUPs where
  possible. (Phase C) full chain ×3 samples ×2 WPs → plots → reviews.

- 2026-07-21 (round 5) — **Phase B (code) DONE + compiles; Phase A skim L1 branch DONE + validated.**
  **Code (all four changes written, ACLiC/NTP compile clean):**
  - #1 (truth-seed revert): `PythiaFullSimExtras.{c,h}` — deleted the reco-seeded
    `ProcessEventFullsimMCTrig` + `IsRealRecoMuon` + the provenance report; `ProcessEventFullsim`
    is again the UNIFIED truth-seeded loop with `store_mc_trigger` conditionals (single-muon
    reco-loose gate, pair `ev_pass_*` + `pass2mu4`); kept one-sided dp/p + the multi-file farm
    check; `Muon.h` `MuonFullsimExtra` gained `bool pass_l1`. Overlay uses the same base loop
    (`CallInitInput<E>` calls each mixin) ⇒ change #1 excludes HIJING muons there.
  - #2 (q·η split (−2.4,−2.0)→(−2.4,−2.2)+(−2.2,−2.0)): `CommonEffcyConfig.h:16`,
    `SingleMuEffcyPtTurnOnFitter.cxx:45-49`, `plot_mc_trig_eff.cxx` kQEtaSuffix/kQEtaRange (11
    bins) + legend pad 11→12. Fitter/MCEffEvaluator auto-follow the config.
  - #3 (Step-2 L1/HLT split): NTP propagates `muon_match_L1MU3V`→`m.pass_l1` (bind-if-present,
    loud warning + throw-on-mixed); `FillMCTrigEffHists` books `numl1`(=pass_l1) and
    `numhlt`(=passmu4&&pass_l1) for Step-1 singles AND Step-2 legs; `plot_mc_trig_eff` Step-2 now
    loops 3 stages into subdirs `step2_dr_binned_singles/{full_chain,L1,HLT}/` (numk/denk pick
    the hist infix; eff(chain)=eff(L1)·eff(HLT|L1); numhlt=chain&&L1 keeps eff≤1).
  - #4 (Step-3 pair-η): `FillMCTrigEffHists` books TH3D `h_mc_dr_{zoom,full}_vs_pt_eta_{num,denom}`
    (dR × coarse pair-pT `pair_pt_coarse_bins` × coarse pair-η `pair_eta_proj_ranges_coarse_incl_gap`);
    `plot_mc_trig_eff` adds `step3_eps_dr_{zoom,full}_pair_eta_pt.png` (subplot=pair-η, line=pair-pT)
    + two tables `step3_plateau_pair_eta_pt.txt` / `step3_plateau_fluctuation_pair_eta_pt.txt`
    (plateau mean±stat-err over ΔR∈[1,4] per (pT,η) cell; fluctuation = stat err + weighted RMS
    scatter). `pair_eta` is an existing reco column → no NTP change for #4.
  **Skim L1 branch (delegated, `_sub_skim_l1branch.md`, then merged here):** added per-muon
  `muon_match_L1MU3V` (vector<bool>) to `SkimCode .../TrigRates.{cxx,h}` mirroring
  `muon_match_mu4roi`: retrieve `LVL1MuonRoIs` (xAOD::MuonRoIContainer), match offline↔RoI ΔR<0.2,
  prescale-free. **KEY FINDING: `MuonRoI::thrValue()` is in GeV not MeV** (header comment wrong);
  MU3V is the lowest L1 muon threshold so EVERY in-time L1 muon RoI passes it ⇒ match any RoI
  (no thr cut). **Validated (300 evt): closure P[L1|chain]=0.979** (≥0.9 ✓), true fractions
  L1 0.763 > mu4roi 0.719 > chain 0.681 (L1 loosest → highest, correct), branch superset +1,
  length-invariant 0/300 mismatches, clean build exit 0. INFO (accepted, matches mu4roi template):
  dφ has no 2π wrap (~3% φ=±π strip; folded into closure). Files: TrigRates.h:280,326;
  TrigRates.cxx:15,98,926,1041,1162-1167,1431-1447.
  **Next:** /review-analysis-code (skim gates the grid); submit grid re-skim (pp-full 6 DSIDs +
  overlay + r17663) with the L1 branch; re-farm pp-full to LGD; re-NTP (truth-seeded + L1);
  full chain ×3 samples ×2 WPs; data-side #2 re-fill (pp24+pbpb23, Stage-5/6 SKIP_CONDOR); /review-plot.

- 2026-07-21 (round 5) — **Data-side #2 re-fill DONE + validated** (delegated; `_sub_data_refill.md`).
  All 4 DATA T&P sets regenerated via the RDF class + fitter entry points directly (NOT the full
  pipeline — the combined `muon_pairs` NTUPs were current, so re-ran only Stage-5 fill + Stage-6 fit,
  reading existing NTUPs; no re-hadd, no condor, no re-skim; 3 threads, no OOM):
  - pp24 Tight/Medium: `RDFBasedHillingPP pp(24); pp.trigger_mode=1; pp.isTight=T/F; pp.Run();`
    → `single_muon_trig_effcy_pT_fitting("" / "_medium_wp")`.
  - pbpb23 Tight/Medium: `RDFBasedHistFillingPbPb pbpb(23); pbpb.isTight=T/F; pbpb.Run();`
    → `single_muon_trig_effcy_pT_fitting_PbPb(23,"" / "_medium_wp")`.
  `isTight` is the WP switch (false → `_medium_wp` filenames; the 2 Medium fit files were NEW).
  **11-bin split confirmed** in every fit file: pp24 44 TF1s each (11 q·η × 2 trig × 2 sign), pbpb23
  264 each (× 6 ctr); new `minus2_40_TO_minus2_20` + `minus2_20_TO_minus2_00` present, old
  `minus2_40_TO_minus2_00` = 0. **Sanity (pp24 μ⁺ Tight 2mu4 data): forward (−2.4,−2.2)=0.785 <
  (−2.2,−2.0)=0.854** — the split separates the lower-eff more-forward slice ⇒ the forward region
  DOES have q·η sub-structure (change #2 motivation confirmed on data). 6 pre-split files backed up
  `.bak_pre_qeta_split_20260721`. **Crossx now stale on the forward-bin trig-eff weight (out of scope, noted).**

- 2026-07-21 (round 5) — **r17663 dry-run: all 4 code changes VALIDATED at runtime (NTP/Fill/Fit/hist);
  one plot ORDERING dependency found (not a bug)** (delegated; `_sub_r17663_chain.md`). r17663 v2 NTUP
  downloaded (`grid_monitor.sh --mode no_overlay 51643610`, 9995 entries, `muon_match_L1MU3V` PRESENT).
  Re-NTP (truth-seeded + L1) rc=0: SF report present, **reco-provenance report GONE** (#1 confirmed),
  no L1-absent warning (`has_l1_match=true`). Fill/Fit rc=0: **#2** 22 fits/WP incl. both new forward
  bins; **#3** numl1/numhlt non-empty, muplus dr≥1: eff(L1)=0.721 ≥ eff(chain)=0.700,
  eff(HLT\|L1)=0.914 (≤1), eff(L1)·eff(HLT\|L1)=0.659≈0.700 (~4% ΔR-window gap, as designed);
  **#4** 4 TH3D populated. ε_dR plateau (ΔR∈[1,4]) 0.871 (T)/0.871 (M) (round-4 was 0.853, consistent).
  **⚠ ORDERING CONSTRAINT (record for post-compaction): `plot_mc_trig_eff("noovl")` throws
  (missing `g_mc_pt_vs_q_eta_*_minus2_40_TO_minus2_20`) because its comparison series
  (`qeta_black_mc_file`=pp_full fit, `cmp_fit_file`=overlay fit) are still round-4 10-bin.** ⇒ the
  noovl plot MUST run LAST, AFTER pp_full + overlay are processed with round-5 (their `single_mu_effcy_pT_fit_mc{,_medium_wp}.root` become 11-bin). "pp"/"pp_full"/"overlay" plots have NO
  comparison series (empty cmp_fit_file) ⇒ self-contained, no cross-dependency. **Processing order:
  overlay + pp_full (any order) → re-run plot_mc_trig_eff("noovl",{T,M}) → /review-plot.** r17663's
  NTP+Fill+Fit outputs (incl. its own 11-bin fit) are DONE and on disk; only its final plot is pending.

- 2026-07-22 (round 5) — **Grid re-skim DONE (13/13); download/farm DONE; re-NTP IN PROGRESS.**
  Overlay: 6 v2 NTUPs downloaded+verified (10000 ent + `muon_match_L1MU3V` each), 132.5 GB reclaimed
  (v1 `.bak_20260722` deleted post-verify). pp-full LGD **farm rebuilt for v2** (75 symlinks,
  15/25/17/12/3/3 parts; `fullsim_pp24_full_to_lgd.sh` edited to v2 task IDs + VER_TAG; stale v1
  symlinks removed first; verified reads v2 with has_L1=1, 150000 ent/part). **Now running:**
  pp-full re-NTP (`ppfull_rentp.sh`: both mc_trig scripts over the farm, ~2-3 h, truth-seeded+L1);
  overlay chain (subagent, NTP slice ~3/6 → then Fill/Fit/plot ×2WP). r17663 chain done except its
  final plot. **Next:** pp-full Fill/Fit/plot ×2WP → noovl plot ×2WP (needs pp_full+overlay 11-bin
  fits) → /review-plot on all 3 samples.

- 2026-07-22 (round 5) — **OVERLAY chain DONE + validated (both WPs).** Re-NTP truth-seeded+L1:
  **#1 CONFIRMED** — single-muon tree **95389** (down from reco-seeded 99763; the ~4374 real HIJING
  muons now EXCLUDED). Full chain rc=0, all 4 changes produced outputs:
  - #2: 22 fits/WP (11 q·η bins), 0 failed.
  - #3: `step2_dr_binned_singles/{full_chain,L1,HLT}/` all populated. **μ⁺ ΔR≥1: eff(L1)=0.788 ≥
    eff(chain)=0.644, eff(HLT|L1)=0.795 (≤1), L1·(HLT|L1)=0.627 ≈ chain 0.644** (~2.6% ΔR-window
    gap) — L1/HLT decomposition sound. Medium consistent.
  - #4: `step3_eps_dr_{zoom,full}_pair_eta_pt.png` + plateau/fluctuation tables (VERY noisy — overlay
    60k ev, 0–5% only → sparse (pT,η) cells; only pp-full will populate them well).
  - **ε_ΔR^cross plateau (ΔR[1,4]): Tight 0.8795 / Medium 0.8697** (shifted UP from reco-seeded
    round-3 0.8429/0.8352 — a real consequence of #1: the pair sample is now Pythia-only). Plateau
    normalization still required before it dresses the PbPb union cross term.

- 2026-07-22 (round 5) — **pp-full re-NTP DONE (truth-seeded+L1), final processing running.**
  Re-NTP over the LGD farm (9.87M events) rc=0 both scripts: SF report 19.65M reco-matched muons,
  no reco-provenance report (truth-seeded), no L1-ABSENT warning (L1 bound); outputs 4.6 GB pairs /
  1.0 GB singles. For pp truth-seed=reco-seed (no HIJING) ⇒ change #1 is a no-op for pp (as
  expected). **Now running (`ppfull_noovl_fullchain.sh`):** pp_full Fill/Fit/plot ×2WP → then
  plot_mc_trig_eff("noovl",×2WP) (pp_full+overlay 11-bin fits now exist, resolving the dry-run
  ordering blocker) → completes ALL 3 samples' plots → then /review-plot.

- 2026-07-22 (round 5) — **ALL 3 SAMPLES' PLOTS DONE (both WPs); noovl three-way plot resolved.**
  pp_full full chain rc=0: Step-1 P(mu4) μ⁺ 0.772 (T)/0.761 (M); 22 fits/WP; L1/HLT + pair-η all
  produced. noovl plot now COMPLETES (pp_full+overlay 11-bin fits exist → the dry-run ordering
  blocker cleared); its three-way q·η panels + L1/HLT subdirs produced, both WPs.
  **Headline deliverables (round 5):**
  | quantity (ΔR∈[1,4]) | pp FULL (T/M) | overlay (T/M) | r17663 (T/M) |
  |---|---|---|---|
  | ε_ΔR^2mu4 / ^cross plateau | 0.9855 / 0.9875 | 0.8795 / 0.8697 | (low-stat, recorded) |
  - **#3 L1/HLT (μ⁺, ΔR≥1):** pp eff(L1)=0.775 ≥ eff(chain)=0.771, eff(HLT|L1)=0.948 (≤1);
    overlay eff(L1)=0.788 ≥ eff(chain)=0.644, eff(HLT|L1)=0.795. Product ≤ eff(chain) by
    construction (numhlt=chain∧L1); gap = 1−P[L1|chain] (~5% pp, ~3% overlay, ΔR-window).
  - **#4 pp-full pair-η PLATEAU table (the well-populated deliverable):** plateau **~0.97–1.00** for
    low-pair-pT cells across all pair-η; deviation from 1 ~1–3%, comparable to the per-cell RMS
    scatter ⇒ **the manual plateau→1 normalization is validated and the systematic sized**; mild
    pair-pT dependence (larger at high pT), mild pair-η dependence. Tables at
    `pp_trigger_efficiency/mc_based/step3_dr_correction/step3_plateau{,_fluctuation}_pair_eta_pt.txt`.
  - **#1:** pp truth-seed=reco-seed (no-op); overlay HIJING excluded (singles 99763→95389).
  - **#2:** 11 q·η bins (both new forward bins) in all fits/panels, all samples.
  **Remaining: /review-plot on all 3 samples → commit → done.**

- 2026-07-28 (round 6) — **Step 4 (§3.4): single-leg ΔR correction ε_ΔR^single — CODE + FILLS DONE;
  /review-analysis-code PASS iter 1.** Resolves Remaining Work 1 (§2 union reformulated to
  `ε_ΔR^single·(ε₁+ε₂) − ε₁ε₂·ε_ΔR^cross`).
  - **Code:** `do_step4` mode in `FillMCTrigEffHists.cxx` (4th param, mutually exclusive with
    do_step3; leg-level inverse weighting — denom = all selected legs [MC weight]; num = leg's OWN
    `lg_passmu4` weighted `weight/ε_MC(pT,q·η)`; role-swap + SS/OS summed; writes a separate
    `..._step4.root` with `h_mc_single_dr_{zoom,full}[_vs_pt_eta]_{num,denom}`). Step-4 rendering
    block in `plot_mc_trig_eff.cxx` → new `step4_dr_correction_singles/` dir (ε_single(ΔR) zoom+full
    with plateau line + pair-pT slices + pair-η panels + plateau/fluctuation tables). Both ACLiC clean.
  - **/review-analysis-code PASS iter 1** (2 INFO, both applied; log
    `review-analysis-code-20260728-150308-mc-trig-eff-step4-singles.md`): exact single-leg analog of
    step3, single ε factor (not ε₁ε₂), no trigger req on denom/partner, MC-derived ε in weights,
    provenance clean; overlay numbers verified <0.1%.
  - **Fills (both WPs; inputs = round-5 §3.1 fits + pair trees, NO re-skim/NTP/fit):**
    | ε_single plateau (ΔR∈[1,4]) | small-ΔR rise ε(<0.05)/(0.05,0.1) |
    |---|---|
    | **pp_full T 0.9930 / M 0.9942** | 1.282 / 1.213 (T) |
    | **overlay T 0.9513 / M 0.9435** | 1.359 / 1.165 (T) |
    **VALIDATION (§3.4): pp FULL-sample plateau = 0.993 ≈ 1** — the inverse-weighting is
    self-consistent (plateau=1 by construction with good fits + high stats); the overlay's 0.951 is
    the test-sample fit-quality offset (plateau-normalized away downstream, like ε_ΔR^cross). Small-ΔR
    rise: pp plateau-normalized ≈ 1.29 first bin, ~1.2 over [0,0.12] = **R4's pp 1.20**; overlay ≈ 1.43
    first bin, ~1.29 over [0,0.12] = **R4's overlay 1.34** ⇒ the continuous inverse-weighted curve
    reproduces R4's coarse saturation. WP-consistent (like ε_ΔR^cross). **pp is validation only — NOT
    applied to pp 2mu4 (§2).** Remaining: plots verify + /review-plot → commit.

- 2026-07-28 (round 6) — **Step 4 PLOTS DONE + /review-plot PASS iter 1 → ROUND 6 COMPLETE.**
  New `step4_dr_correction_singles/` plot dir (Tight root + `medium/` subdir) for overlay
  (`pbpb_trigger_efficiency/mc_based/`) + pp_full (`pp_trigger_efficiency/mc_based/`); 5 PNGs + 2
  tables × 4 dirs. step1/2/3 dirs byte-untouched (verified mtimes). Plot plateaus (1D deliverable):
  overlay **0.9508±0.0162 (T) / 0.9433±0.0153 (M)**; pp_full **0.9930±0.0008 (T) / 0.9942±0.0007 (M)**.
  - **pp_full per-cell plateau table = the validation:** every (pair-pT,pair-η) cell 0.94–1.01 (most
    0.98–1.00, tightest at low pT / central η; a few low-stat forward/high-pT corners 0.94±0.06–0.12,
    all within ~1σ of 1) ⇒ inverse-weighting is **kinematics-independent, plateau=1 everywhere, no
    per-cell bug**. Overlay per-cell table + pair-η panels are noisy (test sample, 36 cells, 0–5% only)
    — honest low-stat scatter with error bars + off-scale arrows, NOT a bug; the robust deliverable is
    the 1D plateau. Full overlay production (RW7) will populate the cells.
  - pp plots annotated red "VALIDATION only -- NOT applied to pp 2mu4"; overlay gray "dresses the
    union linear terms (ε₁+ε₂)". (pp note was shortened after an initial right-edge clip; re-verified.)
  - **/review-analysis-code PASS iter 1** (2 INFO applied) + **/review-plot PASS iter 1** (1 INFO, no
    fix; all plateau/small-ΔR numbers verified MATCH, C1/C2/C3 clean). Logs:
    `review-analysis-code-20260728-150308-mc-trig-eff-step4-singles.md`,
    `review-plot-20260728-152100-mc-trig-step4-singles-plots.md`.
  - **Remaining Work 1 (§2 union-weight decision) RESOLVED:** union reformulated
    `P = ε_ΔR^single(ΔR)·(ε₁+ε₂) − ε₁ε₂·ε_ΔR^cross(ΔR)`; ε_ΔR^single delivered (plateau-normalize like
    ε_ΔR^cross before wiring into crossx — RW6 unchanged). pp/2mu4 untouched. Next: commit.

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

### R7b. Round-3 MC rerun bookkeeping (2026-07-16)

- NTP ×4 rc=0, log-grepped clean. Provenance counters unchanged (pp REAL 477 135 — the dp/p
  cut acts downstream of the provenance report). Tree entry counts unchanged (pp singles
  475 108, overlay 99 763): the trigger-mode tree gate is `reco_match && pt>3 && |η|<2.6` (R5);
  the generic cuts incl. dp/p live in the stored `pass_tight`/`pass_medium` booleans, so the
  one-sided revert changes branch VALUES, not entry counts.
- RDF chain (Fill→Fit→Step3, 2 samples × 2 WPs) rc=0. One fit now returns **status 1**:
  `f_mc_pt_vs_q_eta_muplus_minus2_40_TO_minus2_00` (pp Tight, χ²/ndf 58.4/36) — this is
  exactly the R3/R8 anomalous SATURATED bin (flat ≈0.9, no turn-on shape ⇒ degenerate erf
  parameters; it was already the worst bin in rounds 1–2 at χ²/ndf 56.9). The curve is
  numerically indistinguishable from round-2's status-0 fit (Eval differs <0.001 at
  4.5–55 GeV) ⇒ benign, accepted with this note; the bin's absolute MC turn-on is invalidated
  by R8 anyway.

### R8. Forward-endcap L1 anomaly: production-configuration origin (2026-07-16, round 3, delegated investigation)

**The question (user):** same bin (q·η∈(−2.4,−2), pT 4–6) — pp24 fullsim MC ≫ data with a
saturated ~0.9 flat "turn-on" (data rises 0.45→0.92), while the HIJING overlay (PbPb23, 0–5%)
is BELOW its data. Opposite signs cannot be a physical pp-vs-PbPb difference (Run-2: HI ≈ pp,
KB `atlas_run2_muon_trigger.md`). Code bug or simulation issue?

**AMI tag chains (fetched via pyami, both verified reasonable for their campaigns):**
| item | pp24 test `e8599_s4521_s4483_r16578` | overlay test `e8599_s4614_r17618_r15970` |
|---|---|---|
| sim | AthSimulation 24.0.90, FullG4MT_QS | Athena 23.0.56, FullG4MT_QS, **fixed vtx (−0.6,−0.4,−3.3), zero smear** |
| recon release | Athena **24.0.95** | Athena **24.0.58** |
| conditions tag | **OFLCOND-MC23-SDR-RUN3-09** | **OFLCOND-MC23-SDR-RUN3-05** |
| conditions run | autoConfig, runNumber=**801170** | forced ConditionsRunNumber=**460000** |
| overlay | pp pileup (Py8+Epos minbias) | HIJING PbPb UCC ip 0–5 |
| L1 / HLT menu | `Physics_HI_run3_v1` / `PhysicsP1_pp_lowMu_run3_v1` | `MC_HI_run3_v1` / `Dev_HI_run3_v1` |

Same evgen (e8599, Py8.308), same geometry (ATLAS-R3S-2021-03-02-00). The four candidate
discriminants (release, conditions tag, conditions run number, L1 menu) **all flip together**.

**Hypotheses killed on the NTUPs (NTP conventions mirrored exactly, AMI-weighted, pT 4–6):**
1. *Matching artifact:* `_0_01` and `_V3` identical to the bare 0.02 match to the 4th decimal
   in every bin, both samples (`_V2` is a dead R21-only branch, always false in R25 skims —
   `TrigRates.cxx:1345-1365`).
2. *Vertex-z:* pp bin-A ε flat (0.889–0.921) across vz∈[−60,60] mm; overlay fixed-z (rms 0)
   confirmed — irrelevant, the pp saturation exists at every vz including vz≈−3.3.
3. *Single-side hardware:* μ⁺ side C 0.9003 vs μ⁻ side A 0.8967 — sides identical. The
   organizing variable is **q·η (bending direction)**: forward q·η<0 = 0.90 vs q·η>0 = 0.51
   (2.0–2.2) / **0.23 (2.2–2.4)** — MC hugely exaggerates a mild charge×bending asymmetry that
   data also shows (0.600 vs 0.366, R3).
4. **L1 isolated (⚠ "decisive" HEDGED 2026-07-20 — see R10):** in single-selected-muon pp events, the L1_MU3V item alone is
   ~saturated — L1TBP = 0.968 (loose) / 0.949 (strict) in bin A at pT 4–6, where data's full
   chain is 0.60 ⇒ the over-efficiency is dominantly **L1-level**, HLT nearly fully efficient
   on top. (Overlay L1TBP ≈ 0.98 everywhere from PbPb ambient activity — uninformative there.)
   Overlay per-muon ε in the same bin = 0.505 vs pp 0.899 under identical skim + NTP code.
   **⚠ ADDENDUM 2026-07-20 (R10):** the inference "L1 item saturated ⇒ HLT nearly fully
   efficient on top" is NOT safe. R10 measures, in this same bin, an L1_MU3V item at ~0.96 in
   r17663 *alongside* a per-muon chain match of only 0.560 — a direct counterexample to that
   step. R8's evidence is MC-vs-data and R10's is MC-vs-MC, so this is a tension rather than a
   refutation, but item 4 is no longer "decisive" and the L1 attribution below must be read
   through R10's correction.

**VERDICT: NOT a code bug — a production-configuration difference of the pp24 r16578 tag
relative to both 2024 pp data and the HI-conditions simulation (~95% confidence).**
**⚠ the words "endcap L1 configuration" in the original wording of this verdict OVERSTATE the
localization — see R10:** the effect is endcap-specific and confined to pT ≲ 6 GeV, but our
observables cannot separate the L1 threshold/coincidence from the HLT hypo or the chain
matching. Read this verdict as "the r16578 tag's trigger configuration", step unlocalized.
KB limit stated: the KB does not document the Run-3 endcap L1 inner-station (NSW/EI/Tile)
coincidence or TGC coincidence-window LUT configuration in MC — hence a question, not an
assertion. Analysis impact unchanged: ε^nc stays data-driven, ε_ΔR is self-normalized; only
the Step-1 absolute forward turn-on validation is invalidated (+ caution for forward small-ΔR).

**r17663 contingency (user-suggested; VERIFIED in AMI/rucio, NOT submitted):**
`mc23_5p36TeV.802781.Py8EG_A14_pp_hQCD_DiMu_pTH8_14.merge.AOD.e8599_s4614_r17663_r15970`
(tid50446417, 10 files/10 000 ev/5.0 GB, VALID) = r17618 minus HIJING (PileUp=False), all four
candidate discriminants identical to r17618; identical evgen 802781 ⇒ same AMI weight, none
new needed. Processing = ONE cheap grid skim (clone `grid_sub_r17662_signalonly.sh` →
`_r17663` suffix, mode `ppmcfullsim_hioverlay24`) + the diagnostic battery. **Outcome tree
(bin A, pT 4–6):** ε≈0.9 (pp-like) ⇒ NOT the r-tag conditions — deeper L1-simulation issue
common to quiet MC; ε≈0.5 (overlay-like) ⇒ cause IS the r16578 conditions/config — pins the
trigger-group question, HIJING occupancy exonerated.
**→ REALIZED 2026-07-20: the second branch fired (ε = 0.5672 ± 0.0202). See R10 — with two
refinements to this pre-registered wording: occupancy is exonerated as the DOMINANT cause,
not as a contributor at all (it is the whole of the barrel deficit, and a 2.7σ residual
remains in bin A), and the question is pinned to the r16578 CONFIGURATION but NOT to L1.**

**Draft question for the trigger group:**
> In our mc23_5p36TeV 5.36 TeV productions we see opposite-sign data/MC L1 muon endcap
> behaviour between two reco tags. pp-reference fullsim **r16578** (DSIDs 802758–802781,
> `e8599_s4521_s4483_r16578`; Athena 24.0.95; OFLCOND-MC23-SDR-RUN3-09; autoConfiguration,
> runNumber=801170; L1 `Physics_HI_run3_v1`, HLT `PhysicsP1_pp_lowMu_run3_v1` via
> Campaigns.MC23ppReferenceRun2024): per-muon HLT_mu4_L1MU3V efficiency for offline Tight
> muons with |q·η|∈(2.0,2.4), pT 4–6 GeV is ~0.90 and FLAT (the L1_MU3V item alone ~0.95),
> while 2024 pp 5.36 TeV data shows a turn-on from ~0.45; the excess is confined to |η|>2.0
> and strongly q·η-asymmetric (negative q·η saturated, positive q·η suppressed). The
> HI-conditions counterpart **r17618/r17663** (`e8599_s4614_r17618_r15970`; Athena 24.0.58;
> OFLCOND-MC23-SDR-RUN3-05; ConditionsRunNumber=460000; L1 `MC_HI_run3_v1`, HLT
> `Dev_HI_run3_v1`) shows a normal turn-on in the same bin.
> *(⚠ 2026-07-20: before sending, widen this per R10 / Remaining Work 3 — ask which STEP
> (L1 threshold-coincidence vs HLT hypo vs chain matching) differs, rather than presupposing
> L1; note that an RoI-seeded HLT CB muon is already reconstructed in ~90% of the failing
> cases, and that r17663 shows the discriminant is the r16578 configuration, not occupancy.)*
> (1) Which endcap L1 (TGC Sector Logic) configuration does the r16578 trigger simulation
> apply — the big-wheel coincidence-window LUTs and the inner-station coincidence
> (NSW / EI / Tile) — and is the NSW coincidence ENABLED there?
> (2) Is that the configuration applied online in the 2024 5.36 TeV pp reference run, and if
> not, which CW/coincidence set was?
> (3) Which of {Athena release 24.0.95 vs 24.0.58, conditions tag RUN3-09 vs RUN3-05,
> conditions run number 801170 vs 460000, L1 menu Physics_HI_run3_v1 vs MC_HI_run3_v1} drives
> the TGC/NSW coincidence configuration difference between these tags?

### R10. r17663 no-overlay VERDICT (2026-07-20, round 4) — the R8 outcome tree is RESOLVED

**The discriminator fired: r17663 is OVERLAY-LIKE ⇒ the DOMINANT cause of the forward anomaly
is the r16578 tag's conditions/configuration, NOT HIJING occupancy.** ("Dominant", not "sole":
occupancy is real and large elsewhere — it is the whole of the overlay's *barrel* deficit
(0.859 → 0.670) — and r17663 sits 2.7σ *above* the overlay in the forward bin, a residual
consistent with a sub-dominant occupancy component riding on top of the configuration effect.)

Sample: `mc23_5p36TeV.802781...e8599_s4614_r17663_r15970`, skimmed as jediTaskID 51600634
(9 995 events, 224 branches, trigger present), full chain run into its own identity
(`r17663_no_overlay`). All stages rc=0; the one caught runtime message in the NTP logs
("first bin of the same-sign cut-acceptance histogram has ZERO bin content ... return without
normalizing") is the pre-existing, benign `store_mc_trigger` case — the cut-acceptance path is
bypassed in trigger mode — and appears identically in the pp24 and overlay runs. Provenance confirms it is a QUIET sample: **92.2% of reco muons are
real** (6.9% fake, 0.8% hadronic) — essentially pp24's 91.3%, nothing like the overlay's 36%.

**Three-way MC visualization (2026-07-21, user request):** the r17663 Step-1 q·η-bin panels
(`step1_eff_pt_in_q_eta_bins_mu{plus,minus}.png`, both WPs) were reworked into a pure
MC-vs-MC-vs-MC comparison — **red = r17663 (Pb+Pb23 cond., no overlay, + fit); black = pp24
conditions fullsim MC (FULL sample; replaces the pp24 *data* T&P that was there before);
blue = r17618 (Pb+Pb23 cond., with HIJING overlay), no fit.** The ratio pad is now
**r17663 / pp24-cond. MC** (= the R10 comparison, visualized). The pp24-cond. black uses the
FULL sample: bin-A pT4–6 ε = **0.8954** (test sample 0.8977 — consistent to 0.3%), so the
saturated-black-vs-turning-on-red-and-blue picture is unchanged. (The pp-full **Medium** fit
did not exist — the sibling pp-full session had produced only Tight — so it was produced here
from the sibling's complete/stable pp-full trees; new files, nothing clobbered.) Only these 4
PNGs changed; the pp24 and overlay plot sets are untouched. `/review-plot` PASS (below).

**Per-muon mu4 efficiency (Step-1 hists, charge-summed, MC-weighted):**

| q·η bin | pT | **r17663 (no ovl, HI cond.)** | pp24 fullsim (r16578) | HIJING overlay (r17618) |
|---|---|---|---|---|
| **A (−2.4,−2.0)** | **4–6** | **0.5672 ± 0.0202** | **0.8977 ± 0.0039** | **0.4997 ± 0.0144** |
| A (−2.4,−2.0) | 8–15 | 0.9020 ± 0.0416 | 0.9236 ± 0.0059 | 0.8253 ± 0.0190 |
| B (2.0,2.2) | 4–6 | 0.3831 ± 0.0277 | 0.5114 ± 0.0089 | 0.3662 ± 0.0200 |
| BAR (−0.5,−0.1) | 4–6 | 0.8590 ± 0.0124 | 0.8503 ± 0.0038 | 0.6700 ± 0.0113 |

In the anomalous forward bin r17663 sits **16.1σ from pp24 and 2.7σ from the overlay**
(both errors added in quadrature: 0.3305/0.0206 and 0.0675/0.0248) — it
reproduces the overlay's turn-on, NOT pp24's saturation. Per the R8 outcome tree this pins the
cause to the r16578-vs-r17618/r17663 difference (Athena 24.0.95 vs 24.0.58 / RUN3-09 vs -05 /
conditions run 801170 vs 460000 / L1 `Physics_HI_run3_v1` vs `MC_HI_run3_v1`), all of which
r17663 shares with r17618.

**Supporting observables — and an explicit limit on what they localize.** Two further
per-muon quantities, and the barrel as a control:

| observable (bin A, pT 4–6) | pp24 | r17663 | overlay |
|---|---|---|---|
| RoI-seeded HLT CB muon near the offline muon (`muon_match_mu4roi`, ΔR<0.1, pT>4) | 0.908 | **0.902** | 0.727 |
| mu4 chain per-muon match (`muon_b_HLT_mu4_L1MU3V`, ΔR<0.02) | 0.898 | **0.560** | 0.505 |
| barrel (−0.5,−0.1) chain match | 0.850 | **0.858** | 0.670 |

**Barrel control (solid):** r17663 agrees with pp24 (0.67σ) and is far from the overlay
(11.3σ) ⇒ the overlay's barrel deficit is HIJING **occupancy**, and there is **no conditions
effect in the barrel**. The r16578-vs-r17618/r17663 effect is ENDCAP-specific.

**⚠ WHAT THIS DOES *NOT* SHOW (corrected 2026-07-20 after /review-analysis-code CRITICAL).**
An earlier draft of this section read the `muon_match_mu4roi` row as "HLT muon reconstruction,
independent of L1" and concluded the deficit "is at the L1/RoI step ⇒ an endcap L1 coincidence
configuration difference". **That inference is backwards and is withdrawn.**
`HLT_MuonsCB_RoI` (`SkimCode .../TrigRates.cxx:95`, comment at `:1563`) is the **RoI-SEEDED**
HLT container — HLT RoI reconstruction only runs *inside an L1 muon RoI*. So r17663's 0.902
means an L1 RoI **was** present and an HLT CB muon **was** reconstructed there ~90% of the
time, which is evidence *against* an L1-RoI deficit, not for one. The ~44% that fail the chain
match are therefore **predominantly downstream of a successfully reconstructed RoI-seeded HLT
CB muon** — at least **77%** of them, from the marginals alone (chain-fail 0.440, no-RoI-muon
0.098 ⇒ failures-with-an-RoI-muon ≥ 0.342; the two are marginals with *different* cones,
0.02 vs 0.1, and the joint distribution was not measured, so this is a bound, not an equality).
The candidates are then the chain hypo decision, the chain's L1 threshold/quality requirement,
or the tighter ΔR<0.02 chain match — and our observables **cannot separate L1 from the HLT
hypo**. Two more reasons not to
localize it: (i) event-level `L1_MU3V` TBP is ~0.96 in r17663 as in pp24 (though this sample is
DiMu-filtered, so the *event* item can fire on the OTHER muon and is not a per-muon
discriminator either way); (ii) the split is not clean outside bin A — in bin B (2.0,2.2) the
RoI match is *also* depressed in r17663 (0.543 vs pp24 0.595) while the chain match falls
further (0.378 vs 0.507), so "RoI = occupancy only" does not hold generally.
**The KB has no coverage of Run-3 endcap L1 MC configuration, so no microscopic mechanism is
asserted here** — that is the question for the trigger group, not our answer.
**Slice-mixture is not a confounder:** r17663 is single-slice (pTH8_14) while pp24 is 24
beam×slice configs, but round 1 already established the pp saturation is **present in all 24
slices**, so it is not a pT-hat-mixture artifact; and at fixed (pT, q·η) the L1 response is a
detector-configuration effect.

**Consequence for the trigger-group question (R8):** question (3) is now ANSWERED by
elimination on our side — the discriminant is the r16578 tag's **configuration**, not the
environment (HIJING occupancy) — so the question narrows to (1)+(2), and should be **widened
from "which L1 endcap configuration" to "which step"**: we can say the effect is
endcap-specific, confined to pT ≲ 6 GeV, and present with an RoI-seeded HLT CB muon already
reconstructed at the muon — but we CANNOT say whether it is the L1 threshold/coincidence, the
HLT hypo, or the chain matching. Ask them where it is; do not tell them.

**Analysis impact: UNCHANGED.** ε^nc stays data-driven and ε_ΔR is a self-normalized MC ratio,
so no deliverable moves; what this closes is the *interpretation* of the Step-1 forward
validation. r17663's own ε_ΔR plateau (ΔR∈[1,4]) = **0.8534 ± 0.0355** (Tight) /
0.8535 ± 0.0331 (Medium) — low-statistics, single slice; recorded, not used.
(WP attribution verified per-log and against the canvases + an independent recomputation
from `h_mc_dr_full_num/denom` — a first pass had the two WPs swapped.)

### R9. Round-3 headline numbers (2026-07-16, one-sided dp/p everywhere; ALL results current)

| quantity (ΔR∈[1,4] plateau) | Tight | Medium | round-2 (T / M) |
|---|---|---|---|
| pp ε_ΔR^2mu4 | **0.9569 ± 0.0094** | **0.9576 ± 0.0090** | 0.9583 / 0.9589 |
| overlay ε_ΔR^cross | **0.8429 ± 0.0232** | **0.8352 ± 0.0219** | 0.8426 / 0.8407 |

All shifts ≪ 1σ (pp −0.0014 both WPs; overlay +0.0003 T / −0.0055 M), consistent with the
~2.6% re-admitted negative-tail muons; Tight/Medium remain compatible (WP-independence holds).
Data T&P yields moved UP ×1.025–1.043 (reviewer-verified vs the pre-fix snapshot), the
expected sign and size. Data references now carry the one-sided cut (pp24 + pbpb23 pipelines
+ Medium refills rerun); 48 PNGs regenerated — **Tight in the step dirs, Medium in
`step*/medium/` subdirs** (same filenames).

### R11. Run 2 precedent: how the two reference analyses derive ε_trig (2026-07-28, literature)

Answers the question "data-driven, MC, or MC-corrected-to-data?" for the two Run 2 analyses
ours descends from. **They differ, and neither reweights MC events.**

**(a) HF-muon R_AA (ANA-HION-2019-58-INT1 §4.5) — MC central value × data scale factor.**
- pp (Eqs. 15–16): `ε_trig^pp(pT,η) = ε_mu4^MC(pT,η) · SF_trig^pp(pT,η)`, with
  `SF = ε_mu4^pp,data / ε_mu4^MC`. Central value from **Υ(nS)→μμ Pythia8 simulation**, binned
  finely in (pT, q·η) to capture kinematic structure (Fig. 26); SF from **J/ψ→μμ tag-and-probe
  in pp data**, binned **coarsely** (detector geometry) to damp the statistical fluctuations of
  the data measurement (Figs. 27, 38). This split (fine MC shape, coarse data normalization) is
  the stated *reason* for the hybrid.
- Pb+Pb (Eqs. 17–18): `ε_trig^PbPb = ε_mu4^MC · SF_trig^pp · A_trig(pT,η,ΣE_T^FCal,year)`, with
  the extra term **purely data-driven**, `A_trig = ε_mu4^PbPb,data / ε_mu4^pp,data` (same J/ψ T&P),
  Fermi-fitted to smooth bin-to-bin fluctuations; its 68% CI is the systematic. The note's
  explicit reason for not getting the PbPb SF from MC: *"As trigger simulation is not available
  in data overlay, the trigger efficiency scale factor in Pb+Pb cannot be determined by comparing
  Pb+Pb data and data overlay"* (§4.5.2, l. 519–520).
- Systematics (§5.1.2) are assigned to the **SF** (signal/background fit model, ΔR(HLT,offline)
  0.01→0.02, tag quality Medium→Tight, overlap removal) and to the A_trig fit — not to ε_MC itself.
- Note the inversion vs their reco eff, which is **data-driven** (MCP 13 TeV T&P map applied
  directly, §4.6.1 Fig. 30; PbPb ID eff from PbPb J/ψ T&P, Eq. 20).

**(b) Back-to-back dimuon (ATL-COM-PHYS-2021-1094 `Corrections.tex`; identical to HION-2020-10
§3.2.1) — 100% data-driven, no MC anywhere in ε_trig.**
- ε_mu4 from the **MinBias stream**: fraction of offline muons with an HLT muon (RoI) within
  ΔR<0.1, vs (pT, q·η); Fermi+linear fit for pT>8 GeV, raw data-point interpolation below.
- ε_mu4noL1 and ε_(mu4∩mu4noL1) from the **HardProbes single-muon stream** by tag-and-probe on the
  second muon of a pair whose first muon fired mu4/mu6/mu8.
- Centrality dependence as a multiplicative data-measured `CentDep(pT, centrality)` factor.
- The factorization `ε(2mu4)=ε_mu4(a)·ε_mu4(b)` is validated **in data** (HP-stream ε_mu4 vs
  MinBias ε_mu4, Fig. `can_eff1D_Mu4_RefAll`), not against MC. MC (STARlight+HIJING overlay)
  enters only the **reco** efficiency, plus MCP scale factors.
- (The v_n note HION-2019-11 has no efficiency section at all — the EP flow observable is
  self-normalizing.)

**(c) Closure tests: NEITHER analysis closure-tests its efficiencies.**
- HF-muon note: the word "closure" does not occur. Validation = the data/MC efficiency overlay
  that *defines* the SF (Fig. 27) + re-running the analysis under efficiency variations (§5.1.2,
  Figs. 39–40, a systematic, not a closure).
- Dimuon note: `McClosure.tex` exists but the section is renamed *"MC comparisons between truth
  and reconstructed quantities"* — truth-match probability and reco−truth pT/η/φ residuals only,
  no efficiency applied. Their real check is **data-internal**: ε_mu4 from the MinBias stream
  vs from the HardProbes stream agree, which validates both the HP method and the 2mu4
  factorization `Prob(mu4(b)) = Prob(2mu4(a,b)|mu4(a))`. **That check is inclusive in ΔR and
  therefore blind to exactly the correlation this doc measures.** Plus a ⟨w⟩≈2.3–2.4 sanity plot.
- HION-2020-10 has one genuine closure test, but of the **unfolding** (§3.10.2, Fig. 59: split MC,
  unfold one half with the other's response, converges by iteration 5) — not of the efficiencies.
- Structural reason the trigger side could not be closure-tested in Run 2: their overlay MC has
  **no trigger simulation** (R_AA note §4.5.2). That constraint is exactly what the `_July2026`
  skims removed for us (`mc_trigger_info_skim.md`), so an MC trigger closure is available to this
  analysis and was not available to them.

**Bearing on this doc.** Our design (§1–§2) — data-derived ε^nc, MC used only for the ΔR
**ratios** — is a third pattern, and is the conservative one: it never imports an MC absolute
normalization, so the known per-leg MC L1 over-efficiency (≈1.13, barrel-only ~1.2 per R3)
cancels. The HF-muon precedent shows the alternative (MC absolute value dressed by a coarse
data/MC SF) is an accepted ATLAS practice if we ever want MC's fine (pT,q·η) granularity; the
dimuon precedent shows a fully data-driven ε_trig is also accepted. **Neither analysis reweights
MC events to force MC/data efficiency agreement** — the correction is always a multiplicative
factor on the efficiency map applied to *data* yields.

### R12. ΔR-correction error bars: over-stated by 1.5–3.2× (2026-08-03, round 7) — FIXED

**User question:** "How are the error bars for the ΔR corrections determined? For the highest
pair-pT bin (40.6–120 GeV) they are quite large — understandable from statistics — but the
bin-to-bin variations don't look large enough to be compatible with them. Is this a bug?"

**Answer: YES, it was a bug — in the UNCERTAINTY only. Central values were never affected.**

**(a) What the code did.** `FillMCTrigEffHists.cxx` books the denominator with weight `w` (the MC
event weight) and the numerator with `w/ε` (Step 3: `w/(ε₁ε₂)`), Sumw2 on, so per ΔR bin
`D = Σw, e_D = √(Σw²)` and `N = Σ_fired w/ε, e_N = √(Σ_fired w²/ε²)`. The ratio was then a plain,
option-less `TH1::Divide` (5 sites: `plot_mc_trig_eff.cxx:977, 1057, 1155` Step 3 and `1330, 1397`
Step 4, pre-round-7 numbering), i.e.
`e_R = R·√((e_N/N)² + (e_D/D)²)` — verified against the file contents to 3.6e-16.

**(b) Why it is wrong.** That formula assumes the numerator and denominator are INDEPENDENT. They
are not: the numerator is a re-weighted **subset** of the denominator. Conditioning on the MC
sample (D fixed, only the Bernoulli trigger decisions fluctuate) gives
`Var(N) = Σ a_i² p_i(1−p_i)`, `a_i = w_i/ε_i`, `p_i = ε_i·R`, so the code over-states the error by
exactly `√((1+ε)/(1−ε))`. The old comment ("inverse weights > 1 → Bayes invalid") is only half
right: it correctly rules out the *Bayesian* `TGraphAsymmErrors::Divide`, but not the correlated
propagation.

**(c) Evidence.** The measured `e_code/e_correct` matches `√((1+ε)/(1−ε))` to 1–2% in **12/12**
cases (pp_full + overlay × Step 3/4 × 4 pair-pT bins). χ²/ndf of a constant fit over flat
windows: **mean 0.25 with the code errors, 1.0 with the corrected ones** — the smoking gun.
Worst case pp_full Step 4, pT_pair 50–150: **3.22×** over-stated. Internal consistency:
ε_eff(Step 3) ≈ ε_eff(Step 4)² (0.826² = 0.682 vs 0.701), as required for a pair vs a leg.
Ruled out as causes: MC event-weight spread (N_eff/N_raw = 0.863 over the whole pp_full pair
tree — costs only 14% and inflates the true scatter as much as the quoted error), 1/ε weight
spread (N_eff,num ≈ 0.73 N_eff,den), shared events between ΔR bins (none — each pair fills
exactly one bin), rebinning/smoothing (none).

**(d) Why it is worst in the 40.6–120 GeV bin** — two compounding effects: ε rises with pair pT
along the mu4 turn-on (ε_leg 0.69→0.83, ε_pair 0.48→0.70), so the inflation factor itself grows
1.7→3.2; and N_eff is smallest there, so the (already inflated) error is largest.

**(e) Fix applied (round 7).** The fill books, on the numerator node, `errA = Σ w²/ε²` and
`errB = Σ w²/ε`; the ratio error is then `√(A − R·B)/D`. For **Step 4** the two legs of a pair
land in the SAME ΔR bin and their trigger decisions are correlated (that correlation is what
Step 3 measures), so the exact cross term is booked too — once per pair, at `leg == 1`:
`covP = Σ_{both fired} 2w²/(ε₁ε₂)`, `covQ = Σ_{all pairs} 2w²`, giving
`Var = A − R·B + covP − R²·covQ`. Unweighted limit: A = B = N ⇒ `e_R = √(R(1−R)/D)`, the textbook
binomial, as it must be. Consumed by `SetConditionalRatioErrors()` in `plot_mc_trig_eff.cxx`.
- **Boundary guard.** When every effective entry in a bin fired, the conditional variance
  genuinely vanishes (the k = n binomial artefact) and `var` comes out 0, slightly negative, or a
  catastrophic cancellation of terms orders of magnitude larger (observed in the near-empty
  high-pair-pT cells of the 10 k-event overlay: A = 1.1e-06 vs var = 1e-22, which produced a
  bogus `25.56 ± 0.0000` plateau cell). Detected by comparing `var` with the SCALE of its terms
  and replaced by the "1/n rule" with the effective denominator count `n_eff = (D/e_D)²`.
- **Measured effect** (Tight): pp_full Step 3 ±0.0012 → **±0.0008**, Step 4 ±0.0008 → **±0.0004**;
  overlay Step 3 ±0.0256 → **±0.0183**, Step 4 ±0.0162 → **±0.0085**. Central values unchanged.

**(f) PHYSICS CONSEQUENCE — the plateau offset is now significant.** The normalization systematic
is sized by comparing |plateau − 1| with its stat error. With the old errors pp_full Step 4,
pT_pair 50–150 read 0.9846 ± 0.0209 (0.7σ from 1, "consistent"); with the correct error it is
0.9840 ± 0.0065 (**2.5σ**). Inclusively, pp_full Step 4 is now 0.9911 ± 0.0004 — **22σ from 1**.
The plateau≠1 offset is therefore a REAL effect to be normalized away and assigned a systematic,
not a statistical fluctuation. **This enlarges the plateau-normalization systematic and it must
be re-derived** (RW 6).

### R13. Truth fiducial as the default sample selection (2026-08-03, round 7, user)

`truth pT > 4 GeV && |truth η| < 2.4` is now applied to every MC muon in Steps 1–4, all samples,
on top of the data-like reco cuts (both legs for pairs). Implemented purely at the RDF stage —
`truth_pt`/`truth_eta` were already stored on both the single-muon and pair trees — so **no
NTuple-processing re-run was needed** for it.
- **Impact (pp24 FULL, Tight):** keeps **98.49%** of the weighted reco-selected sample;
  integrated ε(mu4) 0.7682 → **0.7722**. 13 649 705 selected muons remain.
- **Headline shifts** (Tight / Medium), vs the round-6 values:
  | quantity (ΔR∈[1,4] plateau) | round 6 | round 7 (truth fiducial + correct errors) |
  |---|---|---|
  | pp_full ε_ΔR^2mu4 | 0.9855 / 0.9875 | **0.9817 ± 0.0008 / 0.9826 ± 0.0008** |
  | pp_full ε_ΔR^single | 0.993 / 0.994 | **0.9911 ± 0.0004 / 0.9918 ± 0.0003** |
  | overlay ε_ΔR^cross | 0.8795 / 0.8697 | **0.8734 ± 0.0183 / 0.8614 ± 0.0176** |
  | overlay ε_ΔR^single | 0.951 / 0.943 | **0.9488 ± 0.0085 / 0.9399 ± 0.0082** |
  | noovl ε_ΔR^cross / ε_ΔR^single | — | 0.8521 ± 0.0237 / 0.9258 ± 0.0113 (Tight) |
  All shifts are small (≤0.006) — the truth fiducial removes a 1.5% population, it does not
  reshape the correction.
- **DOCUMENTED ASYMMETRY (needs human awareness):** the cut has **no data analogue** — the data
  tag-and-probe denominator cannot be truth-gated. The Step-1 MC/data comparison therefore
  carries one more MC-only selection than before. It does not affect the deliverables (the ΔR
  **ratios**, which are MC-internal), but any quoted Step-1 MC/data ratio must state it (cf. R5).

### R15. ΔR-correction fits + plateau ROOT file (2026-08-04, round 7, delegated)

**Integrated flow, no hardcoded plateaus.** The Step-3/Step-4 blocks of `plot_mc_trig_eff.cxx`
now write `<mc_dir>/dr_correction_plateaus_<label>[_medium_wp].root` (Step 3 RECREATE, Step 4
UPDATE): per step, TH2Ds of plateau value / stat error / weighted RMS / n_bins over
(pair pT × pair η) with **axes cloned from the source TH3D**, a 1-bin inclusive plateau, and a
`prov_stepN` TNamed provenance stamp. `fit_dr_corrections.cxx` reads that file — never a
.txt/.md — checks the axis edges *and* the provenance stamp (all samples share a binning, so
only the stamp can catch a cross-sample mix-up), warns if the plateau file is older than the
histograms, runs the guard, divides each cell by **its own** plateau, and fits. Driver:
`pipelines/run_dr_correction_fits.sh` (knobs `SAMPLES WPS STEPS METHODS SKIP_MEASURE SKIP_FIT
STRICT_GUARD`; exit 3 = artefact failure, 2 = guard failure).

**Guard outcome — see `docs/systematic_uncertainties.md` §1a.** Post-veto, two pp24 cells exceed
`|plateau−1| > 0.1`, both in pT_pair[50,150): Step 3 × η_pair[1.0,1.5) = 0.8573 ± 0.0350 and
Step 4 × η_pair[0.5,1.0) = 0.8955 ± 0.0092. **User decision 2026-08-04: accept and normalize;
carry |plateau−1| there as a systematic.** (Pre-veto there were three; the η_pair[2.0,2.4)
failure was removed by the forward veto.)

**Fit functions tried** (inclusive cell, Tight, χ²/ndf with the round-7 conditional errors):

| | polyu_fixedRp | powerlaw_fixedRp | powerlaw_floatRp | expo | interp |
|---|---|---|---|---|---|
| pp_full step3 | **5.89** | 18.31 | 11.70 | 7.79 | exact |
| pp_full step4 | 25.18 | 47.07 | **15.68** | 15.39 | exact |
| overlay step3 | 1.343 | **1.272** | 1.319 | 1.355 | exact |
| overlay step4 | 1.205 | **1.138** | 1.205 | 1.231 | exact |

Flat-beyond-R_p audit (the §-requirement that Step 4 is flat for ΔR ≳ 0.3 and Step 3 for
ΔR ≳ 0.5): `polyu_fixedRp`, `powerlaw_fixedRp`, `interp` → **0 violations, worst |f−1| = 0**
(exact by construction); `powerlaw_floatRp` → violations with worst |f−1| ≈ 0.10–0.13; `expo`
→ many violations, worst 0.112 (pp) / 16.6 (overlay) — it never returns to 1.

**Recommendation:** nominal `polyu_fixedRp`, `f = 1 + u²(a₂+a₃u+a₄u²)` with
`u = max(0, 1−ΔR/R_p)` — exactly 1 and C¹ beyond R_p, and the only parametric form that
reproduces the measured non-monotonic small-ΔR shape. Systematic variant `powerlaw_fixedRp`.
**Reject `expo` and `powerlaw_floatRp`: they violate the required flatness.**
- ⚠ **Honest caveat:** on the pp FULL sample **no** 2–4-parameter form reaches χ²/ndf ≈ 1 — with
  per-mille errors the data resolve real structure the smooth forms cannot follow (this is the
  same residual ~0.5% non-flatness of the "plateau" noted in R12(f)). `interp` is exact by
  construction and may be the better choice for pp; on the overlay every method gives ≈1.1–1.4.
  **A human should pick between `polyu_fixedRp` and `interp` for pp.**

**Persistence verified.** All functions are **TFormula-string** TF1s (never C++ lambdas): re-read
in a *fresh* ROOT session with no macro loaded, `f(0)=1.0361`, `f(0.5)=f(10)=f(50)=1.0` — nothing
collapses to 0 outside the fit range (the `project_tf1_eval_out_of_range` trap). 0 persistence
failures across all files. Consumers must still clamp at the point of use.

**Physics worth a human's eye:** the small-ΔR correction has **opposite sign** in the two
systems — pp `ε_ΔR^2mu4` **dips to ≈0.81** (2mu4 needs two L1 RoIs and close muons share one)
while overlay `ε_ΔR^cross` **rises to 2.22 ± 0.19** (mu4 cross term). Do not assume one sign.

**⚠ Overlay per-cell plateaus are NOT measurable today:** 33 of 36 overlay Step-3 cells have
|plateau−1| > 0.1, spanning 0.04 to 16.1 — 10 000 events give too few pairs per (pair pT, pair η)
cell in ΔR ∈ [1,4]. **Only the inclusive overlay curve is usable** (plateau 0.8681 ± 0.0189);
per-cell PbPb corrections need the full-statistics overlay (RW 7).

### R16. Corrected-MC study: does correcting MC to data change the ΔR corrections? (2026-08-04)

`SF(pT, q·η) = ε_data/ε_MC` from the fitted TF1s, evaluated continuously; **every MC muon that
FIRED carries the extra weight SF, the denominator is untouched**, so
`ε_corr = Σ_fired w·SF / Σ_all w ≈ ε_data`. Step 3/4 corrected numerator weight `w·SF/ε_corr`
with ε_corr from the **corrected** fit file (self-consistency).

**Q1 — do the two methods give the same single-muon efficiency?** Corrected MC reproduces the
data turn-on to **mean |⟨ε_corr/ε_data⟩ − 1| = 0.51% (pp24 FULL)** and **5.6% (overlay)**. Worst
point: 0.69 at pT = 4 GeV in q·η ∈ (−2.4,−2.2) μ⁻ (pp) — i.e. exactly the anomalous
forward-negative turn-on where ε_MC is saturated and SF is extreme, the region the forward veto
now removes from Steps 2–4.

**Q2 — does correcting the MC change the correction terms?** **No, to the precision that
matters.** Algebraically it should not: if ε_corr were exactly ε_data then `SF/ε_corr = 1/ε_MC`,
the original weight, and the corrections would be identical bin-by-bin. Measured
(pull = (corr−orig)/σ_orig, σ_orig the conditional error; the two series share the same events
so a quadrature error would be meaningless):

| sample | Step 3 max \| mean per-cell median | Step 4 max \| mean per-cell median |
|---|---|---|
| pp_full | 3.20σ \| **0.047σ** | 3.49σ \| **0.049σ** |
| overlay | 28.1σ \| **0.184σ** | 51.6σ \| **0.210σ** |

⇒ typical displacement is a few hundredths of a σ. The isolated large pulls are all in the
noise-dominated overlay cells (the same ones whose plateaus swing 0.04–16.1) where the
conditional error collapses at the k=n boundary, so a tiny difference becomes a huge pull — an
artefact of a vanishing denominator, not physics. **This validates the nominal procedure: the ΔR
correction is insensitive to the overall single-muon efficiency normalization and therefore
transfers from the MC world to the data world.**

**Two real bugs found and fixed inside the new corrected code path** (they never touched the
nominal path): (i) `TGraphAsymmErrors::BayesDivide` is undefined when the corrected numerator
exceeds the denominator (SF > 1) — it returns an EMPTY graph, and `TGraph::Fit` on an empty graph
is a **no-op that still returns status 0**, leaving the TF1 at its initial parameters (a silent,
badly wrong turn-on of 0.4500 at pT = 4). Replaced by an explicit conditional-error estimator, and
the fitter now counts `Npts == 0` / `ndf <= 0` as FAILED. (ii) An `eff = 0 ± 0` point (denominator
> 0, no fired muons) is infinitely constraining and gave χ²/ndf = 30627; the "1/n rule" now covers
**both** binomial boundaries.

**Validation:** an `sf_closure` mode forces SF ≡ 1, under which the corrected chain must reproduce
the nominal one exactly — verified bin-by-bin, worst relative difference **4.6e-16** (Step 3) /
**3.9e-15** (Step 4). That closure also caught two live concurrency hazards during the work (a
half-written input NTP file, and a mid-task selection change), each of which would otherwise have
been reported as physics.

## Remaining Work

**Blocking / needs user decision:**
1. ~~The §2 assumption-(i) violation (R4)~~ **RESOLVED 2026-07-28 (round 6, user).** The PbPb union
   weight now dresses its linear terms with a ΔR-dependent ratio correction
   `ε_ΔR^single(ΔR)=ε_single(ΔR)/ε_single(ΔR>1)`, measured by leg-level inverse weighting
   (**Step 4 / §3.4**). §2 reformulated to `P = ε_ΔR^single·(ε₁+ε₂) − ε₁ε₂·ε_ΔR^cross`. Delivered for
   overlay (T 0.951 / M 0.943 plateau, small-ΔR rise → R4's 1.34); pp_full ran as validation (plateau
   0.993 ≈ 1, per-cell 0.94–1.01 ⇒ machinery self-consistent) and is **NOT applied to pp** (2mu4
   product absorbs it). Like ε_ΔR^cross, **plateau-normalize ε_ΔR^single before wiring into the PbPb
   crossx union weight (RW 6)** — that application step is the remaining downstream work.
2. ~~Merge held~~ — **DONE 2026-07-14, merge commit `8edc4fb`** (46 commits, no conflicts).
   The gate cleared when the sibling session committed its `isTestSample` refactor, which ties
   the fullsim input dir to its isospin treatment in ONE switch and thereby fixes the latent
   bug (`FullSimSampleUsesFourBeams(pp, is_test=true)` → 4 beams, as the pp24 test sample
   requires). **Verified before merging:** the committed runners reproduce this session's
   results BIT-IDENTICALLY (475 108 tree entries; anomaly-bin 8136 denom / 7289 num; all 6
   per-slice w_factors equal), and the plot macro rebuilds and reproduces the pp plateau
   (0.9588 ± 0.0094) from master. Follow-up work continues from master.

**Open questions / follow-ups:**
2b. **FUTURE TODO (user, 2026-07-16, not urgent — PbPb-only, qualitatively similar either way):
   decide whether real HIJING muons stay INCLUDED in the MC trig-eff muon sample.** Current
   choice (D4/§3.0) includes them (like the template-fit real/hadronic/fake axis); the
   alternative excludes them (consistent with reco-efficiency's Pythia-only construction).
   Affects only the overlay (pp provably identical either way — set-equality check, T5).
3. **Trigger-group question — ready to send, now SHARPENED by the r17663 result (R10).**
   The r17663 skim was run (2026-07-20) and **resolved the outcome tree: the dominant cause is
   the r16578 configuration, not HIJING occupancy**; its part (3) is answered by elimination.
   Include the R10 evidence *as corrected there*: the effect is **endcap-specific** (the barrel
   is clean), confined to **pT ≲ 6 GeV**, and present even though an **RoI-seeded HLT CB muon is
   already reconstructed** at the offline muon in ~90% of cases — but **the failing step is NOT
   localized** (our observables cannot separate the L1 threshold/coincidence from the HLT hypo
   or the chain matching). Ask them where it is; do not assert L1. Only the sending remains —
   it needs the user.
3b. **Full data-side dp/p cascade (U1 blast radius):** all data NTP outputs (pp24 nominal +
   pbpb23/24/25 nominal & trig-eff) and downstream (crossx, R_AA, template fits) still carry
   the old `fabs(dp/p)`; only the pp24+pbpb23 trig-eff references are being rerun in round 3
   (user decision 2026-07-16). Schedule the rest as one batch.
3c. **Minor hardening (INFO from the round-4 review):** `FitMCSinglesEffcy` writes
   `single_mu_effcy_pT_fit_mc[_medium_wp].root` with a basename that is IDENTICAL across pp24,
   overlay and noovl — only the directory separates them. Nothing was clobbered (dirs verified
   distinct), but it is the one place in the new wiring where a wrong `dir` would silently
   overwrite a sibling instead of erroring. Append `cfg.label` when that file is next touched.
4. **Double-matching cross-check for R4 — DONE incidentally in R8** (2026-07-16): `_0_01` and
   `_V3` per-muon matches are identical to the bare 0.02 match to the 4th decimal in every
   bin, both samples ⇒ double-matching excluded outright (no NTP flag needed; checked at the
   raw-NTUP level with NTP conventions mirrored).
5. **Quantify the Step-1 denominator asymmetry** (MC = real muons only; data = includes fakes)
   before quoting the Step-1 MC/data ratio as a number in the note (R5).
6. **Plateau-normalize ε_ΔR** before it is applied to crossx (unchanged precondition);
   current values = **R9** (round 3, Tight): pp 0.9569 ± 0.0094, overlay 0.8429 ± 0.0232.
   **Applies equally to the round-6 ε_ΔR^single** (overlay plateau 0.951 T / 0.943 M) — divide by the
   large-ΔR plateau so ε_ΔR^single(ΔR>1)=1 before it dresses the PbPb union linear terms (§2).
7. Re-run downstream when the full-stat productions arrive (all-centrality PbPb24-conditions
   overlay; lifts D2). **D3 is now LIFTED** — all 24 pp slices are in.

8. **★ NEXT TO-DO — MC CLOSURE TEST (round 8). Written 2026-08-03; deliberately NOT executed:
   the user will first check the round-7 results (the corrected-MC study and the ΔR-correction
   fits), because the closure design depends on both.**

   **What it is.** The MC sample is the only place where BOTH the unbiased denominator (every
   event stored, `StoreAllEvents`) and the per-muon trigger decision exist. So the whole
   correction chain can be closed on itself: take the **triggered** MC pairs, apply the analysis's
   per-pair trigger weight, and check that the result reproduces the **all-pairs** (no trigger
   requirement) MC yield — differentially, not just inclusively.
   - pp / 2mu4: `w = 1 / (ε₁ ε₂ · ε_ΔR^2mu4(ΔR))`, numerator condition = pair passes 2mu4.
   - PbPb / mu4 union: `w = 1 / (ε_ΔR^single(ΔR)·(ε₁+ε₂) − ε₁ε₂·ε_ΔR^cross(ΔR))` (§2),
     numerator condition = at least one leg mu4-matched.
   Closure variable = (weighted triggered yield) / (all-pairs yield), which must be **1** within
   uncertainties in every bin. Bin it in **ΔR, pair pT and pair η** — the inclusive ratio can
   close by construction while the differential one does not, and it is the differential
   behaviour that the correction is for. This test was **not available to either Run 2 reference
   analysis** (their overlay MC had no trigger simulation — R11), so there is no precedent to
   inherit; it is a genuine addition.

   **Two design choices that ROUND 7 MUST SETTLE FIRST (hence the wait):**
   (a) *Which ε goes in the weight.* Using the MC ε_MC with the MC-derived ΔR corrections is
   self-consistent and should close trivially — a code check, not a physics check. The physics
   question is whether the chain still closes when the **data-derived** ε^nc is used with the
   MC-derived ΔR corrections, i.e. the configuration the analysis actually applies. The
   corrected-MC study (contract item 5 / RW 9) decides how meaningful that variant is: if
   correcting the MC to the data efficiency provably leaves the ΔR corrections unchanged, the
   data-ε closure test is the right one and is not circular.
   (b) *Fitted vs binned ΔR correction.* Closure should use the **fitted/interpolated**,
   plateau-normalized corrections (contract item 6) — that is what the analysis will apply — so
   the test also validates the fit, not just the binned ratios. Which functional form is nominal
   is decided by item 6.

   **Acceptance:** closure consistent with 1 within the (round-7, conditional) errors in every
   (ΔR, pair pT, pair η) bin with meaningful statistics; any structured deviation is a systematic
   on the trigger correction and goes into `docs/systematic_uncertainties.md` §1a.

9. **Corrected-MC study (contract item 5, round 7)** — whether re-weighting each MC muon by
   `SF = ε_data/ε_MC(pT, q·η)` changes the ΔR correction terms. Feeds MC-closure design choice
   (a) above.

## Latest Stage

**2026-08-03 (round 7) — PAUSED mid-round at user request (tmux restart). RESUME HERE.**
Commits `f990173`..`b3b2a99` on master; working tree clean at the pause.

**DONE (committed):**
1. **Truth fiducial default** (§3.0c'), full product set regenerated — contract item 1. ✅
2. **ΔR-correction error bars fixed** (R12) — contract item 4 answered. ✅
3. **§3.5 sanity check** — `n_vtx` NTP propagation (both NTP passes re-run, entry counts
   identical: pp 19 563 129, overlay 95 389), `do_sanity` fill mode + plot block, run for
   pp_full / overlay / noovl × both WPs. **VERDICT: the forward anomaly is REAL** — see R14
   below. Contract item 2 ✅ except the sanity PNGs, which are not yet rendered (the plot macro
   was locked by the parallel fit task).
4. **Forward low-pT veto** `pT > 7 || q·η > −2` applied to Steps 2/3/4 and all fills re-run —
   contract item 3. ✅ (plots not yet re-rendered)
5. **`docs/systematic_uncertainties.md`** created — contract item 7. ✅
6. **MC-closure test written as Remaining Work 8**, NOT executed — contract item 8. ✅
7. **ΔR-correction fit stage delivered** (delegated): plateau ROOT file written by the measuring
   step, guard, fits, plots, driver — contract item 6. Files committed EXCEPT the
   `plot_mc_trig_eff.cxx` hunks (see below).

**REMAINING — pick up in this order:**
- **(a) MERGE `plot_mc_trig_eff.cxx`.** The fit task's version is staged at
  `/tmp/plot_mc_trig_eff_FROM_T5.cxx`; its six hunks (which write
  `<mc_dir>/dr_correction_plateaus_<label>[_medium_wp].root` from the Step-3/Step-4 blocks) are
  listed verbatim in `docs/tracking/_sub_drfit_1.md`. They must be merged into the committed
  version, which additionally carries the Step-1 sanity block. **Do not overwrite; do not
  `git merge` the worktree branch — its base commit is stale.**
- **(b) Re-render all plots** (`pipelines/run_mc_trigeff_round7.sh`, plot stage) so the sanity
  PNGs and the post-veto Step-2/3/4 PNGs exist, then re-run
  `pipelines/run_dr_correction_fits.sh` on the post-veto inputs.
- **(c) ⚠ USER DECISION — the full-sample plateau guard FAILS.** `|plateau − 1| > 0.1` in three
  pp_full cells, all in the top pair-pT bin (50–150 GeV): Step 3 η_pair[1.0,1.5) = 0.8597 ±
  0.0346 (worst) and η_pair[2.0,2.4) = 1.1082 ± 0.0119; Step 4 η_pair[0.5,1.0) = 0.8954 ±
  0.0091. The guard behaves exactly as specified (throws); what to DO about those cells is a
  physics decision. Note these values predate the veto refill and must be re-measured first.
- **(d) The corrected-MC study (contract item 5) was still running when the session paused** and
  was killed by the restart. **Its work is SAFE — staged into master by `8b5d269`, do not go
  looking in the worktree:**
  - `docs/tracking/_sub_corrmc_1.md` — its scratch doc, the source of truth for its findings.
    **Read this first**; do not assume the run finished.
  - `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_corrected.cxx` and
    `pipelines/run_mc_trigeff_corrected.sh` — usable as-is.
  - `docs/tracking/_staged_corrmc/{FillMCTrigEffHists.cxx,FitMCSinglesEffcy.cxx}.t4` — its
    versions of the two SHARED files. They are master@`3be0c5a` + its corrected-mode edits, so
    they LACK the forward veto and the vertex-comment fix. **3-way merge them; never copy over.**
  - Its ROOT outputs already exist on disk: `mc_trig_eff_hists_pp24_full[_medium_wp]_corrected*`
    (incl. `_step3`, `_step4` and a `_sfclosure` variant), `single_mu_effcy_pT_fit_mc_corrected*`,
    and the six `step{1,3,4}_corrected_mc/` plot directories under both trigger-efficiency plot
    roots. **They were produced BEFORE the forward veto landed ⇒ stale; regenerate.**
- **(e)** Then `/review-analysis-code` + `/review-plot`, `/wrap-up`, final summary.

### R14. Sanity-check VERDICT: the forward MC anomaly is REAL (2026-08-03, round 7)

pp24 FULL, Tight, μ⁺+μ⁻ summed, ε(mu4) per variant (`orig` = round-7 selection,
`vtx` = +1 reconstructed vertex, `ptm` = +|ΔpT|/pT^truth < 0.10, `both` = +both):

| q·η bin | pT 4–4.5 | 4.5–5 | 5–6 | 6–8 | 8–15 | 15–60 |
|---|---|---|---|---|---|---|
| (−2.4,−2.2) orig → both | 0.9013 → 0.9021 | 0.9095 → 0.9108 | 0.876 → 0.871 | 0.865 → 0.859 | 0.909 → 0.908 | 0.922 → 0.917 |
| (−2.2,−2.0) orig → both | 0.8971 → 0.8993 | 0.9113 → 0.9078 | 0.902 → 0.900 | 0.916 → 0.917 | 0.949 → 0.950 | 0.951 → 0.953 |
| (−2.0,−1.6) orig → both | 0.9011 → 0.9020 | 0.9264 → 0.9269 | 0.942 → 0.941 | 0.958 → 0.958 | 0.967 → 0.969 | 0.966 → 0.969 |
| **(+2.0,+2.2)** orig → both | **0.2610 → 0.2694** | 0.4818 → 0.4945 | 0.722 → 0.733 | 0.864 → 0.877 | 0.930 → 0.930 | 0.948 → 0.940 |

**Both split forward-negative bins are flat at ≈0.90 from the very first pT bin — no turn-on at
all — and neither extra requirement moves them by more than 0.006** (data in the same region
rises 0.45 → 0.92). The **mirror bin (+2.0,+2.2) in the same events shows a textbook turn-on
0.26 → 0.95**, which rules out any global reconstruction or selection artefact. ⇒ Neither
pile-up muons nor badly measured muons explain it; the anomaly is a genuine property of the
r16578 trigger configuration (consistent with R8/R10) ⇒ the §3.5 decision rule fires and the
forward low-pT veto is applied to Steps 2–4.

Pass fractions (pp24 FULL, Tight): `n_vtx == 1` keeps **5.65%** weighted (758 883 raw muons —
ample); `|ΔpT|/pT^truth < 0.10` keeps **99.82%**; both **5.64%**. Integrated ε moves only
0.7759 → 0.7737 (μ⁺). Overlay and noovl: the vertex requirement is a **no-op** (they have
exactly one track-bearing vertex in 100% of events — no pile-up) and the pT match keeps 99.8%,
changing ε by 0.0002.

---

**2026-07-28 (round 6) — COMPLETE. Step 4: single-leg ΔR correction ε_ΔR^single (§3.4) delivered;
Remaining Work 1 RESOLVED; both reviews PASS.** §2 union reformulated
(`P = ε_ΔR^single(ΔR)·(ε₁+ε₂) − ε₁ε₂·ε_ΔR^cross(ΔR)`); §3.4 = leg-level inverse weighting
(num = leg's own mu4 match ÷ ε_MC, denom = all legs, role-swap + SS/OS). `do_step4` in
FillMCTrigEffHists + Step-4 block in plot_mc_trig_eff → new `step4_dr_correction_singles/` dir
(step1/2/3 untouched), both WPs. **Deliverable (overlay/PbPb): plateau 0.951 T / 0.943 M, small-ΔR
rise → R4's 1.34. Validation (pp_full, full stats, NOT applied to pp): plateau 0.993 ≈ 1, per-cell
0.94–1.01 ⇒ inverse-weighting self-consistent.** /review-analysis-code PASS iter 1 + /review-plot
PASS iter 1 (all numbers MATCH, C1/C2/C3 clean). **Downstream remaining (RW6): plateau-normalize
ε_ΔR^single, then wire it into the PbPb crossx union weight** (application step, separate task).

---

**2026-07-22 (round 5) — COMPLETE. All 4 user changes delivered + validated + committed; both
reviews PASS.** #1 truth-seed revert (overlay HIJING excluded 99763→95389; pp no-op), #2 q·η
(−2.4,−2.0) split (11 bins, MC + data T&P), #3 Step-2 L1/HLT split via a re-skim (new per-muon
`muon_match_L1MU3V`; eff(L1)≥eff(chain), eff(HLT\|L1)≤1, closure ~0.95–0.98), #4 Step-3 pair-η
panels + plateau/fluctuation tables (pp plateau ~0.97–1.00 validates plateau→1). All 3 samples
grid-re-skimmed (13 tasks) → pp-full LGD farm v2 + overlay/r17663 download → re-NTP (truth-seeded+L1)
→ Fill/Fit/plot ×3×2WP → noovl three-way plot. /review-analysis-code PASS (`review-...-round5...md`);
/review-plot PASS iter 4 (pair-η super-title placement fixed). Commits `dea8c26`..`ec8ef7d` +
`fullsim_pp24_full_to_lgd.sh`→v2. **Carry-over (NOT round-5, pre-existing): §2 union-weight decision
(Remaining Work 1), trigger-group question send (RW 3), full data-side dp/p cascade (RW 3b);
plateau-normalize ε_ΔR before crossx (RW 6); crossx now also stale on the q·η-split forward-bin
trig-eff weight (#2 blast radius).**

---

**2026-07-22 (round 5, earlier) — GRID RE-SKIM COMPLETE (13/13 done), download/farm was IN PROGRESS.**
- All 13 grid re-skim tasks **succeeded** (verified via BigPanDA; poller exited 12:45). r17663 already
  downloaded + through NTP/Fill/Fit (dry-run; only its plot pending, needs pp_full+overlay 11-bin fits).
- **Download/farm launched (background):** pp-full LGD farm `fullsim_pp24_full_to_lgd.sh --no-devslice`
  (edited to v2: TASKS→51643327.. + `outds_for` VER_TAG v2; 75 stale v1 symlinks removed first,
  record backed up) → rucio add-rule the 6 v2 datasets to LGD + rebuild the symlink farm the NTP
  globs. overlay `grid_monitor.sh --mode overlay <6 v2 ids>` → 132 GB download+hadd (renames v1→.bak;
  DELETE the .bak after verify to reclaim ~132 GB; quota headroom ~175 GB soft was checked).
- **Overlay DONE downloading + verified (14:08):** all 6 v2 NTUPs 10000 entries + `muon_match_L1MU3V`
  present; **132.5 GB reclaimed** (deleted the 6 v1 `.bak_20260722` after verify). Overlay chain
  (re-NTP truth-seeded+L1 → Fill/Fit/plot ×2WP) delegated + running. **pp-full farm building**
  (v2 LGD rules replicate fast BNL→BNL; pTH14_24 already 25 symlinks).
- **Then:** pp_full re-NTP+full chain when farm done → re-run plot_mc_trig_eff("noovl",{T,M})
  (needs pp_full+overlay 11-bin fits) → /review-plot on all 3 samples.

**2026-07-21 (round 5) — Phase A/B (code + skim + submit + data-refill + dry-run). Four user changes
(#1 truth-seed revert, #2 q·η split, #3 Step-2 L1/HLT split via re-skim, #4 Step-3 pair-η).**
- **Phase B (all code) DONE + compiles clean + /review-analysis-code PASS iter 1** (0C/0W/1 INFO;
  log `review-analysis-code-20260721-214317-round5-mc-trig-eff.md`). The one INFO (L1-match φ has
  no 2π wrap, mirroring `muon_match_mu4roi`) was FIXED in the skim (φ folded into (−π,π]) since the
  L1 match is compared against the wrapping TDT chain match.
- **Phase A (grid re-skim for the L1 branch) SUBMITTED.** Skim rebuilt with `muon_match_L1MU3V`
  (φ-fix pp closure 0.9786→**0.9829**; HI test PASSED — `LVL1MuonRoIs` present in the r17618
  overlay AOD, no CHECK failure). **jediTaskIDs:**
  - **pp-full** `FullJuly2026.v2` (6 DSIDs 803015-020): **51643327, 51643336, 51643344, 51643353,
    51643363, 51643375**.
  - **overlay** `July2026.v2` (6 DSIDs 802776-781, r17618): 51643516, 51643528, 51643539,
    51643567, 51643586, 51643596.
  - **r17663** `July2026.v2` (802781 r17663): **51643610**.
  All 13 rc=0, zero errors. HI validation: P[L1\|chain]=0.948 (overlay) / 0.983 (pp), ordering
  L1>mu4roi>chain holds; lower absolute HI fractions expected (busier Dev_HI + occupancy).
  **IN FLIGHT (2026-07-21 ~22:20):** (a) grid task poller running (`/tmp/claude-101379/round5_task_poll.sh`,
  BigPanDA API, 13 tasks, re-invokes on all-terminal); (b) data-side #2 re-fill delegated + running
  (pp24+pbpb23 T&P graphs/fits, both WPs, 11-bin q·η split, Stage-5/6 SKIP_CONDOR, ≤3 threads,
  backs up 10-bin outputs to `.bak_pre_qeta_split_20260721`).
  **REMAINING when grid done:** re-farm pp-full to LGD (v2 datasets; the farm/NTP glob points at
  `pythia_fullsim_full_sample/`; watch quota — overlay download is ~132 GB) + download/hadd
  overlay (into `..._hijing_overlay_test_sample/`, `grid_monitor.sh --mode overlay <6 ids>`) +
  r17663 (`--mode no_overlay 51643610`) → re-NTP the mc_trig scripts (truth-seeded + L1) for all 3
  samples → Fill→Fit→Fill(step3)→plot_mc_trig_eff ×{pp_full,overlay,noovl}×{Tight,Medium} →
  /review-plot. **pp plots use the FULL sample.** The mc_trig NTP now REQUIRES nothing new but
  the L1 branch is bind-if-present (a v1 NTUP would warn + give empty L1/HLT hists).
  **STORAGE (checked 2026-07-21): data-fileset headroom ~175 GB to soft / ~445 GB to hard (halved).
  Overlay v2 download ~132 GB fits — rename v1→.bak (free), download v2, verify, then DELETE v1.bak
  promptly (per grid_monitor .bak rule). pp-full → LGD farm (no GPFS).** r17663 already downloaded
  (166.9 MB, done). **PROCESSING ORDER (from the dry-run): overlay + pp_full first → noovl plot LAST
  (its comparison series needs their 11-bin fits) → /review-plot.**
- Blast radius noted (out of scope, Remaining Work): #2 shifts the data single-muon trig-eff fits
  → crossx would need a rerun to stay consistent.

---

**2026-07-20 (round 4) — COMPLETE. r17663 no-overlay discriminator: submitted, run
end-to-end, verdict delivered (R10); both reviews PASS.**
- Grid skim jediTaskID **51600634** → 9 995 events, trigger present. Full chain run: NTP ×2 →
  Fill/Fit/Step3 × 2 WPs → **24 PNGs** in `plots/r17663_no_overlay_trigger_efficiency/mc_based/`
  (Tight in step dirs, Medium in `medium/`). **Separation verified: ZERO files under
  `pp_trigger_efficiency/` or `pbpb_trigger_efficiency/` modified.** New sample identity
  (`FullSimSampleType::noovl`) committed in `9916bb4`.
- **RESULT (R10): r17663 is OVERLAY-LIKE** — bin A (q·η∈(−2.4,−2), pT 4–6) ε = 0.5672 ± 0.0202
  vs pp24 0.8977 ± 0.0039 (16.1σ) and overlay 0.4997 ± 0.0144 (2.7σ) ⇒ **the dominant cause of
  the forward anomaly is the r16578 tag's trigger CONFIGURATION, not HIJING occupancy**; the
  barrel control (r17663 vs pp24 0.67σ, vs overlay 11.3σ) shows the overlay's barrel deficit IS
  occupancy and that the configuration effect is endcap-specific.
- **The failing step is deliberately NOT localized** (R10's "WHAT THIS DOES NOT SHOW"): an
  earlier draft blamed L1, which was backwards — `HLT_MuonsCB_RoI` is RoI-SEEDED, so r17663's
  0.902 argues *against* an L1-RoI deficit. R8's item-4 "decisive" and its verdict wording were
  hedged accordingly.
- Reviews: **/review-analysis-code PASS iter 4** (iter-1 CRITICAL = the L1 overclaim; iter-2/3
  WARNINGs = the claim surviving in Remaining Work, an unhedged "is downstream" → ≥77% bound,
  a stale Latest Stage, R8 internal inconsistency — all fixed). **/review-plot PASS iter 2**
  (iter-1 WARNING = a Tight/Medium plateau label swap in the doc; plots were always correct).
- **Carries forward (needs the user):** send the trigger-group question — now narrowed to the
  r16578 configuration and widened from "which L1 config" to "which STEP" (Remaining Work 3).

---

**2026-07-16 (round 3) — COMPLETE (except two user decisions pending).** All 5 items done:
U1 one-sided dp/p fixed everywhere + full trig-eff rerun (MC both samples/WPs + data pp24 &
pbpb23 refs incl. Medium refills; R9 headline table); U2 answered (truth fiducial gate =
`PassCuts_PythiaCore`, removed in `store_mc_trigger` only; `require_signal_cuts` untouched);
U3(a) HIJING-inclusion TODO recorded (Remaining Work 2b), U3(b) forward anomaly root-caused
to pp r16578 L1 endcap production config (R8; trigger-group question drafted); U4 Medium →
`medium/` subdirs (48 PNGs); U5 WP wiring verified NO BUG (plots identical because the mu4
response is nearly WP-blind — 94–95%/88–90% Tight∩Medium overlap, |Δε| ≤ 0.011).
Reviews: /review-analysis-code APPROVED iter 3; /review-plot APPROVED iter 1.
**Awaiting user:** (i) submit the r17663 no-overlay skim? (ii) send the R8 trigger-group
question? Plus the standing items: §2 assumption-(i) union-weight decision (Remaining Work 1)
and the full data-side dp/p cascade (Remaining Work 3b).

---

**Original round-3 plan (written before work):**

- **U1 — dp/p cut polarity (user: `fabs` IS WRONG; the cut is ONE-SIDED, `dp/p < thrsh`, data AND MC).**
  Archaeology DONE: data has used `fabs` since the repo's initial commit (`eac35a3`, 2022-11-02,
  `MuonNTupleFirstPass.C:28`) through every restructure to `DimuonDataAlgCoreT.c:599` — no recent
  commit introduced it; it was always wrong. MC was one-sided all along until round-2 D5
  (`221f49a`) copied data's `fabs` into `PythiaFullSimExtras.c:160`. Fix: one-sided in BOTH
  `DimuonDataAlgCoreT.c:599` and `PythiaFullSimExtras.c:160`; rewrite D5 (the "fix" direction
  reverses: data moves to MC's convention). `PowhegFullSimExtras.c:39` already one-sided.
  Blast radius: DATA-side generic muon cut change ⇒ every data NTP output (pp24 + pbpb23/24/25)
  and everything downstream (T&P, fits, crossx, template fits) is stale — rerun scope needs a
  user decision; MC trig-eff chain rerun is local/cheap and proceeds now.
- **U2 — explain the removed "truth fiducial gate" (answer only, no code change):** the old
  trigger-mode pair path passed pairs through `PassCuts_PythiaCore` (`PythiaAlgCoreT.c:741`,
  TRUTH pT>4 & truth |η|<2.4 on both legs); removed ONLY inside `store_mc_trigger`
  (`ProcessEventFullsimMCTrig`); nominal fullsim path untouched; the RDF-level
  `require_signal_cuts` reco-eff mode is a different layer, NOT touched.
- **U3 — (a) future-TODO: decide HIJING-muon inclusion for MC trig-eff (template-fit-like
  include vs reco-eff-like exclude); (b) DELEGATED /review-investigation: pp Step-1 forward
  bin (q·η −2.4..−2, pT 4–6) MC≫data flat-at-~0.9 vs overlay SAME bin MC<data — code bug vs
  simulation; includes AMI-tag comparison (pp e8599_s4162?_r16578 vs overlay e8599_s4614_r17618
  — verify actual tags) and the r17663 no-overlay 10k sample as contingency.
- **U4 — plots: Medium WP → `medium/` subdirectory per step** (currently `_medium_wp` filename
  suffix in the same dir).
- **U5 — DELEGATED verification: are the Tight plots actually Tight** (pp + PbPb; code path +
  input-file histogram bin contents; user observes Tight ≈ Medium visually).

---

**2026-07-14 (round 2) — COMPLETE. Reco-seeded truth-matched-real muons + both WPs; both
reviews PASS (/review-analysis-code iter 2, /review-plot iter 3).** Four result sets delivered
({pp, overlay} × {Tight, Medium}), 48 PNGs. Headlines: ε_ΔR^2mu4 = 0.9583 ± 0.0095 (T) /
0.9589 ± 0.0091 (M); ε_ΔR^cross = 0.8426 ± 0.0234 (T) / 0.8407 ± 0.0224 (M) — **WP-independent**
(pp χ²/ndf 0.01, overlay 0.03), so the correlation procedure is valid at both WPs.
**Two things carry forward:**
1. **(unchanged, still the open user decision)** the §2 assumption-(i) violation — the singles
   efficiency IS ΔR-dependent, so the **PbPb mu4 UNION** weight is at risk (pp/2mu4 product form
   is safe). Round 2 *strengthened* this: the small-ΔR enhancement survives on a truth-matched
   REAL-muon sample, so it is genuine L1 correlation, not fake/hadronic contamination.
2. **NEW blast radius (D5):** the `fabs(Δp/p)` fix is in shared fullsim code ⇒ **reco-efficiency,
   detector-response and template-fit MC are now STALE** (~2.6% of reco muons). Not rerun here
   (out of scope); must be rerun before those results are used. `PowhegFullSimExtras.c:25`
   carries the same bug (POWHEG fullsim obsolete — left, noted).

---

**Original round-2 PLAN (written before work), per §3.0 / D4 / D5:**
- **T1 — NTP (`store_mc_trigger` only).** Bind `muon_truth_id`, `muon_truth_IsPrimary`. Add
  `IsRealRecoMuon(reco_ind)`. Rebuild the muon list RECO-SEEDED over real muons; build pair
  trees from all (i<j) selected reco muons (not Pythia truth pairs). Fix the `fabs(Δp/p)` bug
  (D5). Report fake/hadronic/real fractions + how many real HIJING muons are gained.
  → /review-analysis-code (include §3.0 + D4 + D5 in the prompt).
- **T2 — Data Medium reference.** Rerun the data T&P (P2) + trig-eff fit with `isTight=false`
  for pp24 + pbpb23 → `*_medium_wp` outputs. (No data re-skim: the data trees are Medium
  supersets carrying a per-pair Tight flag.)
- **T3 — RDF/fit/plot both-WP wiring.** `FillMCTrigEffHists` / `FitMCSinglesEffcy` already take
  `use_tight_wp`; `plot_mc_trig_eff` must pick the **WP-matched data file** (`_medium_wp`) —
  currently it hardcodes the Tight data path. 
- **T4 — Rerun everything:** NTP ×4 (pp = all 24 beam×slice configs; overlay), then
  Fill→Fit→Step3→Plot for **2 samples × 2 WPs** = 4 full result sets.
- **T5 — Verify + review + docs + git.** /review-plot on the 4 plot sets; confirm the ΔR
  correlation procedure (Steps 2–3) is valid at BOTH WPs.

---

**2026-07-14 (round 1) — Steps S10–S15 ALL COMPLETE. Branch merged to master (`8edc4fb`).**
All six requested items done and reviewed (/review-plot APPROVED iter 3; the q·η investigation
returned NOT-a-code-bug with a KB-grounded evidence chain; merge verified bit-identical against
the committed runners). **One thing now sits with the user and blocks nothing else: the §2
assumption-(i) violation** — the singles efficiency IS ΔR-dependent (R4), which puts the
**PbPb mu4 union** weight at risk while leaving **pp/2mu4 safe by construction**. The Physics
Procedure §2 is deliberately unedited pending that decision (Remaining Work item 1).

**Next actionable step (pp, unblocked):** plateau-normalize ε_ΔR^2mu4 (0.9588 ± 0.0094) and
wire it into the pp crossx weight — roadmap Q4. The PbPb side should wait for the item-1
decision.

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
