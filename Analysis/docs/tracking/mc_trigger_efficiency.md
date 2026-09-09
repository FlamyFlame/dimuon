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

## Autonomy Contract (round 13 — **DONE 2026-09-08**; all 7 Done items met, both reviews PASS)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `fit_dr_corrections.cxx`: the **Step-3 `expo`** fit carries the two user-required
     parameter restrictions (small-ΔR plateau BELOW the large-ΔR plateau; turning point at
     ΔR > 0 / zero slope at ΔR = 0). Step 4 keeps its present limits — its physics is the
     OPPOSITE sign (R4: ε_single is ENHANCED at small ΔR).
  2. pp24 **DATA** single-muon mu4 T&P rerun on the new gap cut + the new coarse q·η top bin
     (`pipeline_pp_trig_eff.sh` stages 5–8, `SKIP_CONDOR=1`), both WPs.
  3. **MC pp_full** Steps 1–4 refilled + turn-ons refitted + replotted on the new
     `single_mu_fiducial_gap_cuts` AND the new pair-level `|η^pair| < 2.2`, both WPs
     (`run_mc_trigeff_round7.sh`, `SAMPLES="pp_full"`).
  4. **Step-3 ΔR-correction fits** re-measured and refitted for pp_full, both WPs, all
     methods × all sign series × all plateau modes, with the new restrictions
     (`run_dr_correction_fits.sh`, `STEPS=3`, no skips), plots regenerated.
  5. The most-forward pair-η cells are `(−2.2,−2.0)` / `(2.0,2.2)` in the 9-bin views and
     `2.0 ≤ |η^pair| < 2.2` in the 3-bin folded views — VERIFIED on the produced artefacts,
     not assumed from the source.
  6. `/review-analysis-code` on the fit-code change and `/review-plot` on the regenerated
     Step-1 + Step-3 plot sets, both PASS/APPROVED.
  7. Doc + INDEX updated, work committed.
- **NOT in scope (user, explicit):** propagating ε_ΔR to the pp24 cross-section; the HIJING
  overlay and `noovl` MC; PbPb data T&P; Step-4 fits; the MC closure.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Autonomy Contract (round 9 — **DONE 2026-08-11**; all 7 Done items met, both reviews PASS)

**Origin:** user request 2026-08-10 (session "MC-based trigger efficiency followup"). Focus is
**pp**; PbPb (HIJING overlay) code is updated for consistency only and is **NOT rerun**.
**r17663 (`noovl`) is explicitly OUT OF SCOPE for every item** — it was a one-time Step-1
cross-check and must not acquire Steps 2/3/4 or analysis-level decisions.

- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a
  passing small test, or one pipeline stage is NOT a stopping point.
- Done = all of the following, for **pp (pp_full)** at **both WPs** (Tight nominal + Medium):
  1. **Step-1 2D efficiency plots.** MC single-muon mu4 efficiency as a **2D map in
     (pT, q·η)**, on the SAME binning as the data single-muon mu4 efficiency 2D plots:
     (a) one PNG with **2 subplots — μ⁺ left, μ⁻ right**; (b) one PNG with **μ⁺ and μ⁻
     combined**. Written under `step1_singles_data_mc/`.
  2. **Step-3 ΔR-correction distribution plots restructured.** `step3_dr_correction/` gains
     **three subdirectories** — ΔR ∈ [0,1], [0,2], and full range — each containing **one PNG
     per pair-pT bin**, with **one subplot per pair-η bin**, and the **plateau value of that
     (pair pT, pair η) cell** (the very value used to normalize before the ΔR fit) drawn as a
     **horizontal line**. The **pair-pT/pair-η-INTEGRATED Step-3 plots are removed** (user: not
     helpful for the analysis).
  3. **Sign-separated Step-3 and Step-4 fits.** Under **each fit-function subdirectory** of
     `step3_dr_fit/` and `step4_dr_fit/` there are now subdirectories separating
     **sign-integrated** (`sign_intgr/`, = today's single-marker + single-fit plots, MOVED
     there unchanged) from **sign-separated** (`sign_sepr/`, same-sign and opposite-sign
     markers AND their two fitted curves overlaid on each subplot; e.g. red / dark red for
     opposite sign, blue / dark blue for same sign). The sign-integrated results are KEPT —
     the sign-separated ones are NEW fits, measured and fitted independently, including an
     independent plateau per sign where a plateau is needed.
  4. **MC statistics tables** in a new directory `mc_statistics_pt8bins/` under
     `pp_trigger_efficiency/`, on the canonical 8-bin log pair-pT binning (8→150 GeV).
     **Three sets × 2 CSV files = 6 CSVs**, one file for raw muon-pair COUNTS and one for
     SUMMED CROSS SECTION in each set:
     (i) same-sign and opposite-sign as 2 columns × 8 pair-pT rows (pair-η integrated);
     (ii) same-sign only, pair-pT as columns × pair-η as rows;
     (iii) opposite-sign only, pair-pT as columns × pair-η as rows.
  5. **PbPb (HIJING overlay) plotting/filling code updated for consistency** with every code
     change above, compiled and shown to be consistent — but **NOT rerun** (user: optional).
  6. **r17663 / `noovl` untouched** by every item above.
  7. Reviews (`/review-plot`, and `/review-analysis-code` for the C++/RDF changes), tracking
     doc + INDEX updated, commits, and a final summary listing **every new plot and file path
     produced this round**.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

**Carried forward, NOT part of round 9's Done:** the R23 OPEN finding (top-pair-pT plateaus
20–42 % high — ε_MC pT clamp vs coarse-q·η parameterisation, C4 investigation) remains open and
blocks *use* of the per-cell high-pair-pT corrections; round 9 produces plots and tables on top
of the existing round-8 measurement and does not resolve it.

## Autonomy Contract (round 8 — ACTIVE, opened 2026-08-04; re-read on every compaction)

**Origin:** advisor feedback relayed by the user 2026-08-04. Three physics points, then a
concrete work list. Round 7 is CLOSED except for its two reviews and wrap-up.

- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a
  passing small test, or one pipeline stage is NOT a stopping point.
- Done = all of the following:
  1. **Backup** of the current (fine q·η + gap-fallback) results: the whole
     `~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/` directories copied to
     `{pp,pbpb}_trigger_efficiency_fine_q_eta_bins_w_gap/`.
  2. **Canonical coarse pair-pT binning = 8 LOGARITHMIC bins, 8 → 150 GeV**, in
     `ParamsSet.h`; the two existing 4-bin definitions (`pair_pt_coarse_bins` 8–120 and
     `pair_pt_coarse_bins_pt150` 8–150) are **DELETED**, `N_COARSE_PAIR_PT_BINS = 8`.
     *Physics:* pair pT is the key observable; 4 bins smear its dependence out of the
     efficiency correction, which propagates into the pair-pT spectrum and R_AA.
     **Crossx consumers (`RDFBasedHistFillingPP.cxx`, `RDFBasedHistFillingPbPb.cxx`) and the
     plots they feed are an explicit FUTURE TO-DO — postponed, not done here.**
  3. **Run-3 coarse q·η binning** `{-2.4,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.2}` (10 bins)
     becomes `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap` (today's version merges
     (−0.5,0.5) into one bin — that is the only difference).
  4. **Gap fiducial cut applied** to every trigger-efficiency measurement, from
     `ParamsSet::single_mu_fiducial_gap_cuts` (already declared, **currently read by nothing**):
     reject a muon whose `q·η` lies in any of `(−1.20,−1.05)`, `(−0.06,0.06)`, `(2.20,2.40)`;
     reject a PAIR if EITHER leg is rejected. Applied to the single-muon tree AND the pair
     trees via a gap-cut mode in the ntuple-processing code (**default true**, distinct output
     suffix, nominal never clobbered), per the NTuple-Processing Provenance rule.
     **Data ntuple processing + nominal signal selection are explicitly OUT of scope
     (larger blast radius, needs a human decision) — future to-do.**
  5. **Single-muon efficiencies remade** (data AND MC) on the coarse q·η binning. Because the
     coarse binning INCLUDES the gaps, **the gap-region 2D fallback disappears** — the
     single-muon efficiency retrieval must be adjusted accordingly, with the fine-q·η path kept
     as an opt-in legacy mode. `mu4_mu4noL1` is not maintained (memory
     `feedback_mu4_mu4noL1_not_maintained`) ⇒ mu4 only.
  6. **Forward-edge decision plot:** single-muon mu4 efficiency in the most positive q·η bin,
     pT dependence with upper edges **(2,2.2), (2,2.25), (2,2.3), (2,2.4)** overlaid; 4 subplots
     — μ⁺ left / μ⁻ right, data top / MC bottom. Purpose: decide whether the forward cut must
     stay at 2.2 or can loosen to 2.25 / 2.3. ((2,2.4) is known to be bad.)
     *(Historical record of the REQUEST. 2.25 is not a bin boundary on the fine q·η axis and
     resolves to 2.26 — see the R22 correction box.)*
  7. **All remaining MC trigger-efficiency plots remade** on the 8 pair-pT bins, keeping the
     gap cut and the round-7 forward veto (`pT > 7 || q·η > −2`) on Steps 2/3/4.
  8. **Blast-radius flag:** everything outside trigger efficiency that the q·η-binning and
     gap-cut changes affect is enumerated in the doc as a future to-do (not executed).
  9. Reviews (`/review-analysis-code`, `/review-plot`), docs, INDEX, commit.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

**Point 3 of the advisor feedback (recorded, changes the R16 method):** weighting each MC muon
by `SF = ε_data/ε_MC` is **not a reasonable procedure**. The analysis proceeds with
**data-driven single-muon mu4 efficiencies + MC-derived ΔR corrections**. If the SF is to be
tested at all, it is applied as a direct multiplication `ε_MC × SF`, never as a per-muon weight.
⇒ the round-7 corrected-MC study (R16) keeps its *conclusion* (the ΔR corrections are
insensitive to the single-muon normalization) but its per-muon-SF **implementation is
deprecated**; do not extend it.

**⚠ Cross-doc coupling (read before touching the gap cut):** `muon_gap_cuts_acceptance.md`
(ACTIVE, a SIBLING SESSION is writing it) owns the derivation of these gap windows — its F6/F7
measure the structure, F8 declares the vector. `pp_trig_eff_highpt_jump.md` (ACTIVE, blocked on
a user decision) documents the live `w_trig = 0` bug whose candidate fix (d) *is* this fiducial
cut. Do not edit either doc from here; cross-reference them.

## Autonomy Contract (round 7 — CLOSING 2026-08-04; re-read on every compaction)
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
  - `ε_ΔR^single(ΔR) = ε_single(ΔR) / ε_single(plateau window)` is the ΔR correction on the
    single-leg (marginal) trigger probability — the deliverable of **Step 4 (§3.4)**,
    plateau-normalized so it is 1 for well-separated muons; and
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

with εᵢ = ε^nc(pTᵢ, q·ηᵢ), the **data-derived** single-muon mu4 efficiency of an isolated muon,
measured at ΔR > 0.8 so that it is already an isolated-muon value and the linear terms reduce to
εᵢ at large ΔR as required. **Consistency caveat (2026-08-04):** the MC ratios are now normalized
over ΔR ∈ [2, 3.5] (§3.3), which is NOT the ΔR > 0.8 region in which ε^nc is defined, so the two
references do not coincide by construction. Measured size of the mismatch on pp24 Tight: the
inclusive ε_ΔR^cross averaged over ΔR ∈ [1,2] is 0.97599/0.97159 = **1.0045**, i.e. 0.45 %
(up to ~3 % in individual cells). That residual is covered by the plateau-window systematic
(§3.3); it is not assumed away. MC supplies only the ΔR **ratios** ε_ΔR^single and ε_ΔR^cross — the
per-leg L1 over-efficiency cancels in each (§1, §4). The two assumptions this doc
measures/tests: (i) ~~the singles terms ε₁, ε₂ carry **no** ΔR dependence~~ **VIOLATED (R4) and
now CORRECTED**: the singles ARE ΔR-dependent, so the union's linear terms are dressed by
`ε_ΔR^single(ΔR)` (Step 4) rather than assumed flat; (ii) the correlation is a function of ΔR
alone (checked via kinematic binning).

**Decision note (2026-07-28, user — resolves Remaining Work 1).** The RW-1 union-weight
question is settled in favour of the "ΔR-dependent ratio correction on the linear terms" option:
`ε_ΔR^single(ΔR) = ε_single(ΔR)/ε_single(plateau window)`, measured by inverse weighting as Step 3
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

**THE PLATEAU WINDOW (authoritative; code constant in
`Analysis/Utilities/MCTrigEffPlateauWindow.h`, never retyped).** ΔR ∈ **[2, 3.5]** since
2026-08-04 (user decision; before that [1,4]). It applies identically to §3.3 and §3.4. Both
edges are cut, each for an independently measured reason:
- **Lower edge 1.0 → 2.0 — per-cell structure.** In several (pair pT, pair η) cells the
  ΔR ∈ [1,2] half sits significantly above the far half (worst: pT_pair[8,14) × η_pair[−0.5,0.5),
  +0.054 at 7.9σ). Cutting it improves the mean per-cell constant-fit χ²/ndf from 1.64 to 1.27
  (Step 3) and 1.66 to 1.34 (Step 4). Cutting only the upper edge does **not** fix this.
- **Upper edge 4.0 → 3.5 — the tail is an ARTEFACT, and it is diagnostic 1 firing.** The
  inclusive pp24 curve is flat from ΔR ≈ 0.6 to ≈ 3.4 (0.9702–0.9792) and then falls
  monotonically: 0.9571 (3.62), 0.9521 (3.88), 0.9341, 0.9262, 0.9156, 0.8965, 0.668, 0.390.
  The cause is **geometric**: with Δφ ≤ π, ΔR > 3.5 forces |Δη| > 1.54 and ΔR > 4 forces
  |Δη| > 2.47, so large-ΔR pairs push **both** legs into the endcaps — precisely the r16578
  forward-endcap region (R8/R10/R14) where ε_MC is badly parameterized. Run 2's analogous L1
  close-by-RoI correction stays at unity out to large ΔR and shows no such fall.
  Inclusive constant-fit χ²/ndf (Step 3 / Step 4): [1,4] 3.67/3.85, [2,4] 5.06/5.37,
  **[2,3.5] 1.08/0.76**.
- **Cost:** median per-cell relative stat error 0.95% → 1.40% (Step 3), 0.46% → 0.66% (Step 4).
  All cells stay usable (`fit_ok = 1` in 36/36).
- **The retired [1,4] window is kept ONLY as the normalization systematic**
  `|plateau[2,3.5] − plateau[1,4]|`, written per cell into the plateau ROOT file and reported in
  `plateau_guard_report.txt`. It is an uncertainty — nothing may normalize a curve by it.
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
  out; plateau-normalized → `ε_ΔR^single(ΔR) = ε_single(ΔR)/ε_single(plateau window)`. This
  dresses the union's linear terms (§2).
- **Method (inverse weighting — the leg-level analog of §3.3):**
  - **Object:** each muon **leg** of every MC reco pair (role-swap: both legs of a pair are
    probes; SS + OS trees summed — the trigger response is a per-muon detector property, blind to
    the pair charge product), binned by the **pair ΔR** (fine bins, as Step 3).
  - **Denominator:** all selected legs, unit weight × MC event weight. **NO trigger requirement
    on the leg or on its partner** — the partner only defines ΔR (§4; identical to §3.2).
  - **Numerator:** legs whose **own** per-muon mu4 match fires (never the event-level chain,
    never the partner's decision), each weighted `1/ε_MC(pTleg, q·ηleg)` with ε_MC the **MC-derived**
    §3.1 fit, TF1-Eval'd continuously (clamp to fit range, floor, cap — as §3.3).
  - **Ratio vs ΔR = ε_single(ΔR)**; plateau-normalize over the plateau window → ε_ΔR^single(ΔR).
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

### R19. Plateau window → ΔR ∈ [2, 3.5]; round-7 closed (2026-08-04)

**Two user decisions (AskUserQuestion) and one mid-review correction.**

**D-P1 — crossx binning: DEFER AND BUNDLE.** R17 recorded "crossx outputs are now stale"; that
was too broad. Measured: the ONLY crossx consumers of `pair_pt_coarse_bins` are the two
extended-mass 0–20 GeV histograms `h2d_crossx_minv_0_20_vs_pair_pt_coarse_{op,ss}_dsigma`
(`RDFBasedHistFillingPP.cxx:552-560`, `PbPb.cxx:1069-1082`) — the THStack / control-region study
inputs. The **nominal dσ/dp_T spectra were never affected** (they use `pT_bins_120`/`ptb150`).
Decision: refill those once, together with the pp + PbPb 23/24/25 refill that the sibling
session's `w_trig = 0` gap bug (`pp_trig_eff_highpt_jump.md`, `9f97818`) forces anyway.

**D-P2 — plateau window.** First measured near-half-vs-far-half and proposed [2,4]. The
`/review-analysis-code` reviewer found — and an independent recomputation confirmed — that this
was **half the story and [2,4] was inclusively WORSE than the [1,4] it replaced**. The inclusive
pp24 Tight curve is flat from ΔR ≈ 0.6 to ≈ 3.4 and then falls monotonically: 0.9571 (3.62),
0.9521 (3.88), 0.9341, 0.9262, 0.9156, 0.8965, 0.668, 0.390. Constant-fit χ²/ndf (S3 / S4):
[1,4] 3.67/3.85, **[2,4] 5.06/5.37**, [1,3.5] 0.90/0.74, **[2,3.5] 1.08/0.76**.
Cause is GEOMETRIC — Δφ ≤ π ⇒ ΔR > 3.5 forces |Δη| > 1.54 and ΔR > 4 forces |Δη| > 2.47, pushing
BOTH legs into the r16578 forward-endcap region (R8/R10/R14) where ε_MC is badly parameterized;
i.e. §3.3 diagnostic 1 firing. Run 2's analogous L1 close-by-RoI correction stays at unity out to
large ΔR. **User decision: nominal ΔR ∈ [2, 3.5]**, each edge cut for its own measured reason
(lower = per-cell structure, mean per-cell χ²/ndf 1.64→1.27 S3 and 1.66→1.34 S4; upper = the
artefact tail). Systematic = `|plateau[2,3.5] − plateau[1,4]|` per cell.

**Final numbers (Tight, canonical pair-pT binning as of round 7):**
| | inclusive plateau | flagged cells | failing | verdict |
|---|---|---|---|---|
| pp_full Step 3 | **0.9731 ± 0.0009** | 1 (pT_pair[41,120)×η_pair[1.0,1.5) = 0.8850 ± 0.0367) | 0 | PASS |
| pp_full Step 4 | **0.9867 ± 0.0004** | **0** | 0 | PASS |
| overlay Step 3 | 0.8669 ± 0.0218 | 4 (TEST sample, exempt) | — | reported |
0 unmeasurable cells anywhere. Window systematic: median ≈ 0.007, max 0.171.

**Code:** new `Analysis/Utilities/MCTrigEffPlateauWindow.h` is the single source of truth for both
windows, included by `plot_mc_trig_eff.cxx`, `plot_mc_trig_eff_corrected.cxx`,
`fit_dr_corrections.cxx` and `FillMCTrigEffHists.cxx`; **three silently-diverged literal copies
deleted** (the corrected macro and two sites in the filler still held [1,4]). Plateau ROOT file
gains `h_stepN_plateau_syst` / `h_stepN_plateau_syst_inclusive` (−1 = not evaluable, never 0);
guard report and Table B gained the window systematic. `FitMCSinglesEffcy` PNGs now carry the WP
token (Medium had been overwriting Tight — a deterministic clobber, not a race).

**Review:** 1 CRITICAL (the window, above) + 7 WARNINGs; every reported number verified MATCH by
independent recomputation. One correction to my own text: "per-cell errors grow ~2–5×" was
overstated (median 1.31×). Two findings — `plot_dr_correction_fits.cxx` draws only ΔR < 2 (so no
plateau-defining bin is on the canvas) and its hand-mirrored fit constants — live in a file a
SIBLING SESSION holds uncommitted; reported, not edited. **Carry both into round 8.**

### R20. Round 8 opened — advisor feedback (2026-08-04)

Backup of the pre-change state: `~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency_
fine_q_eta_bins_w_gap/` (311 + 526 PNGs, verified equal counts).

**Canonical binnings changed (done, smoke-tested):**
- `ParamsSet::pair_pt_coarse_bins` = **8 log bins 8→150 GeV**, generated by `fillLogBinningArray`
  (never retyped): 8, 11.54, 16.65, 24.01, 34.64, 49.97, 72.08, 104.0, 150.
  `N_COARSE_PAIR_PT_BINS = 8`; **`pair_pt_coarse_bins_pt150` DELETED** (it had no readers).
  ⚠ The coarse axis is no longer derived from `pT_bins_120` — **no interior edge coincides with
  the fine axis**; nothing may assume coarse ⊂ fine.
- `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap` = 10 contiguous bins
  `{-2.4,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.2}` (only change: (−0.5,0.5) split at 0).
- New `ParamsSet::PassSingleMuFiducialGap(eta, charge)` + `FiducialGapCutExpr(q_eta_expr)` —
  ONE definition for MC and data, deliberately NOT an extension of
  `PassSingleMuonGapCut`/`MuPairPassGapCut` (different object; those define every existing
  `_wgapcut` histogram and redefining them would silently change already-produced outputs).
  Verified: q·η −1.1/0.0/2.3 rejected; −2.3/1.1/−1.25/2.1 kept; charge folding correct.

**User decisions this round:** data tag-and-probe gap cut applies to the **PROBE ONLY**
(ε^nc is a per-muon efficiency evaluated only for non-gap muons; the tag is left uncut, and this
is the like-for-like match to MC Step 1, which has no tag). Per-cell ΔR corrections: **NO
fallback of any kind** — unmeasurable cells are **not plotted**, and the count of unmeasurable
cells is reported explicitly; a silent fallback would be misleading, and any future fallback must
be requested explicitly and tested against the original. Step-3/Step-4 pair-pT overlays must be
**split into two PNGs** (bottom four and top four pair-pT bins) — 8 lines in one subplot is
unreadable.

**Scout findings that change the work (three independent read-only sweeps, cross-checked):**
- **FIVE hardcoded copies** of the q·η bin list, not three: `SingleMuEffcyPtTurnOnFitter.cxx:45-51`,
  `plot_mc_trig_eff.cxx:409-415` **and `:416-418` (numeric)**, `plot_mc_trig_eff_corrected.cxx:424-430`
  **and `:431-433` (numeric)**.
- **`useCoarseQEtaBin` already exists** (`RDFBasedHistFillingData.h:175`) and the graph PRODUCER
  honours it, but the FITTER and the READER (`RDFBasedHistFillingData.cxx:628`) both hardcode
  *fine* — that asymmetry is the actual bug to fix. Also `PP.cxx:296`/`PbPb.cxx:642` skip the
  per-charge/per-centrality graphs in coarse mode; trigger efficiency needs them.
- **FATAL at 8 bins:** `plot_mc_trig_eff.cxx:1412-1413` `scol`/`smark` are 4-element vectors
  indexed WITHOUT modulo at `:1459-1461, :1476-1477` (index runs to 7) ⇒ undefined behaviour.
- Removing the gap fallback **fixes the `w_trig = 0` sentinel bug** of `pp_trig_eff_highpt_jump.md`.
- **Latent bug (verified):** `RDFBasedHistFillingData.h:206-207` — `output_generic_hists` and
  `output_gapcut_hists` have **no initializer** and are read uninitialized at
  `RDFBasedHistFillingData.cxx:188` by the trig-eff pipelines.
- **A live second-binning bug, independent of this work:** the data-area extended-mass macros
  (`plot_extmass.C`, `thstack_extmass.C`) label panels 8/15/27/50/150 while the histograms they
  project were filled with 8/13.75/23.63/40.62/120. They also index `PTLAB[5]`/`px[4][4]`/
  `colc[4]` by pair-pT bin ⇒ **out of bounds at 8 bins**.
- ⚠ **The crossx fillers have NO opt-out** — they read the same vector, so any recompile switches
  their axis to 8 bins. "Postponed" means *do not rerun crossx*, not *code unaffected*.
- ⚠ Editing `q_eta_proj_ranges_coarse_incl_gap` **also changes the single-muon RECO-efficiency
  plots** (`plot_single_muon_reco_effcy.cxx:103-105`, `..._r17618_vs_r17662.cxx:171-173`) — the
  "anything other than trigger efficiency" flag the user asked for. Do NOT touch the `_run2`
  vectors: `RDFBasedHistFillingData.cxx:743-748` key-matches against them.

Full inventories: scratch docs `_sub_qeta_binning_1.md`, `_sub_gapcut_wiring_2.md`,
`_sub_pairpt_binning_3.md` (Explore agents are read-only and could not write their own; the
orchestrator persisted their returned text verbatim).

**R20b. Round-8 wiring as implemented (2026-08-04).**

*Gap cut* — `ParamsSet::PassSingleMuFiducialGap(eta,charge)` + `FiducialGapCutExpr(q_eta_expr)`
(string builder, so the windows are never retyped into a JIT filter). Applied at the RDF stage,
NOT in the ntuple processing: `eta`/`charge` are already branches on both trees, every other
muon-level fiducial cut in this chain already lives in the RDF macro, and the DATA gap cut cannot
go in the data NTP (out of scope) — putting the MC cut in the MC NTP would split one fiducial
definition across two stages and make MC/data asymmetric. Sites:
`FillMCTrigEffHists.cxx` `sel_single` (Steps 1 + sanity), `sel_pair_legs` (Steps 2 + 4), and the
**separate Step-3 selection string** (Step 3 does not use `sel_pair_full` — the one easy to miss);
data probe-side in `RDFBasedHistFillingPP.cxx` / `PbPb.cxx` right after `q_eta2nd` is defined.
Output token `_qeta_fid`, so no gap-cut output can overwrite a no-gap-cut one.

*Coarse q·η* — `useCoarseQEtaBin` now defaults **true** (fine = legacy opt-in). The
producer/reader asymmetry is fixed: `EvaluateSingleMuonEffcyPtFitted` reads `q_eta_proj_ranges`
(set from the same flag) instead of the hardcoded fine list, via new file-scope mirrors
`s_q_eta_proj_ranges` / `s_use_coarse_q_eta` (the evaluator is a STATIC member and cannot see
instance state). **The 2D fallback is now confined to the legacy fine path**; on the nominal
coarse path a missing turn-on THROWS instead of returning the −1 sentinel — which is the
`w_trig = 0` pair-dropping bug of `pp_trig_eff_highpt_jump.md`, now structurally impossible.
`MakeAndWriteSingleMuonTrigEffPtGraphs` branched on `useCoarseQEtaBin` and produced NO per-charge
graphs in coarse mode; it now branches on `isForSoumya` (obsolete) so per-charge graphs are always
made — otherwise every efficiency lookup would have silently missed.
**All five hardcoded copies of the q·η bin list are gone** (data fitter, plot macro ×2 incl. the
numeric ones); each now derives from `CommonEffcyConfig`. Legend pads derive from the bin count
instead of a literal `cd(12)`.

*Plots* — the 4-element `scol`/`smark` and `ptcol`/`ptmark` palettes (indexed WITHOUT modulo up
to 7 at 8 bins ⇒ undefined behaviour) extended to 8 distinct entries and modulo-guarded. Step-3
and Step-4 pair-pT overlays **split into `_lowpt` / `_highpt` PNGs** (user: 8 overlaid series in
one pad are unreadable), y-range computed per half. Pair-pT labels `%.0f` → `%.1f` (the new edges
11.54 / 16.65 / 24.01 … are misstated by up to 0.5 GeV at `%.0f`).

*Latent bug fixed* — `RDFBasedHistFillingData.h` `output_generic_hists` / `output_gapcut_hists`
were UNINITIALISED and read by the trig-eff pipelines; now `true` / `false`.

*Runs so far:* pp24 data trig-eff refilled on the coarse binning + probe gap cut (60 per-charge
`_py_*_divided` graphs) and refitted (40 TF1s). MC chain (fill → fit → Steps 3/4 → sanity, 3
samples × 2 WPs) launched.

**⚠ OPERATIONAL TRAP HIT (2026-08-04):** the PbPb23 data trig-eff fill with
`EnableImplicitMT(8)` was **SIGKILLed (exit 9), almost certainly OOM** — PbPb books the
centrality-binned families on top of everything else. It left a **450-byte output file** and the
wrapping shell still reported exit 0, so an exit-code check would have called it a success and the
next stage would have fitted an empty file. This is the standing "validate ARTEFACTS, not exit
codes" rule, and it now has a second failure mode: **check the output SIZE, not just existence.**
Re-run with `EnableImplicitMT(2)`.

**Round-8 run status at last checkpoint:** pp24 data (fill + fit) DONE; MC fill/fit/Steps 3-4/
sanity DONE for all 3 samples × 2 WPs; pp_full plot set DONE on the 8 pair-pT bins (the
`_lowpt`/`_highpt` split PNGs are present and the plateau table shows the new edges
`pTpair[8.0,11.5) … [104.0,150.0)`). Overlay/noovl plots and the ΔR fits are BLOCKED on the
PbPb23 data reference. **Confirmed cost of the finer binning, exactly as predicted:** the pp_full
Step-3 top cell reads `1.0879 ± 1.0879` — a 100 % relative error, i.e. genuinely unmeasurable.
Per the user's instruction such cells must NOT be plotted and their COUNT must be reported; no
fallback of any kind. That reporting is still TO DO.

### R21. ★ THE GAP HYPOTHESIS IS CONFIRMED — plateaus move to ≈1 (2026-08-04, Tight)

The advisor's hypothesis was that the plateaus sat away from 1 because the **gap regions** were
handled by an unfitted 2D fallback with large binning fluctuations. Cutting those muons and
fitting a CONTIGUOUS coarse q·η binning does exactly what was predicted:

| inclusive plateau | round 7 (fine q·η + gap fallback) | round 8 (coarse q·η + gap cut) |
|---|---|---|
| pp_full Step 3 | 0.9731 ± 0.0009 | **0.9945 ± 0.0009** |
| pp_full Step 4 | 0.9867 ± 0.0004 | **0.9974 ± 0.0004** |
| overlay Step 3 | 0.8669 ± 0.0218 | **0.9614 ± 0.0249** |
| overlay Step 4 | 0.9477 ± 0.0090 | **0.9924 ± 0.0114** |

**⚠ STALE NUMBERS (noted 2026-08-11):** the round-8 column above was measured with the forward q·η
edge at 2.20. After R22's edge moved to 2.30 the same quantities read **0.9921 ± 0.0009** (Step 3)
and **0.9962** (Step 4). The *conclusion* of R21 — that the gap handling was what held the plateaus
away from 1 — is unaffected; only the third-decimal values are superseded.

**Medium WP reproduces it** (pp Step 3 0.9948 ± 0.0008, Step 4 0.9977 ± 0.0004; same failing/
flagged counts 12/5 and 2/7) ⇒ the improvement is WP-independent, as a detector-geometry
effect must be.

The overlay Step-3 move (0.867 → 0.961) is the striking one: the long-standing ~13 % offset was
**mostly the gap-region fallback**, not a genuine failure of the inverse-weighting closure. The
plateau-normalization systematic (§1a) shrinks accordingly and must be re-derived on these values.

**The cost is per-cell, and it is real.** With 8 × 9 = 72 cells (Tight):

| | unmeasurable | flagged (0.10–0.15) | FAILING (>0.15) | verdict |
|---|---|---|---|---|
| pp_full Step 3 | 0 | 4 | **12** | FAIL (FULL sample) |
| pp_full Step 4 | 0 | 7 | **2** | FAIL (FULL sample) |
| overlay Step 3 | 18 | 1 | 45 | reported, exempt (TEST) |
| overlay Step 4 | 15 | 5 | 34 | reported, exempt (TEST) |

The failures concentrate in `pT_pair ≳ 50 GeV`, where the old single [41,120) column is now split
four ways and pp Pythia runs out of yield (e.g. Step 3 `pT_pair[72.1,104.0) × η_pair[−2.4,−2.0)`
= 1.1060 ± 0.0059). The guard did its job: it reported every cell and exited non-zero, and the
fits exist only under an explicit override, for inspection. **Open decision for the user:**
whether 8 bins is right, or the top edge should come down from 150 GeV — the inclusive and
mid-pT corrections are excellent either way.

### R22. Forward q·η edge scan — 2.30 is safe, 2.40 is not (2026-08-04)

Plateau ε (pT > 8 GeV), all four panels agreeing (`forward_qeta_edge_scan/`):

| upper edge | data μ⁺ | data μ⁻ | MC μ⁺ | MC μ⁻ | probes (data μ⁺) |
|---|---|---|---|---|---|
| 2.20 | 0.9073 | 0.9177 | 0.9339 | 0.9306 | 3 378 |
| 2.26 | 0.9004 | 0.9107 | 0.9268 | 0.9245 | 4 318 |
| **2.30** | 0.9026 | 0.9101 | 0.9281 | 0.9251 | **4 898** |
| 2.40 | 0.8492 | 0.8563 | 0.8851 | 0.8804 | 6 220 |

2.20 → 2.30 costs ~0.5 % in plateau efficiency and recovers **45 % more probes**; 2.40 costs ~6 %.
⇒ the forward window can be loosened to `{2.30, 2.40}`. **APPLIED** — `ParamsSet.h` now carries
`{2.30, 2.40}`.

> **Correction (2026-08-13, R31).** The row originally headed **2.25** was mislabelled: the fine
> q·η axis is 0.02-wide above 2.20 (edges … 2.24, **2.26** …), so 2.25 is *not* a bin boundary and
> the projection necessarily stopped at **2.26**. The efficiencies were always right — only the
> stated edge was wrong (and the probe counts were rounded to 3 380 / 4 320 / 4 900 where the
> integrals are 3 378 / 4 318 / 4 898). Header and counts corrected above; the pp figure was
> regenerated on the fixed
> labelling (2026-08-13), and every label is now written from the achieved bin edge so it can no
> longer disagree with the data. **Central values unchanged.**

*New macro* `plotting_codes/trig_effcy/mc_based/plot_forward_qeta_edge_scan.cxx` — the forward-edge
decision plot (2×2: μ⁺/μ⁻ × data/MC, upper edges REQUESTED 2.20/2.25/2.30/2.40 — the drawn ladder
resolves to 2.20/**2.26**/2.30/2.40, see the correction box above). **It requires NO-GAP-CUT
inputs** (the nominal chain removes the forward window outright, which would draw four identical
curves); it asserts the q·η > 2.2 region is populated and refuses to run otherwise. It needs one MC
pass with the gap cut disabled (`_nogapcut` output) — **produced 2026-08-04**,
`pythia_fullsim_full_sample/mc_trig_eff_hists_pp24_full_nogapcut.root`, and it is what the MC
panels of the pp figure are drawn from.

### R23. Both reviewers run — one OPEN physics finding that REFUTES an earlier claim (2026-08-06)

`/review-plot` FAIL (11 CRITICAL, 11 WARNING) and `/review-analysis-code` FAIL (3 CRITICAL,
14 WARNING). Every reported number verified MATCH in both. Logs:
`.claude/logs/review-plot-20260806-073032-round8-mc-trigeff-gapcut.md`,
`.claude/logs/review-analysis-code-20260806-205802-round8-gapcut-binning.md`.

**★ OPEN, BLOCKING BEFORE THESE CORRECTIONS ARE USED — the top pair-pT plateaus are a real
effect, NOT a statistics effect.** This doc previously attributed the 12 failing pp Step-3 cells
to the top bins "running out of yield" at 8 bins. **That is wrong.** Measured:
`pT_pair[104,150) × η_pair[0.5,1.0)` = **1.2198 ± 0.0037**, weighted RMS 0.0063 over 4 ΔR bins —
**59σ from 1**; also 1.2020 ± 0.0111 (18σ) and 1.1922 ± 0.0088 (22σ). Twelve cells across the top
three pair-pT bins are 20–42 % high, coherently across ΔR, with per-cell scatter far below the
offset. Round 7 had **0** failing pp cells and verdict PASS; round 8 has 12. Two candidate
mechanisms, both testable:
- **(a) the ε_MC pT clamp.** `MCEffEvaluator::Eval` clamps pT into the TF1 fit range [4, 60] GeV.
  In a 104–150 GeV *pair* bin both legs sit near 50–75 GeV, so ε is read at the 60 GeV edge;
  ε₁ε₂ biased low pushes the inverse-weighted ratio UP. Same family as the out-of-range TF1 trap
  of `pp_trig_eff_highpt_jump.md`. Test: extend the single-muon turn-on fit range above 60 GeV.
- **(b) the new COARSE q·η bins** are up to 0.5 wide and straddle the barrel/endcap transition,
  so ⟨ε(bin)⟩ is a poorer parameterisation than the retired fine bins — and the quality of that
  parameterisation is exactly what §3.3 diagnostic 1 measures.
⇒ **C4 investigation required.** Do NOT adopt the per-cell high-pair-pT corrections until it is
resolved; the inclusive and mid-pT corrections are unaffected.

**Fixed this round (both reviews):** the `interp` branch set `fit_ok = 1` unconditionally,
publishing 20 overlay cells (e.g. plateau 0.0078 ± 0.0051, a ×128 inflation) that the guard called
unmeasurable and the plot refused to draw — the consumed artefact disagreed with both the guard and
the figures; `fit_ok` never tested the plateau ERROR, so the 1.0879 ± 1.0879 single-ΔR-bin cell was
published as usable; the y-range loop screened differently from the drawing code, so undrawn cells
pushed every PNG to the 3.0 cap; the plateau-window label printed `[2,4]` for a `[2,3.5]` window
(`%.0f` rounding onto the RETIRED edge); fit headlines misstated every pair-pT edge; the Step-3
palette drew 8 series in 4 colours; the inclusive canvas headline collided with its panel label;
the `interp` annotation was prose with an undefined `R_p`; 8 stale PNGs; and the nested `medium/`
subdirectories (a LOCAL `wp_dir` in `plot_mc_trig_eff.cxx` bypassed the retired helper).
**New (user):** Step-3/4 fit plots draw fit-domain points black and **excluded points blue**.

**Open WARNINGs carried forward:** the 2.30 gap-cut ↔ coarse-q·η-top-edge tie is comment-only and
only throws when a muon lands in the mismatch window (needs a startup assert); the tag-and-probe
filename is retyped in 8 places and two of them load the COARSE file for a fallback only the
LEGACY path uses; `MCTRIGEFF_PAIRPT_4BIN=""` flips the C++ and shell tests apart (the shell would
then `rm -f` the NOMINAL fit file); the PNG-count check still expects 5 not 9; `PlotSubdir()` is
dead and contradicts the adopted layout; `run_pipeline2_pbpb.C` defaults to the legacy fine path,
re-arming the −1 sentinel; **the coarse q·η change silently re-bins the single-muon RECO-efficiency
plots (9→10 bins, forward edge 2.2→2.3) under unchanged filenames — a blast-radius item missing
from the Contract**; `plot_muon_q_eta_spectrum.cxx` still asserts the removed fallback/sentinel
behaviour on a figure; and stale "2.2"/"NOT YET WIRED IN" comments contradict the live code on the
BLOCKING gap-cut vector.

**★ PHYSICS ITEM FOR THE BUNDLED CROSSX RERUN:** ε^nc is now measured on non-gap PROBES, but the
signal selection still has NO gap cut, and the coarse q·η bins STRADDLE the gaps. A bin-averaged
turn-on measured without gap muons will therefore be applied to signal muons that ARE gap muons,
whose true efficiency is lower ⇒ under-correction. Either apply the gap cut to the signal selection
in the same rerun, or measure ε^nc without the probe cut in the straddling bins. State the choice
in the Physics Procedure.

### R24b. The `_pt4bin` variant trees are STALE as of round 9 (2026-08-11)

`{pp,pbpb}_trigger_efficiency/mc_based_pt4bin{,_medium}/` were last written on the round-8 code.
They have none of the round-9 structure — no `sign_intgr/` / `sign_sepr/`, no
`step3_dr_correction/dr_*` subdirectories — and their `fit_report*.txt` still carry the `%.0f`
pair-pT labels that round 9 replaced everywhere else (24 files). They were deliberately NOT
regenerated: round 9 is scoped to the canonical 8-bin axis, and rewriting only their reports would
leave a half-migrated tree, which is worse than a uniformly old one. If the 4-bin comparison is
still wanted, it needs one full `MCTRIGEFF_PAIRPT_4BIN=1` pass of the chain, not a patch.
(The `*_fine_q_eta_bins_w_gap/` trees are the round-8 backup snapshot and are frozen by design.)

### R26. ★ OPEN — the NOMINAL `expo` form cannot describe the measured shape in the most
### populated cell (2026-08-11, plot review; pre-existing, exposed by round 9)

**Not a round-9 regression — a property of the fit form chosen in `440e4a0`** ("expo is NOMINAL").

In `pT_pair ∈ [8.0, 11.5) GeV × η_pair ∈ [−0.5, 0.5)` — the single largest cell of the sample (that
pair-pT bin holds 144 721 same-sign + 776 604 opposite-sign pairs, and (−0.5, 0.5) is the widest
pair-η bin) — the sign-integrated `expo` fit has **χ²/ndf = 48.18**. The measured plateau-normalized
points **dip** to 0.908 at ΔR ≈ 0.225 and then **rise above 1** to 1.064–1.114 across
ΔR ∈ [0.425, 0.975], while the fitted curve is pinned at 1.000 over that whole region
(A = +0.0535, λ = 0.0269 — the fit has collapsed onto its own asymptote). The monotone form
`1 + A·exp[−(ΔR/λ)^p]` **cannot represent a dip-then-overshoot shape at all**, so the delivered
correction is low by 8–11 % over more than half the fit domain in the most populated cell. Because
`usable` in `fit_dr_corrections.cxx` has **no χ² term**, the cell carries `fit_ok = 1` and every
consumer will apply it.

**It is systematic, not one bad cell — and NEITHER parametric form escapes it.**

**⚠ CORRECTION (2026-08-11, caught by the plot reviewer).** An earlier version of this entry
carried a χ²/ndf table claiming `polyu_fixedRp` had max 1.8 / 2.8 against `expo`'s 48.2 / 34.0, and
recommended switching the nominal on that basis. **That table was wrong** — it came from a parse
that mis-read the column on every row carrying a trailing `AT LIMIT:` annotation, so those rows
entered as zero. The corrected numbers, taken from the reports' OWN summary lines and confirmed by
an independent re-parse (Tight, sign-integrated, `pp_full`):

| step / method | mean | median | max | cells > 5 | cells > 10 | n |
|---|---|---|---|---|---|---|
| Step 3 / **expo (NOMINAL)** | 3.695 | 1.896 | 48.2 | 7 | 6 | 71 |
| Step 3 / polyu_fixedRp | 3.645 | 1.778 | 35.3 | **10** | 7 | 71 |
| Step 4 / **expo (NOMINAL)** | 12.465 | 2.749 | 641.1 | 13 | 6 | 72 |
| Step 4 / polyu_fixedRp | **23.411** | 2.716 | **932.4** | **19** | 10 | 72 |

⇒ **`polyu_fixedRp` is NOT globally better, and the earlier recommendation to adopt it is
WITHDRAWN.** It edges `expo` on the Step-3 mean/median/max and on the worst cell (35.3 vs 48.2 —
both still bad fits), but it is clearly WORSE on the Step-4 aggregate (mean 23.4 vs 12.5, max 932
vs 641) and has more badly-fitted cells than `expo` in both steps. The one place its advantage is
real and verified is the **at-limit count: 0 against `expo`'s 34/73 in Step 4** (Step 3: 0 vs 7).

The finding that survives, and it is the more important one: **both parametric forms fail badly in
a substantial minority of cells** — 7–19 cells of ~71 above χ²/ndf 5, with Step-4 tails reaching
the hundreds — while `usable` has no χ² term, so every one of them is published with `fit_ok = 1`.
The problem is the parametric family, not the choice between these two members of it.
(`interp` has χ²/ndf ≈ 0 by construction — it passes through every point — so it is not a
fit-quality statement, but it also has the fewest rejected same-sign cells, 11/72 vs 25 and 31.)

**Three options, all physics decisions ⇒ USER (raised 2026-08-11):**
1. **Add a χ²/ndf screen to `usable`** in `fit_dr_corrections.cxx`, so a badly-fitted cell is
   rejected on the canvas and for consumers exactly as an unphysical `f(0) < 0` cell already is.
   This is the change most clearly supported by the evidence, and it is method-agnostic. Cost: it
   removes the correction in those cells rather than fixing it — including the most populated one.
2. **Use `interp` for the per-cell corrections** (no functional-form bias by construction; also the
   best same-sign coverage), keeping a parametric form only for the inclusive summary curve.
   Cost: an interpolation carries no smoothing, so it propagates per-bin statistical noise.
3. **Adopt a form that can describe a dip-then-overshoot.** Most work, and only worth it if the
   overshoot is physical rather than a residual of the plateau normalization — which is itself
   worth establishing first.
**Nothing changed here** — the nominal fit function is not something to swap autonomously, and on
the corrected numbers there is no obvious swap to make.

### R27. The no-plateau-correction Step-3 fit variant (2026-08-11, round 10, user request)

**What it is.** `step3_dr_fit/` now carries two complete parallel trees,
`plateau_corrected/` (the previous results, relocated one level deeper, verified **byte-identical**
— 408 files compared, 0 differing) and `no_plateau_correction/`, each holding
`<method>/{sign_intgr,sign_sepr}/`. Step 4 is corrected-only (user) but moved to
`step4_dr_fit/plateau_corrected/` so the two steps keep the same shape. 486 PNG + 126 TXT.

**The fit.** Raw, un-normalized ε_ΔR fitted with a **free additive baseline**
`f(ΔR) = C + A·exp[−(ΔR/λ)^p]` (`polyu_fixedRp` likewise; `interp` pins its flat branch to the last
measured knot). **Fit domain unchanged and confirmed in code: ΔR ∈ [0, 1] only** — the `zoom`
histogram, 20 bins of 0.05, centres 0.025…0.975; the wide-bin histogram (which carries the [2, 3.5]
window) is read by the PLOT stage only. **So C is fixed exclusively by data in [0, 1] and no
plateau-window bin enters any fit** — which is the whole point.

**The plateau guard does not apply here**, and that is a consequence, not a choice: nothing is
normalized by the plateau, so `|plateau − 1| > 0.15` cannot disqualify a cell and an unmeasurable
far-ΔR plateau does not prevent a fit. Cells fitted, Step 3, sign-integrated, 72 cells:

| WP / method | corrected | no correction | only corrected | only no-correction |
|---|---|---|---|---|
| Tight `expo` | 59 | **68** | 1 | 10 |
| Tight `polyu_fixedRp` | 59 | **66** | 0 | 7 |
| Tight `interp` | 59 | **70** | 0 | 11 |
| Medium `expo` | 60 | **67** | 2 | 9 |
| Medium `polyu_fixedRp` | 60 | **66** | 0 | 6 |
| Medium `interp` | 60 | **70** | 0 | 10 |

The "only no-correction" cells are exactly those the corrected mode threw out on its plateau. The
three "only corrected" cells fail the positivity screen — the un-normalized fit dips below zero on
[0, R_p] (worst: Medium `pT_pair[50,72.1) × η_pair[1.5,2.0)`, min f = −0.671).

**Inclusive cell, `expo`, sign-integrated:**

| | A | λ | p | C | χ² / ndf |
|---|---|---|---|---|---|
| Tight, corrected | −0.1881 ± 0.0024 | 0.2530 ± 0.0026 | 2.481 ± 0.083 | — | 174.22 / 17 = 10.25 |
| Tight, no correction | −0.1848 ± 0.0030 | 0.2513 ± 0.0031 | 2.520 ± 0.093 | **0.9906 ± 0.0015** | 173.18 / 16 = 10.82 |
| Medium, corrected | −0.1853 ± 0.0023 | 0.2594 ± 0.0025 | 2.536 ± 0.083 | — | 191.52 / 17 = 11.27 |
| Medium, no correction | −0.1832 ± 0.0028 | 0.2588 ± 0.0030 | 2.551 ± 0.091 | **0.9920 ± 0.0014** | 191.36 / 16 = 11.96 |

Nested-model sanity check passes: **χ² DROPS when C is freed** (174.22 → 173.18), so the extra
parameter is not fighting the data; χ²/ndf rises only because ndf goes 17 → 16.

**★ The inclusive cell is NOT where the two methods disagree — the per-cell ones are.** Inclusively
C reproduces the unused [2, 3.5] plateau to 0.15 % (0.9906 vs 0.9921). Per cell they part company
exactly where the user predicted: `pT_pair[8.0,11.5) × η_pair[−2.4,−2.0)` gives **C = 0.849 ± 0.013
against a plateau of 0.997** — a 15 % difference in the normalization applied to the small-ΔR region
of interest, in a gap-adjacent forward η bin. That is the quantitative form of the concern that
motivated this variant, and it is why the comparison has to be made per cell.

**Presentation.** `C (fitted plateau) = v ± e` is printed on every panel and on the inclusive
canvas, per sign in `sign_sepr`, with the `(at limit)` marker where it binds; the drawn equation and
a header line define it. The [2, 3.5] value is deliberately **not** drawn on a no-correction canvas,
so two different numbers called "plateau" never appear together; it stays in the report. The y-axis
title differs between the modes (ratio-to-plateau vs the measured efficiency itself).

**Open trade-off to be aware of:** the per-pair-pT canvases share one y range across the directory
(user request, round 10 item 1), and in the no-correction tree that range is stretched by a few
wild high-pair-pT cells, so a typical panel uses ~20 % of its frame. A percentile-based shared range
would fix it without giving up comparability — **not done, user's call.**

**Next (NOT done, user's stated intent):** compare MC closure with and against the plateau
correction. That is Remaining Work 8's closure test, now with two correction variants to run it on.

### R28. sign_sepr marker/colour scheme + the fit-rejection standard (2026-08-12, round 10 item 3, user)

**Change (plot only, no physics).** In `plot_dr_correction_fits.cxx` the sign-separated series
were `kBlue+3` / `kRed+3` markers on `kBlue` / `kRed` curves, with SQUARES for same sign and
CIRCLES for opposite sign — two very dark, similarly-valued hues that the user could not tell apart
where the series overlap at small ΔR, which is the whole point of the comparison. Now
(`plot_dr_correction_fits.cxx:280–283`):

| series | markers + error-bar lines | fitted curve | marker shape (in fit domain / outside) |
|---|---|---|---|
| same sign | `kBlue+1` | `kBlue` | triangle 22 / 26 |
| opposite sign | `kRed+1` | `kRed` | circle 20 / 24 |

The shape now carries the sign as well as the hue (greyscale- and print-safe). `sign_intgr` is
untouched. Regenerated with `SKIP_MEASURE=1 SKIP_FIT=1 SAMPLES=pp_full STEPS="3 4" WPS="tight medium"`
(plot stage only — the reviewer verified from the filesystem that every `fit_report*.txt` still
carries its 08-11 mtime, so no fit and no measurement was re-run): Step 3 = 12 `sign_sepr`
directories (2 WP × 2 plateau modes × 3 methods) × 18 PNG = 216, plus 108 `sign_intgr` PNG;
Step 4 = **162** PNG. No artefact failures. The pipeline still exits 2 on the **pre-existing** R23
plateau-guard FAIL (`pp_full/tight/step3`, `pp_full/medium/step3`) — unchanged by this edit.

**One defect found and fixed on the way (pre-existing, both steps).** The defining equation carried
a negative kern around the conditioning bar, `" |#kern[-0.20]{#DeltaR})"`. A kern is an
**absolute**-font-size nudge and this macro draws two canvas sizes, so on the 900 px INCLUSIVE
canvases it pulled the `Δ` on top of the `|` and `P(pair passes 2mu4 |ΔR)` rendered as one garbled
glyph, while the 1556 px 3×3 canvases looked fine — i.e. it was invisible to anyone who checked
only the per-pair-pT figures. Both kerns removed (Step 3, and the Step-4 analog
`"|#kern[-0.15]{leg}"`), with the reason recorded in the code so it is not "tightened" back. Step 4
was regenerated for that reason alone. This defeated the definition of the analysis-invented symbol
ε_ΔR on every inclusive canvas of every method since the equation was introduced.

**`/review-plot` APPROVED at iteration 2** (log
`.claude/logs/review-plot-20260812-175808-sign-sepr-marker-colours.md`). The reviewer confirmed
separability at full resolution in the dense small-ΔR overlap, that the error bars really carry
`c_mark` (what looks black in a downscaled view is dense blue), that filled-inside/open-outside
survives the shape change for both series, and re-verified the sign mapping at source
(`FillMCTrigEffHists.cxx:714-715`, `sign1 → ss_ → "same sign"`).

**Left alone deliberately, flagged by the reviewer:** the sibling trees
(`pbpb_trigger_efficiency/*`, `*_pt4bin*`, `r17663_no_overlay_*`, `*_fine_q_eta_bins_w_gap/*`) still
carry BOTH the old kern and the old square/dark-shade markers — and they are also still on the
pre-round-10 layout (no `plateau_corrected/` level), so they need a re-run, not a restyle.

**The fit-rejection standard (user question, answered from the code).** One gate,
`h_step3_fit_ok`; consumers must require `== 1`. A cell is rejected if ANY of:
1. **unmeasurable cell** — `DrCorrPlateauUsable`: `plateau > 0`, `≥ 0.5`, `err < plateau`, and
   ≥1 plateau bin (`dr_correction_sample_cfg.h:42`, `fit_dr_corrections.cxx:754`). No fit is even
   attempted.
2. **too few informative points** — points with non-zero error in ΔR ∈ [0,1] must be
   `≥ nfree + 2` (`:791`).
3. **fit did not converge** — `TFitResult::IsValid()` and `ndf > 0` (`:945`).
4. **unphysical** — `f(ΔR) ≤ 0` anywhere on 201 samples over [0, R_p] (`:961`); a trigger
   correction can never be ≤ 0.
5. **plateau off unity** — `|plateau − 1| > 0.15` (`kPlateauGuardTol`, `:976`). **Dropped in
   `no_plateau_correction`**, where nothing is divided by the plateau.
`interp` has no parameters, so only 1, 4 (`interp_physical`) and 5 apply (`:860–875`).

NOT part of the gate, and worth restating: **χ²/ndf never enters `usable`** — this is the open
finding **R26**. Distinct from all of the above is the *sample-level* plateau guard
(`|plateau−1| > 0.10` FLAGGED, `> 0.15` FAIL) which is fatal only for the sign-integrated series of
a FULL production and is what makes the pipeline exit non-zero today.

Current Step-3 `expo`, pp24 FULL, Tight, of 72 cells: plateau_corrected 32 rejected (same sign) /
17 (opposite sign); no_plateau_correction 27 / 4 — i.e. criterion 5 is what rejects most same-sign
cells in the corrected tree, exactly the pathology the no-correction variant was built for.

### R29. Round 11 — the gap-window rerun, and what it moved (2026-08-13)

**Trigger.** `single_mu_fiducial_gap_cuts` central crack window `{-0.06,+0.06}` → `{-0.10,+0.06}`
(commit `1f143e6`; derivation and cost in `muon_gap_cuts_acceptance.md` F11a — the depletion is
one-sided, N(neg)/N(mirror pos) = 0.509 pp / 0.578 PbPb in (−0.10, 0.00) and ≈1 elsewhere, so the
symmetric window cut healthy acceptance on the + side and left the depleted slice in). The window
is compiled into the RDF/JIT selection at FILL time, so this needed a REFILL, not a replot.

**What was rerun, in order** (each stage's inputs strictly older than its outputs):

| # | stage | driver | outputs |
|---|---|---|---|
| 1 | data T&P, pp24, Tight | `SKIP_CONDOR=1 pipeline_pp_trig_eff.sh` | fill 02:15→02:17, fit + plots 02:17 |
| 2 | data T&P, PbPb 23/24/25, Tight (+ the inverse-weighted data ΔR stage) | `SKIP_CONDOR=1 pipeline_pbpb_trig_eff.sh` | fills 02:26 / 02:30 / 02:36, fits 02:37, plots → 02:38 |
| 3 | data T&P, **Medium** WP, pp24 + PbPb23 | **`run_data_trigeff_medium_wp.sh` (NEW, see below)** | 02:39 / 02:44 |
| 4 | MC chain: refill Steps 1–4, refit turn-ons, replot, 2D maps, statistics tables — 3 samples × 2 WPs | `run_mc_trigeff_round7.sh` | 02:44→03:01 |
| 5 | ΔR corrections: guard + refit + replot (both plateau modes, 3 methods, 3 sign series, both views) | `SKIP_MEASURE=1 run_dr_correction_fits.sh` | 03:0x→03:2x, 54 nominal fit files |

The NTuple processing was NOT rerun and did not need to be — the cut is applied at the RDF stage
(D6). Verified directly on the refilled pp24 data file: the probe q·η spectrum is **exactly zero**
in [−0.100, +0.060) and populated on both sides (1972 entries in [−0.110,−0.100), 2228 in
[0.060,0.070)), i.e. the asymmetric window is live and is the ONLY thing that changed.

**A hole this exposed, now fixed: the data pipelines never ran the Medium WP.** Both
`pipeline_{pp,pbpb}_trig_eff.sh` leave `isTight` at its default `true`, so every rerun refreshed
the Tight tag-and-probe file and silently left the `_medium_wp` one behind — and that file is not
a spare: `plot_mc_trig_eff.cxx` with `use_tight_wp = false` reads it as the data reference the
Medium MC is compared against. Before this round the Medium data files dated from 2026-08-04, so a
Medium comparison would have put **new-window MC against old-window data**, with no error and no
warning. New driver `pipelines/run_data_trigeff_medium_wp.sh` runs exactly the two WP-dependent
stages (RDF fill with `isTight=false`, turn-on fit with `wp_suffix="_medium_wp"`) for pp24 and
PbPb23 — the two files anything consumes. **PbPb 24/25 Medium T&P files have never existed and
still do not**; nothing reads them today, so they were not invented here.

**What it moved (pp24 FULL, Tight, Step-3 `expo`, plateau-corrected), old → new:**

| quantity | before (2026-08-11) | after |
|---|---|---|
| plateau-guard FAILING cells (\|plateau−1\| > 0.15) | 12 | **11** |
| rejected cells, same sign / opposite sign (of 72) | 32 / 17 | **32 / 17** (unchanged) |
| inclusive χ²/ndf, sign-integrated | 10.248 | **10.278** |
| inclusive χ²/ndf, same sign / opposite sign | 1.308 / 9.803 | **1.323 / 9.860** |

Small, as expected for a window that removes ~1.6 % of pp muons in a narrow q·η slice, and in the
harmless direction (one fewer cell fails the guard). **R23 is NOT resolved by this** — the guard
still returns `verdict: FAIL` on the top pair-pT cells, so the pipeline still exits non-zero, and
R25/R26 are untouched: the sign-dependence and the χ²-blind `fit_ok` gate are properties of the
fit, not of the fiducial region.

**Blast radius, stated and checked:** the gap cut is still NOT in the signal selection
(`muon_gap_cuts_acceptance.md` F10.3), so **crossx and R_AA did not move** and
`signal_selection_change_impact.md` was not triggered. `_pt4bin` was deliberately excluded (already
STALE per R24b).

### R30. Post-rerun plot review: two header-text defects fixed, four findings left open (2026-08-13)

`/review-plot` ran once over the FINAL state (layout restructure + gap-cut rerun together) and
**APPROVED at iteration 3**. The physics half is the part worth recording: C1/C2 clean on the
refilled results (every departure maps to an already-open finding — R4's forward rise, R25's sign
split, R23's top-pair-pT plateaus); the widened window is visible as an empty q·η column from
−0.10 to +0.06 measured off the rendered Step-1 2D map and, correctly, invisible in the pair-η
binning whose edges straddle it; the dependency chain `ParamsSet.h → .so → data T&P (both WPs) →
turn-on fits → MC hists → plateaus → 108 fit ROOTs → 810 PNGs` was verified node by node to be in
order with **zero stale products**; and C3 against the Run 2 muon-trigger KB entry is consistent in
both the ρ_ΔR turn-on scale and the HLT_mu4 barrel/endcap magnitudes (including the *sign* of the
MC-over-data barrel excess, which is our separately measured 1.21/leg, not a new effect). The two
dR views were shown **pixel-identical in their annotation strips** for all 192 nine-panel pairs and
all 24 inclusive pairs.

**Fixed (plot-macro only — no refill, no refit, so nothing the closure thread reads moved):**
1. The header's last line pitch was 18 px for a 17 px font, so `eq_line1`'s descenders sat on
   `eq_line2`'s cap height — the `a_2`/`a_3` subscripts of the `polyu_fixedRp` equation printed on
   top of `u ≡ max(0, 1 − ΔR/R_p)`, and `interp`'s line 2 on line 3. Pre-existing, on ~2/3 of the
   fit canvases, and the fit equation is a MANDATORY on-canvas element. Pitch → 22 px with
   everything below it moved down by the same 4 px.
2. **A regression of round 10 that the reviewer did NOT catch, found by direct inspection:** on the
   900×826 INCLUSIVE canvases the panel label's superscripts (`p_T^{pair}`, `η^{pair}`) were
   CLIPPED by the pad edge. Moving the label into the strip put its baseline at NDC 0.965, where
   the baseline fits but the glyph *ascent* does not; the shorter 3×3 sub-pads happened to survive,
   which is exactly why a "label present and legible?" check passed it. Label → 0.940, top margin
   tracking it. *Lesson: position text from its tallest glyph, not its baseline — and crop-and-zoom
   the rendered PNG, because a superscript with its top sliced off still reads as present.*

**Left open deliberately (all INFO):**
- An **unconstrained** fit published with `fit_ok = 1` (`A = −0.153 ± 1.7`, `C = 0.706 ± 2.3`, i.e.
  300 % errors) whose **χ²/ndf is 1.45** — evidence for the open **R26** decision that a χ² screen
  alone would not catch this class; a relative-parameter-error screen would.
- ~~`forward_qeta_edge_scan/` still shows the **old symmetric window** (2026-08-04). Regenerating it
  needs a dedicated `MCTRIGEFF_NO_GAPCUT=1` fill — user's call.~~ **Resolved for PbPb by R31
  (2026-08-13): the PbPb figure is now DATA-ONLY, so no `_nogapcut` MC fill is needed there. The pp
  figure is untouched and still carries the 2026-08-04 MC panels.**
- `plot_mc_trig_eff.cxx` defines ε_ΔR as "P(both μ fire | ΔR)" where the fit plots say "P(pair
  passes 2mu4 | ΔR)". The latter is correct; fixing it rewrites the plateau files the concurrent
  closure thread now reads, so it is NOT done here.
- Pre-restructure orphan directories in the PbPb fit-plot trees (housekeeping, PbPb not refilled).

### R31. PbPb forward q·η edge scan is now DATA-ONLY: by-year + by-centrality (2026-08-13, user)

**What changed and why.** The PbPb `forward_qeta_edge_scan/` figure was a 2 × 2 with *data on top,
HIJING overlay on the bottom*. The overlay is still a **test sample**, and in a single forward q·η
slice (2.0 < q·η < 2.4) it has far too few muons to say anything about where the edge belongs — the
bottom row carried no information and invited over-reading. **User decision: drop the MC entirely
from the PbPb figure** and spend the panels on the two axes PbPb actually has and pp does not.
`single_mu_eff_forward_qeta_edge_scan.png` was **deleted**; two PNGs replace it (same directory,
`…/plots/pbpb_trigger_efficiency/mc_based/forward_qeta_edge_scan/`):

| PNG | grid | content |
|---|---|---|
| `…_by_year.png` | 3 rows × 2 cols | centrality-integrated; **μ⁺ left, μ⁻ right**; rows = PbPb **2023 / 2024 / 2025** |
| `…_centrality.png` | 3 rows × 2 cols | **2023+2024+2025 summed AND μ⁺ + μ⁻ summed**; one panel per centrality bin |

Centrality bins are read from **`ParamsSet::ctrbins = {0,5,10,20,30,50,80}`** — the histogram key
`_ctr<lo>_<hi>_` *and* the panel label are both built from that one vector, so they cannot drift
apart (CLAUDE.md BLOCKING binning rule). Its 6 bins are exactly the 3 × 2 grid the user asked for.
**pp is untouched** (still the 2 × 2 data/MC canvas), and no other PbPb trigger-efficiency plot was
regenerated.

**Inputs.** The un-suffixed no-gap-cut tag-and-probe outputs
`pbpb_20{23,24,25}/histograms_real_pairs_pbpb_20YY_single_mu4_fine_q_eta_bin.root`
(the gap-cut variant carries `_qeta_fid`; a gap-cut input would have nothing above q·η = 2.2 and
the macro's `AssertForwardRegionPopulated` guard refuses to run on one). Numerator / denominator =
`_2mu4_sepr` / `_mu4_sepr`, i.e. the `{"_2mu4","_mu4"}` entry of `trigs_pair` for `trigger_mode=1`.

**Two candidate-edge ladders, and why.** The requested ladder was 2.20 / 2.25 / 2.30 / 2.40, but
neither figure could draw it as written, for two independent binning reasons — both verified
against the ROOT files by the reviewer:

1. **2.25 is not a bin boundary.** The fine q·η axis is 0.10 wide up to 2.20 and 0.02 wide above
   it (edges … 2.20, 2.22, 2.24, **2.26** …), so a projection asked to stop below 2.25 stops at
   2.26. The by-year ladder is therefore **2.20 / 2.26 / 2.30 / 2.40**.
2. **The per-centrality q·η binning is not uniform in centrality.** It is registered per bin in
   `hist_binning_map` (`binName + "_ctr<lo>_<hi>"`) and the peripheral bins are deliberately
   coarser: **184** bins for 0–5 / 5–10 / 10–20 / 20–30 %, **102** for 30–50 %, **61** for
   50–80 %. 2.30 exists on the 184-bin grid but **not** on the 61-bin grid, whose neighbouring
   boundaries are 2.28 and 2.36. Scanning "2.30" in all six panels would have meant [2.0,2.30) in
   four panels, [2.0,2.32) in 30–50 % and [2.0,2.36) in 50–80 % — one legend entry over three
   different physical windows, exactly the CLAUDE.md §Binnings failure. The centrality figure
   therefore scans the boundaries **common to all six panels**, computed at run time by
   `CommonQEtaBoundaries`: **2.20 / 2.28 / 2.36 / 2.40**.

Every legend entry is written from the edge the axis actually delivered (`ResolveQEtaRange`), so
no label can state a q·η the drawn data does not have. **The adopted edge 2.30 is drawn as a curve
in the by-year figure; in the centrality figure it is bracketed by the 2.28 and 2.36 curves**, and
both canvases carry its value in the headline, read from
`ParamsSet::single_mu_fiducial_gap_cuts` rather than retyped.

**Numbers — plateau ε (p_T > 8 GeV) and probe count D, per candidate upper edge.**
All 48 values below were independently recomputed from the three `_fine_q_eta_bin` files by the
reviewer: **48/48 MATCH**.

*By year (all centralities), ladder 2.20 / 2.26 / 2.30 / 2.40:*

| | [2.0,2.20) | [2.0,2.26) | [2.0,2.30) | [2.0,2.40) |
|---|---|---|---|---|
| 2023 μ⁺ | 0.8514 (2571) | 0.8534 (3295) | **0.8565 (3763)** | 0.8133 (4767) |
| 2023 μ⁻ | 0.8649 (2605) | 0.8638 (3327) | **0.8691 (3758)** | 0.8212 (4765) |
| 2024 μ⁺ | 0.9181 (1783) | 0.9083 (2280) | **0.9065 (2567)** | 0.8545 (3195) |
| 2024 μ⁻ | 0.9053 (1796) | 0.8961 (2300) | **0.8986 (2623)** | 0.8435 (3271) |
| 2025 μ⁺ | 0.9037 (5045) | 0.9018 (6488) | **0.9050 (7409)** | 0.8575 (9324) |
| 2025 μ⁻ | 0.9102 (5102) | 0.9037 (6570) | **0.9054 (7417)** | 0.8493 (9351) |

*By centrality (23+24+25, μ⁺+μ⁻), common ladder 2.20 / 2.28 / 2.36 / 2.40:*

| centrality | [2.0,2.20) | [2.0,2.28) | [2.0,2.36) | [2.0,2.40) |
|---|---|---|---|---|
| 0–5 % | 0.8678 (3723) | **0.8665 (5100)** | 0.8558 (6291) | 0.8171 (6846) |
| 5–10 % | 0.8830 (3394) | **0.8845 (4641)** | 0.8729 (5713) | 0.8254 (6279) |
| 10–20 % | 0.9018 (5174) | **0.8983 (7136)** | 0.8926 (8787) | 0.8510 (9549) |
| 20–30 % | 0.9052 (3122) | **0.8989 (4271)** | 0.8921 (5262) | 0.8530 (5694) |
| 30–50 % | 0.9171 (2653) | **0.9115 (3627)** | 0.9087 (4426) | 0.8676 (4782) |
| 50–80 % | 0.9032 (682) | **0.8966 (928)** | 0.8921 (1159) | 0.8565 (1254) |

**Physics read.** The adopted forward edge **2.30 is confirmed in PbPb on both new axes**, and the
conclusion is the one R22 reached in pp:
All Δε and probe-gain numbers below are quoted **relative to the tightest edge, 2.20**, and are
plateau quantities (p_T > 8 GeV), i.e. read straight off the two tables above.
- **Loosening from 2.20 is free.** By year, 2.20 → 2.30 changes the plateau by |Δε| ≤ **0.0116**
  (max 2024 μ⁺; *positive*, i.e. an improvement, in both 2023 charges) while the probe sample grows
  by **+44.0 – +46.9 %**. By centrality, 2.20 → 2.28 costs |Δε| ≤ **0.0066** in every bin (again
  positive in 5–10 %) for **+36.1 – +37.9 %** probes, and even 2.36 costs only **0.0084 – 0.0131**
  for **+66.8 – +69.9 %**.
- **2.40 is not free.** It costs **0.038 – 0.064** by year and **0.047 – 0.058** by centrality —
  a 4–6 % absolute loss, at least 5× and up to ~40× the 2.28/2.30 step, and visible as a red
  turn-on *shifted* to higher p_T rather than merely scaled down (muons falling off the endcap
  trigger acceptance, not a flat inefficiency).
- **The effect is charge-blind** (μ⁺ and μ⁻ agree to within **0.0135** at every edge in every year,
  i.e. ~1–1.5 %) and **occupancy-blind**: the 2.28-vs-2.20 step does not grow toward central
  collisions — it is smallest in 0–5 % (−0.0013) and largest in 50–80 % (−0.0066) — so the
  loosened window is as safe in central as in peripheral events.
- A separate, edge-independent trend is visible: the plateau itself rises from 0.8678 ± 0.0055
  (0–5 %) to 0.9171 ± 0.0054 (30–50 %), i.e. Δ = 0.049 ± 0.008, a **~6σ** effect — not a
  fluctuation. It is not monotone in occupancy (50–80 % sits at 0.9032 ± 0.0113, *below* 30–50 %).
  It is present at *every* candidate edge, so it factorizes out of the edge comparison and does not
  bear on the decision. Note it is **not** covered by the Run 2 reference:
  `atlas_run2_muon_trigger.md` §11 records *no significant* central-vs-peripheral difference for
  HLT_mu4 in Pb+Pb, but that statement is η-integrated whereas this is the outermost 0.4 units of
  q·η, so the two are not in contradiction — flagged here as unexplained rather than expected.
- 2023 sits ~0.04 below 2024/2025 at every edge; this is a year-level offset, not an edge-scan
  effect (the *shape* of the edge dependence is identical in all three years). See the caveat
  below before reading it as physical.

**⚠ Caveat on the inputs (stated, not fixed).** The three `_fine_q_eta_bin` files were filled at
different times: **2023 on 2026-07-21, 2024 and 2025 on 2026-07-08** — i.e. before the one-sided
Δp/p fix `71fcf1c` (2026-07-16) reached them, whereas 2023 is after it. The muon-pair trees
themselves are all current (regenerated 2026-08-13 02:19), so a refill would make the three years
consistent. **It was NOT done**: `TrigEffPlotterPbPb::configureDataFiles` reads these exact
un-suffixed filenames, so refilling them would silently change the input of the
`pbpb_trigger_efficiency/mu4/no_corr/` plot set, which the user explicitly fenced off for this
task. The effect on this figure is expected to be small (Δp/p is a probe-quality cut and largely
cancels in an efficiency ratio) and it cannot fake the 2.40 cliff, but the ~0.04 year-to-year
offset above should not be read as purely physical until the refill is done.

**Note on the output directory.** Both PNGs stay in
`plots/pbpb_trigger_efficiency/**mc_based**/forward_qeta_edge_scan/`. That parent name is now a
misnomer for a data-only figure, but it is kept deliberately: renaming would split the figure from
its pp twin and orphan the `mc_based_medium` sibling. **The PbPb content of that directory is
data-only.**

*Macro* `plotting_codes/trig_effcy/mc_based/plot_forward_qeta_edge_scan.cxx` — PbPb branch
rewritten (`PlotPbPbByYear` / `PlotPbPbByCentrality`), pp branch factored into `PlotPP`, plus three
shared fixes that also touch pp: `ResolveQEtaRange` (labels written from the achieved bin edge),
`CommonQEtaBoundaries` (the centrality common grid), and the adopted-edge value on the headline
read from `ParamsSet`. **pp was regenerated on the fixed labelling** — its efficiencies are
unchanged, only the "2.25" legend entry became the truthful "2.26" (see the R22 correction box).
Run: `root -l -b -q 'plot_forward_qeta_edge_scan.cxx+("pbpb")'` / `...+("pp")`.

**Reviewed** `/review-plot`, log `.claude/logs/review-plot-20260813-183530-pbpb-forward-qeta-edge-scan.md`.

### R32. Round 12 — the THIRD Step-3 plateau mode: no plateau correction with the LAST TWO
### pair-pT BINS MERGED (2026-08-17, user request; pp + overlay, both WPs, all methods, all signs)

**What it is.** `step3_dr_fit/` now carries THREE parallel trees. Beside `plateau_corrected/`
(nominal) and `no_plateau_correction/` (R27) there is
**`no_plateau_correction_last2ptbins_merged/`**, holding the SAME raw fit — `f(ΔR) = C + A·exp[−(ΔR/λ)^p]`
with the baseline C free and fitted on ΔR ∈ [0,1] alone — measured on a pair-pT axis whose **last two
cells are merged**: `p_T^pair ∈ [72.1,104)` + `[104,150)` → **one cell [72.1,150) GeV**, so the grid
is **7 × 9 = 63 cells** instead of 8 × 9 = 72. Everything else is untouched (same samples, same gap
cut, same pair-η binning, same plateau window, same methods, same sign series, same two ΔR views,
same styling).

**Why (physics).** The top two cells of the 8-bin log axis run past where pp Pythia has yield. They
are where the plateau guard fails (R24/R27), where the ΔR fits are noise-dominated, and where the
closure collapses (`mc_trig_eff_closure.md` R1b, η^pair ∈ [−2.4,−2.0) × p_T^pair ∈ [72,104)).
Merging them buys statistics in exactly one place and leaves every measuring cell alone.

**It is NOT a new binning** (.claude/CLAUDE.md §Binnings). `ParamsSet::pair_pt_coarse_bins` and the
FILLED histograms are unchanged; a merged cell is the two source bins **projected together**
(num/denom/errA/errB summed *before* the ratio, `DrCellRatioRange`), which is numerically identical
to having filled a 7-bin axis. The grouping is derived in ONE place from the source histogram's own
axis (`dr_correction_pt_groups.h`) — no edge value is typed in a macro, a report or a plot — and the
variant is opt-in and suffixed (`_nocorr_ptmerge`), so it can neither overwrite nor be mistaken for
the other two trees. **8-bin axis only:** with `MCTRIGEFF_PAIRPT_4BIN` set the C++ throws and the
driver skips it with a note (the `_pt4bin` variant already merges the top cells by construction).

**The plateau in this mode.** As in `nocorr`, nothing is normalized by it and the guard's fatal tier
is off. But the plateau map on disk describes the UN-merged grid, so for the merged cells it is
**re-measured** from the `full_vs_pt_eta_*` histograms in the same file with the SAME estimator and
the SAME window — `PlateauFromRatio` was moved out of `plot_mc_trig_eff.cxx` into the shared
`dr_correction_plateau.h` rather than copied (R17's lesson). It is reported, and applied to nothing.

**★ WHAT THE MERGE ACTUALLY BUYS — the delivered correction, not the χ².** Census from
`DrCorrectionEvaluator` (pp24 FULL, opposite sign, the series the analysis applies):

| WP / method | cells | fitted | raw-bin placeholder | **no correction at all** | delivered ε_ΔR span |
|---|---|---|---|---|---|
| Tight `expo` un-merged | 72 | 68 | 2 | **2** | 0.061 – 1.354 |
| Tight `expo` **MERGED** | 63 | 60 | 3 | **0** | 0.061 – 1.318 |
| Tight `polyu_fixedRp` un-merged | 72 | 65 | 4 | **3** | 0.078 – **2.736** |
| Tight `polyu_fixedRp` **MERGED** | 63 | 62 | 1 | **0** | 0.059 – 1.240 |
| Medium `expo` un-merged | 72 | 70 | 0 | **2** | 0.015 – 1.354 |
| Medium `expo` **MERGED** | 63 | 61 | 2 | **0** | 0.015 – 1.219 |
| Medium `polyu_fixedRp` un-merged | 72 | 66 | 3 | **3** | 0.065 – **2.749** |
| Medium `polyu_fixedRp` **MERGED** | 63 | 62 | 1 | **0** | 0.063 – 1.243 |

Two things follow, and they are the point of the variant:
1. **No dead cell is left.** In the un-merged mode 2–3 cells per configuration deliver **no
   correction at all** (ε_ΔR ≡ 1 by default) — all of them in the top two pair-pT bins. Merged,
   every one of the 63 cells delivers a measured correction.
2. **The ε_ΔR = 2.74 fit artefact is GONE.** The `polyu_fixedRp` cell with C = 0.053 ± 3.18 that
   inflated f/C to 2.74 — flagged as an artefact in `mc_trig_eff_closure.md` R2 and named in this
   doc's open "sanity bound on the delivered ε_ΔR" — does not survive the merge: the merged maximum
   is 1.24. A 2mu4 close-by correction is a LOSS, so a delivered value above ~1 was never a
   measurement; the merged set no longer produces one above 1.32.

**What it does NOT buy: fit quality.** χ²/ndf does not improve, and should not be expected to — with
twice the pairs the error bars shrink, so the same shape mismatch costs more χ². pp24 Tight,
opposite sign, `expo`: merged mean 2.03 / median 1.59 / max 8.80 over 63 cells, against un-merged
1.81 / 1.51 / 6.07 over 70. **R26 is neither fixed nor worsened by this change** — it is the same
parametric family on better-populated cells.

**Inclusive cell unchanged, as it must be** — it integrates over pair pT, so the merge cannot touch
it: χ²/ndf = 10.319 (Tight, opposite sign, `expo`), identical to the un-merged `nocorr` value.

**★ THE pp24 CROSSX DEFAULT (user, 2026-08-17; TEMPORARY).** `expo` + **opposite sign** + this merged
mode is the variant the pp24 cross-section will apply, named once in
`dr_correction_sample_cfg.h` as `DrCorrCrossxMethod()` / `DrCorrCrossxSign()` / `DrCorrCrossxMode()`
so no consumer retypes it. The applied form is unchanged from `mc_trig_eff_closure.md` §3.2:
`ε_ΔR(ΔR) = f(ΔR)/C` for ΔR < 1, `= 1` for ΔR ≥ 1. TEMPORARY because R26 is still OPEN.
**Scope was limited to producing the variant (user): the crossx chain was NOT touched here**, and
`DrCorrectionEvaluator::Load()` still DEFAULTS to the un-merged `nocorr`, so the MC-closure thread's
results do not move.

**Artefacts.** 576 PNGs + reports across four trees (pp / PbPb × Tight / Medium), 8 canvases per
(method, sign mode, ΔR view) = 7 pair-pT + inclusive; 18 pp + 18 overlay fit ROOT files
`dr_correction_fits_<label>{,_medium_wp}_step3_<method>{,_ss,_os}_nocorr_ptmerge.root`. Driver run
`SAMPLES="pp_full overlay" WPS="tight medium" STEPS="3" PLATEAU_MODES="nocorr_ptmerge"
SKIP_MEASURE=1 ./run_dr_correction_fits.sh` → exit 0, "all stages complete", zero artefact failures,
zero read-back persistence failures. Step 4 remains corrected-mode-only (user).

**★ TWO OPEN ITEMS THE ROUND-12 PLOT REVIEW RAISED — BOTH NEED A USER DECISION, BOTH ARE R26.**
Neither was changed autonomously: each would alter the published `fit_ok` artefact for all three
plateau modes and therefore the closure thread's inputs as well.
1. **The new default series contains a cell whose fit is accepted but unusable.** Merged top cell
   × η^pair ∈ [−2.4,−2.0), opposite sign, `expo`: **A = +0.5855 — a POSITIVE amplitude, i.e. a
   close-by ENHANCEMENT**, which a 2mu4 *product* trigger cannot physically be; C = 0.4900 ± 0.7162
   (a **146 %** relative error), `p` pinned on its limit, χ²/ndf = 8.80, f(0)/C = **2.19**, and
   `fit_ok = 1`. Identical at Medium. **It is NOT delivered as 2.19** — `DrCorrectionEvaluator`'s
   own baseline screen `DrCorrPlateauUsable(C, err)` rejects it on both criteria and falls back to
   the raw bins, which deliver **1.168** at ΔR → 0 (in the un-merged mode that cell delivered no
   correction at all). But two things remain true: the CANVAS draws a curve the analysis will never
   use (R28 defined ONE gate, `h_stepN_fit_ok`; `dr_correction_apply.h` adds a second screen the
   canvas does not know about), and 1.168 > 1 is still not a close-by loss. **A screen on `A > 0`
   (equivalently `f(0) > C`), or the χ² term R26 is about, would remove this class of cell.**
2. **The Pb+Pb overlay gains exactly one unconstrained fit.** In the merged top cell ×
   η^pair ∈ [−0.5,0.5) the overlay now has 7 points and is fitted with A = 0.586 ± 5.5,
   λ = 0.0371 ± 2.1, p = 3.42 ± 4.4, C = 0.661 ± 0.15 — **every parameter error ≥ its parameter** —
   and it is drawn with the same authority as the pp FULL-sample panels. In the un-merged tree all
   18 overlay top-p_T cells were "too few points → no fit". The `k < nfree + 2` threshold admits
   it. Raising that threshold is again a change to the acceptance rule.

**What the merge actually buys, measured (plot review).** The statistics gain is **modest**:
N_eff rises by only **1.06–1.48×** depending on η^pair (the `[104,150)` bin's filled ΔR bins are a
SUBSET of `[72.1,104)`'s — the filled-bin count is unchanged in all nine cells), so error bars
shrink **3–18 %**, not dramatically. The real gain is **coverage**, and it is confined to the
deliverable's top row (pp24 Tight, opposite sign, `expo`, ΔR → 0):

| η^pair bin | un-merged `[104,150)` | merged `[72.1,150)` |
|---|---|---|
| [−2.4,−2.0) | **no correction (ε_ΔR ≡ 1)** | raw placeholder, 1.168 |
| [−2.0,−1.5) … [0.5,1.0) | fitted, 0.123 – 0.382 | fitted, 0.116 – 0.417 |
| [1.0,1.5) | raw placeholder | raw placeholder |
| [1.5,2.0) | fitted, 0.400 | fitted, 0.399 |
| [2.0,2.4) | **no correction (ε_ΔR ≡ 1)** | **fitted, 0.304** |

One caveat the review raised and this doc records rather than dismisses: because the two source
bins can carry genuinely different ε_ΔR shapes, a merged cell describes a **mixture**. In the
SIGN-INTEGRATED series that shows up as two cells whose merged fit is *worse* than the un-merged
`[72.1,104)` one (η^pair ∈ [1.0,1.5) becomes rejected, f(0) = −0.95; η^pair ∈ [−1.5,−1.0) falls
into a cusp solution). **Neither is present in the opposite-sign deliverable** (checked cell by
cell, table above) — but the sign-integrated tree carries them and they belong to any future
decision about widening the merge.

**Remaining known artefact of the crossx deliverable (from the round-12 code review, INFO).** One
merged cell still delivers ε_ΔR above 1: `p_T^pair ∈ [11.5,16.6) × η^pair ∈ [−2.4,−2.0)`, Tight,
opposite sign, `expo`, where the fit returns **A = +0.2838** (a positive amplitude = a close-by
*enhancement*, which a 2mu4 product cannot physically be), **λ = 0.0208 pinned on its 0.02 lower
limit**, C = 0.8922, χ²/ndf = 4.65 ⇒ f(0)/C = **1.318**. It is **pre-existing and untouched by the merge**: the
un-merged `nocorr` tree has this identical cell at the identical 1.318 (that pair-pT bin is nowhere
near the merged one). The un-merged tree's own maximum, **1.354, sat in a different cell** —
`[72.1,104) × [−2.4,−2.0)` — and that cell is absorbed by the merge. The code's own alarm only fires
above 1.5, so both are published silently. It belongs to **OPEN R26**: a screen on
`A > 0` (equivalently `f(0) > C`) would remove this class of cell, and that is a fit-form decision,
not something to patch here. Consumers should treat any delivered value above 1 as an artefact.

**★ WHAT THE MERGE COSTS AND WHAT IT DOES NOT BUY (round-12 plot review, measured).** Two honest
qualifications, both measured from the source 3D histograms rather than assumed:
- **The statistics gain is modest.** Pooling `[104,150)` into `[72.1,104)` raises the effective
  entry count in the fit domain by only **1.06×–1.48×** depending on pair-η — error bars shrink
  **3–18 %**, not dramatically — and the number of FILLED ΔR bins is unchanged in all nine cells
  (13/17/20/20/20/20/19/17/12 before and after), because `[104,150)`'s filled bins are a subset of
  its neighbour's. The merge's real benefit is **coverage**, not precision (below).
- **In the sign-INTEGRATED series it costs one cell**: `η^pair ∈ [1.0,1.5)` had `fit_ok = 1` at
  `[72.1,104)` and the merged cell is rejected (f(0) = −0.95), so nothing above 72.1 GeV is fitted
  there. **This does NOT affect the deliverable**: in the opposite-sign series that cell was already
  on the raw-bin placeholder before the merge and still is.

**Per-cell status of the DELIVERABLE (pp24 Tight, opposite sign, `expo`), top pair-pT row** — the
row the merge acts on, at ΔR → 0:

| η^pair bin | un-merged `[104,150)` | merged `[72.1,150)` |
|---|---|---|
| [−2.4,−2.0) | **no correction (ε_ΔR ≡ 1)** | raw placeholder, **1.168** |
| [−2.0,−1.5) | fit 0.268 | fit 0.417 |
| [−1.5,−1.0) | fit 0.382 | fit 0.221 |
| [−1.0,−0.5) | fit 0.123 | fit 0.150 |
| [−0.5,0.5) | fit 0.198 | fit 0.211 |
| [0.5,1.0) | fit 0.167 | fit 0.116 |
| [1.0,1.5) | raw 0.188 | raw 0.436 |
| [1.5,2.0) | fit 0.400 | fit 0.399 |
| [2.0,2.4) | **no correction (ε_ΔR ≡ 1)** | fit 0.304 |

⇒ the two dead cells are gone and `[2.0,2.4)` gains a real fit; the price is that `[−2.4,−2.0)`
now delivers **1.168 from the raw-bin placeholder where it previously delivered nothing** — above 1,
so an artefact of that forward cell rather than a measurement (same forward region as the r16578
anomaly, R8/R10/R14).

**★ OPEN, ESCALATED TO THE USER (round-12 reviews) — R26 is now biting the DEFAULT series.**
1. **`fit_ok` still carries no χ²/shape term**, and in the merged top cell × `η^pair ∈ [−2.4,−2.0)`
   the `expo` fit returns **A = +0.5855 — a POSITIVE amplitude, i.e. a close-by *enhancement***,
   which a 2mu4 product cannot physically be — with **C = 0.4900 ± 0.7162 (146 % relative error)**
   and `p` pinned at its limit, χ²/ndf = 8.80, yet `fit_ok = 1` and the canvas draws it. It is only
   the CONSUMER's extra baseline screen (`DrCorrPlateauUsable(C)` in `dr_correction_apply.h`) that
   keeps f(0)/C = 2.19 out of the analysis. Two consequences: the canvas and the consumer disagree
   about which cells are usable (R28 defines ONE gate; `apply.h` adds a second the canvas does not
   know about), and the fallback still delivers 1.168 > 1. A screen on **A > 0** (equivalently
   f(0) > C) or on χ² inside `usable` would remove this class of cell — but it changes the published
   artefact for ALL three plateau modes and the MC-closure thread's inputs, so it is a **USER
   DECISION**, deliberately not taken here.
2. **One overlay panel gains an unconstrained fit.** In the Pb+Pb top cell the merge turns one
   "too few points" cell (`[72.1,150) × [−0.5,0.5)`, 7 points) into a fitted one whose **three
   SHAPE parameters are entirely unconstrained** — A = 0.586 ± 5.5, λ = 0.037 ± 2.1, p = 3.42 ± 4.4,
   each error larger than its own parameter — while only the baseline C = 0.661 ± 0.15 is measured
   (values from the sign-integrated series; opposite sign is numerically equivalent, C identical).
   It is drawn with the same authority as a pp FULL-sample panel. The overlay is a
   10 000-event TEST sample and is not a deliverable, but the acceptance rule
   (`k ≥ nfree + 2`) is what admits it — also a **USER DECISION**.

**Review round (both mandated reviewers run 2026-08-17/18).** Four defects were found and FIXED,
three of them PRE-EXISTING and inherited by the new tree rather than caused by it:
1. **The read-back report named the WRONG fit file in every non-nominal plateau mode**
   (`plot_dr_correction_fits.cxx` called `DrCorrFitFile(...)` without the mode, which defaults to
   the plateau-corrected name). So the persistence audit of the `nocorr` tree — and now of the
   *deliverable* merged tree — could not be traced to the artefact it audited. Fixed; the
   `no_plateau_correction/` and merged trees' reports were regenerated and verified on disk to name
   their own file. `plateau_corrected/` needed none — for that mode the fixed call produces exactly
   the name the old defaulted call already produced.
2. **Every pp24 FULL report described the sample as "a 10k-event TEST sample"** in the line
   carrying the headline inclusive χ²/ndf. The parenthetical is now conditional on
   `is_full_sample`: pp24 reports drop it, the genuinely 10 000-event overlay keeps it.
3. The merged-mode plateau re-measurement fetched the Step-4 covariance terms with a tolerant
   `Get` while the four beside them threw; now uniformly `GetObj<TH3D>` (unreachable in step 3,
   where `has_cov` is false, but it would have silently under-stated merged Step-4 plateau errors).
4. `FillMCTrigEffClosure.cxx` duplicated the `"nocorr"` token as a literal in its provenance stamp
   while `Load()` had just gained a defaulted mode argument; it now carries an explicit
   `kDrPlateauMode` constant passed to both, so the stamp cannot drift from what was loaded.
Only the MERGED tree was refitted for (2). The `plateau_corrected/` and `no_plateau_correction/`
fit files were deliberately **NOT** regenerated: a sibling session's MC-closure results consume
them and rewriting their mtimes would make that session's freshness gate call its own outputs
stale. **Those two trees therefore still carry the old report text**; the code is fixed and they
will pick it up at their next legitimate regeneration.

**Regression: the shared-header refactor moved NOTHING.** `PlateauFromRatio` was relocated and
`DrCellRatio` became a thin wrapper over the new explicit-range form, both of which are on the
NOMINAL path. Re-running pp24 Tight `expo` in the pre-existing `corr` (sign-integrated) and `nocorr`
(opposite-sign) modes after the change reproduced their `fit_report*.txt` **byte-identically**
(`diff` clean, both).

### D6: the gap cut lives at the RDF stage, not in the ntuple processing (2026-08-05)
Round-8 Contract item 4 specified "a gap-cut mode in the ntuple-processing code". It was
implemented at the RDF stage instead. **Provenance rule satisfied** (the cut is applied to
ntuple-processing OUTPUT branches — `charge`, `eta` — never re-derived from raw NTUPs), and three
independent reasons favour it: every other muon-level fiducial cut in this chain already lives
there; the DATA gap cut cannot go in the data NTP (out of scope), so an NTP-side MC cut would make
MC and data asymmetric; and it avoids re-running 8 grid-scale NTP jobs. Recorded here because the
reviewer correctly flagged the deviation as undocumented.

## Results & Observations

### R34. Round 14 — the Step-3 POLYNOMIAL fit constrained at ΔR = 0 (2026-09-08, pp24 only)

**User request:** *"For the step 3 polynomial fit, require that when dR = 0, the efficiency cannot
have value exceeding plateau."* This is requirement (1) of R33(b) — the small-ΔR value must lie
below the large-ΔR plateau — now imposed on the BACKUP form `polyu_fixedRp` as well as on the
nominal `expo`. Scope mirrors round 13: **pp_full only, Step 3 only**; the HIJING overlay,
`noovl`, PbPb, the Step-4 fits and the `_pt4bin` variant stay knowingly STALE.

**(a) Why it needed a reparametrization, not a limit.** `f = C + u²(a₂ + a₃u + a₄u²)` with
`u ≡ max(0, 1 − ΔR/R_p)` and `u(0) = 1`, so `f(0) = C + (a₂+a₃+a₄)`. The requirement `f(0) ≤ C`
is the LINEAR constraint `a₂ + a₃ + a₄ ≤ 0` — a bound on a COMBINATION of three parameters, which
MINUIT cannot express as a `SetParLimits` on any one of them. (This is the one structural
difference from `expo`, where the same requirement was already a bound on the single parameter
`A`.) The fit is therefore written with the constrained combination AS parameter `p0`:

    A ≡ f(0) − C = a₂ + a₃ + a₄ ,   a₂ = A − a₃ − a₄
    f = C + u² [ A + a₃(u − 1) + a₄(u² − 1) ]

The function FAMILY is untouched — expanding gives back `C + a₂u² + a₃u³ + a₄u⁴`, and `a₃`/`a₄`
still ARE the `u³`/`u⁴` coefficients — so the fit's flexibility is unchanged; only which of its
directions can be bounded. The requirement is then the single limit **`A ∈ [−50, 0]`** (was
`[−50, 50]`). Same symbol, same meaning as `expo`'s `A`. Verified as an exact identity before any
production run: over 15 (C, a₂, a₃, a₄) sets and ΔR ∈ [0, 2], max |old − new| = **2.2e-16**,
`f(0) − C = A` to 1e-12, and `f(R_p+) = C` exactly.

**(b) What is deliberately NOT constrained.** Only the ΔR = 0 value, which is what was asked. The
polynomial may still be non-monotonic and may still rise above `C` BETWEEN 0 and `R_p` — that
flexibility is the whole reason this backup form exists next to the monotone `expo` (which, by
`A ≤ 0`, cannot exceed `C` anywhere). No turning-point / zero-slope condition is imposed here:
R33(b)'s condition (2) is specific to the exponential's `p`, and the quartic has no analogous
single parameter.

**Step 4 is EXCLUDED**, by the same physics as `expo` (R4/§3.4: the single-leg correction is an
ENHANCEMENT at small ΔR, `A > 0`). `restrict_shape = (step == 3)` is now defined ONCE and shared
by both restricted methods, so they cannot drift apart on which step is restricted.

**(c) What was rerun.** `run_dr_correction_fits.sh` with `SAMPLES=pp_full STEPS=3 SKIP_MEASURE=1
METHODS=polyu_fixedRp` — 2 WPs × 5 plateau modes × 3 sign series = **30 fit + plot
configurations** (the plateaus are an unchanged input, hence `SKIP_MEASURE`). Exit code 2 is the
PRE-EXISTING plateau-guard failure of the `plateau_corrected` mode, unchanged by this work.
Then the MC closure — see (e).

**(d) The restriction holds.** Over all 30 regenerated reports (the FINAL set, 2026-09-08
23:52–23:57): **1482 fitted cells, ZERO cells with `A > 0`, and `f(0) ≤ C` in 1482 / 1482.** `A`
reproduces `f(0) − C` to print precision (1e-4, the 4-dp `C=` column; ~1e-16 in
`plateau_corrected` where `C = 1` exactly), and the driver's flatness audit reports worst
`|f − C|` beyond `R_p` = **0.00e+00** in all 30 blocks.

**⚠ THE INPUT MOVED BETWEEN THIS ROUND'S TWO RUNS, AND THE TIGHT FIT SET IS NOW INTERNALLY
INCONSISTENT ACROSS METHODS.** The CONCURRENT workstream of (h) **refilled the Tight Step-3
histogram at 2026-09-08 18:48 and its plateau map at 18:51** — on its NEW 9 GeV pair-pT axis —
between this round's first run (16:48) and its second (23:52, done only to pick up the corrected
report note of (g)). Read off the artefacts, first pair-pT cell of the Tight `nocorr` reports:

| Tight, Step 3 | `expo` | `interp` | `polyu_fixedRp` |
|---|---|---|---|
| first pair-pT cell | `[8.0,11.5)` | `[8.0,11.5)` | **`[9.0,12.8)`** |
| written | 09-08 00:53 | 09-08 00:56 | 09-08 23:53 |

So **the Tight polyu fits are on the 9 GeV axis while the Tight `expo`/`interp` fits are still on
8 GeV**, i.e. the three tiers of the delivered `expo → polyu → interp` cascade no longer describe
the same cells for Tight. **Medium is self-consistent** (all three methods at 8 GeV; its Step-3
input is still 2026-09-07 21:21), so the Tight and Medium halves of the table in (e) also come
from different fills. Nothing is published from this state — `DrCorrectionCascadeEvaluator`'s
canonical-edge check and `combine_dr_correction_fits`'s point-by-point raw cross-check both THROW
— and the fix is the one (h) already names: refill Step 3 on the 9 GeV axis and refit ALL THREE
methods together. **The restriction itself held identically in both runs** (0 violations,
`f(0) ≤ C` everywhere, on either axis); only the cell count (1483 → 1482) and the rail counts
moved. Every number quoted here is read off the artefacts now on disk.

**(e) How often the limit actually binds — and why the `AT LIMIT` flag over-counts here by ~2×.**
`ParAtLimit` flags a parameter within `1e-3 × (hi − lo)` of a limit. For polyu that range is
`[−50, 0]`, so the tolerance is **0.05 — ten times `expo`'s** (range `[−5,0]` → 0.005) and
comparable to a typical fitted `|A|` of 0.1–0.8. **171 cells are FLAGGED but only 95 are actually
railed** (`|A| < 1e-4`); the rest are ordinary fits sitting just inside the tolerance. The flag
count is therefore NOT comparable with `expo`'s without correcting for that factor, and the fit
report now says so explicitly. Both counts, per series:

| WP | plateau mode | fitted (intgr/os/ss) | flagged `AT LIMIT: A` | railed \|A\|<1e-4 |
|---|---|---|---|---|
| T | plateau_corrected | 71 / 71 / 65 | 10 / 11 / 7 | 8 / 9 / 4 |
| T | nocorr | 71 / 71 / 65 | 7 / 8 / 3 | 4 / 4 / 2 |
| T | nocorr_ptmerge | 64 / 64 / 62 | 7 / 7 / 2 | 4 / 3 / 1 |
| T | nocorr_etamerge | 24 / 24 / 23 | 3 / 3 / 1 | 3 / 3 / 0 |
| T | nocorr_etamerge_ptmerge | 22 / 22 / 21 | 3 / 3 / 1 | 3 / 3 / 0 |
| M | plateau_corrected | 71 / 71 / 65 | 10 / 10 / 4 | 6 / 6 / 3 |
| M | nocorr | 71 / 71 / 65 | 10 / 10 / 6 | 4 / 3 / 4 |
| M | nocorr_ptmerge | 64 / 64 / 62 | 9 / 9 / 5 | 3 / 2 / 3 |
| M | nocorr_etamerge | 24 / 24 / 24 | 4 / 4 / 3 | 2 / 2 / 1 |
| M | nocorr_etamerge_ptmerge | 22 / 22 / 22 | 4 / 4 / 3 | 2 / 2 / 1 |

**95 of 1482 (6.4 %) genuinely rail** — the same order as `expo` once the tolerance difference is
accounted for, and the expected consequence of closing a half-space on a flexible quartic in
low-statistics forward / high-pT cells. A railed cell delivers `ε_ΔR(0) = f(0)/C = 1`, i.e. no
correction at ΔR = 0, and must be read as CONSTRAINED, not measured.

**(f) TWO MEASURED SIDE-EFFECTS of the endpoint-only form of the requirement. Both were put to
the user, who confirmed ΔR = 0 only, as asked** — so they are recorded here as known properties of
the delivered fit, not as defects.

1. *The `>plateau` pathology is partly DISPLACED, not removed.* Only `f(0)` is pinned; the quartic
   may still rise above `C` between 0 and `R_p`, and it does. Max `f/C` on `(0, R_p]` exceeds 1 by
   >1 % in **701 cells (47 %)**, >5 % in **278 (19 %)**, >20 % in **90 (6.1 %)**, >100 % in **32**;
   worst **11.47** (Medium, `nocorr`, pT [50.0,72.1) × η [−2.2,−2.0), same sign). And the railed
   cells do it MORE often than the rest — **25 % vs 18 % above 1.05** — i.e. the constraint pushes
   some of the overshoot from ΔR = 0 to ΔR ≈ 0.1–0.3. This is exactly the flexibility that makes
   `polyu` the backup form next to the monotone `expo` (which, by `A ≤ 0`, cannot exceed `C`
   anywhere); the stronger requirement `ε_ΔR ≤ 1` over the whole domain is NOT expressible as a box
   limit (it is a bound on `max_u`) and would need a penalty term or a post-fit `usable` screen.
   **Not implemented — the user chose the ΔR = 0 requirement as stated.**
2. *In the `nocorr*` family the fit has a second way to satisfy the constraint: collapse `C`.*
   `C` is free there, so `f(0) ≤ C` can be met by lowering the baseline instead of the endpoint.
   **11 of the 95 railed cells** came out with `C < 0.5` or `σ_C > 0.3·C`; worst
   `C = 0.095 ± 0.110` against a MEASURED plateau of 1.134 (same cell as above), with
   `a₃ = −36.7, a₄ = +19.2`. **9 of the 11** are caught downstream by
   `DrCorrPlateauUsable` (`dr_correction_sample_cfg.h`: `C ≥ 0.5` and `σ_C < C`, the same screen a
   bad measured plateau gets) and fall through to the `interp` tier, so the cascade ROUTING moves.
   **The other 2 CLEAR the screen and ARE delivered at the polyu tier**: `C = 0.7805 ± 0.5334`
   (Tight, `nocorr`, same sign, pT [105.5,150.0) × η [0.5,1.0)) and `C = 0.7776 ± 0.5335` (Medium,
   the same cell family), both with `σ_C/C ≈ 0.68` and max `f/C` ≈ 1.90 / 1.94 — exactly the case
   `dr_correction_apply.h`'s own comment already warns about ("a C that only just clears the
   screen … can still put f/C above 2"). So it is NOT true that nothing wrong can be published;
   what is true is that both survivors are same-sign, top-pair-pT cells. The producer-side remedy
   (adding a `C` sanity term to `fit_ok` itself) was offered and NOT chosen.

**(g) Two review fixes carried in the same change** (`/review-analysis-code`, iteration 1 FAIL →
both addressed, reports regenerated):
- the polyu report note now carries the `ParAtLimit`-tolerance caveat of (e), which the `expo`
  note already had — without it the note asserted that every flagged cell is a constrained value,
  false for half of them;
- `plot_dr_correction_fits.cxx::LoadFunc` now THROWS if a `polyu_fixedRp` TF1's parameter 0 is not
  named `A`. `MethodFormulaTex()` is keyed on the METHOD NAME alone, and the reparametrization
  left the parameter COUNT unchanged and the persisted TF1 still `Eval`-correct — so replotting an
  older fit file (the `_pt4bin` ones on disk, 2026-08-07/13, are exactly this) would have drawn
  the NEW equation over OLD `(a₂, a₃, a₄)` parameters and nothing would have noticed. The
  parameter NAME is the discriminator and it survives write/read.

**(h) NOT RERUN — the MC closure, and why (user decision).** Plan item 6 could not run: a
CONCURRENT session in this same checkout has (uncommitted, ~16:19 on 2026-09-08) adopted the muon
`pT > 4.5 GeV` cut and moved the canonical pair-pT axes **8 → 9 GeV**
(`ParamsSet::signal_pair_pt_min`, `pair_pt_coarse_bins`, `pT_bins_150`, …). Every Step-3 histogram
and fit file on disk is still on the 8 GeV axis, so `DrCorrectionCascadeEvaluator` correctly throws
*"pair-pT edge 0 is 8.000000 in the fit file but 9.000000 canonically — stale fit file"*. **This
staleness predates and is independent of this change** — the round-13 fit files fail identically.
Refilling Step-3 on the 9 GeV axis is that other workstream's blast radius and would write to the
same files and output paths it is actively producing, so the user directed this round to **stop at
the fits**. The four-approach closure figures of `mc_trigeff_dr_binning_approaches.md` (2026-09-08
01:33) are therefore **STALE with respect to the polyu tier of the cascade** and must be
regenerated by whoever refills Step-3 on the new binning — recorded there too.

*Each Step-3 fit artefact is consistent with the histogram it was fitted to — but as (d) records,
those histograms are no longer the same across methods for Tight (polyu on the 9 GeV refill of
18:48; `expo`/`interp` on the 8 GeV fill of 00:49–00:59). `dr_correction_cell_groups.h` reads only
`ParamsSet::N_COARSE_PAIR_PT_BINS`, still 8, so the CELL COUNT is unaffected by the concurrent
edits — only the EDGES moved, and only for the Tight polyu set.*

### R33. Round 13 — the new gap-cut set + the pair-level |η^pair| < 2.2, and the Step-3 `expo`
### shape restriction (2026-09-07/08, pp24 only)

**Scope, confirmed with the user before any work (AskUserQuestion, 2026-09-07): pp_full only on
the MC side, pp24 only on the data side.** The HIJING overlay, `noovl`, PbPb data, the Step-4
fits and the MC closure are OUT of scope and are knowingly left STALE. ε_ΔR was NOT propagated to
the pp24 cross-section (user: future work, other agents in flight).

**(a) What moved in the input.** `muon_gap_cuts_acceptance.md` F17 changed the cut set and added a
pair-level window; that session committed the code (`f7c5abe`…`8da992e`) but ran nothing. This
round is the rerun. Live in every number below:
`single_mu_fiducial_gap_cuts` = **{{−1.30,−1.05}, {−0.10,+0.06}, {2.20,2.40}}**; NEW
`ParamsSet::pair_eta_fiducial_max = 2.2` on `kGapLeg` (Steps 2/4) and `kGapPair` (Step 3) but NOT
`kGapSingle` (Step 1 fills a single-muon map — there is no pair to cut on); coarse q·η top bin
{2.0,2.3} → **{2.0,2.2}** (forced by the `SetIOPaths` startup throw); coarse pair-η outer bins
±(2.0,2.4) → **±(2.0,2.2)**.

*A concurrency check that mattered.* The sibling session edited `FillMCTrigEffHists.cxx` at
21:27:49, i.e. AFTER this round's Step-1–4 fills (21:17–21:26) and while the plot stage ran. The
edit refactored `kGapPair` from an inline string to `MCTrigEffPairSel::FiducialGapCut()`. Both
forms were compiled standalone and print the **byte-identical** selection string, so the filled
histograms are valid against HEAD. Nothing else the sibling touched feeds this chain.

**(b) The Step-3 `expo` shape restriction (user).** `f(ΔR) = C + A·exp[−(ΔR/λ)^p]`. The two
requirements map onto the parameters exactly, and **they are not nested — they are the same
constraint**:
1. *small-ΔR plateau below the large-ΔR plateau* ⟺ `f(0) = C+A < C` ⟺ **`A < 0`**, which also makes
   `f` monotone increasing and so removes the "decreasing with ΔR" failure at the same time.
2. *turning point at ΔR > 0*, with *"slope at ΔR = 0 should be zero"* to be tried first. For this
   form `f'' = 0` at `ΔR_infl = λ((p−1)/p)^{1/p}`, real and positive **iff `p > 1`**; and
   `f'(0) = 0` **iff `p > 1`** as well (`p = 1` → finite slope `|A|/λ`; `p < 1` → infinite slope at
   0 and no positive inflection, which IS the "concave rise, turning point < 0" failure). So the
   stronger requirement and its stated fallback collapse to the single condition `p > 1`, and
   **there is no weaker variant to relax to** if too many cells are rejected. Derived analytically
   and confirmed numerically at p = 0.5/0.8/1.0/1.001/1.5/3.0.
Both are STRICT inequalities that MINUIT cannot express, so the imposed limits are their closures
`A ∈ [−5,0]`, `p ∈ [1,8]`. A cell railed at `A = 0` or `p = 1` is the BOUNDARY case — a constrained
value, not a measurement — and both the fit report and the ROOT `provenance` now say so, with both
ends of each rail spelled out (`p = 8` and `A = −5` are the pre-existing limits and mean the
opposite things).
Imposed as parameter limits, not as a post-hoc `fit_ok` screen: a screened-out cell delivers NO
correction (it falls through to the backup method / raw-bin placeholder), whereas a constrained fit
delivers the best physically admissible fit — and "restrict the fitting parameters" is what was
asked.
**Step 4 is deliberately EXCLUDED** (`restrict_shape = (step == 3)`): its single-leg correction is
physically an ENHANCEMENT at small ΔR (R4/§3.4, ε(ΔR<0.12)/ε(ΔR>1) ≈ 1.20 pp), i.e. `A > 0`.

**(c) What was rerun.** Data pp24 T&P Tight (`pipeline_pp_trig_eff.sh` stages 5–8) + Medium
(`run_data_trigeff_medium_wp.sh`); MC `pp_full` Steps 1–4 + turn-on refits + sanity + plots + 2D
maps + statistics tables (`run_mc_trigeff_round7.sh`); Step-3 ΔR fits, all 3 methods × 3 sign
series × 5 plateau modes, both WPs (`run_dr_correction_fits.sh STEPS=3`). Data top q·η bin is now
`2_00_TO_2_20` in the turn-on fit file.

**(d) Headline numbers.**
- Step-1 MC single-muon ε(mu4), pp24 FULL Tight, on the new cuts: **0.8025 (μ⁺) / 0.7961 (μ⁻)**
  (was 0.7722 summed on the old cuts — the rise is expected, the widened crack and barrel/endcap
  windows remove low-efficiency muons).
- Step-1 MC/data ratio vs pT (μ⁺): 1.206 (4.3 GeV), 1.115, 1.103, 1.114, 1.129, 1.081 (40 GeV) —
  consistent with the known ~1.13 per-leg L1 MC over-efficiency (§3.1).
- **Every pair-η cell edge is now what the user asked for**: the 9-bin views run
  `[−2.2,−2.0) … [2.0,2.2)` and the folded views `[0,1), [1,2), [2,2.2)`. Verified on the produced
  artefacts, not assumed: **zero occurrences of 2.3 or 2.4 in any pair-η context** anywhere in the
  source, the reports or the canvases. (A `−2.4` on a SINGLE-MUON q·η axis is correct and stays —
  the forward window is one-sided.) Nothing is retyped: `MakeDrEtaGroups` reads the outer edge off
  the filled axis.
- **The restriction holds identically**: over all 30 Step-3 `expo` reports (both WPs × 5 plateau
  modes × 3 sign series), 1542 rows = 1453 carrying fitted parameters (1452 converged + 1
  non-converged) + 30 `inclusive` + 59 with no fit (46 too-few-points + 13 no-plateau).
  **Zero violations of `A ≤ 0` and zero of `p ≥ 1`.** A spans [−5.0000, −0.0000], p spans
  [1.0000, 8.0000].
- **How much the restriction actually bound** — cells railed on the NEW limits vs the pre-existing
  ones (`expo`, per-cell rows; `ptmerge` = last two pair-pT bins merged, `etamerge` = folded |η|):

  | WP | plateau mode | series | fitted | p@1 NEW | p@8 old | A@0 NEW | A@−5 old | fit_ok=0 |
  |---|---|---|---|---|---|---|---|---|
  | T | plateau_corrected | intgr | 70 | 8 | 8 | 0 | 0 | 13 |
  | T | plateau_corrected | os | 70 | 10 | 7 | 0 | 0 | 16 |
  | T | plateau_corrected | ss | 64 | 11 | 12 | 1 | 0 | 25 |
  | T | nocorr | intgr | 70 | 4 | 8 | 0 | 0 | 4 |
  | T | nocorr | os | 70 | 6 | 10 | 0 | 0 | 4 |
  | T | nocorr | ss | 64 | 9 | 14 | 2 | 0 | 19 |
  | T | nocorr_ptmerge | intgr | 63 | 5 | 9 | 0 | 0 | 2 |
  | T | **nocorr_ptmerge** | **os** | **63** | **7** | **10** | **0** | **0** | **1** |
  | T | nocorr_ptmerge | ss | 61 | 8 | 12 | 1 | 1 | 12 |
  | T | nocorr_etamerge | intgr | 23 | 2 | 3 | 0 | 0 | 1 |
  | T | nocorr_etamerge | os | 23 | 1 | 3 | 0 | 0 | 1 |
  | T | nocorr_etamerge | ss | 22 | 1 | 3 | 1 | 0 | 3 |
  | T | nocorr_etamerge_ptmerge | intgr | 21 | 2 | 3 | 0 | 0 | 0 |
  | T | nocorr_etamerge_ptmerge | os | 21 | 1 | 3 | 0 | 0 | 0 |
  | T | nocorr_etamerge_ptmerge | ss | 21 | 1 | 2 | 1 | 0 | 1 |
  | M | plateau_corrected | intgr | 70 | 9 | 9 | 0 | 0 | 13 |
  | M | plateau_corrected | os | 70 | 11 | 8 | 0 | 0 | 15 |
  | M | plateau_corrected | ss | 64 | 11 | 12 | 2 | 0 | 25 |
  | M | nocorr | intgr | 70 | 5 | 8 | 0 | 0 | 4 |
  | M | nocorr | os | 70 | 5 | 10 | 0 | 0 | 4 |
  | M | nocorr | ss | 64 | 8 | 14 | 3 | 0 | 19 |
  | M | nocorr_ptmerge | intgr | 63 | 5 | 8 | 0 | 0 | 2 |
  | M | nocorr_ptmerge | os | 63 | 6 | 10 | 0 | 0 | 1 |
  | M | nocorr_ptmerge | ss | 61 | 8 | 13 | 2 | 1 | 11 |
  | M | nocorr_etamerge | intgr | 23 | 2 | 3 | 0 | 0 | 1 |
  | M | nocorr_etamerge | os | 23 | 1 | 3 | 0 | 0 | 1 |
  | M | nocorr_etamerge | ss | 23 | 1 | 4 | 1 | 0 | 2 |
  | M | nocorr_etamerge_ptmerge | intgr | 21 | 2 | 3 | 0 | 0 | 0 |
  | M | nocorr_etamerge_ptmerge | os | 21 | 1 | 3 | 0 | 0 | 0 |
  | M | nocorr_etamerge_ptmerge | ss | 21 | 0 | 2 | 1 | 0 | 1 |

  Read this as the measurement of how often the unconstrained fit wanted a flipped shape: `p@1` is
  a cell whose free fit wanted `p ≤ 1` (the concave-rise pathology), `A@0` a cell whose free fit
  wanted a small-ΔR ENHANCEMENT. **Same-sign is by far the worst** (9–11 of ~64 at `p = 1`, and
  every `A@0` cell but none of the sign-integrated or opposite-sign ones), which is consistent with
  R25: the same-sign small-ΔR correlation is the real one and the `expo` form struggles with it.
  The pre-existing `p = 8` rail (the fit wanting a step function, R26) is untouched and comparably
  common — it is a different pathology and the report now tells the reader to separate them.
- **★ The restriction closes R32 item 0 for fitted cells.** With `A ≤ 0`, `ε_ΔR = f(ΔR)/C ≤ 1`
  identically, so the "ε_ΔR = 1.168 / 2.19 at ΔR → 0" artefact — a close-by 2mu4 correction that
  came out as a GAIN — is now **structurally impossible for a fitted cell**. Measured on the
  crossx-consumed series (Tight, `expo`, opposite sign, `nocorr_ptmerge`): `f(0)` spans
  **−0.4626 … 0.9098** over 63 cells with **no cell above 1**, and the delivered
  `ε_ΔR(0) = f(0)/C` spans **0.113 … 0.944** over the 62 usable cells. Only **1** cell is rejected
  (was 3 before). It is NOT closed for the **raw-bin placeholder** path in `dr_correction_apply.h`,
  which is unconstrained and can still deliver > 1 — that remains open.
- ε_ΔR(0) deepens monotonically with pair pT: ≈0.86–0.92 in `p_T^pair ∈ [8,11.5)` GeV → ≈0.14–0.53
  in `[72.1,150)`, which is the expected close-by-L1 behaviour.

**(e) Reviews.** `/review-analysis-code` **APPROVED at iteration 4**, `/review-plot` **APPROVED at
iteration 2**. Between them they found and got fixed: a false statement about the `p = 1` boundary
published in all 30 fit reports (at `p = 1` neither requirement actually holds — the limit is the
CLOSURE of `p > 1`, and the report now says so); a rail description that documented only one end of
the `A` limit; a **§Binnings violation** — a NEWLY retyped `2.0 ≤ |η| < 2.4` in a header comment
while the produced axis was `[2, 2.2)`, which I had also *reported as fixed when it was not*; a
missing shape-restriction stamp in the ROOT `provenance` (the file the crossx evaluator actually
consumes); a clipped Step-1 data legend losing its closing bracket in 6 PNGs; a missing ε_ΔR
defining equation on 8 canvases; and **a refuted explanation left standing in this doc** (see the
corrected `raw_joint_trigger_probability` entry above — the OS-only resonance veto cannot be the
cause, because this sample is the `_no_data_resonance_cuts_` tree and there is no `minv` cut
anywhere in the MC trig-eff selection). No fitted number moved across any of the amendments:
90/90 fit reports were byte-identical in their data rows after each re-emit.

**(f) STALE, knowingly — do not read these as current.**
1. **Step-4 fits.** The Step-4 *histograms* were refilled on the new cuts (21:15–21:26) but the
   Step-4 *fits* (`step4_dr_fit/`, 2026-08-13) were NOT re-run — user: "ignore step4 fit for now".
   The Step-4 subtree is therefore internally MIXED: `step4_dr_correction_singles/` (measurement
   plots) is fresh, `step4_dr_fit/` is not.
2. **The `_pt4bin` variant** (R24b) — still stale, and now also carries UNRESTRICTED Step-3 `expo`
   fits with no shape stamp. Distinct `FileSuffix()`, no nominal consumer reads it.
3. **The MC closure** (`mc_trig_eff_closure.md`) — its inputs all moved; it must be regenerated.
4. **The HIJING overlay, `noovl`, and PbPb data T&P** — out of scope this round.
5. `ε_ΔR` is **not** propagated to the pp24 cross-section (user's explicit instruction).

**(g) Carried open, unchanged.** R23 (top-pair-pT plateaus 20–42 % high; 12 guard failures in the
Tight sign-integrated `corr` mode this round), R25 (the correction is not charge-blind), R26 (both
parametric forms fail in a minority of cells; the `p = 8` rail is its signature), R27 (the two
unquoted correlated uncertainties, and the plateau-error `k = n` binomial fallback), and the
`ParAtLimit` bare-bool limitation (the summary aggregate merges the `p = 1` and `p = 8` rails; the
per-cell `p` column and the new report text separate them).

### R24. Round 9 — sign-separated corrections, per-cell ΔR distributions, 2D singles map,
### and the statistics of the 8-bin pair-pT axis (2026-08-10/11, pp only)

**Scope note.** pp (`pp_full`) only; the HIJING overlay was updated for CODE consistency and its
PLOT stage re-run (same round-8 ROOT inputs → plateaus verified bit-identical, 288/288 cells), but
it was **not refilled**, so it has no per-sign histograms. `noovl` (r17663) was not touched at all,
per the user's standing instruction that it is a one-time Step-1 cross-check only.

**(a) Step-1: MC single-muon mu4 efficiency as a 2D (q·η, pT) map.** New macro
`plotting_codes/trig_effcy/mc_based/plot_mc_singles_2d_effcy.cxx`. Two canvases per WP — μ⁺/μ⁻ side
by side, and the two charges combined (numerators and denominators **added before dividing**, not
an average of two efficiencies). The histograms already existed (`h_mc_pt_vs_q_eta_{num,denom}_*`,
booked since round 7); no refill was needed. Binning verified **bit-identical** to the data
tag-and-probe 2D map — 184 q·η × 41 pT bins, max |edge difference| = 0 on both axes — so the MC and
data maps can be read against each other directly. Plateau (pT > 8 GeV, all q·η), Tight:
μ⁺ 0.9138, μ⁻ 0.9155, combined 0.9147; Medium 0.9115 / 0.9135 / 0.9125.
- *Confirmed not a bug:* the pT axis carries a **duplicate 8.0 GeV edge** (`pT_bins_8` ⊕ `pT_bins_60`
  concatenated without dropping the shared edge) ⇒ one permanently empty zero-width bin. It is
  inherited deliberately from the data construction (`RDFBasedHistFillingData.cxx`), and both sides
  have it, so the two maps stay aligned. 1144 of 7544 cells are undefined = 24 fully-empty q·η
  columns (the three gap-cut windows, measured at `[−1.200,−1.060)`, `[−0.060,0.060)`,
  `[2.300,2.400)` — exactly `single_mu_fiducial_gap_cuts` with the new 2.30 forward edge) plus the
  zero-width pT row: 40·24 + 184 = 1144 exactly. Structural, not statistical.

**(b) Step-3 ΔR-correction distributions, per cell, with the plateau drawn.**
`step3_dr_correction/` now has three subdirectories — `dr_0_to_1/`, `dr_0_to_2/`,
`dr_full_range/` — each holding **one PNG per pair-pT bin with one subplot per pair-η bin**, the
per-cell plateau drawn as a horizontal line (solid over ΔR ∈ [2, 3.5] where the window is on the
canvas, dashed across the pad otherwise). Three ranges because no single x range shows both the
small-ΔR rise (below ΔR ≈ 0.3) and the plateau window legibly. The 0–2 view concatenates the fine
0.05-wide bins below 1 with the wide bins above, which is the SAME point set the fit canvases draw
— the two figure sets are directly comparable, and no third ΔR binning was invented.
- **The two pair-pT/pair-η INTEGRATED canvases were dropped** (user: not useful — the correction is
  applied per cell, never inclusively). The inclusive plateau is still MEASURED and still written to
  the plateau ROOT file, which the fit stage and the guard both read; only the two figures are gone.
  Per the user's decision the pair-η-integrated pair-pT-slice canvases and the 9-panel overlays STAY.
- Rendering fixes found only by looking at the rendered PNGs: the super-title ran into the legend
  (the ε symbol was repeated in the title although every y axis carries it — removed, legend moved
  right), and the y range was being set from value ± error, which let a handful of wide-ΔR bins with
  errors of order 1 stretch every panel to the 3.0 cap and flatten the structure. Central values
  only now, as the fit canvases already did.
- Cells with no measurable plateau carry no plateau line and are counted in the log: pp 1/72 (Tight)
  and 0/72 (Medium); overlay 36/72 and 38/72 — the expected cost of a 10 000-event TEST sample on 72
  cells, reported rather than papered over.

**(c) Sign-separated Step-3 and Step-4 (`FillMCTrigEffHists.cxx`, `plot_mc_trig_eff.cxx`).**
Every Step-3/Step-4 histogram is now booked **twice**: sign-integrated (unchanged, still nominal)
and with a per-sign prefix `h_mc_dr_{ss,os}_` / `h_mc_single_dr_{ss,os}_`, one pair tree each
(`muon_pair_tree_sign1` = same sign, `sign2` = opposite sign — the split is on the TRUTH charges).
Closure checked: `ss + os = integrated` exactly (denominator 6.983715 + 30.225879 = 37.209594 nb).
Each sign gets its **own** plateau, written into the same plateau ROOT file under sign-tagged keys
(`h_step3_ss_plateau`, `prov_step3_os`, …) — normalizing a same-sign curve by an
opposite-sign-dominated plateau would import the very charge dependence the split exists to test.
Samples filled before this round degrade to a printed note instead of throwing (verified on the
overlay). **Inclusive plateaus, pp_full:**

| | Tight same sign | Tight opposite sign | Medium same sign | Medium opposite sign |
|---|---|---|---|---|
| Step 3 (ε_ΔR^2mu4) | 0.9889 ± 0.0016 | 0.9933 ± 0.0010 | 0.9900 ± 0.0015 | 0.9935 ± 0.0010 |
| Step 4 (ε_ΔR^single) | 0.9950 ± 0.0007 | 0.9967 ± 0.0005 | 0.9956 ± 0.0007 | 0.9970 ± 0.0005 |

The **plateaus** agree to ≈0.4 % (Step 3) and ≈0.2 % (Step 4). **That is NOT the same as
charge-blindness, and an earlier version of this entry wrongly concluded that it was.** §3.3 says
explicitly that "the physical content is the small-ΔR shape relative to the plateau" — and that is
exactly where the two signs part company. See R25.

### R25. ★ OPEN — the Step-3 (both-legs) ΔR correction is NOT charge-blind at small ΔR
### (2026-08-11, pp_full Tight; refutes the §3.4 "blind to the pair charge product" rationale for Step 3)

Plateau-normalized inclusive measurements, read directly from the Step-3 histograms (no fit):

| ΔR | same sign | opposite sign | same/opposite |
|---|---|---|---|
| 0.025 | 0.3570 | 0.8356 | **0.427** |
| 0.075 | 0.4794 | 0.8020 | **0.598** |
| 0.125 | 0.7186 | 0.8477 | 0.848 |
| 0.175 | 0.8602 | 0.8755 | 0.982 |
| 0.225 | 0.9422 | 0.9106 | 1.035 |
| 0.275 | 0.9528 | 0.9419 | 1.012 |

So the divergence is a **factor 2.3 in the first ΔR bin**, dies out by ΔR ≈ 0.175, and is present in
the MEASURED points — it is not a fit artefact. The `expo` fits reproduce it: same sign
f(0) = 0.304 (A = −0.696 ± 0.035, λ = 0.137 ± 0.006), opposite sign f(0) = 0.816
(A = −0.184 ± 0.002, λ = 0.258 ± 0.003); Medium behaves identically.
**Step 4 (single leg) shows nothing of the kind** — the two signs agree at every ΔR
(f(0) 1.047 same sign vs 1.148 opposite sign). So the effect lives ONLY in the two-leg/joint term.

**Consequence for the analysis (not yet acted on).** The nominal ε_ΔR is sign-INTEGRATED and is
therefore ≈ the opposite-sign curve, because opposite-sign pairs are 87 % of the sample by count
(2 169 349 vs 323 349). Applying it to same-sign pairs **under-corrects them by up to a factor 2.3
below ΔR ≈ 0.15**. That is not academic: the ΔR > 0.05 cut was removed from the signal selection
(`remove_dr_cut_signal_selection.md`), so small-ΔR pairs are in the signal region, and the
background subtraction is **OS − SS**, i.e. the same-sign spectrum enters the result directly.

**TWO COMPETING EXPLANATIONS — both plausible, NOT resolved, do not adopt either yet.**
1. **A genuine L1 close-by-RoI effect.** In the toroid, two same-charge muons bend the same way, so
   a close same-sign pair stays close in the muon spectrometer and its two L1 RoIs merge into one —
   the pair then fails a 2-of-2 trigger. Opposite-charge muons bend apart and separate their RoIs.
   This predicts exactly what is seen: an effect in the JOINT term only (Step 3) and none in the
   single-leg marginal (Step 4), dying out once the pair is wider than an RoI.
2. **An ε_MC parameterisation artefact in the inverse weight** (raised by the plot reviewer, and it
   is a good point). For a close SAME-sign pair q·η₁ ≈ q·η₂, so the 1/(ε₁ε₂) weight **squares** any
   mis-parameterisation of ε_MC in that q·η region; for a close OPPOSITE-sign pair
   q·η₁ ≈ −q·η₂ and the two errors partly **cancel**. Given **R23 is still open on exactly that
   question** (the top-pair-pT plateaus 20–42 % high, ε_MC pT clamp vs coarse q·η binning), this
   cannot be waved away.

**★ THE TEST WAS RUN, AND EXPLANATION 2 IS REFUTED (2026-08-11).** A weight-only numerator
(`numraw`, the SAME trigger-passing node weighted by the plain MC weight) was booked beside the
inverse-weighted one and Step 3 refilled, giving the RAW joint trigger probability
P(both fire | ΔR) with **no ε_MC anywhere in it**. pp_full, Tight:

| ΔR | P_raw same sign | P_raw opposite sign | raw ratio | ε_ΔR ratio (inverse-weighted) |
|---|---|---|---|---|
| 0.025 | 0.2375 | 0.5564 | **0.427** | **0.425** |
| 0.075 | 0.3135 | 0.5489 | 0.571 | 0.595 |
| 0.125 | 0.4714 | 0.5762 | 0.818 | 0.844 |
| 0.175 | 0.5595 | 0.5789 | 0.967 | 0.978 |
| 0.225 | 0.6316 | 0.6147 | 1.027 | 1.030 |
| 0.325 | 0.6420 | 0.6503 | 0.987 | 1.025 |

**Both columns are ratios of the UN-plateau-normalized curves** (that is the like-for-like
comparison: the raw probability has no plateau normalization to apply). They agree bin by bin to
**≲4 %** — largest deviation 4.2 % at ΔR = 0.075, 3.8 % at ΔR = 0.325. For reference the
plateau-normalized ε_ΔR ratio, which is what the ratio canvases draw, is 0.4273 / 0.5977 / 0.8477 /
0.9825 / 1.0347 / 1.0296 — the same story. **Two conclusions, and they point the same way:**
- The split is **already present in the raw trigger probability**, so the 1/(ε₁ε₂) weight neither
  creates it nor materially changes it ⇒ **explanation 2 (ε_MC mis-parameterisation squaring for
  close same-sign pairs) is ruled out.**
- The inverse weighting is precisely what divides out the (pT, q·η) kinematic dependence, and it
  does **not** remove the split ⇒ a different (pT, q·η) MIX for close same-sign vs opposite-sign
  pairs is ruled out as well.

**Unexplained, noticed on the new raw-probability figure (2026-08-11) — STILL OPEN; the original
candidate explanation was REFUTED 2026-09-08 (round-13 plot review).** The opposite-sign raw
probability **dips around ΔR ≈ 0.65–0.75** while the same-sign one stays flat, pushing P_SS/P_OS up
to ≈1.4 there. Re-measured on the round-13 outputs: P_OS falls from ≈0.66 at ΔR ≈ 0.38 to **0.455 at
ΔR ≈ 0.72** and recovers to ≈0.61 by ΔR ≈ 0.83, against P_SS flat at ≈0.61–0.63; the error bars are
≈±0.01, far smaller than the excursion, so it is not a fluctuation.

- **REFUTED: the OS-only ntuple-stage resonance veto.** The original entry proposed
  `project_os_resonance_veto` (the OS pair tree is resonance-vetoed at the ntuple stage, the SS tree
  is not). That cannot be the cause **for this sample**: the MC trig-eff fill reads the
  `..._no_data_resonance_cuts_mc_trig_...` trees (`FillMCTrigEffHists.cxx:165–166`), i.e. the
  variant produced *without* the veto, and neither `FillMCTrigEffHists.cxx` nor
  `MCTrigEffPairSel::Step3PairSelection()` (`Utilities/MCTrigEffPairSelection.h:62–70`) contains any
  `minv` cut at all — verified by grep, 2026-09-08. Both samples therefore keep their resonances.
- **Candidate that survives, NOT yet tested:** precisely *because* the resonances are kept, the OS
  sample contains J/ψ (prompt or from B — the argument rests only on an OS-only narrow resonance
  at m = 3.10 GeV, not on where it came from) while the SS sample cannot. For a two-body decay
  ΔR ≈ 2m/p_T^pair, so a J/ψ (m = 3.10 GeV) lands at ΔR ≈ 0.65–0.75 for
  p_T^pair ≈ 8–9.5 GeV — the most populated pair-pT cell — and each leg is then ≈4–5 GeV, i.e. sitting
  in the mu4 turn-on where the per-leg efficiency is lowest. That would depress the OS **raw** joint
  probability in exactly this ΔR window, and would be divided out by the 1/(ε₁ε₂) weighting, which is
  consistent with ε_ΔR showing no corresponding feature. Testable by re-filling the raw-probability
  figure with the J/ψ mass band excluded, or by splitting it in pair pT.
- It sits far above the ΔR ≲ 0.15 region this entry is about and does not affect the conclusion
  below, but it is a real asymmetry between the two samples and **should be understood before the
  per-sign corrections are used**.
Figure: `step3_dr_correction/step3_raw_joint_trigger_probability_by_sign.png`.

⇒ **a genuine charge-dependent two-body trigger correlation**, confined to the JOINT term (Step 3
only, nothing in Step 4's single-leg marginal) and to ΔR ≲ 0.15. That is the signature of L1
close-by-RoI merging: same-charge muons bend the same way in the toroid, so a close same-sign pair
stays close in the muon spectrometer and its two RoIs merge into one, failing a 2-of-2 trigger;
opposite-charge muons bend apart. Explanation 1 stands. (Independent support, RECOVERY SCALE ONLY: the
opposite-sign curve recovers by ΔR ≈ 0.3, the same scale as the Run 2 close-by-RoI correction ρ_ΔR,
which was measured on opposite-charge dimuons. **The magnitude may NOT be compared** — that KB entry
states explicitly "do not infer the size of our ε_dR from these Run 2 ρ_ΔR curves — different
system, energy, and trigger" (`.claude/kb/physics/detector/atlas_run2_muon_trigger.md` §12.4). No
Run 2 same-sign analog exists, so that branch is RUN2-CROSSCHECK UNVERIFIED.)

**OPEN — USER DECISION, with a rerun blast radius.** The measurement is settled; how the analysis
USES it is not. The nominal ε_ΔR is sign-integrated and ≈ the opposite-sign curve, so same-sign
pairs are under-corrected by up to ×2.3 below ΔR ≈ 0.15. The natural fix is to apply the
**sign-dependent** ε_ΔR (same-sign correction to the same-sign spectrum, opposite-sign to the
opposite-sign one) in the pp 2mu4 weight — the per-sign fits already exist, at both working points
and for all three fit functions. **But coverage is NOT free, and the user should know that before
deciding:** the same-sign curve falls to f(0) ≈ 0.1–0.3, which a form pinned to 1 at large ΔR can
only reach with a large negative amplitude, so it frequently goes unphysical and is rejected.
Same-sign cells with `fit_ok = 0` (Step 3, corrected — see the ⚠ below):

| WP | sign-integrated | same sign, `expo` | same sign, `polyu_fixedRp` | same sign, `interp` |
|---|---|---|---|---|
| Tight | **13**/72 | **32**/72 | **38**/72 | **15**/72 |
| Medium | **12**/72 | **31**/72 | **37**/72 | **15**/72 |

**Only the interpolation comes close to covering the same-sign cells at the nominal rate**, and even
it rejects 15 against the nominal 13.

**⚠ CORRECTED (2026-08-11, code review).** An earlier version of this entry quoted 25 / 31 / 11
against 12, taken from the fit reports' `cells marked UNUSABLE` line. **That line under-counted:**
`n_unusable` was incremented only for cells that reached the end of the fit loop, while the two
early-exit paths (`no plateau -> cell skipped`, `too few points -> no fit`) persisted `fit_ok = 0`
without counting it. The reported number therefore disagreed with the `h_stepN_fit_ok` map it
described — worst in the same-sign series, i.e. exactly the one this decision is about. Fixed in
`fit_dr_corrections.cxx`; every report now equals a direct bin count of its own map. Step 4 moved
the same way (Tight same sign: expo 5 → 12, `polyu_fixedRp` 4 → 11, `interp` 2 → 6).

**A second thing on the table before deciding — a definitional mismatch (code review, 2026-08-11).**
The MC same-sign/opposite-sign split is on the **TRUTH** charges (`MuonPairMC.h`
`truth_same_sign = (m1.truth_charge == m2.truth_charge)`, selected by `GetMuPairIsSameSign`), but
the consumer this proposal has in mind — applying the same-sign correction to the same-sign spectrum
in the OS − SS subtraction — classifies DATA pairs by **RECO** charge. The leakage is small
(muon charge mis-identification is ≲0.1 %, i.e. ≲1 % relative on the 13 % same-sign sample) and it
does not change the conclusion, but it is a real definitional difference and it should be on the
table rather than discovered later. That changes the trigger-corrected yields and therefore
**requires a crossx refill + full replot** (`signal_selection_change_impact.md`), which is exactly
why it is not done here. It should be bundled with the refill the `w_trig = 0` gap bug already
forces (`pp_trig_eff_highpt_jump.md`). **Raised to the user 2026-08-11.**

**(c2) THREE pre-existing bugs found while doing (c), all fixed.**
- `plot_dr_correction_fits.cxx` built `hist_path` **without** `MCTrigEffPairPt::FileSuffix()` — the
  only stage in the chain that did — so under `MCTRIGEFF_PAIRPT_4BIN=1` it paired 8-bin histograms
  with the 4-bin fit file. The fit stage has explicit cell-count and edge guards for exactly this;
  the plot stage had none.
- **`root_file_has_objects()` in `run_dr_correction_fits.sh` was a silent no-op.** It fed a heredoc
  to `root -l -b -q`, and `root -q` with no macro argument **quits before reading stdin**, so the
  function returned success for every input — *every* "file missing/incomplete" artefact check in
  that pipeline had been passing unconditionally. Verified both ways before and after the fix. This
  is the standing "validate ARTEFACTS, not exit codes" rule failing at the validator itself.
- The driver's PNG-count check still expected the stale literal **5** (from "4 pair-pT bins +
  inclusive"); it is now derived from the plateau map's x-axis.

**(d) Statistics of the 8-bin pair-pT axis** — new macro `write_mc_pair_statistics_tables.cxx`,
reading three TH2Ds per sign (`h_mc_paircount / pairsumw / pairsumw2 _vs_pt_eta_{ss,os}`) that
`FillMCTrigEffHists.cxx` now books **on the Step-3 denominator node**, so the tables describe
exactly the selection the ΔR correction is measured under and cannot drift from it. Six CSVs per WP
in `mc_statistics_pt8bins[_medium]/`. **Key numbers (Tight, inside the binned range):**

| pair pT [GeV] | same-sign pairs | σ_SS [nb] | opposite-sign pairs | σ_OS [nb] |
|---|---|---|---|---|
| 8.0–11.5 | 144 721 | 1.2349 | 776 604 | 7.9344 |
| 11.5–16.6 | 99 382 | 0.6734 | 773 805 | 6.5984 |
| 16.6–24.0 | 46 722 | 0.2365 | 383 473 | 2.5811 |
| 24.0–34.6 | 19 953 | 0.06766 | 150 976 | 0.7644 |
| 34.6–50.0 | 8 455 | 0.01699 | 55 792 | 0.18666 |
| 50.0–72.1 | 3 199 | 0.003650 | 20 623 | 0.038618 |
| 72.1–104.0 | 792 | 0.000580 | 6 788 | 0.007115 |
| 104.0–150.0 | **125** | 0.0000812 | **1 288** | 0.000959 |
| total | 323 349 | 2.2338 | 2 169 349 | 18.1117 |

(Medium: SS 352 665 / 2.4626 nb, OS 2 385 330 / 20.1386 nb. Cross sections are in **nb** — AMI
cross sections are nb, ×1000 to compare with pp data in pb⁻¹.)
- **★ This directly sizes the R23 open question.** The top pair-pT bin holds **125 same-sign and
  1 288 opposite-sign pairs, spread over 9 pair-η bins** — of order 10–150 pairs per cell. Per-cell
  corrections in the top bins are statistics-starved regardless of what the R23 investigation
  concludes about the ε_MC pT clamp vs the coarse q·η parameterisation.
- **Also measured, and new:** only **40.77 %** of selected same-sign pairs and **66.61 %** of
  opposite-sign pairs fall inside `pair pT ∈ [8, 150] GeV`; the rest sit **below 8 GeV** and are
  outside the coarse axis entirely. The ΔR correction is therefore measured on well under half the
  same-sign sample. Whether the coarse pair-pT axis should extend below 8 GeV is a binning question
  and hence a **user decision** — flagged here, not acted on.


### R27. AUDIT of the ΔR-correction error bars — what they ARE, and the three things they are NOT
### (2026-08-11, user question; pp24 FULL, Tight, current round-9 outputs)

**Question.** How are the error bars on the inverse-weighted ΔR corrections computed, why not
`TH1::Divide`, and are they right — in particular, how much of the offset-from-1 and of the
bin-to-bin scatter is statistics, and how much is the ε_MC uncertainty propagated through the
inverse weighting?

**(a) What is computed (one implementation, `plotting_codes/trig_effcy/mc_based/dr_correction_ratio.h`
`SetConditionalRatioErrors`).** With `D = Σ_all w`, `N = Σ_fired a_i`, `a_i = w_i/(ε₁ε₂)` (Step 3)
or `w_i/ε_i` (Step 4), `R = N/D`, the sample is held FIXED and only the Bernoulli trigger decisions
fluctuate: `Var(N) = Σ a_i² p_i(1−p_i) + 2Σ_pairs a₁a₂(p₁₂ − p₁p₂)`, `p_i = ε_i·R`, estimated from
histograms booked in `FillMCTrigEffHists.cxx` as `Var = A − R·B (+ covP − R²·covQ)` with
`A = Σ_fired a²`, `B = Σ_fired a²·ε`, and `e_R = √Var/D`. `A − R·B` is exactly unbiased for
`Σ a²p(1−p)` under `p_i = ε_i R`. Unweighted limit A = B = N ⇒ `e_R = √(R(1−R)/D)` — textbook
binomial. covP/covQ exist only for Step 4 (both legs of a pair share one ΔR bin); Step 3 has one
Bernoulli trial per entry, so they are null.

**(b) Why not `TH1::Divide`.** Option-less Divide propagates `e_R = R√((e_N/N)²+(e_D/D)²)`, which
assumes numerator and denominator INDEPENDENT — but the numerator is a re-weighted SUBSET of the
denominator. Ratio of the two, in the uniform limit, is `√((1+p)/(1−p))`, p the effective pair
efficiency. Option `"B"` / `TGraphAsymmErrors::Divide` is also unusable, for a different reason:
the entries are not counts (weights 1/(ε₁ε₂) > 1, so N > D is allowed and the binomial/Bayes
estimators are undefined — R16 already hit exactly that as a silent empty-graph bug).
**Re-measured today on the current round-9 pp files (inclusive Step 3, `h_mc_dr_full`):**
`σ_Divide/σ_cond = 1.51–1.87`, and over the plateau window [2, 3.5] the constant-fit
**χ²/ndf = 1.13 with the conditional errors vs 0.39 with Divide** (Step 4: 0.73 vs 0.13).
The conditional bars are the ones the scatter actually supports.

**(c) Scale of the statistical bar — it is NOT √(R(1−R)/n).** Because R ≈ 1 while the per-pair
firing probability is p ≈ ε₁ε₂ ≈ 0.5, the fluctuating quantity is the trigger decision, not the
ratio: `σ_R ≈ √(R(1−εR)/(ε·n_eff))`. At ΔR = 0.125 (inclusive, `h_mc_dr_full`) that is
**0.00117 measured vs 0.00046** for the naive √(R(1−R)/n_eff) — inverse weighting costs a factor
≈ 1/√ε ≈ 1.4 on top, plus the MC event-weight spread. n_eff there is 5.9e5 pairs.

**(d) Current inclusive numbers (round 9, plateau window [2, 3.5], pp24 FULL Tight):**
Step 3 plateau **0.99208 ± 0.00087 (9.1σ from 1)**, Step 4 **0.99621 ± 0.00040 (9.5σ)**. (R12(f)'s
"0.9911 ± 0.0004, 22σ" was Step 4 on the retired [1,4] window, 2026-08-03 — superseded, same
conclusion: the offset is real, not a fluctuation.) The offset is a fit-quality offset in ε_MC, and
that is what §3.3 diagnostic 2 says it is.

**(e) ★ THE MISSING TERM — ε_MC uncertainty is NOT in any bar.** The bars condition on the ε map
being EXACT. Measured size of the ε_MC uncertainty itself (from
`single_mu_effcy_pT_fit_mc.root`): fit-parameter statistical error **0.15–1.0 % absolute**;
fit-FORM residual (point − fit)/fit rms **0.7–2.7 % in the central q·η bins**, **6.5–8.2 %** in
q·η ∈ (1.0,1.5) and (2.0,2.3), **21 %** in the forward (−2.4,−2.0) bin that `kVetoFwdLowPt` removes.
So per leg the real ε uncertainty is the FORM one, ~1–2 %, not the fit-stat one.
- **Why it mostly cancels.** `R ∝ Σ_fired w/(ε₁ε₂)`, so ε → ε(1+δ) gives
  `δR/R = −⟨δ₁+δ₂⟩_fired(ΔR)`. After plateau normalization only the DIFFERENCE survives:
  `δ(ε_ΔR)/ε_ΔR = −[⟨δ₁+δ₂⟩(ΔR) − ⟨δ₁+δ₂⟩(plateau)]`. A global ε normalization error cancels
  EXACTLY; only the ΔR-dependent kinematic-mix difference is left.
- **Measured leverage** (nominal vs the corrected-MC files, i.e. ε_MC → ε_corr ≈ ε_data, a real
  10–25 % per-leg perturbation with a realistic barrel-vs-endcap shape): the RAW plateau moves
  **0.99208 → 0.97193 (−2.0 %)**, but the plateau-NORMALIZED curve moves by only **+0.019 at
  ΔR = 0.025, +0.015 at ΔR = 0.125, < 0.004 above ΔR ≈ 0.5**. Leverage ≈ **0.1 at ΔR → 0, ≈ 0 in
  the plateau** — this is the quantitative version of R16's "insensitive to the overall ε
  normalization", now resolved in ΔR instead of as a per-cell median.
- **But it is not negligible where the physics is.** At ΔR = 0.125 that shift is **12.6 σ_stat**;
  at ΔR = 0.025 it is 5.9 σ_stat. Folding the leverage 0.1 onto the ~1–2 % real form uncertainty
  gives **~0.1–0.2 % at ΔR → 0**, i.e. 0.3–0.6× the statistical bar there — comparable, one-sided,
  and fully correlated across ΔR bins. **⇒ an ε_MC-parameterization systematic on ε_ΔR should be
  quoted; today none is.** It is NOT the same thing as the plateau-window systematic
  `|plateau[2,3.5] − plateau[1,4]|`, which sizes *where the reference is taken*, not *how well
  ε_MC(pT, q·η) is parameterized*.

**(f) Second missing term — the plateau normalization error is dropped.**
`fit_dr_corrections.cxx:657` does `r->Scale(1.0/plateau)`, which scales contents AND errors by the
same factor, so the plateau's own uncertainty never enters the normalized points, the fitted
parameters, or `f_at_0`. For the per-cell fits this is the LARGER of the two normalization
effects: `dr_correction_plateaus_pp24_full.root` `h_step3_plateau` has **median per-cell relative
error 1.53 % over 72 cells** (inclusive: 0.09 %). Dividing content and error together is the right
treatment for a FIT (a correlated normalization must not be added per point), but the resulting
±1.5 % must be carried separately as a correlated systematic on every per-cell correction, and is
not. This is the same quantity the round-10 `no_plateau_correction` variant exists to sidestep.

**(g) Third — event-level correlation, checked and NEGLIGIBLE.** pp's Step-3 numerator condition is
the EVENT-level `pass2mu4`, so two selected pairs in one event share ONE Bernoulli trial, which the
diagonal `A − R·B` does not model (the code comment flags it as unmodelled). Measured on the pp24
FULL pair trees (Tight, pT > 4, |η| < 2.4): **2.52 % of opposite-sign and 0.46 % of same-sign
selected pairs come from events with more than one selected pair** (1.013 / 1.003 pairs per event).
Worst-case error underestimate √(1 + 0.025) = **1.3 %** — below the 1–2 % agreement of (b). Not
worth modelling.

**Verdict.** The bars are CORRECT as what they claim to be — conditional (binomial) statistical
errors on an inverse-weighted efficiency ratio, with the sample held fixed — and the plateau
χ²/ndf ≈ 1.13 confirms the bin-to-bin scatter matches them. They are NOT the total uncertainty on
the deliverable: the 9σ offset from 1 is a genuine ε_MC parameterization effect (correctly
normalized away), and the two correlated normalization terms in (e) and (f) are real, of order
0.1–0.2 % (shape) and 1.5 % (per-cell plateau), and are currently unquoted. **New Remaining-Work
item; needs a user decision on how to size them** (see Remaining Work 9).

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

**Recommendation (REVISED 2026-08-04 by the user — supersedes the first pass):**
**nominal `expo`**, `f(ΔR) = 1 + A·exp[−(ΔR/λ)^p]`; **backup `polyu_fixedRp`**; cross-check
`interp`. **The two power laws are DROPPED.**
- **Why the power laws are out — they are NOT SMOOTH.** `f = 1 + A·u^n` is continuous in *value*
  at R_p but its slope `−A n u^(n−1)/R_p` **diverges for n < 1** — a cusp. The fitted n rails to
  its lower limit 0.2 in a large fraction of cells: pp Step 3 **11/37**, Step 4 9/37; overlay
  Step 3 **19/37**, Step 4 17/37. That is the typical case, not an edge case. *(The first pass
  passed them on a "flat beyond R_p" audit that only checked the VALUE beyond R_p, never the
  SLOPE — the audit missed exactly this.)*
- **Why `expo` is back in, as nominal.** It is smooth everywhere and → 1 asymptotically. It was
  first rejected for not being *exactly* 1 beyond R_p; **that was wrong** — "flat for ΔR ≳ 0.5"
  is an ESTIMATE, not a strict bound, and a cell whose behaviour departs wildly from it is
  something to investigate, not grounds to reject the functional form. Measured residual
  |f−1| at ΔR = 1 over pp cells: **median 0.0000, 75th pct 0.0000, 90th pct 0.0143 (Step 3) /
  0.0000 (Step 4)**, max 0.2051 / 0.0332 in one or two outlier cells. The large overlay
  residuals (75th pct 0.94) sit in cells already marked `fit_ok = 0`.
- **Why `polyu_fixedRp` is only the backup.** It is C¹ at R_p by construction and follows the
  non-monotonic small-ΔR shape, but its higher-order terms can equally absorb shapes that are
  **procedure artefacts rather than physics** — so it serves as the cross-check on `expo`, not
  the default.
- χ²/ndf, inclusive cell (Tight): `expo` pp 7.80 (S3) / 15.39 (S4), overlay 1.355 / 1.231;
  `polyu_fixedRp` pp 5.90 / 25.17, overlay 1.343 / 1.205. On pp no 3–4-parameter form reaches
  ≈1 (see the caveat below); on the overlay everything is ≈1.1–1.4.
- **Rejected-method outputs REMOVED (2026-08-04, by the user):** the
  `step{3,4}_dr_fit/powerlaw_{fixedRp,floatRp}/` directories (8 dirs, 89 MB) and their 16
  `dr_correction_fits_*powerlaw*.root` files are gone; verified that only `expo`,
  `polyu_fixedRp` and `interp` survive (24 fit ROOT files, 6 plateau files intact). Regenerable
  if ever needed via
  `METHODS="powerlaw_fixedRp powerlaw_floatRp" SKIP_MEASURE=1 bash pipelines/run_dr_correction_fits.sh`.

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

### R17. TWO COEXISTING pair-pT BINNINGS — found and eliminated (2026-08-04, user)

**The user spotted, from the plot text, that the pair-pT ranges in the Step-3 slices panel
(40.6–120 GeV) did not match those in the plateau tables (50–150 GeV).** They did not: two
pair-pT binnings had been running side by side inside the same plot directory.

**History (git, not inference):**
| what | when | binning |
|---|---|---|
| `ParamsSet::pair_pt_coarse_bins` created (`3bfc9cc`, "single source of truth") | 2026-07-01 | `{8, 15, 27, 50, 150}` — **only version in history; never edited** |
| MC trig-eff Step-3 slices panel + its fine 2D (`3e6f75e`, `3c2c846`) | 2026-07-10 | `pT_bins_120` grouped → 8–13.75 / 13.75–23.63 / 23.63–40.62 / 40.62–120 |
| Coarse 3D enters the fill (`68be85c` round-5 #4; `8c9b674` Step 4) | 2026-07-22 / 07-28 | `{8, 15, 27, 50, 150}` |

So **nothing ever changed a binning** — the split was created on 2026-07-22 when the pair-η
panels and plateau tables were added reading `ParamsSet`, while the older slices panel kept its
own grouping. From then until 2026-08-04 the plateau tables and the panel beside them described
**different cells**. Every plot rendered and every number came out: a silent failure.
The user had also previously asked for `pair_pt_coarse_bins` to be
`{8, 13.8, 23.6, 40.6, 120}`; **that was never carried out** (git shows one version only).

**Resolution (user, 2026-08-04):**
- `ParamsSet::pair_pt_coarse_bins` is now the **pT_bins_120 group edges
  `{8, 13.7502, 23.6334, 40.6205, 120}`**, *derived* from `pT_bins_120` (indices 0/3/6/9/15) in
  the constructor rather than retyped, so the coarse and fine axes cannot drift apart. pp data
  barely reaches beyond 120 GeV, so 120 is the physical top of the axis.
- `ParamsSet::pair_pt_coarse_bins_pt150 = {8, 15, 27, 50, 150}` survives as an **opt-in** variant.
- The inconsistent fine 2D `h_mc_dr_zoom_vs_pair_pt_*` is **deleted**; every pair-pT view
  (Step-3 slices, Step-3/4 pair-η panels, plateau tables, fits) now projects the **same 3D**.
- ⚠ **Blast radius, NOT yet rerun:** `pair_pt_coarse_bins` is also read by the crossx hist filling
  (`RDFBasedHistFillingPP.cxx:555`, `RDFBasedHistFillingPbPb.cxx:1072`), so **crossx outputs are
  now stale**. Awaiting a user decision: rerun crossx on the canonical binning, or pin those two
  call sites to `_pt150`.
- Rule recorded in `.claude/CLAUDE.md` §Binnings (BLOCKING) and in `.claude/conventions/
  atlas-plotting.md`.

**Effect on the results.** The inclusive numbers are unchanged (the inclusive cell integrates
over pair-pT): pp_full Step 3 = 0.9723 ± 0.0008, Step 4 = 0.9862 ± 0.0004 (Tight); the fit-method
χ²/ndf ranking of R15 is unchanged. The **per-cell** picture improves — on the canonical binning
**no pp24 cell fails** the 0.15 tier (two are flagged in the 0.10–0.15 band, both Step 3 in the
top pair-pT bin: η_pair[−2.4,−2.0) = 1.1490 ± 0.0105 and η_pair[1.0,1.5) = 0.8648 ± 0.0277;
Step 4 has none), and `h_stepN_fit_ok = 0` marks **0** pp24 cells vs **31/36** (Step 3) and
**13/36** (Step 4) overlay cells.

### R18. Publication standard for the plots (2026-08-04, user)

Two standing rules, now enforced by `/review-plot` (`.claude/conventions/atlas-plotting.md`
criteria **P1** and **P2**):
- **P1 — no illustrative text on a plot.** Anything a physics professor being shown the figure
  should not see must not be drawn on it: explanatory/tutorial sentences, pointers to files or
  code ("see fit_report.txt", ROOT/histogram names, internal mode identifiers such as
  `polyu_fixedRp` or `do_step4`), drawing asides ("(points only)", "(read back)"), hedging prose.
  **If a number is worth showing, draw the number.** That material belongs in this doc and in the
  summary to the user. All such text was removed from the three MC trig-eff plot macros.
- **P2 — a fit must always carry its exact equation**, placed above the parameter values with the
  same symbols and any auxiliary variable defined; the legend entry is simply `fit`; internal
  method names never appear on the canvas; fixed parameters are marked `(fixed)` rather than
  printed as `± 0`.

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

10. **★ OPEN (R27, 2026-08-11) — TWO UNQUOTED CORRELATED UNCERTAINTIES ON ε_ΔR. Needs a user
    decision on how to size them.** The published bars are conditional STATISTICAL errors only
    (and are correct as such: plateau χ²/ndf 1.13). Missing:
    (i) **ε_MC parameterization** — the ε map is treated as exact. Its real (form) uncertainty is
    ~1–2 % per leg; measured leverage onto the plateau-normalized correction is ≈ 0.1 at ΔR → 0 and
    ≈ 0 in the plateau ⇒ ~0.1–0.2 % at small ΔR, one-sided and fully ΔR-correlated, i.e. 0.3–0.6×
    the statistical bar exactly where the physics lives. Distinct from the plateau-WINDOW
    systematic. A ready-made evaluation exists: the nominal-vs-corrected-MC difference of the
    normalized curves (R27(e)).
    (ii) **plateau normalization** — `fit_dr_corrections.cxx:657` scales contents and errors
    together, so the plateau's own error (median **1.53 %** per cell, 72 cells, pp Step 3) never
    reaches the fitted parameters or `f_at_0`. Correct for the fit, but it must travel as a
    correlated systematic and today does not. Interacts with the round-10
    `no_plateau_correction` variant, which removes the term rather than sizing it.

## Latest Stage

**2026-09-08/09 (round 14) — ✅ DONE for the FITS; the MC closure is deliberately NOT rerun (user
decision, see below). The Step-3 `polyu_fixedRp` (polynomial) fit is constrained so that at
ΔR = 0 the efficiency cannot exceed the plateau. Full results and numbers in R34.**

Delivered: the polynomial reparametrized in `A ≡ f(0) − C` so the user's requirement is the single
fit limit `A ≤ 0` (R34(a)), Step 3 only; the 30 Step-3 polyu fit + plot configurations regenerated
(2 WPs × 5 plateau modes × 3 sign series) — **1482 fitted cells, zero `A > 0`, `f(0) ≤ C` in all of
them**, flat at `C` beyond `R_p` to 0.00e+00 in all 30 blocks; and two `/review-analysis-code`
fixes (R34(g)): the report's `AT LIMIT` tolerance caveat, and a parametrization guard in the
plotter that refuses to draw the new equation over an old fit file.

**Three things the user should know.**
1. **Two side-effects were measured and put to the user, who chose the requirement as stated
   ("ΔR = 0 only")** — R34(f): the polynomial may still exceed the plateau BETWEEN 0 and `R_p`
   (47 % of cells by >1 %, worst 11.5), and in the `nocorr` family 11 of the 95 railed cells met
   the constraint by collapsing the free baseline `C` instead — 9 of those are caught by the
   existing `DrCorrPlateauUsable` screen and rerouted to the `interp` tier, but **2 clear it**
   (`C ≈ 0.78 ± 0.53`, `σ_C/C ≈ 0.68`) and are delivered with `f/C` up to 1.94.
2. **The MC closure could NOT be rerun** — R34(h): a concurrent, uncommitted workstream moved the
   canonical pair-pT axis 8 → 9 GeV, so `DrCorrectionCascadeEvaluator` throws on every Step-3 fit
   file on disk. The staleness predates this change. The user directed this round to stop at the
   fits. **`mc_trigeff_dr_binning_approaches.md`'s four-approach χ²/ndof figures are now stale for
   two independent reasons and must not be used for the open approach decision** — flagged at the
   top of that doc's Latest Stage.
3. **The Tight Step-3 fit set is now internally INCONSISTENT across methods** — R34(d). That same
   workstream refilled the Tight Step-3 histogram on its new 9 GeV axis at 2026-09-08 18:48,
   between this round's two runs, so the Tight `polyu` fits are on 9 GeV while the Tight
   `expo`/`interp` fits are still on 8 GeV; Medium is self-consistent at 8 GeV. Every consumer
   THROWS rather than silently mixing them, and the fix is the refill + refit of all three methods
   that item 2 already requires.

*The plan, as written before the work, follows.*

---

**2026-09-08 (round 14) — plan. Extend the Step-3 shape restriction to the `polyu_fixedRp`
(polynomial) fit: at ΔR = 0 the efficiency may not exceed the plateau.**

*User request (verbatim intent): "For the step 3 polynomial fit, require that when dR = 0, the
efficiency cannot have value exceeding plateau."* This is requirement (1) of R33(b) — the
small-ΔR plateau must lie BELOW the large-ΔR one — now imposed on the BACKUP form as well as on
the nominal `expo`. Nothing else about the polynomial is constrained (no monotonicity, no
turning-point condition): the user asked for the ΔR = 0 value only.

**How it maps onto the parameters.** `f = C + u²(a₂ + a₃u + a₄u²)`, `u ≡ max(0, 1 − ΔR/R_p)`,
`u(0) = 1`, so `f(0) = C + (a₂ + a₃ + a₄)`. The requirement is therefore the LINEAR constraint
`a₂ + a₃ + a₄ ≤ 0`, which is NOT a box constraint on any single existing parameter and so cannot
be imposed as a MINUIT limit on `a₂`. It is imposed by REPARAMETRIZING the same function family
so that the constrained quantity IS a parameter:

    A ≡ f(0) − C = a₂ + a₃ + a₄     (so a₂ = A − a₃ − a₄)
    f = C + u²·[ A + a₃(u − 1) + a₄(u² − 1) ]

`a₃` and `a₄` keep their exact meaning (still the coefficients of `u³` and `u⁴`); only the
leading coefficient is traded for `A`, which is the SAME symbol and the SAME meaning `expo`
already uses (`f(0) − C`). Then the restriction is the closure `A ≤ 0`, i.e. limits `[−50, 0]`
in place of the pre-existing `[−50, 50]`. A cell railed at `A = 0` is the boundary case — the
unconstrained fit wanted a small-ΔR ENHANCEMENT — and delivers a curve that is flat at `C` at
ΔR = 0; it is flagged `AT LIMIT: A` per cell exactly as `expo`'s is.

**Step 4 is EXCLUDED, for the same physics reason as `expo`** (R33(b), R4/§3.4): the single-leg
correction is physically an ENHANCEMENT at small ΔR, so `A > 0` there. `restrict_shape =
(step == 3)`.

**Plan (per §3.3 Step 3; reviewer: `/review-analysis-code` for the fit code, `/review-plot` for
the regenerated plot sets):**
1. `fit_dr_corrections.cxx` — reparametrized `polyu_fixedRp` formula (BOTH plateau families),
   parameter names, seed (with the same inside-the-limits clamp `expo` got), the `A ≤ 0` limit
   for step 3 only, and the fit-report + ROOT-`provenance` shape notes extended to the method.
2. The two places that DUPLICATE the formula string / parameter names, which must move with it:
   `combine_dr_correction_fits.cxx::NocorrMethodFormulas()` and
   `plot_dr_correction_fits.cxx::MethodFormulaTex()`.
3. The two closure figures that print the polynomial's equation as prose:
   `plot_mc_trig_eff_closure.cxx`, `plot_mc_trig_eff_closure_compare.cxx`.
4. Standalone check that the reparametrization is an identity on the OLD fitted parameters
   (same curve, same `f(0)`) before any production rerun.
5. Rerun `run_dr_correction_fits.sh` for `pp_full`, STEPS=3, both WPs, all methods × 3 sign
   series × 5 plateau modes (`SKIP_MEASURE=1` — the plateaus are an input and did not change);
   plots regenerated by the same driver.
6. **Blast radius beyond this doc:** the delivered ε_ΔR is the `expo → polyu → interp` cascade
   (`mc_trigeff_dr_binning_approaches.md` §PP-3), so every cell where `expo` is NOT accepted and
   `polyu` IS now changes. That doc's four-approach MC closure χ²/ndof — currently the evidence
   in front of the user for an OPEN approach decision — therefore goes stale and is rerun
   (`run_mc_trigeff_closure.sh`, four modes, both WPs). The pp24 cross-section is unaffected:
   ε_ΔR is still not propagated to it, and `DrCorrCrossxMethod()` is `expo`.
7. Report how often the new limit binds (cells railed at `A = 0`, per WP × mode × sign), doc +
   INDEX, commit.

*The round-13 entry, now closed, follows.*

---

**2026-09-08 (round 13) — ✅ DONE. Rerun on the NEW fiducial cut set + the NEW pair-level
`|η^pair| < 2.2`, and the Step-3 `expo` fit CONSTRAINED so it cannot come out flipped. Both
reviews PASS (`/review-analysis-code` iter 4, `/review-plot` iter 2); full results and numbers in
R33.**

Delivered, pp24 only (user-confirmed scope): the pp24 DATA single-muon mu4 T&P refit at both WPs on
the new probe gap cut and the new coarse q·η top bin `[2.0,2.2)`; the MC `pp_full` Steps 1–4 refill
+ turn-on refit + full plot set at both WPs; and the Step-3 ΔR-correction fits re-measured and
refitted for 3 methods × 3 sign series × 5 plateau modes at both WPs. Step-1 MC ε(mu4) =
**0.8025 (μ⁺) / 0.7961 (μ⁻)**. Every pair-η cell is now `[−2.2,−2.0) … [2.0,2.2)` (folded:
`[0,1), [1,2), [2,2.2)`) — verified on the artefacts, zero `2.3`/`2.4` anywhere in a pair-η context.
The `expo` restriction `A ≤ 0, p ≥ 1` holds in **all 1453 fitted cells with zero violations**, and
it makes `ε_ΔR ≤ 1` structurally, closing R32 item 0 for fitted cells (the crossx-consumed series
now delivers `ε_ΔR(0) ∈ [0.113, 0.944]` with 1 rejected cell instead of 3).

**Two things the user should know.** (1) The two stated fit requirements are **the same
constraint** for this functional form (`p > 1`), so the "impose the stronger one first, relax to
the weaker if too many fits are rejected" path does not exist — see R33(b). (2) **Knowingly left
STALE** and listed in R33(f): the Step-4 *fits* (their histograms were refilled, so that subtree is
now mixed), the `_pt4bin` variant, the MC closure, the HIJING overlay / `noovl` / PbPb data, and —
per the explicit instruction — the pp24 cross-section, to which ε_ΔR was NOT propagated.

*The original plan, as written before the work, follows.*

---

**2026-09-07 (round 13) — plan. Rerun on the NEW fiducial cut set + the NEW pair-level
`|η^pair| < 2.2`, and CONSTRAIN the Step-3 `expo` fit so it cannot come out flipped.**

*Plan, written before the work (per-step protocol). Scope confirmed with the user
(AskUserQuestion, 2026-09-07): **pp_full only** on the MC side, **pp24 only** on the data side.
The HIJING overlay, `noovl`, PbPb data, Step-4 fits and the MC closure are OUT of scope and are
knowingly left STALE. ε_ΔR is NOT propagated to the pp24 cross-section (user: future work, other
agents in flight).*

**Why (the input change).** `muon_gap_cuts_acceptance.md` F17 (2026-09-07, user) changed the
fiducial cut set and added a pair-level window; the code was edited there but **nothing was
rerun**, so every artefact in this doc is stale. Specifically:
- `ParamsSet::single_mu_fiducial_gap_cuts` `{{-1.20,-1.05},{-0.10,+0.06},{2.30,2.40}}` →
  **`{{-1.30,-1.05},{-0.10,+0.06},{2.20,2.40}}`**.
- NEW `ParamsSet::pair_eta_fiducial_max = 2.2` (symmetric, strict `<`), joined to `kGapLeg` and
  `kGapPair` in `FillMCTrigEffHists.cxx` and to `MCTrigEffPairSelection.h::FiducialGapCut()`.
  **NOT** on `kGapSingle` — Step 1 fills a single-muon map, there is no pair to cut on (this is
  exactly the user's "the |pair η| cut doesn't apply to step 1, but it does to steps 3/4").
- FORCED: `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap` top bin `{2.0,2.3}` →
  **`{2.0,2.2}`** (the `RDFBasedHistFillingData::SetIOPaths` startup throw slaves it to the
  forward window's lower edge).
- USER DECISION: `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` outer bins
  `±(2.0,2.4)` → **`±(2.0,2.2)`**, so the 9 panels exactly tile the surviving region. This is
  what makes the most-forward Step-3 cells `(−2.2,−2.0)` and `(2.0,2.2)`, and the top folded
  cell `2.0 ≤ |η^pair| < 2.2`, as the user asked. **No canonical binning is retyped anywhere:**
  `MakeDrEtaGroups` reads `eta_max` off the source histogram's own axis, so the folded top group
  follows the pair-η vector automatically (`dr_correction_cell_groups.h`); the pair-pT axis is
  untouched (`.claude/CLAUDE.md` §Binnings).

**The fit restriction (user, §3.3 deliverable).** The Step-3 `expo` form is
`f(ΔR) = C + A·exp[−(ΔR/λ)^p]` (`C ≡ 1` in the plateau-corrected mode, free in the `nocorr`
family). Some cells fit "flipped": decreasing with ΔR, or rising concavely with the turning point
below 0. Mapping the user's two requirements onto the parameters:
1. *"the small-ΔR plateau must lie BELOW the large-ΔR plateau"* — `f(0) = C + A`, `f(∞) = C`, so
   this is exactly **`A < 0`**. It also makes `f` monotonically increasing
   (`f' = −A(p/λ)(ΔR/λ)^{p−1}e^{−u} > 0`), which kills the "decreasing with ΔR" failure too.
2. *"the turning point must be > 0"*, with *"in principle the slope at ΔR = 0 should be zero"*
   tried first — for this form the two are the SAME condition. `f'' = 0` at
   `ΔR_infl = λ·((p−1)/p)^{1/p}`, which is real and positive **iff `p > 1`**; and
   `f'(0) = 0 ⟺ p > 1` as well (`p = 1` gives slope `|A|/λ`, `p < 1` gives an infinite slope and
   no positive inflection — precisely the "concave rise, turning point < 0" failure). So the
   stronger requirement and its fallback collapse to one constraint, **`p ≥ 1`**, and there is no
   weaker variant to fall back to if too many cells are rejected. Implemented as the parameter
   limit, so MINUIT searches only the physical region rather than the fit being screened after
   the fact.
   *Why a limit and not a post-hoc `fit_ok` screen:* a screened-out cell delivers NO correction
   (it falls through to the backup method / raw-bin placeholder); a constrained fit delivers the
   best fit that is physically admissible. The user asked to "restrict the fitting parameters".
3. **Step 4 is deliberately EXCLUDED.** Its physics is the opposite sign: R4/§3.4 measured the
   single-leg efficiency to be *enhanced* at small ΔR (`ε(ΔR<0.12)/ε(ΔR>1) ≈ 1.20 pp`), i.e.
   `A > 0`. Imposing `A < 0` there would be wrong. The restriction is gated on `step == 3`.

**Execution order (data first — the MC Step-1 canvases overlay the data turn-on, and the coarse
q·η binning moved on BOTH sides):**
1. `fit_dr_corrections.cxx` — the Step-3 `expo` limits; `/review-analysis-code`.
2. **Data pp24** — `SKIP_CONDOR=1 pipeline_pp_trig_eff.sh` stages 5–8, Tight + Medium.
3. **MC pp_full** — `SAMPLES="pp_full" run_mc_trigeff_round7.sh` (Steps 1–4 refill, turn-on
   refit, sanity, replot, 2D maps, statistics tables), both WPs.
4. **Step-3 ΔR corrections** — `SAMPLES="pp_full" STEPS=3 run_dr_correction_fits.sh`, no skips,
   all methods × signs × plateau modes.
5. Verify the pair-η cell edges on the produced artefacts; compare the fit reports against the
   round-11/12 values (how many cells the restriction moved / rejected); `/review-plot`.

*Nothing else in this doc changes: no sample, no weight, no fit method, no plateau window, no
pair-pT binning.*

---

**2026-08-18 (round 12) — ✅ DONE. The third Step-3 plateau mode is built, produced, reviewed and
committed; results and numbers in R32.** All 11 plan steps below are complete: code + full
production (pp24 + overlay × both WPs × 3 methods × 3 sign series, step 3), the byte-identical
regression check on the two pre-existing modes, `/review-analysis-code` and the plot review (four
defects found, three of them pre-existing, all fixed and re-verified; APPROVED at iteration 2),
docs and INDEX, commits `22b8454`, `b853b81`, `6841045`, `6c11caa`. **The artefacts are FINAL and
usable — the hand-off block above is live.**

**Carried forward as OPEN, needing a user decision (both are R26; see R32 for the numbers):** a
positive-amplitude / χ² screen for `fit_ok`, and the minimum-points threshold that admits one
unconstrained overlay fit. Neither was changed here: each would alter the published acceptance rule
for all three plateau modes and therefore the MC-closure thread's inputs.

> ### 📌 HAND-OFF TO THE pp24-CROSSX THREAD (`pp24_crossx_rerun_2026_08.md`) — READ THIS FIRST
>
> **The 7-pair-pT-cell ΔR correction your Objective item 3 asks for EXISTS and is usable as of
> 2026-08-17 23:29.** It is the new `nocorr_ptmerge` plateau mode (R32). Nothing about the
> definition, the fit domain or the applied form changed — only the pair-pT cell grid.
>
> **What to load** (named in code so you never retype it —
> `plotting_codes/trig_effcy/mc_based/dr_correction_sample_cfg.h`):
>
> | | value | accessor |
> |---|---|---|
> | method | `expo` | `DrCorrCrossxMethod()` |
> | sign series | `os` (opposite sign) | `DrCorrCrossxSign()` |
> | plateau mode | `nocorr_ptmerge` | `DrCorrCrossxMode()` |
>
> **Files** (pp24 Pythia fullsim FULL, `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/`),
> written **2026-08-18 09:46–09:53**, all AFTER every input they consume. (They were first produced
> 2026-08-17 23:16–23:24 and REFITTED on 2026-08-18 for a report-TEXT correction only — every fitted
> number is unchanged and was re-verified against the pre-refit values.)
> - Tight  : `dr_correction_fits_pp24_full_step3_expo_os_nocorr_ptmerge.root`
> - Medium : `dr_correction_fits_pp24_full_medium_wp_step3_expo_os_nocorr_ptmerge.root`
> - the `polyu_fixedRp` and `interp` siblings exist under the same naming (swap the method token),
>   both signs and both WPs — 18 files in total.
> Do NOT build the path by hand: `DrCorrFitFile(cfg, use_tight_wp, 3, method, sign, mode)`.
>
> **How to apply it — use the existing consumer, do not re-derive.**
> `plotting_codes/trig_effcy/mc_based/dr_correction_apply.h` (`DrCorrectionEvaluator`) already
> implements the agreed form, `ε_ΔR(ΔR) = f(ΔR)/C` for ΔR < 1 and `= 1` for ΔR ≥ 1
> (`mc_trig_eff_closure.md` §3.2), with the baseline screen on C, the TF1 range clamp, the
> floor/cap counters and the per-cell census printout. `Eval(dR, pair_pt, pair_eta)` looks the cell
> up on the FIT FILE's own axes, so the 7-cell grid needs nothing from you.
> **Call it as** `ev.Load(cfg, use_tight_wp, DrCorrCrossxMethod(), DrCorrCrossxSign(), DrCorrCrossxMode());`
> — the `plateau_mode` argument is NEW (2026-08-17) and **defaults to the un-merged `nocorr`**, which
> is deliberate: the MC-closure thread consumes that default and must not move. Pass the mode
> explicitly or you will silently get the 8-cell correction.
>
> **⚠ FOUR THINGS YOU MUST KNOW**
> 0. **One forward cell of the top pair-pT row delivers ε_ΔR = 1.168 at ΔR → 0**
>    (`η^pair ∈ [−2.4,−2.0)`, from the raw-bin placeholder — its `expo`/`polyu` fits are rejected by
>    the baseline screen). A close-by 2mu4 correction is a LOSS, so a value above 1 is an artefact,
>    not a measurement; it is ~1 % of pairs but it is in your weight. It belongs to OPEN R26 and is
>    **escalated to the user** (R32) — do not paper over it, and state it wherever the crossx result
>    is quoted until the user decides.
> 1. **The `expo → polyu → raw` cascade is NOT in `dr_correction_apply.h`** — that class does ONE
>    method → raw-bin placeholder. You have already built the cascade on top of the accessors above
>    in `Analysis/Utilities/DrCorrectionCrossxEvaluator.h` (seen here 2026-08-17, with its own
>    `DrCorrCrossxBackupMethod() = "polyu_fixedRp"`), so nothing is missing — but note that
>    `dr_correction_apply.h` itself was **edited on 2026-08-17** by this thread (the new
>    `plateau_mode` argument, the grouped raw-bin fallback), so re-read it before changing it.
>    Sizing for your routing counters, from the census in R32 (pp24 Tight, opposite sign, merged):
>    `expo` delivers a fit in 60 of 63 cells with 3 on the raw-bin placeholder; `polyu_fixedRp`
>    delivers 62 of 63 with 1. **No cell is left without a correction in either** — the 2–3 dead
>    cells of the un-merged mode are gone, which is most of what the merge bought you. Expect your
>    route census to come out roughly 60 primary / ~2 backup / ~1 raw, not 63/0/0.
> 2. **The raw-bin placeholder is still a flagged TEMPORARY stop-gap** (`mc_trig_eff_closure.md`
>    §3.3, this doc's OPEN R26). It is fine to run on, but it must be stated wherever the result is.
> 3. **The merged mode is defined for the 8-bin nominal pair-pT axis ONLY.** With
>    `MCTRIGEFF_PAIRPT_4BIN` set it throws by design.
>
> **Sanity numbers to check yourself against** (R32). Over the cells delivered from a FIT, ε_ΔR
> spans **0.061 – 1.318** (Tight `expo`, opposite sign) — note this extremum, printed by
> `DrCorrectionEvaluator`, covers the FITTED cells only, **not** the raw-bin placeholder ones.
> Per-cell values at ΔR → 0 in the merged top row `p_T^pair ∈ [72.1,150)`, opposite sign, `expo`:
> η^pair bins 2–6, 8, 9 deliver 0.116 – 0.417 from a fit (a close-by LOSS, as the physics requires);
> bin 7 `[1.0,1.5)` delivers 0.436 from the raw-bin placeholder; **bin 1 `[−2.4,−2.0)` delivers
> 1.168 from the raw-bin placeholder** — above 1, i.e. an artefact of that forward cell, not a
> measurement (in the un-merged mode that cell delivered no correction at all, ε_ΔR ≡ 1).
> A 2mu4 close-by correction is a LOSS, so treat any delivered value above 1 as an artefact and say
> so wherever the result is quoted; the un-merged `polyu` maximum of 2.74 flagged in
> `mc_trig_eff_closure.md` R2 does not survive the merge (merged max 1.24).
>
> **Concurrency.** This thread wrote only: the new `no_plateau_correction_last2ptbins_merged/` plot
> trees, the `*_nocorr_ptmerge.root` fit files, and the source files listed in R32
> (`dr_correction_*.h`, `fit_dr_corrections.cxx`, `plot_dr_correction_fits.cxx`,
> `plot_mc_trig_eff.cxx`, `run_dr_correction_fits.sh`). It did **not** touch
> `RDFBasedHistFilling*`, `MCTrigEffPairSelection.h`, the crossx code, or any pre-existing
> trigger-efficiency output. We share the ACLiC build directory
> `plotting_codes/trig_effcy/mc_based/*.so` — recompile before you run.

*Plan, written before the work (per-step protocol).*

**What the user asked for.** Beside `plateau_corrected/` and `no_plateau_correction/`, add a third
top-level directory under `step3_dr_fit/` carrying the **no-plateau-correction fit with the last
two pair-pT bins combined**. It applies to the **8-bin pair-pT version ONLY** (confirmed with the
user 2026-08-17; the `_pt4bin` variant already has a single top cell). The variant with
**opposite-sign pairs + the `expo` form** becomes the **temporary default the pp24 crossx
application will use**; the other methods (`polyu_fixedRp`, `interp`) and the other sign series are
still produced, with every other setting unchanged.

**Why it is physics-motivated.** On the 8-bin log axis the top two cells — `p_T^pair ∈ [72.1,104)`
and `[104,150)` GeV — run past where pp Pythia has yield: they hold the plateau-guard failures
(R24/R27), they are where the closure collapsed (`mc_trig_eff_closure.md` R1b), and their ΔR fits
are noise-dominated. Merging them buys back statistics in exactly one place without touching the
cells that carry the measurement.

**Scope (confirmed with the user 2026-08-17): PRODUCE THE VARIANT ONLY.** The crossx chain is NOT
touched — ε_ΔR is not yet applied anywhere in the crossx hist filling (only the MC closure consumes
it), so wiring it in is a separate step with its own blast radius. Nothing in
`signal_selection_change_impact.md` is triggered.

**Binning statement (.claude/CLAUDE.md §Binnings).** No canonical binning is changed:
`ParamsSet::pair_pt_coarse_bins` and the FILLED histograms keep their 8 bins. The merge is a
**grouping applied at the fit/plot stage** — bins 7 and 8 are projected together — which is
numerically identical to filling a 7-bin axis, and it is **opt-in and suffixed**
(`_nocorr_ptmerge` / `no_plateau_correction_last2ptbins_merged/`), never a silent second default.
The grouping is derived in ONE place from the source histogram's own axis, so no edge value is
retyped anywhere.

**Implementation plan (per §3.3):**
1. `dr_correction_sample_cfg.h` — third plateau-mode token `nocorr_ptmerge` + the two predicates
   (`no plateau` / `merge last two pT bins`) every stage asks instead of comparing strings;
   directory + file-name token; the named constants recording the pp24-crossx default
   (method `expo`, sign `os`, mode `nocorr_ptmerge`).
2. NEW `dr_correction_pt_groups.h` — the pair-pT grouping (group → source-bin range, group edges,
   merged-axis map booking), built from the source axis, with a hard throw if the merge is asked
   for on the 4-bin variant.
3. `dr_correction_ratio.h` — `DrCellRatio` gains an explicit-bin-range form; the existing
   signature stays a thin wrapper so every current call is byte-identical.
4. NEW `dr_correction_plateau.h` — `PlateauFromRatio` MOVED out of `plot_mc_trig_eff.cxx` so the
   merged cell's plateau is measured by the SAME code that measures every other cell's (it is
   reported-only in a no-plateau-correction mode, but a second copy of a weighted mean is how two
   plateau definitions start to drift — R17).
5. `fit_dr_corrections.cxx` — mode plumbing via the predicates; in the merged mode the plateau /
   nbins / window-systematic maps are re-measured on the MERGED cells from the full-range 3D
   histograms and the fit runs on the 7-cell axis; report + provenance state the merge.
6. `plot_dr_correction_fits.cxx` — same predicates; the measured points are re-projected with the
   grouped bin range; a guard that the fit file's pair-pT axis matches the mode.
7. `dr_correction_apply.h` — optional `plateau_mode` argument (**default unchanged**, so the
   closure thread's results do not move) + grouped raw-bin fallback; header documents the
   temporary pp24-crossx default.
8. `pipelines/run_dr_correction_fits.sh` — the third mode in the defaults and in the three mirror
   functions; skipped with a printed note when `MCTRIGEFF_PAIRPT_4BIN` is set.
9. Run: `pp_full` + `overlay`, both WPs, step 3, all methods, all sign series, mode
   `nocorr_ptmerge` (`SKIP_MEASURE=1` — the plateau files are the round-11 ones and are current).
10. Regression check: re-run one `nocorr` and one `corr` series and diff the fit reports against
    the pre-change copies — the shared-header refactor must not move a single number.
11. `/review-analysis-code` + `/review-plot`, doc + INDEX update, commits.

---


**2026-08-13 (side task) — ✅ DONE. PbPb forward q·η edge scan re-cut as DATA-ONLY (R31):
`single_mu_eff_forward_qeta_edge_scan.png` deleted, replaced by `…_by_year.png` (μ⁺/μ⁻ ×
2023/2024/2025, centrality-integrated) and `…_centrality.png` (23+24+25, μ⁺+μ⁻, one panel per
`ParamsSet::ctrbins` bin). The adopted forward edge **2.30 is confirmed in PbPb**: relative to
2.20 it costs |Δε| ≤ 0.0116 (by year, at 2.30) for **+44.0–46.9 %** probes and ≤ 0.0066 (by
centrality, at 2.28) for **+36.1–37.9 %**, while 2.40 costs **0.038–0.064** (by year) /
**0.047–0.058** (by centrality) — charge-blind (≤ 0.0135) and centrality-blind.
pp untouched; no other PbPb trigger-efficiency plot regenerated. One stated caveat: the 2024/2025
`_fine_q_eta_bin` inputs predate the Δp/p fix `71fcf1c` (see R31).**

**2026-08-13 (round 11) — ✅ DONE, 03:2x. The gap-window rerun is COMPLETE. Every MC and data
trigger-efficiency product in this doc is current with `single_mu_fiducial_gap_cuts =
{{-1.20,-1.05}, {-0.10,+0.06}, {2.30,2.40}}`.**

> ### 📌 HAND-OFF TO THE MC-CLOSURE THREAD (`mc_trig_eff_closure.md`) — READ THIS FIRST
>
> **The rerun is finished; the inputs you consume are final as of 2026-08-13 03:2x.** Use them
> as they now stand:
>
> | what you read | path (pp24 FULL) | written |
> |---|---|---|
> | MC trig-eff hists, Step 3 / Step 4 | `pythia_fullsim_full_sample/mc_trig_eff_hists_pp24_full{,_medium_wp}_step{3,4}.root` | 02:51 / 02:54 |
> | MC single-muon turn-on fits | `…/single_mu_effcy_pT_fit_mc{,_medium_wp}.root` | 02:50 |
> | ΔR plateaus | `…/dr_correction_plateaus_pp24_full{,_medium_wp}.root` | 02:57 / 02:59 |
> | **ΔR-correction fits (54 nominal files)** | `…/dr_correction_fits_pp24_full*<step><method><sign><pmode>.root` | **03:13–03:2x** |
> | data T&P reference, Tight / Medium | `dimuon_data/pp_2024/…_qeta_fid{,_medium_wp}.root` + `…/single_mu_effcy_pT_fit{,_medium_wp}.root` | 02:17 / 02:39 |
> | data T&P, PbPb 23/24/25 Tight; PbPb23 Medium | `dimuon_data/pbpb_20{23,24,25}/…_qeta_fid*.root` + their fits | 02:26–02:38 / 02:44 |
>
> **⚠ Your outputs from before ~03:25 are STALE — please regenerate.** Our two runs overlapped:
> `mc_trig_eff_closure_pp24_full.root` (02:53:04) and `…_medium_wp.root` (02:58:04) were written
> **while** this rerun was in flight, so they combine the new Step-3/4 hists with the **old
> (2026-08-11, pre-widening) ΔR-correction fits** — the fits were only rewritten at 03:13–03:2x.
> A closure built on that mixture is not a closure of anything: the correction applied and the
> efficiency it was derived from used different fiducial regions.
>
> **Concurrency, going forward.** We share `pythia_fullsim_full_sample/` and the RDF/ACLiC
> build directory, and none of it is under git, so a collision is silent. This thread is now
> **idle** on those paths — nothing further will be written here without saying so in this doc
> first. Your `_pt4bin` Step-3 refill (02:59:56) did not collide with anything of ours (the
> nominal token is empty), but note `_pt4bin` remains STALE for everything else (**R24b**) — it
> was deliberately NOT part of this rerun.

*Plan, written before the work (per-step protocol).*

**Why.** `ParamsSet::single_mu_fiducial_gap_cuts` central crack window was widened
`{-0.06, +0.06}` → `{-0.10, +0.06}` (commit `1f143e6`, user decision recorded in
`muon_gap_cuts_acceptance.md` F11a; the measured depletion is one-sided, N(neg)/N(mirror pos) =
0.509 pp / 0.578 PbPb in (−0.10, 0.00) and ≈1 everywhere else, so the symmetric window was
cutting healthy acceptance on the positive side and leaving the depleted slice in). That vector is
LIVE in the trigger-efficiency chain — `FillMCTrigEffHists.cxx:641–663` (MC, Steps 1–4, single /
leg / pair forms) and `RDFBasedHistFillingPP.cxx:172` + `RDFBasedHistFillingPbPb.cxx:342` (data
tag-and-probe, PROBE leg only) — and the windows are built into the RDF/JIT selection strings at
FILL time, so the numbers are baked into the filled histograms: a replot cannot pick the change up,
the hists must be REFILLED. F11b flagged this and left the rerun to this doc's session; the user
has now asked for it, explicitly including the data single-muon mu4 efficiency.

**Blast radius (stated before running).** The gap cut is NOT yet in the signal selection (deferred,
`muon_gap_cuts_acceptance.md` F10.3), so **crossx and R_AA are NOT affected** and
`signal_selection_change_impact.md` is not triggered. Stale, and only: the data single-muon
efficiency (pp24, PbPb 23/24/25) and everything downstream of it, and the whole MC trig-eff chain.
The NTuple processing is untouched — the cut lives at the RDF stage (D6) — so no condor, no
re-skim. `CommonEffcyConfig` consistency is unaffected: only the FORWARD window's edge (2.30) is
coupled to the coarse q·η binning; the central window sits inside existing bin edges.

**Order (data first, because the MC Step-1 canvases overlay the data turn-on):**
1. **Data, pp24** — `SKIP_CONDOR=1 pipeline_pp_trig_eff.sh` (stages 5–8: RDF fill → turn-on fit →
   validate → Pipeline-2 plots).
2. **Data, PbPb 23/24/25** — `SKIP_CONDOR=1 pipeline_pbpb_trig_eff.sh` (stages 5–10; 9–10 are the
   inverse-weighted data ΔR corrections, which consume the refitted ε_single and therefore move
   too).
3. **MC chain** — `run_mc_trigeff_round7.sh` (refill Steps 1–4, refit turn-ons, replot) for
   `pp_full` + the two test samples, both working points.
4. **ΔR corrections** — `run_dr_correction_fits.sh` with NO skips (measure → guard+fit → plot),
   which also re-lays the round-10 two-view tree on the new numbers.
5. Re-verify: plateau guard verdict, fit-report rejection counts, and the Step-1 data/MC comparison,
   against the pre-change values recorded in R24–R28.

*Nothing else in this doc changes: no binning, no sample, no weight, no fit method, no plateau
window.*

---


**2026-08-11 (round 10) — IN PROGRESS. Two user follow-ups on the round-9 Step-3 outputs.**

*Plan, written before the work (per-step protocol). pp only; Step 4 is explicitly out of scope
for the new variant.*

1. **One y range per ΔR-range subdirectory** (not per PNG). ✅ **DONE** — `plot_mc_trig_eff.cxx`
   now computes the range once over every (pair pT, pair η) cell of a view and applies it to all
   eight pair-pT canvases; it is printed to the log. Tight: `dr_0_to_1` [0.039, 1.815],
   `dr_0_to_2` [0.000, 2.545], `dr_full_range` [0.000, 2.584]. Regenerated pp × both WPs,
   overwriting in place. *Motivation:* with a per-file range, flipping through the eight files
   compared curves drawn on different axes — a cell could look deeper than one in another file
   purely from the scaling. The Step-3 FIT canvases already shared a range across their pair-pT
   PNGs (`compute_range` in `plot_dr_correction_fits.cxx`), verified, so no change was needed
   there; the standalone inclusive canvas keeps its own range (round-9 review finding: the shared
   one squeezed its data into <10 % of the frame).

2. **A no-plateau-correction Step-3 fit variant** — ✅ **DONE (see R27 for the results).**
   *User motivation, recorded verbatim:* the inverse-weighted ΔR distribution has "complicated
   structures ... even in the full range, especially in the gap-enclosing pair-η bins, which means
   the plateau determined in a large ΔR region might not be accurate for small ΔR (region of
   interest)". This is the same pathology R23 and R26 circle: in `pT_pair [104, 150)` the per-cell
   plateaus sit at **1.15–1.32** while the small-ΔR points sit at **0.2–0.5**, so a far-region
   normalization moves the region of interest by 20–30 %.
   *Two user decisions (AskUserQuestion, 2026-08-11):*
   - **(D-R1) The no-correction fit floats its own asymptote:** `f(ΔR) = C + A·exp[−(ΔR/λ)^p]`
     with **C free**, fitted to the **raw, un-normalized** ε_ΔR. The baseline is then determined by
     the ΔR < 1 data itself and the [2, 3.5] window **never enters the fit**. `polyu_fixedRp` gets
     the same treatment (its leading 1 becomes a free constant); `interp` pins its flat region to
     the last measured knot instead of to 1.
   - **(D-R2) Layout:** the plateau mode is the TOP level —
     `step3_dr_fit/{plateau_corrected,no_plateau_correction}/<method>/{sign_intgr,sign_sepr}/`, so
     each mode is a complete parallel tree comparable as a unit (which is what the future MC-closure
     comparison needs). Step 4 stays corrected-only per the user, but its tree moves to
     `step4_dr_fit/plateau_corrected/...` so the two steps do not end up with different shapes.
   *Follow-up user requirement:* **C must be printed on every subplot**, explicitly labelled as the
   fitted plateau, per sign in `sign_sepr`, with the `(at limit)` marker if it hits a bound — and it
   must never be confusable with the `plateau = …` line the corrected mode prints from the [2, 3.5]
   window, which the no-correction mode does not use.
   *Consequence settled from the stated purpose (not a free choice):* **the plateau guard does not
   apply in no-correction mode.** Nothing is normalized by the plateau, so `|plateau − 1| > 0.15`
   cannot disqualify a cell and a cell whose far-ΔR plateau is unmeasurable is still fittable. The
   variant should therefore cover MORE cells than the corrected one — itself a useful number.
   *Fit domain, confirmed in code and reported to the user:* ΔR ∈ **[0, 1]** only — `kFitLo/kFitHi`,
   the `zoom` histogram alone (20 bins of 0.05, centres 0.025…0.975), TF1 built on that range and
   fitted with `"QRNS"`. Points drawn between ΔR = 1 and 2 come from the wide-bin histogram and are
   excluded. So C is fixed by the ΔR ≈ 0.4–1.0 points; where a cell has not flattened by ΔR = 1, C
   absorbs the residual slope and becomes correlated with A — visible as a large error on C.

3. **`sign_sepr` marker/colour scheme** (user, 2026-08-12). ✅ **DONE — see R28.** Same sign =
   blue triangles (`kBlue+1` markers and error bars) on a `kBlue` fit; opposite sign = red circles
   (`kRed+1`) on a `kRed` fit. All 12 `sign_sepr` directories regenerated (plot stage only; no fit
   or measurement was re-run). R28 also records, in full, the fit-rejection standard the user asked
   for — the five `fit_ok` criteria, what `interp` uses, and the fact that χ²/ndf is NOT among them
   (R26).

*Still open and unchanged:* R25 (sign-dependent correction — user decision, crossx blast radius),
R26 (both parametric forms fail in a minority of cells), R23, and the earlier carry-overs.

---


**2026-08-11 (round 9) — DONE. Both reviews PASS. Two OPEN findings handed to the user.**

Delivered, pp only (`pp_full`, both working points): the Step-1 2D (q·η, pT) single-muon efficiency
maps; `step3_dr_correction/` restructured into three ΔR ranges × one PNG per pair-pT bin × one
subplot per pair-η bin with the per-cell plateau drawn; same-sign / opposite-sign separated Step-3
and Step-4 fits under `sign_intgr/` + `sign_sepr/` (with per-sign ratio canvases); and the MC
pair-count / cross-section CSV tables on the 8-bin pair-pT axis. PbPb code updated and its plot
stage re-run for consistency (plateaus verified bit-identical, 288/288 cells) but **not refilled**,
so it has no per-sign histograms — every per-sign path degrades to a printed note, verified.
**r17663 (`noovl`) untouched throughout.**

`/review-plot` **APPROVED at iteration 4**; `/review-analysis-code` **PASS at iteration 2**. Every
numerical claim in R24–R26 was independently reproduced by the reviewers. Between them the two
loops found and fixed, in the round-9 work: canvases drawing fits the analysis' own `fit_ok` screen
rejects (32/72 same-sign Step-3 `expo` cells) with full parameters; measured points silently dropped
above the y cap (15 hidden on the overlay, one at 11.19 in a panel that IS drawn); a legend
overrunning the canvas; inclusive canvases squeezed into <10 % of their frame; a `#kern` with empty
braces that shifted nothing; and the raw-probability figure drawn with the sign colours inverted.
Pre-existing bugs fixed along the way: a missing pair-pT file-suffix in the fit plotter; a
**`root_file_has_objects()` that was a silent no-op**, so every artefact check in that pipeline had
been passing unconditionally; a stale PNG-count literal; and `n_unusable` disagreeing with the
`fit_ok` map it described.

**★ TWO THINGS WAIT ON THE USER — neither is acted on:**
1. **R25** — the Step-3 (both-legs) ΔR correction is **not charge-blind**: same-sign/opposite-sign
   = 0.43 at ΔR = 0.025, dying out by ΔR ≈ 0.175, absent from Step 4. Proven to be in the RAW
   trigger probability (`numraw`), which rules out both an ε_MC-weight artefact and a kinematic-mix
   artefact. The sign-integrated nominal ≈ the opposite-sign curve (87 % of pairs), so it
   under-corrects same-sign pairs by up to ×2.3 at small ΔR — and the background subtraction is
   OS − SS. Applying the sign-dependent correction needs a crossx refill; coverage cost and the
   TRUTH-vs-RECO charge caveat are tabulated in R25.
2. **R26** — **both** parametric fit forms fail badly in a substantial minority of cells (7–19 of
   ~71 above χ²/ndf 5; the densest cell reaches 48.2 with a dip-then-overshoot no monotone form can
   represent), while `usable` has no χ² term so all of them are published with `fit_ok = 1`. Most
   likely resolution is a χ² screen in `usable` or per-cell interpolation — **not** the
   `polyu_fixedRp` swap an earlier, wrongly-parsed version of R26 recommended (withdrawn).

Also open and unchanged from earlier rounds: R23 (top-pair-pT plateaus 20–42 % high), the §2
plateau-normalization before crossx, the trigger-group send, and the full data-side dp/p cascade.

---

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
