# Adoption of muon pT > 4.5 GeV, gap window (−1.25,−1.05), pair pT > 9 GeV — full rerun

Mode: **IMPLEMENTATION**. Opened 2026-09-08.

## Objective

Adopt three user-decided selection changes across the whole analysis, data **and** MC, and
rerun every stage they invalidate:

1. **Muon pT threshold 4.0 → 4.5 GeV** — reconstructed pT in data, truth pT in MC, and
   reconstructed pT inside the trigger-efficiency measurement. This executes the decision
   that `muon_pt45_cut_diagnostic.md` was PARKED awaiting (that doc is hereby unblocked and
   its temporary diagnostic code is to be deleted).
2. **Fiducial gap window `(−1.30, −1.05)` → `(−1.25, −1.05)`** in
   `ParamsSet::single_mu_fiducial_gap_cuts`. The other two windows
   (`(−0.10,+0.06)`, `(2.20,2.40)`) and the pair-level `|η^pair| < 2.2` are unchanged.
3. **Signal-region pair pT `> 8` → `> 9` GeV**, in every signal selection — data and every
   truth/reco MC analog.

Together these also carry two user-decided **binning** moves and a full **Pb+Pb migration**
(see Design Decisions).

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. **Cuts** in place at every site: data reco `pt > 4.5`; MC truth `truth_pt > 4.5`;
     signal-region `pair_pt > 9` in data and in every truth/reco analog; gap window
     `(−1.25,−1.05)`.
  2. **Trigger efficiency** selects on reco `pt > 4.5` AND truth `truth_pt > 4.5` AND the
     FULL nominal reconstructed-muon cut set (quality/WP, |η|, one-sided Δp/p, |d0|,
     |z0 sinθ|, trigger matching as applicable) — every gap against the nominal NTuple
     selection either closed or explicitly recorded as unclosable with the reason.
  3. **Binnings** moved (user decisions D2/D3): `ParamsSet::pTbins` and
     `single_mu_pt_coarse_bins` low edge 4 → 4.5; turn-on fit range `[4,60]` → `[4.5,60]`
     with the log/linear pivot at 4.5 (data AND MC fitters); `pair_pt_coarse_bins`,
     `pT_bins_120`, `pair_pt_coarse_bins_pt150` and the 4-bin variant low edge 8 → 9,
     same bin counts, same log-spacing rule.
  4. **Pb+Pb migrated** (user decision D1): `RDFBasedHistFillingPbPb` and the overlay
     per-centrality nodes adopt `single_mu_fiducial_gap_cuts` + `pair_eta_fiducial_max`,
     and `ParamsSet::N_PAIR_ETA_CROSSX_BINS` (48) replaces the retyped `44, −2.4, 2.4`;
     the `EvaluateSingleMuonEffcyPtFitted` `−1.0f` sentinel bug of
     `pp_trig_eff_highpt_jump.md` is resolved as a by-product and that doc closed.
  5. **NTuple Condor rerun**: pp24 both modes, Pb+Pb 2023/2024/2025, Pythia fullsim pp24
     (full `_pdf`), HIJING overlay, Pythia truth, POWHEG truth (+ POWHEG fullsim if any
     live product needs it).
  6. **Derived corrections rebuilt**: data mu4 turn-on refits (pp + Pb+Pb, Tight AND
     Medium), MC trig-eff Steps 1–4, ΔR-correction fits, MC closure,
     `pair_reco_eff_pp24_full.root`.
  7. **Final results regenerated**: pp24 crossx, Pb+Pb crossx (23+24+25 combined), R_AA,
     MC-vs-data comparison, signal acceptance + cutflow, crossx sanity/stage plots.
  8. Temporary `_diag_mupt45` / `_logbins_from9` diagnostic code DELETED
     (`muon_pt45_cut_diagnostic.md` Remaining Work).
  9. Docs updated — this doc, `analysis_overview.md` §2, `signal_selection_change_impact.md`,
     `muon_gap_cuts_acceptance.md`, `muon_pt45_cut_diagnostic.md`,
     `pp24_all_vertex_pairs.md`, `pp_trig_eff_highpt_jump.md`, `ami_weights.md` if any
     sample moves, `INDEX.md` — and the work committed (one commit per logical change).
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

**(a) Muon pT > 4.5 GeV.** `muon_gap_cuts_acceptance.md` F14 measured the single-muon q·η
shape in four pT slices and found that muons with reconstructed pT ∈ [4, 4.5) GeV carry
deep, narrow depletions at q·η ≈ +0.6 and +1.1–1.3 that shrink monotonically with pT and are
gone above 6 GeV — i.e. exactly where the adopted fiducial windows have no coverage. Their
reconstruction and trigger efficiencies are therefore both the smallest and the most
rapidly varying in q·η, so an efficiency correction applied to them carries the largest
systematic per corrected pair. Raising the threshold removes that population instead of
modelling it. The statistics price was measured in `muon_pt45_cut_diagnostic.md`:
**pp24 −23.25 %, Pb+Pb 23+24+25 −34.53 %** of raw signal-region OS pairs. The user has
decided that price is worth paying.

**(b) Gap window (−1.25, −1.05).** The barrel/endcap dip measured in F6 has its minimum in
[−1.12, −1.10) at 32 % of plateau and half-depth points at ≈ −1.17 and ≈ −1.08. The window
adopted on 2026-09-07 ran to −1.30, i.e. ≈0.13 beyond the measured half-depth on the
negative side, where the yield is only 7–33 % depleted. Pulling the lower edge back to
−1.25 recovers that largely-good acceptance while still covering the dip. Same physics
argument as F6/F11a, applied to the outer edge.

**(c) Pair pT > 9 GeV.** Raising the single-muon threshold to 4.5 GeV raises the kinematic
floor of the pair: two muons at 4.5 GeV cannot make an arbitrarily soft pair, and the region
just above 8 GeV becomes a sculpted corner of phase space whose acceptance is set by the
single-muon cut rather than by a smooth detector response. Moving the signal cut to 9 GeV
removes that corner rather than correcting it.

### 2. Top-level equation

Unchanged in form (`pp24_crossx_rerun_2026_08.md` §2). What changes is the domain of the
sum and the definition of the per-pair weight's inputs:

    dsigma/dX = (1 / L_int) * SUM_{pairs in signal region} w_reco(pair) * w_trig(pair)

Signal region, **pp24 (reco), after this change**:

    minv in (1.08, 2.9)  &&  pair_pt > 9
      &&  BOTH muons: pt > 4.5 && ParamsSet::PassSingleMuFiducialGap(eta, charge)
      &&  |pair_eta| < ParamsSet::pair_eta_fiducial_max

with `single_mu_fiducial_gap_cuts = {(−1.25,−1.05), (−0.10,+0.06), (2.20,2.40)}`, rejected
on CLOSED intervals, built in FLOAT. **Pb+Pb adopts the identical form in this change**
(decision D1), so pp and Pb+Pb share a signal region again and R_AA is formable.

Truth analog: the same with `truth_*` variables plus `from_same_b`, and `truth_pt > 4.5`.

### 3. Step-by-step method

**(a) Where the muon pT cut lives.** It is an **NTuple-processing** cut, not an RDF signal
cut — `DimuonDataAlgCoreT.c` for data, the truth-pT gates in `PythiaAlgCoreT.c` /
`PowhegAlgCoreT.c` and the reco/truth gates in `PythiaFullSimExtras.c` /
`PowhegFullSimExtras.c` for MC. Raising it therefore requires a **full NTuple rerun**
(Condor) in data and MC, and it breaks the §1 boundary of
`signal_selection_change_impact.md` exactly as `pp24_all_vertex_pairs.md` did.

**(b) Where the pair-pT cut lives.** At the **RDF** stage, in every signal-region filter
(`signal_cuts` and the truth analogs). It needs no NTuple rerun of its own — but it rides
along with (a).

**(c) Where the gap window lives.** One vector in `ParamsSet.h`, consumed at the RDF stage
via `FiducialGapCutExpr` / `PassSingleMuFiducialGap` at the 13 live sites inventoried in
`muon_gap_cuts_acceptance.md` F17. Changing the value is a one-line edit with the full
downstream blast radius and no NTuple consequence.

**(d) The trigger-efficiency selection must be the nominal reconstructed-muon selection.**
The MC trigger efficiency measures the probability that a *nominally selected* muon pair
fires the trigger. Its muon population must therefore be defined by the SAME reconstructed
cuts the nominal analysis applies — not a reduced subset. Any cut present in nominal NTuple
processing but absent from the trigger-efficiency selection means ε is measured on a
different (looser) population than the one it is later applied to, and the mismatch does not
cancel: it is a population-definition error, not a ratio bias. Both the reco threshold and
the truth fiducial move to 4.5 GeV (user decision D4).

**(e) Binning low edges follow the cuts** (user decisions D2/D3). An axis whose first bin is
partially empty reports a differential quantity diluted against its own bin width while its
label claims the full range. Bin COUNTS and the log-spacing rule are preserved; only the low
edge moves.

**(f) Efficiencies are measured and applied on the same fiducial region.** The Pb+Pb
migration (D1) exists so that the Pb+Pb trigger efficiency (already measured WITH the probe
fiducial cut) and the Pb+Pb signal selection stop describing different regions —
`muon_gap_cuts_acceptance.md` F17 open defect 1 — and so that R_AA stops mixing fiducial
regions (open defect 2).

### 4. Negative constraints (what the code must NOT do)

1. Do **NOT** apply the 4.5 GeV cut only at the RDF stage. That was the *diagnostic's*
   deliberate shortcut (`muon_pt45_cut_diagnostic.md` Design Decisions) to avoid a Condor
   rerun; the adopted change is an NTuple-stage change, in data and in MC.
2. Do **NOT** leave any pT threshold at 4 in a selection that mirrors another one at 4.5.
   The data reco cut, the MC truth cut, the MC reco cut, the trig-eff reco cut and the
   trig-eff truth fiducial must all move together.
3. Do **NOT** re-invent or silently re-derive any binning. Only the low edges named in
   Done item 3 move, and they move in `ParamsSet.h` alone; every consumer reads them.
4. Do **NOT** change the crack window `(−0.10,+0.06)`, the forward window `(2.20,2.40)`,
   the pair-level `|η^pair| < 2.2`, `minv ∈ (1.08,2.9)`, the one-sided Δp/p cut, or the
   all-vertex pp impact-parameter procedure. None of them is in scope.
5. Do **NOT** reuse another production's AMI weights when any MC sample is re-processed
   (`ami_weights.md`, BLOCKING). The productions are unchanged here, so the existing
   `ami_info/` applies — but the `expected_ami_dsids` guard must still pass, not be bypassed.
6. Do **NOT** form R_AA until pp and Pb+Pb are BOTH on the new signal region. With D1 they
   will be, but not until the Pb+Pb crossx refill completes.

## Context

- **Unblocks** `muon_pt45_cut_diagnostic.md` (PARKED since 2026-09-07 pending exactly this
  decision).
- **Subsumes** the unfinished MC half of `pp24_all_vertex_pairs.md` (its steps 6, 7a–c, 8
  were started and never completed; nothing was running when this task opened). The
  all-vertex pp selection itself is DONE and stays.
- **Subsumes** the "nothing has been rerun" state of `muon_gap_cuts_acceptance.md` F17.
- **Closes** `pp_trig_eff_highpt_jump.md` item (2) by its own option (d).
- Blast radius reference: `signal_selection_change_impact.md` (which must itself be amended
  — its §1 NTuple-unchanged boundary does not hold for the pT cut).

## Scope

IN: every data and MC NTuple stage; every RDF hist-filling stage; the data and MC trigger
efficiencies; the ΔR corrections and closure; the pp24 pair reco-efficiency; pp24 and Pb+Pb
crossx; R_AA; MC-vs-data comparison; signal acceptance and cutflow; the Pb+Pb fiducial /
pair-eta-axis migration; all affected docs.

OUT: template fits themselves (their MC inputs are regenerated, but no fit method changes);
unfolding; ε_acc construction (still unbuilt); the F14 open question of extra q·η windows at
+0.6 / +1.1–1.3; the single-value pair trigger efficiency study (not wired into crossx).

## Design Decisions

- **D1 (2026-09-08, user).** Pb+Pb gets the FULL migration, not just an NTuple rerun:
  fiducial + pair-level gap cuts and `N_PAIR_ETA_CROSSX_BINS` (48), then crossx + R_AA.
  Reason: it is the only option that leaves pp and Pb+Pb on one signal region, and it
  resolves two standing blockers as a by-product.
- **D2 (2026-09-08, user).** Single-muon pT axes move to 4.5: `ParamsSet::pTbins`,
  `single_mu_pt_coarse_bins`, and the turn-on fit range/pivot in BOTH fitters.
- **D3 (2026-09-08, user).** All canonical pair-pT axes move to start at 9 GeV, same bin
  counts and same log-spacing rule.
- **D4 (2026-09-08, user).** In the trigger efficiency BOTH the reco threshold and the truth
  fiducial move to 4.5 GeV; and the trigger-efficiency selection must apply the full nominal
  reconstructed-muon cut set, not just WP + pT + |η|.
  **Audit outcome (`_sub_trigeff_muon_cuts_audit.md`): the second half is MOSTLY already
  satisfied and the premise was wrong in one respect** — `wp_col` is not a bare quality bit but
  `pass_tight | pass_medium`, per-muon booleans computed in NTuple processing by
  `PythiaFullSimExtras::PassMuonMediumCuts`, a cut-by-cut mirror of the data's
  `PassCuts_DataCore`. So the trig-eff selection already carries combined, WP, IDCuts,
  MuonCuts, |η|<2.4, pT, one-sided Δp/p, |d0| and |z0 sinθ|; the literal `pt > 4 && |eta| < 2.4`
  in the RDF strings is a redundant restatement. **ONE genuine gap, pp only, undocumented:**
  the pp SAME-VERTEX pair requirement lives only in `pair_pass_{medium,tight}` (via
  `ip_pair_ok`), and Steps 2/3/4 + closure + pair-eff use the PER-MUON flags on a pair, so they
  never apply it. Being closed in this change. Effect expected sub-0.1 % on the ratio; it
  matters as provenance — ε_ΔR must be measured on the population it is applied to. Deliberate,
  documented omissions kept: trigger match (the measurand), resonance veto (inputs are the
  `_no_data_resonance_cuts` trees), photoproduction veto, PbPb event cleaning.
- **D5 (2026-09-08, user).** The default fine pair-pT crossx axis becomes
  `pT_bins_150` = **16 log bins 9→150**, chosen so every coarse edge IS a fine edge
  (coarse k = fine 2k, exactly two fine bins per coarse cell) — the nesting is the reason for
  16 rather than 15. `pT_bins_120` becomes the **opt-in `_pt_120` alternative view**, 16 log
  bins 9→120, which deliberately does NOT nest and must never bin a correction. The old
  `_pt_150` plot directories are cleaned up (the nominal now reaches 150).
  For the record, the previous axes were **15** log bins: 8→120 and 8→150.
- **D7 (2026-09-08, agent, flagged to user).** The truth/reco MC analogs that were still on the
  RETIRED standalone `q·η < 2.2` are migrated onto the fiducial + pair-level windows in the same
  step: `PythiaTruth` (template `kin_cuts` + truth signal acceptance), `PowhegTruth` (signal
  acceptance), `PythiaFullsimOverlay` (per-centrality `pass_signal_{truth,reco}`),
  `PowhegFullsim` (legacy `pass_signal_{truth,reco}`), and the acceptance CUTFLOW macro.
  **Rationale, and why this is a completion of D1 rather than new scope:** D1's stated purpose is
  that pp and Pb+Pb describe ONE signal region. If the truth analogs keep the retired cut, the
  truth acceptance and ε_acc describe a different region from the reco selection they are supposed
  to correct — `muon_gap_cuts_acceptance.md` F17 open defect 4, which says in as many words that
  the truth acceptance "cannot yet supply ε_acc for the current fiducial region". Reversible if
  the user disagrees; nothing else depends on it.
- **D8 (2026-09-08, user).** The DATA muon-pT cut is **mode-dependent**: the trigger-efficiency
  NTuple mode (`trigger_effcy_calc`) keeps **pT > 4.0**, the nominal analysis mode uses **4.5**.
  **Reason:** ε^nc(pT, q·η) is a PER-MUON efficiency, so the population it is MEASURED on need
  not equal the population it is APPLIED to, and it is only ever evaluated above 4.5. Cutting
  probes at 4.5 would delete ~40 % of the 4→6 GeV mu4 rise and leave the fitted turn-on midpoint
  OUTSIDE the fit range in 63–80 % of q·η cells (see R3).
- **D9 (2026-09-08, user).** MC abandons the 4.0–4.5 population **entirely** — no loose flags, no
  second set of branches, no mode dependence. **User's reasoning, checked and agreed:** the MC
  trigger efficiency is NOT tag-and-probe (it is truth-seeded and reco-matched), so there is no
  probe population to widen; MC statistics are ample at low pT and scarce at high pT, so nothing
  is bought; and the binnings start at 4.5 / 9 to match data. Checked for a counter-argument and
  found none that survives: the MC deliverable is ε_ΔR, which is **pair-level and must be measured
  on the analysis population anyway**, and Step 1's singles map is a closure / MC-vs-data
  comparison object rather than the applied correction (the applied one is the data T&P).
  An earlier proposal to loosen `pass_medium`/`pass_tight` to 4.0 was REJECTED by the user on the
  grounds that a future consumer would then inherit an MC/data mismatch BY DEFAULT and only memory
  would prevent it — a correct objection, and the reason the loose population is gone rather than
  renamed.
- **D10 (2026-09-08, user).** `ParamsSet::pT_bins_8` — the low half of the trigger-efficiency
  `pt2nd` MEASUREMENT axis, distinct from the analysis axes — becomes **two segments with 4.5
  forced as an interior edge**: 2 log bins 4.0→4.5 plus 18 log bins 4.5→8.0 (still 20 bins).
  Data fills the whole axis and keeps its turn-on leverage; MC starts exactly ON a bin edge, so it
  has no partially-filled first bin and no x-shift where the turn-on is steepest. Above 4.5 the
  two maps stay bit-identical, so the Step-1 MC/data ratio is untouched.
- **D11 (2026-09-08, user, follows D8/D9/D10).** Fit ranges DIVERGE on purpose: the DATA fitter
  `SingleMuEffcyPtTurnOnFitter` keeps **[4, 60] with the pivot at 4.0**; the MC twin
  `FitMCSinglesEffcy` uses **[4.5, 60] with the pivot at 4.5**. The MC fitter's "mirrors the data
  fitter EXACTLY" header is amended to name this as the one deliberate difference, with the reason,
  so it does not read as drift.
- **D6 (2026-09-08, user).** The FINE q·η axis is left alone. The new window edge −1.25 is not
  a fine-axis edge (`makeEtaTrigEffcyBinning` runs …−1.26, −1.24… in 0.02 steps), so it splits
  a bin. Accepted: the cut is applied to the muon's own q·η value and the correction uses the
  COARSE contiguous bins, so only the shaded band on the fine-axis diagnostics and the 2D
  (q·η, pT) map cell edge are affected — not any corrected number.

## Implementation Plan

1. Blast-radius enumeration by independent subagents; cross-check. (in flight)
2. Audit the trigger-efficiency muon selection against nominal NTuple processing (D4). (in flight)
3. `ParamsSet.h` + `CommonEffcyConfig.h`: gap window, pT axes, pair-pT axes. → `/review-analysis-code`
4. NTuple-processing pT cuts (data + all MC). → `/review-analysis-code`
5. RDF signal selections: pair pT > 9 everywhere; Pb+Pb fiducial + pair-level cut + 48-bin
   pair-eta axis; overlay per-centrality nodes. → `/review-analysis-code`
6. Trigger-efficiency selection (D4). → `/review-analysis-code`
7. Turn-on fitters: range + pivot. → `/review-analysis-code`
8. Delete the temporary `_diag_mupt45` / `_logbins_from9` diagnostic code.
9. Compile everything; small local tests before any Condor submission.
10. Condor: data pp24 (2 modes) + Pb+Pb 23/24/25; MC fullsim pp24, overlay, Pythia truth,
    POWHEG truth.
11. Turn-on refits (pp + Pb+Pb, Tight + Medium).
12. MC trig-eff Steps 1–4; ΔR fits; closure.
13. `pair_reco_eff_pp24_full.root` rebuild.
14. crossx pp + Pb+Pb; R_AA; MC-data comparison; acceptance + cutflow; sanity/stage plots.
    → `/review-plot`
15. Docs + commits.

## Progress Log

- 2026-09-08 Step 0: doc triage (`INDEX.md`, `signal_selection_change_impact.md`,
  `muon_pt45_cut_diagnostic.md`, `muon_gap_cuts_acceptance.md`, `pp24_all_vertex_pairs.md`);
  verified NOTHING was running (0 Condor jobs, no background processes), so the unfinished
  MC half of `pp24_all_vertex_pairs.md` is quiescent and safe to subsume. Four user
  decisions taken (D1–D4). The concurrent thread's uncommitted dr-correction / pair-eff work
  was committed unchanged first (`3a97a3f`) so it cannot be tangled with this rerun.

- 2026-09-08 Steps 3–8 DONE (code edits). Delegated to three parallel subagents on DISJOINT file
  sets (NTuple processing / trig-eff selection / turn-on fitters) plus the orchestrator's own
  work on `ParamsSet.h`, the Pb+Pb migration, the truth analogs and the diagnostic deletion.
  Their scratch docs (`_sub_edit_ntuple_ptcut.md`, `_sub_edit_trigeff_selection.md`,
  `_sub_edit_turnon_fit_range.md`, `_sub_mupt45_rerun_map.md`, `_sub_mupt45_sites_B.md`,
  `_sub_trigeff_muon_cuts_audit.md`) are merged into R1/R2/R3 and the Design Decisions above,
  then deleted. Details worth keeping that are not already recorded elsewhere in this doc:
  * **NTuple thresholds, 7 live sites:** `DimuonDataAlgCoreT.c` (data reco, now mode-dependent
    per D8), `PythiaAlgCoreT.c` + `PowhegAlgCoreT.c` (truth pair cuts),
    `PythiaFullSimExtras.c` ×2 (reco flag + truth-tree gate), `PowhegFullSimExtras.c` ×2.
    A pre-existing convention asymmetry was PRESERVED rather than silently normalised: pair
    cuts are inclusive (`reject < 4.5`, i.e. keep ≥ 4.5) while single-muon truth gates are
    exclusive (`> 4.5`). The `store_mc_trigger` loose gate stays at `pt > 3.0` — deliberately
    below the fiducial so a truth-pT gate cannot sculpt the reco-pT turn-on; it now has 1.5 GeV
    of margin instead of 1.0.
  * **Dead NTuple code left at 4.0 on purpose**, with liveness established by evidence (only
    three entry headers are `.L`-loaded by live run/`.sub` scripts): `PythiaNTupleFirstPass*`,
    `original_no_template_class/*`, `maybe_old_unsure_pbpb/*`.
  * **The trig-eff same-vertex fix REPLACES the per-muon WP terms with a `pair_wp` alias**
    (`pair_pass_tight`/`pair_pass_medium`) rather than ANDing them on. Because
    `pair_pass_X = m1.pass_X && m2.pass_X && ip_pair_ok`, it strictly implies both per-muon
    flags, so ANDing would be redundant AND would mask a future divergence instead of surfacing
    it. Step 1 keeps the per-muon flag — it is a single-muon map with no pair. Measured effect
    on `pp_full`: **795 / 1 698 987 = 0.047 %** of pairs removed, consistent with the audit's
    sub-0.1 % expectation; an exact no-op on the overlay and `noovl` (`ip_pair_ok ≡ true`).
  * **The turn-on pivot is no longer a literal anywhere.** Both fitters now compose it from
    `pT_min`, which makes the pre-existing comment's claim ("the 4.0 pivot IS pT_min")
    structurally true instead of merely asserted. Two further copies were found and fixed: the
    drawn LaTeX equation hardcoded `(p_{T}-4)/4` (the published figure would have misstated the
    function it was drawn from), and `FitMCSinglesEffcy.cxx` had a hardcoded
    `SetLimits(4.0, 60.0)`. Fit parameter inits/limits were deliberately NOT moved — they are
    physics priors on where the mu4 turn-on sits, which is a trigger property independent of
    where the offline analysis cuts.

## Results & Observations

### R1 — Recon findings that change the plan (three independent subagents, 2026-09-08)

Scratch docs `_sub_mupt45_rerun_map.md`, `_sub_mupt45_sites_B.md`,
`_sub_trigeff_muon_cuts_audit.md` (merged here, then deleted).

**Containment — no re-skim needed.** `SkimCode/.../TrigRates.cxx` stores every `xAOD::Muon`
with **no pT selection anywhere**; the two 4-GeV-ish numbers in it are an L1-RoI matching
*flag* (`pt>4000 → match_mu4roi`) and the `truth_mupair_pass` branch (`pt>=3800`), neither of
which filters. Raising the analysis cut to 4.5 GeV is fully contained in the existing skims.

**The gap-window change is binning-neutral and is a ONE-LINE edit.** Only the *forward*
window edge is slaved to a binning edge (the `SetIOPaths` startup throw); −1.30/−1.25 both lie
strictly inside coarse q·η bin (−1.5,−1.0). All 12 live sites read the vector. The pT change
strictly contains this change's blast radius, so both go in one pass.

**★ Hazard 1 — the stale-file guards would NOT have caught the muon-pT move.**
`PairRecoEffEvaluator`, `DrCorrectionCrossxEvaluator`, `PairTrigEffEvaluator` and
`plot_pp24_fullsim_pair_reco_eff` all `CheckCanonicalBinning` against pair-pT / pair-η / ΔR
AXIS EDGES. A muon-pT threshold change moves none of them, so every stale ε_reco / ε_ΔR /
ε_pair file would have loaded clean and been applied to a 4.5 GeV sample. It is only the
**pair-pT 8→9** half of this change that trips all four, loudly, at `Load()` time (before the
event loop, so genuinely loud). The provenance `TNamed` that `build_pp24_fullsim_pair_reco_eff.C`
writes names the gap cut and the axes but **not** the muon-pT threshold, and
`PairRecoEffEvaluator` never reads it.

**★ Hazard 2 — silently-stale artefacts (nothing checks them at all).** The data mu4 turn-on
TF1 files (`.../trg_effcy_pT_fitting_to_{erf,fermi}_plus_log/single_mu_effcy_pT_fit[_medium_wp].root`)
are loaded BY NAME and a stale file reweights every pair with no complaint; likewise the T&P
graph inputs, the MC singles map, the Step-3/4 plateau file (its guard is a *value* check,
|plateau−1|≤0.1, not a provenance check), and every NTuple tree — **no output filename anywhere
carries a muon-pT token**, so 4.0 and 4.5 outputs overwrite each other. ⇒ back up before
overwriting, and validate stages by freshness + a required histogram, never by file-exists.

**★ Hazard 3 — a Pb+Pb crossx output is ALREADY destroyed.**
`pbpb_2024/histograms_real_pairs_pbpb_2024_single_mu4_no_trg_plots_nominal.root` is **851 bytes**
(2026-09-06) — the swallowed-RDF-exception failure mode, from the `pp_trig_eff_highpt_jump.md`
throw. 2023 (3.6 MB) and 2025 (4.3 MB) are last-good 2026-07-08. All three are refilled by this
task, but they must be backed up first.

**★ Hazard 4 — a second, independent Pb+Pb blocker, newly measured.** All three Pb+Pb turn-on
fit files contain 240 TF1s with complete coverage, but their top bin is named `_2_00_TO_2_30`
while the reader now builds `_2_00_TO_2_20` (the coarse q·η edge moved 2026-09-07; pp24 was
refit, Pb+Pb was not). Every Pb+Pb muon with q·η ∈ [2.0,2.2) therefore throws — deterministic
and independent of the pre-existing sentinel bug. The Pb+Pb T&P refit in this task cures it.

**POWHEG fullsim re-enters the blast radius.** `pp24_all_vertex_pairs.md` carved it out because
its only live node is truth-level and an impact-parameter/reco change cannot reach it. That
argument does **not** survive the pT change, which IS a truth-pT cut — so
`run_powheg_fullsim_wtruth_{bb,cc}.sub` (11+11 jobs) must run for the MC-vs-data comparison curve.

**Interrupted predecessor.** The pp24 fullsim mc-trig pair file written at 15:48 by the
abandoned `pp24_all_vertex_pairs.md` step 7a reports "probably not closed / recovered 6 keys",
i.e. a truncated file from a killed local job. It is regenerated here regardless.

**Retyped binning copies that will NOT follow the `ParamsSet` edit** (the exact Binnings failure
mode, live today): `SingleBAnalysis/SingleBAnalysisBase.cxx` and
`plotting_codes/single_b_analysis/plot_dr_vs_pair_pt_diagnostic.cxx` each re-derive the fine
pair-pT axis as `(15, 8.0, 120.0)`; `RDFBasedHistFillingPbPb.cxx` retypes `44, -2.4, 2.4`. All
three are fixed in this change.

**Additional pT axes not in the original Done list, folded in under D2/D3** (same rule — the low
edge follows the cut that defines it, bin counts unchanged): `pT_bins_8` 20 log 4→8 ⇒ 4.5→8
(this is the low half of the trigger-efficiency `pt2nd` axis — leaving it at 4 would have put
the graph point at a bin CENTRE of 4.517 while the survivors' mean is ≈4.55, an x-shift exactly
where the turn-on is steepest); `pT_bins_40` 18 log 4→40 ⇒ 4.5→40; `pT_bins_80` 12 log 8→80 ⇒
9→80. `pT_bins_60` (8→60) is the single-muon HIGH half and 8 is not a cut there — unchanged.

### R3 — ★ The turn-on degeneracy risk that changed the plan (measured, 2026-09-08)

Measured READ-ONLY from the 380 stored nominal turn-on fits (`_sub_edit_turnon_fit_range.md`),
before anything was rerun. This is why D8/D10/D11 exist.

1. **~40 % of the 4→6 GeV mu4 rise lies inside the 0.5 GeV that a 4.5 probe cut would remove**
   (pp 0.400; PbPb 23/24/25 0.428/0.399/0.421; per-cell range 0.13–0.99). The removed segment is
   worth 15–19 % of plateau height on average, up to 53–73 %.
2. **In 63–80 % of q·η cells the fitted turn-on MIDPOINT already sits below 4.5 GeV** (pp 16/20),
   and in 31–41 % it is more than one full width below. Fitting above 4.5 would therefore
   determine the turn-on location by extrapolation, with (mean, sigma) degenerate along a flat
   interior valley — limits ([0,10] × [0.01,50]) are vastly wider than the fitted values
   (2.76–5.10, 0.60–1.55), so there is no wall anywhere near it.
3. **The existing warning provably cannot catch this.** Across all 380 nominal fits the
   at-a-fit-limit `*` flag fires only on `corrCoef` at its 0 boundary and fires on a SHAPE
   parameter **zero times, in any file**. χ²/ndf is blind by construction: a degenerate fit
   reproduces the retained points exactly. The MC twin's `status/GetN/GetNDF` test is equally
   blind.
4. Note the axis was **rebinned, not truncated** — only 3 of 40 points are lost — so the damage
   would have been loss of *leverage on the turn-on location*, not of statistics.

**Resolution (user):** data probes stay at 4.0 (D8), MC drops the population entirely (D9), the
measurement axis forces 4.5 to be an edge so MC has no partial bin (D10), and the two fit ranges
diverge on purpose (D11). A diagnostic remains worth printing after the refit even so:
**σ_ε(4.5)/ε(4.5)** from the fit covariance (`"QR"` → `"QRS"`, then `GetConfidenceIntervals` at
`pT_min`), flag > 0.10, stop > 0.20; companions `|ρ(mean,sigma)| > 0.95` and the census
`d = (pT_min − loc)/width > 1`.

### R4 — ★ The textual-vs-inherited rule, and which on-disk MC artefacts are invalid right now

Established 2026-09-08/09 with the `polyn fit restriction` and `additional fullsim statistics`
sessions, after this doc's first characterisation ("the whole 17:11→18:52 chain is a chimera")
was shown to be too broad. **The correct boundary is mechanical, not per-pipeline:**

> **A fill that re-applies the cut TEXTUALLY is immune to a stale tree; a fill that INHERITS the
> cut from a persisted flag is not.**

Because the NTuple selection is a per-pair cut, a 4.0 GeV tree is a strict SUPERSET of a 4.5 GeV
one: `{quality, pt>4.0} ∩ {pt>4.5} == {quality, pt>4.5}`. So any selection that re-states
`pt > 4.5` reproduces the 4.5 population exactly, even from an old tree. A selection that only
says `pair_pass_tight` gets whatever threshold was compiled into that branch.

**Measured partition** (grep for an explicit `pt > 4.5` / `truth_pt > 4.5` term, following the
shared header):

| Consumer | Explicit muon-pT term? | Verdict |
|---|---|---|
| `FillMCTrigEffHists.cxx` (Steps 1–4) | yes, 12 sites, on every leg | **IMMUNE** |
| `FillMCTrigEffClosure.cxx` | yes, via `MCTrigEffPairSel::Step3PairSelection()` / `SingleBSignalCutsReco()` | **IMMUNE** |
| `FillMCTrigEffPairEff.cxx` | yes, same shared header | **IMMUNE** |
| `RDFBasedHistFillingPythiaFullsim.cxx` | **none anywhere in the file** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPowhegFullsim.cxx` | **none** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPythiaFullsimOverlay.cxx` | **none** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPythiaTruth.cxx` (truth signal acceptance + template `kin_cuts`) | **none, reco OR truth** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPowhegTruth.cxx` (signal acceptance) | **none** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPP.cxx` (T&P probe leg, generic family) | **none** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPbPb.cxx` (T&P probe leg, generic family) | **none** | **STALE until the NTuple rerun** |
| `RDFBasedHistFillingPowhegFullsimSingleMuon.cxx` | **none** | **STALE until the NTuple rerun** |

**CORRECTED 2026-09-09 (iteration-3 review):** this table first listed only the three fullsim
classes. It is **eight**. The truth signal-acceptance and template `kin_cuts` outputs carry no
muon-pT term of any kind, reco or truth, so they are exactly as stale as the fullsim ones — and
the earlier wording ("the artefacts that are genuinely invalid on disk right now are the ones the
three fullsim classes produce") under-reported the blast radius.
Also tightened: the closure and pair-eff rows are immune **entirely** through
`Step3PairSelection()`. `SingleBSignalCutsReco()` carries NO muon-pT term at all, so naming it as
a second source of immunity would invite a future filter reorder that silently removes it.

⇒ The artefacts that are genuinely invalid on disk right now are the ones the three fullsim
classes produce: **ε_reco (`pair_reco_eff_pp24_full.root`, written 18:52 today), the detector
response, the template-fit MC, and the overlay per-centrality ε_reco** — not the trigger-efficiency
chain, which this doc initially and wrongly included.

**★ The sharpest form of it** (found by `polyn fit restriction`): in
`RDFBasedHistFillingPythiaFullsim.cxx` the `pass_signal_truth` / `pass_signal_reco` columns are
RDF-**Defined**, built from the LIVE `ParamsSet::SignalPairPtCutExpr` and the live gap expression,
while the muon-pT term sits in the persisted `pair_pass_*`. So a fill from an old tree writes
**one selection string carrying two vintages**: new 9 GeV pair-pT and new (−1.25,−1.05) gap, old
4.0 GeV muon pT. The `_pass_{medium,tight}_sig_*` variants are therefore MORE internally mixed
than the plain `_weighted` ones, not less.

**Why no guard caught it, and the standing recommendation.** Every stale-file guard in the chain
(`PairRecoEffEvaluator`, both `DrCorrection*Evaluator`s, `PairTrigEffEvaluator`,
`PairEtaPanelBins`) compares AXIS EDGES. A muon-pT threshold change moves no axis, so none of them
can fire on it — R1 Hazard 1 predicted exactly this before it was observed. The durable fix is to
stamp the muon-pT threshold into artefact provenance and have those guards compare it; that was
put to the user on 2026-09-08 and NOT adopted in this round, so the exposure remains for the next
threshold change. All of these artefacts are regenerated by this task's rerun, so nothing here
requires a code change — only that nothing consumes them in the meantime.

### R2 — Step 3 DONE: `ParamsSet.h` (the keystone)

Gap window `{-1.30,-1.05}` → `{-1.25,-1.05}`; `pTbins` low edge 4 → 4.5; `pT_bins_40`,
`pT_bins_8` → 4.5; `pT_bins_80` → 9; `pT_bins_120` → **16 log 9→120** (alternative);
`pT_bins_150` → **16 log 9→150** (new default); `pair_pt_coarse_bins` and its 4-bin variant →
9→150; `single_mu_pt_coarse_bins` → {4.5, 8, 14, 25, 100}; `pt_titles[0]` relabelled.
Comment blocks updated for truthfulness, including marking the measured gap-cut cost and ε_acc
tables as STALE IN TWO WAYS (superseded window AND measured at pT > 4).

## Remaining Work

Implementation Plan steps 9–15: finish compiling, local smoke tests,
`/review-analysis-code`, then the Condor reruns and the whole downstream chain.

Carried forward, NOT part of this task:
- `plot_reco_distr_singleb_vs_op_pp24.C` is left on the retired `q·η < 2.2` **and** on a
  `dr > 0.05` cut that was removed from the analysis on 2026-06-22. Doubly stale, unmaintained,
  loose top-level macro; recorded rather than fixed.
- The legacy pre-RDF `SingleBAnalysis/SingleBAnalysisBase.cxx` still retypes the fine pair-pT
  axis as `(15, 8.0, 120.0)` and its own `signal_cuts`. Not in the active chain (see
  `signal_selection_change_impact.md` §2 "Legacy / retired"); left, and listed here so it is not
  revived unnoticed.
- Dead NTuple copies (`PythiaNTupleFirstPass*`, `original_no_template_class/*`,
  `maybe_old_unsure_pbpb/*`) still carry the 4.0 GeV threshold — reviving any of them revives it.

## Latest Stage

**CODE EDITS COMPLETE; COMPILING. Nothing has been rerun and no Condor job submitted.**

Done so far (steps 1–8):
- `ParamsSet.h` keystone (R2) + a NEW `signal_pair_pt_min` / `SignalPairPtCutExpr` single source
  of truth replacing the pair-pT threshold that had been RETYPED at ~10 sites, and the
  two-segment `pT_bins_8` of D10.
- NTuple processing: all 7 live pT thresholds → 4.5, plus the D8 mode-dependent form (the
  trigger-efficiency mode keeps 4.0). Dead copies deliberately left at 4.0 and inventoried.
- Pb+Pb MIGRATED (D1): fiducial + pair-level gap cuts adopted in the signal region and both
  template blocks; the retyped `44, -2.4, 2.4` pair-η axis replaced by
  `ParamsSet::N_PAIR_ETA_CROSSX_BINS` at all 10 sites.
- Truth/fullsim/overlay analogs migrated off the retired `q·η < 2.2` (D7); the acceptance
  cutflow now READS its cut list from `ParamsSet` instead of retyping it.
- Trigger-efficiency selection: thresholds moved and the one genuine gap the audit found (the pp
  same-vertex pair requirement, missing from Steps 2/3/4 + closure + pair-eff) closed via
  `pair_pass_*`, behind a throwing pre-loop guard. `SingleBSignalCutsReco` now reads the pair-pT
  threshold from `ParamsSet` — it had been a silent mirror (its guard compares only the `minv`
  half).
- Fit ranges split per D11; the MC fitter's "mirrors the data fitter EXACTLY" header amended.
- Temporary `_diag_mupt45` / `_logbins_from9` diagnostic code DELETED (3 files + 384 lines of
  plotter methods + the PP.cxx block), zero residual references.

**Verified, not assumed** — a standalone ROOT test of the new axes prints:
gap windows `(-1.25,-1.05) (-0.10,0.06) (2.20,2.40)`; `signal_pair_pt_min = 9`;
`pTbins[0] = single_mu_pt_coarse_bins[0] = 4.5`; `pT_bins_8` 20 bins 4.0→8.0 **with 4.5 an edge
at index 2**; `pT_bins_150` 16 bins 9→150 and `pair_pt_coarse_bins` 8 bins 9→150 with
**coarse[k] == fine[2k] for all k** (new coarse edges 9, 12.79, 18.18, 25.85, 36.74, 52.23,
74.24, 105.53, 150); `_pt_120` 16 bins 9→120. ALL CHECKS PASSED.

**D5 IMPLEMENTATION INVERTED (2026-09-08), on the strength of the coverage-gap report.**
The first approach kept every histogram name and just re-pointed the plotters at the `_pt_150`
family. The gap report showed that family is INCOMPLETE: no `_pt_150` twin exists for the minv
and ΔR 2D histograms (pp and Pb+Pb, plus their `_counts`), for `_no_trig_corr`, for the four
correction-stage suffixes, for the same-sign crossx, for the 3D minv — **or for the two R_AA
global 3Ds**, which would have left R_AA on the 9→120 alternative while the cross-section sat on
9→150. Two coexisting pair-pT binnings is precisely the failure the Binnings rule was written
after. Booking ~15 new twins would also have cost memory in every job.
**Inverted instead:** the UNSUFFIXED family — which is already the complete nominal set — is
rebooked on `pT_bins_150`, and the partial `_pt_150` family is rebooked on `pT_bins_120` with its
name token renamed to `_pt_120`. The coverage gap disappears because the default now uses the
family that already covers everything, and no new histogram is created. Partial coverage is
acceptable for an opt-in alternative view. The 1D `dsigma` family and the template-fit
`pair_pt_log_150` histograms were ALREADY on `pT_bins_150`, so they are correct untouched and
their names stay truthful.

**Three bugs caught in my own edits during this step, by verifying rather than assuming:**
(i) a scripted replacement left the Pb+Pb **R_AA global 3D** with its bin COUNT from
`pT_bins_150` and its EDGES from `pT_bins_120` — a silently malformed axis on the R_AA input;
(ii) the Pb+Pb per-centrality nominal block (the family BOTH the TAA-weighted and the counts
plots are drawn from) was left on the 120 axis; (iii) 10 further retyped `44`-bin pair-η axes in
the `44, eta_edges.data()` spelling, which the earlier `44, -2.4, 2.4` regex could not see.
All fixed; a sweep now confirms every `npt`/`ptbins` pair agrees on its vector and no bare `44`
remains in Pb+Pb code. Earlier in the same step a brace-matcher that ignored string literals cut
three method signatures and left their bodies — the file was restored from git and the deletion
redone with a literal-aware matcher.

Compile status: `RDFBasedHistFillingPP`, `RDFBasedHistFillingPbPb`, the three NTuple analysis
class headers (so the D8 mode-dependent cut compiles), `PythiaTruth`, `PowhegTruth`,
`PythiaFullsim`, `PowhegFullsim`, `PythiaFullsimOverlay`, `FillMCTrigEffHists` and
`FillMCTrigEffClosure` all compile with **no errors**; the remaining three classes and a PP/PbPb
recompile after the inversion are in flight.
