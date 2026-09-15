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
  3. **Binnings** moved (user decisions D2/D3/D10): `ParamsSet::pTbins`,
     `single_mu_pt_coarse_bins` and `pT_bins_40` low edge 4 → 4.5; `pT_bins_8` becomes TWO
     SEGMENTS so 4.5 is a forced interior edge; `pair_pt_coarse_bins` (+ 4-bin variant),
     `pT_bins_80`, `pT_bins_120` and `pT_bins_150` low edge 8 → 9, with `pT_bins_150` = 16
     log bins 9→150 as the DEFAULT fine crossx axis and `pT_bins_120` = 16 log bins 9→120
     as the opt-in `_pt_120` view.
     ⚠ **SUPERSEDED, do not act on the original wording of this item:** it said the turn-on
     fit range moves to `[4.5,60]` in BOTH fitters. **D11 reversed that.** The DATA fitter
     `SingleMuEffcyPtTurnOnFitter` deliberately STAYS at `[4,60]` with the pivot at 4.0,
     because the trigger-efficiency NTuple mode keeps 4.0 GeV probes (D8); only the MC twin
     `FitMCSinglesEffcy` uses `[4.5,60]`. Do NOT "fix" the data fitter to 4.5 — R3 measures
     what that would cost.
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
- **D12 (2026-09-09, user).** A PRE-EXISTING normalization hole, surfaced by the
  `additional fullsim statistics` session and confirmed independently: `PythiaAlgCoreT.c`
  built the weight as `ami_w * nom_ratio / N_beam` while the event loop runs over
  `N_proc = min(N_beam, nevents_max)`, so a run truncated by `nevents_max` wrote an absolute
  cross-section low by exactly `N_proc/N_beam`. Undetectable: the weight, the NTUP chain
  entry count and AMI `totalEvents` all return `N_beam` and AGREE with each other. Reachable:
  the pipeline smoke test passes `NEVENTS_MAX` through the same run script while
  `extra_output_suffix` stays `_full`, and `--dry-run` does not change the output path, so a
  2000-event file can land on the nominal filename and pass Stage 4's exists+non-empty check.
  Dates to `6f78163` (2026-04-16); NOT introduced by this task.
  **Fixed on both counts.** (a) `meta_tree_out` now records `N_proc` and `N_beam` per
  (kn, beam) and is Filled on the FULLSIM path — it used to be Filled only under
  `getIsPrivate()`, so every fullsim file carried an empty meta tree; a truncated chain also
  prints a loud warning. (b) `N_proc` is hoisted ABOVE the weight and used as the denominator
  in the fullsim AND the non-private-truth branch, and `p->crossx` multiplies back the SAME
  count (`* N_to_process`, was `* N_beam`) — fixing only the weight would have reintroduced
  the factor in `crossx`. **A no-op for every nominal run** (`N_proc == N_beam` when
  `nevents_max` is unset), so nothing this rerun produces can move because of it.
  NOT adopted: giving `--dry-run` its own output suffix (offered, user chose the narrower fix).
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

- 2026-09-09 Steps 9–15 (partial). Everything COMPILES (21 targets: 3 NTuple analysis headers,
  11 RDF/trig-eff classes ACLiC, 7 interpreted macros). `/review-analysis-code` ran THREE
  iterations — see R5. D12 implemented (both halves). **COMMITTED in 7 logical commits:**
  `06f9bfb` ParamsSet/binnings · `4506fd9` NTuple pT + D8 + D12 · `0054c90` trig-eff 4.5 +
  same-vertex · `139262d` fitters/D11 · `e5e5312` Pb+Pb migration + pair pT 9 + default axis ·
  `01031f0` plots/pipelines (R_AA top bin, Pb+Pb serialization, stale-macro gate) ·
  `c9ab2df` docs. Working tree clean apart from `.claude/logs/tracking.jsonl`.
  **NOTHING HAS BEEN RERUN. No Condor job has been submitted. No output file overwritten.**

- 2026-09-09 22:30-23:00 **REVIEW ITERATION 4 — the decision to run it was correct.** Two independent
  read-only reviewers on the previously unverified `01031f0` + `c9ab2df` batch, split by scope
  (A = `pipelines/` + run/submit scripts; B = `plotting_codes/` + `RAA_plotting.cxx` + docs).
  **Both returned FAIL.** Combined: **6 CRITICAL, 31 WARNING, 6 INFO** — see R7. Three of the
  CRITICALs would have silently corrupted or wasted this rerun, and one is a plain build break in a
  file that had been reported as compiling. Merged findings in R7; the Phase-1-reaching half is
  fixed and committed (`abb9103`, `227041f`) BEFORE any job was submitted.

- 2026-09-09 22:54 **PHASE 1 LAUNCHED.** All RDF classes were pre-compiled serially first
  (13 targets; `RDFBasedHistFillingPowhegFullsimSingleMuon` was MISSING a `.so` and was built now),
  so every `.so` is newer than every header and the parallel pipelines only LOAD — nothing compiles
  concurrently. That matters because `RDFBasedHistFillingPP.cxx` and `...PbPb.cxx` both include the
  unguarded `RDFBasedHistFillingData.cxx`.
  Five drivers, on DISJOINT output trees, each harness-owned so its exit notifies:
  | Job | Driver | Covers | Writes to |
  |---|---|---|---|
  | A | `ENABLE_MC_TRIG_EFF=1 pipeline_pythia_fullsim_pp.sh full` | 1g (local, 6-8 h) | `pythia_fullsim_full_sample/` |
  | B | `pipeline_pp_trig_eff.sh` | 1b + 2a(pp, Tight) | `pp_2024/` |
  | C | `SKIP_EVSEL=1 run_pbpb_all.sh` | 1c + 2a(Pb+Pb, Tight+Medium) + 3b | `pbpb_2023/24/25/` |
  | D | `pipeline_pythia_truth.sh` x3 modes | 1d | pythia truth dirs |
  | F | `pipeline_pythia_fullsim_overlay.sh hijing` | 1f | overlay dir |
  Plus Phase 1e submitted directly (no pipeline exists for it): clusters **2420** `powheg_truth_bb`,
  **2421** `powheg_truth_cc`, **2422** `powheg_fullsim_wtruth_bb`, **2423** `powheg_fullsim_wtruth_cc`.
  Job A Stage 0 preflight PASSED: all 6 pT-hat slices and exactly 6 AMI files present, so the
  `expected_ami_dsids` guard is armed on the FULL `_pdf` production.
  Logs: `.claude/logs/pipeline-runs/mupt45_rerun_20260909_225405/`. A persistent monitor tails all
  of them for stage transitions AND failure signatures.

- 2026-09-10 **PHASE 1 COMPLETE; PHASE 2 COMPLETE except POWHEG; PHASE 3 IN FLIGHT.** Done since the
  launch entry above:
  * **Phase 1 finished** — all 97 Condor jobs plus the ~7 h local fullsim pass. Jobs A/B/D/F exited
    0; C's trig-eff half exited 0 and its crossx half failed (R8) and was fixed and rerun.
  * **★ Independent verification that the whole point of the rerun took effect** (the pipelines'
    own checks were inert for all of Phase 1 — see R9 — so this was done from outside):
    every trig-eff-mode tree sits at EXACTLY 4.0 GeV on purpose and every Pb+Pb nominal tree at
    EXACTLY 4.5, with zero violations either way. **D8's mode-dependence works end to end.**
    | tree | OS pairs | min muon pT | below threshold |
    |---|---|---|---|
    | pp24 trig-eff | 779 603 | 4.0000 | 0 below 4.0 |
    | pbpb 23 / 24 / 25 trig-eff | 799 732 / 565 081 / 1 602 734 | 4.0000 | 0 below 4.0 |
    | pbpb 23 / 24 / 25 nominal | 434 966 / 306 628 / 869 072 | 4.5000 | 0 below 4.5 |
  * **Phase 2a DONE** — data mu4 turn-on refits, pp AND Pb+Pb, **Tight AND Medium**, all fresh.
    The Pb+Pb refit cured R1 Hazard 4: the fits now key their top q*eta bin `_2_00_TO_2_20`, the way
    the reader builds it. `pbpb_2024` and `pbpb_2025` Medium fits were CREATED (neither existed).
  * **Phase 2b DONE for pp** — job A's Stage 10 ran the MC trig-eff chain end to end
    (FillMCTrigEffHists step1 → FitMCSinglesEffcy → step3 → plot_mc_trig_eff) on the fresh sample.
  * **Phase 2e DONE** — `pair_reco_eff_pp24_full.root` rebuilt (14:10, was 2026-09-08 18:52).
    Sanity-checked: every efficiency in [0,1], 143 filled 3D cells, no unphysical value.
  * **Phase 3b DONE** — Pb+Pb crossx completed for 2023/2024/2025 after the R8 fix (3.2 / 2.8 /
    3.9 MB). **`pbpb_2024` is a real file again for the first time since 2026-09-06**, when it
    became an 851-byte corpse. The Stage-7 sanity plot drew the **Pb+Pb panels** successfully,
    which independently confirms the D1 48-bin pair-eta migration is live on the refilled
    histograms — the `PairEtaPanels::Bins` guard passes where it used to throw.
  * **Phase 3a IN FLIGHT** — `pipeline_pp_crossx.sh` (cluster 2429): pp24 nominal NTuple → hadd →
    crossx → plots, with `INCLUDE_PBPB_SANITY=true` now that Pb+Pb is refilled (the flag flip the
    A7 fix was built for). Run with `SKIP_MC_DATA_COMPR=1` because POWHEG was not ready — see R10.
  * **Peers pinged** (R6 discharged): `polyn fit restriction` and `additional fullsim statistics`
    both told the fresh trees and Step-3 histograms exist. The latter replied confirming and
    FIXING both defects reported to it, and flagged one stale retyped edge in
    `Utilities/PairTrigEffEvaluator.h` (mine), now de-numeralised.

- 2026-09-15 **Old-binning plot cleanup (user request).** Plot families whose filenames carry the
  pair-pT edges are NOT overwritten by the rerun (new edges = new names), so old 8-GeV-axis PNGs
  sat next to the new 9-GeV ones. Every plot directory under `dimuon_data/plots`, the fullsim /
  overlay / truth / POWHEG sample `plots/` and `mc_trig_eff_fit_plots/` was scanned for files
  OLDER than the Phase-1 launch (2026-09-09 22:54) coexisting with NEWER ones. **752 files moved
  (not copied) into a per-directory `old_4GeV_cut/` subdirectory** (102 subdirs), list at
  `dimuon_data/plots/old_4GeV_cut_moved_files_20260915.txt`:
  * 732 old-edge pair-pT PNGs (`pairpt_8.0_11.5 … 104.0_150.0` / `72.1_150.0`) in
    `pp_trigger_efficiency/mc_based{,_medium}/step3_dr_correction/dr_*` (48) and
    `pp_trigger_efficiency/mc_based/step3_dr_fit/{5 variants}/{expo,interp,polyu_fixedRp}/sign_*`
    (684). Verified: zero moved files carry a new edge, zero are newer than the launch.
  * 20 other old pT-binned artefacts in the same mixed dirs: `step1_eff_2d_pt_vs_q_eta_charge_*`
    (4, optional 2D macro not rerun), `dr_correction_data_and_fits_tight_os.root` (2, 09-03,
    `combine_dr_correction_fits` not rerun; no reader in code), the Aug `_corrected*` MC turn-on
    fit PNGs in `pythia_fullsim_full_sample/mc_trig_eff_fit_plots` (8), and the June `mu4noL1`
    single-muon plots in `pp_trigger_efficiency/mu4/no_corr/pp24/single_muon_effcy` (6).
  Left in place: non-pT files (`README.md`, `single_b_analysis/nevt_*`, 2025 `pythia_private_sample`
  files). No IntNote manifest entry or live code path references a moved file.
  **NOT touched — old-ONLY directories on the 8-GeV axis with no new sibling** (stale, not
  confusing in the same way; user to decide): `pp_trigger_efficiency/mc_based_medium/step3_dr_fit/*`
  and `mc_based{,_medium}/step4_dr_fit/plateau_corrected` (Medium ΔR fits / Step-4 fits not refit),
  ALL of `pbpb_trigger_efficiency/mc_based*` and `r17663_no_overlay_trigger_efficiency/*` (the
  Pb+Pb / no-overlay MC trig-eff chain was not rerun at all), `*_fine_q_eta_bins_w_gap/mc_based`,
  the `mc_based_pt4bin*` variants, the already-named `mc_based.bak_testsample_*` snapshots, and the
  `_pt_150` directories of Open item 2.

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

### R5 — Review outcome: three iterations, and what the pattern says

`/review-analysis-code`, log
`.claude/logs/review-analysis-code-20260909-000600-mupt45-gap125-pairpt9-adoption.md`.

| Iter | CRITICAL | WARNING | Character |
|---|---|---|---|
| 1 | 1 | 9 | R_AA silently dropping its top pair-pT bin, + stale text |
| 2 | 0 | 8 | Stale text ONLY — no code or physics defect |
| 3 | 1 (PRE-EXISTING) | 20 | First sweep of `plotting_codes/` + `pipelines/` |

**Iteration 1's CRITICAL:** `RAA_plotting.cxx` grouped the fine pair-pT axis with a hardcoded
15-bin map while the axis moved to 16, so bin 16 (125.8–150 GeV) was dropped from R_AA mode 3
and all three labels were wrong. Its two existing guards compare pp against Pb+Pb and both axes
moved together, so neither could fire. Fixed at root: the grouping is DERIVED from `ParamsSet`,
split over the COARSE cells and expanded 2:1 so every group boundary is a `pair_pt_coarse_bins`
edge — each R_AA group is a whole number of the cells the corrections were measured in
(3+3+2 cells → 9–25.8 / 25.8–74.2 / 74.2–150 GeV) — plus a nesting assert and a coverage assert
that fires if the grouping does not tile the histogram read from disk.

**Iteration 3's CRITICAL is D12**, pre-existing and now fixed.

**The dominant defect class, by a wide margin, is a stale CLAIM rather than wrong code** — a
comment, a doc, a drawn label, or a string emitted at runtime. Two would have reached output: the
R_AA figure legend published `p_T^pair > 8 GeV` on a final-results plot, and the ΔR evaluators
printed "pair pT < 8 GeV" into every job log. Everything of that kind is now COMPOSED from
`ParamsSet` rather than retyped.

**★ The lesson that matters for a future round: my own amendments were not reliably complete.**
Iteration 3 found FOUR iteration-2 fixes half-applied — most tellingly, a printed string left
untouched while the comment beside it was corrected. Amendments must be re-verified, not assumed;
a 4th iteration has NOT been run and the last batch of fixes is unverified by a reviewer.

Two operational defects found in iteration 3 that would have broken THIS rerun:
`run_pbpb_all.sh` launched trig-eff and crossx concurrently although crossx consumes trig-eff's
turn-on fits (now serialized), and the Pb+Pb combiner silently dropped years missing a histogram
while still labelling the canvas "2023, 2024, 2025 combined" (now probes every year and reports
contributors).

### R7 — Review iteration 4: what the "unverified batch" actually contained

Log `.claude/logs/review-analysis-code-20260909-223311-mupt45-iter4-amendment-verify.md`.
Two reviewers, disjoint scope, both **FAIL**. R5 predicted the amendments were unreliable; they
were. **The single most important number here: a file that this doc recorded as one of "21 compile
targets clean" does not compile.** So "it compiled" was itself an unverified claim.

**Phase-1-reaching — FIXED before any job was submitted (`abb9103`, `227041f`):**

| # | Defect | Why it mattered |
|---|---|---|
| A1 | `pipeline_pythia_fullsim_pp.sh` validated Stages 1-3/5 by existence + entry count, never freshness | Stages 1-3 `cp -a` the OLD file and leave it at the nominal path; a throw inside the event loop is caught by TRint and the job exits 0 ⇒ the whole 6-8 h chain measures ε_reco and the detector response on the pre-change 4.0 GeV population and applies it to new-selection data. Now: `mktemp` stamps + `require_fresh` on every product, a real tree check on the single-muon file, and a FILLED-histogram probe replacing the key count (a key count passes on a partially flushed file) |
| A3 | `run_pbpb_all.sh` `YEARS=(${YEARS:-...})` destroyed the scalar the children read | Arrays are not exported ⇒ every sub-pipeline saw `YEARS` unset and ran all three years. `YEARS=24 SKIP_CONDOR=1 ./run_pbpb_all.sh` — the natural re-entry after one year fails — would resubmit ~24 jobs and overwrite good output. Verified empirically before AND after the fix |
| A9 | Nothing in the Pb+Pb path regenerated the **Medium-WP** turn-on fits | `pipeline_pbpb_trig_eff.sh` leaves `isTight` true, so Medium would have been the ONLY correction left on the superseded gap window and the old `_2_00_TO_2_30` q*eta key. A Medium stage now runs between trig-eff and crossx |
| A12 | The mc_trig products were not in the backup loop | Its comment promised to back up "anything we are about to overwrite" and did not cover `ENABLE_MC_TRIG_EFF=1` — the configuration this rerun uses |
| A11 | `pipeline_powheg_fullsim_single_muon.sh` validated an explicit part list then hadded a **glob** | The glob also matches this repo's own `.bak_<timestamp>.root` names; a double-counted batch in a reco-eff denominator is invisible downstream. Now merges the arrays it validated |
| A5 | `run_pbpb_all.sh` serialization diagnostics were dead code under `set -e` | Safety property was fine (crossx cannot start early); the operator just got a generic ERR line instead of "ABORTING before crossx" |

**★ A fix of mine that was itself wrong, caught by testing rather than assuming** — the R5 lesson,
live: the reviewer's suggested `if ! cmd; then RC=$?; fi` (and my first application of it) captures
the status of the `!` NEGATION, which is **always 0**, so a child exiting 7 reports `rc=0`. The
working form is `cmd || RC=$?` with `RC` preset. Also caught the same way: I guessed the mc_trig
single-muon suffix as `_single_muon_mc_trig` when the real composition is `_mc_trig_single_muon`
(trig_suffix precedes extra_output_suffix), so the freshness gate I had just added would have
aborted Stage 4 of a *correct* run. Both were found by exercising the change against the real
shell / the real filenames.

**NOT acted on, deliberately:** reviewer A re-raised (as CRITICAL) that `--dry-run` writes to the
NOMINAL output filenames. D12 records that giving `--dry-run` its own suffix was offered to the
user, who chose the narrower fix. Left as the user decided; this rerun does not use `--dry-run`.

**Outstanding — NOT Phase-1-reaching, being fixed while Condor runs.** Reviewer A: #4
`run_mc_trigeff_round7.sh` accepts a previous production's file as proof a stage ran and (no
`set -e`) ignores ROOT's exit code entirely (Phase 2b); #6 pp crossx required-histogram probe's
failure message unreachable under `set -e`; #7 the pp pipeline still skips the Pb+Pb sanity panels
citing a 44-bin axis that no longer exists; #8 Pb+Pb crossx Stage 5 lacks the stamp + probe its pp
twin has, on the exact stage that produced the 851-byte corpse; #10 `run_all_crossx.sh` fills pp24
via `test_crossx_pp24.sh`, which omits the `_wgapcut` family the pipeline version writes, and gates
four RDF stages on PNG existence alone. Reviewer B: **B1 build break** (`const int nb = 6;`
swallowed into a `//` comment in `plot_single_muon_reco_effcy_r17618_vs_r17662.cxx:45`); **B2**
Pb+Pb combined crossx records `combined_years_` and never reads it, so the drawn label still claims
all three years while the code's own warning says it must not; **B3** `signal_selection_change_impact.md`
asserts three lines above its own correction that Pb+Pb adopted neither gap cut; **B5** the new
mode-3 coverage assert checks only the bin COUNT and `pT_bins_120` has the same count (16), so the
one realistic wrong axis passes; plus ~20 stale comments/headers/labels and the doc items.

**Two findings that belong to a PEER session, not this one** (`additional fullsim statistics` owns
`mc_pthat_slice_mass_statistics.cxx` and `pthat_slice_mass_stats_sample_request.md`) — flagged to
the user, NOT edited here:
- **B23 is physics-results-bending.** `drop_sigma`'s "above" region is anchored on the data's own
  peak bin while the header and the emitted CSV both say the two regions were "fixed before looking",
  so the region systematically excludes its own largest bin and `nbins_above` varies 6/7/9/15/16
  across rows (different mass ranges compared as if alike). Reviewer measured the effect on the
  shipped histograms: pTH125_300 OS top-3 **+4.5σ → +3.9σ**, pTH70_125 SS top-3 **−4.7σ → −5.6σ**.
  That moves the peer doc's stated decisive number.
- **B22**: that macro ships an `N_beam`-vs-`N_proc` "residual risk" narrative, in the code and in all
  9 emitted CSV headers, that D12 already fixed upstream — harmless today (`N_proc == N_beam` for a
  nominal run) but it would now misdiagnose a truncated file.

### R8 — ★ The Pb+Pb crossx blocker: a SECOND latent throw, exposed by curing the first

**Phase 1 ran; Pb+Pb crossx was the one failure.** Everything else completed: job A (Pythia
fullsim pp24 FULL, ~7 h, Stages 1-10 incl. MC trig-eff), B (pp trig-eff + Tight turn-on refit),
D (Pythia truth x3), F (HIJING overlay), and the Pb+Pb trig-eff + Tight + **Medium** turn-on refits
for all three years.

**What happened.** `pipeline_pbpb_crossx.sh` Stage 5, year 2023, threw inside the RDF event loop:

> `EvaluateSingleMuonEffcyPtFitted: no fitted turn-on for q_eta=-1.969 pt=5.96 (ctr='', _sign2)`

ROOT's TRint caught it, the macro printed `[RUN] pbpb23 FillHistogramsCrossx completed
successfully!`, the job exited 0 — and left an **851-byte, 0-key corpse** at
`pbpb_2023/histograms_real_pairs_pbpb_2023_single_mu4_no_trg_plots_nominal.root`, byte-for-byte the
same failure that destroyed the 2024 file on 2026-09-06. **Nothing was lost**: the pre-change
backup holds the good 3.6 MB 2023 file, and 2025 was never reached.

**★ The repaired validation layer is what caught it.** `validate_root_file_quick` rejected the
corpse (`BAD RDF crossx yr23`) and the pipeline aborted before plotting. With the `-q` bug still in
place that check would have PASSED — it passed on anything — and Phase 3b would have plotted the
previous production's histograms as if they were the new selection.

**Root cause — NOT the q*eta binning, and NOT introduced by this change.** The q*eta bin exists:
the refit file carries `minus2_00_TO_minus1_50`, which contains -1.969, and its top bin is now
`_2_00_TO_2_20`, so **R1 Hazard 4 is cured**. The failure is the CENTRALITY: `ctr=''`.
`RDFBasedHistFillingData::FindCtrSuffix` returns `""` for `centrality >= 80` and for
`centrality < 0`, because the analysis' centrality binning is `ParamsSet::ctrbins =
{0,5,10,20,30,50,80}`. The empty suffix then reaches the trigger evaluator, which correctly
refuses to invent an efficiency and throws.

Measured on the FRESH trees (signal-region OS pairs):

| year | OS pairs | centrality >= 80 | centrality < 0 |
|---|---|---|---|
| 2023 | 434 966 | 438 (0.101 %) | 2 999 |
| 2024 | 306 628 | 317 (0.103 %) | 1 281 |
| 2025 | 869 072 | 805 (0.093 %) | 4 866 |

**Why it is latent, not new.** The surviving population after this change is a strict SUBSET of the
old one (every new cut removes pairs), so these pairs were always there. They never threw before
because a DIFFERENT throw fired earlier: every Pb+Pb muon with q*eta in [2.0,2.2) hit the retired
`_2_00_TO_2_30` fit key (Hazard 4). Curing that let the event loop run far enough to reach this.

**The fix, and why it is provably a no-op on every physics result.** These pairs contribute to
NOTHING as the code stands, established by direct inspection rather than assumption:
1. `FindCtrSuffix` returns `""` ⇒ they enter no `_ctr*` histogram.
2. The ONLY centrality-inclusive histograms in the Pb+Pb nominal output are
   `h3d_{op,ss}_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt`, whose **Z axis is exactly
   [0, 80] with edges {0,5,10,20,30,50,80}** ⇒ these pairs land in its over/underflow.
3. Every consumer excludes that over/underflow: `RAA_plotting.cxx:369,375` projects an explicit
   centrality bin RANGE (`bin_num_first, bin_num_last`), never `0..nbins+1`.
4. Counted on the pre-change 2023 output: **164 keys, 0 of which are centrality-inclusive without a
   centrality axis.** The template-fit output likewise: **60 keys, 0 without centrality binning.**

So excluding them up front removes the throw and changes no filled bin. Implemented in
`RDFBasedHistFillingPbPb::FillHistogramsCrossx` as a single filter on `df_op` and `df_ss`, with the
edges READ from `ParamsSet::ctrbins` (never retyped).

**Not treated as a stop-and-ask** because it is not a physics-results-bending choice: the
measurement's centrality acceptance is already 0-80 % by construction everywhere else, and the
filter only stops the evaluator being called on rows that have no home. Had any
centrality-inclusive histogram lacked a centrality axis, this WOULD have been a user decision.

### R9 — ★ The validation layer had never run, and what that cost

Found while fixing a NARROWER defect reviewer A reported (the pp crossx probe's failure message
being unreachable under `set -e`). The message was unreachable; the deeper truth is **the probe
never ran at all**.

> **`root -l -b -q` with no macro argument QUITS BEFORE READING STDIN.** The heredoc is never
> executed and the command exits 0 **unconditionally**. The `$(...)` capture form is equally dead:
> empty output, status 0.

So `validate_root_file_quick`, the post-hadd non-empty-tree check and the crossx
required-histogram probe **passed on anything** — a zombie file, zero keys, missing trees.
**19 occurrences across 8 pipelines.** The guard against the exact failure mode this whole task is
built around was itself inert.

Verified on the real artefact rather than argued: with the fix,
`validate_root_file_quick` REJECTS the 851-byte 0-key `pbpb_2024` corpse (rc=3) and ACCEPTS the
good 3.6 MB 2023 file. Before the fix, **both passed**.

Note the trap was already found and documented ONCE, in `run_dr_correction_fits.sh` on 2026-08-11
("`root -l -b -q` … QUITS BEFORE READING STDIN … this function returned 0 unconditionally").
Nobody swept the other seven pipelines. The lesson is not that the bug was subtle — it is that a
fix recorded in one file does not propagate itself.

**What it cost, and what it saved.** Cost: the whole of Phase 1 ran with no working artefact
validation, so every output had to be re-verified independently afterwards (done — see the Progress
Log; that is where the muon-pT numbers come from). Saved: the FIRST pipeline to run with the
repaired check, Pb+Pb crossx, immediately caught a swallowed exception (R8) and aborted before
plotting. With the `-q` bug still in place it would have passed and Phase 3b would have plotted the
previous production's histograms as the new selection.

Fixed in two commits because four of the eight pipelines were mid-run when it was found and
editing a script bash is still reading is unsafe: `47e7dea` (four idle) and `f047fd4` (the rest,
after they exited). Zero occurrences remain.

### R10 — POWHEG fullsim had not compiled since the all-vertex carve-out

All 22 `run_powheg_fullsim_wtruth_{bb,cc}` jobs wrote NOTHING, and **Condor recorded "Normal
termination (return value 0)" for every one** — ROOT exits 0 even when `.L` fails to COMPILE. Only
a freshness check on the part files revealed it: 0 of 10 were newer than the submission, while
POWHEG *truth* had all 6 fresh. A file-exists check would have seen ten healthy 100 MB parts from
March and reported success.

Two errors, both pre-existing:
1. `PowhegFullSimExtras.h:49` — `UseAllVertexIP() const` calls `getIsFullsimOverlay()`, which
   `PowhegAlgCoreT.h:89` declared **non-const**. The `PythiaAlgCoreT` twin has been const all along.
2. `PowhegFullSimExtras.h:63,75` — `ResonanceTaggingImpl` takes **4** arguments; both sites passed
   **3**. Never updated when `minv_cuts_to_use` was added.

**Why they survived:** `pp24_all_vertex_pairs.md` deliberately carved POWHEG fullsim out of its
blast radius, so `UseAllVertexIP` was added to a file that **nothing ever built**. R1 predicted this
task would drag POWHEG fullsim back in (a truth-pT cut reaches truth-level nodes); doing so is what
finally compiled it.

**The one judgement call**, flagged rather than buried: the 4th argument is
`self().pmsRef().minv_cuts_v2`, mirroring the Pythia fullsim twin EXACTLY. Determined rather than
chosen — the two pp-conditions fullsim samples are drawn on the SAME MC-vs-data figures and must
carry the same pair selection (this file's own header says so), and it is inert where it matters
since v1 and v2 differ only below 1.06 GeV while the signal region starts at 1.08. Reversible.

### R11 — ★ The MC trig-eff chain ran TIGHT-ONLY, and step4/sanity were stale for BOTH WPs

**Caught by the `polyn fit restriction` peer session, not by me, and not by either reviewer.** A
genuine unmet Done item (item 6, "MC trig-eff Steps 1-4"), and a good argument for the peer pings.

`pipeline_pythia_fullsim_pp.sh` Stage 10 honours `USE_TIGHT_WP`, whose **default is 1**. Job A was
launched as `ENABLE_MC_TRIG_EFF=1 ... full` without setting it, so the whole MC chain
(`FillMCTrigEffHists` step1 → `FitMCSinglesEffcy` → step3 → `plot_mc_trig_eff`) ran at **Tight
only**. Verified independently on disk before acting:

| product | Tight | Medium |
|---|---|---|
| `mc_trig_eff_hists_pp24_full[_medium_wp].root` | 09-10 06:00 ✅ | **09-07 21:20 ✗** |
| `..._step3.root` | 09-10 06:01 ✅ | **09-07 21:21 ✗** |
| `single_mu_effcy_pT_fit_mc[_medium_wp].root` | 09-10 06:00 ✅ | **09-07 21:20 ✗** |
| `..._step4.root` | **09-07 21:24 ✗** | **09-07 21:24 ✗** |
| `..._sanity.root` | **09-07 21:26 ✗** | **09-07 21:26 ✗** |

So the Medium MC sat on 4.0 GeV muons, the retired 8 GeV pair-pT axis AND the superseded
(−1.30,−1.05) gap window simultaneously — which makes the WP systematic meaningless, since
`plot_mc_trig_eff` with `use_tight_wp=false` compares that against the FRESH Medium DATA file
(refit 14:09). **A stale Medium file produces a comparison in which data and MC had different
selections applied, with no error and no warning** — which is exactly what
`run_data_trigeff_medium_wp.sh`'s own header warns about, for the data half of the same problem.

**Wider than reported:** the peer flagged Medium; checking the actual mtimes showed **step4 and
sanity were stale for Tight as well**, because Stage 10 runs only 10a-10d and stops at the plot.

**Concurrency constraint that shaped the fix.** The peer was at that moment READING the Tight
`_step3.root` for its ΔR refit. `run_mc_trigeff_round7.sh` would have rewritten it mid-read, so it
was NOT used. Instead: the full Medium chain (step1 → fit → step3 → step4 → sanity), which touches
only `*_medium_wp*`, plus the two missing Tight products (`_step4`, `_sanity`), which the peer does
not read. Tight `_step3` deliberately left untouched.

**CLOSED 2026-09-10 14:45-14:56.** Every MC trig-eff product is now same-day at BOTH working
points: Medium step1 14:45, fit 14:45, step3 14:46, step4 14:48, sanity 14:50, plots 14:56
(`mc_based_medium/`); Tight step4 14:52 and sanity 14:54 (its step1/fit/step3 stay at 06:00-06:01
from job A, deliberately untouched while the peer was reading them). No stage failed.

**A `USE_TIGHT_WP` default of 1 on a pipeline whose deliverable set includes both WPs is the real
defect** — the same shape as the data-side gap that A9 fixed in `run_pbpb_all.sh`. Recorded under
Remaining Work rather than fixed mid-rerun.

### R12 — ★★ BLOCKING: one ε_ΔR cell inflates the pp24 crossx weight 3-6x. NOT FIXED — user decision.

Reported by the `polyn fit restriction` peer session and **independently confirmed here in full**
(the fit report, the consuming code path, the screen that lets it through, and the yield impact).
It is in the NOMINAL `expo` method, not in that session's polyu change.

**The cell** — series `no_plateau_correction / expo / opposite sign / Tight`, i.e. exactly
`DrCorrCrossxMode()` / `DrCorrCrossxMethod()` / `DrCorrCrossxSign()`, line 94 of
`fit_report_opposite_sign.txt` (regenerated 2026-09-10 14:44):

| pair pT | pair eta | plateau | npts | A | λ | p | χ²/ndf | f(0) | C |
|---|---|---|---|---|---|---|---|---|---|
| [52.2,74.2) | [1.0,1.5) | 0.8225 | 20 | −3.3428 | **2.9998 AT LIMIT** | 1.3904 | 4.072 | 0.6327 | **3.9755 ± 1.1774** |

**Every other cell in that table has C ∈ [0.99, 1.27] with σ_C ≈ 0.003-0.05.** λ railed at its upper
limit 3.0, so over the fit domain the exponential is nearly a straight line, and MINUIT paid for it
by driving the free baseline to **4.8× the cell's own measured plateau**.

**The delivered correction is ε_ΔR = f/C**, so in this cell it is **0.159 at ΔR=0 and 0.323 at
ΔR=1** — it never approaches 1, which it must by construction.

**It IS in the cross-section.** `RDFBasedHistFillingPP.cxx:387-392`:
`eps_trig^pair = eps^nc_1 · eps^nc_2 · eps_dR`, then `w_trig = 1/eps_trig`. (This supersedes the
"not propagated" state recorded earlier.)

**Measured impact on the figures produced today** (pp24 OS, minv ∈ (1.08,2.9), pair pT > 9):

| quantity | value |
|---|---|
| pairs in the affected cell | **205** |
| ...as a fraction of the whole signal region | 0.032 % |
| ...as a fraction of the [52.2,74.2) GeV pT slice | **15.6 %** |
| ...with ΔR < 1, i.e. actually corrected | **205 (100 %)** |
| weight inflation per affected pair | **3.1× - 6.3×** |

⇒ globally negligible, **locally severe**: the (52-74 GeV, η 1.0-1.5) panel is inflated 3-6×
outright, the η-integrated pair-pT spectrum in that bin by roughly 1.3-1.8×, and **R_AA inherits it
as a SUPPRESSION** there because pp is the denominator.

**Why nothing caught it — three independent screens, all blind to this failure:**
1. `fit_ok` carries **no χ² term**, so χ²/ndf = 4.072 passes.
2. `DrCorrPlateauUsable(plateau, err)` bounds the baseline only **from BELOW**
   (`plateau >= 0.5 && err < plateau`). C = 3.98 with σ_C = 1.18 satisfies both. Its own comment
   shows why: it was written against a *collapsing* baseline ("dividing by 0.01 inflates the curve
   100x") and is structurally blind to a runaway one.
3. `DrCorrectionCrossxEvaluator::Eval` applies **no cap** to the returned f/C.
Plus: the 2026-09-07 shape restriction constrained A and p but left λ ∈ [0.02, 3.0] free, and
railing λ is the mechanism — a cell with no plateau within reach buys a flat-looking fit by
inflating C.

**RESOLVED 2026-09-10 (user chose: bound C against the cell's own measured plateau).**
Implemented as `DrCorrBaselineConsistent(C, measured_plateau)` in `dr_correction_sample_cfg.h`, a
two-sided factor-1.5 window, wired into BOTH acceptance sites of `dr_correction_apply.h`. The
per-cell measured plateau was already persisted as `h_step3_plateau` in every fit file, so no new
input was needed. Chosen over the two alternatives on measurement, not preference:

| candidate | why not |
|---|---|
| χ² term in `fit_ok` | **poor discriminator here**: 5 cells exceed the bad cell's χ²/ndf = 4.072 and **4 of them have a perfectly good baseline** (C/plateau 0.75-1.05). The worst χ² in the table, 13.986, has C/plateau = 0.75 and is fine. Any threshold catching the bad cell rejects several good ones |
| disqualify `AT LIMIT: lambda` | clean today — exactly 1 cell has it, the bad one — but rests on a single instance. (16 cells have `AT LIMIT: p` with C ∈ [0.947,1.195], all healthy, so p-at-limit must NOT disqualify) |

**Verified end to end, not argued:**
- The cascade census moves **70 expo / 0 polyu / 0 interp / 2 raw → 61 expo / 1 polyu / 1 raw.**
  The polynomial tier — which exists for exactly this and was being reached in **zero** cells
  because the expo screen accepted almost everything — now takes the bad cell.
- Exactly **two** cells are rejected, both predicted: ratio 4.834 and 1.991.
- In the pathological cell ε_ΔR now runs **0.56 at ΔR=0.05 → 1.0 by ΔR=0.5**, i.e. it tends to 1 as
  it must, matching its healthy neighbours; the weight factor falls from **3.1-6.3× to 1.0-2.0×**.

**STILL OPEN, and explicitly NOT covered by this fix** — a bound on the DELIVERED correction f/C is
a different question from a bound on C. Cells remain whose C is entirely normal (C/plateau ≈ 1.03)
but whose f/C dips low at small ΔR; the extremum over the accepted cells is 0.0363 at
(pair-pT bin 7, pair-eta bin 7). **Important nuance measured here:** that 0.0363 is the value at
ΔR **exactly 0**, which is kinematically unreachable (ΔR ≳ 2m/pT ≈ 0.014 at 150 GeV). At the
smallest *reachable* ΔR the same cell gives ≈0.27, and it rises to 0.98 by ΔR = 0.3 — the right
shape. So the delivered-magnitude question is real but much less severe than the raw extremum
suggests, and it is entangled with reviewer finding R13-3 (the barrel high-pT correction). Not
resolved here.

**Superseded status line:** The remedy is a physics choice with three candidate
forms (add a χ² term to `fit_ok`; bound C from ABOVE against the cell's own measured plateau, which
is already in the file; or treat `AT LIMIT: lambda` as disqualifying — any of which routes this cell
to the polyu or interp tier). Under the Autonomy Contract this is a physics-results-bending
ambiguity ⇒ **stop and ask**. **Every pp24 crossx and R_AA figure produced on 2026-09-10 carries
this defect** and must be regenerated after the fix.

### R13 — `/review-plot`: two reviewers, both FAIL, and what remains open

Disjoint scopes, read-only, physics-results criteria C1-C4 applied. Both independently reached the
ε_ΔV pathology of R12 from the FIGURES, without being told, which is a good sign for the criteria.

**Confirmed the R12 defect was visible in the final results** (reviewer 1): η-integrated, the pp24
bin [52.23,62.27] sat **×1.36, +7.2σ** above a weighted local power law through its neighbours
(every other bin within ±2.4σ). Split by η panel it was **×3.50 / +11.1σ** in η ∈ [1.0,1.5) — the
R12 cell exactly — and it propagated into R_AA as a **coherent dip in ALL SIX centrality bins**.
Reviewer 2 saw the same edge independently in the MC/data ratio: flat at 1.30-1.45 for all ten bins
from 9 to 52.23 GeV, then 0.92 in [52.23,62.27] — an **8σ** step, with both generators moving
together, i.e. the discontinuity is in the DATA denominator.

**★ MEASURED AFTER THE FIX (same local-power-law test the reviewer used, on the refilled
histogram):** the R12 cell's panel is fully corrected and the residual is entirely the second cell.

| panel | bin 11 [52.23,62.27) pull BEFORE | AFTER |
|---|---|---|
| η ∈ [1.0,1.5) — the R12 cell | **+11.1σ** (×3.50) | **+0.9σ** — every bin now within ±0.9σ |
| η-integrated | **+7.2σ** (×1.36) | **+3.0σ** |
| η ∈ [0.5,1.0) — the SECOND cell | +4.1σ | **+3.2σ** (unchanged, exactly as predicted) |

So R12 accounted for the dominant part of the discontinuity, and what remains is item 1 below,
whose C/plateau is a perfectly normal 1.098 — the screen correctly leaves it alone, which is why
it must be diagnosed separately rather than by loosening this one.

**★ STILL OPEN — findings the R12 fix does NOT explain, and that need their own investigation:**
1. **A second bad cell.** With η ∈ [1.0,1.5) removed entirely, the [52.23,62.27] bin is *still*
   ×1.16 / +3.4σ high, carried by **η ∈ [0.5,1.0) at ×1.53 / +4.1σ** in the same pT bin (its
   neighbour is ×0.80/−2.6σ, so the step across 52.23 GeV in that panel is ×1.9). Its C/plateau is
   normal, so R12's screen does not touch it.
2. **A systematic barrel high-pT blow-up.** The pp total trigger correction reaches 8.2/7.0/9.0 in
   the three central η panels at 74-105 GeV and 9.6/7.7/10.4 at 105-150, against 1.9-2.4 in the
   forward panels. Barrel/endcap is ~1.7 at low pT (consistent with RPC-vs-TGC L1 acceptance) and
   grows to ~5 at high pT, which geometric acceptance does not explain. **Pb+Pb is immune** (flat
   1.55-2.36 everywhere — it applies the bare trigger union, no ε_ΔR), so R_AA inherits this
   entirely from the pp denominator.
3. **Step-3 ε_ΔR plateaus sit at 1.10-1.27 in the two highest pair-pT columns** (all 16 cells),
   where they are 1±0.05 in the 56 lower-pT cells. Since |value−1| feeds the systematic, that is a
   10-27 % systematic confined to the top pT cells — the same region as 1 and 2.

These three are one region and plausibly one cause. **Not investigated here** — they are
pre-existing, they are not what this task set out to change, and the honest position is that the
pp24 high-pair-pT cross-section and R_AA above ~50 GeV remain PROVISIONAL until they are settled.

**Fixed from the review:** the Tight Step-4 plots were drawn at 06:03 from the pre-rerun Step-4
histogram (Stage 10 stops at the plot, so the histogram I regenerated at 14:52 was never redrawn) —
they carried the RETIRED 8→150 coarse binning while the Medium twin, run later, was on the
canonical 9→150. Two coexisting pair-pT binnings in one plot set, the exact BLOCKING failure.
Redrawn. Also redrawn: two `sanity_check_crossx` figures from 2026-08-24 still showing ±2.4 η panels.

**Reviewer findings deliberately NOT acted on in this task** (cosmetic or pre-existing; recorded so
they are not lost): legend text overrunning the frame in `DrawPairPtByEtaWithDrLines`; R_AA drawn on
a LINEAR x-axis although the pair-pT axis is log-binned; R_AA markers all black because
`SetMarkerColor` is never called; Pb+Pb `counts/` figures labelling a pair count `N_{events}`;
no WP config var in `RAA_plotting.cxx` / `plot_crossx_trig_corr_sanity.C`; `label_line2_ = "tight WP"`
set but never drawn on the Pb+Pb figures; 48 retired-binning PNGs from 09-08 still on disk; ~20
blank canvases in the pp and Pb+Pb pair-trigger sets; the `mc_data_compr` READMEs quoting
pre-adoption numbers; the Pb+Pb turn-on canvases labelled "2mu4" when Pb+Pb runs single mu4.

**Both reviewers independently PASSED** the things this task actually changed: axis provenance
(nominal 150 vs opt-in 120 never crossed, 2:1 nesting exact), no surviving "> 8 GeV" or "4 GeV"
label anywhere, the R_AA group labels matching the axis they were drawn from, the Pb+Pb 48-bin
migration visible on disk with the gap bands at η≈0 and ±1.15, the per-histogram Pb+Pb year label
being honest, PNG-only, `Scale(N,"width")` on every differential quantity, and — notably — that the
old `w_trig = 0` sentinel pathology is GONE: 1/ε ≥ 1 in **every** bin of pp and Pb+Pb.
Run 2 cross-check: pp barrel plateau 0.68-0.76 / endcap 0.84-0.95 against Run 2 mu4 ≈0.70/≈0.90,
and the Step-1 MC/data ≈1.10-1.15 reproduces the known Run 2 barrel data/MC deficit — consistent.

### R6 — Cross-session state (three peers share this checkout)

Four peer sessions exist; three were engaged and all cleared this rerun. **They share the working
tree**, so commit by explicit path and never `git add -A`.

- **`polyn fit restriction`** — owns the Step-3 `polyu_fixedRp` reparametrization (`d23fb81`):
  `A ≡ f(0) − C` is now parameter 0 with `A ∈ [−50, 0]`, Step 3 ONLY, Step 4 unchanged.
  **Binding constraints on this rerun:** (i) `plot_dr_correction_fits.cxx::LoadFunc` THROWS if a
  polyu TF1's parameter 0 is not named `A` — so **refit polyu wherever you replot it, and never
  replot a polyu file you did not just produce**; that throw means the file predates `d23fb81`,
  not that the code is broken; (ii) **refit expo + polyu + interp TOGETHER**, never one alone;
  (iii) the four-approach χ²/ndof ranking in `mc_trigeff_dr_binning_approaches.md` is UNUSABLE
  until after this rerun, and `run_mc_trigeff_closure.sh` currently dies in Stage 2 on the binning
  guard. Correct order afterwards: refill Step 3 → refit all three methods → rerun closure →
  re-read the ranking. Wants a ping when the fresh trees and Step-3 histograms exist.
- **`additional fullsim statistics`** — read-only on
  `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root`; writes only to a new
  `plots/pp_trigger_efficiency/mc_pthat_slice_mass_stats{,_medium}/`. Sizing a NEW MC production
  request; its current numbers are on the SUPERSEDED selection and it will redo on the fresh
  trees. Found D12. Wants the same ping.
- **`pT > 4.5GeV cut`** — the diagnostic thread this task executes the adopt branch of. Confirmed
  consistent, closed its doc (`85ac4cb`), idle. Nothing outstanding.

**A disagreement worth preserving, because the resolution is a reusable rule.** This doc first
claimed the whole 17:11→18:52 chain was a chimera. The peer disproved that for the trig-eff chain
by measurement (`mc_trig_eff_hists_pp24_full.root` 18:46: `h_mc_pt_denom_mu±` axis **4.50**,
**underflow = 0**; the 09-07 Medium control: axis **4.00**) — `FillMCTrigEffHists` re-applies
`pt > 4.5` textually, and a 4.0 GeV NTuple is a strict SUPERSET, so the re-applied cut reproduces
the 4.5 population exactly. The claim survives only for the ε_reco side. That is the R4 rule, and
the peer's own per-file grep formulation had to be corrected too: the test is **"does the
selection string reaching the Filter contain the cut, AFTER resolving shared helpers"** — a
per-file grep reports the closure and pair-eff as stale when they are immune via
`MCTrigEffPairSelection.h`.

**Unexplained and worth knowing:** something ran the complete 17:11→18:52 MC chain on 2026-09-08,
ending in `pair_reco_eff_pp24_full.root` at 18:52, using the then-partially-edited `ParamsSet.h`.
All three peers and this session deny it; `condor_q` was empty and no process survived. Best
hypothesis: the abandoned step 7a of `pp24_all_vertex_pairs.md`, whose Latest Stage still reads
"MC half NOW RUNNING" and whose chain matches exactly. Everything it wrote is regenerated here.

### R2 — Step 3 DONE: `ParamsSet.h` (the keystone)

Gap window `{-1.30,-1.05}` → `{-1.25,-1.05}`; `pTbins` low edge 4 → 4.5; `pT_bins_40`,
`pT_bins_8` → 4.5; `pT_bins_80` → 9; `pT_bins_120` → **16 log 9→120** (alternative);
`pT_bins_150` → **16 log 9→150** (new default); `pair_pt_coarse_bins` and its 4-bin variant →
9→150; `single_mu_pt_coarse_bins` → {4.5, 8, 14, 25, 100}; `pt_titles[0]` relabelled.
Comment blocks updated for truthfulness, including marking the measured gap-cut cost and ε_acc
tables as STALE IN TWO WAYS (superseded window AND measured at pT > 4).

## Remaining Work

### THE RERUN — status 2026-09-10

**Phase 1 (NTuple) — DONE.** 97 Condor jobs + the ~7 h local fullsim pass, plus a 22-job POWHEG
fullsim resubmission after R10. All eight data trees independently verified for the muon-pT cut
(see Progress Log): nominal at exactly 4.5, trig-eff mode at exactly 4.0 by design, zero violations.

**Phase 2 (derived corrections) — DONE.**
- 2a data turn-on refits: pp + Pb+Pb, Tight AND Medium. Cured R1 Hazard 4.
- 2b MC trig-eff Steps 1-4 + sanity, **both WPs** (Tight-only gap closed — R11).
- 2c ΔR-correction fits: **owned by the `polyn fit restriction` peer session**, running now, all
  three methods together (expo / polyu_fixedRp / interp) per R6.
- 2d MC closure: follows the peer's refit.
- 2e `pair_reco_eff_pp24_full.root` rebuilt and sanity-checked (all values in [0,1], 143 filled
  3D cells).

**Phase 3 (final results) — DONE except the review.**
- 3a pp24 crossx + pp24 nominal NTuple: DONE. 3b Pb+Pb crossx 23/24/25: DONE (after R8).
- 3c R_AA modes 1/2/3 on run_year_trigger_mode 6: DONE. Both existing guards passed, which
  independently confirms pp and Pb+Pb are on the SAME 48-bin pair-eta axis and the same fine
  pair-pT axis — the whole point of D1.
- MC-vs-data comparison (3 macros, 15 PNGs), signal acceptance + cutflow (Pythia + POWHEG, nominal
  + `_pt_120`), crossx sanity plots incl. the Pb+Pb panels: DONE.
- **3e `/review-plot` IN FLIGHT** — two read-only reviewers, disjoint scopes (crossx + R_AA;
  MC-data + acceptance + trig-eff), applying physics-results criteria C1-C4.
- `/sync-note-figures` + `/check-note-sync`: NOT started.

### Open items for the user

1. **POWHEG `cc` truth is now available for the first time.** All 6 `cc` parts were produced this
   round (488 513 / 390 784 / 391 307 / 488 205 / 448 918 / 233 562 OS entries); previously parts 2
   and 6 were missing, which is why `PlotMCDataComprBaseClass.h` states the curve is bb-ONLY and the
   legend must not imply otherwise. **I did NOT merge it**, so `RDFBasedHistFillingPowheg` logged
   `Skipping missing input file: ...cc_truth.root` and the bb-only configuration is preserved
   EXACTLY. This matters because that class loops over `{bb, cc}` and includes whatever EXISTS — so
   merging the cc file would silently turn the POWHEG curve into bb+cc while its legend still read
   "POWHEG bb". Adding cc would change a final-results figure and needs a decision. (Impact on the
   comparison curve itself is likely small: its histograms filter on `from_same_b`, which cc pairs
   fail; the flavour-binned families would change a lot.)
2. **Old `_pt_150` plot directories.** Five plainly superseded (D5 authorises removal now that
   Phase 3 has regenerated the nominal dirs) and five DATED SNAPSHOTS the subagent flagged as
   "keep or decide separately". 7.1 MB total, PNGs only. Not deleted.
3. **`USE_TIGHT_WP` defaults to 1 in `pipeline_pythia_fullsim_pp.sh`** while its deliverable set
   includes both working points — the defect behind R11, and the same shape as the data-side gap
   A9 fixed in `run_pbpb_all.sh`. Not changed mid-rerun.

### Constraints that bind the rerun (do not rediscover these)

- **Validate every stage by FRESHNESS + a required histogram, never by file-exists** — and note
  R9: the helper that was supposed to do this had never actually run.
- **`root -l -b -q` with a HEREDOC executes NOTHING and exits 0.** Never use `-q` when the macro
  comes from stdin. 19 sites fixed; zero remain.
- **ROOT exits 0 even when `.L` fails to COMPILE** (R10) — a Condor "Normal termination (return
  value 0)" is not evidence a job did anything.
- **polyu:** refit wherever you replot; never replot a polyu file you did not just produce.
- **Never refit only one ΔR method.**
- Backups: `~/usatlasdata/dimuon_data/pre_mupt45_backup_20260908/` (166 files, 21.7 GB). They
  earned their keep: the good pre-change Pb+Pb 2023 crossx file was needed when R8 left a corpse.

### Carried forward, NOT part of this task

- `SingleBAnalysis/SingleBAnalysisBase.cxx` (legacy pre-RDF) still retypes the fine pair-pT axis
  and its own `signal_cuts`; not in the active chain.
- Dead NTuple copies (`PythiaNTupleFirstPass*`, `original_no_template_class/*`,
  `maybe_old_unsure_pbpb/*`) still carry 4.0 GeV — reviving any revives it.
- `plot_reco_distr_singleb_vs_op_pp24.C` keeps its retired selection but is now gated behind
  `RUN_RECO_DISTR=1` with a STALE banner instead of running as a nominal pipeline stage.
- **Open INFO items needing a user decision** (deliberately not acted on): `pT_bins_80` is reused
  as a SINGLE-muon axis in `var1D_pythia_truth.json` / `var1D_powheg_truth.json` while now
  starting at 9; retyped reco-eff ΔR slices `{8,12,20,∞}`; the
  `trigger_effcy_calc`/`pbpb_run3_mu4_force_nominal` coupling (verified correct in all eight live
  Pb+Pb scripts, but the flag silently gained a selection meaning); and **the muon-pT provenance
  stamp** — proposed and NOT adopted, so the next threshold change has the same silent exposure
  R4 describes.

## Latest Stage

**THE RERUN IS COMPLETE. Phases 1-3 done and regenerated on the corrected ε_ΔR. Two physics
questions remain OPEN and are the user's to settle; until they are, the pp24 cross-section and R_AA
ABOVE ~50 GeV are PROVISIONAL.**

State as of 2026-09-10 16:12:
- **Phase 1** (97 Condor jobs + the ~7 h local fullsim pass + a 22-job POWHEG resubmit) — done, and
  the cut verified directly in all eight trees, not inferred: nominal at exactly 4.5, trig-eff mode
  at exactly 4.0 by design (D8), zero violations either way.
- **Phase 2** — data turn-on refits pp + Pb+Pb, Tight AND Medium; MC trig-eff Steps 1-4 + sanity at
  BOTH working points; `pair_reco_eff_pp24_full.root` rebuilt and sanity-checked;
  `pair_trig_eff_pp24_full.root` regenerated on the 9 GeV axis.
- **Phase 3** — pp24 and Pb+Pb cross-sections, R_AA (3 modes), MC-vs-data (3 macros), signal
  acceptance + cutflow, all sanity plots. All regenerated AFTER the ε_ΔR fix.
- **`/review-plot` ran** (2 reviewers, both FAIL) and everything actionable within this task's scope
  is fixed; the residue is listed in R13 as recorded-not-acted-on.
- **ΔR fits / MC closure** remain the `polyn fit restriction` peer's; it has been unblocked twice
  (fresh trees, then `pair_trig_eff_pp24_full.root`) and owes the four-approach χ²/ndof ranking.

### What a resuming agent must NOT assume

1. **That the high-pair-pT pp results are final.** R13 lists three OPEN findings in one region — a
   second bad ε_ΔR cell in η ∈ [0.5,1.0) (+3.2σ at 52.23 GeV), a barrel high-pT correction blow-up
   that geometric acceptance does not explain, and Step-3 plateaus at 1.10-1.27 in the top two pT
   columns. Plausibly one cause. Pb+Pb is immune, so R_AA inherits all of it from the pp denominator.
2. **That "it compiled" or "the job exited 0" means anything**, unless a freshness check says so.
   This rerun found: a validation layer that had never executed (R9), 22 Condor jobs reporting
   normal termination while writing nothing (R10), and an RDF event loop that threw, printed
   "completed successfully", exited 0 and left an 851-byte corpse (R8).
3. **That a reviewer's or peer's premise is true.** Three were checked and wrong: a peer reported
   the pair-eff sources as uncommitted (all three are clean at HEAD); a reviewer's suggested
   rc-capture idiom silently discards the exit code; and my own guessed filename matched no file.
   Every one was caught by testing rather than by reading.

### Open, needing the user

1. **A bound on the DELIVERED ε_ΔR (f/C)**, as opposed to the bound on C that was adopted — the two
   remaining findings above may or may not be one problem with it. See R12/R13.
2. **POWHEG `cc` truth is now producible** for the first time; deliberately NOT merged, so the curve
   stays bb-only. Merging it would silently change a final figure under an unchanged legend.
3. The old `_pt_150` plot directories; the `USE_TIGHT_WP` default that caused R11.
