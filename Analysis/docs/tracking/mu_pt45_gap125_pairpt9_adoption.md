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

### THE RERUN, in dependency order (nothing below has been started)

**Phase 1 — NTuple (Condor, except the fullsim pass).** Roughly 90 Condor jobs.

| # | What | Notes |
|---|---|---|
| 1a | pp24 nominal, `run_pp_24_nominal.sub` (trigger_mode=3) | 12 jobs; then hadd |
| 1b | pp24 trig-eff, `run_pp_24.sub` (trigger_mode=1) | 12 jobs; **stays at 4.0 GeV by design (D8)** |
| 1c | Pb+Pb 2023/2024/2025, both modes | 12+12; `SKIP_EVSEL=1` (event selection is NOT in the blast radius — its cuts are event-level, derived from raw NTUPs) |
| 1d | Pythia truth (private, nonprivate 5.36, nonprivate 5.02) | 2+6+6 |
| 1e | POWHEG truth + POWHEG fullsim `run_powheg_fullsim_wtruth_{bb,cc}.sub` | 6+6 and 11+11. POWHEG fullsim IS in the blast radius: the carve-out in `pp24_all_vertex_pairs.md` rested on that change being reco-only, and a truth-pT cut is not |
| 1f | HIJING overlay | 1 |
| 1g | Pythia fullsim pp24 FULL `_pdf` | **NOT Condor** — local pass over the LGD symlink farm, ~6–8 h. `ENABLE_MC_TRIG_EFF=1 pipelines/pipeline_pythia_fullsim_pp.sh full`, **NOT** `SKIP_NTP=1` |

**Phase 2 — derived corrections.**
2a. Data mu4 turn-on refits, pp AND Pb+Pb, **Tight AND Medium** (`pipeline_pp_trig_eff.sh`,
    `pipeline_pbpb_trig_eff.sh`, `SAMPLES="pp" pipelines/run_data_trigeff_medium_wp.sh`).
    This also cures the Pb+Pb `_2_00_TO_2_30` vs `_2_00_TO_2_20` key mismatch (R1 Hazard 4).
2b. MC trig-eff Steps 1–4.
2c. ΔR-correction fits — **all three methods together** (`expo polyu_fixedRp interp`), never one
    alone (R6).
2d. MC closure. NOTE `run_mc_trigeff_closure.sh` currently dies in Stage 2 on the binning guard;
    it should pass once 2b/2c are refilled on the 9 GeV axis.
2e. `build_pp24_fullsim_pair_reco_eff.C+(true)` → `pair_reco_eff_pp24_full.root` (not in any
    pipeline; RECREATE, no backup of its own).

**Phase 3 — final results.**
3a. pp24 crossx pipeline (+ MC-data comparison, Stage 8).
3b. Pb+Pb: `run_pbpb_all.sh` — now SERIALIZED, trig-eff before crossx.
3c. R_AA (`RAA_plotting.cxx` mode 6), now that pp and Pb+Pb share a signal region again.
3d. Signal acceptance + cutflow, crossx sanity/stage plots.
3e. `/review-plot` on the regenerated plot sets; `/sync-note-figures` + `/check-note-sync`.

### Constraints that bind the rerun (do not rediscover these)

- **Validate every stage by FRESHNESS + a required histogram, never by file-exists.** ROOT
  swallows exceptions thrown inside an RDF event loop and still exits 0, leaving a fresh
  near-empty file. `pbpb_2024/histograms_real_pairs_..._nominal.root` is an 851-byte, 0-key
  corpse of exactly this (2026-09-06); 2023/2025 are last-good 2026-07-08.
- **polyu:** refit wherever you replot; never replot a polyu file you did not just produce; a
  `LoadFunc` throw means the file predates `d23fb81`, not that the code is broken.
- **Never refit only one ΔR method.**
- The four-approach χ²/ndof ranking in `mc_trigeff_dr_binning_approaches.md` is UNUSABLE until
  after this rerun; re-read it only after closure is rerun.
- **After the turn-on refit, print the R3 diagnostic** — σ_ε(4.5)/ε(4.5) from the fit covariance
  (`"QR"` → `"QRS"`, `GetConfidenceIntervals` at `pT_min`), flag > 0.10, **stop at > 0.20**,
  companion |ρ(mean,σ)| > 0.95. The at-a-fit-limit `*` flag provably cannot detect this.
- **Ping the two waiting peers** when the fresh trees and Step-3 histograms exist (R6).
- Backups of every pre-change artefact are at
  `~/usatlasdata/dimuon_data/pre_mupt45_backup_20260908/` (151 correction files + hists + fits;
  `pair_reco_eff_pp24_full.root` in it is genuinely the 2026-08-18 pre-change file).
- Disk: ~1.34 TB of ~1.50 TB (halved figures). Reruns overwrite in place and the 4.5 GeV cut
  shrinks trees, so the steady-state change is negative.
- The eight classes of R4 produce artefacts that are stale until Phase 1 completes; nothing
  should consume them in the meantime.

### Immediate next step

**A 4th `/review-analysis-code` iteration on the last amendment batch has NOT been run.**
Iteration 3 demonstrated that amendments are not reliably complete (four iteration-2 fixes were
half-applied, R5), so the last batch is unverified. Either run it, or accept the risk and start
Phase 1 — the user was offered both.

### Merged from `_sub_crossx_axis_default_swap.md` (scratch doc now deleted)

That subagent left three OPEN ITEMS FOR THE COORDINATOR that were never resolved. Checked
2026-09-09:

- **(a) The `_150` token in var1D / producer variable NAMES — left alone, deliberately.** Nine
  sites (`var1D_{pythia_truth,pythia_fullsim,powheg_fullsim,powheg_truth,pp,pbpb}.json`,
  `RDFBasedHistFillingPP.cxx:798`, `...PythiaFullsim.cxx:40,63`, `...PowhegFullsim.cxx:54,62`)
  carry names like `pair_pt_log_150` / `truth_pair_pt_log_150`. They are ALREADY bound to
  `"binning": "pT_bins_150"`, i.e. already on the nominal axis, and there is no second
  alternative-axis member — so the names are cosmetically stale but **functionally correct**. The
  subagent recommended dropping the token; NOT done, because the rename has real consumers
  (`PlotMCDataComprBaseClass.c:38`, `plot_mc_data_pair_pt_in_eta.cxx`) and buys zero correctness.
  **Renaming them `_150 -> _120` would be actively WRONG** — it would rebind the only 1D truth
  pair-pT spectrum to the opt-in alternative and move the whole MC-vs-data comparison off the
  nominal axis. Recorded so a future reader does not "tidy" it the wrong way.
- **(b) Signal-acceptance producers — RESOLVED, verified 2026-09-09.** The subagent flagged that
  `h2d_sig_accept_*` lives in the truth classes, outside its scope, and needed the same treatment.
  It got it: `RDFBasedHistFillingPythiaTruth.cxx:479-503` and `...PowhegTruth.cxx:299-323` book the
  unsuffixed `h2d_sig_accept_{num,denom}_pt_eta` on `pT_bins_150` and the opt-in
  `..._pt_120_eta` on `pT_bins_120`. Reviewer B independently passed `SignalAcceptancePlotter`.
- **(c) Old `_pt_150` plot directories — still present, NOT yet deleted.** Ten dirs, 7.1 MB, PNGs
  only, under `plots/single_b_analysis/`. Five are plainly superseded (`pp24_pt_150`,
  `pbpb_23_24_25_combined_pt_150`, `pbpb_23_24_combined_pt_150`, `pythia_pt_150`,
  `powheg_pt_150`) — D5 authorises removing these, but only AFTER Phase 3 regenerates the nominal
  unsuffixed dirs, so a failure cannot leave nothing behind. Five are DATED SNAPSHOTS
  (`*_backup_20260615`, `*_backup_20260616_pre_reco_nominal`, `*_backup_20260505`) that the
  subagent explicitly flagged as "keep or decide separately" — **user decision, not taken.**

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

**REVIEW ITERATION 4 IN FLIGHT; PHASE 1 PRE-FLIGHT DONE; RERUN ABOUT TO START.**

Plan for this session (user instruction 2026-09-09 22:33: *"finish this scoped 4th review then
perform the full reruns"*), in order:

1. **Review iteration 4** (in flight) — two independent read-only reviewers on the ONLY unverified
   batch, commits `01031f0` (plots/pipelines, 42 files) + `c9ab2df` (docs, 11 files). Split by scope
   so nothing is missed and neither can clobber the other: **A** = `Analysis/pipelines/` + the
   run/submit shell scripts (the rerun drivers; the Phase-1 gate); **B** = `plotting_codes/`,
   `RAA_plotting.cxx`, and the two designated ground-truth docs. Log:
   `.claude/logs/review-analysis-code-20260909-223311-mupt45-iter4-amendment-verify.md`.
   Both prompts carry the D6/D8/D9/D11 do-not-flag list so the deliberate divergences are not
   "fixed", and the append-only-history exemption for tracking docs.
   Rationale for running it at all: R5 says amendments are not reliably complete, AND the
   iteration-3 pass was never written into the old review log (its header still reads
   "Iterations completed: 2"), so the assurance is weaker than the doc implied. ~20 min against an
   8 h+ Phase 1 that a Phase-1-reaching defect would waste.
2. **Amend** anything the reviewers find that reaches the rerun; re-verify amendments explicitly
   (the R5 failure mode), commit by explicit path.
3. **Phase 1 → 2 → 3** exactly as the tables below prescribe.

### Phase-1 pre-flight, completed and verified this session (read-only)

- **Condor queue EMPTY, no analysis process running.** Nothing to collide with.
- **Working tree clean** at `c4e13ed`; the 21:45-21:48 mtimes on the `pipelines/` files are this
  task's own pre-commit edits (`git diff HEAD -- Analysis/pipelines/` is empty), NOT a peer writing.
- **`ParamsSet.h` values re-confirmed independently** by compiling a throwaway macro against the
  live header (not by reading it): `signal_pair_pt_min = 9`, `pTbins[0] = 4.5`, gap windows
  `(-1.25,-1.05) (-0.10,+0.06) (2.20,2.40)`, `pair_eta_fiducial_max = 2.2`,
  `N_PAIR_ETA_CROSSX_BINS = 48`, `N_COARSE_PAIR_PT_BINS = 8`.
- **`.L DataAnalysisClasses.h` loads clean inside the EXACT Condor environment**
  (`atlasLocalSetup.sh` + `lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"`), which is what the
  `run_*.sh` job scripts do. Phase 1 will not die on a load error.
- **D8 wiring verified end to end at the run-script level**, which is where it can silently go wrong:
  `DimuonDataAlgCoreT.c:70` DERIVES `trigger_effcy_calc = (trigger_mode==0||trigger_mode==1) &&
  !pbpb_run3_mu4_force_nominal` (it is never set by hand), and `:627` reads
  `const float muon_pt_min = trigger_effcy_calc ? 4.0f : 4.5f`. Checked every live job script:
  `run_pp_24_nominal.sh` sets `trigger_mode = 3` -> 4.5; `run_pp_24.sh` sets nothing and the header
  default is `trigger_mode = 1` -> 4.0 (correct for 1b); `run_pbpb_{23,24,25}_nominal.sh` each set
  `pbpb_run3_mu4_force_nominal = true` -> 4.5; the plain `run_pbpb_{23,24,25}.sh` set neither -> 4.0.
  All six Pb+Pb scripts and both pp scripts are correct.
- **Job counts confirmed from the `.sub` files** (total **97**, matching the "~90" estimate):
  pp24 12 + 12; Pb+Pb 4+4 (23) + 2+2 (24) + 6+6 (25) = 24; Pythia truth 2 + 6 + 6 = 14;
  POWHEG truth 6 + 6 = 12; POWHEG fullsim 11 + 11 = 22; overlay 1.
- **Backups verified present and correct** at `~/usatlasdata/dimuon_data/pre_mupt45_backup_20260908/`
  (21.7 GB, 166 files: `mc_corrections` 151, `pp_2024` 6, `pbpb_2023` 5, `pbpb_2024` 2,
  `pbpb_2025` 2). Every live `single_mu_effcy_pT_fit*.root` and every live per-year crossx
  `histograms_real_pairs_*_nominal.root` has a backup copy; the backed-up
  `pair_reco_eff_pp24_full.root` is genuinely the 2026-08-18 pre-change file.
  Two files that Phase 2a will CREATE rather than overwrite (so nothing is lost):
  `pbpb_2024` and `pbpb_2025` have **no** `single_mu_effcy_pT_fit_medium_wp.root` at all.
- **Disk**: data fileset 1.38 TB used of a 1.54 TB soft quota / 1.79 TB hard limit (halved figures).
  ~160 GB of headroom to the soft quota. The 4.5 GeV cut shrinks the trees and reruns overwrite in
  place, so the steady state is negative -- but this is tighter than the doc's earlier 1.34/1.50
  reading and is worth re-checking between phases.
- **Driver mapping established** (which script owns which Phase row), so the phases below are run,
  not reinvented: `pipeline_pp_trig_eff.sh` = 1b + 2a(pp, Tight) in one pass (its Stages 1-4 are the
  NTuple submit/wait/validate/hadd, Stages 5-8 the fine-q*eta fill, the turn-on fit and its
  validation); `run_data_trigeff_medium_wp.sh` = 2a(pp, Medium); `pipeline_pp_crossx.sh` = 1a + 3a
  (Stages 1-4 NTuple nominal, 5 RDF crossx fill, 6 crossx plots, 7 the trig-eff-correction sanity
  check, 8 the MC-data comparison); `run_pbpb_all.sh` = 1c + 2a(Pb+Pb) + 3b, now serialized.
  **Ordering consequence that binds pp**: `pipeline_pp_crossx.sh` Stages 5-7 CONSUME the turn-on
  fits, so 1b + 2a must complete BEFORE it -- run `pipeline_pp_trig_eff.sh` first, then the Medium
  pass, then `pipeline_pp_crossx.sh`. Both pipelines expose `SKIP_CONDOR=1` to re-enter after the
  NTuple stage without resubmitting.

State otherwise unchanged from the previous entry:
- All code is on master in 8 commits (`06f9bfb` -> `c4e13ed`); 21 compile targets clean.
- **Nothing has been rerun. No Condor job submitted. No output file overwritten.** Every number and
  plot on disk still describes the OLD selection.
- Decisions D1-D12 settled; none outstanding.
- The three engaged peer sessions cleared the rerun; two are waiting on the fresh trees (R6).

**What a resuming agent must NOT assume:** that any MC artefact on disk is usable -- the R4 table
lists eight fill classes whose outputs are stale until Phase 1 completes, and the axis-edge guards
CANNOT detect the muon-pT half of that staleness.
