# ΔR-Correction Cell-Binning Approaches — 4-way comparison + MC closure on the crossx binning

**Mode:** Implementation. **Created:** 2026-08-24. **Session:** "MC trigger efficiency 4 approaches & closure update".
**Parents:** `mc_trigger_efficiency.md` (ACTIVE, Step 3 / R27 / R32) and `mc_trig_eff_closure.md`
(ACTIVE, the closure machinery this doc re-binns and extends). Nothing here re-derives Step 3's
measurement; it re-*groups* the cells it is fitted in, and re-bins the closure's presentation.

---

## Objective

The Step-3 ΔR correction ε_ΔR is fitted per (pair p_T, pair η) cell, and the top pair-p_T cells are
statistics-starved: their fits are noise-dominated, some are rejected outright, and the closure
collapses there (`mc_trig_eff_closure.md` R5). One remedy is already in the repo — merge the last
two pair-p_T bins (`nocorr_ptmerge`, parent R32). This doc adds **two more ways of buying the same
statistics** and compares all four on an equal footing:

| # | approach | ε_ΔR fit cells | mode token | plot subdirectory |
|---|---|---|---|---|
| **A** | no combination (the un-merged reference) | 8 pair-p_T × 9 pair-η = 72 | `nocorr` | `no_plateau_correction/` |
| **B** | last two pair-p_T bins merged | 7 × 9 = 63 | `nocorr_ptmerge` | `no_plateau_correction_last2ptbins_merged/` |
| **C** | pair-η merged into 3 | 8 × 3 = 24 | `nocorr_etamerge` | `no_plateau_correction_paireta_merged/` |
| **D** | pair-η merged into 3 AND last two pair-p_T merged | 7 × 3 = 21 | `nocorr_etamerge_ptmerge` | `no_plateau_correction_paireta_merged_last2ptbins_merged/` |

The three merged pair-η cells are **negative-η endcap (−2.4, −1.0)**, **barrel (−1.0, +1.0)**,
**positive-η endcap (+1.0, +2.4)** — the physical detector regions, not an arbitrary regrouping.

Each approach gets its own Step-3 fit tree (`expo` + `polyu_fixedRp` + `interp`, no-plateau-correction,
opposite sign) and its own MC closure. The four closures are then **overlaid in one figure** so the
approach that corrects the pp24 cross-section most accurately can be chosen from evidence.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. Step-3 no-plateau-correction fit trees exist for **all four** approaches, opposite sign, for
     `expo`, `polyu_fixedRp` and `interp`, at both working points, with their fit-plot sets in the
     four subdirectories named in the table above.
  2. **MC closure re-run for all four approaches**, on the **pp24 crossx binning** (§PP-1):
     `ParamsSet::pT_bins_150` (15 log bins, 8–150 GeV) × the 9
     `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` panels — the SAME binning in every
     approach, whatever its correction cells are. One subdirectory per approach, both sample
     versions (all opposite-sign / single-b signal region), both WPs.
  3. The corrected series of every closure is the **expo → polyu_fixedRp → interp** cascade (§PP-3).
  4. Approach A additionally keeps the per-fit-form figures (separate `expo` and `polyu` lines) in
     its own **separate subdirectory**.
  5. The **5-line comparison figure** (no-trigger + the 4 approaches' corrected series, with the
     closure ratio pad), in its own new subdirectory, both sample versions, colours kRed / kBlue /
     kGreen+2 / kMagenta for A/B/C/D. Its per-approach values MUST equal the ones in that
     approach's own closure subdirectory.
  6. `/review-analysis-code` on the C++/RDF, `/review-plot` on the figures; this doc + INDEX
     updated; commits.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

ε_ΔR^2mu4(ΔR) is the probability that a *close-by* pair fires 2mu4 relative to two independent
legs. It is measured by inverse weighting (parent §3.3) and fitted **per (pair p_T, pair η) cell**,
because the correlation is not the same in the barrel as in the endcaps or at low vs high pair p_T.
Finer cells therefore mean a more faithful correction — until the cells run out of pairs, at which
point the fit stops measuring the shape and starts fitting noise, and the *delivered* correction
becomes worse than a coarser one would have been. Parent R32 measured this trade-off for the pair-p_T
axis: merging the top two bins removes every dead cell and every ε_ΔR > 1.3 artefact, at the cost of
describing a mixture of two cells.

**Pair η is the second axis with the same trade-off, and it has never been tested.** The 9-bin
pair-η grid is inherited from the cross-section's presentation binning, not chosen for ε_ΔR's
statistics. Grouping it into the three physical detector regions — negative endcap, barrel, positive
endcap — is the pair-η analogue of R32's pair-p_T merge: it triples the pairs per fit inside each
region while keeping the one distinction that is physically motivated (the endcap L1 trigger
geometry differs from the barrel's, and the forward-negative r16578 anomaly, parent R8/R10/R14,
lives in one endcap only, which is why the two endcaps are NOT merged with each other).

Which of the four is right is an empirical question about the delivered correction, and the MC
closure is the instrument that answers it.

### 2. Top-level equation

Unchanged from `mc_trig_eff_closure.md` §2 — the per-pair pp24 weight is

```
ε_trig^pair(pair) = ε(p_T1, q·η1) · ε(p_T2, q·η2) · ε_ΔR^2mu4(ΔR ; cell)
```

and the closure variable is

```
C(p_T^pair, η^pair) = Σ_{pairs passing 2mu4} w_MC / ε_trig^pair
                      ─────────────────────────────────────────
                          Σ_{all pairs} w_MC
```

with the applied single-muon ε = **ε_MC** (closure doc D4: the test is self-contained, so its
residual is the ΔR correction alone). The four approaches differ ONLY in **which cell** supplies
ε_ΔR for a given pair — i.e. in how the (pair p_T, pair η) plane is partitioned before the fit.

### 3. Step-by-step method

#### §PP-1 The closure is binned on the pp24 CROSS-SECTION's binning — always, in every approach

The closure histograms and figures use

- pair p_T: **`ParamsSet::pT_bins_150`** — 15 logarithmic bins, 8 → 150 GeV (the crossx "8–150 GeV
  version", `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts` / `RDFBasedHistFillingPP.cxx`);
- pair η: the **9** `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` panels — the same
  ranges `SingleBCrossxPlotterBase::DrawPairPtByEta` slices the cross-section in.

**This is INDEPENDENT of the ε_ΔR cell grid**, and deliberately so (user, 2026-08-24). The question
the comparison answers is *"which approach gives the pp24 cross-section the most accurate efficiency
correction?"*, so all four must be read out in the cross-section's own cells; a closure drawn on each
approach's own correction cells would put four different x-axes side by side and make the comparison
meaningless. It supersedes closure-doc **D2** (see D8 below).

The correction itself is still looked up **per pair**, in its own approach's cell grid — a pair's
ε_ΔR does not know which histogram bin it will later be filled into, so nothing about the physics
changes when the presentation binning does.

#### §PP-2 The pair-η grouping

Group boundaries: the **interior** edges −1.0 and +1.0, which MUST already exist in
`pair_eta_proj_ranges_coarse_incl_gap`; the outer edges are read from the source axis. A boundary
that is not an existing edge is a THROW, never a rebin. Consequences, exactly as for the pair-p_T
merge (parent R32):

- **It is NOT a new binning** (.claude/CLAUDE.md §Binnings). The FILLED histograms keep the 9
  pair-η bins; a merged cell is the three source bins **projected together** — num / denom / errA /
  errB summed *before* the ratio — which is numerically identical to having filled a 3-bin axis.
- The grouping is derived in ONE place from the source histogram's own axis
  (`dr_correction_cell_groups.h`); no edge value is retyped in a macro, a report or a plot.
- The variant is opt-in and suffixed (`_nocorr_etamerge`), so it can neither overwrite nor be
  mistaken for the un-merged trees.
- A merged cell describes a **mixture** of its three source cells. Inside an endcap group the three
  source bins span 1.4/0.5/0.4 units of η and the forward-negative anomaly sits in only one of them,
  so the mixture is not uniform — this is a cost of the approach, to be read off the closure, not
  assumed away.

#### §PP-3 The delivered ε_ΔR: the expo → polyu → interp cascade

Per (pair p_T, pair η) fit cell, in order:

1. the **exponential** fit `f = C + A·exp[−(ΔR/λ)^p]`, if accepted;
2. else the **polynomial** fit `f = C + u²(a₂ + a₃u + a₄u²)`, if accepted;
3. else the **interpolation** — linear through the measured points below R_p, flat at the last
   measured knot above it — if accepted;
4. else no correction (ε_ΔR ≡ 1), counted and printed.

"Accepted" = the producer's own `h_step3_fit_ok == 1` **and** the baseline screen
`DrCorrPlateauUsable(C, σ_C)`, because C is what the cell is divided by. In every tier the delivered
correction is

```
ε_ΔR(ΔR) = f(ΔR) / C     for ΔR < 1 ;      ε_ΔR(ΔR) = 1     for ΔR ≥ 1
```

with C the tier's own baseline: the fitted free asymptote for tiers 1–2, and the **flat-branch
value** (the last measured knot, `fit_dr_corrections.cxx`'s `flat_val`) for the interpolation. The
interpolation therefore enters the cascade on exactly the same footing as the two fits, not as a
different kind of object.

**Tier 3 is the `interp` fit, NOT the raw measured bins** (user, 2026-08-24). It replaces the
raw-bin placeholder that `dr_correction_apply.h` has carried since the closure was first built —
that placeholder is flagged TEMPORARY there, and per-cell interpolation is one of the three
resolutions parent-doc R26 names. The raw-bin branch survives only as tier 4's implementation
detail for a cell where even the interpolation fails its positivity screen; it is counted and
printed whenever it is reached.

#### §PP-4 What is plotted

Per approach, per WP, per sample version — the two versions are unchanged from the closure doc
(**all opposite-sign pairs**, and **data-like single-b signal pairs**):

- **upper pad** dσ/dp_T^pair (log–log), one panel per pair-η bin (9 → 3×3): the no-trigger-requirement
  series, and the cascade-corrected 2mu4 series.
- **lower pad** the closure ratio C, with a line at 1, conditional (binomial-correct) errors
  (`SetConditionalRatioErrors`, closure doc §3.4).

Approach A additionally keeps its **per-fit-form** figures — separate `expo` and `polyu` lines —
in a subdirectory of its own, because they answer a different question (how much the two parametric
forms disagree, parent R26) than the cascade figure does.

**The comparison figure**: the same layout, with **five** series — no trigger requirement (black),
and the cascade-corrected series of A (kRed), B (kBlue), C (kGreen+2), D (kMagenta). Each
approach's points MUST be numerically identical to the ones in its own subdirectory: same fill,
same cascade, only overlaid.

### 4. Negative constraints

- **The pp24 cross-section application is NOT changed by this work** (user, 2026-08-24).
  `DrCorrCrossxMethod/Sign/Mode` and `DrCorrectionCrossxEvaluator` stay exactly as they are —
  including its `expo → polyu → raw bins` cascade — until an approach is chosen. Consequence,
  stated rather than hidden: approach B's closure no longer tests the crossx evaluator
  byte-for-byte (closure doc D5 is relaxed; see D9).
- **Do NOT change any binning.** The pair-p_T correction cells come from `MCTrigEffPairPt::Edges`,
  the pair-η ones from `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`, the closure's
  presentation axis from `ParamsSet::pT_bins_150` — all read, never retyped. The merges are
  projections of those axes, not new edges.
- **The two endcaps are never merged with each other.** They differ physically (parent R8/R10/R14:
  the r16578 anomaly is forward-*negative* only), and the closure's principal finding (closure doc
  R5) lives in one of them.
- **No trigger requirement anywhere in the denominator**, and **do NOT apply ε_ΔR^single (Step 4)**
  to the pp weight — both inherited unchanged from the closure doc §4.
- **Do NOT re-derive the sample selection from raw NTUPs.** Everything reads the ntuple-processing
  output through `MCTrigEffPairSel::Step3PairSelection()`.

---

## Scope

**In:** pp24 fullsim FULL sample (`pp_full`), 2mu4, opposite sign, Step 3 only, the
no-plateau-correction mode family; both WPs; the fit stage, its plot stage, the closure fill, the
closure plot, the new comparison plot, and the two drivers.

**Out:** Pb+Pb / the HIJING overlay; same-sign closure; the `corr` (plateau-corrected) mode for the
new approaches; changing what the pp24 cross-section applies; resolving parent R26.

## Design Decisions

### D8: the closure's presentation binning is the CROSS-SECTION's, not the correction's (user, 2026-08-24)
**Old (closure doc D2):** the closure's pair-p_T axis IS the correction's own cell axis.
**New:** it is `ParamsSet::pT_bins_150` × the 9 crossx pair-η panels, in every approach.
**Reason:** the deliverable of this thread is a *choice between four corrections*, judged by how
accurately each corrects the pp24 cross-section. Four figures on four different x-axes cannot be
compared, and none of them would be on the axis the decision is about. D2's concern — never mixing
two binnings of one quantity — is still respected: the closure figure now consistently uses ONE
binning (the cross-section's) for the plotted quantity, and the correction cells are not plotted at
all. What D2 forbade was *deriving* a 1-D view by regrouping a different axis of the same object;
here the correction is applied per pair and the histogram is filled in the crossx cells directly.

### D9: closure doc D5 is relaxed — the closure uses a general cascade evaluator, not the crossx one (2026-08-24)
`DrCorrectionCrossxEvaluator` hard-codes the crossx defaults and the `expo → polyu → raw` tier list.
The four approaches need `expo → polyu → interp` on four different modes, so the closure now uses a
new `Utilities/DrCorrectionCascadeEvaluator.h`. `DrCorrectionCrossxEvaluator` is left untouched (the
user's decision to leave the cross-section alone), which means two cascade implementations coexist.
**Mitigation:** a MIRROR NOTICE in both files naming the other, and Remaining Work 1 — when an
approach is chosen, `DrCorrectionCrossxEvaluator` is to be *replaced* by the general class
configured with the winning mode and tier list, not edited to match it.

### D10: the pair-η group boundaries are validated against the source axis, never invented (2026-08-24)
`MakeDrEtaGroups` takes the interior boundaries {−1.0, +1.0} and *looks them up* in the histogram's
own pair-η axis, throwing if either is not an existing edge. The outer edges come from the axis.
This is the same discipline `MakeDrPtGroups` uses (which reads every edge from the axis and merges
by index) applied to a grouping that cannot be expressed by an index rule alone.

---

## Implementation Plan

1. [x] `dr_correction_pt_groups.h` → **`dr_correction_cell_groups.h`**: generalise `DrPtGroups` to
   `DrAxisGroups`, add `MakeDrEtaGroups` (per §PP-2 / D10), and give `DrGroupCellRatio` /
   `BookDrGroupMap` BOTH axes' groupings. → `/review-analysis-code`.
2. [x] `dr_correction_sample_cfg.h`: the two new mode tokens, their directories, and
   `DrCorrModeMergeEta()`. → same review.
3. [x] `fit_dr_corrections.cxx`: thread the pair-η grouping through the cell loop, the plateau
   re-measurement (merged cells have no entry in the on-disk map) and the reports (per §PP-2).
4. [x] `plot_dr_correction_fits.cxx`: same grouping; 1×3 panel grid for 3 pair-η cells.
5. [x] `dr_correction_apply.h`: pair-η grouping in the raw fallback; **`interp` support** — the
   `gknots` TGraph tier with C = the flat-branch value (per §PP-3).
6. [x] NEW `Utilities/DrCorrectionCascadeEvaluator.h`: the expo → polyu → interp → none cascade for
   an arbitrary mode (per §PP-3 / D9).
7. [x] `FillMCTrigEffClosure.cxx`: the crossx presentation binning (§PP-1), four modes, the cascade
   series everywhere + the two per-form series for approach A.
8. [x] `plot_mc_trig_eff_closure.cxx`: mode subdirectories; approach A's per-form figures into
   their own subdirectory. → `/review-plot`.
9. [x] NEW `plot_mc_trig_eff_closure_compare.cxx`: the 5-line figure (§PP-4). → `/review-plot`.
10. [x] Drivers: `run_dr_correction_fits.sh` (new modes) and `run_mc_trigeff_closure.sh` (four
    modes, `interp` in METHODS, the comparison stage).
11. [ ] Run everything, sanity-check, record numbers here; reviews; commit.

## Progress Log

*(append-only; newest entries at the END)*

- 2026-08-24 — **Steps 1, 2, 5, 6, 7, 8, 9, 10 written and COMPILING** (steps 3, 4 and the fit
  driver were delegated to a subagent working on `fit_dr_corrections.cxx`,
  `plot_dr_correction_fits.cxx` and `run_dr_correction_fits.sh` only — disjoint file ownership;
  the orchestrator owns git and every other file). What was written:
  * `plotting_codes/trig_effcy/mc_based/dr_correction_pt_groups.h` → **`dr_correction_cell_groups.h`**
    (`git mv`). `DrPtGroups` generalised to `DrAxisGroups`; new `MakeDrIdentityGroups`,
    `MakeDrEtaGroups` (interior boundaries {−1.0,+1.0} **looked up** in the source axis, THROW if
    either is not an existing edge — D10), `DrGroupsDescribe`; `DrGroupCellRatio` and
    `BookDrGroupMap` now take BOTH axes' groupings. Verified on the canonical axes: 3 η groups =
    source bins 1–3 / 4–6 / 7–9, edges −2.4 / −1.0 / +1.0 / +2.4; 7 pair-pT groups
    8 / 11.54 / 16.65 / 24.01 / 34.64 / 49.97 / 72.08 / 150.
  * `dr_correction_sample_cfg.h`: tokens `nocorr_etamerge` / `nocorr_etamerge_ptmerge`, their
    directories, `DrCorrModeMergeEta()`; `DrCorrModeNoPlateau` now answers for every `nocorr*`
    mode and `DrCorrModeMergeLastTwoPt` for both pT-merging ones.
  * `dr_correction_apply.h`: **`interp` support** — the `gknots_step3_pt<i>_eta<j>` TGraph tier
    with C = the flat-branch value `gk->Eval(1.0)`, screened by `DrCorrPlateauUsable(C, 0)`; the
    `method == "interp"` THROW is gone; `Cell::HasCurve()` added; the raw-bin fallback now projects
    BOTH groupings.
  * NEW `Utilities/DrCorrectionCascadeEvaluator.h` — the N-tier cascade for any mode/sign, with a
    canonical-binning guard that groups `ParamsSet::pair_pt_coarse_bins` /
    `pair_eta_proj_ranges_coarse_incl_gap` through the SAME helpers the fit stage uses. Carries the
    MIRROR NOTICE for `DrCorrectionCrossxEvaluator` (D9).
  * `FillMCTrigEffClosure.cxx`: crossx presentation axes (D8) decoupled from the correction cells,
    with a guard that the two cover the same region; four modes; the cascade series in every mode
    (key `cascade`) plus the two per-form series in `nocorr` only; provenance stamp now records
    both grids.
  * `plot_mc_trig_eff_closure.cxx`: restructured into `DrawClosureSet`, called once per figure set;
    `nocorr` gets a second set in `separate_fit_forms/`. The canvas headline reads the cell counts
    out of the file's own provenance rather than re-deriving them.
  * NEW `plot_mc_trig_eff_closure_compare.cxx`: the 5-line overlay; the four no-trigger
    denominators are compared bin by bin and it THROWS if they differ (that identity is what makes
    the overlay legitimate).
  * `pipelines/run_mc_trigeff_closure.sh`: four modes, `interp` in METHODS, per-form PNG
    validation, Stage 4 (the overlay), 4-bin upstream regen looped over the modes it supports.
  All five compiled objects rebuilt with ACLiC and verified newer than their sources.

- 2026-08-24 — Doc created. Triage: read `INDEX.md`, `mc_trig_eff_closure.md` (full),
  `mc_trigger_efficiency.md` §Physics Procedure + R27 + R32; skipped the other ACTIVE docs
  (`intnote_trig_reco_eff_sections.md`, `pp24_crossx_rerun_2026_08.md`, `muon_gap_cuts_acceptance.md`,
  `pp_trig_eff_highpt_jump.md`, `raa_from_rdf_crossx.md`, the template-fit family,
  `analysis_roadmap_2026_06.md`, `analysis_status_summary.md`, `academic_writing_workflow.md`,
  `pythia_fullsim_pp24_full_sample_skim.md`, `tight_wp_default_change.md`) — out of scope.
  Two blocking ambiguities resolved with the user BEFORE any code: (i) the third cascade tier is
  the `interp` fit, replacing the raw-bin placeholder (§PP-3); (ii) the pp24 cross-section
  application is NOT changed by this work (§4 / D9). Plan recorded before any edit.

## Results & Observations

### R1. Steps 3+4 — the pair-η grouping in the fit and its plot stage (2026-08-24, delegated)

*(merged from the subagent's scratch doc `_sub_etamerge_fitstage_1.md`, now deleted; every number
below was measured, not asserted.)*

**Regression on the un-merged path — clean.** `nocorr` / `expo` / pp_full / Tight / step 3 /
opposite sign: `plateau_guard_report_opposite_sign.txt` **byte-identical** to the pre-work copy;
`fit_report_opposite_sign.txt` identical in **every one of its 72 cell rows**, every plateau, every
fitted parameter, every χ²/ndf, the unusable count and the at-limit list. It differs in **one line
only**, and that line is a **stale-baseline artefact rather than a regression**: the "(the only
statistically meaningful curve for a 10k-event TEST sample)" parenthetical on the inclusive χ²/ndf
became sample-conditional in commit `22b8454` (2026-08-18 10:10), while the baseline on disk was
written 2026-08-17 23:13 by the previous binary. pp24 is a FULL production, so omitting it is the
correct current output. The plot path is unchanged too — `readback_check.txt` byte-identical,
9 PNGs at 1556×1540 px.

**The new modes, measured** (pp_full, Tight, step 3, `expo`, opposite sign):

| mode | `h_step3_fit_ok` | cells | accepted (`fit_ok = 1`) | unusable | χ²/ndf over the cells |
|---|---|---|---|---|---|
| `nocorr` (reference) | 8 × 9 | 72 | 68 | 4 | mean 1.809, median 1.511 (n = 70) |
| `nocorr_etamerge` | **8 × 3** | **24** | **23** | 1 | mean 1.852, median 1.421 (n = 24) |
| `nocorr_etamerge_ptmerge` | **7 × 3** | **21** | **20** | 1 | mean 2.007, median 1.529 (n = 21) |

The 68/72 of the reference reproduces `mc_trigger_efficiency.md` R27's independently recorded
"Tight `expo`, no correction: 68" — an outside check that the un-merged path did not move.
**χ²/ndf does not improve with the merge, and should not be expected to** — the same point R32
makes for the pair-p_T merge: with three times the pairs the error bars shrink, so the same shape
mismatch costs more χ². What the merge buys is coverage, and that is what the closure measures.

The merged pair-η axis reads **[−2.4, −1.0)**, **[−1.0, +1.0)**, **[+1.0, +2.4)** in the fit file
and in the report labels — the three physical detector regions, from the source axis's own edges.
`polyu_fixedRp` (24 cells, mean 1.913) and `interp` run clean on the new modes as well.

**Panel grid.** 3 pair-η cells are drawn **1 × 3** (1556 × 600 px), not the 2 × 2 the bare
`ceil(sqrt(N))` rule would give; 9 cells stay 3 × 3 (1556 × 1540 px). The header strip needed no
change and the reason is now in the code: the canvas is sized `520·ncol × (470·nrow + header)`, so
a panel is 520 × 470 px in every layout and the strip, specified in canvas pixels, keeps its
absolute geometry (checked numerically and by opening the PNG).

**One process note worth keeping.** An intermediate edit left an unbalanced `)`; ACLiC reported
`expected ':'` **and left the previous `.so` in place**. It was caught only by the
check-the-log-AND-the-timestamp rule the drivers already encode, and every measurement above was
redone with the final binary. This is the "stale .so" trap in its natural habitat.

## Remaining Work

1. **Replace `DrCorrectionCrossxEvaluator` with the general cascade class** once an approach is
   chosen (D9) — two cascade implementations must not coexist longer than this comparison.

## Latest Stage

**2026-08-24 — closure side written and compiling; waiting on the fit stage.** Steps 1, 2, 5–10 are
done (see the Progress Log entry). Steps 3 and 4 — the pair-η grouping inside
`fit_dr_corrections.cxx` and `plot_dr_correction_fits.cxx`, plus `run_dr_correction_fits.sh` — are
with a delegated subagent (scratch doc `_sub_etamerge_fitstage_1.md`), which is also running the
byte-identical regression of the existing `nocorr` fit report.

**Next, in order:** (1) merge the subagent's work and re-verify the regression myself; (2) run
`run_dr_correction_fits.sh` for the two NEW modes, opposite sign, all three methods, both WPs;
(3) run `run_mc_trigeff_closure.sh` for all four approaches; (4) `/review-analysis-code` +
`/review-plot`; (5) record numbers here and commit.
