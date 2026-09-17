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

**SUPERSEDED 2026-09-03 (D11):** the pair-η merge is now the SIGN-INDEPENDENT |η^pair| grouping
**|η| < 1.0**, **1.0 ≤ |η| < 2.0**, **2.0 ≤ |η| < 2.4** — see D11 below. Every "negative-eta
endcap / positive-eta endcap" reference in this doc from this point back describes the RETIRED
2026-08-24 grouping; R1/R2's numbers were measured against it and are historical, not current.

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

#### §PP-2 The pair-η grouping (SUPERSEDED 2026-09-03 — see D11 for the current definition)

**Original (2026-08-24), retired:** group boundaries the **interior** edges −1.0 and +1.0
(negative-eta endcap / barrel / positive-eta endcap), which MUST already exist in
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

**Current (2026-09-03), per D11:** group boundaries the **interior |η| edges** 1.0 and 2.0 (a
sign-independent fold: |η| < 1.0, 1.0 ≤ |η| < 2.0, 2.0 ≤ |η| < 2.4), each of which MUST be an
existing edge on **both sides of 0** of `pair_eta_proj_ranges_coarse_incl_gap` (the source axis
must be symmetric about 0 for a sign-independent fold to be well defined). The barrel group
(|η| < 1.0) straddles zero and is therefore one contiguous source range; each forward group is the
**union** of a negative-eta and a positive-eta source sub-range — the one grouping on this axis a
single contiguous range cannot describe, so `DrAxisGroups` (`dr_correction_cell_groups.h`) now
carries a **list** of source sub-ranges per group rather than one, and
`DrCellRatioMultiRange` (`dr_correction_ratio.h`) sums every sub-range's num/denom/errA/errB before
the ratio, generalizing `DrCellRatioRange`. Same non-new-binning guarantees as above; same opt-in
suffixing. **`DrCorrectionEvaluator::Eval` (`dr_correction_apply.h`) looks a pair's cell up by
`|pair_eta|` whenever the loaded mode folds eta** (`eta_folded` member, set from
`DrCorrModeMergeEta`) — the folded output axis runs 0 → η_max, so a signed lookup would send every
negative-η pair to the underflow bin and silently return `eps_dR = 1` for half the sample; this was
fixed in the same edit that introduced the fold, even though nothing downstream calls `Eval` with
a folded mode yet (§4 below).

#### §PP-3 The delivered ε_ΔR: the expo → polyu → interp cascade

Per (pair p_T, pair η) fit cell, in order:

1. the **exponential** fit `f = C + A·exp[−(ΔR/λ)^p]`, if accepted;
2. else the **polynomial** fit `f = C + u²[A + a₃(u − 1) + a₄(u² − 1)]`, if accepted — the
   same quartic as `C + a₂u² + a₃u³ + a₄u⁴` with `a₂ = A − a₃ − a₄`, written in
   `A ≡ f(0) − C` since 2026-09-08 so that the user's Step-3 requirement "at ΔR = 0 the
   efficiency cannot exceed the plateau" is the single fit limit `A ≤ 0`
   (`mc_trigger_efficiency.md` R34);
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

#### §PP-5 The pair-η-summed non-closure chi^2/ndof (user, 2026-08-25; RENORMALIZED 2026-09-03)

A second figure in the same `approach_comparison/` subdirectory compresses the nine pair-η panels
of one approach into a single curve, so the four can be read against each other bin by bin:

```
chi^2(p_T^pair)      = Sum over pair-eta bins of ( C(p_T^pair, eta^pair) - 1 )^2 / sigma_C^2
chi^2/ndof(p_T^pair) = chi^2(p_T^pair) / n_dof(p_T^pair)
```

over exactly the C (and its conditional/binomial-correct sigma_C) the 9-panel figure's lower pads
draw, on the same crossx presentation binning (§PP-1). One curve per approach, four overlaid in one
panel with a TLegend, plus a dashed reference line at chi^2/ndof = 1 (perfect closure).

**RENORMALIZED (user, 2026-09-03): it is now a PROPER chi^2** — the 2026-08-25 original divided
neither term by its uncertainty ("a distance from closure, NOT a chi^2 ... the user's
definition"), and R2 below showed why that reading is not meaningful: the un-normalized sum was
dominated by whichever cells happened to be statistics-starved, on equal footing with cells that
are genuinely mis-corrected (the reviewer's own check on it, "bin 14 alone is 6.998 ± 10.8 — a
~0.6σ effect," makes the point directly). Each term is now divided by `sigma_C^2` before summing,
so `chi^2/ndof` says how SIGNIFICANT the deviation from closure is, comparably across pair-p_T
bins with different `n_dof`, not merely how far the raw ratios sit from 1. The bar on each point is
NOT a propagated point error — a chi^2 statistic is judged against its own null distribution, not
against itself with a Gaussian uncertainty — it is the null-hypothesis (perfect-closure) width
`sqrt(2/n_dof(p_T))`, i.e. how much `chi^2/ndof` would fluctuate around 1 by statistics alone if
the correction closed exactly.

Two constraints on how it is formed (unchanged from 2026-08-25):

- **A pair-η cell enters the sum only where the no-trigger DENOMINATOR is non-empty AND
  sigma_C > 0.** An empty cell has C = 0 and sigma_C = 0 by construction (`dr_correction_ratio.h`
  returns both as 0 for D <= 0) and is excluded, avoiding a division by zero. A cell with a filled
  denominator but nothing passing the trigger has C = 0 and sigma_C = max(0,1)/n_eff > 0 (the k=0
  boundary case in `SetConditionalRatioErrors`), so it correctly contributes a large but finite
  term — real total non-closure, not a spurious infinity.
- **The four sums run over the SAME cells**, because the four denominators are already checked
  bin by bin to be identical; the same cells are skipped in every approach, so `n_dof(p_T)` does
  not depend on the approach.

`chi^2` (and its `n_dof`) are pure sums, **not a density** — never width-scaled, and their
`Integral` is taken without the `"width"` option. Do NOT combine pair-η bins' statistics (numerator
and denominator) before forming C and only then compare that combined ratio to 1: pair-η bins can
carry opposite-sign non-closure (the r16578 forward anomaly, parent R8/R10/R14, is not charge/η
blind), and pre-combining would let such structured, physically real non-closure cancel in the
combined ratio — exactly the failure §PP-1's "read out in the cross-section's own cells" and D2's
never-regroup-to-hide-structure principle exist to prevent. Summing independent per-η
`(C-1)^2/sigma_C^2` terms is the standard way to build a multi-bin chi^2 across independent cells
and preserves that structure; only the missing `/sigma_C^2` normalization was the defect.

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
- **SUPERSEDED 2026-09-03 (D11):** this doc's original constraint was "the two endcaps are never
  merged with each other" (parent R8/R10/R14: the r16578 anomaly is forward-*negative* only). D11
  reverses this — the pair-η merge now explicitly DOES fold the negative- and positive-eta endcap
  bins together by |η|, because the measured dR correction barely depends on sign. The r16578
  anomaly is a single-muon-leg reconstruction effect (parent R8/R10/R14), a different quantity from
  the pair-level dR correlation this doc fits; the two findings are not in tension.
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
**Superseded 2026-09-03 by D11's boundary set — the LOOKUP discipline itself is unchanged.**

### D11: the pair-η merge is now SIGN-INDEPENDENT (|η| fold), not the signed 3-region split (user, 2026-09-03)
**Old (D10 / 2026-08-24):** negative-eta endcap (−2.4,−1.0) / barrel (−1.0,1.0) / positive-eta
endcap (1.0,2.4) — a signed 3-way split, chosen because the r16578 forward anomaly is
negative-eta-only (parent R8/R10/R14) and the two endcaps were deliberately kept apart.
**New:** |η^pair| < 1.0 (barrel) / 1.0 ≤ |η| < 2.0 / 2.0 ≤ |η| < 2.4 — a sign-independent fold.
**Reason (user, 2026-09-03):** having now measured the dR correction in both groupings' cells, the
correlation barely depends on the SIGN of pair η, while the >2.0-vs-<2.0 split WITHIN the endcap is
a much bigger effect than any negative/positive asymmetry — the boundary that matters is how
forward a pair is, not which side of η = 0 it is on. This is the **new default** merged pair-η
grouping for the dR correction; the un-merged 9-bin grid (`nocorr`, `nocorr_ptmerge`) is untouched
and remains the alternative.
**Consequences:**
- `MakeDrEtaGroups`'s boundary set moves from `DrEtaMergeInteriorBoundaries()` = {−1.0, +1.0}
  (signed, contiguous groups) to `DrEtaAbsMergeInteriorBoundaries()` = {1.0, 2.0} (|η|, folded
  groups) — see the current §PP-2 above and `dr_correction_cell_groups.h`.
- `DrAxisGroups` generalises from one `[lo,hi]` per group to a **list** of source sub-ranges per
  group (`DrCellRatioMultiRange` in `dr_correction_ratio.h` sums every sub-range before the
  ratio) — the only structural change this decision required, because a forward |η| group is the
  union of a negative- and a positive-eta source sub-range.
- **R1/R2's numbers (below) describe the RETIRED signed grouping and are historical.** They are
  kept, not deleted (tracking-doc discipline), but must not be read as current.
- `DrCorrectionEvaluator::Eval` was fixed to look a pair's cell up by `|pair_eta|` for a folded
  mode (previously it always used the signed value, which would have silently broken any consumer
  of a merge-eta mode — caught and fixed in this edit even though nothing calls `Eval` with a
  merge-eta mode yet, §4).
- **Overwrites, does not preserve, the old grouping's outputs**: `nocorr_etamerge` and
  `nocorr_etamerge_ptmerge` fit files and plots (user: "overwrite previous results" — this is not a
  third alternative alongside the signed split, it replaces it).
- **MC closure re-run DONE (2026-09-03, same day, following request):** `FillMCTrigEffClosure.cxx`'s
  coverage guard (`.edges.front()/.back()` against the signed presentation axis) DID throw for a
  merge-eta mode exactly as flagged, since `MakeDrEtaGroups(...).edges` is |η|-based [0, 2.4] when
  folded, not signed [−2.4, 2.4]. Fixed by comparing in |η| space for a folded mode (presentation
  |η| extent is [0, `eta_edges.back()`], since the presentation axis is symmetric about 0 — checked
  in-code, not assumed). **A second, more serious bug was found while fixing the first**: the
  RDF sample Filter two lines below the guard applied the SAME (folded, |η|-space) bounds directly
  to the SIGNED `pair_eta` branch — `pair_eta >= 0.0 && pair_eta < 2.4` — which would have silently
  DROPPED EVERY NEGATIVE-η pair from the merge-eta closure sample (a ~50% sample loss, no crash,
  no warning: exactly the class of silent binning bug `.claude/CLAUDE.md` "Binnings" warns about).
  Fixed by mapping the folded cell extent back to signed `pair_eta` bounds
  `[-eta_hi, +eta_hi]` before filtering; verified post-fix that the merge-eta selected-pair survival
  fraction (66.63%) now matches the un-merged reference exactly, confirming no pairs are being lost
  to sign. Re-ran `FillMCTrigEffClosure` for `nocorr_etamerge` / `nocorr_etamerge_ptmerge`, both
  WPs (`nocorr` / `nocorr_ptmerge` untouched — their grouping is unaffected by D11) and
  `plot_mc_trig_eff_closure_compare.cxx`, both WPs. **Result is a genuine, sizeable surprise: the
  new fold grouping closes MUCH WORSE than the retired signed grouping** — see R4. The pp24 crossx
  application (`DrCorrCrossxMode()` is `nocorr_ptmerge`, which never merges eta) is unaffected
  either way. **Propagation to the crossx / MC-data comparison plots still awaits explicit user
  approval** — R4 argues AGAINST adopting the fold grouping, not for it. This code fix did **NOT**
  go through a formal `/review-analysis-code` pass (self-verified only, per the user's standing
  instruction this session to skip formal review and look directly) — flag for a review pass before
  this closure result is used to make the grouping decision.

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
12. [x] `plot_mc_trig_eff_closure_compare.cxx`: the **summed squared non-closure** panel
    (user, 2026-08-25; §PP-5). → `/review-plot`.
13. [x] `plot_mc_trig_eff_closure_compare.cxx`: **renormalize** the non-closure panel into a proper
    `chi^2/ndof` (user, 2026-09-03; §PP-5). `/review-plot` did NOT complete (reviewer subagent
    stopped by the user; done from the executor's own verification per explicit user request).

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

- 2026-08-25 — **Step 12: the summed squared non-closure panel** (§PP-5), added to
  `plot_mc_trig_eff_closure_compare.cxx` and written into the SAME
  `closure/approach_comparison/` subdirectory as the 9-panel overlay:
  `closure_compare_nonclosure_squared_{all_opposite_sign,single_b_signal_cuts}.png`, both WPs
  (4 PNGs). `pipelines/run_mc_trigeff_closure.sh` Stage 4 now deletes and validates all four.
  Layout: `#Sigma_{#eta^{pair}}` rather than `#sum` (the big-operator glyph is drawn ~3x the text
  size and collides with the axis labels at any offset that still fits); log y taken only when the
  dynamic range exceeds 100 (it does: 632 all-OS, 727 signal). `/review-plot` PASS at iteration 1,
  every reported number independently re-derived from the TH2Ds by the reviewer (MATCH).
  See R2 for the numbers.

- 2026-09-03 — **Step 13: RENORMALIZED the non-closure panel into a proper χ²/ndof** (§PP-5), per
  the advisor's flag that the 2026-08-25 unnormalized D² was not statistically meaningful. Each
  term `(C_η−1)²` is now divided by `σ_C,η²` (the same conditional/binomial-correct error already
  used for the ratio pads); plotted as `χ²/ndof(p_T) = [Σ_η (C−1)²/σ_C²] / n_dof` with `n_dof` the
  existing per-bin contributing-cell count (`ncell[bx]`, identical across all four approaches — the
  denominators are enforced identical). Point error bar is the null-hypothesis width `√(2/ndof)`,
  not a propagated point error. PNGs renamed `closure_compare_nonclosure_squared_*` →
  `closure_compare_nonclosure_chi2_*` (macro `file_chi2` field + `run_mc_trigeff_closure.sh`
  `CMP_PNGS`/comment/log-grep, all updated together). Deliberately reused the EXISTING closure ROOT
  files on disk (2026-08-24/25, verified by mtime to predate D11 — this fix does not touch the
  pair-η-grouping work; C/D's numbers below are still on the retired signed grouping, per D11's own
  note). Compiled clean (ACLiC, no new warnings), rerun both WPs. **Formal `/review-plot` did not
  complete — the reviewer subagent was stopped by the user mid-run, who then explicitly asked to
  mark this done from the executor's own compile/rerun/visual verification instead** (this is a
  deviation from the usual binding review-on-plot-change rule, made explicitly by the user, not a
  skipped step). See R3 for the numbers.

- 2026-09-03 — **D11's deferred closure re-run, executed on request.** Attempting it hit the
  coverage-guard throw D11 had flagged and left unfixed; fixed it to compare in |η| space for a
  folded mode. While fixing it, found and fixed a SECOND, more serious bug in the same function:
  the RDF sample `Filter` applied the folded (|η|-space, [0, 2.4]) cell bounds directly to the
  SIGNED `pair_eta` branch, which would have silently kept only `pair_eta >= 0` — dropping every
  negative-η pair from a merge-eta closure fill (no crash, no warning). Fixed by mapping the folded
  extent back to signed bounds before filtering; verified the merge-eta selected-pair survival
  fraction (66.63%) now matches the un-merged reference exactly. Re-ran `FillMCTrigEffClosure` for
  `nocorr_etamerge`/`nocorr_etamerge_ptmerge`, both WPs, and `plot_mc_trig_eff_closure_compare.cxx`,
  both WPs (`nocorr`/`nocorr_ptmerge` untouched — D11 does not affect their grouping). Compiled
  clean. **Not formally reviewed** (self-verified, per the user's standing instruction this session
  to skip formal review) — flag for `/review-analysis-code` before the result below drives a
  grouping decision. See R4: the new fold grouping's χ²/ndof is 4-8x WORSE than both the retired
  signed grouping (R3) and the un-merged/pT-merged approaches, in every version and WP — a genuine
  reversal of the earlier ranking, not a bug artefact (checked: the fold's own selected-pair count
  and cell coverage match the un-merged reference).


- 2026-09-08 — **★ R4 IS RETRACTED: the fold's poor closure was a BUG in the cascade evaluator, not
  the grouping.** Found by `/review-analysis-code` while reviewing the sibling thread
  `mc_trigeff_single_value_pair_eff.md`. `DrCorrectionCascadeEvaluator::Eval` looked the pair-eta
  cell up with the **signed** `pair_eta` on the **folded** (|eta|, 0 -> 2.2) Y axis, so in every
  `*_etamerge*` mode every pair with `pair_eta < 0` fell into the underflow, was counted "outside
  the cell grid" and received **`eps_dR = 1` — no correction at all**. Half the sample, silently:
  the etamerge closure log read `628453 OUTSIDE the cell grid (50.67%)` against `0 OUTSIDE` for an
  un-folded mode. This is the SAME bug D11 fixed in `DrCorrectionEvaluator::Eval` and in the
  closure's sample Filter — but the CASCADE class had reimplemented those two lines and was missed,
  and the cascade class is the one that builds the delivered series every closure figure draws.
  `DrCorrectionCrossxEvaluator` carried it too (latent: `DrCorrCrossxMode()` is un-folded).
  **Fix:** one shared `DrCorrectionEvaluator::CellLookupEta(pair_eta, folded)` in
  `dr_correction_apply.h`, called from all four sites, so the fold cannot be forgotten a fifth time.
  **Rerun (2026-09-08):** all four modes x both WPs refilled and replotted, plus the 4-approach
  overlay and the chi^2/ndof panel.
  **THE RANKING REVERSES.** chi^2/ndof summed over the 15 crossx pair-p_T bins, SIGNAL version,
  ndof = 133 — both working points, each row read from its OWN log
  (`pipelines/logs_closure/compare_{tight,medium}.log:21`), against R4's same-version rows:

  | | A `nocorr` | B `nocorr_ptmerge` | C `nocorr_etamerge` | D `nocorr_etamerge_ptmerge` |
  |---|---|---|---|---|
  | **Tight**, R4 (bug present) | 8.44 | 8.28 | **61.33** | **61.56** |
  | **Tight**, now | 6.557 (chi2 872.1) | 6.503 (864.9) | **5.029 (668.8)** | 5.129 (682.2) |
  | **Medium**, R4 (bug present) | 7.94 | 8.04 | **64.56** | **64.77** |
  | **Medium**, now | 6.321 (840.7) | 6.266 (833.4) | **5.222 (694.5)** | 5.332 (709.2) |

  **⚠ THE "now" ROWS CARRY TWO CHANGES, NOT ONE.** The fold fix acts ONLY on C and D — A and B are
  un-folded, so `CellLookupEta` is a no-op for them — yet A moves 8.44 -> 6.557 as well. That part
  is the 2026-09-07 gap-cut rerun (`muon_gap_cuts_acceptance.md` F17/F18): every closure file was
  refilled today against Step-3 fits regenerated on the new fiducial windows. **The fold fix's own
  effect is the C/D column**, where 61 -> 5 cannot be attributed to a gap-cut change that moved A
  and B by 20 %. Read the table by column, not by row.
  The |eta| fold is now the BEST of the four by this metric, not the worst. **R4's conclusion
  ("argues AGAINST adopting the fold on closure grounds") no longer holds and must be re-judged by
  the user against the regenerated figures.** The sibling doc's independent measurement agrees:
  above 50 GeV approach D's inclusive closure moved 0.7433 -> 0.9362 and its "sign asymmetry"
  (0.44 at eta^pair in [-1,-0.5) against 0.94 at the mirror) vanished entirely — that pattern was
  the bug's signature. See `mc_trigeff_single_value_pair_eff.md` R5 for the full evidence.

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

### R2. The summed squared non-closure, measured (2026-08-25)

D^2 summed over all 15 crossx pair-p_T bins — **lower is better**, but read the split below before
ranking anything:

| version | A 8x9 | B 7x9 | C 8x3 | D 7x3 |
|---|---|---|---|---|
| Tight, all opposite-sign | 3.666 | 2.175 | 2.086 | **2.066** |
| Tight, single-b signal | 6.028 | 5.955 | 10.499 | **4.693** |
| Medium, all opposite-sign | 3.383 | 2.196 | 2.079 | **2.058** |
| Medium, single-b signal | 5.082 | 5.748 | 10.195 | **4.598** |

**The grand total is not the whole story, and the code now prints the split.** The pair-p_T merge
fuses the top TWO of the 8 eps_dR cells, so below that edge A and B coincide bin by bin and so do
C and D: every difference the p_T merge makes lives in the crossx bins above it. The split bin is
FOUND, not typed (the first bin where a merged-p_T approach departs from its un-merged twin) and
comes out as **bin 12, p_T^pair = 68.65 GeV**:

| version | range | A | B | C | D |
|---|---|---|---|---|---|
| Tight, all OS | bins 1-11 (< 68.65 GeV) | 0.321 | 0.321 | **0.186** | **0.186** |
| Tight, all OS | bins 12-15 | 3.345 | **1.854** | 1.900 | 1.880 |
| Tight, signal | bins 1-11 | **0.473** | **0.473** | 0.714 | 0.714 |
| Tight, signal | bins 12-15 | 5.555 | 5.482 | 9.785 | **3.978** |
| Medium, all OS | bins 1-11 | 0.337 | 0.337 | **0.186** | **0.186** |
| Medium, all OS | bins 12-15 | 3.046 | **1.859** | 1.893 | 1.872 |
| Medium, signal | bins 1-11 | **0.466** | **0.466** | 0.709 | 0.709 |
| Medium, signal | bins 12-15 | 4.616 | 5.283 | 9.486 | **3.890** |

Three readings, all of which the choice of approach has to face:

1. **~90 % of every total sits in the top four bins** (Tight all-OS approach A: 3.345 of 3.666),
   which are exactly the statistics-starved ones. The reviewer's re-derivation: Tight signal,
   approach C, **bin 14 alone is 6.998 +- 10.8** of that column's 10.499 — the single bin that
   makes C look worst is a ~0.6 sigma effect. D^2 also carries a positive noise bias
   E[D^2] = Sum dev_true^2 + Sum sigma_C^2, largest precisely there. A ranking read off the grand
   total is a ranking of those four bins.
2. **Below the merge edge the pair-eta merge flips sign between the two sample versions**: it
   HELPS on all opposite-sign pairs (0.186 vs 0.321) and HURTS in the single-b signal region
   (0.714 vs 0.473), at both WPs. Same correction, different population — the signal region is the
   one the cross-section uses.
3. **Above the edge the p_T merge is what pays**, and it pays in every version (A -> B: 3.345 ->
   1.854 all-OS; the eta merge alone, C, is the worst column in the signal region).

Nothing here is yet a decision — recorded as measurement.

### R3. The non-closure χ²/ndof, renormalized (2026-09-03)

Same inputs as R2 (the retired signed pair-η grouping, historical — D11 postdates this fill), same
four approaches, but each `(C-1)^2` term is now divided by `σ_C^2` before summing, per §PP-5's
2026-09-03 renormalization.

χ²/ndof summed over all 15 crossx pair-p_T bins (`n_dof` shared by all four approaches, since their
denominators are identical):

| version | n_dof | A `nocorr` | B `nocorr_ptmerge` | C `nocorr_etamerge` | D `nocorr_etamerge_ptmerge` |
|---|---|---|---|---|---|
| Tight, all opposite-sign | 135 | 15.474 | 15.333 | 12.129 | 12.185 |
| Tight, single-b signal | 133 | 8.441 | 8.285 | 9.336 | 9.414 |
| Medium, all opposite-sign | 135 | 14.909 | 14.953 | 12.696 | 12.745 |
| Medium, single-b signal | 133 | 7.940 | 8.043 | 9.555 | 9.638 |

Split at the pair-p_T merge edge (bin 12, 68.65 GeV — same bin R2 found):

| version | range | n_dof | A | B | C | D |
|---|---|---|---|---|---|---|
| Tight, all OS | bins 1-11 | 99 | 20.151 | 20.151 | 15.789 | 15.789 |
| Tight, all OS | bins 12-15 | 36 | 2.609 | 2.081 | 2.063 | 2.275 |
| Tight, signal | bins 1-11 | 99 | 10.375 | 10.375 | 11.390 | 11.390 |
| Tight, signal | bins 12-15 | 34 | 2.808 | 2.199 | 3.357 | 3.660 |
| Medium, all OS | bins 1-11 | 99 | 19.618 | 19.618 | 16.569 | 16.569 |
| Medium, all OS | bins 12-15 | 36 | 1.958 | 2.122 | 2.045 | 2.229 |
| Medium, signal | bins 1-11 | 99 | 10.053 | 10.053 | 11.710 | 11.710 |
| Medium, signal | bins 12-15 | 34 | 1.788 | 2.191 | 3.278 | 3.604 |

**This changes the reading, not just the number.** R2's un-normalized D² put ~90% of the total in
the top 4 (statistics-starved) bins and concluded the ranking was mostly a ranking of noise there.
Normalized, it is the OPPOSITE: **the well-measured low-p_T bins (1-11, huge MC statistics) carry
χ²/ndof ~ 10-20 — far from the closure reference of 1 — while the statistics-starved high-p_T bins
(12-15) sit at χ²/ndof ~ 2-4, close to consistent with 1 given their own bar `√(2/ndof)`.** The
residual non-closure is small in absolute size (the "inclusive closure" ratios all sit at
0.988-0.996, i.e. within 0.4-1.2% of 1) but, with pp24's full-sample statistics, that small
residual is highly statistically significant almost everywhere on the axis — not a noise artefact
confined to a few bins. This is a genuinely different conclusion from R2's and belongs in front of
whichever approach gets chosen next (Remaining Work 1 / R2's "nothing here is yet a decision").
The four-approach ranking itself is qualitatively similar to R2 (C/D beat A/B on all-OS; A/B beat
C/D on signal), so which grouping to adopt is still an open call, now on more solid statistical
footing.

### R4. χ²/ndof re-measured on the NEW |η|-fold pair-η grouping (2026-09-03, D11)

Same χ²/ndof machinery as R3, same `nocorr`/`nocorr_ptmerge` inputs (A/B, unaffected by D11), but
C/D (`nocorr_etamerge`/`nocorr_etamerge_ptmerge`) now read the FRESH closure fill on the sign-
independent |η| fold (barrel <1.0 / [1.0,2.0) / [2.0,2.4)), after the two bugs above were fixed.

Inclusive closure (integral over all cells — a weak, non-differential check, quoted for context):

| version | A `nocorr` | B `nocorr_ptmerge` | C `nocorr_etamerge` (NEW fold) | D `nocorr_etamerge_ptmerge` (NEW fold) |
|---|---|---|---|---|
| Tight, all-OS | 0.99400 | 0.99399 | 0.96937 | 0.96937 |
| Tight, signal | 0.98879 | 0.98878 | 0.95469 | 0.95468 |
| Medium, all-OS | 0.99552 | 0.99552 | 0.97007 | 0.97007 |
| Medium, signal | 0.99051 | 0.99049 | 0.95520 | 0.95520 |

χ²/ndof summed over all 15 crossx pair-p_T bins:

| version | n_dof | A | B | C (NEW fold) | D (NEW fold) |
|---|---|---|---|---|---|
| Tight, all-OS | 135 | 15.47 | 15.33 | **93.15** | **93.28** |
| Tight, signal | 133 | 8.44 | 8.28 | **61.33** | **61.56** |
| Medium, all-OS | 135 | 14.91 | 14.95 | **96.92** | **97.07** |
| Medium, signal | 133 | 7.94 | 8.04 | **64.56** | **64.77** |

Split at the pair-p_T merge edge (Tight, all-OS shown; the other three versions are the same
pattern): bins 1-11 (ndof=99) go from C/D's **15.79** under the retired signed grouping (R3) to
**112.8** under the new fold; bins 12-15 (ndof=36) go from **2.06-2.27** to **39.1-39.6**. The
degradation is dominated by the well-measured low-p_T bins, exactly where R3 already found the
un-normalized D² had been hiding the real (small-but-significant) signed-grouping non-closure —
the fold grouping is worse almost everywhere on the axis, not just in a few cells.

**This REVERSES the R2/R3 ranking.** Under the retired signed grouping, C/D were competitive with
or better than A/B (R3: 12.1-12.2 vs 15.3-15.5 all-OS). Under the new |η| fold, C/D are 4-8x WORSE
than A/B in every version and WP. The physics reason the fold was adopted (D11: the dR correction
barely depends on the SIGN of pair η, measured from the Step-3 FIT quality/χ²-per-dof) does not
carry over to the CLOSURE test, which is a stricter, more differential check of the delivered
correction against the trigger-weighted spectrum — the two questions are related but not the same,
and this result says they disagree here. **Recorded as measurement, not a decision**: whether the
fold grouping should still be adopted (e.g. for reasons independent of this closure metric) is for
the user to weigh against R4, not something this doc resolves on its own.

## Remaining Work

1. **Replace `DrCorrectionCrossxEvaluator` with the general cascade class** once an approach is
   chosen (D9) — two cascade implementations must not coexist longer than this comparison.

## Latest Stage

**2026-09-17 — the pp24 cross-section ADOPTED approach C's |η^pair| fold (`nocorr_etamerge`, 8 × 3
cells, un-merged pair-pT) for coarse pair-pT bins 1–6, with the last two bins served by the
single-value pair efficiency instead of any ΔR fit (user decision;
`docs/tracking/pp24_trig_eff_hybrid_application.md`). The four-approach closure ranking is still to
be regenerated on the 9 GeV axis as recorded below; it no longer gates the application.**

**2026-09-08 — ⚠ THE FOUR-APPROACH CLOSURE FIGURES ARE NOW STALE (two independent reasons).**

1. **The polyu tier of the delivered cascade changed.** `mc_trigger_efficiency.md` R34 constrained
   the Step-3 `polyu_fixedRp` fit so `f(0) ≤ C` (the user's "at ΔR = 0 the efficiency cannot exceed
   the plateau"), and reparametrized it in `A ≡ f(0) − C`. §PP-3 tier 2 therefore moved in every
   cell where `expo` is rejected and `polyu` accepted, and 11 cells additionally have a COLLAPSED
   fitted baseline `C` — 9 rerouted to `interp` by `DrCorrPlateauUsable`, 2 delivered at the polyu
   tier with `f/C` up to 1.94 (`mc_trigger_efficiency.md` R34(f)(2)). All 30
   Step-3 polyu fit files were regenerated 2026-09-08; the approach-comparison figures are from
   01:33 that morning.
2. **The canonical pair-pT axis moved 8 → 9 GeV** (a CONCURRENT, uncommitted workstream adopting
   muon `pT > 4.5 GeV`: `ParamsSet::signal_pair_pt_min`, `pair_pt_coarse_bins`, `pT_bins_150`).
   `run_mc_trigeff_closure.sh` now dies in Stage 2 with `DrCorrectionCascadeEvaluator: pair-pT edge
   0 is 8.000000 in the fit file but 9.000000 canonically -- stale fit file`. **The closure cannot
   be regenerated at all until the Step-3 histograms are refilled on the 9 GeV axis.**

**Consequence for the OPEN decision.** The χ²/ndof ranking below (Tight signal 6.557 A / 6.503 B /
5.029 C / 5.129 D) was the evidence the approach choice awaits judgement on. Do NOT decide on it:
it predates both changes. The correct sequence is — refill Step-3 on the new 9 GeV axis → refit all
three methods → rerun this closure → re-read the ranking.

*The previous Latest Stage follows.*

---


**2026-09-08 — R4 RETRACTED, all four approaches refilled and replotted on the fixed cascade
lookup (see the last Progress Log entry).** Nothing is in flight in this doc. The open item is a
USER JUDGEMENT: with the bug gone the |eta| fold has the BEST chi^2/ndof of the four
(Tight, signal: 5.029 / 5.129 against 6.557 / 6.503; Medium: 5.222 / 5.332 against 6.321 /
6.266), the reverse of what R4 concluded, so the choice of approach for
the pp24 cross-section is reopened on the regenerated evidence. The cross-section application is
still deliberately unchanged (`DrCorrCrossxMode()` = `nocorr_ptmerge`).
