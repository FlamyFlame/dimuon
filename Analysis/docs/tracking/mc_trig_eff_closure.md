# MC Trigger-Efficiency Closure (pp24, 2mu4)

**Mode:** Implementation. **Created:** 2026-08-13. **Session:** "MC trigger efficiency closure test".
**Parent:** `mc_trigger_efficiency.md` (ACTIVE) — this doc executes its **Remaining Work 8**
(the MC closure test), which that doc wrote in round 7 and deliberately did not run.
**Why a separate doc:** a sibling session owns `mc_trigger_efficiency.md` and has uncommitted
round-11 edits in it; two sessions appending to one file would race. Everything this doc measures
is a *consumer* of that doc's Step-3 deliverable — nothing here re-derives it.

---

## Objective

Close the pp24 2mu4 trigger-correction chain on itself, using the MC sample as the only place
where an **unbiased denominator** (every event stored, `StoreAllEvents`) and the **per-pair
trigger decision** both exist.

Take the MC pairs that **pass 2mu4**, weight each by the inverse of the per-pair trigger
probability the analysis actually applies, and check that the result reproduces the **all-pairs**
(no trigger requirement) MC spectrum — **differentially in pair p_T, inside pair-η bins**.

Pb+Pb is explicitly OUT of scope (user): mu4 is a union weight needing the Step-4 single-leg
correction as well, and the full HIJING overlay production is still in flight.

## Autonomy Contract (**DONE** 2026-08-18 — kept for the record, no longer governs)

> Round-1 `Done` was met on 2026-08-13. The round-2 request (apply ε_MC; rerun on current inputs;
> two ε_ΔR variants, one subdirectory each) came later and is met too: 12 PNGs in 6 mode
> subdirectories, both reviews run, doc + INDEX updated. **A post-compaction agent should read
> Latest Stage, not the round-1 `Done` list below** — item 2 there describes the superseded 4-PNG
> / three-line round-1 deliverable.

- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a
  passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. Closure **fill** code (RDF over the MC OS pair tree, applying the per-pair weight) and
     closure **plot** code exist, compile, and are driven by one pipeline script.
  2. **Four deliverable PNGs** at the nominal Tight working point — 2 sample versions
     (all opposite-sign pairs / pairs passing the data-like single-b signal cuts) × 2 pair-p_T
     binnings (canonical 8-bin, 4-bin variant) — each with **one subplot per pair-η bin** and
     **three lines**: no trigger requirement, 2mu4 corrected with the `expo` fit, 2mu4 corrected
     with the `polyu` fit. The Medium working point is produced too (WP registry).
  3. The **4-bin upstream chain** (Step-3 fill → plateau → opposite-sign no-plateau-correction
     fits) regenerated on the current code, because the existing `_pt4bin` outputs predate the
     per-sign and no-plateau-correction work (parent doc R24b).
  4. Reviews (`/review-analysis-code` for the C++/RDF, `/review-plot` for the figures), this doc
     + INDEX updated, commits.
  5. ~~**Re-run on the final inputs** once the sibling session's hand-off note lands.~~ **DONE
     2026-08-13**: the hand-off landed (commit `d043677`, parent doc R29 + Latest Stage); every
     round-1 result was produced at 04:16–04:19, after all of its inputs, and the driver's
     freshness gate enforces it. *(Round-2 results were produced 2026-08-18 11:51–11:55; see the
     round-2 provenance note in Results.)*
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The analysis corrects each surviving muon pair by `1/ε_trig^pair`. For pp24 / 2mu4 that is

```
ε_trig^pair(pair) = ε^nc(p_T1, q·η1) · ε^nc(p_T2, q·η2) · ε_ΔR^2mu4(ΔR)
```

(this is the **analysis's** weight — the closure applies ε_MC in place of ε^nc since round 2, §2/D4)
with ε^nc the **data-derived** single-muon tag-and-probe efficiency and ε_ΔR^2mu4 the
**MC-derived** two-body ΔR correction (parent doc §3.3, Step 3). Every ingredient has been
measured, but the *product* has never been tested against a sample where the true answer is
known. MC is the only sample where it can be: the trigger decision is recorded without selecting
on it, so the "what the yield would have been without the trigger" reference exists.

A closure that holds bin-by-bin validates the factorized form, the ε^nc parameterization, the ΔR
correction *and* its fit. A structured deviation is a systematic on the trigger correction.

### 2. Top-level equation

Closure variable, per (pair p_T, pair η) bin:

```
C(p_T^pair, η^pair) = Σ_{pairs passing 2mu4} w_MC / ε_trig^pair
                      ─────────────────────────────────────────
                          Σ_{all pairs} w_MC
```

`C = 1` in every bin ⇔ the correction closes. `w_MC` is the per-pair MC weight
(σ_slice·ε_filt·r_isospin/N_slice). The plotted lines are **3 in the `nocorr` variant** (the denominator plus the numerator for each
of the two ΔR-fit forms) and **2 in `nocorr_ptmerge`** (the denominator plus the single cascade
series) — see §3.2 — all as dσ/dp_T^pair; the ratio panel is `C`.

**Which ε goes in the weight — ROUND 2 (user, 2026-08-18): the MC-derived ε_MC**, so that the test
is SELF-CONTAINED. Round 1 plotted the as-applied mixed configuration (data ε^nc × MC ε_ΔR) and
measured 1.270; R1 established that the entire offset is ⟨r₁·r₂⟩, the MC/data single-muon
over-efficiency, and NOT the ΔR correction. A closure of that configuration therefore cannot come
out at 1 on an MC sample by construction, and it confounds the two effects. Weighting with ε_MC —
the very efficiency the Step-3 inverse weighting divided by — removes the data/MC difference from
the test entirely, so what remains is exactly the question the closure exists to answer: **does the
ΔR correction, its fit, the cell lookup and the inverse weighting close?**

**The two ε roles are SWAPPED relative to round 1.** ε_MC is now the plotted/applied series; the
data ε^nc numerator is still filled and still printed per pair-η bin, as the DIAGNOSTIC. Their
ratio is unchanged in meaning and is still the measurement R1 reported:

```
C(ε_MC)                 → does the ΔR correction + its fit close?   ← THE PLOTTED SERIES (round 2)
C(ε_data) / C(ε_MC)     → ⟨r₁·r₂⟩,  r_j = ε_MC(leg j)/ε_data(leg j),
                          weighted by w_MC/(ε_MC,1·ε_MC,2·ε_ΔR)     ← printed diagnostic
```

**⚠ It is the mean of the PRODUCT of the two legs' ratios, NOT ⟨r²⟩.** Reading a per-leg number as
`√⟨r₁r₂⟩` is valid only where both legs sample the same q·η regime. For an **opposite-sign** pair
the legs sit in **mirrored** q·η bins (q·η₁ ≈ −q·η₂), and the parent doc's R8/R10/R14 record a
strong q·η asymmetry (the forward-**negative** MC turn-on is saturated, its mirror is not). So the
square-root reading holds in the barrel panels and **not** in the endcap ones.

Both numerators are written to the output ROOT file and printed per pair-η bin; the FIGURE carries
the ε_MC series only.

### 3. Step-by-step method

#### §3.1 The pair sample

**Base selection = the parent doc's Step-3 selection, unchanged** (that is the sample the ΔR
correction was measured on, so the correction is only defined there):

- truth-seeded, reco-matched Pythia muons (parent §3.0, round-5 revert);
- both legs: nominal WP bit, `p_T > 4 GeV`, `|η| < 2.4`;
- both legs: **truth fiducial** `truth p_T > 4 GeV`, `|truth η| < 2.4` (parent §3.0(c'));
- both legs: **forward low-p_T veto** `p_T > 7 || q·η > −2` (parent §3.5 decision rule);
- both legs: **fiducial gap cut** from `ParamsSet::single_mu_fiducial_gap_cuts`
  (currently `(−1.20,−1.05)`, `(−0.10,+0.06)`, `(2.30,2.40)` in q·η);
- **opposite sign only** (`muon_pair_tree_sign2`). Same-sign is excluded by the user because the
  no-plateau-correction fit fails in many same-sign cells and the single-b analysis is an
  opposite-sign measurement.
- `pair p_T` inside the coarse axis `[8, 150] GeV` — outside it no ΔR-correction cell exists.

**Two VERSIONS of the sample:**
- **(A) all opposite-sign pairs** — the base selection alone.
- **(B) data-like single-b signal pairs** — base selection **plus the reco signal cuts**, byte-identical
  to `RDFBasedHistFillingPP.cxx`'s `signal_cuts`:
  `minv > 1.08 && minv < 2.9 && pair_pt > 8 &&` the `ParamsSet::single_mu_fiducial_gap_cuts`
  expression (`ParamsSet::FiducialGapCutExpr`) on **both** muons.
  **Updated 2026-08-17** (`pp24_crossx_rerun_2026_08.md`): the per-muon one-sided `q·η < 2.2` was
  REPLACED by the detector-gap fiducial cut on both muons, in the data analysis and here at once —
  both sides read the same `ParamsSet` vector, so they cannot drift. *(The round-1 text said the
  gap cut was "not yet applied to the data signal selection"; that is no longer true.)*
  The gap cut therefore appears twice — once in the base selection, once here — which is a
  no-op intersection, not a second cut.

**Two SUB-VERSIONS** per version: the ΔR correction taken from the canonical **8**-bin pair-p_T
cells and from the **4**-bin variant (`MCTRIGEFF_PAIRPT_4BIN=1`). The plotted pair-p_T axis is the
SAME binning as the correction's cells in each case — a 1D view of the cells the correction was
measured in (CLAUDE.md §Binnings item 2).

#### §3.2 The per-pair trigger weight

**ROUND 2 — TWO VARIANTS OF ε_ΔR, one output subdirectory each (user, 2026-08-18):**

| variant | ε_ΔR cells | how ε_ΔR is chosen per cell | plotted lines |
|---|---|---|---|
| `nocorr` | canonical 8 pair-pT × 9 pair-η | the `expo` fit and the `polyu_fixedRp` fit, drawn as two SEPARATE series | 3 (no-trigger + 2 corrected) |
| `nocorr_ptmerge` | 7 pair-pT (last two merged, [72.1,150) GeV) × 9 pair-η | ONE series, the **cross-method cascade the pp24 cross-section applies**: `expo` fit → `polyu_fixedRp` fit where the exponential is rejected → the RAW measured bins where both are | 2 (no-trigger + 1 corrected) |

The `nocorr_ptmerge` series is built by `Utilities/DrCorrectionCrossxEvaluator.h` — the SAME class
the cross-section fills with, reading the SAME named defaults (`DrCorrCrossxMethod/Sign/Mode`), so
the closure tests what is actually applied rather than a re-implementation of it.


For a pair with legs (p_T1, q·η1), (p_T2, q·η2), separation ΔR, in cell (i_ptpair, i_ηpair):

```
ε_trig^pair = ε_MC(p_T1, q·η1) · ε_MC(p_T2, q·η2) · ε_ΔR(ΔR ; cell)          ← APPLIED (round 2)
```

- **ε_MC (APPLIED, §2/D4)** — the MC direct conditional probability P[mu4 | reco muon] from
  `single_mu_effcy_pT_fit_mc{WP}.root` (`SingleMuEffEvaluator::Src::kMCDirect`), looked up per
  (charge, coarse q·η bin) and **evaluated continuously at the exact p_T** (`TF1::Eval`, never
  resample-to-nearest), with the p_T clamped into the fit range, capped at 1 and floored at 0.02.
- **ε^nc_data — DIAGNOSTIC NUMERATOR ONLY (§2/D4), never plotted** — the pp24 tag-and-probe
  turn-on TF1s, the identical lookup the analysis and `FillMCTrigEffHists.cxx`'s `DataEffEvaluator`
  use. Its numerator is booked off the same node with the same ε_ΔR, so the ratio of the two is
  ⟨r₁·r₂⟩ and nothing else. Both lookups go through the SAME `SingleMuEffEvaluator` struct, hence
  the same clamp/cap/floor guards — a guard difference would show up in their ratio as if it were
  physics. *(Round 1 applied this one; do not re-swap them from the round-1 text.)*
- **ε_ΔR** — from the **no-plateau-correction** Step-3 fit of the **opposite-sign** series, then
  divided by its own fitted baseline `C`:

  ```
  ε_ΔR(ΔR) = f_cell(ΔR) / C_cell     for ΔR < 1
  ε_ΔR(ΔR) = 1                        for ΔR ≥ 1
  ```

  `f(ΔR) = C + A·exp[−(ΔR/λ)^p]` (`expo`) or the `polyu_fixedRp` analogue; `C` is the fit's own
  free asymptote, determined by the ΔR < 1 data alone (parent doc R27). **Why this variant:** the
  measured ε_ΔR carries structure out to large ΔR, worst in the gap-enclosing pair-η bins, so the
  far-region [2, 3.5] plateau is not necessarily the right baseline for the small-ΔR region the
  correction is about. Dividing by `C` afterwards puts the plateau at 1, which is what the
  factorized weight assumes, and leaves **no step at ΔR = 1** because `f(1)/C → 1` there.
- **Fit domain is ΔR ∈ [0, 1] only** — beyond it no correction is assumed (`ε_ΔR ≡ 1`).

#### §3.3 Cells with no usable fit — TEMPORARY PLACEHOLDER

> **⚠ THIS IS A TEMPORARY PLACEHOLDER, NOT THE FINAL SOLUTION.**
> The parent doc's open finding **R26** (both parametric forms fail badly in a substantial
> minority of cells, and `fit_ok` carries no χ² term) is unresolved and is a *user decision*.
> Until it is taken, a cell whose opposite-sign fit is rejected (`h_step3_fit_ok == 0`) is filled
> with the **raw binned** measured ε_ΔR:
> ```
> ε_ΔR(ΔR) = ε_ΔR^raw(bin containing ΔR) / C_raw ,  ΔR < 1 ;   1 otherwise
> C_raw = unweighted mean of the raw bins with a finite error over ΔR ∈ [0.5, 1.0]
> ```
> `C_raw` exists for the same reason the fitted `C` does — without it the raw values carry the
> cell's own normalization offset and the correction would step at ΔR = 1. A raw bin with no
> information (zero denominator) and a cell with no usable `C_raw` get `ε_ΔR = 1`; both are
> counted and printed. **This is a stop-gap so the closure can be produced at all; it must be
> replaced once R26 is decided** (χ² screen in `usable`, per-cell interpolation, or a fit form
> that can describe the measured shape).

#### §3.4 What is plotted

One PNG per (version, pair-p_T binning, working point, **ε_ΔR variant**). **One subplot per pair-η
bin** of `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` (9 bins → 3×3). Each subplot:

- **upper pad** — dσ/dp_T^pair (log–log). `nocorr`: three series — (1) all pairs, no trigger
  requirement; (2) pairs passing 2mu4, corrected with the `expo` ΔR fit; (3) same, `polyu` ΔR fit.
  `nocorr_ptmerge`: two series — (1) and (2') the single cascade-corrected series of §3.2.
- **lower pad** — the closure ratio `C` (each corrected series over series 1), with a line at 1.

The two variants go to **separate subdirectories** named by the repo's own mode helper
`DrCorrPlateauModeDir`: `closure/no_plateau_correction/` and
`closure/no_plateau_correction_last2ptbins_merged/`. The pair-p_T axis of each figure IS that
variant's own cell axis (D2), so the merged variant is drawn on 7 bins.

Statistical errors on `C` are the **conditional (binomial-correct)** form, not `TH1::Divide`'s
independent propagation: the numerator is a re-weighted SUBSET of the denominator (parent doc
R12 / `dr_correction_ratio.h`). With `a_i = w_i/p̂_i`, `p̂_i = ε_trig,i` the **predicted** per-pair
probability and `π_i` the **true** one,

```
Var(N) = Σ_all a²·π(1−π)  →  estimated from the fired pairs as  Σ_fired a²(1−π)
π_i ≈ R·p̂_i  (first order; R = the bin's closure ratio)
⇒  Var = A − R·B ,   A = Σ_fired a² ,   B = Σ_fired a²·p̂ ,   σ_C = √Var / D
```

**The `R` factor is not droppable here.** `Var = A − B` is the special case `R = 1`, i.e. it
assumes the closure already holds — the very thing being tested. Since `B ≤ A`, dropping `R` shifts
the variance by `(R−1)·B`, so it **over**states σ where `R > 1` and **under**states it where
`R < 1` — the second direction being the dangerous one, and common in the statistics-poor
high-pair-pT cells.
**Round-2 numbers (applied ε_MC, `nocorr`/`expo`/all-OS/Tight):** inclusive `R = 0.994`, but per
cell `R` spans **0.521–1.382**, and dropping the factor overstates σ by up to **1.11×** and
understates it by up to **2.98×**; **0 of 72** cells reach the binomial `k = n` boundary and 0 have
`A − B ≤ 0`.
*(The round-1 values this paragraph originally quoted, on the data-ε configuration: inclusive offset
~27 %, overstatement up to 1.32× at `R = 1.52`, understatement up to 2.46× at `R = 0.82`, 1 cell of
72 at the `k = n` boundary where a naive implementation would set the error to 0 and both consumers
would drop the point.)* The computation is therefore delegated to
`SetConditionalRatioErrors` in `dr_correction_ratio.h` — the repo's **single** implementation,
which carries that boundary fallback (`max(R,1)/n_eff`).

### 4. Negative constraints

- **NO trigger requirement anywhere in the denominator.** The denominator is every selected pair.
- **The APPLIED ε is ε_MC (§2 / D4, round 2); the data ε^nc numerator is the DIAGNOSTIC and must
  NOT be plotted.** *(Round-1 history, superseded: the constraint here used to read the opposite —
  "do NOT weight with ε_MC" — because round 1 plotted the as-applied mixed configuration. R1 showed
  that configuration cannot close at 1 on an MC sample by construction, its offset being ⟨r₁·r₂⟩,
  so the user swapped the roles. Do not re-swap them from the round-1 text.)*
- **Do NOT apply `ε_ΔR^single` (Step 4) to the pp weight** — the 2mu4 product already absorbs the
  single-leg ΔR dependence; adding it would double-count (parent §4).
- **Do NOT re-derive the sample selection from raw NTUPs.** The closure reads the ntuple-processing
  output (`*_mc_trig*.root`) through the same RDF selection strings as `FillMCTrigEffHists.cxx`.
- **Do NOT change any binning.** The pair-p_T cells come from `MCTrigEffPairPt::Edges`, the pair-η
  cells from `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` — read, never retyped.
- **Do NOT touch** `mc_trigger_efficiency.md`, `run_dr_correction_fits.sh`,
  `plot_dr_correction_fits.cxx` or the 8-bin outputs: the sibling session owns them.

---

## Scope

**In:** pp24 fullsim FULL sample (`pp_full`), 2mu4, opposite sign; the closure fill + plot code;
the 4-bin upstream regeneration needed to have a 4-bin opposite-sign no-plateau-correction fit.

**Out:** Pb+Pb / the HIJING overlay (user); same-sign closure (user); resolving R25 (sign-dependent
correction) or R26 (fit-form failure); applying anything to crossx.

## Design Decisions

### D1: opposite-sign fits, not the sign-integrated ones (2026-08-13)
The closure sample is opposite-sign only, so the correction that belongs on it is the
**opposite-sign** Step-3 series (`_os`). The sign-integrated series is numerically close (87 % of
pairs are opposite sign, parent R25) but it is not the same measurement, and the user's stated
reason for restricting the closure to opposite sign — the same-sign fits fail — only makes sense
for the sign-separated fits.

### D2: the pair-p_T plot axis IS the correction's cell axis (2026-08-13)
Plotting the closure on a finer pair-p_T axis than the cells the correction was measured in would
mix two binnings of the same quantity — exactly the failure CLAUDE.md §Binnings item 2 exists to
prevent. So the 8-bin sub-version plots on the 8 canonical log bins and the 4-bin sub-version on
the 4-bin variant, both read from `MCTrigEffPairPt::Edges`.

### D3b: the shared headers were extracted but `FillMCTrigEffHists.cxx` was NOT migrated onto them (2026-08-13)
`MCTrigEffPairSelection.h` (the Step-3 pair selection) and `SingleMuEffEvaluator.h` (the fitted
single-muon turn-on lookup) were extracted **verbatim** from `FillMCTrigEffHists.cxx`, which still
carries its own inline copies. Normally the extraction would be completed by pointing that file at
the headers — one definition, no drift.
**Why it was not done:** a concurrent session is re-running the trigger-efficiency chain and
executes `FillMCTrigEffHists.cxx` at unpredictable times; editing it mid-run would break that run.
**Why the duplication is the highest-consequence kind here:** if the two copies ever drift, the ΔR
correction is applied to a *different* sample from the one it was measured on — silently, since
every histogram still fills. Mitigation in place: a MIRROR NOTICE in both headers naming the
original file and lines, and the extracted string was verified **byte-identical** at extraction
time (a small ROOT program compared the concatenated selection strings). See Remaining Work 3.

### D4: ROUND 2 — the applied single-muon efficiency is ε_MC, not the data ε^nc (user, 2026-08-18)
The closure is now **self-contained**: every ingredient of the weight comes from the same MC
sample, so the residual is the ΔR correction alone. Round 1's mixed configuration is not
discarded — its numerator is still filled and printed as the diagnostic, and its ratio to the
ε_MC series is still ⟨r₁·r₂⟩ (R1). Physics Procedure §2 was updated with this change, which is a
user instruction, not an agent choice.
**Consequence for reading the numbers:** the round-1 headline (1.270) and the round-2 headline are
different measurements, not a revision of one another — 1.270 was the as-applied data-ε closure,
round 2 plots what round 1 called the ε_MC diagnostic.

### D5: the merged variant is built by the CROSS-SECTION's own evaluator, not a copy (2026-08-18)
`nocorr_ptmerge` is what pp24 crossx applies, and the point of closing on it is to test *that*
object. So the closure instantiates `Utilities/DrCorrectionCrossxEvaluator.h` itself — same class,
same `DrCorrCrossxMethod/Sign/Mode` defaults, same expo→polyu→raw cascade and the same binning
guard. A second implementation of the cascade inside the closure could drift from the cross-section
silently (every histogram would still fill), which is the failure mode D3b already warns about for
the pair selection. The macro additionally THROWS if the requested mode is `nocorr_ptmerge` while
`DrCorrCrossxMode()` has moved elsewhere, so the two cannot diverge unnoticed.

### D3: forward low-p_T veto kept in the closure sample (2026-08-13)
`p_T > 7 || q·η > −2` is part of the sample the ΔR correction was measured on. Dropping it for the
closure would apply the correction to a population it was not measured on, and would re-admit the
forward-endcap region where ε_MC is badly parameterized (parent R8/R10/R14).

---

## Implementation Plan

1. [x] **Closure fill** `RDFBasedHistFilling/FillMCTrigEffClosure.cxx` (per §3.1–§3.3):
   RDF over `muon_pair_tree_sign2`, two versions, two ΔR-fit methods; writes
   `mc_trig_eff_closure_<label><wp><ptbin>.root` with (pair p_T × pair η) TH2Ds
   `den / num_<method> / numA_<method> / numB_<method>`. → `/review-analysis-code`.
2. [x] **Closure plot** `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_closure.cxx`
   (per §3.4) → `<out_base>closure/`. → `/review-plot`.
3. [x] **Driver** `pipelines/run_mc_trigeff_closure.sh`, including the 4-bin upstream
   regeneration (Step-3 fill → `plot_mc_trig_eff` plateau → `fit_dr_corrections` os/nocorr).
4. [x] **Run on current inputs**, sanity-check, record numbers here.
5. [x] **Re-run after the sibling's gap-cut hand-off** — hand-off landed (`d043677`); all four
   fills + eight PNGs regenerated at 04:16–04:19 on the final inputs.

### Round 2 (2026-08-18, user request)

6. [x] **Fill**: applied ε → `Src::kMCDirect` (data ε^nc demoted to the diagnostic numerator);
   `plateau_mode` becomes an argument; `nocorr_ptmerge` served by `DrCorrectionCrossxEvaluator`;
   output file gains `DrCorrPlateauModeTag`. → `/review-analysis-code`.
7. [x] **Plot**: `plateau_mode` argument → mode subdirectory + mode-dependent series list
   (3 lines / 2 lines); header text names ε_MC and the cascade. → `/review-plot`.
8. [x] **Driver**: loop the modes (8-bin: both; 4-bin: `nocorr` only — the merged mode is defined
   for the canonical axis only and THROWS otherwise), freshness-check each mode's own fit files.
9. [x] **Run** on the current inputs, sanity-check, record numbers.

## Progress Log

*(append-only; newest entries at the END)*

- 2026-08-13 — Doc created. Physics Procedure written from the user's specification plus the
  parent doc's §2/§3.3/§4 and R27 (no-plateau-correction fit variant). Plan recorded before any
  code was written.
- 2026-08-13 — Steps 1-3 DONE. New files, all of them NEW (nothing owned by the sibling session
  was touched):
  * `Analysis/Utilities/MCTrigEffPairSelection.h` — the Step-3 pair selection + the data-like
    single-b reco signal cuts as one builder.
  * `Analysis/Utilities/SingleMuEffEvaluator.h` — one evaluator for BOTH fitted turn-on sets
    (data tag-and-probe / MC direct), so the two carry identical clamp/cap/floor guards.
  * `Analysis/plotting_codes/trig_effcy/mc_based/dr_correction_apply.h` — the consumer side of
    `fit_dr_corrections`: `eps_dR(dR) = f(dR)/C` for dR < 1, 1 above, with the §3.3 placeholder.
  * `Analysis/RDFBasedHistFilling/FillMCTrigEffClosure.cxx`
  * `Analysis/plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_closure.cxx`
  * `Analysis/pipelines/run_mc_trigeff_closure.sh`
  Two defects found and fixed during the trial run: (i) pair eta is the pseudorapidity of the SUM
  4-vector and is NOT bounded by the two legs' eta, so 1.1 % of selected pairs sat outside the
  |eta^pair| < 2.4 cell grid and were silently uncorrected — the sample now carries an explicit
  pair-eta bound; (ii) every panel after the first was created INSIDE the previous panel's ratio
  pad (a new TPad is adopted by gPad), so only one of nine panels appeared — `c.cd()` per panel.
- 2026-08-13 — Step 4 DONE: the 4-bin upstream was regenerated (its `_pt4bin` outputs predated
  both the per-sign booking and the no-plateau-correction variant, parent doc R24b) and all four
  deliverable PNGs + their Medium counterparts were produced. Results in R1-R3 below.

- 2026-08-18 — **ROUND 2 DONE** (user request: apply ε_MC; rerun on current inputs; two ε_ΔR
  variants in one subdirectory each). Plan items 6–9 complete.
  * `FillMCTrigEffClosure.cxx`: applied ε → `Src::kMCDirect`; the data ε^nc numerator demoted to
    the diagnostic; new `plateau_mode` argument; new `SeriesEval` list (one entry per corrected
    numerator) so neither variant needs its own code path; `nocorr_ptmerge` served by
    `Utilities/DrCorrectionCrossxEvaluator.h` itself (D5); histograms renamed by which ε went into
    them (`num_epsmc_*` applied / `num_epsdata_*` diagnostic) so files written on either side of
    the D4 swap can never be confused; output file gains `DrCorrPlateauModeTag`.
  * `Utilities/DrCorrectionCrossxEvaluator.h` (shared with the pp24 crossx fill): per-cell route
    census exposed as members and stamped into the closure's provenance; `CheckCanonicalBinning`
    now edge-checks the BACKUP fit file's axes too (it previously compared only bin counts, so a
    stale backup with the same cell count but different edges would have corrected pairs with
    another cell's curve, silently). Both changes are additive — verified they cannot alter or
    newly throw for the cross-section (both WPs' backup fit files carry the canonical merged edges).
  * `plot_mc_trig_eff_closure.cxx`: mode argument → mode subdirectory + mode-dependent series list;
    header names the MC turn-on and, for the merged variant, the cascade.
  * `pipelines/run_mc_trigeff_closure.sh`: MODES loop, per-mode fit-file pre-flight and freshness
    gate, per-mode plot subdirectory; 4-bin × `nocorr_ptmerge` skipped loudly (undefined by design).
  * Reviews: `/review-analysis-code` ran 3 iterations. All findings were doc-consistency or
    guard-completeness, none changed a number; every reported value was independently re-read from
    the ROOT files and matched. The recurring finding is worth recording as a lesson: **after a
    role swap like D4, the stale statements are spread across the whole Physics Procedure** — §4,
    then §3.1(B), then §3.2's own equation each had to be caught in a separate pass.
  * `/review-plot` PASS (2 runs — the first reviewer died on an API error mid-pass and was
    relaunched). It confirmed the layout matches the request exactly, that no internal method name
    (`expo`, `polyu_fixedRp`, `crossx_cascade`) reaches any canvas, that the ratio pads' percentile
    range always contains 1 with every off-scale point counted (printed count exact in all 12
    figures, and 0 cells with a zero numerator that could hide an uncounted one), and that all 108
    spectra fall smoothly with no unexplained discontinuity. It ruled the merged figure's line
    defining which `f` each cell uses to be a legitimate definition rather than prose — without it
    the symbol `f` would be undefined on that canvas.
  * **Superseded round-1 artefacts QUARANTINED, not deleted** (they carried the same file names one
    directory above the new ones): 8 PNGs → `<plot dir>/closure/superseded_round1_data_eps/`,
    4 untagged ROOT files → `<mc dir>/superseded_round1_data_eps/`, 8 round-1 logs →
    `pipelines/logs_closure/superseded_round1_data_eps/`. A glob-based figure sweep now sees 20
    closure PNGs, 8 of them superseded — delete the quarantine directories once the user confirms.

## Results & Observations

> **R1–R3 below are ROUND 1** — the AS-APPLIED configuration (data ε^nc × MC ε_ΔR), which round 2
> demoted to the diagnostic. The round-2 results (applied ε_MC, both ε_ΔR variants) are **R4–R6**.
> Where round 2 re-measured a round-1 number on the current inputs, R4 says so.

*(**Round 1**; the round-2 provenance is one paragraph below.) Measured on the **FINAL
post-gap-cut inputs**. The sibling session's round-11 rerun completed at
03:2x on 2026-08-13 and its hand-off note (`mc_trigger_efficiency.md` Latest Stage / R29) declares
those products final; the numbers below come from closure fills at **04:16–04:19**, i.e. after
every one of them. Verified input mtimes: MC Step-3 hists 02:51 (8-bin) / 02:59–03:03 (4-bin),
MC single-muon turn-ons 02:50, data tag-and-probe turn-ons 02:17 (Tight) / 02:39 (Medium),
opposite-sign no-plateau-correction ΔR fits 03:05–03:18. The driver enforces this: `val_fresh`
refuses a closure output older than any input it consumed.*

***Round-2 provenance:*** *fills 2026-08-18 11:51–11:55, on ΔR fits of 2026-08-17 23:13
(`expo_os_nocorr`, Tight) / 2026-08-13 03:17 (Medium) and 2026-08-18 09:46–09:52
(`_nocorr_ptmerge`, both WPs), MC single-muon turn-ons 08-13 02:50, Step-3 hists 08-13 02:51. Each
output's `provenance` TNamed carries its own inputs' mtimes, and `val_fresh` refuses an output
older than any of them.*

*(The hand-off also flags that closure files written at 02:53/02:58 — during the overlap — mixed
new Step-3 hists with pre-widening ΔR fits. Those were superseded by three later full re-runs;
nothing from that window survives.)*

### R1. ★ THE HEADLINE — the chain closes to ~1 % on eps_MC, and the ~27 % offset in the
### as-applied (data-eps) configuration is the known MC/data single-muon over-efficiency

Inclusive closure, pp24 FULL, opposite sign, Tight, canonical 8 pair-pT bins:

| sample version | uncorrected 2mu4 / all | corrected, ε^nc **data** | corrected, **ε_MC** (diagnostic) |
|---|---|---|---|
| all opposite-sign | 0.6116 | **1.2699** (exponential) / 1.2696 (polynomial) | **0.9940** / 0.9936 |
| single-b signal region | 0.5870 | **1.2678** / 1.2678 | **0.9924** / 0.9923 |

Medium is the same picture (1.2862 / 0.9955 and 1.2858 / 0.9939); the 4-bin variant agrees to
better than 0.5 % (1.2687 / 0.9927 and 1.2629 / 0.9881).

*(Source, as quarantined on 2026-08-18:
`pipelines/logs_closure/superseded_round1_data_eps/fill_{8,4}bin_{tight,medium}.log` and
`<mc dir>/superseded_round1_data_eps/mc_trig_eff_closure_pp24_full{_medium_wp}{_pt4bin}.root`,
fills run 2026-08-13 04:16–04:19 (the "03:35–03:39" written here originally was the earlier,
superseded pass). An earlier
version of this table was written from a fill made against ΔR fit files that the sibling session
rewrote minutes later — the values drifted by 0.06–0.35 %. The driver now refuses to accept a
closure output older than any of its efficiency inputs, and the output's `provenance` string
records each fit file's path, mtime and size.)*

**Read it as two separate statements, because a single closure number confounds them:**

1. **The dR correction, its fit, the inverse weighting and the cell lookup are self-consistent to
   ~1 %.** Weighting with eps_MC — the very efficiency the Step-3 inverse weighting divided by —
   the closure is 0.993 inclusively and 0.94-1.03 across the nine pair-eta bins. **This is the
   statement the user's question ("test the Step-3 dR fit") actually asks for, and it PASSES.**
2. **The as-applied configuration (data ε^nc × MC ΔR correction) over-corrects by ~27 %**, and
   that offset is NOT the ΔR correction. Per pair-η bin, the ratio of the two columns is
   ⟨r₁·r₂⟩ (r_j = ε_MC/ε_data of leg j — **not** ⟨r²⟩; see §2):

   | η^pair | data-ε closure | ε_MC closure | ratio = ⟨r₁·r₂⟩ |
   |---|---|---|---|
   | [-2.4,-2.0) | 1.114 | **0.940** | 1.185 |
   | [-2.0,-1.5) | 1.120 | 0.994 | 1.127 |
   | [-1.5,-1.0) | 1.104 | 0.980 | 1.127 |
   | [-1.0,-0.5) | 1.281 | 0.984 | 1.302 |
   | [-0.5,+0.5) | 1.517 | 1.033 | 1.469 |
   | [+0.5,+1.0) | 1.464 | 0.992 | 1.476 |
   | [+1.0,+1.5) | 1.097 | 0.975 | 1.125 |
   | [+1.5,+2.0) | 1.142 | 0.989 | 1.154 |
   | [+2.0,+2.4) | 1.208 | 0.977 | 1.237 |

   The structure is **barrel-large, endcap-small**, which is the same structure as the parent
   doc's independently measured MC/data single-muon over-efficiency (**barrel-only ~1.2, endcap
   ~1.01**, R3/R8). In the **barrel** both legs of an opposite-sign pair sample the same q·η
   regime, so the per-leg reading √1.47 = **1.21** is valid and matches the parent doc's ~1.2. In
   the **endcap** panels the two legs sit in *mirrored* q·η bins, one of which carries the
   forward-negative anomaly and the other not, so ⟨r₁r₂⟩ = 1.13 is a geometric mean across two
   different regimes and **√1.13 = 1.06 is NOT a per-leg statement** (see §2). Either way, the
   closure figure is measuring the MC trigger simulation being more efficient than the detector,
   not a defect of the correction.

**Consequence.** A closure test of the AS-APPLIED weight cannot come out at 1 on an MC sample
whose own trigger efficiency differs from data's — the deviation is ⟨r₁·r₂⟩ by construction. What the test CAN establish, and does, is (a) the correction machinery is
self-consistent at the 1 % level and (b) a quantitative, pair-eta-differential map of the residual
mis-closure the analysis inherits from the MC/data efficiency difference. Whether that residual
belongs in `docs/systematic_uncertainties.md` as a trigger systematic is a **user decision**.

### R1b. ★ THE CLOSURE TEST'S PRINCIPAL FINDING — the correction over-predicts the 2mu4 firing
### probability by ~2× in the forward-NEGATIVE cell at high pair p_T, and it is the ε map, not ε_ΔR

Found by the plot review (2026-08-13). In `η^pair ∈ [−2.4,−2.0) × p_T^pair ∈ [72, 104) GeV` the
closure collapses to **0.535 ± 0.071** (exponential) / 0.737 ± 0.102 (polynomial), against
1.118 ± 0.062 and 0.852 ± 0.158 in the neighbouring p_T bins — 6.5σ from 1 and ~6σ from its own
neighbour, far outside its error bar, so C1's "not just statistics" clause applies.

**It is NOT the MC/data over-efficiency of R1.** The ε_MC diagnostic in the same cell reads
**0.521**, i.e. the failure survives when the sample's *own* efficiency is used. Whatever is wrong
is inside the MC-internal chain, not in the data-vs-MC difference.

**The MIRROR bin settles what it is.** Compare `η^pair ∈ [−2.4,−2.0)` with `[+2.0,+2.4)` in the
same p_T bin — same acceptance, same statistics, opposite q·η sign:

| | N pairs | ⟨ΔR⟩ | ΔR < 1 | fired 2mu4 | C(data ε) | C(ε_MC) |
|---|---|---|---|---|---|---|
| `η^pair ∈ [−2.4,−2.0)` | 141 | 0.335 | 91.5 % | **46.8 %** | **0.535** | **0.521** |
| `η^pair ∈ [+2.0,+2.4)` | 138 | 0.230 | 94.2 % | 58.0 % | 1.100 | 1.081 |

The two cells are kinematically and statistically alike, the ΔR correction is fully active in both
(>90 % of pairs at ΔR < 1), and the measured firing fractions differ by only 1.24× — yet the
*corrected* results differ by 2.1×. So the correction is over-predicting the firing probability in
the forward-negative cell: with C(ε_MC) = 0.52 at a measured 46.8 % firing rate, the ε map is
asserting a pair probability near 0.9 where reality delivers 0.47.

**Attribution: the documented forward-negative r16578 anomaly, at high p_T where the veto does not
reach.** Parent doc R8/R10/R14: in `q·η ∈ (−2.4,−2.0)` the MC single-muon efficiency is *saturated*
at ~0.90 from the first p_T bin while data rises 0.45 → 0.92, and the sanity check ruled out bad
muons/pile-up — it tracks the r16578 trigger **configuration**. The `p_T > 7 || q·η > −2` veto
(§3.1, parent §3.5) removes only the *low*-p_T forward-negative corner; these legs are high-p_T and
are exactly what the veto deliberately keeps. The mirror bin, which has no such anomaly, closes
normally. This is the ε_MC(p_T, q·η) map being wrong there, **not** a defect of the Step-3 ΔR fit.

**Status: OPEN, and it is a genuine result of the test** — this is the kind of failure an MC
closure exists to expose, and it was invisible in every plot the chain produced before. **Do not
act on it from these numbers**: N = 141 pairs, and the whole chain is re-run on the post-gap-cut
inputs (Remaining Work 2). Re-check it there first; if it survives, it belongs with the parent
doc's forward-endcap thread (R8/R10/R14) and argues that the forward veto's p_T threshold, or the
forward-negative q·η bin itself, needs revisiting for the high-pair-p_T region.

Related and consistent: the ε_MC closure drifts upward with pair p_T in the barrel panels
(`η^pair ∈ [+0.5,+1.0)`: 0.931 → 1.224 → 1.382 over the top three p_T bins), which R2's generic
"statistics-starved" wording does not explain either. Same re-check applies.

### R2. Differential behaviour

Both fit forms agree closely **inclusively** — 1.26989 vs 1.26957, i.e. the choice between the
exponential and the polynomial ΔR-correction form moves the integrated closure by well under 1 % — consistent with the parent doc's R26
finding that the two forms differ mainly in how badly they fit, not in the delivered correction.
**Differentially they do NOT always agree**, and that is what the figure shows: 17 cells differ by
more than 1.5σ, the worst reaching 57 % relative (`pT_pair` bin 8 × `η_pair ∈ [-1.0,-0.5)`:
polynomial 0.795 ± 0.142 against exponential 1.867 ± 0.395). Those disagreements are a direct
read-out of R26 — where the two forms fit the same cell differently, the delivered correction
differs — and they are the reason the closure is drawn for both forms rather than one.
The remaining visible outliers sit in the highest pair-pT bin, where the cells are
statistics-starved (parent doc R23 has those cells' plateaus 20-42 % high and open) — but note that
"statistics" was NOT the explanation for the worst of them; see the fitted-baseline finding
immediately below.

**A pathological fitted baseline was found and screened out (2026-08-13, review iteration 2).**
In `polyu_fixedRp`, cell `pT_pair` bin 8 × `η_pair ∈ [1.5, 2.0)`, the free baseline came out
**C = 0.0534 ± 3.18** — a normalization 20× below unity, with an error 60× the value — yet the cell
carried `fit_ok = 1` (parent doc R26: `usable` has no χ² term). Dividing by it inflated
`f(ΔR)/C` to **21.1** over ΔR < 1, which `kMaxCorr = 5` then truncated, and drove the delivered
closure point in that cell to **0.193 ± 0.031** — 26σ from 1, an artefact of the fit, not a
measurement. **Fixed** by applying the repo's own unusable-normalization screen
`DrCorrPlateauUsable(C, σ_C)` to the fitted baseline — the same test, for the same stated reason,
that `dr_correction_sample_cfg.h` applies to a measured plateau ("a plateau far below 1 is not a
normalization offset but an empty cell; dividing by 0.01 inflates the curve 100× rather than
normalizing it"). The cell now falls through to the raw-bin placeholder and that point reads
**0.737**. Rejections by this screen, per (method, binning), **`nocorr` mode**: `expo` **0** in both 8-bin
runs, **1** in each 4-bin run; `polyu_fixedRp` **1** in every run. *(The merged mode rejects one
`expo` cell as well — the forward-negative top cell; see R5.)* Every rejection is printed by name and the count
travels in the output file's `provenance`. Inclusive numbers move by < 1e-4.

**★ OPEN — the screen bounds the DENOMINATOR, not the delivered correction.** `DrCorrPlateauUsable`
asks "is C a sane normalization?"; it does not ask "is `f/C` a sane correction?". A cell that only
just clears it still delivers a large correction: `polyu_fixedRp`, `pT_pair` bin 8 ×
`η_pair ∈ [-1.0,-0.5)`, has **C = 0.553 ± 0.339** (a 61 % relative error) against **C = 1.200 ±
0.022** for `expo` in the *same* cell, and its ε_ΔR reaches **2.74**. Delivered ε_ΔR now spans
**[0.061, 1.354] for `expo`** and **[0.078, 2.736] for `polyu_fixedRp`** (8-bin Tight; the extrema
and their cells are printed on every run). A close-by 2mu4 correction is a **loss** — the two muons
share an L1 RoI — so ε_ΔR > 1 at small ΔR is a fit artefact rather than a measurement, and that
cell is the 2.6σ split between the two methods visible in the top pair-pT bin of the
`-1.0 < η^pair < -0.5` panel (polynomial 0.795 ± 0.142 vs exponential 1.867 ± 0.395).
**Nothing is rejected on this**: a hard bound on the correction itself is a physics choice
entangled with the parent doc's OPEN R26 (both fit forms fail in a minority of cells; `usable`
carries no χ² term), so it is measured, printed and raised here rather than chosen silently.

**One named exception in the ε_MC diagnostic: `η^pair ∈ [-2.4,-2.0)` closes at 0.940**, a −6 %
deficit well outside the ~1 % of every other bin. It is **not** simply forward statistics: the
mirror bin `[+2.0,+2.4)` closes at 0.977. The natural candidate is the documented
forward-**negative** q·η anomaly of the r16578 production (parent doc R8/R10/R14), which the
`p_T > 7 || q·η > −2` veto only removes below 7 GeV, so the high-p_T forward-negative legs the veto
deliberately keeps are still in this bin. **Left OPEN and named here rather than attributed** — it
would need the same investigation as R23, and it is the one place where the "self-consistent to
~1 %" headline does not hold.

### R3. Coverage and the placeholder

- **66.6 % of selected opposite-sign pairs (Tight) fall inside the correction cells** (pair pT in
  [8,150), |eta^pair| < 2.4). The rest are below 8 GeV — consistent with the parent doc's R27
  (66.6 % of opposite-sign pairs inside the 8-150 GeV axis) — and no correction is defined there.
- **The §3.3 raw-bin placeholder is rarely reached**: 8-bin Tight, 2 of 72 cells for the
  exponential form and 4 of 72 for the polynomial, covering 0.075 % / 0.020 % of the evaluated
  pairs; cells with neither a fit nor a usable raw curve: 2 (21 pairs) for `expo`, 3 (54 pairs)
  for `polyu_fixedRp`, ε_ΔR = 1 there. 4-bin: 2 of 36 (`expo`) / 1 of 36 (`polyu_fixedRp`).
  The per-cell census is stamped into each output file's `provenance` string, so a consumer that
  opens only the ROOT file can still see how much of the map is a placeholder.
  It is nevertheless **a stop-gap and must be replaced once parent-doc R26 is decided.**
- Guard counters clean: 0 efficiency floorings, 0 pairs outside the cell grid. After the
  fitted-baseline screen the polynomial form no longer reaches the ε_ΔR cap of 5. The predicted
  per-pair probability `p = ε₁ε₂ε_ΔR` exceeds 1 for 17 (exponential) / 37 (polynomial) of 1.34 M
  triggered pairs — counted and printed, not clamped, because clamping `p` would silently alter
  the weight as well as the error term.

### R4. ★ ROUND 2 — the SELF-CONTAINED closure holds at ~1 %, and the two ε_ΔR variants agree

Applied ε = ε_MC (§2/D4). Inclusive C = Σ num / Σ den over all cells, pp24 FULL, opposite sign:

| ε_ΔR variant | sample version | WP | **C (applied, ε_MC)** | C (ε^nc data diagnostic) |
|---|---|---|---|---|
| `nocorr`, exponential | all OS | Tight | **0.9940** | 1.2699 |
| `nocorr`, polynomial | all OS | Tight | 0.9936 | 1.2696 |
| `nocorr_ptmerge`, cascade | all OS | Tight | **0.9940** | 1.2699 |
| `nocorr`, exponential | signal | Tight | **0.9888** | 1.2630 |
| `nocorr_ptmerge`, cascade | signal | Tight | 0.9888 | 1.2630 |
| `nocorr`, exponential | all OS | Medium | 0.9955 | 1.2862 |
| `nocorr_ptmerge`, cascade | all OS | Medium | 0.9955 | 1.2862 |
| `nocorr`, exponential (4-bin) | all OS | Tight | 0.9927 | 1.2687 |

Uncorrected 2mu4/all = 0.6116 (all OS) / 0.5849 (signal), Tight.

**Three statements the table supports.**
1. **The correction machinery closes.** With the sample's own ε the chain reproduces the
   no-trigger spectrum to 0.6 % inclusively (1.1 % in the signal region), across both WPs and both
   pair-pT binnings. This is round 1's ε_MC diagnostic promoted to the figure. For **`all_os` it is
   numerically unchanged** from round 1 (0.993978 both rounds — which also shows the 2026-08-17
   re-fit of the `expo`/`os`/`nocorr` file was numerically identical). The **`signal` version DID
   move**, in every configuration, because its selection changed (§3.1(B)): applied closure
   Tight 8-bin 0.99244 → 0.98877, Medium 8-bin 0.99387 → 0.99051, Tight 4-bin 0.98809 → 0.98438,
   and the uncorrected 2mu4 fraction 0.5870 → 0.5849.
2. **The merged variant the cross-section applies is not a different answer — with one measured
   exception.** `nocorr_ptmerge` agrees with the un-merged `expo` series to 1.3×10⁻⁵ inclusively,
   and **53 of the 54 cells the merge does not touch are bitwise identical**. The exception is
   `p_T^pair ∈ [50, 72) × η^pair ∈ [−1.5,−1.0)`, where they differ by **9.7 %** (1.010 un-merged vs
   1.108 merged). It is **not a merge effect**: `expo` is rejected in that cell in BOTH modes
   (`h_step3_fit_ok = 0`), so the un-merged `expo` series falls to the raw-bin placeholder there
   while the cascade routes the cell to the accepted `polyu_fixedRp` fit. It is a direct read-out of
   what the cross-method cascade buys over a single-form series — the same R26 story.
3. **The data-ε diagnostic reproduces R1** (1.2699 all-OS). Its **signal-region** value moved
   1.2678 → 1.2630, for the same selection change as in point 1. The all-OS values are untouched,
   as expected — that version never had the signal cuts.

**Per pair-η bin** the closure spans **0.940–1.033** for `all_os` and **0.868–1.033** for the
signal region (both Tight, `nocorr`/`expo`). The signal region's top pair-pT cells run further from
1 than the all-OS ones — `C = 2.058 ± 0.429` at η^pair ∈ [+0.5,+1.0) and `1.715 ± 0.597` at
[−1.0,−0.5), against 1.382 and 1.347 in `all_os` — the same cells and the same direction as the
barrel high-pT drift of R5, with error bars that span the excursion. Any statement about the
residual must therefore say which sample version it refers to.

**Caveat on "closes to ~1 %".** That is an INCLUSIVE statement. Per cell the mis-closure is often
many σ, because the statistical errors are small: 18.1σ at η^pair ∈ [−0.5,0.5) × p_T ∈ [11.5,16.6)
(C = 1.034 ± 0.002) and ~16σ in two more cells. If this residual is ever written up as a trigger
systematic it must be quoted **per cell**, not as the inclusive 1 %.

### R5. The forward-negative top cell (R1b) SURVIVES round 2 — and in the merged variant it is the raw-bin placeholder

| cell (Tight, all OS) | C (ε_MC) | mirror bin | fired/all |
|---|---|---|---|
| `nocorr`: η^pair ∈ [−2.4,−2.0) × p_T ∈ [72.1,104) | **0.5213** | 1.0807 (η ∈ [+2.0,+2.4)) | 0.49 |
| `nocorr_ptmerge`: η^pair ∈ [−2.4,−2.0) × p_T ∈ [72.1,150) | **0.5341** | 1.1278 | 0.53 |

Round 1 measured 0.521 for the same cell in its ε_MC diagnostic, so the finding is **unchanged on
the current inputs** and is now the *plotted* configuration rather than a printed diagnostic.
Two things round 2 adds:

- **Merging does not rescue it.** The neighbouring [104,150) cell closes at 0.84, and the merged
  cell lands at 0.534 — so this is not a top-bin statistics artefact.
- **In the merged variant that cell is the ONE cell the cross-section's cascade routes to the
  RAW-BIN PLACEHOLDER**: both fits are rejected there (`expo` C = 0.490 ± 0.716, `polyu_fixedRp`
  C = 0.4905 ± 0.400, both failing `DrCorrPlateauUsable`), and the placeholder delivers ε_ΔR up to
  **2.157** at ΔR = 0.355 — where the physics requires ε_ΔR ≤ 1, a close-by 2mu4 correction being a
  LOSS. In the un-merged variant the same cell IS fitted, but its delivered ε_ΔR still reaches 1.354.
- **How much of the collapse is ε_ΔR — measured, not asserted.** With a firing fraction of 0.4925
  and C = 0.5213, the implied mean predicted pair probability is **0.945**, while ε_MC·ε_MC ≤
  0.90² = 0.81 there (the parent doc's forward-negative saturation, R8) ⇒ the mean applied ε_ΔR is
  **≳ 1.17**, above the physics bound and contributing. It is **not the dominant term**, though:
  rebuilding the merged cell from the *fitted* un-merged numerators gives C = 0.5605 against the
  placeholder-corrected 0.5341, so the placeholder-vs-fit swap is a **−4.7 %** effect, and the cell
  collapses just as hard (0.5213) in the un-merged variant where it *is* fitted. The ε_MC map
  carries roughly two thirds of the suppression (≈1.6× against ≈1.17×).
  So the round-1 attribution (the forward-negative ε_MC map, parent R8/R10/R14) **remains the
  leading term but is not the whole story**: a delivered ε_ΔR > 1 in that very cell contributes, and
  it is the cell the cross-section corrects with a placeholder. R1b therefore touches Remaining
  Work 0 / parent R26 as well as the ε-map thread.

**This is the cell the pp24 cross-section is currently applying that placeholder to** — the same
cell its own thread flagged (`pp24_crossx_rerun_2026_08.md`: "one forward ε_ΔR cell delivers > 1").

### R6. Guards, coverage and the placeholder census (round 2, Tight)

- Cell census, 8-bin `nocorr`: `expo` 68 fitted / 2 raw-bin placeholder / 2 with no correction;
  `polyu_fixedRp` 65 / 4 / 3. Merged `nocorr_ptmerge`: the cascade routes **60 cells to the
  exponential, 2 to the polynomial, 1 to the raw bins** — identical to what the pp24 crossx
  pipeline reported, which is the cross-check the merged variant exists for. The census now travels
  inside each output file's `provenance` string, not only in the log.
- Placeholder reach per EVALUATION: 0.075 % (`nocorr`/`expo`), 0.020 % (`nocorr`/`polyu`),
  **0.006 % (cascade, from its own per-pair census)**. The 0.013 % its sub-evaluator printed was an
  artefact: the load-time scan that bounds the delivered ε_ΔR landed in the same counters (99
  raw-bin + 101 empty-bin diagnostic calls on top of a real 71). **Fixed 2026-08-18** — the scan now
  snapshots and restores the sub-evaluators' counters, and the cascade prints its own per-pair
  share. ~22 % of evaluations are at ΔR ≥ 1, where no correction is applied by
  construction; 78 % come from a fit.
- 0 pairs outside the cell grid, 0 ε_ΔR floorings. `p = ε_MC·ε_MC·ε_ΔR > 1` in 49 (`nocorr`/`expo`)
  and 70 (cascade) of 1.34 M triggered pairs — counted and printed, never clamped, because clamping
  `p` would silently alter the weight as well as the Bernoulli error term.
- 66.63 % of selected opposite-sign pairs fall inside the correction cells (the rest are below
  8 GeV pair p_T, where no correction is defined) — unchanged from R3.

## Remaining Work

0. **★ NEEDS A USER DECISION — a sanity bound on the DELIVERED ε_ΔR** (see R2). The baseline screen
   catches a catastrophic `C = 0.053` but not a marginal `C = 0.553 ± 0.339` that still yields
   ε_ΔR = 2.74 where the physics says the correction must be ≤ 1. Candidate screens: (a) reject any
   cell whose ε_ΔR exceeds 1 anywhere on ΔR < 1; (b) reject on the baseline's relative error
   σ_C/C; (c) leave it and let the R26 resolution fix it. Coupled to R26, hence not chosen here.
1. The §3.3 raw-binned fallback is a **temporary placeholder** — replace once parent doc R26 is
   decided by the user.
2. ~~Re-run on the post-gap-cut-rerun inputs.~~ **DONE 2026-08-13** — see the note at the head of
   Results. The `_pt4bin` variant remains outside the parent doc's rerun by design (its R24b), but
   this thread regenerated the 4-bin Step-3 fill, plateaus and opposite-sign no-plateau-correction
   fits itself (02:59–03:05), all on top of the final 02:50 MC turn-ons, so the 4-bin sub-version
   is internally consistent too.
3. **Migrate `FillMCTrigEffHists.cxx` onto `Utilities/MCTrigEffPairSelection.h` and
   `Utilities/SingleMuEffEvaluator.h`** once the sibling session releases that file (D3b). Until
   then two byte-identical copies of the Step-3 pair selection and of the turn-on lookup coexist,
   and nothing enforces that they stay identical. Cheap interim guard worth adding at the same
   time: have `FillMCTrigEffHists.cxx` assert its locally built `sel` equals
   `MCTrigEffPairSel::Step3PairSelection()`.

## Latest Stage

**2026-08-18 — ROUND 2 DONE.** Applied ε = **ε_MC** (self-contained closure, D4); rerun on the
current inputs; **two ε_ΔR variants, one subdirectory each**. `/review-analysis-code` PASS after 3
iterations (all findings were doc/guard consistency, none moved a number; every value independently
re-read from the ROOT files). `/review-plot` run on the figures.

Deliverables — 12 PNGs, all under
`/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/`:

| pair-pT binning | WP | directory | variant subdirectories |
|---|---|---|---|
| 8 (canonical) | Tight | `mc_based/closure/` | `no_plateau_correction/` (3 lines) + `no_plateau_correction_last2ptbins_merged/` (2 lines) |
| 8 (canonical) | Medium | `mc_based_medium/closure/` | both |
| 4 (variant) | Tight | `mc_based_pt4bin/closure/` | `no_plateau_correction/` only |
| 4 (variant) | Medium | `mc_based_pt4bin_medium/closure/` | `no_plateau_correction/` only |

Each subdirectory holds `closure_pair_pt_all_opposite_sign.png` and
`closure_pair_pt_single_b_signal_cuts.png`. Histograms:
`<mc dir>/mc_trig_eff_closure_pp24_full{_medium_wp}{_pt4bin}{_nocorr|_nocorr_ptmerge}.root`.
Driver: `pipelines/run_mc_trigeff_closure.sh` (`MODES`, `WPS`, `PTBINS`, `REGEN_4BIN`, `SKIP_FILL`).

**Headline: C(ε_MC) = 0.994 (all OS) / 0.989 (signal region), Tight; the merged cascade variant
agrees with the un-merged exponential one to 1×10⁻⁵.** The correction machinery closes.

**Open for the user:**
1. **(R5, was R1b) The forward-negative top cell.** η^pair ∈ [−2.4,−2.0) × p_T^pair > 72 GeV closes
   at **0.52–0.53** against a mirror bin at 1.08–1.13, and in the merged variant it is the ONE cell
   the cross-section's cascade corrects with the **raw-bin placeholder**, delivering ε_ΔR up to
   **2.157** where the physics requires ≤ 1. Round 1 attributed the collapse to the forward-negative
   ε_MC map; round 2 shows a delivered ε_ΔR > 1 in that same cell contributes directly. Coupled to
   Remaining Work 0 and parent-doc R26.
2. **(R4 caveat) Per-cell vs inclusive.** "Closes to ~1 %" is inclusive; individual cells sit at
   16–18σ from 1. A trigger systematic built from this must be per cell.
3. **Quarantined round-1 outputs** (8 PNGs, 4 ROOT files, 8 logs, all under
   `superseded_round1_data_eps/`) are kept, not deleted — say the word and they go.
