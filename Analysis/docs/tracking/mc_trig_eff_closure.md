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

## Autonomy Contract (ACTIVE — re-read on every compaction)

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
     result here was produced at 04:16–04:19, after all of its inputs, and the driver's freshness
     gate enforces it.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The analysis corrects each surviving muon pair by `1/ε_trig^pair`. For pp24 / 2mu4 that is

```
ε_trig^pair(pair) = ε^nc(p_T1, q·η1) · ε^nc(p_T2, q·η2) · ε_ΔR^2mu4(ΔR)
```

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
(σ_slice·ε_filt·r_isospin/N_slice). The three plotted lines are the denominator (line 1) and the
numerator for two ΔR-fit forms (lines 2, 3), all as dσ/dp_T^pair; the ratio panel is `C`.

**Which ε goes in the weight (parent doc RW-8 design choice (a), settled by the user 2026-08-13):
the DATA-derived ε^nc**, not ε_MC. That is the configuration the analysis actually applies, so it
is the physics test rather than a code check.

**The ε_MC variant is nevertheless measured, as a DIAGNOSTIC (not a plotted series).** It is the
self-consistency check: ε_MC is the very efficiency the Step-3 inverse weighting divided by, so
that variant closes at 1 if and only if the ΔR correction and its fit are self-consistent. The two
together separate what a single closure number confounds:

```
C(ε_MC)                 → does the ΔR correction + its fit close?
C(ε_data) / C(ε_MC)     → ⟨r₁·r₂⟩,  r_j = ε_MC(leg j)/ε_data(leg j),
                          weighted by w_MC/(ε_MC,1·ε_MC,2·ε_ΔR)
```

**⚠ It is the mean of the PRODUCT of the two legs' ratios, NOT ⟨r²⟩.** Reading a per-leg number as
`√⟨r₁r₂⟩` is valid only where both legs sample the same q·η regime. For an **opposite-sign** pair
the legs sit in **mirrored** q·η bins (q·η₁ ≈ −q·η₂), and the parent doc's R8/R10/R14 record a
strong q·η asymmetry (the forward-**negative** MC turn-on is saturated, its mirror is not). So the
square-root reading holds in the barrel panels and **not** in the endcap ones.

Without it a reader cannot tell a broken ΔR correction from the MC trigger simulation simply being
more efficient than the detector — and R1 shows that is exactly which of the two the data-ε
closure is dominated by. The diagnostic is written to the output ROOT file and printed per pair-η
bin; the FIGURE carries only the three series the user specified.

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
  `minv > 1.08 && minv < 2.9 && pair_pt > 8 && m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2`.
  The gap cut is already in the base selection; in the *data* analysis it is not yet applied to
  the signal selection (a future to-do owned by `muon_gap_cuts_acceptance.md`, explicitly out of
  scope here).

**Two SUB-VERSIONS** per version: the ΔR correction taken from the canonical **8**-bin pair-p_T
cells and from the **4**-bin variant (`MCTRIGEFF_PAIRPT_4BIN=1`). The plotted pair-p_T axis is the
SAME binning as the correction's cells in each case — a 1D view of the cells the correction was
measured in (CLAUDE.md §Binnings item 2).

#### §3.2 The per-pair trigger weight

For a pair with legs (p_T1, q·η1), (p_T2, q·η2), separation ΔR, in cell (i_ptpair, i_ηpair):

```
ε_trig^pair = ε^nc_data(p_T1, q·η1) · ε^nc_data(p_T2, q·η2) · ε_ΔR(ΔR ; cell)
```

- **ε^nc_data** — the pp24 tag-and-probe turn-on TF1s, looked up per (charge, coarse q·η bin) and
  **evaluated continuously at the exact p_T** (`TF1::Eval`, never resample-to-nearest), with the
  p_T clamped into the fit range, capped at 1 and floored at 0.02. This is the identical lookup
  the analysis and `FillMCTrigEffHists.cxx`'s `DataEffEvaluator` use.
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

One PNG per (version, pair-p_T binning, working point). **One subplot per pair-η bin** of
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` (9 bins → 3×3). Each subplot:

- **upper pad** — dσ/dp_T^pair (log–log), three series: (1) all pairs, no trigger requirement;
  (2) pairs passing 2mu4, corrected with the `expo` ΔR fit; (3) same, `polyu` ΔR fit.
- **lower pad** — the closure ratio `C` (series 2 and 3 over series 1), with a line at 1.

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
assumes the closure already holds — the very thing being tested, and it misses by ~27 %. Since
`B ≤ A`, dropping `R` shifts the variance by `(R−1)·B`, so it **over**states σ where `R > 1`
(measured: up to 1.32× at `R = 1.52`) and **under**states it where `R < 1` (up to 2.46× at
`R = 0.82`) — the second direction is the dangerous one and it is not rare, since `R < 1` in the
statistics-poor high-pair-pT cells. It is the *correct* form that can reach the binomial `k = n`
boundary (`Var → 0`; 1 cell of 72 here, `R = 1.098`), and a naive implementation would set that
error to 0 and both consumers would drop the point. The computation is therefore delegated to
`SetConditionalRatioErrors` in `dr_correction_ratio.h` — the repo's **single** implementation,
which carries that boundary fallback (`max(R,1)/n_eff`).

### 4. Negative constraints

- **NO trigger requirement anywhere in the denominator.** The denominator is every selected pair.
- **Do NOT weight with ε_MC.** That is the self-consistency check, not the closure the analysis
  needs (§2). ε_MC appears in this measurement only where it already lived: inside the Step-3
  inverse weighting that produced ε_ΔR.
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

## Results & Observations

*Measured on the **FINAL post-gap-cut inputs**. The sibling session's round-11 rerun completed at
03:2x on 2026-08-13 and its hand-off note (`mc_trigger_efficiency.md` Latest Stage / R29) declares
those products final; the numbers below come from closure fills at **04:16–04:19**, i.e. after
every one of them. Verified input mtimes: MC Step-3 hists 02:51 (8-bin) / 02:59–03:03 (4-bin),
MC single-muon turn-ons 02:50, data tag-and-probe turn-ons 02:17 (Tight) / 02:39 (Medium),
opposite-sign no-plateau-correction ΔR fits 03:05–03:18. The driver enforces this: `val_fresh`
refuses a closure output older than any input it consumed.*

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

*(Source: `pipelines/logs_closure/fill_{8,4}bin_{tight,medium}.log` and
`mc_trig_eff_closure_pp24_full{_medium_wp}{_pt4bin}.root`, run 2026-08-13 03:35–03:39. An earlier
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
**0.737**. Rejections by this screen, per (method, binning): `expo` **0** in both 8-bin runs, **1** in each
4-bin run; `polyu_fixedRp` **1** in every run. Every rejection is printed by name and the count
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

**2026-08-13 — ALL Implementation Plan steps DONE, on the FINAL post-gap-cut inputs. Both reviews
run; two findings escalated to the user (below).**

Delivered (all under `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/`):

| pair-pT binning | WP | directory | PNGs |
|---|---|---|---|
| 8 (canonical) | Tight | `mc_based/closure/` | `closure_pair_pt_all_opposite_sign.png`, `closure_pair_pt_single_b_signal_cuts.png` |
| 4 (variant) | Tight | `mc_based_pt4bin/closure/` | same two |
| 8 (canonical) | Medium | `mc_based_medium/closure/` | same two |
| 4 (variant) | Medium | `mc_based_pt4bin_medium/closure/` | same two |

Histograms: `<mc dir>/mc_trig_eff_closure_pp24_full{_medium_wp}{_pt4bin}.root`.
Driver: `pipelines/run_mc_trigeff_closure.sh` (`REGEN_4BIN=1` rebuilds the 4-bin upstream).

**Open for the user (R1b):** the closure exposes a ~2× over-prediction in
`η^pair ∈ [−2.4,−2.0) × p_T^pair ∈ [72,104) GeV` that survives the ε_MC diagnostic and is absent
from the mirror bin — attributed to the forward-negative r16578 ε-map anomaly above the
`p_T > 7` veto. Re-check on the final inputs before acting.

**Open for the user (R1):** the as-applied closure sits at ~1.27, entirely accounted for by the
MC/data single-muon over-efficiency ⟨r₁·r₂⟩ (**not** ⟨r²⟩ — the two legs of an opposite-sign pair
sit in mirrored q·η bins); the ε_MC diagnostic closes at 0.994. Decide (a) whether the eps_MC
closure should also become a plotted series rather than a printed diagnostic, and (b) whether the
residual mis-closure map of R1 becomes a trigger systematic.
