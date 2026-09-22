# Sizing a new high-pT MC request: mass filter and pair statistics in the top two pT-hat slices

**Mode:** Implementation. **Created:** 2026-09-09.
**Session:** "additional fullsim statistics".
**Parents:** `mc_trigeff_single_value_pair_eff.md` (ACTIVE — the single-value pair efficiency this
request is meant to feed), `mc_trigger_efficiency.md` (ACTIVE — Step 3, the ΔR correction),
`mc_trigeff_dr_binning_approaches.md` (ACTIVE — the cell groupings),
`mu_pt45_gap125_pairpt9_adoption.md` (ACTIVE — the selection every input here is measured under),
`ami_weights.md` (the DSID/σ registry the denominator is proved against).

---

## Objective

Provide the numbers needed to **specify a new Pythia fullsim MC request** in the two highest
pT-hat slices (**70–125 GeV** and **125–300 GeV**), with a **loose dimuon mass filter** applied at
generation:

1. Where does the **back-to-back mass peak** sit, per slice? — the mass distribution of the
   selected pairs, so the filter can be placed below it.
2. How many **pairs per event** survive `m < 10 GeV` and the signal window
   `1.08 < m < 2.9 GeV`, per slice and per pair sign?
3. The same, resolved in the **3 highest coarse pair-pT cells × 3 |η^pair| groups** — the cells
   where both the ΔR correction and the single-value pair efficiency are statistics-starved, and
   the reason the request exists at all.

Deliverables are a plot set and CSV tables under a new directory in the pp trigger-efficiency
plot tree.

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The pp24 per-pair 2mu4 trigger weight is factorized as
`ε(p_T1,q·η_1) · ε(p_T2,q·η_2) · ε_ΔR(ΔR ; cell)`. In the top coarse pair-p_T cells the existing
fullsim sample cannot constrain the ΔR shape: fits are rejected, and one cell is corrected by a
raw-bin placeholder at ε_ΔR = 2.157 where the physics requires ≤ 1
(`mc_trigger_efficiency.md` R26/R32, `mc_trig_eff_closure.md` R5). The single-value pair
efficiency (`mc_trigeff_single_value_pair_eff.md`) removes the fit but still needs enough pairs
per cell to give a usable binomial error. **Both procedures are limited by the same thing: pairs
in the top pair-p_T cells.** That is what the new request must buy.

Requesting raw events in those slices is inefficient: most generated events give no low-mass
dimuon at all. A **generator-level dimuon mass filter** raises the useful fraction per event.
The filter has two competing requirements, and the whole point of this study is to show that
both can be met at once:

- **Low enough to remove the back-to-back pairs.** At high pT-hat a large population of pairs
  comes from two muons recoiling against each other from *separate* heavy-flavour quarks. Their
  opening angle is ΔR ≈ π, so their invariant mass is large (m ≈ 2√(p_T1 p_T2) for a
  back-to-back, central pair). They are **different physics** from the single-b signal — a
  different ΔR mixture, a different leg-p_T sharing, hence potentially a different trigger
  efficiency. Averaging them into the same (p_T^pair, η^pair) cell would bias the measured
  efficiency of the cell towards a population the analysis does not use.
- **Loose enough not to bias the result.** The filter must sit well above the signal window
  (1.08–2.9 GeV) so that the measured efficiency is not conditioned on a cut the analysis itself
  does not apply, and it must leave enough of the ΔR range populated for the **ΔR correction**
  procedure to remain measurable — not only the single-value procedure. `m < 10 GeV` is the
  proposed value; this doc measures whether it satisfies both.

### 2. Top-level quantities

Per pT-hat slice `s`, per pair sign, per mass range `R`, and optionally per
(pair-p_T, |η^pair|) cell:

```
N(R)  = number of selected pairs with m in R          [raw, UNWEIGHTED]
r(R)  = N(R) / N_events(s)                            [pairs per generated event]
```

`r` is the deliverable: the expected pair count of a **new** request of `N` events in that slice
is `r · N`.

**⚠ WHAT `N` COUNTS — an open question that must be settled before the request goes out.** `r` is
measured on a sample with **no generator-level mass filter**, so `r · N` is the yield from `N`
events generated **before** any such filter. If the new request's event count is quoted **after**
the mass filter — the usual convention, and what AMI `totalEvents` records — the yield is larger,
by `1 / ε_massfilter`. **ε_massfilter is a generator-level number and is NOT measured anywhere in
this doc**, so as it stands this study sizes the *unfiltered* request exactly and the *filtered*
one only up to that factor. Raised by `/review-analysis-code`; see Remaining Work.

**Why raw and unweighted.** This is a sample-size question, not a cross-section question. Within
ONE pT-hat slice the per-pair MC weight is a single constant `σ_s · ε_filt,s / N_s`, so the
weighted and unweighted numbers carry identical information, and only the raw one answers
"how many pairs will N events give me". These are explicitly **not** yields: slices are combined
with very different weights, and no number in this doc may be summed across slices.

**Why the denominator must be right, and where it comes from.** Every percentage here is
`N / N_events`, so a wrong `N_events` silently rescales the entire request. It is **read from the
pair file**, not inferred: `meta_tree_out` carries `nproc_kin<K>_beam<B>`, the number of events the
ntuple processing **actually looped over** for that slice, alongside `nbeam_kin<K>_beam<B>` (the
events available). A pair file with no usable meta tree is **refused**
(see R0 for why that guard exists).

Three cross-checks on the number read from the file, all required to pass:

- **(a)** `nbeam_kin<K>_beam0` equals the entries of the NTUP chain on disk — catches a pair file
  and an NTUP farm that are out of step, or a partially symlinked farm;
- **(b)** `σ·ε_filt / w = N_proc`, where `w` is the constant per-pair weight of the slice and
  `σ·ε_filt` comes from that slice's own **DSID-guarded** AMI file (`ami_weights.md` — the AMI
  files are named by beam + slice only, so another production's file opens silently). Since the
  2026-09-09 upstream fix the weight divides by `N_proc`, so this is exact;
- **(c)** the AMI production record's `totalEvents`, which lives off this machine entirely.
  Tolerance 2 % (= 6 400 of 320 000). The NTUP farm is **3 part-files of ~107 k events** each
  (AMI's `nFiles: 32` counts AODs, not NTUP parts), so a dropped part is a 33 % deficit and throws,
  while the one event genuinely absent from each `_pdf` slice is 0.0003 % and passes.

A truncated run is now **reported**, not hidden: if `N_proc < N_beam` the macro says so loudly and
notes that the per-event rates remain correct (they divide by `N_proc`) while the absolute counts
are those of a partial run.

**One thing upstream does NOT persist, checked rather than assumed.** `PythiaAlgCoreT` has an
in-memory `meta_fullsim_truncated` flag (`PythiaAlgCoreT.h:141`) but **never `Branch()`es it**, so
it does not reach the pair file: `meta_tree_out` carries only the `nentries_`, `nproc_` and
`nbeam_` families (54 branches, enumerated). **`nproc != nbeam` is therefore the only file-level
truncation detector**, and it is the one this macro uses.

### 3. Step-by-step method

**§3a The pair sample.** `MCTrigEffPairSel::Step3PairSelection(true)`, read from the per-kn pair
trees `muon_pair_tree_kin{4,5}_sign{1,2}` that the ntuple processing already writes
(`PythiaAlgCoreT::fill_kn_trees_fullsim`). This is **byte-identical to the population the ΔR
correction (§3.3 of `mc_trigger_efficiency.md`) and the single-value pair efficiency are measured
on**, which is what makes the counts a statement about the statistics *those two procedures* would
gain. Explicitly: both legs' pair-level WP flag (quality bits + one-sided Δp/p < 0.12 + d0/z0 +
the pp same-vertex pair requirement), p_T > 4.5 GeV, |η| < 2.4, the truth fiducial
(truth p_T > 4.5, |truth η| < 2.4), the forward low-p_T veto (p_T > 7 || q·η > −2), the per-leg
fiducial gap windows and the pair-level |η^pair| < 2.2.

**§3b What is NOT required.** *No trigger requirement* — these are the Step-3 denominator pairs,
and the request is about how many pairs exist to measure an efficiency *from*, not how many fire.
*No signal-region pair-p_T cut* and *no resonance veto* (the MC pair file carries none).

**§3c The mass ranges.** `m < 10 GeV` (the proposed filter, the one number this study introduces)
and the signal window `1.08 < m < 2.9 GeV`, taken from `PairTrigEff::Window("sig")` and guarded
against the signal region itself by `PairTrigEff::CheckSignalWindowMirror()`. **The two ranges
overlap**: the signal window is a subset of `m < 10`, and the tables must be read that way.

**§3d The cells.** Pair p_T = the 3 highest cells of `ParamsSet::pair_pt_coarse_bins` (the
canonical 8 log bins 9→150 GeV, so the 3 highest are **[52.23, 74.24), [74.24, 105.53),
[105.53, 150)** GeV), located by `PairTrigEff::FirstDeliveredPtEdge()` so the selection of "the
3 highest" follows a binning change instead of being retyped. |η^pair| = the 3 sign-independent
groups `PairTrigEff::AbsEtaGroups()` builds from
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` — **|η| < 1.0 / 1.0–2.0 / 2.0–2.2**.
Both are the same axes the single-value pair efficiency is measured on, by construction.

**§3e The mass figure.** Raw pair counts on an axis uniform in log(m), so the drawn height is
dN/dln(m) up to one constant and the position of a peak is not a binning artefact. Drawn twice:
over all selected pairs, and restricted to the 3 highest pair-p_T cells (the region the request
is about). Both signs on the same axes, with the proposed 10 GeV filter and the signal window
marked.

### 4. Negative constraints

- **These are not yields and not cross sections.** Nothing here may be summed across pT-hat
  slices, and no number here belongs in a σ or R_AA. The per-slice weights differ by a
  factor of **7.4** (`ami_weights.md` table B: 1.4588 nb vs 0.1978 nb).
- **The counts do not include a trigger requirement.** They are the denominator population. A
  reader wanting "how many pairs will fire 2mu4" must multiply by the pair efficiency, which is
  what `mc_trigeff_single_value_pair_eff.md` measures — it is not in this doc.
- **The percentages are per-event rates, not efficiencies.** A cell can exceed 100 % in principle
  (more than one selected pair per event); it is a rate, not a probability.
- **This is not the muon-p_T-4.5 acceptance study.** The −23.25 % cost of the 4.0 → 4.5 move is
  recorded in `muon_pt45_cut_diagnostic.md` / `mu_pt45_gap125_pairpt9_adoption.md` and is
  *already spent* here: the Step-3 selection applies p_T > 4.5 literally, so this population is
  the post-cut one.

---

## Context

- User request (2026-09-09): request additional high-p_T statistics in the two highest pT-hat
  slices, with a loose mass cut (e.g. m < 10 GeV), one sign-independent sample and one filtering
  on same-sign pairs. Deliverables as listed in Objective.
- Reference sample size in the request: **350 k events** per slice. The existing slices hold
  **319 999** events each, so the counts are reported on the existing sample with the percentage
  as the **per-event rate** (user decision, 2026-09-09) — multiply by 350 000 for the projection.
- **User request (2026-09-17) — a SAME-SIGN-ONLY request, wider in pT-hat and lower in
  pair p_T.** The same-sign MC statistics are poor down to pair p_T ≈ 25 GeV, so the second
  request is now: same-sign-filtered, a wider range of pT-hat slices, a much lower pair-p_T cut
  (~20–25 GeV) and a loose mass cut whose value may need re-deciding (the back-to-back peak
  moves DOWN in mass at lower pair p_T). Preparation asked for: (1) the inclusive mass figure for
  the top FOUR slices (add 24–40 and 40–70); (2) the same with a pair-p_T cut of 20–150 GeV in
  place of the 52.23–150 top-3 cells; (3) the fullsim per-slice statistics plots
  (`plot_pythia_fullsim_kn_pt_crossx.cxx`, existing statistics, NOT projected) for same-sign
  pairs with m < 10 GeV in place of the single-b signal pairs. Nothing existing is overwritten.
- The three highest pair-p_T bins used to be 50–72 / 72–104 / 104–150 GeV. Since the
  2026-09-08 adoption of the 9 GeV pair-p_T floor
  (`mu_pt45_gap125_pairpt9_adoption.md`), the coarse axis is 8 log bins **9→150** and the three
  highest are **52.23–74.24 / 74.24–105.53 / 105.53–150** GeV. Read from `ParamsSet`, never
  retyped.

## Scope

**In:** the two highest pT-hat slices of the pp24-conditions Pythia fullsim FULL sample; the
Step-3 trigger-efficiency pair population; both muon working points (Tight nominal, Medium
alongside); the mass figure and the nine CSV tables per WP (2 mass-window counts, 2 mass compositions,
4 cell tables, 1 spectrum-landmarks file).

**Out:** Pb+Pb / the HIJING overlay (the per-pair weight there is the mu4 union, a different
formula); the trigger decision and any efficiency; the actual text of the production request;
any change to the ntuple processing or to a canonical binning.

## Design Decisions

**D1 — A standalone macro, not an extension of `FillMCTrigEffPairEff.cxx`.**
The single-value fill is a *measurement* producer wired into the closure comparison; this is a
*sizing* study with a different denominator (events, not pairs), a different weighting (raw, not
σ-weighted) and a per-pT-hat-slice split the measurement deliberately integrates over. Bending it
to carry both would put a sample-size question inside a physics-result file. The macro reuses the
measurement's selection and cell headers verbatim (`MCTrigEffPairSel`, `PairTrigEff`), so the two
cannot describe different populations.
*Rejected alternative:* extending the `book_kn_stats` block of `FillMCTrigEffHists.cxx` (added
2026-09-03 for the earlier per-slice statistics request). It already splits by kn and by sign
under the same Step-3 selection, but has **no mass axis at all**, and adding one would have put
mass-windowed statistics bookkeeping into the ΔR-correction fill, whose output file every
downstream fit stage opens.

**D2 — Raw, unweighted counts; the percentage is the per-event rate.** See §2. Physics reason:
within one pT-hat slice the weight is one constant, so weighting adds nothing, while the raw
count is the only form that answers the request.

**D3 — The denominator is READ from the file, then cross-checked three ways.** See §2.
**REVISED TWICE.** (i) The original decision claimed "proved twice … the two independent routes
fail in different ways"; that was wrong — routes (a) and (b) both reduced to `N_beam`, because the
weight was itself built from `N_beam`. (ii) A third route (AMI `totalEvents`) was added, but none
of the three could see `nevents_max` truncation, which was then documented as a residual risk and
escalated. That risk is now **closed**: the upstream fix records `N_proc` in `meta_tree_out` and
divides the weight by it, and this macro reads `N_proc` as the authoritative denominator and
refuses a pair file that lacks it (R0). Truncation is detected as `nproc != nbeam`; the upstream
`meta_fullsim_truncated` flag is not persisted to the file (R0). Physics reason unchanged throughout: a wrong `N_events`
rescales the whole request.

**D4 — A local log mass axis for the figure only.** No canonical minv binning spans the required
range: `hist_binning_map["minv_log_bins_{ss,op}"]` stops at 60 GeV (below the back-to-back region
the figure exists to show) and is sign-dependent, so SS and OS could not share an axis;
`ParamsSet::hq_minvBins` starts at 2.41 GeV (above the signal window the same figure must show).
The axis is **80 bins uniform in log(m), 0.211 → 400 GeV**. **REVISED 2026-09-09** from 70 bins /
0.2 → 220 GeV: 220 GeV truncated the spectrum (the highest selected pair mass is **358.96 GeV**,
and 30 pairs — 17 OS + 9 SS in pTH125_300, 2 + 2 in pTH70_125 — sat in overflow, silently absent
from the drawn curve); the low edge moved 0.2 → 0.211 so that bin 1 =
[0.21100, 0.23188) is only **1.66 %** dead below the 2m_μ = 0.21133 GeV kinematic threshold
instead of 58 %. It **still straddles** that threshold — no log axis edge lands on it and none
can — it simply no longer reads as a dip. Overflow and underflow are now verified
**0** in all 16 histograms at both working points.

**No table may depend on this axis** — every counted number comes from an exact
`minv > lo && minv < hi` filter. **This rule was VIOLATED once and the violation is why the
`mass_composition_top3_*.csv` tables exist** (D7): R2's first fractions table was integrated off
display-bin edges instead. **No log display axis has an edge at 2.9, 4 or 10** — on the axis in
use when that error was made they fell at 2.8602 / 4.1822 / 9.8358, and they move again every time
the axis does — so the bands were split at the wrong masses, and the table was wrong by
`/review-analysis-code` independently, 2026-09-09.

**D7 — The mass-composition table is booked as exact counts, not read off the figure.** Four
exclusive, exhaustive bands of the top-3 pair-pT population (`< 2.9`, `2.9–4`, `4–10`, `≥ 10`),
each an RDF `Count()` on an exact `minv` filter, plus a J/ψ-window sub-count. Every boundary is
read from a canonical source — `PairTrigEff::Window("sig").hi`, `Window("wide").hi`,
`kLooseMassMax`, and `ParamsSet::minv_cuts` for the J/ψ window — and a guard throws if they ever
stop being ordered. Physics reason: this is the table the *filter decision* is argued from, and
the failure it replaces was a silent 6-point error in exactly that number.

**D5 — Log y on the mass figure.** Permitted because the x axis is log-binned and the spectrum
falls by several decades; on a linear y the low-mass continuum would compress the back-to-back
region the figure exists to show into a flat line. Documented in the code, per the log-scale rule.

**D6 — An input-stability guard.** The pair file is regenerated by the ntuple processing, which
another session runs on the same checkout. The macro stats the file before and after the event
loop and **throws before writing anything** if size or mtime moved: a mid-read regeneration would
otherwise leave every histogram filled with fewer pairs and no error from RDF — indistinguishable
from the honest answer to the question being asked.

**D8 — The 4-slice / 20 GeV figures are a FIGURE-ONLY variant of the same macro, never a second
copy of the selection.** `mc_pthat_slice_mass_statistics` gains two arguments,
`n_top_slices` (2 = the nominal request study; 4 adds pTH24_40 = DSID 803017 and pTH40_70 =
803018, both from `ami_weights.md` table B and both DSID-guarded like the first two) and
`figure_pair_pt_lo` (NaN = the canonical 3-highest-cells edge from `ParamsSet`; 20 = the user's
figure cut). Any non-nominal value switches the macro to figure-only: the nine CSV tables are
shaped for the 3 highest coarse cells of the FIRST request and are NOT re-emitted, and every
output name carries the variant suffix (`_top4`, `_pairpt20`), so the 2026-09-10 deliverables
are untouched. **20 GeV is a cut on the figure, not a binning**: it is not an edge of
`pair_pt_coarse_bins` (nearest edges 18.18 / 25.84) and is labelled as the range it is.
Same population, same axis, same denominator checks as D1–D6. The proposed filter value is
moved to ONE header (`Utilities/MCRequestLooseMassFilter.h`) so that the fullsim statistics
macro (D9) and this one cannot drift apart on it.

**D9 — The same-sign statistics plots are a selection VARIANT of `plot_pythia_fullsim_kn_pt_crossx.cxx`,
writing to its own subdirectory.** Nominal = opposite-sign tree, `from_same_b` (the single-b
signal). Variant `ss_loose_mass` = same-sign tree (`muon_pair_tree_kin<k>_sign1`, split by
`truth_same_sign` upstream), truth pair p_T under `truth_minv < 10` and reco pair p_T under
`minv < 10 && pair_pass_<WP>` — each variable cut at its own level, because the generator filter
of the request is a truth-mass cut while the analysis consumes reco mass. Weighted (dσ/dp_T in
nb/GeV, as the nominal plots) — here the slices ARE compared across each other, unlike D2, so
the σ-weighting is required. The projected-statistics function is guarded against the variant
(it hard-codes the signal selection). Outputs → `plots/ss_mass_lt10/`, same file names.

## Implementation Plan

1. **DONE** — Doc triage; establish the canonical cells and the slice/DSID map; resolve the
   350 k-vs-320 k denominator, the selection and the WP with the user.
2. **DONE** — `plotting_codes/trig_effcy/mc_based/mc_pthat_slice_mass_statistics.cxx`
   (implements §3a–§3e). Reviewer: `/review-analysis-code` for the RDF/selection half and
   `/review-plot` for the figure; include §2, §3a–§3e and §4 in the task prompts.
3. **DONE** — run for Tight and Medium; numbers in Results & Observations (R1–R6).
4. **DONE** — `/review-analysis-code` (RDF/selection/denominator, §2/§3a–§3d/§4) and
   `/review-plot` (the mass figure, §3e/D4/D5); committed.
5. **DONE 2026-09-17 (no reviewer pass — user: urgent, reviews explicitly waived)** — same-sign-request preparation, three figure sets (§3e, D8):
   (a) `pair_mass_by_pthat_slice_top4.png` — 4 slices (kin 2–5), all selected pairs;
   (b) `pair_mass_by_pthat_slice_top4_pairpt20.png` — 4 slices, 20 < p_T^pair < 150 GeV;
   (c) `<full sample>/plots/ss_mass_lt10/{truth,reco}_pair_pt_kn{,_stat_error,_err_frac}{,_150GeV}.png`
   — the per-slice dσ/dp_T statistics plots for SAME-SIGN pairs with m < 10 GeV.
   Reviewers: `/review-analysis-code` (selection, slice/DSID map, denominator guards on the two
   added slices) + `/review-plot` (§3e, D4, D5, D8) — one combined pass on both macros.

## Progress Log

- **2026-09-09, step 1.** Doc triage: read `INDEX.md`,
  `mc_trigeff_single_value_pair_eff.md`, `ami_weights.md`,
  `pythia_fullsim_pp24_full_sample_skim.md`, `Utilities/MCTrigEffPairSelection.h`.
  Established: the two highest slices are `ikin` 4 (pTH70_125, DSID 803019) and 5 (pTH125_300,
  DSID 803015); each holds **319 999** events on disk (32 files × 10 k per `_pdf` dataset, one
  event short), **not** the 350 k of the request; the canonical coarse pair-p_T axis is already
  9→150 GeV so the 3 highest cells are [52.2273, 74.2385), [74.2385, 105.5262), [105.5262, 150);
  the |η^pair| groups are [0,1), [1,2), [2,2.2).
  User decisions: raw counts with the percentage as the per-event rate; Step-3 selection; both WPs.
- **2026-09-09, cross-session.** Session "analysis-a8" (`mu_pt45_gap125_pairpt9_adoption.md`) is
  regenerating this sample. Established and communicated that the Step-3 population on the
  pre-rerun tree is **provably identical** to the post-rerun one: the persisted flag encodes
  p_T > 4.0, the selection ANDs p_T > 4.5 literally, and
  `{quality, p_T>4.0} ∩ {p_T>4.5} = {quality, p_T>4.5}`; the gap windows and every pair-p_T axis
  are read from `ParamsSet` at selection time. The immunity covers the *selection*, not the
  *sample*, so the numbers will be re-derived on the regenerated trees before the request goes out.

- **2026-09-09, step 2 — the macro.** Wrote
  `plotting_codes/trig_effcy/mc_based/mc_pthat_slice_mass_statistics.cxx` (implements §3a–§3e).
  Design D1: standalone rather than an extension of `FillMCTrigEffPairEff.cxx` or of the
  `book_kn_stats` block of `FillMCTrigEffHists.cxx` (the latter already splits by kn and sign
  under the same Step-3 selection but has no mass axis at all). Compiles with ACLiC.
- **2026-09-09, step 3 — the runs.** Tight and Medium on the 2026-09-08 17:11 pair file. Both
  passed the D6 input-stability guard and the R1 denominator check. Outputs in
  `plots/pp_trigger_efficiency/mc_pthat_slice_mass_stats{,_medium}/`. Numbers → R1–R6.
- **2026-09-09, step 4, iteration 1 — `/review-analysis-code` + `/review-plot`, both FAIL.**
  Code: 2 WARNING (legend on data; "denominator proved twice" is not two independent routes) +
  6 INFO. Plot: 3 WARNING (legend; the "clear valley at 8–15 GeV" claim is not in the histograms;
  low-statistics SS drawn as HIST with no errors) + 4 INFO. **The code reviewer independently
  re-derived every reported number with a from-scratch retyped selection — all ~30 MATCH**,
  including both N_slice routes, both σ·ε_filt against `ami_weights.md` table B, all four Table-1
  rows, both 3×3 grids cell by cell, and zero TH2D over/underflow. Sign convention independently
  confirmed by the J/ψ spike appearing only in `sign2`.
  *Amendments:* legend → upper right, cap → ×60, one y range shared by both panels; top-3 figure
  → `E1` with Poisson bars; `mass_top` gained its upper pair-pT bound; mass axis → 0.2–400 GeV,
  80 bins; `ConstantSliceWeight` on both sign trees; third denominator route (AMI `totalEvents`)
  added; CSV headers rewritten to stop claiming "the events actually processed" and to name the
  `nevents_max` residual risk; display ranges reset before `Write()`; `PlainText` advance fixed.
  R2's "valley" claim verified against the histograms by the executor and **withdrawn** (R0/R2).
  Upstream `meta_tree_out` fix raised with the session owning `PythiaAlgCoreT`.
- **2026-09-09, step 4, iteration 2 — both reviewers FAIL again, and both found the same defect
  independently.** R2's replacement fractions table had been integrated off the **display
  histogram's** bin edges (2.86 / 4.16 / 9.84 instead of 2.9 / 4 / 10), wrong by up to **6.2
  percentage points** — the headline same-sign number read 48.4 % where the truth is 54.55 %.
  This violated D4's own "no table depends on the display axis". Both reviewers' exact-cut values
  agreed with each other.
  *Amendments:* **D7** — four exclusive exact-cut mass bands plus a J/ψ sub-count booked in the
  macro and emitted as `mass_composition_top3_<slice>.csv`, every boundary read from
  `PairTrigEff::Windows()`, `kLooseMassMax` and `ParamsSet::minv_cuts`, with an ordering guard;
  R2's table and conclusion rebuilt from that CSV; the denominator is now stated once per
  sentence (it had switched mid-paragraph); the top-3 canvas label now carries **both** canonical
  pair-pT edges (it promised `> 52.2 GeV` for a histogram bounded at 150); D3 and D4 rewritten
  (both still asserted claims that had been withdrawn); the pTH70_125 top-3 peak/minimum row
  re-labelled as unresolvable Poisson noise (~15 ± 4 per bin) rather than a measurement; per-figure
  headroom; axis low edge → 0.211 so bin 1's dead fraction below 2m_μ falls from 58 % to 1.66 % (it still straddles); `totalEvents `
  key given the trailing space its neighbours use; a retyped `1.08-2.9` in a CSV comment replaced
  by `ranges[1].csv_text`; the NTUP-farm granularity statement corrected (3 parts × ~107 k, 2 % =
  6 400 events, not "32 × 10 k").

- **2026-09-09, step 4, iteration 3 — both reviewers FAIL again, same class of defect a THIRD
  time.** R2's peak/minimum table had been measured on the *previous* display axis (0.2 → 400) and
  was never re-derived when the low edge moved to 0.211 in the iteration-2 batch; every bin range
  it quoted no longer existed in the shipped histograms. Additionally: the code still said the
  denominator was "established twice"; D4 attributed the old axis's bin edges to the current axis
  and overstated the 2m_μ fix ("no longer straddles" — it does, just 1.6 % dead instead of 58 %);
  the peak-moves-up-with-pT-hat claim was a **tie** on the current axis; and the pTH125_300 top-3
  "minimum bin" was itself a 2.5 σ dip on a plateau — the same single-bin error the doc had already
  corrected for the other slice. Every count, rate, band, grid and denominator route was
  independently re-derived for the third time and **matched exactly**; only derived prose was wrong.
  *Amendments:* **the landmarks are now EMITTED** — `mass_spectrum_landmarks.csv`, computed from
  each run's own histograms with derived search regions, one row per slice × sign × scope, with a
  standing rule in R2 that no spectrum landmark may be hand-transcribed again. Two statistical
  fixes came with it: the single-bin peak-vs-min comparison is look-elsewhere biased (the extremes
  of a flat distribution are far apart by construction), so it is emitted only under a column named
  `..._LOOK_ELSEWHERE_BIASED` and is not quoted; and the verdict is now **`drop_sigma`**, the
  significance of the step in mean density across the filter over two regions fixed before looking.
  Also added: a checked assertion that the mass bands really are exhaustive (they sum to the top-3
  total); the sign→per-kn-tree index now derived from `PairTrigEff::PairSign::tree` instead of a
  token ternary that mapped any unexpected token to sign2; (a)/(b)/(c) denominator labels made
  consistent between code, doc and CSV headers; D4's two misstatements fixed.
  **This produced the study's sharpest result (R2):** the density step across 10 GeV is +4.5 σ /
  +2.6 σ for OS but **0.0 σ / −4.7 σ for SS** — for same-sign pairs at high pair p_T the filter
  cuts the middle of a flat continuum.

- **2026-09-10, step 5 — rerun on the REGENERATED 4.5 GeV trees, and two defects fixed.**
  Input `muon_pairs_..._mc_trig_full.root` mtime **2026-09-10 04:34:21**, 3 107 190 105 B (the
  file SHRANK from 4 636 631 852 B — the tighter muon cut writes fewer pairs).
  **Every Table-1, 3×3 and composition number is BYTE-IDENTICAL to the pre-rerun run.** That is
  the predicted result and it closes the open item: the Step-3 population is invariant under the
  muon-p_T change because `{quality, p_T>4.0} ∩ {p_T>4.5} = {quality, p_T>4.5}`, and the gap
  windows and pair-p_T axes were always read from `ParamsSet` at selection time.
  *Defect 1 (CRITICAL, `/review-analysis-code` + the peer session independently):* `drop_sigma`'s
  "above" window ended at the DATA's own peak bin while claiming to be "fixed before looking" —
  biased positive by construction, with the compared mass range varying row to row. Fixed to a
  pre-registered log-symmetric window derived from canonical mass scales ([4,10) vs [10,25) GeV,
  R = filter/template-fit-top = 2.5), with `drop_sigma_wide` emitted beside it as a stated
  window-dependence check. **The headline "0.0 σ — perfectly flat" for pTH125_300 SS was an
  artefact and is withdrawn**; the corrected values are +5.2/+2.9 σ (OS) vs +1.8/−1.6 σ (SS).
  *Defect 2:* the macro's `N_proc` narrative was stale after the upstream fix — rewritten, and
  `N_proc` is now read from `meta_tree_out` as the authoritative denominator (R0, now CLOSED).
  Also fixed from the iteration-4 review (**but see the step-6 entry: two of these were logged
  as fixed while the edit that would have made them had aborted, and were only actually applied in
  step 6**): D4's false "no longer straddles 2m_μ" claim; the
  unsupported "minimum lies above 10 GeV wherever resolvable" headline (deleted — it rested on a
  look-elsewhere-biased single bin); "within ~12 %" (OS only; SS differs by 40–58 %); the Medium
  range (1.040–1.076, not 1.04–1.07); "more than an order of magnitude" → ×7.4; the retired axis's
  4.164 → 4.1822; and the ⚠ open question of whether a requested `N` is counted before or after
  the mass filter, now BLOCKING in Remaining Work.

- **2026-09-10, step 6 — iteration-5 review, four CRITICALs, and a false-verification of my own.**
  *(1) The worst finding was not a number but a claim about a check.* R0 said
  `meta_fullsim_truncated` had been **verified unset in the regenerated file**. That branch is not
  in the file at all — upstream sets the flag in memory and never `Branch()`es it. The original
  check had run a `TTree::Scan` that printed `Bad numerical expression: "meta_fullsim_truncated"`
  and an empty column, and I read that as "it is a bool" instead of "it does not exist". The
  macro's read of it was correspondingly a silent no-op. **A verification that cannot have happened
  is worse than no verification.** Both the claim and the dead read are removed; `nproc != nbeam`
  is now stated as the only file-level truncation detector, and the file's 54 branches are
  enumerated as evidence.
  *(2) Two "fixes" recorded in the step-5 log had never been applied* — D4's straddle claim and the
  retired-axis edge `4.164`. The edit scripts that would have made them aborted on an unmatched
  pattern and wrote nothing, and the log recorded the intent rather than the outcome. Both are now
  applied and re-checked in place, and every edit in this step was verified against the file
  afterwards rather than assumed.
  *(3) §2 still carried a pre-R0-close sentence* defining `N_events` as the NTUP-chain entry count,
  spliced onto the ⚠ block and contradicting the paragraph below it. Deleted.
  *(4) `drop_sigma` vs `drop_sigma_wide` disagree in SIGN on all four `all selected pairs` rows*,
  which by the CSV's own stated rule makes those rows unquotable — and nothing marked them. A
  `window_dependence_ok` column now does, R2 states the failure explicitly, and the mechanism is
  recorded: inclusively the back-to-back peak sits INSIDE the narrow [10,25) window, whereas in the
  top-3 cells it is at 45–50 GeV and outside it. The accompanying honest caveat — that the top-3
  step therefore measures the shoulder, not the rise, and depends on the peak being outside the
  window — is now in R2 rather than left for a reader to notice.
  Also: Scope said six CSVs where the macro writes nine; the top-3 canvas label rounded 52.2273
  down to `52.2` with `%.1f` (now `%.2f`); Remaining Work still listed the completed regenerated-tree
  re-derivation; the INDEX line was dated 2026-09-09 and still advertised R0 as an open gap.

## Progress Log (continued)

- **2026-09-17, step 5 — same-sign-request preparation figures (user waived the reviewers:
  "urgent, do not run code/plot review").** New shared header
  `Utilities/MCRequestLooseMassFilter.h` (`MCRequest::kLooseMassMax = 10`), consumed by both
  macros. `mc_pthat_slice_mass_statistics.cxx`: 4-slice table (pTH24_40 = 803017, pTH40_70 =
  803018 added, DSID-guarded), `n_top_slices` + `figure_pair_pt_lo` arguments, figure-only mode
  with variant suffix, N-panel layout (2×2 for 4). All four slices passed the denominator
  cross-checks: N_proc = N_beam = NTUP chain = 2 399 995 / 1 199 997 / 319 999 / 319 999, weight
  implies N to 0.1 event, AMI totalEvents 2.4M / 1.2M / 320k / 320k. Outputs (Tight, in
  `mc_pthat_slice_mass_stats/`, nothing pre-existing touched):
  `pair_mass_by_pthat_slice_top4.png`, `pair_mass_by_pthat_slice_top4_pairpt20.png` (20 <
  p_T^pair < 150 GeV, user's figure cut — NOT a canonical edge), plus the by-product
  `pair_mass_by_pthat_slice_top3_pairpt_top4.png` and `mc_pthat_slice_mass_stats_top4{,_pairpt20}.root`.
  `plot_pythia_fullsim_kn_pt_crossx.cxx`: `ss_loose_mass` argument (D9) → 12 PNGs in
  `pythia_fullsim_full_sample/plots/ss_mass_lt10/` (`{truth,reco}_pair_pt_kn{,_stat_error,_err_frac}{,_150GeV}.png`),
  FULL sample, Tight; header text shrunk 0.038 → 0.032 for the longer label only.
  Input files: `..._mc_trig_full.root` (2026-09-10 04:34) and `..._full.root` (2026-09-10 01:21).
  **R7 observations** below.

## Results & Observations

All numbers below are **Tight** WP, raw pair counts, `N_slice = 319 999` events per slice.
Medium is 4–8 % higher inclusively and 2–6 % higher in the 3×3 cells (R6), and is in
`mc_pthat_slice_mass_stats_medium/`.
Percentages are **per-event rates**: multiply by the requested event count for the projection.

### R0 — Upstream gap found by review: **RAISED, FIXED UPSTREAM, AND NOW CONSUMED (CLOSED)**

`/review-analysis-code` established that the denominator's two original routes both reduced to
`N_beam` and neither could see `nevents_max` truncation, so a smoke-test pair file — which lands
on this same `_full` path — would have made every check agree while every per-event rate was
understated. Nothing in the pair file recorded `N_proc`.

**Resolution.** Raised with the session executing `mu_pt45_gap125_pairpt9_adoption.md` **before**
it submitted the NTuple rerun, so the fix could land in the regenerated trees at no extra cost.
It did (their D12): `PythiaAlgCoreT.c:839` now divides the weight by `N_proc`, and `meta_tree_out`
records `nproc_kin<K>_beam<B>` and `nbeam_kin<K>_beam<B>` on the fullsim path. **Verified in the
regenerated file**: `meta_tree_out` has 1 entry and 54 branches, with
`nproc_kin4_beam0 = nbeam_kin4_beam0 = 319 999` and the same for kin5.

**CORRECTION (2026-09-10, `/review-plot` + `/review-analysis-code` iteration 5).** An earlier
version of this paragraph also claimed `meta_fullsim_truncated` was verified unset **in the file**.
It is not in the file at all — upstream sets the flag in memory and never Branches it. The original
check had run a `TTree::Scan` that printed `Bad numerical expression: "meta_fullsim_truncated"` and
an empty column, and that output was misread as "the branch is a bool" rather than "the branch does
not exist". **A verification that cannot have happened is worse than no verification**, and the
macro's read of that branch was correspondingly a silent no-op. Both are removed:
`nproc != nbeam` is the only truncation detector and the code no longer pretends otherwise.

**This macro now consumes it.** `N_slice` is **read** from `nproc_kin<K>_beam0`, not inferred, and
the macro **refuses to run** on a pair file whose meta tree is missing or empty — with an error
naming the upstream fix — rather than silently mis-normalising. The three checks are now
cross-checks on a number read from the file: `nbeam` == NTUP chain entries; `σ·ε_filt / w` ==
`N_proc` (exact since the fix); AMI `totalEvents`. A truncated file is reported loudly instead of
passing. *(Note: fullsim files produced before 2026-09-09 carry an empty meta tree and will
trigger the refusal — that is the intended behaviour.)*

### R1 — The denominator, verified on the regenerated trees

`N_slice` is **read** from the pair file's `meta_tree_out`: `nproc_kin4_beam0` =
`nproc_kin5_beam0` = **319 999**, equal to `nbeam_kin{4,5}_beam0` — a complete, untruncated run.
(`meta_fullsim_truncated` is **not** a branch of this file; see R0.) All three cross-checks pass on both slices:

- **(a)** `nbeam` = 319 999 = the entries of the NTUP chain on disk;
- **(b)** `σ·ε_filt / w` = 319 999.01 (pTH70_125: 1.45881 nb / 4.55879e-06; pTH125_300:
  0.197783 nb / 6.18073e-07) — agrees with `N_proc` to 0.01 events;
- **(c)** AMI `totalEvents` = 320 000; the NTUPs hold one event fewer (0.0003 %, far inside the
  2 % tolerance).

Both σ·ε_filt match `ami_weights.md` table B exactly (1.4588, 0.1978 nb) and both AMI files passed
the **DSID guard** (803019, 803015).

### R2 — Where the back-to-back peak is (the filter question)

`pair_mass_by_pthat_slice.png` / `pair_mass_by_pthat_slice_top3_pairpt.png`.

- **Inclusively**, both slices show a broad back-to-back peak well above the proposed filter. The
  binning-robust locator is the **geometric-mean mass above 10 GeV**: **24.1 GeV** (pTH70_125 OS)
  and **27.4 GeV** (pTH125_300 OS) — it does move up with pT-hat, but only mildly. *(The
  tallest-bin comparison does NOT show this: on the current axis the two inclusive OS maxima are
  a near-tie, so quoting peak bins to show a shift would be reading a coin flip. Superseded claim,
  see the correction note.)*
- **In the 3 highest pair-p_T cells** — the cells the request exists for — the back-to-back
  population sits far higher: geometric-mean mass above the filter **34.1 GeV** (pTH70_125 OS) and
  **35.5 GeV** (pTH125_300 OS), with the tallest bin at [45.68, 50.20) = 174 pairs for
  pTH125_300 OS.

**CORRECTION, and the standing rule that came out of it (2026-09-09).** Three review iterations in
a row caught the *same class* of error here: a peak / minimum / density number transcribed into
this doc from one display axis and left behind when the axis changed (70 bins 0.2–220 → 80 bins
0.2–400 → 80 bins 0.211–400; every bin edge moves each time). Two specific claims were withdrawn:
"a clear valley opens at m ≈ 8–15 GeV" (there is no valley) and "the peak moves up with pT-hat"
as read from peak-bin locations (a tie).

**The structural fix, and the rule:** these landmarks are now **emitted by the macro** into
`mass_spectrum_landmarks.csv`, recomputed from the histograms of each run, with search regions
derived rather than hand-picked. **No spectrum landmark may be hand-transcribed into this doc
again — quote that CSV.** And two statistical rules that fell out of the same reviews:

1. **A single bin is a Poisson draw, not a measurement.** Comparing the tallest and shallowest bin
   of a ~20-bin scan is look-elsewhere biased — the extremes of a *flat* distribution are far apart
   by construction — so that comparison is emitted only under a column explicitly named
   `..._LOOK_ELSEWHERE_BIASED` and is **not quoted**.
2. **The verdict is the density step across the filter**, `drop_sigma` = (mean/bin below −
   mean/bin above) / combined error, over two regions fixed before looking. Positive = the density
   *falls* as the filter is crossed, i.e. the filter sits above the single-b population's edge and
   below the rise to the back-to-back peak — **not in a gap between them**.

**The measured density step across the proposed 10 GeV filter.**

**CORRECTED 2026-09-10 — the first version of this number was computed on a data-chosen window.**
`drop_sigma`'s "above" region originally ended at the *tallest bin above the filter*, while the
code and this doc both claimed the regions were "fixed before looking". That excluded the largest
bin by construction (biasing the step positive) and made the compared mass range differ from row
to row. Flagged as CRITICAL by `/review-analysis-code` and independently by the session running the
p_T-4.5 adoption. **The withdrawn claim is the one that had been promoted to the INDEX headline:
"0.0 σ — perfectly flat across the filter" for pTH125_300 same sign. It was an artefact of where
the window happened to end.**

The window is now pre-registered *and derived*: log-symmetric about the cut,
`[filter/R, filter)` vs `[filter, filter·R)` with **R = filter / template-fit-top = 10/4 = 2.5**,
i.e. **[4, 10) vs [10, 25) GeV**. Neither edge depends on the data. A second column,
`drop_sigma_wide` over [10, 400) GeV, is emitted beside it as an explicit window-dependence check.

Tight, 3 highest pair-p_T cells; pairs per bin on the log-uniform axis (= dN/dln m); source
`mass_spectrum_landmarks.csv`:

| population | [4, 10) | [10, 25) | **step** | [10, 400) | step (wide) |
|---|---|---|---|---|---|
| pTH125_300 **OS** | 144.8 ± 4.3 | 115.7 ± 3.6 | **+5.2 σ** | 83.7 ± 1.5 | +13.6 σ |
| pTH70_125 **OS** | 18.4 ± 1.5 | 12.8 ± 1.2 | **+2.9 σ** | 9.1 ± 0.5 | +5.9 σ |
| pTH125_300 **SS** | 59.9 ± 2.7 | 53.2 ± 2.4 | **+1.8 σ** | 47.2 ± 1.1 | +4.3 σ |
| pTH70_125 **SS** | 3.2 ± 0.6 | 4.8 ± 0.7 | **−1.6 σ** | 5.3 ± 0.4 | −2.8 σ |

**How to read the two columns, and where the check FAILS.** The **magnitude** is not expected to
agree and is not quotable: the wide window reaches into the steeply falling high-mass tail, so it
always gives the larger number. Only the **sign** is the robustness claim, and it agrees in all
four rows of the table above.

**It does NOT agree in the four `all selected pairs` rows** (−25.7 vs +19.2, −17.5 vs +40.1,
−15.6 vs +22.8, −3.2 vs +45.8). By the rule stated in the CSV header, **the inclusive rows' step is
not quotable and nothing in this doc quotes it**; the CSV now carries a `window_dependence_ok`
column marking them so. The mechanism is understood and is not a defect: **inclusively the
back-to-back peak sits at 18–26 GeV, INSIDE the narrow [10, 25) window**, so that window is
dominated by the peak itself. In the 3 highest pair-p_T cells the peak has moved to 45–50 GeV,
**outside** the narrow window, so there the narrow window measures the shoulder just above the
filter — which is precisely the quantity the filter question asks about.

**The honest caveat that follows.** Because the top-3 narrow window sits below the back-to-back
peak, `drop_sigma` there measures the 10–24 GeV shoulder and not the rise to the peak. That is the
right quantity for "is the filter at an edge?", but it does mean the positive step depends on the
peak being outside the window — which is uncomfortably adjacent to the data-dependence the
pre-registration was introduced to remove. It is mitigated by the window being fixed by canonical
mass scales rather than by the data, and by the wide column agreeing in sign; it is not eliminated.

**What it says, and it is weaker than the first version claimed.** For **opposite sign** the
density genuinely falls as the filter is crossed — +5.2 σ and +2.9 σ — so 10 GeV does sit above
the single-b population's edge. For **same sign** the step is **marginal at best**: +1.8 σ in the
harder slice, and in the softer slice it goes the *wrong way* (−1.6 σ, the density is higher above
the cut than below it). So the sign asymmetry that matters for the request survives the
correction — same-sign pairs at high pair p_T show little or no edge at 10 GeV — but the earlier
"perfectly flat" phrasing overstated how clean that null was.

- **What survives the proposed value, in the 3 highest pair-p_T cells.**
  Source: `mass_composition_top3_<slice>.csv`, **exact `minv` cuts** (D7).

  **CORRECTION (2026-09-09, second round).** The first version of this table was integrated off
  the display histogram's bin edges (2.86 / 4.16 / 9.84 instead of 2.9 / 4 / 10) and was wrong by
  up to **6.2 percentage points** — the headline same-sign number read 48.4 % where the truth is
  54.55 %. Both reviewers caught it independently and their exact values agree with the CSV the
  macro now produces. This violated D4's own "no table depends on the display axis"; D7 is the
  structural fix.

  Counts, and % of all selected pairs in those cells:

  | | m < 2.9 | 2.9–4 | 4–10 | m ≥ 10 (removed) | all |
  |---|---|---|---|---|---|
  | pTH125_300 OS | 10 546 (58.27 %) | 2 874 (15.88 %) | 1 392 (7.69 %) | 3 288 (**18.17 %**) | 18 100 |
  | pTH70_125 OS | 2 307 (58.96 %) | 1 065 (27.22 %) | 185 (4.73 %) | 356 (**9.10 %**) | 3 913 |
  | pTH125_300 SS | 287 (10.13 %) | 128 (4.52 %) | 571 (20.16 %) | 1 847 (**65.20 %**) | 2 833 |
  | pTH70_125 SS | 21 (7.72 %) | 9 (3.31 %) | 36 (13.24 %) | 206 (**75.74 %**) | 272 |

  Composition of **what the filter keeps** (the first three bands; these totals are exactly the
  3×3 `m < 10` totals of R5, which is the consistency check between the two tables):

  | | kept | below 2.9 | below 4 | on the 4–10 continuum |
  |---|---|---|---|---|
  | pTH125_300 OS | 14 812 | 71.20 % | 90.60 % | **9.40 %** |
  | pTH70_125 OS | 3 557 | 64.86 % | 94.80 % | **5.20 %** |
  | pTH125_300 SS | 986 | 29.11 % | 42.09 % | **57.91 %** |
  | pTH70_125 SS | 66 | 31.82 % | 45.45 % | **54.55 %** |

  **Caveat on the 2.9–4 band — it is not continuum.** This file carries no resonance veto (§3b),
  and inside the analysis's own J/ψ veto window `ParamsSet::minv_cuts` = [2.9, 3.3) GeV sit
  **84.6 %** (pTH125_300 OS) and **87.3 %** (pTH70_125 OS) of that band — 2 432 of 2 874 and 930
  of 1 065 pairs. So the opposite-sign "90–95 % below 4 GeV" is roughly two-thirds signal-window
  pairs plus a large unvetoed J/ψ, not a smooth low-mass continuum. Same-sign is different again
  (38.3 % / 55.6 % of a much smaller band), as it must be — a J/ψ cannot decay to a same-sign
  pair, so what is there is combinatoric.

- **Conclusion on the proposed value — it meets requirement (2) cleanly and requirement (1) only
  for opposite sign.**
  - *Not biasing (requirement 2): satisfied.* 10 GeV stands ~3.4× above the signal window's
    2.9 GeV edge and ~2.5× above the 4 GeV top of the template-fit window, so neither measurement
    is conditioned on a cut the analysis does not apply, and the ΔR range those windows populate
    is untouched.
  - *Removing the back-to-back pairs (requirement 1): only partly.* Every percentage below is
    **of what the filter keeps** (the `m < 10` sample), stated once so it cannot be confused with
    the fractions-of-all-pairs column above.
    For **OS it works**: **90.6 %** (pTH125_300) and **94.8 %** (pTH70_125) of the kept sample is
    below 4 GeV, and only **9.40 %** and **5.20 %** sits on the 4–10 GeV continuum.
    For **SS it does not**: **57.91 %** (pTH125_300) and **54.55 %** (pTH70_125) of the surviving
    same-sign pairs sit at 4–10 GeV — i.e. more than half of what the filter keeps is on the
    back-to-back continuum — and only 29–32 % are below 2.9 GeV. The same-sign sample at high
    pair p_T is back-to-back-dominated to begin with (65–76 % of *all* its pairs are above
    10 GeV), and a 10 GeV filter does not change that character. The density step says the same, though
    only weakly: **+1.8 σ** (pTH125_300) and **−1.6 σ** (pTH70_125) on the pre-registered window —
    marginal in one slice and the wrong way in the other, against **+5.2 σ / +2.9 σ** for opposite
    sign. For same sign there is little or no edge at 10 GeV to cut on.
  - **This is a decision point for the user, not one to settle here** (§Remaining Work): a tighter
    filter on the same-sign-filtered request (m < 4–5 GeV) would buy far more usable SS pairs per
    event, at the cost of the mass headroom above the template-fit window. It does **not** affect
    the sign-independent request, where 10 GeV is fine.

### R3 — The mass filter removes MORE than half of the same-sign pairs

Fraction of all selected pairs kept by `m < 10 GeV`:

| slice | same sign | opposite sign |
|---|---|---|
| pTH70_125 | 11 662 / 28 508 = **40.9 %** | 83 263 / 109 839 = **75.8 %** |
| pTH125_300 | 16 320 / 36 858 = **44.3 %** | 92 945 / 122 406 = **75.9 %** |

The majority of same-sign pairs in these slices **are** back-to-back. This is the quantitative
form of the physics motivation: without the filter, a same-sign efficiency cell would be
dominated by a population the analysis never uses.

### R4 — Table 1: mass ranges, both signs, no pair-p_T/η restriction

`mass_window_counts_<slice>.csv`. Count (per-event rate); the two rows OVERLAP.

**pTH70_125**

| mass range | same sign | opposite sign |
|---|---|---|
| m < 10 GeV | 11 662 (3.6444 %) | 83 263 (26.0198 %) |
| 1.08 < m < 2.9 GeV | 3 132 (0.9788 %) | 36 091 (11.2785 %) |

**pTH125_300**

| mass range | same sign | opposite sign |
|---|---|---|
| m < 10 GeV | 16 320 (5.1000 %) | 92 945 (29.0454 %) |
| 1.08 < m < 2.9 GeV | 4 946 (1.5456 %) | 37 923 (11.8510 %) |

Inclusively the two slices are within ~12 % of each other **for opposite sign** (+11.6 % at
m < 10 GeV, +5.1 % in the signal window); **same sign differs far more** (+40 % and +58 %). **The slices only separate in the top
pair-p_T cells** (R5) — which is exactly why an inclusive rate is the wrong number to size this
request on.

### R5 — Tables 2–5: the 3 highest pair-p_T cells × 3 |η^pair| groups

`cells_<slice>_<range>.csv`. Count (per-event rate). Cells are
[52.23, 74.24) / [74.24, 105.53) / [105.53, 150) GeV × |η^pair| < 1 / 1–2 / 2–2.2.

**pTH125_300, m < 10 GeV** — opposite sign:

| p_T^pair | \|η\| < 1 | 1–2 | 2–2.2 |
|---|---|---|---|
| 52.2–74.2 | 6 361 (1.9878 %) | 3 007 (0.9397 %) | 314 (0.0981 %) |
| 74.2–105.5 | 2 831 (0.8847 %) | 1 254 (0.3919 %) | 126 (0.0394 %) |
| 105.5–150 | 643 (0.2009 %) | 257 (0.0803 %) | 19 (0.0059 %) |

same sign: 527/264/13, 115/47/0, 17/3/0 — total **986 (0.3081 %)** vs OS **14 812 (4.6288 %)**.

**pTH70_125, m < 10 GeV** — opposite sign 1 830/1 049/157, 288/180/22, 17/14/0 (total
**3 557**, 1.1116 %); same sign 34/27/3, 2/0/0, 0/0/0 (total **66**, 0.0206 %).

Signal-window (1.08–2.9 GeV) 3×3 totals: pTH70_125 SS **19** / OS **1 485**;
pTH125_300 SS **251** / OS **6 968**.

**Three things this table decides:**

1. **pTH125_300 is the slice to request.** At equal event count it gives **4.2× more** OS pairs
   in the 3×3 cells than pTH70_125 (14 812 vs 3 557) and **15× more** SS (986 vs 66).
2. **pTH70_125 cannot populate the top cell at all.** Its [105.5, 150) row is 17/14/0 OS and
   0/0/0 SS — after the signal-window cut, 7/4/0 and 0/0/0. No amount of the 70–125 slice fixes
   the highest pair-p_T cell; only 125–300 does.
3. **The same-sign request must be filtered on same-sign, not merely larger.** In the 3×3 cells
   SS is **1.9 %** (pTH70_125) and **6.7 %** (pTH125_300) of OS. An unfiltered request would need
   ~15× the events to give the SS cells the statistics the OS cells already have — which is the
   direct justification for the second, same-sign-filtered sample.

### R6 — Working-point dependence

Medium/Tight ratio of the 3×3 totals: pTH125_300 OS 15 095/14 812 = 1.019, SS 1 011/986 = 1.025;
pTH70_125 OS 3 630/3 557 = 1.021, SS 70/66 = 1.061. Inclusively **1.040–1.076** (the maximum is
pTH70_125 SS in the signal window, 3 132 → 3 371). The WP choice does not change any conclusion
above.

### R7 — Same-sign request preparation (2026-09-17, figures only; landmarks NOT emitted for these)

- **Inclusive, 4 slices** (`pair_mass_by_pthat_slice_top4.png`): the same-sign back-to-back
  peak moves DOWN with pT-hat — by eye at ≈ 15–20 GeV (24–40), ≈ 20 GeV (40–70), ≈ 25 GeV
  (70–125, 125–300); in the two added slices the SS continuum rises monotonically from 1 GeV
  up to the peak with NO shoulder at 10 GeV. Same-sign is ~10× below opposite-sign below 3 GeV
  in every slice.
- **20 < p_T^pair < 150 GeV** (`..._top4_pairpt20.png`): the SS distribution is essentially
  FLAT from ~3 to ~30 GeV in all four slices (24–40: ~100/bin; 40–70: ~300/bin; 70–125:
  ~250/bin; 125–300: ~500/bin), then falls — so at this pair-pT cut 10 GeV again sits in the
  middle of a plateau, not at an edge. Per the standing rule these are eyeball readings of the
  figure, not emitted landmarks; the mass-cut decision for the SS request is the user's (see
  Remaining Work).
- **SS m < 10 GeV statistics** (`ss_mass_lt10/`): the reco SS rel. stat. error reaches ~10 %
  per fine bin at ≈ 30 GeV (kn3) / ≈ 40 GeV (kn4) / ≈ 60 GeV (kn5) and the error-fraction map
  shows 125–300 dominating above ~80 GeV, 70–125 at 40–80 GeV, 40–70 at 25–40 GeV — i.e.
  ALL FOUR slices matter for a request reaching down to 20–25 GeV, which is the justification
  for widening the slice range beyond the first request's two.

## Remaining Work

- **BLOCKING for the request itself — is the requested `N` counted before or after the mass
  filter?** §2's ⚠ block: `r · N` is exact only if `N` counts events *before* the filter. If it
  counts events *after* (the usual convention), the yield is larger by `1 / ε_massfilter`, which
  is a generator-level number this study does not measure. Settle the convention, and measure or
  obtain ε_massfilter, before quoting any projected yield in the request.
- **USER DECISION — the same-sign-filtered request's mass cut.** R2 establishes that `m < 10 GeV`
  leaves the same-sign sample back-to-back-dominated (54.6–57.9 % of surviving SS pairs at 4–10 GeV,
  only 29–32 % below 2.9 GeV), while it is clean for opposite sign. Whether the SS-filtered
  request should carry a tighter filter (m < 4–5 GeV) is a physics trade-off — more usable SS
  pairs per event vs. mass headroom above the template-fit window — and must be put to the user.
- **Still to state explicitly:** whether the ΔR range surviving the chosen filter supports the
  **ΔR-correction** procedure, not only the single-value one. R2 argues it does on mass-window
  grounds (nothing inside 1.08–4 GeV is touched), but the ΔR distribution of the kept sample in
  the top pair-p_T cells has not been drawn. That is the direct test and it is not yet done.

## Latest Stage

**Step 5 DONE (2026-09-17).** Figures delivered and committed; no reviewer pass (user waived).
Open: the SS request's slice list / pair-pT cut / mass cut are user decisions (Remaining Work).

Plan as executed: (i) shared header
`Utilities/MCRequestLooseMassFilter.h` (kLooseMassMax = 10); (ii)
`mc_pthat_slice_mass_statistics.cxx`: 4-slice table with DSIDs, `n_top_slices` +
`figure_pair_pt_lo` arguments, figure-only mode, 2×2 canvas for N = 4 (subplot rule), variant
suffixes; run `(pp_full, true, 4)` and `(pp_full, true, 4, 20)`; (iii)
`plot_pythia_fullsim_kn_pt_crossx.cxx`: `ss_loose_mass` argument → tree/filter/label/outdir
helpers, run FULL sample Tight; (iv) reviewers; (v) results → R7; commit by explicit path.
Files: the two macros, the new header, this doc, INDEX.md.

Previous stage (step 4, 2026-09-10). Rerun on the regenerated 4.5 GeV trees; all counts unchanged (predicted and now
confirmed). Five review iterations have run (`/review-analysis-code` and `/review-plot`, twice
each in parallel, then a combined adversarial pass). Every count in the deliverable has been
independently re-derived from a retyped selection **four times** and matched exactly every time;
every failure was in *derived statements* — a claimed valley, a table read off display-bin
edges, landmarks left on a superseded axis, and a significance computed on a data-chosen window.
Each has been fixed structurally rather than patched: the composition bands and the spectrum
landmarks are now emitted by the macro, and the density-step window is pre-registered and derived.

Next: iteration-5 review of this state, then commit by explicit path (macro + this doc + the
INDEX line; the working tree is shared, so never `git add -A`).

**Open items** — see Remaining Work: (1) **BLOCKING** — is the requested `N` counted before or
after the mass filter? ε_massfilter is not measured here; (2) USER DECISION — the same-sign
request's mass cut, now supported by the corrected +1.8/−1.6 σ step; (3) the ΔR distribution of
the kept sample in the top cells is still not drawn, and that is the direct test of whether the
filter leaves the ΔR-correction procedure measurable.
