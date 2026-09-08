# Single-value pair 2mu4 efficiency in the top pair-p_T cells (mass-windowed)

**Mode:** Implementation. **Created:** 2026-09-08.
**Session:** "single-value pair 2mu4 efficiency".
**Parents:** `mc_trigger_efficiency.md` (ACTIVE — Step 3, R25/R26/R32),
`mc_trigeff_dr_binning_approaches.md` (ACTIVE — the four cell groupings and the closure this doc
extends), `mc_trig_eff_closure.md` (ACTIVE — the closure machinery),
`muon_gap_cuts_acceptance.md` (ACTIVE — the 2026-09-07 cut set every input here was re-derived on).

---

## Objective

Measure, on the pp24-conditions Pythia fullsim FULL sample, a **single 2mu4 efficiency per
(pair p_T bin, |eta^pair| bin, pair sign) cell inside a dimuon mass window**, for the three
highest pair-p_T cells, as an ALTERNATIVE to the current factorized weight
`eps(p_T1,q*eta_1) * eps(p_T2,q*eta_2) * eps_dR(dR ; cell)` whose dR-shape fit is
statistics-starved there. Deliver the numbers in a ROOT file the crossx pipeline can read, and
compare the two procedures on the MC closure, restricted to p_T^pair > 50 GeV.

## Autonomy Contract (DONE 2026-09-08 — every item met; kept for the record)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. A ROOT file holding eps_2mu4^pair and the data-calibrated K for every
     (pair-p_T coarse bin, |eta^pair| group, sign, mass window, WP) cell, with conditional
     (binomial-correct) errors and raw pair counts, plus a loader header so the crossx pipeline
     can retrieve a cell by (pair p_T, pair eta, sign, window) — NOT wired into the cross-section.
  2. The three requested pair-p_T bins [49.97,72.08), [72.08,103.98), [103.98,150) measured for
     BOTH mass windows (1.08-2.9 and 1.0-4.0 GeV) and BOTH signs; [49.97,72.08) is the control
     region where the dR procedure also works.
  3. The 5-line closure comparison figure on the pp24 crossx binning
     (`ParamsSet::pT_bins_150` x the 9 `pair_eta_proj_ranges_coarse_incl_gap` panels), zoomed to
     p_T^pair > 50 GeV: no-trigger reference + approach B (last two pair-p_T cells merged) +
     approach D (B plus the |eta| fold) + the two single-value mass-window series.
  4. Numbers recorded in this doc (per-cell efficiencies, counts, inclusive closure).
  5. `/review-analysis-code` on the C++/RDF, `/review-plot` on the figure; this doc + INDEX
     updated; commits.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The pp24 per-pair trigger weight is factorized,

```
eps_trig^pair = eps(p_T1, q*eta_1) * eps(p_T2, q*eta_2) * eps_dR^2mu4(dR ; pair-pT cell, pair-eta cell)
```

— two single-muon turn-ons times a two-body **close-by correction fitted as a function of dR** in
each cell. In the top pair-p_T cells the pp24 fullsim sample has too few pairs to constrain that
shape: fits are rejected outright, `polyu` produced an eps_dR = 2.74 artefact, and even after the
last-two-bins merge one cell is corrected by the raw-bin placeholder at eps_dR = 2.157 where the
physics requires <= 1 (parent `mc_trigger_efficiency.md` R26/R32; `mc_trig_eff_closure.md` R5).

**The alternative measures the same physics without the fit and without the factorization.** At
fixed pair p_T and pair eta the dimuon MASS essentially fixes the opening angle
(dR ~ 2m/p_T for a roughly symmetric decay) and constrains the leg-p_T sharing, so a mass window
plus a (p_T, |eta|) cell already pins down most of the pair kinematics the trigger responds to.
Everything the factorized product is built to track — the singles' turn-on and then the L1
close-by-RoI loss that breaks factorization — is contained in ONE directly measured ratio. It is
the statistically cheapest object the same events can produce: one binomial number instead of a
curve.

**What it costs, stated rather than hidden:** all within-cell differential information. The
correction is exact only INTEGRATED over its cell, and it inherits the MC pair spectrum inside
that cell, so a data/MC spectrum difference within a cell biases it. The factorized form is
immune to that because it is evaluated per pair. The closure figure (§PP-4), drawn on the FINER
crossx presentation binning, is precisely the read-out of that within-cell residual.

**Why a mass window, and why two.** The number is only valid for the mass mixture it was measured
on. Window 1 = the signal window, so it is the right number for the signal region. Window 2 is
what a template fit over 1-4 GeV would need and is a DIFFERENT dR mixture at the same pair p_T.
**The difference between the two versions is the measurement of the mass dependence** and decides
whether one efficiency may serve the wider fit window.

**Why the signs are separated.** Parent R25 established that the dR correction is genuinely
charge-dependent at small dR (same-sign/opposite-sign = 0.43 at dR = 0.025, an L1 close-by-RoI
signature present in the RAW trigger probability), and the background estimate is OS - SS. A
same-sign pair must carry its own number.

### 2. Top-level equations

Per cell = (pair-p_T coarse bin) x (|eta^pair| group) x (pair sign), on the pairs of one mass
window, with w = the per-pair MC weight (sigma_slice * eps_filt * r_isospin / N_slice):

```
S0 = sum_{all pairs}         w
S1 = sum_{pairs firing 2mu4} w
S2 = sum_{pairs firing 2mu4} w / [eps_MC(1) * eps_MC(2)]

eps_2mu4^pair(cell) = S1 / S0                (PURE:  replaces the WHOLE weight)
K(cell)             = S2 / S0                (CALIBRATED: multiplies the two single-muon effs)
```

Both are defined so that the corresponding estimator is **unbiased when integrated over the
cell**, which is what fixes the internal weighting:

- applying `w_corr = w / eps_2mu4^pair` to the firing pairs returns `S0` exactly;
- applying `w_corr = w / [eps(1) eps(2) K]` returns `S0` exactly as well, because
  `S2 = sum_pass w/(eps1 eps2)` is the MC estimate of `sum_all w * P_pair/(eps1 eps2)` and the
  unbiased single factor is the `1/(eps1 eps2)`-weighted mean of `P_pair/(eps1 eps2)`.

**Why K exists at all (and why it, not eps^pair, is what the cross-section should apply).**
`eps_2mu4^pair` is a pure-MC number, and the analysis deliberately takes its single-muon
efficiencies from DATA tag-and-probe, using MC only for the correlation. The MC/data pair
over-efficiency has been measured in this very closure: `<r_1 r_2> = 1.27` (all-OS, Tight,
barrel-driven; `mc_trig_eff_closure.md` R1). Applying `eps_2mu4^pair` directly on data would
therefore bias those cross-section bins by ~25 %, and an MC closure CANNOT see it (the closure
applies eps_MC and is self-contained by construction). `K` is the single-number analogue of
eps_dR: it is measured with eps_MC in MC exactly as the Step-3 correction is, and applied on top
of the data eps^nc exactly as the Step-3 correction is. **Both numbers are written; neither is
wired into the cross-section by this doc.** (User, 2026-09-08: proceed on the stated default —
deliver K as the recommended applied form and eps^pair alongside.)

### 3. Step-by-step method

#### §PP-1 The cells

- pair p_T: `ParamsSet::pair_pt_coarse_bins` — the canonical 8 log bins, 8-150 GeV. The
  measurement is made in **all 8** (it costs nothing and the low bins are the sanity check); the
  three the request is about are bins 6/7/8 = **[49.97, 72.08), [72.08, 103.98), [103.98, 150)**
  GeV, of which bin 6 is the **CONTROL REGION** where the dR procedure also has statistics.
- |eta^pair|: the **3 sign-independent groups** `MakeDrEtaGroups(..., true)` builds from
  `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` —
  **|eta| < 1.0 / 1.0 <= |eta| < 2.0 / 2.0 <= |eta| < 2.2** (the outer edge is 2.2 since the
  2026-09-07 pair-level fiducial cut, `muon_gap_cuts_acceptance.md` F17). Same helper as the
  `nocorr_etamerge*` dR modes, so the two groupings cannot drift.
- sign: opposite sign (`muon_pair_tree_sign2`) and same sign (`muon_pair_tree_sign1`), separately.

**No binning is invented here.** Both axes are read from their canonical source; the |eta| fold is
the existing `MakeDrEtaGroups` (.claude/CLAUDE.md §Binnings items 1-2).

#### §PP-2 The pair sample and the mass windows

Base selection = `MCTrigEffPairSel::Step3PairSelection(true)`, i.e. EXACTLY the sample the dR
correction is measured on — truth-seeded reco-matched Pythia legs, both legs' WP bit and
p_T > 4 / |eta| < 2.4, the truth fiducial, the forward low-p_T veto, the per-leg fiducial gap
windows and the pair-level |eta^pair| < 2.2. On top of it, one of two MASS WINDOWS:

| token | window | what it is |
|---|---|---|
| `sig`  | 1.08 < m < 2.9 GeV | the single-b signal window — with the base selection this IS the signal region (`SingleBSignalCutsReco` adds only `pair_pt > 8` and the gap cut, both already in) |
| `wide` | 1.0 < m < 4.0 GeV | the template-fit window, containing the phi, J/psi and psi(2S) mass region |

The MC pair file carries NO resonance veto (`..._no_data_resonance_cuts_mc_trig_full.root`), so
the `wide` window is the honest 1-4 GeV mixture the MC produces. **Caveat to carry forward:** that
mixture is HF continuum as Pythia makes it; if the data's 1-4 GeV mixture is materially more
resonant, the applicable average shifts. A mass-differential eps^pair table is printed as the
diagnostic that would show it.

#### §PP-3 Errors

`eps^pair = S1/S0` and `K = S2/S0` are efficiency-like ratios whose numerator is a re-weighted
SUBSET of the denominator, so their errors are the CONDITIONAL/binomial ones
(`SetConditionalRatioErrors`, `dr_correction_ratio.h`), built from the same
`A = sum a^2`, `B = sum a^2 p` sums the Step-3 correction uses — never `TH1::Divide`'s
independent-error form, which over-states these bars by 1.5-3.2x (parent R12).

#### §PP-4 The comparison figure

The MC closure of `mc_trig_eff_closure.md` §2 with the applied single-muon efficiency eps_MC,
binned on the **pp24 cross-section's** axes (`ParamsSet::pT_bins_150`, 15 log bins 8-150 GeV, x
the 9 pair-eta panels — `mc_trigeff_dr_binning_approaches.md` §PP-1/D8), **x-range restricted to
p_T^pair > 50 GeV**, five series:

| series | weight applied to the firing pairs | colour |
|---|---|---|
| no trigger requirement (the reference) | — | black |
| **B** dR procedure, last two pair-p_T cells merged (`nocorr_ptmerge`) | `w / [eps_MC(1) eps_MC(2) eps_dR]` | kBlue |
| **D** dR procedure, B plus the \|eta\| fold (`nocorr_etamerge_ptmerge`) | same, other cells | kMagenta |
| **single value, `sig` window** | `w / eps_2mu4^pair(cell)` | kOrange+7 |
| **single value, `wide` window** | `w / eps_2mu4^pair(cell)` | kAzure+7 |

A SECOND figure, identical but with the two single-value series in their CALIBRATED form
`w / [eps_MC(1) eps_MC(2) K(cell)]`, is produced beside it — that is the form §2 recommends
applying, and it is the apples-to-apples comparison against B and D (same singles, different pair
correction), while the pure form is the procedure as the user stated it. Neither is a new
measurement: both come from the same fill.

**Two honest statements about what this figure can and cannot decide.**
1. On the `signal` sample version with the `sig` window, the single-value series closes to
   **exactly 1 by construction when integrated over its own coarse cell** (§2). It is NOT
   trivially 1 on the finer crossx presentation bins, so the figure is a genuine test — of the
   RESIDUAL WITHIN-CELL kinematic dependence, not an unbiased head-to-head ranking.
2. The `wide` window carries no such construction on the signal sample: its deviation from 1 IS
   the mass-window bias, and is the unbiased number the figure delivers.

**Coverage.** The single-value correction exists only inside the coarse cells it was measured in.
The presentation bin that STRADDLES the 49.97 GeV cell edge ([46.44, 56.46) of `pT_bins_150`) is
therefore drawn for the dR series only; the single-value series start at the first fully covered
presentation bin. Drawing them in a partially covered bin would show a coverage artefact as if it
were non-closure.

### 4. Negative constraints

- **Nothing here changes what the pp24 cross-section applies.** `DrCorrCrossxMethod/Sign/Mode`
  and `DrCorrectionCrossxEvaluator` are untouched; the new file is produced and readable, not
  wired in.
- **Do NOT change any binning.** Both axes are read from `ParamsSet` / `CommonEffcyConfig` through
  the existing helpers.
- **Do NOT re-derive the sample from raw NTUPs** — everything reads the ntuple-processing output
  through `MCTrigEffPairSel::Step3PairSelection()` (.claude/CLAUDE.md §NTuple-Processing
  Provenance).
- **eps^pair is NOT eps_dR** and must never be routed through `DrCorrectionEvaluator`: it replaces
  the whole weight, carries no dR dependence, and is defined only in the cells it was measured in.
- Pb+Pb is out of scope (the mu4 UNION weight is a different formula).

---

## Scope

**In:** pp24 fullsim FULL sample, 2mu4, both signs, both mass windows, all 8 coarse pair-p_T bins
(the top 3 are the deliverable), the 3 |eta^pair| groups, Tight WP nominal with Medium reachable
by one switch; the new measurement, its ROOT deliverable + loader, the closure extension, and the
zoomed comparison figure.

**Out:** Pb+Pb / the overlay; wiring anything into the cross-section; re-deriving the dR
corrections (the fresh 2026-09-07/08 fit files are used as they are); the data/MC scale factor
that a pure-eps^pair application would need.

## Design Decisions

### D1: the deliverable carries BOTH `eps_2mu4^pair` and `K` (user default accepted, 2026-09-08)
The user's procedure is the pure single number; the analysis's calibration chain needs the K form
(§2). They come from the same three sums, so both are measured, written and plotted; the choice of
what to apply is left to the user with the numbers in hand.

### D2: measured in all 8 coarse pair-p_T bins, delivered for the top 3
Filling the whole canonical axis costs nothing, gives the low-p_T cells as a sanity check against
the dR procedure where the latter is well determined, and avoids inventing a 3-bin axis.

### D3: a separate fill macro + loader, not an extension of the dR-correction chain
`eps^pair` is a different object from `eps_dR` (no dR dependence, whole-weight replacement,
mass-windowed, sign-separated). Routing it through `DrCorrectionEvaluator`/the plateau-mode tokens
would overload a class whose every consumer assumes `f(dR)/C`. New:
`RDFBasedHistFilling/FillMCTrigEffPairEff.cxx` + `Utilities/PairTrigEffEvaluator.h`.

---

## Implementation Plan

1. [x] `Utilities/PairTrigEffEvaluator.h` — the cell definition (mass windows, axes), the writer's
   naming convention, the loader + canonical-binning guard. Per §PP-1/§PP-2.
2. [x] `RDFBasedHistFilling/FillMCTrigEffPairEff.cxx` — the measurement (§2, §PP-2, §PP-3);
   writes the ROOT deliverable + prints the per-cell tables and the mass-differential diagnostic.
3. [x] `FillMCTrigEffClosure.cxx` — the two single-value series (pure and calibrated) as extra
   numerators, per §PP-4. Coverage handling per §PP-4.
4. [x] `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_closure_highpt_compare.cxx` — the
   5-line zoomed figure, both forms. → `/review-plot`.
5. [x] Driver `pipelines/run_mc_trigeff_pair_eff.sh`; run everything; record numbers here.
6. [x] `/review-analysis-code` — **PASS at iteration 3** (2 FAIL rounds first; the first found
   the CRITICAL cascade fold bug, the second the coverage-tolerance one). `/review-plot` — **PASS at
   iteration 2**. INDEX updated; committed as `455fe32`, `7f7e5f2`, `1b249f2`, `849acd8`, `a3f9711`.

## Progress Log

*(append-only; newest entries at the END)*

- 2026-09-08 — **Steps 1-5 DONE: the measurement, the deliverable, the closure extension, the
  comparison figure and the driver are written, compiled and RUN (Tight).** Files:
  * NEW `Utilities/PairTrigEffEvaluator.h` — the cell definition (mass windows; pair-pT axis read
    from `ParamsSet::pair_pt_coarse_bins`; the |eta^pair| axis built by the SAME
    `MakeDrEtaGroups(..., true)` fold the `nocorr_etamerge*` dR modes use), the on-disk naming, and
    the reader (`Covered()` gates `Eval()`; `CheckCanonicalBinning()` re-derives both axes from
    their live sources and throws on a mismatch; `CheckSignalWindowMirror()` asserts the `sig`
    window still appears verbatim in `MCTrigEffPairSel::SingleBSignalCutsReco()`).
  * NEW `RDFBasedHistFilling/FillMCTrigEffPairEff.cxx` — the measurement (§2), 44 TH2D +
    provenance -> `pair_trig_eff_pp24_full.root`, plus the per-cell tables and the
    mass-dependence diagnostic.
  * MODIFIED `RDFBasedHistFilling/FillMCTrigEffClosure.cxx` — four extra numerators
    (`paireff_{sig,wide}` pure, `paireffK_{sig,wide}` calibrated) and the coverage denominators
    `den_paireff_<window>`, booked for the **`signal` version only** (the mass window is part of
    the definition, so the all-OS sample is a different mixture). The dR path is untouched.
  * NEW `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_closure_highpt_compare.cxx` — the
    5-line zoomed figure, one per applied form.
  * NEW `pipelines/run_mc_trigeff_pair_eff.sh` — 4 stages with artefact validation.
  Run: `[00:45:43]..[00:49:19]`, exit 0; logs in `pipelines/logs_pair_eff/`.

- 2026-09-08 — **RESULTS.** See §Results & Observations. Headlines: the single-value procedure
  closes to **1.0000 exactly** against its own coverage denominator (the construction check), and
  per pair-eta panel to **0.88-1.11** above the cell edge — as good as the dR procedure's
  0.93-1.11 with three times fewer cells. The `wide` (1-4 GeV) window applied to signal-region
  pairs closes to **0.8607**, i.e. a **14 % under-correction**: the trigger efficiency's mass
  dependence is real and large. Approach D (the |eta| fold) closes to **0.7433** above 50 GeV, and
  its failure is **entirely in the negative-eta panels** (0.44 at [-1.0,-0.5) against 0.94 at its
  mirror [0.5,1.0)).

- 2026-09-08 — **MEDIUM working point run** (`WPS="medium" ./run_mc_trigeff_pair_eff.sh`,
  `[00:56:54]..[01:00:24]`, exit 0) →
  `pair_trig_eff_pp24_full_medium_wp.root` + `mc_based_medium/closure/single_value_highpt_comparison/`.
  **The result is WP-independent**, as the mu4 response has been throughout this thread: inclusive
  closure above the cell edge 1.0015 (dR p_T-merged) / 0.9393 (|eta| fold, POST fold-bug fix; the pre-fix number was 0.7441) / **1.0000** (single
  value, sig) / **0.8629** (single value, wide), against Tight's 1.0003 / 0.9362 / 1.0000 / 0.8607;
  per panel the two working points agree to 0.026 at worst (the `wide` series in the 0.5 < eta^pair < 1.0 panel, 0.7863 Tight vs 0.8122 Medium). The mass-window bias, the |eta|-fold
  sign asymmetry and the per-panel spread are therefore not WP artefacts.

- 2026-09-08 — **One defect found and fixed in this session, in the printed table only (the figure
  was always right).** The first version of the per-panel table divided the single-value numerator
  by the FULL no-trigger denominator, while that numerator only covers pair p_T >= 49.97 GeV; the
  presentation bin straddling the cell edge holds ~29 % of the yield above it, so the table read
  0.71 where the true closure is 1.0000. Fixed by giving each single-value series its own matched
  `den_paireff_<window>` denominator in the table (the figure instead OMITS the partially covered
  bin, so its five lines keep one common reference). The coverage fraction is now printed.

- 2026-09-08 — **Stage 0 (prerequisite rerun) verified ALREADY DONE, not repeated.** The
  2026-09-07 gap-cut change (`muon_gap_cuts_acceptance.md` F17/F18) had invalidated the whole
  trigger-efficiency chain, but a parallel session had rerun it overnight. Verified on disk:
  `mc_trig_eff_hists_pp24_full[_medium_wp][_step3].root` 2026-09-07 21:20-21:21,
  `single_mu_effcy_pT_fit_mc[_medium_wp].root` 21:20, the pp24 data tag-and-probe
  `single_mu_effcy_pT_fit[_medium_wp].root` 21:14-21:15 (keys now `..._2_00_TO_2_20_divided`), and
  the Step-3 dR fits for all methods x signs x modes x both WPs 2026-09-07 23:52 - 2026-09-08
  00:22. Axis check on the two modes this doc needs:
  `..._expo_os_nocorr_ptmerge.root` X(7) = 8/11.54/16.65/24.01/34.64/49.97/72.08/150, Y(9) =
  **-2.2 ... 2.2**; `..._expo_os_nocorr_etamerge_ptmerge.root` Y(3) = **0 / 1 / 2 / 2.2**. Both
  fresh on the new fiducial region. The closure ROOT files (2026-09-03) are the only stale link
  and are rebuilt by step 3 of the plan anyway.
  **Carried forward as a known, unreviewed input:** the 2026-09-07 `expo` shape restriction
  (A <= 0, p >= 1, Step 3 only) is baked into these fit files and has no `/review-analysis-code`
  log; it is the user's own in-flight change and is used as-is.

## Results & Observations

### R1 — The delivered numbers (pp24 fullsim FULL, Tight, `pair_trig_eff_pp24_full.root`)

Opposite sign, the three delivered pair-p_T cells. `eps^pair` = S1/S0 (replaces the whole weight),
`K` = S2/S0 (multiplies the two single-muon efficiencies). Errors are conditional/binomial.

| p_T^pair [GeV] | \|eta^pair\| | eps^pair (sig) | K (sig) | eps^pair (wide) | K (wide) | raw all / fired (sig) |
|---|---|---|---|---|---|---|
| 49.97-72.08 | 0-1 | 0.3393 ± 0.0100 | 0.4319 ± 0.0130 | 0.4111 ± 0.0083 | 0.5171 ± 0.0105 | 4830 / 1633 |
| 49.97-72.08 | 1-2 | 0.6233 ± 0.0150 | 0.7311 ± 0.0189 | 0.6772 ± 0.0113 | 0.7821 ± 0.0137 | 2142 / 1295 |
| 49.97-72.08 | 2-2.2 | 0.5686 ± 0.0418 | 0.6360 ± 0.0482 | 0.6268 ± 0.0313 | 0.6916 ± 0.0351 | 268 / 156 |
| 72.08-103.98 | 0-1 | 0.2210 ± 0.0148 | 0.2727 ± 0.0182 | 0.2679 ± 0.0133 | 0.3333 ± 0.0168 | 1634 / 395 |
| 72.08-103.98 | 1-2 | 0.5554 ± 0.0319 | 0.6351 ± 0.0365 | 0.6128 ± 0.0236 | 0.6952 ± 0.0267 | 619 / 322 |
| 72.08-103.98 | 2-2.2 | 0.4811 ± 0.0951 | 0.5313 ± 0.1082 | 0.5110 ± 0.0708 | 0.5519 ± 0.0782 | 65 / 31 |
| 103.98-150 | 0-1 | 0.2297 ± 0.0321 | 0.2859 ± 0.0415 | 0.2215 ± 0.0217 | 0.2753 ± 0.0279 | 320 / 71 |
| 103.98-150 | 1-2 | 0.4201 ± 0.0847 | 0.4791 ± 0.0971 | 0.4752 ± 0.0604 | 0.5574 ± 0.0750 | 102 / 36 |
| 103.98-150 | 2-2.2 | 0.8072 ± 0.3150 | 0.8288 ± 0.3195 | 0.8491 ± 0.1526 | 0.8934 ± 0.1616 | **4 / 2** |

Same sign, same cells (`sig` window): 0.3331 ± 0.0429 / 0.5399 ± 0.0841 / 0.6000 ± 0.2191 in
49.97-72.08; 0.3317 ± 0.1318 / 0.5000 ± 0.1768 / **no pairs** in 72.08-103.98; 0.2857 ± 0.1707 /
**no pairs** / **no pairs** in 103.98-150. Raw counts 176/85/5, 43/8/0, 7/0/0.

### R2 — The physics the numbers say

1. **The correction is enormous at high pair p_T, and it is the ΔR correlation.** `K` is the
   single-number analogue of eps_dR: it is 0.92-0.96 in the lowest pair-p_T cells (i.e. eps_dR -> 1,
   as it must) and falls to **0.43** at 50-72 GeV and **0.27** at 72-104 GeV in the barrel. A
   factor 2-4 of the trigger probability at high pair p_T is the close-by-RoI loss, not the
   single-muon turn-ons. Physically expected: high pair p_T means a collimated pair.
2. **The efficiency is strongly pair-eta dependent** at fixed pair p_T — 0.34 (barrel) vs 0.62
   (\|eta\| 1-2) at 50-72 GeV. Folding eta is therefore not free (see R4).
3. **The mass dependence is real and LARGE.** Within 1-4 GeV, splitting at the signal window:

| p_T^pair [GeV] | \|eta^pair\| | eps(1.08-2.9) | eps(outside, in 1-4) | eps(1-4) | wide/sig |
|---|---|---|---|---|---|
| 49.97-72.08 | 0-1 | 0.3393 | 0.5088 | 0.4111 | 1.21 |
| 49.97-72.08 | 1-2 | 0.6233 | 0.7481 | 0.6772 | 1.09 |
| 72.08-103.98 | 0-1 | 0.2210 | 0.3273 | 0.2679 | 1.21 |
| 72.08-103.98 | 1-2 | 0.5554 | 0.6707 | 0.6128 | 1.10 |
| 103.98-150 | 0-1 | 0.2297 | 0.2096 | 0.2215 | 0.96 |

   Higher mass at the same pair p_T means a wider opening angle and less close-by loss, exactly the
   sign observed. **Consequence: the 1-4 GeV window is NOT a usable proxy for the signal window** —
   using it would over-state the efficiency by ~10 % (endcap) to ~21 % (barrel), i.e. under-correct
   the signal yield by the same amount. A template fit over 1-4 GeV needs its efficiency
   differentially in mass, or the fit must be done in the signal window.
4. **Same-sign is not measurable above 72 GeV** (43, 8 and 0 pairs; two cells hold none at all).
   Since the background estimate is OS - SS and R25 showed the correction is genuinely
   charge-dependent, this is a hard limit of the present MC sample, not of the method.
5. **The most forward, highest-p_T cell (103.98-150 x \|eta\| 2-2.2) rests on 4 pairs**
   (eps = 0.807 ± 0.315) and must not be used as delivered.

### R3 — The closure comparison (`closure/single_value_highpt_comparison/`)

Above the 49.97 GeV cell edge, opposite-sign pairs in the single-b signal region, Tight, each
series against its own covered region:

| eta^pair panel | ΔR, p_T-merged (B) | ΔR, + \|eta\| fold (D) | single value, sig | single value, wide |
|---|---|---|---|---|
| [-2.2,-2.0) | 0.9865 | 0.5874 | 0.9385 | 0.8569 |
| [-2.0,-1.5) | 1.0588 | 0.7613 | 1.1087 | 1.0181 |
| [-1.5,-1.0) | 1.1062 | 0.6860 | 0.8656 | 0.7943 |
| [-1.0,-0.5) | 1.0469 | 0.4444 | 0.9425 | 0.7822 |
| [-0.5,0.5) | 0.9491 | 0.6803 | 1.0637 | 0.8814 |
| [0.5,1.0) | 0.9283 | 0.9416 | 0.9443 | 0.7863 |
| [1.0,1.5) | 0.9769 | 0.9555 | 0.8795 | 0.8070 |
| [1.5,2.0) | 1.0971 | 0.9942 | 1.0698 | 0.9821 |
| [2.0,2.2) | 0.9615 | 1.0002 | 1.0522 | 0.9585 |
| **ALL** | **1.0003** | **0.7433** *(bug, see R5)* | **1.0000** | **0.8607** |

(the calibrated K form gives 1.0000 / 0.8699 and a slightly tighter per-panel spread, 0.89-1.04.)

**Readings.**
- **The single-value procedure is competitive.** Its 1.0000 inclusive is by construction, but its
  per-panel spread (0.88-1.11) is the genuine test — the residual within-cell dependence — and it
  is no worse than approach B's (0.93-1.11), with 3 \|eta\| cells instead of 9 and no fit at all.
- **The wide window costs 14 %** inclusively, up to 22 % in the barrel panels — R2 item 3 seen in
  the closure.
- **★ THE APPROACH-D COLUMN ABOVE WAS NOT THE FOLD APPROACH — RETRACTED, see R5.** The first
  reading of this table was that D's poor closure is a *sign asymmetry* (0.44 at [-1,-0.5) against
  0.94 at its mirror). `/review-analysis-code` traced that pattern to a **bug**, not to physics:
  `DrCorrectionCascadeEvaluator::Eval` looked the pair-eta cell up with the SIGNED `pair_eta` on a
  FOLDED (0 -> 2.2) axis, so every negative-eta pair fell in the underflow and got
  `eps_dR = 1`, i.e. no correction at all. The perfectly sign-split pattern is that bug's
  signature. **Both the D column above and the parent doc's R4 are affected.** Fixed and re-run —
  see R5 for the corrected numbers.
- **A genuine C = 0, not a coverage artefact** (checked, so a future reader does not reopen it):
  panel eta^pair in [-1.5,-1.0), presentation bin [123.37, 150) GeV shows C = 0 for **all four**
  corrected series. That bin is FULLY covered -- its cell `[103.98,150) x |eta| 1-2` is delivered
  (102 raw pairs, eps = 0.4201) and `den_paireff = den = 4.33e-06` -- so it simply holds ~7
  effective pairs of which none fired 2mu4. The dR series show the identical zero.
- Coverage: the single-value cells reach **70.94 %** of the no-trigger yield above the first drawn
  presentation bin; the remaining 29 % is the bin straddling the 49.97 GeV cell edge plus the one
  cell the delivery gate refuses (R6), both omitted from the single-value series (doc §PP-4).


### R5 — ★ A CRITICAL BUG FOUND BY `/review-analysis-code`, FIXED, AND EVERYTHING IT TOUCHED RERUN

**The bug (pre-existing, not introduced here).** `DrCorrectionCascadeEvaluator::Eval` looked the
pair-eta cell up with the **signed** `pair_eta` on a **folded** (|eta|, 0 -> 2.2) Y axis. In every
`*_etamerge*` mode that sends every pair with `pair_eta < 0` into the underflow bin, where it is
counted "outside the cell grid" and returned **`eps_dR = 1` — no correction at all**. Half the
sample, silently: `628453 OUTSIDE the cell grid (50.67%)` in the etamerge closure log against
`0 OUTSIDE` for an un-folded mode. It is the SAME bug `mc_trigeff_dr_binning_approaches.md` D11
records and fixed in `DrCorrectionEvaluator::Eval` and in the closure's sample Filter — but the
cascade class had reimplemented the two lines and was missed, and the cascade class is the one the
closure's *delivered* series calls. `DrCorrectionCrossxEvaluator` carried it too, latently: harmless
only because `DrCorrCrossxMode()` is the un-folded `nocorr_ptmerge`, and it would have removed the
trigger correction from half the pp24 cross-section the moment a folded mode was adopted — which is
exactly the open decision in the parent doc's Remaining Work 1.

**The fix.** One shared helper, `DrCorrectionEvaluator::CellLookupEta(pair_eta, folded)` in
`dr_correction_apply.h`, called by all three classes (four call sites). A fold can no longer be
forgotten by the next class that needs a cell index. **Verified**: the etamerge closure now reports
`0 OUTSIDE the cell grid (0%)`.

**RETRACTION.** R3's headline "approach D's failure is a SIGN asymmetry" was **this bug's
signature, not physics**. Corrected numbers (Tight, above the cell edge, PURE form):

| eta^pair panel | ΔR, p_T-merged (B) | ΔR, + \|eta\| fold (D) | single value, sig | single value, wide |
|---|---|---|---|---|
| [-2.2,-2.0) | 0.9865 | 0.9321 | 0.9378 | 0.8555 |
| [-2.0,-1.5) | 1.0588 | 1.0107 | 1.1087 | 1.0181 |
| [-1.5,-1.0) | 1.1062 | 0.9011 | 0.8656 | 0.7943 |
| [-1.0,-0.5) | 1.0469 | 0.9155 | 0.9425 | 0.7822 |
| [-0.5,0.5) | 0.9491 | 0.9040 | 1.0637 | 0.8814 |
| [0.5,1.0) | 0.9283 | 0.9416 | 0.9443 | 0.7863 |
| [1.0,1.5) | 0.9769 | 0.9555 | 0.8795 | 0.8070 |
| [1.5,2.0) | 1.0971 | 0.9942 | 1.0698 | 0.9821 |
| [2.0,2.2) | 0.9615 | 1.0002 | 1.0519 | 0.9581 |
| **ALL** | **1.0003** | **0.9362** | **1.0000** | **0.8607** |

Approach D closes at **0.9362**, not 0.7433, and its per-panel values are 0.90-1.01 on the negative
side against 0.94-1.00 on the positive side — **no sign asymmetry**. D is still the worst of the
three corrected series above 50 GeV (a ~6 % under-correction), but that is now a statement about
the fold's coarseness, not about a broken lookup. **The conclusions of this doc about the two
single-value series and about approach B are UNCHANGED** — the bug touched only the folded modes.

**Blast radius, all rerun (2026-09-08):** the closures for both `*_etamerge*` modes and both WPs,
this doc's comparison figure at both WPs, and — because the parent doc's four-approach comparison
reads the same corrupted files — `mc_trigeff_dr_binning_approaches.md`'s own closure set and
chi^2/ndof figures were regenerated for all four modes at both working points. **That doc's R4
conclusion ("the new fold's chi^2/ndof is 4-8x worse ... argues AGAINST adopting the fold") was
measured with this bug in place and must be re-judged by the user against the regenerated figures.**

### R6 — The delivery gate (added in the same review round)

`PairTrigEffEvaluator` now refuses a cell with fewer than `PairTrigEff::MinCellPairs() = 50` raw
pairs or a value outside (0, 1]. Nothing is removed from the ROOT file — every cell is still
measured, written and printed — only automatic delivery through `Eval()` is gated, and the fill
macro prints every refusal with its reason. Tight, refused:

- **os**, `[103.98, 150) x |eta| [2.0, 2.2)`: eps = 0.8072 ± 0.3150 on **4** raw pairs (`sig`),
  0.8491 ± 0.1526 on 14 (`wide`). This was the one cell that read *higher* than its row-mates when
  the physics requires it to be the lowest — a 2-pair fluctuation.
- **ss**, in the `sig` window: every cell above 72 GeV plus `[49.97,72.08) x |eta| [2.0,2.2)` --
  43, 8, 0, 7, 0, 0 and 5 raw pairs. In the WIDER `wide` window one of them survives,
  `[72.08,103.98) x |eta| < 1`, which clears 50 pairs there. **The same-sign single-value efficiency
  is essentially not measurable above ~72 GeV in this MC sample**, and not at all in the signal
  window.

The threshold is a NEW named choice (50 raw pairs ≈ ±0.07 on a binomial at eps ~ 0.5) and is
flagged for the user in R4.

### R4 — What is NOT settled by this work

- Which of `eps^pair` and `K` the cross-section should apply (§2 / D1). `K` keeps the data
  tag-and-probe calibration; the pure `eps^pair` would need a separate data/MC scale factor, and
  the MC/data pair over-efficiency is measured at `<r_1 r_2> = 1.27`.
- Whether to adopt the single-value procedure at all, and in which cells. Nothing is wired into
  the cross-section.
- The Medium working point HAS been run and is consistent (Progress Log).
- **The delivery-gate threshold `MinCellPairs() = 50`** (R6) is a choice this session made to stop a
  4-pair cell reaching a consumer through the same API as a 4830-pair one. It removes one
  opposite-sign cell and most same-sign cells from automatic delivery. If you want a different
  threshold, or none, it is one constant in `Utilities/PairTrigEffEvaluator.h`.
- **Parent doc `mc_trigeff_dr_binning_approaches.md` R4 must be re-judged** against the regenerated
  four-approach figures (R5).


## Remaining Work

Nothing in this doc's own scope. What is left is USER JUDGEMENT (R4) plus one hand-off:

1. **Which form to apply** — the pure `eps_2mu4^pair` (the procedure as stated, but pure MC, so it
   carries the 1.27 MC/data pair over-efficiency and needs a separate scale factor) or the
   calibrated `K` (multiplies the DATA tag-and-probe singles, the analysis's existing calibration).
   Both are in the ROOT file.
2. **Whether to adopt the single-value procedure at all, and in which cells.** Nothing is wired into
   the cross-section.
3. **The delivery threshold** `MinCellPairs() = 50` (R6) — this session's choice, one constant.
4. **HAND-OFF to `mc_trigeff_dr_binning_approaches.md`:** its R4 is retracted and the four-approach
   ranking reverses, so the choice of dR cell grouping for the pp24 cross-section is reopened on the
   regenerated figures. That decision is that doc's, not this one's.

### R7 — CSV tables (user request, 2026-09-08)

`plotting_codes/trig_effcy/mc_based/write_pair_trig_eff_tables.cxx` reads the deliverable and
writes three CSVs per working point into
`plots/pp_trigger_efficiency/mc_based[_medium]/single_value_pair_eff_tables/`, all with pair-p_T
ROWS and |eta^pair| COLUMNS on the canonical axes:

| file | contents |
|---|---|
| `single_value_pair_eff_opposite_sign[_medium_wp].csv` | eps^pair and K with their conditional errors, BOTH mass windows, plus a per-cell delivery status |
| `single_value_pair_eff_same_sign[_medium_wp].csv` | the same for same-sign pairs |
| `single_value_pair_stats_same_sign[_medium_wp].csv` | the same-sign RAW pair counts in the SIGNAL window (`n_all`, `n_2mu4`), with per-row and per-column totals |

Nothing is recomputed: every value is read from `pair_trig_eff_*.root` as written, and the status
column is obtained by ASKING `PairTrigEffEvaluator` at each cell centre rather than re-implementing
its gate, so a cell marked `delivered` in a CSV is exactly one a consumer can `Eval()`.
Cross-checked against R1 and against the fill log's totals (same-sign 32001 selected / 19177 firing).

The same-sign signal-window statistics make R2 item 4 concrete — 17255 / 13724 / 1022 pairs in the
three |eta^pair| groups summed over all pair p_T, but only **176 / 85 / 5** in 49.97-72.08 GeV,
**43 / 8 / 0** in 72.08-103.98 and **7 / 0 / 0** in 103.98-150.

## Latest Stage

**2026-09-08 — DONE.** Both reviews PASS, everything is committed, and nothing is in flight. The
deliverables are `pair_trig_eff_pp24_full[_medium_wp].root` (+ `Utilities/PairTrigEffEvaluator.h` to
read it) and
`plots/pp_trigger_efficiency/mc_based[_medium]/closure/single_value_highpt_comparison/closure_highpt_single_value_{pure,calibrated}.png`.
The cross-section is deliberately unchanged. Everything further is the user decisions in Remaining
Work.
