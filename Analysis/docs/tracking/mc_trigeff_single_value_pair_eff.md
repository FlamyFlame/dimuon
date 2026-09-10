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

### D4: the data/MC difference is handled by a SCALE FACTOR, not by re-factorizing the weight (user, 2026-09-08)
**Old (D1):** two applied forms were delivered and plotted — the pure `eps_2mu4^pair`, and a
"calibrated" `K = eps^pair / <eps_1 eps_2>` meant to be applied as `eps(1) eps(2) K` on top of the
DATA tag-and-probe singles, so that MC was trusted only for the deviation from factorization.
**New:** the applied object is the **raw** `eps_2mu4^pair`. The data/MC difference will be corrected
by the **product of the two single-muon data/MC scale factors** — a step still under discussion with
colleagues and explicitly a future to-do, NOT part of this doc.
**Reason (user):** `K` re-introduces the factorization the single-value procedure exists to avoid —
writing the weight as `eps(1) eps(2) K` is structurally the old procedure with the dR fit replaced
by a per-cell constant, so it is a third procedure rather than "a single pair-level efficiency".
And the correction is a DATA-application question: any data/MC factor **cancels identically in an
MC closure**, so it cannot belong in this doc's test and would only obscure what the test measures.
**Consequences:** the closure applies the raw value only; the `paireffK_*` closure series and the
`_calibrated` figures are withdrawn (one PNG per comparison per WP now, not two). `K` itself is
still written to `pair_trig_eff_*.root` and to the CSVs as a **DIAGNOSTIC** — the ratio of the
measured pair efficiency to what the factorized product predicts, i.e. how much of the pair
inefficiency is single-muon turn-on and how much is the close-by correlation (0.96 at 8-11 GeV
falling to 0.27 at 72-104 GeV, barrel).

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
  **★ CARRIED FORWARD AS A SECOND UNCOMMITTED, UNREVIEWED INPUT — and the one that defines this
  doc's |eta^pair| axis** (`/review-analysis-code` 2026-09-08, issue 1). `PairTrigEff::AbsEtaGroups()`
  gets the three |eta^pair| cells by calling `MakeDrEtaGroups(..., true)`. That returns the
  sign-independent fold {0, 1, 2, 2.2} only in the version of
  `plotting_codes/trig_effcy/mc_based/dr_correction_cell_groups.h` carrying
  `mc_trigeff_dr_binning_approaches.md` D11 — which is a CONCURRENT session's **uncommitted** work
  (`M` in `git status`, its own review log reading `Iterations completed: 0`). At HEAD the same call
  returns the SIGNED grouping {-2.2, -1.0, 1.0, 2.2}, and since every histogram here is filled with
  `|pair_eta|`, a clean checkout of HEAD would leave the first cell permanently empty, collapse
  |eta| in [1, 2.2) into one cell, and label the cells "-2.2..-1" — a different measurement, with no
  crash. **A clean checkout of HEAD therefore does NOT reproduce the numbers in this doc.**
  Mitigated 2026-09-08 by a guard in `AbsEtaGroups()` that THROWS unless the returned grouping is
  actually a |eta| fold (folded flag set, first edge 0), so the failure is loud rather than silent;
  the dependency itself is the concurrent thread's to commit. §PP-1's "the SAME helper as the
  `nocorr_etamerge*` dR modes" should be read with that caveat.

  **Carried forward as a known, unreviewed input:** the 2026-09-07 `expo` shape restriction
  (A <= 0, p >= 1, Step 3 only) is baked into these fit files and has no `/review-analysis-code`
  log; it is the user's own in-flight change and is used as-is.

- 2026-09-09 — **★ EVERY NUMBER IN THIS DOC IS NOW SUPERSEDED BY A SELECTION CHANGE MADE BY A
  CONCURRENT THREAD. Nothing here was wrong; its sample no longer exists.** State established at
  the start of this session by `git status` / `git log` / file mtimes, before any action:

  **(a) The commit question is SETTLED — nothing of this thread is uncommitted.** All of
  `Utilities/PairTrigEffEvaluator.h`, `RDFBasedHistFilling/FillMCTrigEffPairEff.cxx`,
  `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_closure_highpt_compare.cxx`,
  `pipelines/run_mc_trigeff_pair_eff.sh`, this doc — and, crucially,
  `plotting_codes/trig_effcy/mc_based/dr_correction_cell_groups.h` — are committed
  (`3a97a3f`, plus `660107a` for `write_pair_trig_eff_tables.cxx`). The Stage-0 hazard recorded
  above ("a clean checkout of HEAD does NOT reproduce the numbers in this doc") is therefore
  **RESOLVED**: the D11 |eta| fold is at HEAD.

  **(b) Both guards verified by running them, not by reading them.**
  `PairTrigEff::AbsEtaGroups()` returns `folded=1, edges 0 / 1 / 2 / 2.2` — the fold, as required.
  `PairTrigEffEvaluator::CheckCanonicalBinning()` now REFUSES the delivered file:
  `pair-pT edge 0 of pair_trig_eff_pp24_full.root is 8.000000 but cell mode 'nomerge' of ParamsSet
  says 9.000000 -- stale file`. The guard is doing exactly the job it was written for.

  **(c) WHY it refuses: the analysis moved underneath this measurement.**
  `mu_pt45_gap125_pairpt9_adoption.md` (ACTIVE, concurrent) adopted three user decisions —
  **muon pT 4.0 → 4.5 GeV**, gap window **(−1.30,−1.05) → (−1.25,−1.05)**, signal **pair pT 8 → 9
  GeV** — and with them `ParamsSet::pair_pt_coarse_bins` = 8 log bins **9 → 150**:
  9 / 12.79 / 18.18 / 25.85 / 36.74 / **52.23 / 74.24 / 105.53** / 150. The three cells this doc
  delivers are no longer [49.97,72.08) / [72.08,103.98) / [103.98,150) but
  **[52.23,74.24) / [74.24,105.53) / [105.53,150)**. That thread also closed a genuine gap in
  `MCTrigEffPairSel::Step3PairSelection()` — the pp **same-vertex** pair requirement (`pair_pass_*`),
  which Steps 2/3/4, the closure AND `FillMCTrigEffPairEff.cxx` had all been missing — so the
  measured population changes as well as the axis.

  **(d) The upstream has ALREADY been rerun on the new selection; this thread's outputs have not.**
  On disk in `pythia_fullsim_full_sample/`: the MC pair file
  `muon_pairs_..._mc_trig_full.root` **2026-09-08 17:11**, `mc_trig_eff_hists_pp24_full.root` 18:46,
  `single_mu_effcy_pT_fit_mc.root` 18:47, `..._step3.root` 18:48, the Step-3 ΔR fits 18:51 and a
  further `polyu_fixedRp` sweep 23:52 → 2026-09-09 00:02. Their axes read
  `X: 9 12.79 18.18 25.85 36.74 52.23 74.24 150` — the NEW binning. Against that,
  `pair_trig_eff_pp24_full[_medium_wp].root` (12:35 / 12:39), the four
  `mc_trig_eff_closure_pp24_full*` files (12:37–12:41) and every figure and CSV under
  `closure/single_value_highpt_comparison/` (12:46) are **older than their own input** and carry the
  8 → 150 axis.

  **Consequence:** R1, R2 (item 3's mass table), R3, R5's corrected table, R6's refusal list, R7's
  CSVs and R8's merged cells are all measurements of a sample the analysis no longer uses. They are
  retained as the record of the OLD selection; **none of them may be quoted as current**. The
  procedure, the code, both guards and the physics readings (the ΔR correlation dominates at high
  pair pT; the 1–4 GeV window under-corrects the signal region; the pT merge buys reach cheaply)
  are unaffected in KIND — but every NUMBER needs re-measuring, on ~23 % fewer pp24 pairs.

  **Not started, and deliberately not started:** re-running `run_mc_trigeff_pair_eff.sh` writes
  `pair_trig_eff_*.root`, `mc_trig_eff_closure_*.root` and the plot tree that the concurrent thread
  is itself still writing (its last artefact landed 00:02 this morning) — parallel writers to one
  git-invisible output path clobber each other silently (.claude/CLAUDE.md §Parallel delegation,
  item 2). Held for the user's decision; see Latest Stage.


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


### R8 — The pair-pT-MERGED variant, and the two comparison figures (user, 2026-09-08)

A second way of buying statistics in the top cells, on the pair-pT axis instead of the mass axis:
`ptmerge` combines the last two coarse cells into one **[72.08, 150) GeV** cell, the direct analogue
of the dR correction's own `nocorr_ptmerge` (parent R32). It is **not a new binning** — the merged
cell is the two source cells' num / den / A / B summed BEFORE the ratio — and it is opt-in and
suffixed (`_ptmerge`). Measured for both mass windows (the merge is a projection of sums already
there, so producing `wide` costs nothing); the request's scope, and everything below, is the
SIGNAL window.

**Delivered, opposite sign, Tight, signal window:**

| p_T^pair [GeV] | \|eta^pair\| | eps^pair | K | raw all / fired |
|---|---|---|---|---|
| 49.97-72.08 | 0-1 | 0.3393 ± 0.0100 | 0.4319 ± 0.0130 | 4830 / 1633 |
| 49.97-72.08 | 1-2 | 0.6233 ± 0.0150 | 0.7311 ± 0.0189 | 2142 / 1295 |
| 49.97-72.08 | 2-2.2 | 0.5686 ± 0.0418 | 0.6360 ± 0.0482 | 268 / 156 |
| **72.08-150** | 0-1 | 0.2221 ± 0.0136 | 0.2744 ± 0.0167 | 1954 / 466 |
| **72.08-150** | 1-2 | 0.5400 ± 0.0300 | 0.6173 ± 0.0343 | 721 / 358 |
| **72.08-150** | 2-2.2 | 0.5051 ± 0.0938 | 0.5532 ± 0.1050 | **69 / 33** |

**★ The merge is cheap PER PANEL, and it buys the cell the un-merged version could not deliver.**
- **Panel-integrated closure is unchanged**: above the cell edge the merged series matches the
  un-merged one to **0.0054 at worst** across the nine pair-eta panels (0.9378 -> 0.9432 in the most
  forward-negative one; five panels agree to <= 0.0019), and both close to 1.0000 inclusively.
- **⚠ BUT THAT NUMBER IS PARTLY PROTECTED BY CONSTRUCTION, and it is NOT a per-bin statement**
  (`/review-plot` iteration 3, INFO 3). `eps_merged` is the S0-weighted mean of its two source
  cells, so the corrected yield summed OVER the merged cell is preserved exactly; only the
  redistribution across (pair-eta panel x fine pair-pT bin) can move a panel integral, and just
  **5.2 %** of the yield above 46.4 GeV sits beyond 83.46 GeV. So <= 0.0054 per panel is what the
  construction predicts, not independent evidence. **Per presentation BIN the merge moves the
  corrected yield by up to ~22 %** -- panel [-1.5,-1.0), bin [101.47,123.37): 1.718 -> 1.337; panel
  [-2.0,-1.5), same bin: 0.678 -> 0.555; panel [1.0,1.5), same bin: 0.399 -> 0.326. Quote "the merge
  is nearly free" only at the cell/panel level, never per bin.
- What the merge genuinely trades is RESOLUTION for REACH: it replaces two measured numbers (0.221
  and 0.230 in the barrel, 0.555/0.420, 0.481/0.807) by one, and in exchange every cell is
  measurable.
- **Coverage rises** 70.9413 % -> 70.9702 % of the no-trigger yield above the first drawn bin, and
  the merged series loses only **9** presentation bins to masking against the un-merged series' 11.
  (The figure originally printed their UNION; since 2026-09-08 it reports the two counts separately
  wherever they differ, because that difference IS the comparison the figure is about.)
- **No opposite-sign cell is refused any more.** The forward `[103.98,150) x |eta| 2-2.2` cell that
  rested on 4 raw pairs (R6) becomes part of a 69-pair cell and clears the gate; the un-merged
  version delivers 8 of 9 cells, the merged one 6 of 6.
- **Same sign is helped but not rescued**: `[72.08,150) x |eta| < 1` reaches exactly 50 pairs and is
  delivered (0.3260 ± 0.1178); `|eta| 1-2` still has 8 and `|eta| 2-2.2` none. Refusals fall from 7
  cells to 3, but the 35 % error on the one cell that survives says what it is worth.

**Two figures now, in their own subdirectories under `closure/single_value_highpt_comparison/`:**

| subdirectory | question | series |
|---|---|---|
| `mass_window_compr/` | does the mass window matter? | no-trigger + dR pT-merged + dR +\|eta\| fold + single value 1.08-2.9 + single value 1-4 |
| `pt_merge_compr/` | does merging the top two pair-pT cells cost anything? | the same two dR series + single value signal window on 8 cells + on the merged 7 |

each in the PURE and the CALIBRATED applied form, at both working points. The previous flat-path
PNGs were moved into `mass_window_compr/` (verified replaced, then removed).

**Reading, for the decision in R4:** the merge is the cheaper of the two ways of coping with the
statistics. It halves the numbers needed above 72 GeV, removes every opposite-sign refusal, and
costs <= 0.005 in the panel-integrated closure (with the construction caveat above) — whereas
widening the mass window costs **14 % inclusively, a genuine bias that no construction protects**
(R3). If the single-value procedure is adopted, adopt it merged and in the signal window.

### R9 — Why all four corrected series track each other in the closure figure (user question, 2026-09-09)

**Not a bug.** Verified first that the applied weight is what D4 says: `FillMCTrigEffClosure.cxx`
loads the evaluator with `ApplyForm::kPure` for both windows and both cell modes — the raw
`eps_2mu4^pair`, no K, no scale factor — and the comparison macro draws only those. (The dead
`paireffK_*` histograms are still in the closure ROOT files; nothing reads them.)

**The similarity is by CONSTRUCTION, and the construction covers the dR series too.** Both
corrections are cell-integral-preserving on the very cell they were measured in:

- single value: `sum_fired w / eps^pair = S1 / (S1/S0) = S0` — exactly.
- ΔR procedure: `eps_dR` is the SAME estimator resolved in ΔR instead of integrated over it —
  per ΔR bin j of the cell, `eps_dR,j = [sum_fired w/(eps1 eps2)]_j / [sum_all w]_j`, so
  `sum_fired w/(eps1 eps2 eps_dR) = sum_j S0_j = S0` — exactly, up to the fit smoothing the raw
  ratio. Indeed `K = S2/S0` IS the ΔR-integrated `eps_dR`.

So the two procedures are the same ratio at two resolutions, both pinned to the same cell integral.
They can differ ONLY through within-cell shape — and above 50 GeV both legs sit on the single-muon
plateau (R2 item 1: K falls to 0.27 while the singles' turn-on contributes almost nothing), so the
leg-p_T differential information the ΔR form carries is nearly constant there. **The doc's §PP-4
"honest statement 1" is therefore INCOMPLETE: it says the single-value series closes to 1 by
construction over its own cell, but the same is true of the ΔR series. The figure is much less of a
head-to-head test than it reads as.**

**On top of that, the excursions from 1 are COMMON-MODE, because all four numerators are sums over
the SAME fired pairs.** Only the weight `1/eps_hat` differs; which pairs fired is one shared random
draw. Measured over the 40 populated (pair-eta panel x fine p_T bin) cells above 50 GeV, the
correlation of `(C - 1)` between series is **+0.82** (ΔR p_T-merged vs single value), **+0.94**
(ΔR fold vs single value), **+0.85** (the two ΔR modes). All-panels-summed, per fine bin:

| p_T^pair bin | N_eff | eps | C_B (ΔR) | C_D (fold) | C_sv sig | C_sv wide | sigma_stat |
|---|---|---|---|---|---|---|---|
| 56.46-68.65 | 1491 | 0.419 | 0.9627 | 0.9019 | 0.9504 | 0.8171 | 0.030 |
| 68.65-83.46 | 814 | 0.374 | 1.0208 | 0.9211 | 1.0143 | 0.8689 | 0.045 |
| 83.46-101.47 | 413 | 0.260 | 0.9050 | 0.8014 | 0.8328 | 0.7069 | 0.083 |
| 101.47-123.37 | 166 | 0.325 | 1.3376 | 1.1136 | 1.1495 | 1.1044 | 0.112 |
| 123.37-150.00 | 104 | 0.194 | 0.8193 | 0.7160 | 0.6877 | 0.6844 | 0.200 |

(`N_eff = (sum w)^2 / sum w^2` of the no-trigger denominator; `sigma_stat = sqrt((1-eps)/(eps N_eff))`
is the shared binomial scale. All four series move up together in 101-123 and down together in
83-101 and 123-150, each excursion ~1-1.3 sigma_stat.)

**Consequence for how the figure must be read.** `|C - 1|` is NOT a measure of procedure quality
here: it is dominated by the shared draw, and it is largest exactly where N_eff collapses. The
procedure-dependent information is the RATIO BETWEEN series, in which the shared fluctuation
cancels — and that ratio is not flat:

| p_T^pair bin | C_B / C_sv | C_D / C_sv | C_wide / C_sv |
|---|---|---|---|
| 56.46-68.65 | 1.0129 | 0.9490 | 0.8597 |
| 68.65-83.46 | 1.0064 | 0.9081 | 0.8567 |
| 83.46-101.47 | 1.0866 | 0.9622 | 0.8488 |
| 101.47-123.37 | 1.1637 | 0.9688 | 0.9608 |
| 123.37-150.00 | 1.1915 | 1.0412 | 0.9953 |

**This is the real result the figure carries.** The ΔR procedure and the single value agree to
**~1 %** at 56-83 GeV, where the ΔR fit has statistics — and the ΔR procedure runs **9 / 16 / 19 %
HIGHER** above 83 GeV, i.e. it corrects progressively less, exactly in the region where its fit is
starved and which is the reason this doc exists. The `wide` column is the mass-window bias
(~14-15 %, flat where statistics are good), independent of the shared fluctuation as R2 item 3 said.

**Caveat:** computed on the 2026-09-08 closure files, i.e. the SUPERSEDED selection (see the
2026-09-09 Progress Log entry). The inclusive numbers here start at the 56.46 GeV presentation-bin
edge rather than the coarse cell edge 49.97, so they do not reproduce R5's exact-1.0000 column;
the per-bin structure and the correlations are what the reading rests on and are unaffected.

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

1. ~~Which form to apply~~ — **SETTLED (user, 2026-09-08, D4):** the raw `eps_2mu4^pair`, with the
   data/MC difference to be corrected later by the product of the two single-muon data/MC scale
   factors. That step is under discussion with colleagues and is a future to-do; it does not enter
   the MC closure, where any such factor cancels identically. `K` remains in the file as a
   diagnostic only.
2. **Whether to adopt the single-value procedure at all, and in which cells** — and, if so, whether
   on the merged pair-pT cells, which R8 shows to be nearly free. Nothing is wired into the
   cross-section.
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

Each of those exists for BOTH cell modes (the merged one carries a `_ptmerge` suffix and seven
pair-p_T rows). **Four COMPACT tables also sit beside the figure they belong to**, in
`closure/single_value_highpt_comparison/pt_merge_compr/` — the DELIVERED merged cells only, as a
strict **2 pair-p_T rows x 3 |eta^pair| columns** matrix, signal window (user, 2026-09-08):

| file | blocks, each 2 x 3 |
|---|---|
| `pt_merge_pair_eff_{opposite,same}_sign[_medium_wp].csv` | `eps`, `eps_err`, `K`, `K_err`, `status` |
| `pt_merge_pair_stats_{opposite,same}_sign[_medium_wp].csv` | `n_all`, `n_2mu4` |

Quantities needing more than one number per cell are stacked as further blocks of the SAME shape
rather than widened into extra columns, so every block reads as the matrix it is.

Nothing is recomputed: every value is read from `pair_trig_eff_*.root` as written, and the status
column is obtained by ASKING `PairTrigEffEvaluator` at each cell centre rather than re-implementing
its gate, so a cell marked `delivered` in a CSV is exactly one a consumer can `Eval()`.
Cross-checked against R1 and against the fill log's totals (same-sign 32001 selected / 19177 firing).

The same-sign signal-window statistics make R2 item 4 concrete — 17255 / 13724 / 1022 pairs in the
three |eta^pair| groups summed over all pair p_T, but only **176 / 85 / 5** in 49.97-72.08 GeV,
**43 / 8 / 0** in 72.08-103.98 and **7 / 0 / 0** in 103.98-150.

## Latest Stage

**2026-09-09 — BLOCKED ON A USER DECISION; nothing in flight, nothing running.**

Settled this session: everything is committed (Progress Log (a)); both guards verified by
execution (b); the `dr_correction_cell_groups.h` dependency that made a clean checkout
irreproducible is resolved.

Open, and the reason work stopped: the concurrent `mu_pt45_gap125_pairpt9_adoption.md` thread moved
the muon pT cut to 4.5 GeV, the gap window to (−1.25,−1.05), the signal pair-pT cut to 9 GeV and the
coarse pair-pT axis to 9 → 150, and rebuilt the whole upstream on 2026-09-08 evening. **Every number
this doc delivers is therefore stale** (Progress Log (c)/(d)); the canonical-binning guard already
refuses the delivered ROOT file. Re-measuring is a re-run of all four stages into an output tree
that the other thread is still writing, so it is the user's call, not this doc's.

Also still open and unchanged from 2026-09-08:
- `/review-analysis-code` on the pT-merge increment returned FAIL (2 WARNINGs: the then-uncommitted
  fold dependency — now resolved by (a) — and a stale file header in the plot macro); all four fixes
  were applied but NOT re-reviewed.
- `/review-plot` passed at iteration 3, but the D4 calibrated-form removal and the "N bins omitted"
  label removal landed after that pass, so the current figures are unreviewed. Both reviews are held
  until it is decided whether the figures are regenerated first — reviewing a superseded figure set
  would have to be repeated.
- The closure ROOT files still carry dead `paireffK_*` numerators (replotted without refilling); a
  refill drops them, and a refill is now a full re-measurement rather than a cleanup.
- The user decisions of R4 / Remaining Work 2–4 (adopt the procedure at all and in which cells;
  `MinCellPairs() = 50`; and, in the parent doc, which ΔR cell grouping the cross-section uses now
  that R4's ranking reversed).
