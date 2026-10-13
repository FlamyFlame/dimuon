# 3D bin-by-bin reco correction in (pair pT, η^pair, m_μμ): reco-binned ε + pair-η smoothness test

Mode: **IMPLEMENTATION**. Opened 2026-10-11. Successor of `pair_bin_migration_bin_by_bin_feasibility.md`
(its "awaiting user decision" is resolved here).

## Objective
User request 2026-10-09/11:
1. Record the analysis decision: (i) **3D bin-by-bin correction, not full Bayesian unfolding**, justified
   by the high bin-by-bin efficiency and purity of the migration study; (ii) the 3D variables are
   **(pair pT, η^pair, m_μμ)** — m_μμ REPLACES ΔR.
2. Reco efficiency plots w.r.t. **RECONSTRUCTED** (pair pT, η^pair, m_μμ) — all-OS category only, new
   plots (truth-binned plots untouched).
3. Pair-η smoothness test: correct the pp24 data pair-η distribution with the trigger and the 3D
   reco-binned ε (coarse η^pair axis, exactly as it would be applied) and check the corrected
   distribution — filled in FINE η^pair bins — is smooth across the coarse-cell boundaries (some coarse
   cells contain gap windows, others not).
4. IntNote: bin-by-bin decision (supporting criteria: bin-by-bin efficiency + purity) and the reco-eff
   procedure / results updated to gap cuts + new 3D variables, vs-reco-quantity plots only.
OUT (user): propagation to the crossx — NEXT-STEP TO-DO (two other sessions are hunting bad runs).
5. **User decision 2026-10-13 (supersedes D3 / option (b)):** the ε^reco numerator = reco pairs passing the reco signal
   cuts that are truth-matched AND whose truth passes the truth cuts (truth in T). Truth-matched reco pairs whose truth
   fails the truth cuts (outside the truth fiducial) are FAKES for the template fit, together with the unmatched reco
   pairs. Update the reco efficiency (bin-by-bin) accordingly. OUT (user): crossx propagation — the Pb+Pb 2026 ntuple
   processing must be rerun first (another session's task).

## Physics Procedure

### 1. Motivation
The measured yield is binned in RECO pair kinematics; the result is quoted for the TRUTH population
T (truth signal selections). The migration study (`pair_bin_migration_bin_by_bin_feasibility.md` R1–R3)
found every 3D (pT, η^pair, m_μμ) cell with N_eff ≥ 100 at purity ≥ 0.82 and bin-by-bin efficiency
≥ 0.81 (N[T]-weighted means 0.92), response near-diagonal → a bin-by-bin correction suffices. m_μμ
replaces ΔR: 285/288 cells populated (ΔR: 144/288, because ΔR ≳ 2m/pT is kinematically tied to pT),
ε varies only ≈4–5 % across m bins (smooth), and the out-of-T feed-in is ≈2× smaller (max o 0.090 vs
0.18 in cells N_eff ≥ 100).

### 2. Top-level equation
Per reco 3D cell j = (pT, η^pair, m_μμ) on the canonical axes (`ParamsSet::pair_pt_coarse_bins` ×
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` × `ParamsSet::SignalMinvRecoEffcyBins()`):

    N^T_j = Σ_{data pairs in reco cell j} w_trig / ε^reco_j          (then Σ over m → (pT, η^pair))
    ε^reco_j = N[R_j] / N[T_j]          (= 1/C_j, the bin-by-bin correction factor)

- **T_j**: truth OS pair (any origin), truth signal selections (truth mass window, truth pair pT > 9;
  NTuple truth fiducial muon pT > 4.5, |η| < 2.4), truth (pT, η^pair, m) in cell j. NO truth gap cut.
- **R_j**: truth-matched OS pair (both muons reco-matched) passing WP + the DATA reco signal selection
  incl. gap cuts (q·η windows on both muons, |η^pair| < 2.2), RECO (pT, η^pair, m) in cell j, AND its truth
  passes the truth signal selection (truth in T, any truth cell — migration within T is kept). N[R_j] = Σ_i M_ij.
  User decision 2026-10-13 (D7, supersedes option (b) of 2026-10-11).
- **F_j** (out-of-T feed-in): truth-matched pairs of reco cell j whose truth FAILS the truth signal selection.
  NOT in the numerator: they are FAKES of the template fit, together with the reco pairs that are not truth-matched
  (D7). The corrected yield is therefore N^T_j = (N^data_j − fakes_j incl. F_j)·w_trig / ε^reco_j; MC closure:
  N[R_j]·N[T_j]/N[R_j] = N[T_j]. Until the template fit carries F_j, a correction of the raw data yield returns
  N[T_j]/(1 − o_j), o_j = N[F_j]/(N[R_j] + N[F_j]).
- Output = efficiency-corrected truth distribution | truth signal selections = exactly the numerator
  population of A_sig(pT, η^pair) (`signal_acceptance_mass_window.md`), applied afterwards.
- MC-event weighted; WP Tight nominal, Medium for the WP systematic. Errors: N[R_j] = A + B,
  N[T_j] = A + D with A = N[T_j ∧ R_j], B = N[R_j] − A, D = N[T_j] − A disjoint → linear propagation with
  each piece's Σw².

### 3. Step-by-step method
(a) **ε^reco product** from the migration-study file `reco_eff/pair_bin_migration_pp24_full.root`
(cell index (pT, η^pair, m): `h_cell_minv_T`, `h_cell_minv_R_<wp>`, `h_cell_minv_M_<wp>`), which was
built from the ntuple-processing output with the RDF class's own expression builders and closes
bit-identically on the applied ε_reco legs. 3D, m-integrated 2D (pT, η^pair) and projections are ratios
of summed N[R]/N[T] (never projections of a ratio).
(b) **Plots** vs reco pT / η^pair / m (all OS, Tight + Medium), with the truth-binned ε (same axes) overlaid
for comparison.
(c) **Smoothness test.** pp24 data, the crossx selection and trigger weight EXACTLY (`RDFBasedHistFillingPP`
`signal_cuts`, Tight WP, OS, `AddPairEfficiencyWeightColumns` w_trig), reco weight = 1/ε^reco evaluated
at the pair's reco (pT, η^pair, m) on the COARSE η^pair axis. Data η^pair filled in a FINE, plot-only
binning (uniform 0.1 in [−2.2, 2.2], 44 bins; user 2026-10-11): raw, trigger-corrected, trigger+reco-
corrected. MC closure on the same fine binning: MC R (reco η^pair) corrected by ε^reco vs MC T (truth
η^pair) — isolates the coarse-η averaging effect from data statistics.

(d) **Within-cell η^pair shape systematic** (option A, user 2026-10-13). Keep the coarse cells; vary the MC truth η^pair
shape inside each (pT, coarse η) cell so that MC reco reproduces the trigger-corrected data reco shape in the fine 0.1 bins k,
and recompute the cell ε: ε_k = N[R_k]/N[T_k] (m-integrated), ε_c = Σ_k N[R_k]/Σ_k N[T_k], ε'_c = Σ_k D_k / Σ_k (D_k/ε_k)
(truth weight w_k ∝ (D_k/ε_k)/N[T_k] ⇒ MC reco = D), δ_c = ε_c/ε'_c − 1 = corrected-yield shift; systematic |δ_c|,
symmetric. Stat error by toys (data Poisson on the raw pair count of each fine bin × its mean weight; ε_k Gaussian-
binomial with N_eff of T_k; ε_c recomputed per toy from the same ε_k — fully correlated). δ uses the m-integrated
(pT-summed for merged groups) ε, an approximation to the applied m-resolved correction; D_k includes the data background. pT bins ≥ 36.74 GeV share the δ of their
merged range (fine gap-hole bins have N[R_k] ≈ 0 per bin). Evaluated from the smoothness histograms (no refill).

### 4. Negative constraints
1. ε^reco is evaluated on the CANONICAL coarse η^pair axis — never a finer axis for the test (the test
   asks whether the applied binning is good enough). The fine 0.1 binning is ONLY for the data/MC η^pair
   histograms.
2. Do NOT include out-of-T feed-in F_j in N[R_j] (D7; it is a template-fit fake); do NOT cap ε^reco at 1 (C_j < 1 is allowed:
   in-T migration into a cell can exceed the loss).
3. Do NOT change any cut, the existing truth-binned products, the crossx code or outputs (crossx propagation
   is a separate, later step). No ΔR anywhere in the new correction.
4. Truth-matched pairs whose truth muons fail the NTuple truth fiducial (pT < 4.5 or |η| > 2.4) are absent from the
   truth-seeded tree; under D7 they are out-of-T and belong to the template-fit fakes (no longer a reco-eff blind spot,
   but the fake template must include them).

## Design Decisions
- D1 (user 2026-10-09) bin-by-bin, not Bayesian unfolding; supersedes D-U1 of
  `pair_reco_eff_origin_test_and_views.md` (which required a demonstration first — the migration study).
- D2 (user 2026-10-09) third axis m_μμ replaces ΔR (4 uniform bins, D4 of the origin-test doc).
- D3 (user 2026-10-11) ε^reco numerator = N[R_j] incl. out-of-T feed-in (option b). **SUPERSEDED by D7.**
- D7 (user 2026-10-13) ε^reco numerator = truth-matched ∧ reco signal cuts ∧ truth in T (N[R_j] = Σ_i M_ij); out-of-T
  truth-matched pairs are template-fit fakes with the unmatched pairs. Old: feed-in in the numerator (option b), so the
  efficiency absorbed the out-of-T contamination with the MC rate. New: physics reason — a pair whose truth is outside the
  measured fiducial is not signal; subtracting it as background in the template fit (like the unmatched fakes) keeps
  ε^reco a property of signal pairs only, and the inclusive reco-binned ε equals the truth-binned one (Σ_j N[R_j] = N[T∧R]).
  Applies to the MC-only product; the SF product read by the crossx (`pair_reco_eff_reco_binned_sf_pp24_full.root`, option b)
  and the Pb+Pb placeholder built from it are NOT rebuilt (user: no crossx propagation until the pbpb26 ntuple rerun).
- D4 (user 2026-10-11) smoothness test: fine data η^pair binning (0.1), ε from the coarse 3D.
- D6 (user 2026-10-13) R2 resolved by option A: keep the canonical coarse η^pair cells for ε; within-cell shape
  systematic per §3d (rejected: B, a finer η axis for ε only — not requested, binning rule).
- D5 prior-dependence / closure test of C_j (Remaining Work of the feasibility doc) — NOT requested now;
  listed in Remaining Work.

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Decision → `analysis_overview.md` §4b/c, feasibility + origin docs | §1–2 | — | DONE |
| 2 | ε^reco product + evaluator + plots | §2, §3a–b | /review-analysis-code (§2–4) + /review-plot | DONE; reviews PASS (2026-10-13) |
| 3 | Data η^pair smoothness test + MC closure | §3c, §4.1 | /review-analysis-code + /review-plot | DONE Tight (Medium blocked: stale Medium ε_ΔR fit, pT edge 8 vs 9); reviews PASS; R2 resolved by option A (D6) |
| 3c | Within-cell shape systematic | §3d | /review-analysis-code + /review-plot | DONE; reviews PASS it.2 |
| 4 | IntNote reco section | all | /review-note | DONE; APPROVED it.4 (`.claude/logs/review-note-20261013-002420-reco-eff-bin-by-bin.md`) |
| 5 | Docs, INDEX, commit | — | — | DONE 2026-10-13 |
| 6 | D7 numerator (truth in T): product + evaluator marker, plots, MC closure (R∧T), data test, §3d syst rerun | §2, §3a–d | /review-analysis-code + /review-plot | DONE; both APPROVED it.1 (2026-10-13) |
| 7 | D7 in IntNote + analysis_overview + template-fit docs (decision record) | §2 | /review-note | DONE; note APPROVED it.3 (`.claude/logs/review-note-20261013-165550-reco-eff-d7-numerator.md`) |

## Progress Log
- **Step 1 (2026-10-11).** `docs/analysis_overview.md` §4(c) rewritten (bin-by-bin decision, m_μμ axis, ε^reco_j = N[R_j]/N[T_j] with feed-in, output = A_sig numerator, not yet in crossx); §4b note (superseded, not yet applied), §2 A_sig / §5 table / §4f / §6 wording unfolding → bin-by-bin. Decision stamped in `pair_bin_migration_bin_by_bin_feasibility.md` Latest Stage and `pair_reco_eff_origin_test_and_views.md` D-U1; INDEX lines updated.

- **Steps 2–3 (2026-10-11, subagent; scratch merged below).** Product `reco_eff/pair_reco_eff_reco_binned_pp24_full.root`
  (built from the migration file; input re-verified: refill of N[T_j], N[R_j] from the ntuple identical cell by cell),
  evaluator `Utilities/PairRecoEffRecoBinnedEvaluator.h`, macros `plotting_codes/reco_effcy/reco_binned/`; ε plots
  `plots/pp24_reco_effcy_plots/{tight,medium}/signal_region/all_os/reco_binned/{vs_pair_eta,vs_pair_pt,vs_minv}/`;
  smoothness ROOT `~/usatlasdata/dimuon_data/pp_2024/pair_eta_smoothness/`, plots + `smoothness_values_tightWP.txt`
  in `~/usatlasdata/dimuon_data/plots/single_b_analysis/pp24/pair_eta_smoothness/tight/`. Data subclass of
  RDFBasedHistFillingPP reproduces the crossx (516649 pairs; trig+truth-binned 4167.64 pb, raw 1290.29 pb, identical).
  Fine binning = existing `ParamsSet::pair_eta_measured_bins` (= the approved uniform 0.1 on [−2.2, 2.2]).

- **Step 3b (2026-10-11, user: "make the data pair eta distribution plots but now using the same coarse pair eta binning").**
  `plot_pp24_pair_eta_smoothness.cxx` gained `coarse_eta` mode (fine histograms summed onto the canonical coarse cells; exact,
  coarse edges ⊂ fine edges) → `{data_corrected,data_vs_mc_truth,mc_closure}_coarse_eta_tightWP.png` + `coarse_eta_values_tightWP.txt`
  (same dir). Fine-mode table unchanged (diff empty). "Trigger only" marker 25→21: ROOT 6.34 batch PNG draws open squares
  randomly black (reproduced on a synthetic histogram; styles 21/24/26 unaffected). /review-plot APPROVED it.2
  (`.claude/logs/review-plot-20261011-174244-coarse-eta-smoothness.md`): ratio pads widen per panel (cap 3, arrows beyond), data/MC header states the MC
  normalisation; reviewer flagged the coarse data/MC zig-zag as the R2 item (INFO).

- **Step 3c (2026-10-13, user: "go with A").** `plotting_codes/reco_effcy/reco_binned/syst_pp24_pair_eta_within_cell_shape.cxx`
  (§3d) → `within_cell_shape_syst_tightWP.{png,txt,root}` in the smoothness plot dir (`h_within_cell_shape_delta`,
  `h_within_cell_shape_syst` on the canonical (pT, η^pair) axes). Results in R4.

- **Reviews of steps 2–3c (2026-10-13).** /review-analysis-code (all reco_binned macros + evaluator + syst macro) it.1 FAIL on
  ONE warning: syst toys held ε_c fixed (stat inflated up to 20×) → ε_c recomputed per toy, data toys Poisson on raw counts;
  δ unchanged. /review-plot (reco-binned ε plots, fine smoothness plots, δ plot) it.1 FAIL on ONE warning: D_k/ε_k/ε_cell
  undefined on the δ canvas → 3-line header. Combined re-review it.2 PASS. Provenance verified (data macro = crossx selection,
  516649 pairs, 4167.639 pb identical to the crossx file). INFO kept: δ uses m-integrated ε; D_k includes data background.
- **Step 4 (2026-10-13, user: "then update the internal note").** `IntNotes/tex/reconstruction_efficiency.tex` rewritten
  (new §Bin-by-bin correction with ε^bin/purity criteria + Eq. ε^reco = ε^gen ε^bin/(p(1−o)); definition in reco cells
  incl. feed-in; binning (pT, η, m) + fallback + no cap; results vs reco η/pT/m; Application with \PHbox "not yet applied";
  §gap structure within η^pair interval with MC closure, coarse data plot, δ systematic eq. + table; syst list updated).
  Figures: 6 new (recoeff_pp24_rb_vs_{pair_eta,pair_pt,mass}, _eta_smoothness_{mc_closure,data_coarse}, _within_cell_shape_syst);
  7 truth-binned figs removed (git rm); crossx-stages kept (referenced by trigger_efficiency.tex:604, captioned as currently
  applied). placeholder.md item 3 + analysis_overview §4c updated. /review-note log `.claude/logs/review-note-*-reco-eff-bin-by-bin.md`.

- **Step 6 (2026-10-13, user D7).** `reco_binned_cells.h` numerator N[R_S] = Σ_{j∈S} Σ_i M_ij (+ nF info); evaluator `Value()` marker
  (ValueSF byte-identical); builder: `with_sf` THROWS, + `h_NF_<wp>`; closure filler R := … && inT, + `h_mcF_fine_raw_*`, consistency
  vs Σ_i M_ij and R^all (all max rel diff 0); data series `trig_reco_tb` → `trig_reco_applied` (= crossx w_reco, since a2121af the
  SF reco-binned ε, pre-D7); plot labels. Old MC-only product backed up in the session scratchpad only (reproducible from git).
  Product: inclusive **0.58410 ± 0.00057 Tight / 0.63513 Medium** (= generated-binned TR/T exactly); out-of-T feed-in N[F]/(R+F)
  1.325 % / 1.316 %; 282/288 measured; 3 cells ε > 1 (N_eff ≤ 30); R outside axes 3.1e-7. Projections (Tight): η 0.635 0.789 0.449
  0.656 0.476 0.656 0.446 0.787 0.635; pT 0.534 0.609 0.646 0.667 0.673 0.657 0.667 0.655; m 0.589 0.592 0.589 **0.561** (top mass
  interval now 4–6 % low, as the generated-binned: feed-in from m > 2.9 no longer compensates). Data input REFRESHED by another
  session (new GRL, ntuple 2026-10-13 00:41): **509 919** pairs (was 516 649). Coarse all-pT: trig 2421.7, trig+reco-binned 4158.6,
  trig+applied ε^{reco,SF} 4476.3 pb (applied/rb 1.044–1.150 per cell: SFs ⟨SF1SF2⟩ ≈ 0.915 dominate; D7 alone +1.3 %).
  MC closure per coarse cell ≤ 5.8e-15 (per pT panel), all pT 2.6e-6; data vs MC-corrected within-cell χ² central 329.5/9.
  §3d δ moved ≤ 0.004 (R4 table superseded by `within_cell_shape_syst_tightWP.txt`: central all pT +0.0455(32), merged −0.114(35)).
  Reviews: /review-analysis-code APPROVED it.1 (`.claude/logs/review-analysis-code-20261013-164813-d7-numerator-truth-in-t.md`;
  INFOs: data smoothness still contains F_j → +~1.3 % scale; `fill_pp24_fullsim_pair_reco_sf_cells.C` + ValueSF must change with the
  builder at propagation; ε = 1 ± 0 in two N_eff ≤ 2 cells), /review-plot APPROVED it.1 (`.claude/logs/review-plot-20261013-164944-d7-reco-binned-plots.md`).
  Docs: analysis_overview §2 (fakes) + §4c, `low_mass_dimuon_template_fit.md` Design Decisions (fake definition, NOT implemented),
  placeholder.md item 3.

- **Step 7 (2026-10-13).** IntNote `reconstruction_efficiency.tex`: D7 definition (measured region = truth signal region ∩ grid;
  feed-in = template-fit background), Eq. bbb relation without (1−o), new numbers (0.584/0.635, mass dip 3.5–5.5 % vs 3–4 %
  generated-binned), PHbox: D7 not in the crossx (+1.3 % on propagation; no-SF reference 7.1 % lower), application with pp24
  SF weighting, syst table regenerated, 7 figures synced (incl. the crossx stages figure, old GRL). /review-note APPROVED it.3.

## Results & Observations
### R1 — reco-binned ε^reco (Tight; Medium in brackets) — PRE-D7 numbers (feed-in in the numerator); D7 values in the step-6 log
Inclusive 0.5919 ± 0.0006 (0.6436) vs truth-binned 0.5841 (0.6351); 282/288 cells measured (6 empty: pT 105–150,
|η^pair| > 2.0 → 2D fallback); 5 cells ε > 1 (all N_eff ≤ 30), none in 2D; reco R outside axes 0.00057 %. Agrees with
truth-binned within ≈0.01–0.02 below 52 GeV; top mass bin higher (feed-in from m > 2.9).

### R2 — pair-η smoothness: NOT smooth inside coarse cells that contain gap windows (Tight) — USER DECISION
- MC closure: per coarse η cell corrected-R = T to < 4e-15 (exact, by construction; two cells with an empty 3D cell
  ±1e-5). Within cells, the single-muon gap windows carve dips into η^pair (near ±1.15, ≈0, ±2.1; nearly to zero at
  high pT), so the coarse ε under-corrects the dip and over-corrects the rest (up to 1.5–1.8× truth). Cells
  [−1.5,−1], [−0.5,0.5], [1,1.5] (+ edge cells). Within-cell χ² corrected-MC vs MC truth, all pT: 22557/35.
- Truth-binned ε(pT,η,ΔR) shows the same pattern (same coarse η).
- Data follow the MC-reco carved shape: data(corr)/MC(corr R) normalised, all pT: within-cell χ² 586/35, central
  cell 333/9 (dip shallower in data than MC), [±1,±0.5]-type cells 75–84/4, gap cells [±1.5,±1] 11–12/4.
- Consequence: the coarse-cell yield is correct only to the extent MC models the within-cell reco η^pair shape
  (dip depth/width). Scratch estimate (fine-η MC ε instead): central-cell yield +3.8 % (pT-summed), unstable > 37 GeV.

### R3 — coarse-η^pair view (Tight, all pT; yields per cell)
Corrected data smooth across cells; trigger-only dips in the gap cells removed by ε^reco. trig+reco-binned / trig+truth-binned
= 0.97–0.99 per cell (reco-binned lower by 1–3 %, largest in the edge cells). Data(corrected)/Pythia truth shape (panel-normalised):
edge 0.966/0.964, [∓2,∓1.5] 1.022/1.017, [∓1.5,∓1] 1.017/1.031, [∓1,∓0.5] 0.958/0.974, central 1.007 (±0.004–0.010).

### R4 — within-cell shape systematic δ (Tight; §3d; toy stat in brackets, last digits)
| pT [GeV] | edge −2.2…−2 | −2…−1.5 | −1.5…−1 (gap) | −1…−0.5 | central (gap) | 0.5…1 | 1…1.5 (gap) | 1.5…2 | 2…2.2 |
|---|---|---|---|---|---|---|---|---|---|
| 9–12.79 | −0.007(2) | 0.000(0) | +0.002(2) | +0.005(1) | **+0.033(4)** | +0.004(1) | +0.001(2) | 0.000(0) | −0.011(2) |
| 12.79–18.18 | −0.008(2) | 0.000(0) | −0.013(3) | +0.002(1) | **+0.054(6)** | +0.003(1) | −0.012(3) | 0.000(0) | −0.015(2) |
| 18.18–25.85 | −0.008(3) | 0.000(0) | −0.023(8) | 0.000(1) | **+0.044(15)** | +0.002(1) | −0.023(7) | 0.000(0) | −0.010(3) |
| 25.85–36.74 | −0.015(4) | 0.000(1) | −0.013(25) | −0.003(1) | +0.009(34) | 0.000(1) | −0.042(21) | −0.001(0) | −0.004(4) |
| 36.74–150 (merged) | −0.004(12) | +0.001(2) | −0.096(55) | 0.000(2) | **−0.116(34)** | 0.000(2) | −0.058(64) | 0.000(2) | −0.002(6) |
| all pT | −0.008(1) | 0.000(0) | −0.006(2) | +0.002(0) | **+0.045(3)** | +0.003(0) | −0.007(2) | 0.000(0) | −0.012(1) |
(Stat column regenerated 2026-10-13 after /review-analysis-code it.1: ε_c now recomputed per toy, data Poisson; δ values unchanged.)
- [±2,±1.5] cells: δ ≈ 0 (ε flat inside). [±1,±0.5] cells: δ = +0.2–0.5 % although ε_k is NOT flat inside them (0.58–0.69 all pT) —
  data and MC reco share the ε-correlated part of the shape there; so the R3 data/MC zig-zag in those cells is a truth-η-shape
  difference, not an ε effect.
- Central cell: +3–5 % below 26 GeV (data hole shallower than MC). Merged 36.74–150 GeV: −12 ± 3 %, but driven by two fine bins
  (η = ±0.05 at 36.74–52.23 GeV) with ≈3 and 5 data pairs and ε_k 0.03–0.05 → statistics-dominated, NOT an established sign flip.
- [±1.5,±1] gap cells: 0.1–4 % below 37 GeV, 6–10 % above (stat-dominated). Edge cells: 0.2–1.5 %, significant (e.g. all pT −0.8 ± 0.1 %).

## Remaining Work
- Crossx propagation of D7 (user: after the Pb+Pb 2026 ntuple rerun): rebuild `pair_reco_sf_cells` (R_sf with truth in T),
  update `ValueSF`, unblock + rebuild the SF product, rebuild the Pb+Pb placeholder, rerun crossx; propagate the §3d systematic.
- Template fit: add the out-of-T truth-matched pairs (F) to the fake component (decision recorded in `low_mass_dimuon_template_fit.md`);
  truth-matched pairs from truth muons outside the NTuple fiducial need a reco-seeded source.
- Prior-dependence test of C_j.
- Medium WP smoothness + §3d systematic (blocked: Medium ε_ΔR trigger fit predates the 9 GeV pT edge).
- Single-muon-threshold feed-in (§4.4).

## Latest Stage
**2026-10-13 (later): Steps 6–7 DONE (D7; code/plot reviews APPROVED it.1, IntNote APPROVED it.3).** Open: Remaining Work (D7 crossx propagation after the pbpb26 ntuple rerun; template-fit fake component). NOT touched: crossx
code, SF cells/product, Pb+Pb placeholder, crossx outputs (D7 not propagated, user).

## Appendix — merged subagent scratch doc `_sub_reco_binned_1.md` (verbatim)

# _sub_reco_binned_1 — subagent scratch doc (steps 2–3 of pair_reco_eff_bin_by_bin_reco_binned.md)

Canonical doc: `docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md` (Physics Procedure AUTHORITATIVE).
Subagent: no git, writes only `plotting_codes/reco_effcy/reco_binned/`, `Utilities/PairRecoEffRecoBinnedEvaluator.h`,
outputs listed below. Re-read this doc first after any compaction.

## Mandate / Done / Stop-and-ask
- Mandate: run to DONE; do not pause to confirm progress.
- Done = (A) `reco_eff/pair_reco_eff_reco_binned_pp24_full.root` + evaluator header + Tight/Medium plots
  `all_os/reco_binned/{vs_pair_eta,vs_pair_pt,vs_minv}/` + values txt; (B) data smoothness ROOT + plots +
  MC closure plot + smoothness numbers txt, numbers in this doc.
- Stop-and-ask = ANY physics-results-bending ambiguity (report back to orchestrator, do not decide).

## Context gathered (2026-10-11)
- Migration file `~/usatlasdata/pythia_fullsim_full_sample/reco_eff/pair_bin_migration_pp24_full.root`
  (2026-10-06 20:21) newer than ntuple `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root`
  (2026-09-10 01:21) and the RDF hist file (2026-10-06 11:22) it was closed against. Filler throws on closure
  mismatch, so its existence = closure passed at fill time. Will re-verify TNamed + axes + recheck the closure
  against the hist file numerically in the builder.
- pp24 crossx config (RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh): `RDFBasedHistFillingPP pp(24);
  trigger_mode=3; output_generic_hists=true; output_gapcut_hists=true;` everything else default
  (isTight=true, mindR_trig=0.02 → input `muon_pairs_pp_2024_2mu4_mindR_0_02.root`, V1 resonance-cut).
- Crossx selection (RDFBasedHistFillingPP.cxx:612-616): `SignalMinvCutExpr("minv") && SignalPairPtCutExpr("pair_pt")
  && FiducialGapCutExpr(m1 q·η) && FiducialGapCutExpr(m2 q·η) && PairFiducialEtaCutExpr("pair_eta")`, on df_op after
  ApplyMuonWorkingPointFilter (pair_pass_tight). Weight crossx_weight = weight * pp_crossx_lumi_factor;
  corrected = crossx_weight * w_reco * w_trig. w_trig = effcy_pair>0 ? 1/effcy_pair : 0.
- `RDFBasedHistFillingBaseClass::InitOutput()` RECREATEs `output_file` → subclass MUST override SetIOPathsHook
  to redirect output_file to the smoothness file before InitOutput (else crossx hist file clobbered).
- Data plot area: `~/usatlasdata/dimuon_data/plots/single_b_analysis/pp24/` → new subdir `pair_eta_smoothness/<wp>/`.

## Plan
Step A1: builder `plotting_codes/reco_effcy/reco_binned/build_pp24_fullsim_pair_reco_eff_reco_binned.C`
  reads migration file, verifies TNamed/axes/input; builds 3D eps (TH3D pT×eta×m), 2D m-integrated, inclusive,
  N[T], N[R] TH3Ds, A/B/D error propagation; writes TNamed marker. Reports reco-out-of-axes fraction, eps>1 count.
Step A2: evaluator header. Step A3: plotter (Tight+Medium). 
Step B1: data smoothness macro (subclass RDFBasedHistFillingPP). Step B2: MC closure macro. Step B3: plotter + numbers.

## Progress Log
- **A1+A2 DONE (2026-10-11).** Files: `Utilities/PairRecoEffRecoBinnedEvaluator.h` (evaluator + marker namespace
  `PairRecoEffRecoBinnedDefinition` + canonical edge helpers PtEdges/EtaEdges/MinvEdges + CheckAxis),
  `plotting_codes/reco_effcy/reco_binned/reco_binned_cells.h` (RatioOverCells: A/B/D formula for any cell set; truth-binned
  binomial), `.../build_pp24_fullsim_pair_reco_eff_reco_binned.C`.
  Provenance checks passed: migration TNamed input == current ntuple path; definition TNamed matches; migration mtime
  1791332485 > ntuple 1789017686; 1D axes == canonical; h_cell_minv_T and TR_{tight,medium} re-close vs CURRENT hist file
  3D legs with max rel diff 0.
  Output `~/usatlasdata/pythia_fullsim_full_sample/reco_eff/pair_reco_eff_reco_binned_pp24_full.root`:
  h_eps_reco_binned_<wp> (TH3D), _minv_integrated (TH2D), eps_reco_binned_inclusive_<wp>(+_err), h_NR_<wp>, h_NA_<wp>, h_NT,
  frac_R_outside_axes_<wp>, marker TNamed, input TNamed.
  Numbers: Tight: 282/288 cells measured (3 with N[T]=0, 3 with N[R]=0 & N[T]>0 -- all six in pT 105.5-150, |eta|>2.0);
  5 cells eps>1 (all low N_eff: 52-74/eta 2.0-2.2/m 2.445-2.9 1.085+-0.231 NeffT 12.5; 105-150: 1.500+-0.540 (NeffT 6),
  1.113+-0.276 (30), 1.138+-0.789 (10), 1.125+-0.351 (8)); 0 m-integrated 2D cells >1. 3D eps range [0.25,1.50] (extremes all
  NeffT<=30 in top pT bins). Inclusive eps^reco Tight 0.59194+-0.00058 (N[R]=3.989, N[T]=6.7389, A=3.9362), Medium
  0.64360+-0.00057. Reco pairs outside the axes (reco pT>150): 0.000573 % (Tight), 0.000527 % (Medium) of R.
  Evaluator read-back OK (eps(20,0,1.5)=0.5722).
- **A3 DONE.** Plotter `plotting_codes/reco_effcy/reco_binned/plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx` (arg tight_WP).
  Recomputes every plotted point from the migration file with RatioOverCells (summed counts per cell set, never a projected
  ratio) and THROWS unless its single-cell 3D and m-integrated 2D values equal the product (max |diff| 0, both WPs).
  Outputs per WP under `.../plots/pp24_reco_effcy_plots/<wp>/signal_region/all_os/reco_binned/`:
  `vs_pair_eta/pair_eta_in_pair_pt_bins_<wp>WP.png` (m-integ.), `vs_pair_pt/pair_pt_in_pair_eta_bins_<wp>WP.png` (m-integ.),
  `vs_minv/minv_in_pair_pt_bins_<wp>WP.png` (eta-integ.), `reco_binned_values_<wp>WP.txt`. Truth-binned eps (TR/T, weighted
  binomial) overlaid, ±12 % bin-width offset. Observations: reco- and truth-binned agree within ~0.01-0.02 up to 52 GeV; visible
  differences are the top mass bin (reco-binned higher, out-of-T feed-in from m>2.9, o≈0.046) and the top two pT bins (stat).
- **B1/B2 written, running (2026-10-11).**
  Fine binning: `pair_eta_fine_smoothness` = the EXISTING `ParamsSet::pair_eta_measured_bins` (44 uniform bins on
  [-2.2,2.2], opt-in measured-range axis; edges = crossx fine-axis edges inside |eta|<2.2) -- no new edges invented; checked
  in `pair_eta_smoothness_common.h` against ±pair_eta_fiducial_max, width 0.1, and that every coarse eta edge is a fine edge.
  Data: `fill_pp24_data_pair_eta_smoothness.C(tight_WP)` -- subclass `PairEtaSmoothnessPP : RDFBasedHistFillingPP` overriding
  only SetIOPathsHook (parent's, then output_file := own file -- base InitOutput RECREATEs output_file; var1D json absolute),
  FillHistograms (OpenEffcyPtFitFile, OpenPairEfficiencyInputs, ApplyMuonWorkingPointFilter, df_op.Filter(signal_cuts built with
  the identical ParamsSet calls), AddPairEfficiencyWeightColumns, + w_reco_rb), HistPostProcessDataExtra (census print),
  WriteOutputExtra (+TNamed provenance). Config = run_crossx_hist_filling_pp24.sh (pp(24), trigger_mode=3, defaults).
  No shared code edited. Output `~/usatlasdata/dimuon_data/pp_2024/pair_eta_smoothness/pair_eta_smoothness_pp24_2mu4{,_medium_wp}.root`.
  MC: `fill_pp24_fullsim_pair_eta_smoothness_closure.C` -- filler definitions copied verbatim; also re-fills N[T_j], N[R_j]
  cells and compares to the migration file (consistency with current ntuple). Output
  `.../pair_eta_smoothness/pair_eta_smoothness_mc_closure_pp24_fullsim_full.root`.
- **B1 DONE (Tight).** Data run (~4 min) → `pair_eta_smoothness_pp24_2mu4.root` (45 hists). Crossx hist file untouched (mtime
  2026-10-07). **Exact reproduction of the crossx:** unweighted counts 516649 = h2d_counts (incl. overflow); trig+truth-binned-reco
  yield 4167.64 pb = h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts; raw 1290.29 = ..._no_trig_corr. Sum of the 8 pT panels =
  516646 (3 pairs with reco pT > 150, outside every panel; clamped in both evaluators). Reco-binned evaluator census: 516647 3D,
  2 m-integrated fallback, 0 floored, 6 eps>1. Integrated trig+reco-binned = 4108.15 pb vs trig+truth-binned 4167.64 pb (-1.43 %);
  trig only 2425.67.
- **B1 Medium BLOCKED (upstream, not this code):** the Medium crossx path itself throws in OpenPairEfficiencyInputs:
  `DrCorrectionCrossxEvaluator: pair-pT (expo) edge 0 is 8.000000 in the fit file but 9.000000 canonically -- stale fit file`
  (Medium-WP MC dR-correction fit predates the pT>9 adoption; no Medium crossx output exists either). Needs the Medium dR-correction
  refit (run_dr_correction_fits.sh medium) -- out of scope, report to orchestrator. 428-byte stub output deleted. Medium plotter
  for the data test therefore not run; MC side (product, eps plots, MC closure hists) has Medium.
- **B2 DONE.** MC closure (1.6 min). Migration file consistent with CURRENT ntuple: N[T_j], N[R_j]_{tight,medium} re-filled
  cell by cell, max rel diff 0. Evaluator census (Tight): 530134 evals, 37 pT-clamped, 0 fallback, 78 eps>1.
- **B3 DONE (Tight).** `plot_pp24_pair_eta_smoothness.cxx(tight_WP)` → `~/usatlasdata/dimuon_data/plots/single_b_analysis/pp24/
  pair_eta_smoothness/tight/{data_corrected,data_vs_mc_truth,mc_closure}_tightWP.png, smoothness_values_tightWP.txt`.
  KEY FINDINGS:
  1. Closure: per pT panel, sum(corrR)/sum(T)-1 per coarse cell ≤ 3.6e-15 except the two cells containing an unmeasured 3D cell
     (pT 105.5-150, |eta|>2.0; -0.60, -0.50, ~1 weighted entry); all-pT ≤ 1.0e-5 (37 pT>150 clamped pairs).
  2. Within coarse cells the corrected distribution is NOT smooth, in MC closure as much as in data: the single-muon q·eta gap
     windows ([-1.25,-1.05], [-0.10,0.06], [2.2,2.4]) sculpt eta^pair (dips at eta^pair≈±1.15, ≈0 and edges; near-zero at high pT where
     the muons are collinear) and the coarse eps averages over them → compensating over-correction (peaks to 1.5-1.8×T) in the rest of
     the cell. Cells with gap windows: [-1.5,-1], [-0.5,0.5], [1,1.5] (+ edge cells). MC corrected-R vs T within-cell chi2 summed
     over 8 pT panels 31698/279; data vs MC T 30939/275 (reco-binned) and 30567/275 (truth-binned) -- the fine shape is the coarse-
     averaging artefact, identical for the current truth-binned eps(pT,eta,dR) (pulls nearly identical).
  3. Boundary pulls of data/MC-truth (all pT): -2:-7.5, -1.5:+30.2, -1:+3.3, -0.5:+25.4, +0.5:-22.6, +1:-3.1, +1.5:-30.0, +2:+6.6
     vs reference rms 13.9 for non-boundary adjacent bins -- both dominated by the within-cell sculpting, not a boundary effect.
  4. Data follow MC's sculpted shape: data / (MC corrected reco, normalised) boundary pulls (all pT) -2:-0.4 -1.5:+1.5 -1:+0.6
     -0.5:+4.2 +0.5:-2.7 +1:-0.5 +1.5:-2.1 +2:-0.4 (ref rms 3.05); within-cell chi2 sum over pT panels 1236.6/274; largest misfit
     in the central cell [-0.5,0.5] (332.6/9 all pT): data dip at eta^pair≈0 shallower than MC (gap window [-0.10,0.06] acceptance
     data/MC difference) and [-1,-0.5]/[0.5,1] (75/4, 84/4).
  Interpretation: the coarse 3D correction gives the right yield PER COARSE CELL (closure exact); a fine-eta^pair corrected spectrum
  is meaningless inside gap-containing cells. Implication for the crossx: results must be quoted on the coarse eta^pair cells (or
  eps must resolve the gap structure); and the eta≈0 data/MC dip difference means the MC eps within [-0.5,0.5] is averaged with
  a different within-cell weighting in data than in MC (prior/shape dependence of C_j at the few-% level in that cell) -- for the
  orchestrator/user to judge.
- **B3 add-on (indicative, scratch script, not a product).** Size of the within-cell weighting effect: correct data (trigger-
  weighted, fine eta) with a FINE-eta, m-integrated MC eps (R_raw/T per 0.1 bin, same pT panel) instead of the coarse one; relative
  yield change per coarse cell, all pT summed: -0.75 %, -0.01, -0.13, +0.28, **+3.84 %** ([-0.5,0.5]), +0.31, -1.00, -0.01, -1.18 %.
  Central cell per pT panel: +3.3 % (9-12.8), +5.4 % (12.8-18.2), +4.4 % (18.2-25.9), +0.9 % (25.9-36.7); above 36.7 GeV the fine
  eps hits near-zero bins (collinear pairs inside the 0.16-wide gap window) and the estimate is unstable (meaningless). So the data/MC
  difference of the η≈0 dip ⇒ ~4 % shape (prior) dependence of C_j in the central coarse cell at low pT.

## Status (end of subagent run, 2026-10-11)
Done: A (product, evaluator, Tight+Medium plots), B Tight (data, MC closure, plots, numbers). Open: B Medium data test blocked by the
stale Medium dR-correction fit file (upstream). Reviews (/review-analysis-code, /review-plot) NOT run by the subagent -- orchestrator.
Files: Utilities/PairRecoEffRecoBinnedEvaluator.h; plotting_codes/reco_effcy/reco_binned/{reco_binned_cells.h,
pair_eta_smoothness_common.h, build_pp24_fullsim_pair_reco_eff_reco_binned.C, plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx,
fill_pp24_data_pair_eta_smoothness.C, fill_pp24_fullsim_pair_eta_smoothness_closure.C, plot_pp24_pair_eta_smoothness.cxx}
(+ ACLiC .so/.d/.pcm build artefacts in that dir).
