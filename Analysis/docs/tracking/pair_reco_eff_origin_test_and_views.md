# Pair reco efficiency: origin-independence test + projection views (+ mass-axis alternative)

Mode: **IMPLEMENTATION**. Opened 2026-10-06.

## Objective

User request 2026-10-06 (items 5, 6 of that message):
1. **Test the hypothesis that ε_reco depends on pair kinematics only, not on the physical origin
   of the pair.** Overlay ε_reco for (1) truth OS pairs in the truth signal region WITHOUT
   `from_same_b` and (2) the same WITH `from_same_b` (= the applied `_single_b`), in
   (a) pair-pT dependence, 9 pair-η panels; (b) ΔR dependence, 8 coarse pair-pT panels.
   If they agree within statistics, (1) replaces (2) (more MC statistics).
2. **Projection views** of ε_reco (sample chosen by 1): η^pair in pT panels; pT in η^pair panels;
   ΔR in pT panels and pT in ΔR panels, each η-integrated and per |η^pair| group; pair mass in pT
   panels (needs a new (pT, η^pair, m_μμ) 3D, 4 mass bins — the mass-as-closeness alternative to ΔR).
3. Organise every plot type in its own sub-directory.

Recorded user decisions of the same message (context for future work, not implemented here):
- D-U1: ε_reco stays defined in TRUTH kinematics; detector migration is handled by
  multidimensional Bayesian unfolding (or sequential 1D only if shown order-independent).
  Bin-by-bin correction only if shown that unfolding is not needed — never assumed.
- D-U2: NO same-sign reco efficiency. The template fit stays trigger-corrected only (fakes /
  hadronic / mixed-event have no truth); reco-eff + unfolding apply to the extracted signal yield
  (`low_mass_dimuon_template_fit.md` §3e, unchanged).

## Physics Procedure

### 1. Motivation
ε_reco is a detector property: given the truth pair kinematics, the probability that both muons are
reconstructed, pass the WP and the reco selection should not depend on which process produced them.
Restricting the denominator to `from_same_b` would then only discard statistics. The hypothesis is
tested, not assumed.

### 2. Definition (unchanged from `pair_reco_eff_gap_acceptance.md` §2 except the origin requirement)
    ε_reco^C(cell) = N[ C ∩ T, both reco-matched, pair passes WP, RECO pair passes the data signal
                        selection incl. gap cuts on RECO q·η (both muons) and RECO |η^pair| < 2.2 ]
                   / N[ C ∩ T ]
    T = truth OS pair, truth muons pT > 4.5 & |η| < 2.4 (NTuple truth fiducial),
        truth m_μμ ∈ (1.08, 2.9) GeV, truth pair pT > 9 GeV; NO truth gap cut.
    C = all truth OS pairs (sample 1, `_op`)  or  `from_same_b` OS pairs (sample 2, `_single_b`).
Binned in TRUTH on the canonical axes: `ParamsSet::pair_pt_coarse_bins` ×
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` × `dr_bins_edges_for_reco_effcy`
(mass alternative: ΔR axis replaced by a 4-bin truth m_μμ axis inside (1.08, 2.9)). MC-event
weighted; binomial errors for weighted histograms.

### 3. Step-by-step method
(a) Origin test — existing histograms, no refill: the RDF already fills the coarse 3D
`h_truth_dr_effcy_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse{_op,_single_b}{_pass_signal_truth,
_pass_<wp>_and_signal_truth_and_reco}`. Projections are ratios of projected num/den (never
projections of a ratio). Statistical judgement: `_single_b` ⊂ `_op` (df_single_b =
df_op.Filter(from_same_b)), so the two overlaid curves are ~85 % correlated and their overlap is NOT
a test. The quantitative test compares `_single_b` with its INDEPENDENT complement (`_op` −
`_single_b`; Σw² of the complement = Σw²_op − Σw²_single_b), per 3D cell and per plotted bin:
pulls and χ²/ndf.
(b) Views — projections of the SAME coarse 3D (CLAUDE.md Binnings rule 2); |η^pair| groups from the
canonical fold `MakeDrEtaGroups(axis, true)` (edges {0, 1, 2, 2.2}, derived, never retyped).
(c) Mass alternative — new 3D (truth pT coarse, truth η^pair coarse, truth m_μμ 4 bins) on the same
filters; RDF refill of the pp24 fullsim FULL hist file (stage 5 only).

### 4. Negative constraints
1. Do NOT change the ε_reco definition, any cut, or any canonical binning; do NOT change the applied
   product / crossx in this doc (switching the applied sample is a separate, user-approved step).
2. Do NOT compare `_op` vs `_single_b` error bars as if independent.
3. Do NOT retype pair-pT / η / |η| edges; the 4 mass edges are a USER decision.
4. Do NOT build SS reco-eff (D-U2).

## Design Decisions
- D1 — new standalone plotter `plotting_codes/reco_effcy/plot_pp24_fullsim_pair_reco_eff_views.cxx`
  reading the RDF hist file directly (requires the `pair_reco_eff_definition` marker), so any
  category can be drawn without building a product. Existing applied plotter untouched.
- D2 — output tree `plots/pp24_reco_effcy_plots/<wp>/signal_region/{origin_test, <sample>/vs_<var>/}`.

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Plotter: origin test (2 figures + values/χ² table) | §3a | /review-plot | DONE |
| 2 | Judge agreement → choose sample for views (Remark 2) | §1, §3a | — | USER DECISION (R1); views drawn for both |
| 3 | Plotter: views (η in pT; pT in η; ΔR in pT ×4; pT in ΔR ×4) | §3b | /review-plot | DONE (APPROVED it.2) |
| 4 | Mass 3D: edges from user → RDF + json + refill → mass views | §3c | /review-analysis-code + /review-plot | BLOCKED on mass edges |
| 5 | Docs, INDEX, commit | — | — | |

## Progress Log

- **Steps 1–2 DONE (2026-10-06).** NEW `plotting_codes/reco_effcy/plot_pp24_fullsim_pair_reco_eff_views.cxx`
  (ACLiC OK). Origin test run for Tight + Medium → `~/usatlasdata/pythia_fullsim_full_sample/plots/
  pp24_reco_effcy_plots/{tight,medium}/signal_region/origin_test/{pair_pt_in_pair_eta_bins,
  dr_in_pair_pt_bins}_<wp>WP.png` + `origin_test_values_<wp>WP.txt`. Numbers in R1.

- **Steps 3 + review DONE (2026-10-06).** Views (10 figures per sample × WP) for `_single_b` →
  `<wp>/signal_region/single_b/` and `_op` → `<wp>/signal_region/all_os/`, each with
  `vs_pair_eta/`, `vs_pair_pt/`, `vs_dr/` + `views_dr_values_<wp>WP.txt`. Origin test got a
  single-b/all ratio pad (error propagated through the independent single-b and complement pieces)
  and N_eff-thresholded 3D χ². `/review-plot` APPROVED at it.2
  (`.claude/logs/review-plot-20261006-021319-pair-reco-eff-views.md`); it.1 WARNING = overstated
  origin significance → R1 restated. C3: RUN2-CROSSCHECK UNVERIFIED (Run 2 used single-muon
  MCP/HION-2019-58 values, no pair/gap analogue).

## Results & Observations

### R1 — Origin test: ≈1 % coherent difference at projection level, NONE cell-by-cell in populated cells
(revised after /review-plot it.1, which reproduced every number and flagged the first reading as overstated)
Input `histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` (2026-09-17 20:39, marker OK).
Single-b vs its INDEPENDENT complement (`_op` − `_single_b`, Σw² subtracted):

| | Tight | Medium |
|---|---|---|
| fig A (pT in 9 η^pair panels, ΔR-int.) χ²/ndf, p | 104.7/71, 0.006 | 105.2/71, 0.005 |
| fig B (ΔR in 8 pT panels, η-int.) χ²/ndf, p | 39.1/16, 0.001 | 29.5/16, 0.021 |
| all 3D cells | 182.7/142, 0.012 | 188.4/142, 0.006 |
| 3D cells, N_eff ≥ 100 in both | **84.8/88, 0.58** | **78.4/88, 0.76** |
| 3D cells, N_eff ≥ 1000 in both | **37.7/35, 0.35** | **32.7/35, 0.58** |

- The all-cell significance is driven by low-N_eff cells where the Gaussian weighted-binomial error
  is unreliable (max pull 5.11 at η^pair∈[−2,−1.5], pT 52–74 GeV: complement N_eff = 38, ε = 0.957).
- Projection level: inclusive ε single-b 0.5850 vs complement 0.5792 (N_eff 628k / 115k; +0.0058,
  ≈1 % relative, ≈3.7σ). Fig B leading pull (pT 18–26 GeV, ΔR<0.2, −3.9σ) is partly COMPOSITION:
  reweighting the complement to the single-b η^pair mix reduces it to ≈2.9σ (reviewer). Projections
  integrate over axes whose distributions differ by origin, so they are not an origin test per se.
- Practical size: single-b = 84.9 % of the `_op` denominator weight → `_op` gains only ≈18 % N_eff,
  and `_op` vs `_single_b` differ by ≈0.001 inclusive (0.15 %); ratio panels within ±1 % except
  low-N_eff bins.
- Decision status: Remark 2 ("lift unless significant difference") is NOT clearly decided by this
  evidence → **USER DECISION**. Views drawn for BOTH (`single_b/`, `all_os/`).
## Remaining Work
1. USER: choose `_op` or `_single_b` as the ε_reco sample (R1). Switching the APPLIED product is a
   separate step (builder + crossx rerun).
2. USER: 4 truth-mass bin edges inside (1.08, 2.9) for the mass-axis alternative (Step 4).
   Candidates: uniform {1.08, 1.535, 1.99, 2.445, 2.9} (single-b shares 23/30/28/19 %), or equal
   single-b population ≈{1.08, 1.57, 1.95, 2.34, 2.9} (quantiles from the 0.075-GeV truth histogram).

## Latest Stage
Steps 1–3 DONE and reviewed. Waiting on the user for Remaining Work 1–2 (sample choice, mass edges).
Input `~/usatlasdata/pythia_fullsim_full_sample/histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` (2026-09-17 20:39, marker present).
