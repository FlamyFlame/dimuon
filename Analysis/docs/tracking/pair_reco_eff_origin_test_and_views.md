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
filters; RDF refill of the pp24 fullsim FULL hist file (pipeline stages 4–9 with SKIP_NTP=1; only
stage 5 produces the new histograms, stages 6–9 reproduce unchanged outputs).

### 4. Negative constraints
1. Do NOT change the ε_reco definition, any cut, or any canonical binning. The APPLIED pair
   category changes only by user decision — given 2026-10-06 (D3: `_single_b` → `_op`), implemented
   here in steps 6–8; no other change to the applied product / crossx.
2. Do NOT compare `_op` vs `_single_b` error bars as if independent.
3. Do NOT retype pair-pT / η / |η| edges; the 4 mass edges are a USER decision.
4. Do NOT build SS reco-eff (D-U2).

## Design Decisions
- D1 — new standalone plotter `plotting_codes/reco_effcy/plot_pp24_fullsim_pair_reco_eff_views.cxx`
  reading the RDF hist file directly (requires the `pair_reco_eff_definition` marker), so any
  category can be drawn without building a product. Existing applied plotter untouched.
- D2 — output tree `plots/pp24_reco_effcy_plots/<wp>/signal_region/{origin_test, <sample>/vs_<var>/}`.

- D3 — (user 2026-10-06, after R1) the APPLIED ε_reco uses ALL truth OS pairs (`_op`) in the
  truth signal region, not only `from_same_b`: no cell-level origin dependence in populated cells,
  +18 % N_eff; the ≈0.15 % inclusive shift is far below other ε_reco uncertainties. Old approach:
  `_single_b`. Physics reason: ε_reco is a detector property of the pair kinematics (R1).
- D4 — (user 2026-10-06) mass-axis alternative uses 4 UNIFORM truth-m_μμ bins over the signal window
  (1.08, 2.9): edges {1.08, 1.535, 1.99, 2.445, 2.9}, DERIVED from ParamsSet::signal_minv_min/max
  (single source `ParamsSet::SignalMinvRecoEffcyBins()`), opt-in suffixed axis `truth_minv_effcy`;
  diagnostic only (the applied 3D stays on ΔR).

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Plotter: origin test (2 figures + values/χ² table) | §3a | /review-plot | DONE |
| 2 | Judge agreement → choose sample for views (Remark 2) | §1, §3a | — | DONE — user chose `_op` (D3); views drawn for both |
| 3 | Plotter: views (η in pT; pT in η; ΔR in pT ×4; pT in ΔR ×4) | §3b | /review-plot | DONE (APPROVED it.2) |
| 4 | Mass 3D: edges from user → RDF + json + refill → mass views | §3c | /review-analysis-code + /review-plot | SUPERSEDED by steps 7–8 (edges = D4) |
| 5 | Docs, INDEX, commit | — | — | DONE (6b6c94c, steps 1–3) |
| 6 | APPLIED ε_reco sample `_single_b` → `_op` (D3): builder kBase + category marker in product; applied plotter text from the marker; evaluator/crossx comments; docs | §2 | /review-analysis-code | DONE (/review-analysis-code APPROVED it.3) |
| 7 | Mass axis (D4): `ParamsSet::SignalMinvRecoEffcyBins()` uniform 4 bins in the signal window; json `truth_minv_effcy`; 3D in `reco_effcy_var3Ds` | §3c | /review-analysis-code | DONE (/review-analysis-code APPROVED it.3) |
| 8 | Refill fullsim RDF (pipeline_pythia_fullsim_pp.sh full, SKIP_NTP=1) → rebuild product → applied plots → pp24 crossx (SKIP_CONDOR=1) → R_AA → mass views (`_op`) | §3c, D3 | /review-plot (mass views) | DONE; code /review-analysis-code APPROVED it.3; mass-view /review-plot APPROVED it.4 (`.claude/logs/review-plot-20261006-113938-reco-eff-mass-views.md`) |

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

- **Steps 6–7 code (2026-10-06).** `Utilities/PairRecoEffDefinition.h` (+ CategoryKey/Value/RdfTag =
  `_op`); builder kBase from CategoryRdfTag, writes the category TNamed; `PairRecoEffEvaluator.h`
  and `plot_pp24_fullsim_pair_reco_eff.cxx` REFUSE a product without it (sample text from the marker);
  `ParamsSet::N_SIGNAL_MINV_BINS_RECO_EFFCY = 4`, `SignalMinvRecoEffcyBins()` (uniform, derived from
  signal_minv_min/max); `RDFBasedHistFillingPythiaFullsim.cxx` 3D {pt_coarse, eta_coarse,
  truth_minv_effcy} + binning `minv_bins_reco_effcy`; `var1D_pythia_fullsim.json` entry. The
  overlay class inherits the new 3D (harmless, not rerun). All 6 touched classes compile (ACLiC).
- **Step 8 (2026-10-06).** `SKIP_NTP=1 pipeline_pythia_fullsim_pp.sh full` stages 4–9 OK
  (`pipelines/logs/fullsim_pp_op_minv_20261006.log`; hist file 11:22, backup
  `backup/...full.bak_20261006_112142.root`). All 449 pre-existing histograms reproduce: contents and
  entries bit-identical, errors to 3.5e-12 relative (ImplicitMT summation order). 18 new
  `truth_minv_effcy` 3Ds; mass axis {1.08, 1.535, 1.99, 2.445, 2.9}, no z under/overflow.
  Product rebuilt (old one backed up as `backup/pair_reco_eff_pp24_full.bak_20261006_single_b.root`):
  inclusive ε_reco Tight 0.5850 → **0.5841**, Medium 0.6360 → **0.6351**; 144/288 3D cells no measure
  (was 145 with single-b), 0/72 2D empty. Applied plots redrawn (Tight + Medium). Mass views (`_op`)
  → `<wp>/signal_region/all_os/{vs_minv,vs_pair_pt}/` (8 figures per WP). pp24 crossx rerun:
  `pipelines/logs/pp_crossx_op_reco_eff_20261006.log`.

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
- Decision status: USER chose `_op` (all OS) for the applied ε_reco (2026-10-06, D3).

### R2 — Mass as the closeness variable (Tight, all OS; revised after the mass-view /review-plot it.1)
At fixed pair pT, ε_reco varies only by ≈4–5 % across the four uniform truth-m_μμ bins (η-integrated,
(max−min)/mean: 4.2/4.1/3.8/5.1/4.2/4.5 % in the six populated low/mid pT bins; e.g. 9–12.8 GeV
0.5322/0.5382/0.5391/0.5166), against a factor-2 drop at the largest populated ΔR of a pT bin
(9–12.8 GeV: 0.541/0.545/0.510/0.268). ΔR stays the better third axis for the applied 3D.
Two named features:
- **Top-mass-bin dip (≈4 %, ≈8–9σ against the adjacent bin, in all six populated pT bins and both
  WPs; the two highest pT bins are statistics-limited).** Explained by RECO-mass migration out of the
  window: the numerator requires reco m < 2.9 (and > 1.08), filled at truth m, so only migration OUT is
  possible. Response `h_minv_zoomin_vs_truth_minv_zoomin_single_b_pass_tight` (single-b; mass
  resolution is a muon-level effect, not origin-dependent; 0.075-GeV bins, reco edges at 1.05/1.125
  and 2.85/2.925, none at 1.08 or 2.9): truth 2.775–2.85 → 21.2 % reco ≥ 2.85, truth 2.85–2.925 →
  75.9 % reco ≥ 2.85; truth 1.05–1.125 → 10.8 % reco < 1.05.
  * UPPER bound at the top edge (uniform truth density; reco ≥ 2.85 counted as lost):
    0.759·0.050/0.455 + 0.212·0.075/0.455 = 0.118 → ≲ 12 %.
  * Model ESTIMATE (Gaussian smearing, locally uniform truth density). σ_m is inverted from the
    response rows with the 0.075-GeV RESPONSE bin width: σ_m ≈ 0.020 GeV at 1.1 GeV (10.8 % row),
    ≈ 0.041–0.049 GeV at 2.9 GeV (exact uniform-bin inversion of the 21.2 % / 24.1 % top-edge rows;
    σ_m/m ≈ 1.8 % → 1.4–1.7 %, i.e. NOT constant). The loss over a 0.455-GeV ANALYSIS mass bin is
    ≈ 0.399 σ_m / 0.455 → ≈ 3.6–4.3 % in the top bin (observed ≈ 4 %) and ≈ 1.8 % in the bottom bin
    (observed ≈ 1 % vs the second bin): top bin matches quantitatively, bottom bin consistent within
    ≈1.5σ. A model estimate, not a measurement (no reco edge at 2.9).
- **Forward slope (2.0 < |η^pair| < 2.2): ε falls ≈17–18 % from the first to the last mass bin**
  (Tight 9–12.8 GeV 0.617/0.593/0.576/0.514; 12.8–18.2 GeV 0.711/0.695/0.664/0.581; Medium same).
  ≈4 % of that is the top-bin migration loss above; first→third bin alone is ≈7 % (0.617 → 0.576).
  HYPOTHESIS (unverified): at fixed pT a larger mass means a larger opening angle, pushing one muon
  towards |η| = 2.4 / the reco |η^pair| < 2.2 edge where acceptance and resolution degrade.
### R3 — Applied product: single-b → all OS (D3), 2026-10-06
Inclusive ε_reco Tight 0.5850 → 0.5841, Medium 0.6360 → 0.6351. Per common 3D cell (143):
denominator-weighted mean |new/old − 1| = 0.24 % (Tight) / 0.21 % (Medium); 34 cells change by > 2 %,
all low-N_eff; max 0.50 in the 2-pair cell 2.0<η^pair<2.2, 105.5–150 GeV, ΔR<0.2 (weight 1.2e-6).
Cells with no measure 145 → 144. Crossx evaluator census: 85.589 % 3D map, 14.411 % ΔR-integrated
fallback, 0 inclusive, 0 floored. The pre-rerun crossx file was not backed up, so the crossx change is
bounded by the product comparison (≲ 0.3 %) rather than measured directly. The R_AA rerun also picked
up the Pb+Pb 2026 input (`raa_*_pbpb23_24_25_26_combined_*`; another session's work), so the R_AA
difference is not attributable to D3 alone; the older `raa_*_pbpb23_24_25_combined_*` files are stale.

## Remaining Work
1. IntNote `IntNotes/tex/reconstruction_efficiency.tex` is stale since D3 (2026-10-06): ALL single-b denominator wording and ALL ε_reco numbers in §sec:eff:reco* — e.g. denominator
   definition l.47–48 ("opposite-sign pair of muons from the same b-hadron decay chain"), text l.30, 41,
   100, captions l.172, 204, 219, 243, marker description l.107 (now also the pair-category marker),
   empty-cell count l.137–139 (145 → 144 in total; recompute the note's sub-counts too), inclusive ε l.154–155, 277, 404 (0.585/0.636 → 0.584/0.635),
   per-bin values l.182, 185. Owned by `intnote_trig_reco_eff_update_2026_09.md` (flagged there).

## Latest Stage
Steps 1–8 DONE. Reviews: steps 1–3 /review-plot APPROVED it.2; steps 6–8 /review-analysis-code APPROVED
it.3 (`.claude/logs/review-analysis-code-20261006-113235-reco-eff-op-and-mass-axis.md`); mass views (step 8) /review-plot APPROVED it.4. Applied ε_reco = all OS
(`_op`), product + crossx + MC-vs-data + R_AA regenerated; mass-axis diagnostic 3D + views built.
Open: Remaining Work 1 (IntNote text, owned elsewhere).
