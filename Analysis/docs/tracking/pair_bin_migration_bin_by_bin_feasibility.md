# Pair bin migration: reco efficiency / bin-by-bin efficiency / purity — is bin-by-bin unfolding sufficient?

Mode: **IMPLEMENTATION** (diagnostic study). Opened 2026-10-06.

## Objective

User request 2026-10-06: decide whether a bin-by-bin correction suffices for this analysis or a
full (Bayesian) unfolding is needed (recorded decision D-U1 in
`pair_reco_eff_origin_test_and_views.md`: bin-by-bin only if shown sufficient, never assumed).
Compute bin-by-bin **efficiency, stability and purity** for **pair pT, pair η and pair mass**
(and ΔR, for comparison), on the canonical binnings, and organise the plots in sub-directories.
If the pair mass gives better migration behaviour than ΔR as the third axis, a 3D bin-by-bin
correction in (pair pT, pair η, m_μμ) replacing ΔR is the likely follow-up (user decision, not
taken here).

## Physics Procedure

### 1. Motivation
The measured signal yield is binned in RECO pair kinematics; the cross section is quoted in TRUTH
kinematics. Migration between bins (finite resolution) and across the boundary of the measured
phase space (selection edges) must be corrected. A bin-by-bin correction C_j = N^true_j / N^reco_j
corrects migration only to the extent that the MC truth spectrum shape is right; it is adequate when
migrations are small (high purity AND bin-by-bin efficiency in every bin). This study measures them.

### 2. Definitions
Populations (pp24 Pythia fullsim FULL, all truth OS pairs `_op` = the applied ε_reco category D3,
MC-event weighted; WP = Tight nominal, Medium for the WP systematic):
- **T (truth, measured phase space):** `pass_signal_truth` (truth m_μμ ∈ (1.08, 2.9) GeV, truth
  pair pT > 9 GeV) AND truth (pT, η^pair, m_μμ) inside the canonical axis ranges
  (pT < 150 GeV top edge, |η^pair| < 2.2). NO truth gap cut (gap = acceptance inside ε_reco).
  Truth muons carry the NTuple truth fiducial pT > 4.5, |η| < 2.4 (the tree is truth-seeded).
- **R (reco selected):** `pair_pass_<wp> && pass_signal_reco` (both muons reco-matched, WP incl.
  reco pT > 4.5 / |η| < 2.4 / IP / Δp/p, reco mass window, reco pair pT > 9, gap cuts on reco q·η
  of both muons and reco |η^pair| < 2.2) — exactly the applied ε_reco reco leg, no extra reco
  range cap. A reco value outside an axis (only reco pT > 150 GeV occurs) lands in the overflow and
  counts as migrated OUT of the measured range. Truth may be ANYWHERE (inside or outside T).
For a binning of variable x (1D: one axis; 3D: the cell index of (pT, η^pair, Z), Z = m_μμ or ΔR),
truth bin i, reco bin j:

With M_ij = N[T_i ∧ R_j] (truth bin i in T, reco bin j; reco overflow kept in the row):

    reco efficiency        ε_i     = N[T_i ∧ R] / N[T_i]        (= the applied ε_reco when x = 3D cell)
    bin-by-bin efficiency  ε^bin_i = M_ii / Σ_j M_ij = N[T_i ∧ R_i] / N[T_i ∧ R]
                                     (reconstructed truth-bin-i pairs that stay; ATLAS name — called
                                      "stability" in H1/CMS-TOP-14-013, user 2026-10-08)
    purity                 p_j     = M_jj / Σ_{i∈T} M_ij = N[T_j ∧ R_j] / N[T ∧ R_j]
                                     (MIGRATION-IN only: denominator = truth-matched reco-bin-j pairs
                                      whose truth is in T; user 2026-10-08, D3)
    outside (table only)   o_j     = N[R_j ∧ ¬T] / N[R_j]       (feed-in from outside T; not plotted)
    bin-by-bin correction  C_j     = N[T_j] / N[R_j] = p_j (1 − o_j) / (ε_j ε^bin_j)

Purity/ε^bin separate migration from the reconstruction efficiency (ε) and from the fake/hadronic
background (removed by the Δp/p template fit; absent from the truth-seeded tree anyway). The
out-of-T feed-in o_j is NOT in p_j; it stays in C_j (it is real-muon yield in data) and in the table.

Errors: weighted binomial for subset ratios A ⊂ B: var = [(1−2r) Σw²_A + r² Σw²_B] / B².
1D: "same bin" = same bin of x (the other axes only enter via T/R membership).
3D: "same bin" = same (pT, η^pair, Z) cell — the quantity that matters for a 3D bin-by-bin.

### 3. Step-by-step method
(a) Fill: RDF over the ntuple-processing output `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root`,
tree `muon_pair_tree_sign2`, weight `weight` — exactly the RDF class's `df_op_weighted`. Filters
use the SAME expression builders as `RDFBasedHistFillingPythiaFullsim.cxx`
(`ParamsSet::SignalMinvCutExpr`, `SignalPairPtCutExpr`, `RDFBasedHistFillingPythiaFullsim::RecoGapCutExpr` (inherited static)).
Axes from `ParamsSet::pair_pt_coarse_bins`, `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`,
`ParamsSet::SignalMinvRecoEffcyBins()`, `CommonEffcyConfig::dr_bins_edges_for_reco_effcy`.
Store per WP: N[T](truth x), N[T∧R](truth x), M(truth x, reco x) for T∧R, N[R](reco x),
N[R∧¬T](reco x) — for the 4 1D variables and the two 3D cell indices (full cell×cell matrix kept
for a later Bayesian comparison). Closure check: ε from this file must reproduce the RDF-class
3D ε_reco cell by cell (T∧R ≡ `_op_pass_<wp>_and_signal_truth_and_reco` inside the axis ranges).
(b) Plot: per variable ε / ε^bin / p vs x (+ C_j) and the row-normalised response; per 3D option,
cell purity / bin-by-bin efficiency / ε in pT panels; a summary table (min / weighted mean over cells
with N_eff ≥ 100). o_j is in the table only (D3).

### 4. Negative constraints
1. Do NOT change any cut, the ε_reco definition, or any binning; read every edge from the source
   of truth (CLAUDE.md §Binnings).
2. (Superseded 2026-10-08, D3.) Purity now EXCLUDES the out-of-T feed-in by definition; it must NOT
   then be dropped from C_j, which keeps N[R_j] incl. R_j ∧ ¬T (real-muon yield in data).
3. Known blind spot, NOT to be hidden: reco pairs whose TRUTH muon fails the NTuple truth fiducial
   (truth pT < 4.5 GeV or truth |η| > 2.4) are absent from the truth-seeded tree, so feed-in across
   the single-muon threshold is not counted in o_j / C_j (it is outside T, so not in p_j by definition).
4. This study does NOT decide the unfolding method or replace ΔR — that is the user's decision.

## Design Decisions
- D1 — standalone RDF filler (not a new RDF-class booking + refill): diagnostic study, keeps the
  shared hist file untouched; provenance = ntuple-processing output + the class's own expression
  builders (no retyped cut). Output `reco_eff/pair_bin_migration_pp24_full.root`.
- D1b — R carries NO reco axis-range cap (was: "reco values inside the axes" in the first draft of §2):
  keeps ε_i identical to the applied ε_reco (closure exact); reco overflow (cell-index histograms: underflow) = migration out
  (0.0006 % of Tight R). Flagged by /review-analysis-code it.1.
- D3 (2026-10-08, user) — terminology + purity definition. Old: "stability" s_i; purity
  p_j = N[T_j∧R_j]/N[R_j] (denominator incl. out-of-T feed-in); o_j plotted. New: "bin-by-bin
  efficiency" ε^bin_i (same formula, ATLAS naming); purity p_j = M_jj/Σ_{i∈T} M_ij so it measures
  migration-in only, separated from reco efficiency and from fakes/hadronic background (template
  fit); o_j removed from the plots (kept in the table, and inside C_j). Reference for the stability
  formula: CMS-TOP-14-013 (arXiv:1703.01630) §6; H1 hep-ex/0012053 §4.4 instead folds ε into it.
- D2 — plots under `plots/pp24_reco_effcy_plots/<wp>/signal_region/all_os/bin_migration/`
  `{pair_pt, pair_eta, minv, dr}/` (1D) and `{cells_pt_eta_minv, cells_pt_eta_dr}/` (3D), + summary txt.

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Filler macro + closure vs RDF-class ε_reco | §2, §3a | /review-analysis-code (§2–4) | DONE — APPROVED it.2 (`.claude/logs/review-analysis-code-20261006-201009-pair-bin-migration-filler.md`) |
| 2 | Plotter: 1D ε/s/p/o + response; 3D cell views; summary table | §3b | /review-plot | DONE — APPROVED it.3 |
| 3 | Judge criteria (purity/stability ≳ 0.5–0.8; minv vs ΔR) → report to user | §1 | — | DONE (R1–R3; decision with user) |
| 4 | Docs, INDEX, commit | — | — | DONE |
| 5 | D3: rename stability → bin-by-bin efficiency, purity = M_jj/Σ_{i∈T}M_ij, drop o_j lines; replot | §2 (D3) | /review-analysis-code + /review-plot | DONE — code APPROVED it.4, plot APPROVED it.2 |

## Progress Log

- **Steps 1–2 (2026-10-06).** Filler `plotting_codes/reco_effcy/bin_migration/fill_pp24_fullsim_pair_bin_migration.C`
  → `~/usatlasdata/pythia_fullsim_full_sample/reco_eff/pair_bin_migration_pp24_full.root` (2.6 min, MT).
  R = exactly the applied ε_reco reco leg (no extra reco-range cap; reco pT > 150 → overflow =
  migrated out; 0.0006 % of Tight R). T pairs outside the ΔR axis: 0. **Closure: N[T] and N[T∧R]
  in both 3D flattenings reproduce the RDF-class `truth_{dr,minv}_effcy` 3D legs (both WPs)
  bit-identically (max rel diff 0).** Plotter `plot_pp24_fullsim_pair_bin_migration.cxx` (Tight +
  Medium) → `plots/pp24_reco_effcy_plots/<wp>/signal_region/all_os/bin_migration/` with
  `{pair_pt,pair_eta,minv,dr}/{eff_stab_pur,response}_<wp>WP.png`,
  `cells_pt_eta_{minv,dr}/{purity,stability,efficiency}_<wp>WP.png`,
  `cells_compare/cell_purity_stability_<wp>WP.png`, `bin_migration_values_<wp>WP.txt`.
- **Reviews (2026-10-06).** /review-analysis-code APPROVED it.2 (it.1 WARNING: R-definition text in
  §2 / TNamed contradicted the code → text fixed, D1b). /review-plot APPROVED it.3
  (`.claude/logs/review-plot-20261006-202140-pair-bin-migration-plots.md`): "opposite-sign" header, physical panel
  ranges, o defined in legend, cell points with 68 % Wilson intervals and drawn only if denominator
  N_eff ≥ 10, compare plot = cumulative N[T] fraction below x on documented log y.

- **Step 5 (2026-10-08, D3).** Plotter only (filler unchanged): ε^bin = M_ii/Σ_j M_ij (renamed from
  stability, values unchanged), purity = M_jj/Σ_{i∈T} M_ij (identity Σ_{i∈T} M_ij = N[R_j] − N[Rout_j]
  exact, all variables × WPs), o_j table-only, C_j unchanged. 1D series offset ±12 % bin width; outputs
  renamed `{x}/eff_purity_<wp>WP.png`, `cells_pt_eta_{minv,dr}/{purity,bin_by_bin_efficiency,reco_efficiency}_<wp>WP.png`,
  `cells_compare/cell_purity_bin_by_bin_eff_<wp>WP.png` (old-name PNGs deleted). New purity: pT
  0.975→0.895, m 0.972/0.946/0.936/0.959, ΔR top bin 0.568→0.958. Under D3 ΔR ≥ m_μμ on every
  purity/ε^bin metric; the ΔR edge feed-in now shows only in o_j/C_j (max o for N_eff ≥ 100: ΔR 0.18
  vs m 0.090). Reviews: `.claude/logs/review-analysis-code-20261008-151400-bin-migration-d3-purity.md`
  (APPROVED it.4: doc §3b/R1–R3 fixes, feed-in counts corrected), `.claude/logs/review-plot-20261008-151427-bin-migration-d3-plots.md`
  (APPROVED it.2: compare-plot y title clipped → fixed; marker overlap → offset).

## Results & Observations
(Definitions per §2 after D3, 2026-10-08: ε^bin = bin-by-bin efficiency (formerly "stability",
values unchanged); purity p = M_jj/Σ_{i∈T} M_ij (migration-in only). Old purity = p (1 − o_j).)

### R3 — Verdict on the bin-by-bin criteria (input for the user's decision)
- Migration criterion: MET. Every 1D bin has p, ε^bin ≥ 0.89; every 3D cell with N_eff ≥ 100 in both ratio
  denominators has p ≥ 0.81 and ε^bin ≥ 0.80 (both axis options, both WPs); response matrices near-diagonal (only
  nearest-neighbour off-diagonals, ≤ 8 %).
- NOT yet tested (needed before adopting bin-by-bin): prior (truth-shape) dependence of C_j and a
  closure / Bayesian comparison. Most exposed: pT top bins (largest migration), and every bin with
  large out-of-T feed-in o_j — this feed-in no longer lowers p but sits entirely in C_j and depends on
  the truth spectrum OUTSIDE T (m bin 4: J/ψ-side tail; ΔR kinematic-edge cells).
- Third axis (re-stated under D3): by purity / ε^bin alone ΔR is now at least as good as m_μμ in
  every metric (means 0.961 vs 0.92; min p 0.846 vs 0.816; N[T] fraction with p < 0.8: 0.0000 vs
  0.0002). The ΔR difference that remains is the edge feed-in (o_j 0.25–0.73 in the
  ΔR 0.6–1.0 cells at pT 9–12.8 GeV, o = 0.086 in pT 12.8–18.2 / |η|<0.5 / ΔR 0.4–0.6), which is a
  C_j / prior-dependence issue, not a migration one. In cells with N_eff[T] ≥ 100 the largest feed-in
  is o = 0.18 for ΔR (pT 12.8–18.2, η −2.2…−2.0, ΔR 0.4–0.6) vs o = 0.090 for m_μμ (pT 36.7–52.2,
  η −2.0…−1.5, m 2.445–2.9), i.e. ≈ 2× larger for ΔR; 42 m-cells (Tight) have o > 0.05, 36 of them
  in the top mass bin, mostly forward (J/ψ-side tail). m_μμ: 285/288 cells populated (ΔR: 144/288). The choice therefore rests on the prior-dependence test.

### R1 — 1D (Tight; Medium p and ε^bin within ≈ 0.003 except ΔR top bin ε^bin 0.955→0.938)
| var | ε^bin range | purity range | outside o_j (table only) | note |
|---|---|---|---|---|
| pair pT (8 bins) | 0.978 → 0.917 (falls with pT) | 0.975 → 0.895 | 0.012–0.029 | resolution grows with pT; top bin lowest |
| pair η (9 bins) | ≥ 0.994 | ≥ 0.994 | 0.009–0.022 | η migration negligible |
| m_μμ (4 bins) | 0.975 / 0.947 / 0.936 / 0.954 | 0.972 / 0.946 / 0.936 / 0.959 | 0.015 / 0.001 / 0.001 / 0.046 | window-edge feed-in in edge bins (in C only) |
| ΔR (4 bins) | 0.994 / 0.996 / 0.993 / 0.955 | 0.994 / 0.996 / 0.994 / 0.958 | 0.014 / 0.006 / 0.038 / **0.407** | top bin = kinematic edge, fed from outside T (C = 2.22); ≈ 0.05 % of T |
Mass-window and pT>9 edge migration OUT of the window fails R entirely, so it shows up as an
ε loss (the R2 top-bin dip of `pair_reco_eff_origin_test_and_views.md`), not as an ε^bin loss;
its mirror image (feed-in) is o_j, carried by C_j. Both depend on the truth spectrum beyond the
boundary. Mass bin 4's o = 0.046 is feed-in from truth m > 2.9 GeV (plausibly the J/ψ low-side tail)
→ C there depends on the MC's J/ψ modelling.

### R2 — 3D cells: (pT, η^pair, m_μμ) vs (pT, η^pair, ΔR) (Tight; Medium in brackets)
| | m_μμ cells | ΔR cells |
|---|---|---|
| populated cells (N[T]>0) | 285 / 288 | 144 / 288 |
| N[T]-weighted mean ε^bin / purity | 0.922 / 0.922 (0.920 / 0.920) | 0.961 / 0.961 (0.960 / 0.960) |
| min ε^bin / purity, cells N_eff ≥ 100 (both denominators) | 0.806 / 0.816 (182 cells) | 0.823 / 0.846 (107 cells) |
| N[T] fraction in cells with purity < 0.8 | 0.02 % | 0.00 % |
| largest out-of-T feed-in o_j, cells N_eff[T] ≥ 100 | 0.090 (top mass bin, η −2.0…−1.5) | 0.18 (N_eff 317); 0.086 at N_eff 2584 |
| low-N_eff high-o cells | max 0.33 (one cell: pT 105.5–150, η 1.5–2.0, m 1.08–1.535; N_eff 8) | 0.25–0.73 (ΔR 0.6–1.0, pT 9–12.8; N_eff ≤ 82); 0.33–0.50 (pT 36.7–52.2, ΔR 0.2–0.4; N_eff ≤ 21) |
- ΔR: ΔR is measured precisely → higher p / ε^bin. Cells cut by the kinematic band (ΔR ≈ 2m/pT at
  the window/pT edges) are fed from outside T; under D3 this lowers neither p nor ε^bin but raises C,
  e.g. pT 12.8–18.2, |η|<0.5, ΔR 0.4–0.6: p = 0.896, o = 0.086, C = 2.05.
- m_μμ: uniformly ≈ 0.92 (mass resolution σ_m ≈ 0.02–0.05 GeV vs 0.455 GeV bins → ≈ 5–8 % leaves per
  bin); lowest cells are high-pT (resolution) and forward, e.g. pT 52–74, η 1.0–1.5, m 1.99–2.45:
  p = 0.78 (N_eff[T] 124, ratio-denominator N_eff < 100 → outside the min statistic).
- Low-N_eff cells (pT > 52 GeV, forward) are statistics-limited, not migration-limited.

## Remaining Work
- Prior-dependence test of C_j (reweight MC truth shape) and Bayesian-vs-bin-by-bin comparison —
  only if the user wants to pursue bin-by-bin after this study.
- Quantify the single-muon-threshold feed-in blind spot (§4.3), e.g. from the reco-gated
  `_mc_trig_single_muon` tree.

## Latest Stage
Steps 1–5 DONE. Awaiting user decision: pursue bin-by-bin (→ prior-dependence + closure tests,
Remaining Work) and m_μμ vs ΔR as 3rd axis (now decided by the C_j / feed-in prior dependence, R3).
