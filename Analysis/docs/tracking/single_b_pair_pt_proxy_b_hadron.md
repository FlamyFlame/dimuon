# Single-b: pair pT as a proxy for the b-hadron pT (truth-level study)

## Objective
Quantify how good `pair pT` is as a proxy of the parent b-hadron pT for truth single-b OS pairs
(vs the single-muon pT), with (1) the ratio R = pT / pT(last b hadron) for pair, leading and
subleading muon; (2) the 2D (b-hadron pT, pair pT); (3) two sanity checks of R_pair: sample
(pp24 fullsim truth, pp + proton PDF  vs  central truth-only, 4 isospin + nPDF) and selection
(from_same_b only vs + truth signal cuts). Add (1)+(2) to the internal note (new truth-level
study section).

## Context
- Pre-existing 1D plot: `plotting_codes/pythia_plotting_codes/plot_single_b_pt_ratio_to_last_b_hadron.cxx`
  (Apr 2026): CENTRAL 5.36 TeV truth sample (`pythia_truth_full_sample/pythia_5p36TeV/muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root`),
  NO signal cuts (from_same_b only), predates the muon pT 4.5 adoption, y label "[ub]" is wrong (weights are nb).
  NO 2D plot code or output exists anywhere (grep of Analysis/, IntNotes/, plots/).
- Both truth pair trees (central truth, pp24 fullsim FULL) carry identical columns incl.
  `m{1,2}_last_b_hadron_prt_pt_eta_phi_m`, `truth_pair_pt`, `truth_minv`, `weight`; muon truth
  fiducial pT>4.5, |eta|<2.4 is applied at NTuple stage in BOTH (min pT 4.50/4.50).
- m1 is the leading muon (m1pt<m2pt in 0 pairs); last-b-hadron of m1 vs m2 differ in 191/1.70M (truth), 88/1.03M (fullsim).
- Truth signal cuts = `SignalAcceptanceDefinition::NumExpr()` (from_same_b, pair pT>9, |eta^pair|<2.2, mass window); no gap cuts.

## Latest Stage
Plan: one macro `plot_pair_pt_proxy_b_hadron.cxx` -> 4 PNG (R 3-curve; 2D; sanity-sample; sanity-cuts);
then /review-plot; new note section `IntNotes/tex/truth_level_studies.tex` + manifest entries; /review-note.

## Progress Log

### 2026-10-08 - DONE
- Macro `plotting_codes/pythia_plotting_codes/plot_pair_pt_proxy_b_hadron.cxx` (replaces the deleted Apr-2026
  `plot_single_b_pt_ratio_to_last_b_hadron.cxx`) -> `plots/truth_level_studies/pair_pt_proxy_b_hadron/` (4 PNG; /review-plot APPROVED it. 2).
  Both samples = NTuple-processing output (muon_pair_tree_sign2); selection `SignalAcceptanceDefinition::NumExpr()`; weights nb.
- Q1 answer: old 1D plot used the CENTRAL truth sample but had NO signal cuts, pre-4.5 GeV cut, label [ub] wrong (nb); no 2D code/output existed. Both remade.
- Results (central, signal sel.): <R_pair>=0.556, <R_lead>=0.335, <R_sub>=0.224; <R_pair> 0.539 (9-12.8 GeV) -> 0.678 (105.5-150).
- Sample check: fullsim pp/proton-PDF <R>=0.562 (+1.1%); tilt 0.9->1.1 over R; per-pT-hat-slice <R> agree <=0.33%; reweighting to central slice fractions gives 0.5564 => difference = slice-weight mix (isospin/PDF/normalisation NOT separated).
  Tables: `R_mean_tables.txt` (scripts `check_R_*.C`).
- Cuts check: <R> 0.557 -> 0.556; peak +9%, low-R tail removed, R>=0.9 depleted 20-25%.
- IntNote: new `IntNotes/tex/truth_level_studies.tex` (+ master include, 4 manifest figures); /review-note 2 it., last wording fix un-re-reviewed; builds clean (46 pp).
