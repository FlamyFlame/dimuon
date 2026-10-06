# Signal acceptance A_sig = truth mass-window acceptance in (pair pT, η^pair)

Mode: **IMPLEMENTATION**. Opened 2026-10-05.

## Objective

Redefine and rebuild the truth signal acceptance A_sig (Pythia + POWHEG truth) per the user's
measurement definition (2026-10-05), removing the gap cuts it double counts with ε_reco
(`pair_reco_eff_gap_acceptance.md`), and plot its pair-pT AND pair-η dependence.

## Physics Procedure

### 1. Motivation
The final result is an UNFOLDED cross section of truth single-b OS dimuon pairs, measured in
truth pair p_T > 9 GeV and |truth η^pair| < 2.2 (user, 2026-10-05). These two cuts DEFINE the
measured range (like |y| < 1 in CMS arXiv:2212.01636 §3 p.5, or the p_T/y bins of ATLAS
arXiv:1709.03089 §5.2): they are never acceptance-corrected, and pairs outside them are not
probed (truth |η^pair| ≥ 2.2 reconstructed inside = migration, handled by unfolding).
Every OTHER truth-level cut that a pair inside the measured range can fail, and that the
ε_reco denominator does not already contain, must be corrected by A_sig. ε_reco's truth
denominator (2026-09-17) = single-b OS, both muons truth p_T > 4.5 & |η| < 2.4 (NTuple truth
fiducial), truth mass window, truth pair p_T > 9 — so the mass window is the cut A_sig owns.

### 2. Top-level equation
    dσ/dX (truth X, after unfolding) = N_unfolded(X) / (L · A_sig(X)),   X = (pair p_T, η^pair)

    A_sig(cell) = N[ T_meas  AND  m_μμ^truth ∈ (signal_minv_min, signal_minv_max) ]
                / N[ T_meas ]
    T_meas = truth single-b OS pair (`from_same_b`), both muons truth p_T > 4.5 GeV & |η| < 2.4
             (NTuple-stage truth fiducial — kept fiducial, see §4.1), truth pair p_T > 9 GeV,
             |truth η^pair| < 2.2; ALL masses.
Binned in TRUTH (pair p_T, η^pair) on the cross-section's own 2D axes (`pT_bins_150` ×
`pair_eta_crossx`), MC-event weighted. Interpretation (user): in each 2D bin, the fraction of
truth single-b signal kept by the truth mass window. A_sig · ε_reco = N_reco-selected / N[T_meas]:
no double count, no gap.

### 3. Step-by-step method
(a) RDF truth fill (`RDFBasedHistFillingPythiaTruth.cxx`, `RDFBasedHistFillingPowhegTruth.cxx`,
`FillHistogramsSignalAcceptance`): denominator = `from_same_b && pair p_T > 9 && |η^pair| < 2.2`
(truth, cut expressions from ParamsSet); numerator = denominator && truth mass window. NO gap
cut (single-muon q·η windows) anywhere. Write a TNamed provenance marker naming the definition;
the plotter refuses files without it.
(b) Plots: 2D map; A_sig vs pair p_T in coarse η^pair panels; NEW A_sig vs η^pair in coarse
pair-p_T panels (`pair_pt_coarse_bins`, which are fine-axis edges of `pT_bins_150`); cutflow
diagnostic without gap cuts.
(c) A_sig is applied LAST, to the unfolded yield (not yet wired: unfolding is an identity
placeholder).

### 4. Negative constraints
1. Do NOT include any single-muon gap (q·η) window or |η^pair| gap acceptance in A_sig — it is
   in ε_reco (double count).
2. Do NOT acceptance-correct the pair p_T > 9 or |η^pair| < 2.2 cuts (they define the range;
   apply them to numerator AND denominator).
3. Do NOT remove the NTuple truth fiducial muon cuts (p_T > 4.5, |η| < 2.4) — **TO BE
   CLARIFIED** with collaborators (§4.1); no new MC, no new NTuple flag for now.
4. Do NOT change any binning (crossx fine η^pair axis stays 48 bins over ±2.4; its 2.2–2.4 bins
   are structurally empty).

#### 4.1 TO-BE-CLARIFIED (user, 2026-10-05; checking with collaborators)
Whether A_sig should also extrapolate the single-muon cuts p_T 4.5 → 3.7 GeV and |η| 2.4 → 2.5
(the evgen `xAODMultiMuonFilter` limits; 3.7 GeV is needed for MC efficiency, 2.5 is the ID
coverage). Only the existing filtered MC will be used. Until resolved, muon cuts are fiducial.

## Context
- Old A_sig (plots 2026-09-10): num = from_same_b && mass && pair p_T && gap windows ×2 &&
  |η^pair|<2.2; den = from_same_b (any pair p_T / η^pair) → double counted the gap acceptance.
- Literature: ATLAS W/Z 13 TeV arXiv:1603.09222 §2 Eq.(1) p.2 (A at generator level, separate
  from C); ATLAS quarkonium arXiv:1709.03089 §5.2 pp.5–6 (generator-level A, applied per
  candidate at reco p_T, + FSR correction); CMS arXiv:2212.01636 §3 p.5 (range defined by |y|).

## Design Decisions
- D1 — mass window is the only truth cut A_sig owns (user, 2026-10-05: "extrapolate all mass").
- D2 — muon fiducial cuts kept (user, TO-BE-CLARIFIED §4.1).
- D3 — η-dependence view projects the SAME 2D histograms in `pair_pt_coarse_bins` panels
  (CLAUDE.md Binnings rule 2).

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Producers: new num/den + marker (Pythia + POWHEG truth) | §3a, §4.1–2 | /review-analysis-code | DONE (APPROVED it.2) |
| 2 | Plotter: marker guard, η-dependence view; cutflow kCuts | §3b | /review-analysis-code (plots checked C1/C2) | DONE |
| 3 | Compile; rerun Pythia 5.36 truth RDF + POWHEG truth RDF; replot | §3a–b | /review-plot | Pythia DONE; POWHEG blocked (R2) |
| 4 | Docs (analysis_overview, template-fit §3f pointer, signal_selection_change_impact, powheg.md), INDEX, commit | — | — | DONE |

## Progress Log

- **Steps 1–2 code (2026-10-05).** NEW `Utilities/SignalAcceptanceDefinition.h` (Key/Value marker +
  `DenomExpr()` = `from_same_b && SignalPairPtCutExpr(truth_pair_pt) && PairFiducialEtaCutExpr(truth_pair_eta)`,
  `NumExpr()` = Denom && `SignalMinvCutExpr(truth_minv)`). `RDFBasedHistFillingPythiaTruth.cxx` /
  `RDFBasedHistFillingPowhegTruth.cxx` `FillHistogramsSignalAcceptance` use them (gap cuts removed);
  marker written in `WriteOutputExtra` (new override for POWHEG truth, `RDFBasedHistFillingPowheg.h`).
  `SignalAcceptancePlotter.cxx`: marker guard in `Init`, 2D map and η view restricted to |η^pair|<2.2,
  binomial ("B") ratio errors, NEW `DrawAcceptanceEtaByPt` (A_sig vs η^pair in the 8
  `pair_pt_coarse_bins` panels, edge-alignment checked → `*_sig_accept_pair_eta_in_pt_subplots.png`,
  default pT_bins_150 axis only). `plot_sig_accept_cutflow_above_60GeV.cxx`: kCuts = pair pT>9,
  |η^pair|<2.2, m>1.08, m<2.9 (no gap); prints integrated A_sig. Verified both truth NTuples apply
  truth muon pT>4.5 & |η|<2.4 (`PythiaAlgCoreT.c:906-911`, `PowhegAlgCoreT.c:252-257`; NTuples dated
  2026-09-09/10/11, post 4.5-GeV adoption). Backups:
  `.../pythia_5p36TeV/backup/histograms_pythia_5p36TeV_no_data_resonance_cuts.bak_20261005_sigacc.root`,
  `.../powheg_full_sample/backup/histograms_powheg_truth.bak_20261005_sigacc.root`.

- **Review (2026-10-05):** `/review-analysis-code` APPROVED at iteration 2
  (`.claude/logs/review-analysis-code-20261005-225359-signal-acceptance-mass-window.md`). It.1 WARNING:
  `SetRangeUser(±2.2f)` displayed ±2.3 (float rounding) → fixed with `SetMeasuredEtaRange` (bin
  index); cutflow index → `kNMeasuredRangeCuts`. All numbers independently reproduced. Per coarse
  pair-pT bin A_sig = 0.758, 0.750, 0.738, 0.723, 0.710, 0.711, 0.714, 0.721. C3: RUN2-CROSSCHECK
  UNVERIFIED (no comparable Run 2 quantity). INFO left: zero binomial error at A=0/1 edge cells;
  stored `h2d_sig_accept_pt_eta` uses plain Divide (not plotted with errors).

## Results & Observations

### R1 — Pythia 8.3 (5.36 TeV truth, 2026-10-05 refill)
Integrated A_sig (measured range, all pT) = **0.7493** (num 6.11509 / den 8.16069 in
`h2d_sig_accept_{num,denom}_pt_eta`; cutflow 0.749335), pT>60: 0.7064. Single-b truth mass in the
measured range: m<1.08 → 15.9 % (φ peak [0.98,1.06] 2.4 %), m>2.9 → 9.2 % (J/ψ [3.0,3.2] 3.1 %:
B→J/ψX→μμ counts as single-b, user: "extrapolate all mass"). Shape: ≈0.76 at 9 GeV falling smoothly
to ≈0.71 at ~100 GeV; FLAT in η^pair (old gap dips 0.45/0.53 gone — they were the double-counted
gap acceptance). Plots: `plots/single_b_analysis/pythia/pythia_sig_accept_{pair_pt_vs_pair_eta,
pair_pt_in_eta_subplots,pair_eta_in_pt_subplots,above_60GeV,above_60GeV_accum}.png` (+ `pythia_pt_120/`).

### R2 — ★ POWHEG A_sig INVALID: the POWHEG truth NTuple is resonance-vetoed
POWHEG refill integrated A_sig = 0.9654 (num 2691.44 / den 2787.97) vs Pythia 0.749. Root cause:
`PowhegAlgCoreT.c:327-351` applies `ResonanceTagging` UNCONDITIONALLY (no switch), so
`muon_pairs_powheg_bb_truth.root` OS single-b pairs have ZERO weight for m<1.0 and in [2.9,3.3],
[3.5,3.8] (measured: m<1.08 0.8 %, m>2.9 2.7 %, [3.0,3.2] 0.0 %) — exactly the regions the mass
window removes, so the all-mass denominator is missing them. Pythia reads the `_no_data_resonance_cuts`
NTuple. Fix (provenance rule 1): add a no-resonance-cut flag to POWHEG NTuple processing (distinct
suffix) + rerun bb truth on Condor. POWHEG A_sig plots NOT regenerated (old ones are stale: old gap
definition). AWAITING USER.

## Remaining Work
- §4.1 muon-cut extrapolation decision (collaborators).
- Wire A_sig into the crossx after unfolding exists.

## Latest Stage
Pythia side DONE (code approved, refilled, replotted, docs). BLOCKED on user: POWHEG A_sig (R2 —
add a no-resonance-cut POWHEG truth NTuple mode + Condor rerun?) and §4.1 (collaborators).
