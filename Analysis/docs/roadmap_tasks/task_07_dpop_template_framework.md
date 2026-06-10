# Task 07 — Momentum-imbalance (Δp/p) significance + template-fit framework

**Roadmap:** `analysis_roadmap_2026_06.md` §Q4 / chain [8]
**Can start:** immediately (framework with dummy templates).
**Blocked input from user:** choice of background MC for π/K-decay /
calo-secondary muon templates (Run 2 used Pythia8 jetjet JZ slices with
HIJING overlay — no Run 3 equivalent identified in the repo).
**Reviewer:** `/review-analysis-code`; fits also `/review-plot`.
**Method reference:** kb `run2_dimuon_note.md` §"Background (fake-muon)
estimate"; full detail in
`IntNotesRun2DimuonReference/tex/DpopAna.tex`.

## Objective

Build the fake-muon purity estimate chain (no code exists yet):

1. **Δp/p = (p_ID − p_MS,extrap)/p_ID** per muon. Prerequisite check:
   are p_ID and p_MS-extrapolated branches in the May 2026 skim? Check
   `SkimCode/source/HFtrigValidation/*/Trig*` (memory:
   reference_skimming_source) and the muon_pair trees. If absent, this
   needs a skim-level addition — surface to user immediately.
2. **Significance:** (Δp/p − μ(pT, q·η))/σ(pT, q·η) with μ, σ for true
   muons from MC (Run 2 used the overlay sample) → makes the cut/
   template shape independent of pT, η, centrality.
3. **Pair significance:** quadrature sum of the two muons.
4. **3-template fit** of the data pair-significance: Sig–Sig, Sig–Bkg,
   Bkg–Bkg templates from MC → signal-pair fraction vs centrality.
   Run 2 result was >98% purity (no cut applied) — we expect to quote a
   purity number, not apply a cut.

## Dummy strategy

- μ, σ (true-muon calibration): from bug-fixed overlay test sample
  (task 03) or pp fullsim test sample — limited stats → coarse (pT, q·η)
  bins, labeled DUMMY.
- Sig–Sig template: same samples (true HF muon pairs).
- Bkg templates: **no Run 3 fake-enriched MC yet** — as a stopgap,
  derive a fake-muon Δp/p shape from reco muons with truth_prob ≤ 0.5
  or no truth match in the overlay sample; flag clearly. Replace when
  user designates a background sample.

## Implementation sketch

- New dir `Analysis/DpopAnalysis/` (or extend SingleBAnalysis):
  - calibration producer: μ/σ maps from MC → `dpop_calibration_*.root`
  - data filler: pair-significance hists per centrality / sign (RDF,
    follow `RDFBasedHistFilling` conventions + hist_binning_map rules)
  - template fitter: TFractionFitter or RooFit extended fit; outputs
    signal fraction vs centrality + fit overlay plots (Run 2 note
    Fig. Templatefits style).
- Per CLAUDE.md, this is a multi-editing-cycle implementation → create
  a dedicated implementation tracking doc with a Physics Procedure
  section (source it from the Run 2 note DpopAna section) BEFORE
  coding.

## Verification

- Significance distribution for MC true muons ≈ unit Gaussian,
  centrality/pT/η independent.
- Closure: fit MC pseudo-data (known mixture) recovers fractions.
