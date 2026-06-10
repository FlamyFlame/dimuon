# Task 04 — Dummy reconstruction-efficiency correction files

**Roadmap:** `analysis_roadmap_2026_06.md` §Q4 / chain [4]
**Depends on:** task 03 (bug-fixed efficiencies). **Reviewer:**
`/review-analysis-code`.

## Objective

Produce *placeholder* single-muon reco-efficiency lookup objects
ε_reco(pT, q·η, centrality) in a file format that the crossx weighting
(task 05) can consume, so the full chain runs end-to-end before the
production fullsim/overlay samples exist. Everything produced here must
be clearly labeled DUMMY (file names + plot watermarks).

## Design constraints

- **Mirror the trigger-efficiency lookup pattern:** task 05 will apply
  w⁻¹ = ε_trig_pair · ε_reco(μ₁) · ε_reco(μ₂); reuse the
  `OpenEffcyPtFitFile` / TH2D-fallback machinery conventions from
  `RDFBasedHistFillingData` so the dummy is a drop-in for the real one.
- **Source choice (preferred order):**
  1. Bug-fixed r17618 overlay 60k-evt sample (task 03): single-muon
     ε(pT, q·η) per centrality bin where populated (only ctr 0–5 and
     5–10 have stats; copy the 5–10 shape to peripheral bins, flagged).
  2. pp fullsim test sample for the pp dummy.
  3. If either is unusable: flat ε from the bug-fixed global numbers
     (overlay ~0.90 single-muon medium; pp ~0.90) — crudest fallback.
- Run 2 note method reference (kb `run2_dimuon_note.md` §Corrections):
  efficiencies binned in (pT, q·η, 10–20% centrality), data-point
  interpolation rather than fits.

## Steps

1. Decide format with task 05 in mind (recommend: TH2D per centrality
   bin, `h_reco_eff_{medium,tight}_ctrX_Y` in
   `reco_eff_dummy_{overlay,pp24}.root` under the respective data dirs).
2. Write a small standalone producer macro in
   `Analysis/EfficiencyCorrs/` (e.g. `make_dummy_reco_eff.C`) reading
   the task-03 histogram outputs.
3. Produce validation plots (ε vs pT per q·η slice, per centrality),
   watermark "DUMMY".
4. Document in roadmap ledger + this file's status line.

## Verification

- ε values in [0,1] everywhere; no empty bins in the fiducial region
  (pT 4–40 GeV, |q·η| < 2.4) after fallback filling.
- Lookup round-trip test: sample 100 (pT, q·η, ctr) points, confirm the
  task-05 reader returns the histogram values.

## Replacement plan

When production samples land (ATLHI-666 pp24 fullsim, ATLHI-576-line
overlay production): rerun overlay/pp pipelines, regenerate this file
from full stats with the same producer, drop the DUMMY label. The
task-05 weighting code must not need changes.
