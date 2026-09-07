# Plot Review Log
**Task**: Follow-up to the previously-APPROVED linear-y-scale change
(review-plot-20260906-212055-qeta-linear-yscale.md). Two new user requests:
1. `plot_muon_q_eta_pt_dependence.cxx` (F14 p_T-dependence diagnostic, previously out of
   scope) also switched from log-y to linear-y.
2. The prior fix for the header/axis-multiplier-box collision (forcing full 7-digit y-axis
   labels via `TGaxis::SetMaxDigits(7)`) is REPLACED: the user wants ROOT's native "x10^n"
   scientific-notation axis labels back (not full digits), so the collision is now fixed by
   repositioning the `Header()` text instead of suppressing the notation.

**Log file**: review-plot-20260906-224247-qeta-ptdep-linear-and-header-fix.md
**Started**: 2026-09-06 22:42:47
**Status**: APPROVED at iteration 1
**Iterations completed**: 1
**Max iterations**: 5
**Summary**: plot_muon_q_eta_pt_dependence.cxx converted to linear-y (matching the sibling
macros) and a latent pow(10,...) bug in its DrawGapBands fixed along the way. The prior
digit-forcing fix in plot_muon_q_eta_spectrum.cxx was reverted in favor of ROOT's native
x10^n notation, with the underlying Header/axis-box collision fixed by repositioning Header()
instead. Reviewer independently reran all three macros and matched every number exactly;
zero CRITICAL/WARNING issues.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Numerical verification**: All 8 quantities (Tight/Medium fiducial fail fractions for pp24
and PbPb combined/50-80%, all 4 p_T-slice fractions for Tight, all 4 for Medium, truth ε_acc
pp/PbPb overlay) independently reproduced via recompile+rerun, exact match.

## Execution

Files modified:
- `plotting_codes/single_b_analysis/plot_muon_q_eta_spectrum.cxx` (revert digit-forcing fix,
  fix Header() offsets instead)
- `plotting_codes/single_b_analysis/plot_muon_q_eta_pt_dependence.cxx` (new: log-y -> linear-y)

### plot_muon_q_eta_spectrum.cxx (revision of the already-approved change)
- Removed `#include <TGaxis.h>` and `TGaxis::SetMaxDigits(7)`.
- Reverted `StyleSpectrum`'s y-axis `SetTitleOffset` 3.00 -> 1.70 (original value).
- Reverted `SetLeftMargin` 0.28 -> 0.20 in `MakeSinglePad` and both `p_main`/`p_rat` of
  `MakeRatioPads` (original value).
- Root-cause fix for the actual collision: `Header()`'s near-frame text offsets raised
  0.025->0.055 (empty-sub case) and 0.020->0.050 (sub-line, two-line case), with the
  corresponding main line pushed 0.067->0.095 to keep line spacing. This clears ROOT's
  native "x10^n" multiplier box (now visible again, e.g. "x10^3" over "3500" ticks) without
  needing full-digit labels.
- Verified by cropping/inspecting: old_gap_cuts... wait, old_gap_cuts was deleted per prior
  request; verified instead on new_gap_cuts/medium (highest magnitude, ~3.5M plateau,
  Medium WP PbPb combined) and new_gap_cuts/medium ctr_binned (tightest 2x3 pad grid) --
  both clean, headers/subtitles fully clear of the "x10^3" box on every panel.

### plot_muon_q_eta_pt_dependence.cxx (new conversion)
- Removed `gPad->SetLogy()` (single call site, in the 6-panel draw loop).
- Removed `SetMoreLogLabels()` + its comment from `StyleShape`.
- Replaced `SetCommonLogRange`/`kMaxDecades` (returned an int n_below) with
  `SetCommonLinRange` (floor 0, ceiling kHeadroom=1.15x peak, void return) -- same pattern as
  the sibling macro. Removed `NoteOffScale` function and its two call sites (main per-panel
  call + the canvas-level "N bins below the shared Pb+Pb axis range" block), since nothing
  goes off-scale on a 0-floor linear axis.
- Fixed `DrawGapBands`, which had HARD-CODED `std::pow(10., gPad->GetUymin()/GetUymax())`
  (assumed log unconditionally, unlike the sibling macro's `DrawGapMarkers` which branches on
  `gPad->GetLogy()`) -- now uses the raw NDC-to-user values directly since this macro is
  linear-only (no log/linear toggle argument).
- This macro's content is UNIT-AREA PDFs (Scale to Integral("width")=1), so values are
  O(0.001-0.45) -- never triggers ROOT's axis-multiplier notation, so no Header/margin/offset
  changes were needed here (confirmed by inspection: no "x10^n" box appears on any panel).

Regenerated PNGs (all overwritten in place):
- `muon_gap_cuts/new_gap_cuts/muon_q_eta_spectrum_pp24.png` + 3 siblings (Tight)
- `muon_gap_cuts/new_gap_cuts/medium/` same 4 (Medium)
- `muon_gap_cuts/new_gap_cuts/muon_q_eta_pt_dependence.png` (Tight)
- `muon_gap_cuts/new_gap_cuts/medium/muon_q_eta_pt_dependence_medium_wp.png` (Medium)
- `muon_gap_cuts/truth_q_eta/` all 3 (rerun to confirm unaffected, byte-identical to F12)
`old_gap_cuts/` was NOT regenerated (deleted per prior user request; the macro's
`use_new_gap_cuts=false` code path is untouched and would still work if ever needed again).

Console numbers verified against tracking doc F11a/F14 (unchanged from before -- this is a
y-axis-presentation-only change):
- pp24/PbPb proposed fiducial (Tight): 0.0379 / 0.0560 -- MATCH F11a.
- pp24/PbPb proposed fiducial (Medium): 0.0390 / 0.0647 (50-80%) -- MATCH F11a/F14.
- p_T-dependence per-slice fractions (Tight): pp24 0.0368/0.0329/0.0301/0.0425/0.0379 (all pT)
  -- MATCH F14 exactly.
- Truth ε_acc: 0.9133 (pp) / 0.9121 (PbPb overlay) -- MATCH F12.
