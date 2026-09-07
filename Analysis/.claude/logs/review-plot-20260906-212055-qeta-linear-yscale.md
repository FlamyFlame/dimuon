# Plot Review Log
**Task**: Switch the muon q*eta distribution plots (data spectrum + truth spectrum) from log-y
to LINEAR y scale, overwriting the existing PNGs in place, so the gap depth and steepness are
visible on an absolute scale.

Files to modify:
- Analysis/plotting_codes/single_b_analysis/plot_muon_q_eta_spectrum.cxx
- Analysis/plotting_codes/single_b_analysis/plot_muon_truth_q_eta_spectrum.cxx

Out of scope: plot_muon_q_eta_pt_dependence.cxx (F14, separate diagnostic).

**Log file**: review-plot-20260906-212055-qeta-linear-yscale.md
**Started**: 2026-09-06 21:20:55
**Status**: APPROVED at iteration 1
**Iterations completed**: 1
**Max iterations**: 5
**Summary**: Log-y to linear-y switch on both q*eta spectrum macros (data + truth), 19 PNGs
regenerated, numbers byte-identical to F6/F9/F11a/F12. Reviewer independently reran both
macros and matched every quoted number exactly; only finding was one INFO-level stale code
comment (said margin 0.25, actual 0.28), fixed post-approval (comment-only, no rebuild needed).

## Execution

Files modified:
- `plotting_codes/single_b_analysis/plot_muon_q_eta_spectrum.cxx`
- `plotting_codes/single_b_analysis/plot_muon_truth_q_eta_spectrum.cxx`

Changes (both files):
1. Removed every `gPad->SetLogy()` call and the `SetMoreLogLabels()` line/comment in
   `StyleSpectrum()`.
2. Replaced `SetCommonLogRange`/`kMaxDecades` + `NoteOffScale` (log-floor-with-cap machinery,
   now dead on a linear axis) with `SetCommonLinRange(hs)`: floor pinned at 0, ceiling at
   `kHeadroom=1.15`x the overlaid peak. All call sites updated; `NoteOffScale` calls dropped
   (nothing goes off-scale on a 0-floor linear axis).
3. `DrawGapMarkers`/`AutoLegendBox` unchanged -- they already branch on `gPad->GetLogy()` and
   take the linear path automatically now that `SetLogy()` is never called.

Regression found and fixed during execution (data macro only, `plot_muon_q_eta_spectrum.cxx`):
switching to linear-y made ROOT invoke its "x10^n" axis-multiplier box (needed because raw
counts now show as full numbers up to ~3.5M), which collided with the `Header()`/`CanvasGapKey()`
text placed just above the frame (present on every panel with an empty-sub `Header()` call, and
on the sub-line of two-line headers). Fixed at the root (not by fragile per-panel offset hacks):
- Added `TGaxis::SetMaxDigits(7)` (plus `#include <TGaxis.h>`) so ROOT prints full-digit labels
  instead of the multiplier box, for any value up to 7 digits (max observed ~3.5M).
- That made the full-digit labels (e.g. "3000000") collide with the y-axis title instead, so
  `StyleSpectrum`'s `SetTitleOffset` was raised 1.70 -> 3.00, and `MakeSinglePad`/`MakeRatioPads`
  `SetLeftMargin` raised 0.20 -> 0.28 (both `p_main` AND `p_rat` in `MakeRatioPads`, kept equal
  so the stacked main+ratio pads' x-axes still line up).
- Verified clean (no overlap) across all 4 flag combinations (Tight/Medium x old/new gap cuts)
  on both the 3-column canvas (Plot A/B) and the tighter 2x3 grid (Plot C), by cropping and
  visually inspecting the worst-case (highest-magnitude, PbPb combined Medium WP, ~3.5M plateau)
  and the tightest-pad case (Plot C's 2-column grid).
- `plot_muon_truth_q_eta_spectrum.cxx` values stay at O(10-50) nb (weighted cross-section, not
  raw counts) -- never triggers the digit-count issue, so it needed no margin/offset/TGaxis
  changes.

Regenerated PNGs (all overwritten in place, numbers cross-checked against the tracking doc
before/after -- byte-identical to F6/F9/F11a/F12):
- `muon_gap_cuts/old_gap_cuts/` (4) + `medium/` (4) = 8
- `muon_gap_cuts/new_gap_cuts/` (4) + `medium/` (4) = 8
- `muon_gap_cuts/truth_q_eta/` (3)
Total 19 PNGs. `muon_q_eta_pt_dependence.png` (F14) untouched, per scope.

Console numbers verified against tracking doc:
- pp24 legacy gap-cut fail fraction: 0.0739 (F6/F9) -- MATCH.
- pp24/PbPb proposed fiducial: 0.0379 / 0.0560 Tight, 0.0390/0.0640 Medium (F11a) -- MATCH.
- Truth ε_acc: 0.9133 (pp) / 0.9121 (PbPb overlay) (F12) -- MATCH.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 1 (INFO only)
**Details**:
1. [DOCS/COMMENT] Stale margin value in a code comment (plot_muon_q_eta_spectrum.cxx:248 said
   0.25, actual value 0.28). Severity: INFO. Fixed post-approval (comment text only).
**Numerical verification**: Reviewer independently recompiled and reran BOTH macros end-to-end
(not just cross-referenced doc text) and matched every quoted number exactly: pp24 legacy
0.0739 (672686/9107895); pp24/PbPb proposed fiducial Tight 0.0379/0.0560; per-centrality Tight
0.0555/0.0565/0.0563/0.0564/0.0576; truth ε_acc 0.9133 (pp) / 0.9121 (PbPb overlay) with
per-window breakdowns matching F12 exactly.
