# _sub_iter2fix_1 — round-9 MC trig-eff plots, iteration-2 presentation fixes

Subagent scratch doc (append-only). Orchestrator owns git and all shared docs.

## Mandate
Fix five PRESENTATION findings (G1 CRITICAL, G2/G3 WARNING, G4/G5 INFO) from the
iteration-2 plot review. **No fit function / fit range / parameter limit / plateau
window / usable-fit_ok logic may change — no fitted number may move.**
The sixth finding (nominal `expo` cannot describe the shape) is the orchestrator's
physics decision, explicitly out of scope here.

Editable files (exactly three):
- plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff.cxx
- plotting_codes/trig_effcy/mc_based/plot_dr_correction_fits.cxx
- plotting_codes/trig_effcy/mc_based/fit_dr_corrections.cxx

## Done =
1. G1 legend on the 16 `*_pair_eta_pt.png` canvases readable (2 rows, taller header strip,
   no entry touching a neighbour or the right edge) — verified by LOOKING at the PNG.
2. G2 inclusive canvas + inclusive ratio canvas get their own y range from the inclusive
   points; nine-panel pair-pT files keep the shared directory-wide range (commented why).
3. G3 fit_report pair-pT labels `%.1f`; 36 pp reports regenerated.
4. G4 off-scale note on the fit canvases carries the eta^pair range and the sign.
5. G5 the nine overlay panels share one y range (first pass over all cells), 3.0 cap +
   arrows kept.
6. Regenerated: plot_mc_trig_eff x {pp_full, overlay} x {tight, medium}; the dR fit
   pipeline for pp_full, both WPs, steps 3+4, methods expo/polyu_fixedRp/interp,
   SKIP_MEASURE=1.
7. Fit-report diff vs baseline shows ONLY the %.1f pair-pT label change.

## Step 0 (done) — baseline captured
36 `fit_report*.txt` under
  plots/pp_trigger_efficiency/mc_based{,_medium}/step{3,4}_dr_fit/{expo,polyu_fixedRp,interp}/
copied to
  <scratchpad>/iter2_baseline/  (36 files verified)

## Step 1 (planned) — source edits

## Step 1 (done) — source edits, all three files compile

**G3** `fit_dr_corrections.cxx` ~625: `ptlab` Form `"[%.0f,%.0f)"` -> `"[%.1f,%.1f)"` (+ comment).
Label only; no fitted quantity reads this string.

**G2** `plot_dr_correction_fits.cxx`
- the single `g_ylo/g_yhi` block replaced by a `compute_range(bool incl_only)` lambda; called
  twice -> `grid_rng` (all cells, shared by every nine-panel pair-pT PNG in the method/mode dir)
  and `incl_rng` (the inclusive cell alone). Same screening (`DrCorrPlateauUsable`), same
  central-values-only rule, same annotation-band stretch. Comment states why the two differ.
- inclusive canvas now draws with `incl_rng`.
- ratio: new `ratio_range(bool incl_only)` -> `rgrid` / `rincl`; `r_lo/r_hi` set to `rincl`
  before the inclusive-ratio canvas.
- stdout now prints both ranges.

**G4** `plot_dr_correction_fits.cxx`
- `draw_cell`: new `cell_id` ("#eta^{pair} #in [a,b), ", empty on the inclusive canvas) and
  `sign_id` (`s.legend + ", "` when sepr); both prefixed onto the off-scale entry.
- `draw_ratio_cell`: same `cell_id` prefix (no sign tag -- the figure IS the sign ratio).
- `draw_header` / `ratio_header` gained an `n_show` argument (3 on the 1560 px grid canvas,
  2 on the 900 px inclusive canvas); separator ", " -> ";  ".

**G1 + G5** `plot_mc_trig_eff.cxx`, both nine-panel overlay blocks (Step 3 and Step 4)
- G5: first pass over ALL (pair pT, pair eta) cells -> ONE `ymax` per canvas, from CENTRAL
  values only (round-9 rule); 3.0 cap and per-point arrows kept.
- G1: reserved header strip `kHeaderPx = 128`, grid moved into a `TPad(0,0,1,1-hfrac)`;
  header laid out in canvas PIXELS (`ny(px) = 1 - px/canv_h`): title at 26 px (22 px text),
  legend box 38-100 px with `SetNColumns(4)` (two rows of four, 16 px text), off-scale record
  at 118 px (15 px text).
- new: canvas-level off-scale record (eta cell + pair-pT cell + dR + value, 2 shown + total),
  which these canvases previously did not have at all.

## Step 2 (planned) — regenerate + verify

## Step 2 (done) — two extra fixes found by LOOKING at the rendered PNGs

**2a. Panel label vs off-scale arrow (found in the first G1/G5 render).** On the PbPb
`step3_eps_dr_zoom_pair_eta_pt.png` the in-frame eta label at pad NDC (0.17, 0.86) was struck
through by the off-scale arrows, which are drawn just under the frame top. Fixed the same way
the eps_dR distribution canvases do it: `gPad->SetTopMargin(0.10)` per panel and the label moved
to (0.15, 0.945), i.e. into a reserved strip ABOVE the frame. Both Step 3 and Step 4.

**2b. Minimum-span fallback could CLIP the inclusive data (found from the printed ranges).**
`if (yhi - ylo < 0.2) { ylo = 0.9; yhi = 1.12; }` jumped to a FIXED window. Harmless for the
directory-wide range (its span is always > 0.2) but wrong for the new inclusive range: the
Step-4 inclusive points span 0.99-1.16, span 0.187 < 0.2, so every Step-4 inclusive canvas got
[0.9, 1.12] and the 1.16 maximum -- the peak of the rise the figure exists to show -- would have
been arrowed off scale. Replaced by a widen-around-the-data form (mid +- 0.10), in both
`compute_range` and `ratio_range`. Step-4 inclusive is now [0.9608, 1.1608] (tight) /
[0.9619, 1.1619] (medium).

**2c. Ratio header line spacing.** With the eta cell now in the off-scale string, the
`eta^{pair}` superscript at ny(78) was drawn into the `eps^{2mu4}_{dR}` subscript of the
definition line at ny(56). Off-scale line -> ny(86), `kRatioHeaderPx` 92 -> 104.

**2d. Diagnostic (stdout only, not on canvas).** Each `wrote <png>` line now also reports the
number of points pushed off scale, so an off-scale canvas can be found from the pipeline log
instead of by opening PNGs one at a time.

## Step 3 (done) — regeneration

- `plot_mc_trig_eff.cxx+("pp_full"|"overlay", true|false)` -- 4 runs, all exit 0. 220 PNGs
  rewritten across the 4 output trees.
- `bash pipelines/run_dr_correction_fits.sh` with
  `SAMPLES=pp_full WPS="tight medium" STEPS="3 4" METHODS="expo polyu_fixedRp interp"
  SKIP_MEASURE=1` -- exit 2 (the pre-existing plateau-guard failure, EXPECTED); **no ARTEFACT
  FAILURES**. This run rewrote the 36 fit reports (G3).
- Two later plot-only passes (`SKIP_FIT=1` in addition) for 2b/2c -- the fit ROOT files and the
  reports were not touched by those.
- 324 dR-fit PNGs regenerated under pp_trigger_efficiency/mc_based{,_medium}/step{3,4}_dr_fit/
  (81 per step per WP: 3 methods x [9 sign_intgr + 9 sign_sepr + 9 sign_sepr ratio] + the
  inclusive/inclusive-ratio canvases). 24 inclusive + 12 inclusive-ratio canvases.

## Step 4 (done) — verification

**Fit-report diff (the no-numbers-may-move constraint).** All 36 baseline reports were rewritten
with the 8 old pair-pT tokens mapped to the new ones
([8,12)->[8.0,11.5), [12,17)->[11.5,16.6), [17,24)->[16.6,24.0), [24,35)->[24.0,34.6),
 [35,50)->[34.6,50.0), [50,72)->[50.0,72.1), [72,104)->[72.1,104.0), [104,150)->[104.0,150.0)),
whitespace-normalized, and compared line by line to the regenerated files:
**36 / 36 identical; 0 files differ by anything other than the pair-pT label. 2895 label tokens
rewritten (row labels AND the "parameters pinned ON a fit limit" cell list).** No plateau, no
fitted parameter, no chi2/ndf, no fit_ok count, no cell count moved.

**Rendered PNGs inspected (Read tool):**
| PNG | finding | verdict |
|---|---|---|
| pp mc_based step3_eps_dr_zoom_pair_eta_pt | G1, G5 | legend = 2 rows x 4, no overlap, nothing clipped; all 9 panels share y 0-1.75 |
| pbpb mc_based step3_eps_dr_zoom_pair_eta_pt | G1, G5, 2a | legend clean; shared 0-3.0; labels above the frames, arrows no longer through them; off-scale note "eta in [-2.4,-2.0), pT in [8.0,11.5) GeV, dR = 0.12: 4.9; ... (15 in total)" |
| pbpb mc_based step4_eps_dr_single_zoom_pair_eta_pt | G1, G5 | as above, 4 off-scale points listed |
| pp mc_based step4_eps_dr_single_zoom_pair_eta_pt | G1, G5 | shared 0-2.4, legend clean |
| pp step3 expo sign_intgr inclusive | G2 | frame 0.80-1.15, points fill it, the small-dR rise is the whole figure |
| pp step4 polyu_fixedRp sign_sepr inclusive | G2, 2b | frame 0.96-1.16, the ~15% rise fully visible, no clipping |
| pp step4 polyu_fixedRp sign_sepr inclusive_ratio | G2 | frame 0.88-1.07 (was 0-3.0) |
| pp step3 expo sign_sepr inclusive_ratio | G2 | frame 0.40-1.13 (was 0-3.0) |
| pp step3 expo sign_sepr pairpt_8.0_11.5 | G4 | header: "off scale: eta^{pair} in [2.0,2.4), same sign, dR=0.03: 4.92 +- 4.92" -- panel AND sign named; matches the blue arrow in the [2.0,2.4) panel |
| pp step3 expo sign_sepr pairpt_72.1_104.0_ratio | G4, 2c | 3 entries each named by eta cell + "(4 total)"; header lines no longer overlap |
| pp step4 polyu_fixedRp sign_intgr pairpt_8.0_11.5 | regression check | nominal nine-panel unchanged, shared directory range 0-2.63 as before |

## NOT done / out of scope
- The 12 sign-integrated `fit_report.txt` under `pbpb_trigger_efficiency/mc_based{,_medium}/`
  still carry `%.0f` pair-pT labels: the orchestrator's regeneration instruction scoped the fit
  stage to `SAMPLES=pp_full`. The code fix is global, so one
  `SAMPLES=overlay ... SKIP_MEASURE=1` pass will fix them; it was not run here.
- Finding 6 (the `expo` shape failure) untouched, as instructed. No fit function, fit range,
  parameter limit, plateau window or usable/fit_ok logic was modified in any of the three files.
