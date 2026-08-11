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
