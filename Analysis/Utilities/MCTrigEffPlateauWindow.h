#ifndef MC_TRIG_EFF_PLATEAU_WINDOW_H
#define MC_TRIG_EFF_PLATEAU_WINDOW_H

// =================================================================================================
// THE LARGE-dR PLATEAU WINDOW -- SINGLE SOURCE OF TRUTH
//
// eps_dR^cross (Step 3, mc_trigger_efficiency.md §3.3) and eps_dR^single (Step 4, §3.4) are
// RATIOS normalized to their large-dR plateau: the plateau IS the "well-separated muons"
// reference, so it must be measured over a dR range where the curve is actually FLAT
// (diagnostic 1 of both sections).
//
// Every stage that measures, normalizes by, reports or draws that window MUST read the edges
// from here. They used to be retyped in four places (plot_mc_trig_eff.cxx,
// plot_mc_trig_eff_corrected.cxx, FillMCTrigEffHists.cxx, and the fit-report text in
// fit_dr_corrections.cxx); on 2026-08-04 the nominal window changed and three of those copies
// silently kept the old value -- every plot still rendered and every number still came out,
// which is exactly the failure mode .claude/CLAUDE.md §Binnings exists to prevent.
//
// -------------------------------------------------------------------------------------------
// NOMINAL WINDOW: dR in [2.0, 3.5]  (user decision, 2026-08-04). Both edges are cut, and each
// edge is cut for a DIFFERENT, independently measured reason:
//
//   LOWER edge 1.0 -> 2.0 -- PER-CELL structure. In several (pair pT, pair eta) cells the
//   dR in [1,2] half sits significantly above the far half: worst pT_pair[8,14) x
//   eta_pair[-0.5,0.5) at +0.054 (7.9 sigma). Cutting it improves the mean per-cell
//   constant-fit chi2/ndf from 1.64 to 1.27 (Step 3) and 1.66 to 1.34 (Step 4). Cutting only
//   the upper edge does NOT fix this (per-cell 1.68 / 1.70, no better than the retired window).
//
//   UPPER edge 4.0 -> 3.5 -- the INCLUSIVE tail is an ARTEFACT, not trigger physics. The
//   inclusive pp24 Tight curve is flat from dR ~0.6 to ~3.4 (0.9702-0.9792) and then falls
//   monotonically: 0.9571 (3.62), 0.9521 (3.88), 0.9341, 0.9262, 0.9156, 0.8965, 0.668, 0.390.
//   The cause is geometric: with dPhi <= pi, dR > 3.5 forces |dEta| > 1.54 and dR > 4 forces
//   |dEta| > 2.47, so large-dR pairs push BOTH legs into the endcaps -- precisely where eps_MC
//   is documented as badly parameterized (the r16578 forward-endcap anomaly, R8/R10/R14).
//   §3.3 diagnostic 1 names this exact failure mode ("failure to flatten means residual
//   kinematic mis-parameterization of eps_MC leaking in"). Run 2's analogous L1 close-by-RoI
//   correction stays at unity out to large dR and shows no such fall.
//   Inclusive constant-fit chi2/ndf, pp24 Tight Step 3 / Step 4:
//       [1,4] 3.67 / 3.85      [2,4] 5.06 / 5.37      [2,3.5] 1.08 / 0.76
//   -- i.e. the retired [1,4] and the interim [2,4] were both dominated by that tail.
//
//   COST: per-cell statistical errors grow. Median relative stat error on the plateau,
//   pp24 Tight: Step 3 0.95% -> 1.40%, Step 4 0.46% -> 0.66%; worst cell 3.6% -> 6.6%
//   (Step 3). All 36 cells stay statistically meaningful.
// -------------------------------------------------------------------------------------------
namespace MCTrigEffPlateau {

constexpr double kLo = 2.0;
constexpr double kHi = 3.5;

// The RETIRED [1,4] window, kept ONLY to size the plateau-window normalization systematic
// |plateau[nominal] - plateau[retired]|, written per cell into the plateau ROOT file and
// reported in plateau_guard_report.txt (docs/systematic_uncertainties.md §1a).
// It is an UNCERTAINTY: nothing may ever normalize a curve by it, and it is never a second
// nominal.
constexpr double kSystLo = 1.0;
constexpr double kSystHi = 4.0;

}  // namespace MCTrigEffPlateau

#endif  // MC_TRIG_EFF_PLATEAU_WINDOW_H
