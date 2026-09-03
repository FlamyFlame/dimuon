// PtHatKn45ProjectedStats.h
//
// SINGLE SOURCE OF TRUTH for the "projected statistics" decision aid (2026-09-03 user request;
// docs/tracking/pythia_pp24_pthat_stats_projection.md): what would the pp24 Pythia FullSim FULL
// sample look like if the two highest pT-hat slices (kn4 = 70-125 GeV, kn5 = 125-300 GeV) had
// N_TARGET events instead of the N_CURRENT they were actually produced with. Shared by the
// projected crossx plot (plot_pythia_fullsim_kn_pt_crossx.cxx) and the projected SS/OS pair
// statistics CSVs (write_mc_pair_statistics_tables_projected.cxx) so the two consumers can never
// disagree about the scale factor.
//
// N_CURRENT is MEASURED, not assumed: both slices load exactly 319999 events, confirmed in the
// FULL-sample production log
// (/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/trigeff_rerun_tight_20260720_221052.log,
// "Fullsim pTH70_125/pTH125_300 beam=pp N=319999/319999" -- the InitInputFullsim/ProcessDataHook
// event count the per-pair weight `w = sigma_slice*genFiltEff*r_isospin/N_slice` is built from).
// Matches the user's "320k".
//
// PHYSICS (see the tracking doc's Physics Procedure §1): w scales as 1/N while the pair count
// scales as N, so sum(w) -- a cross section -- is an UNBIASED estimator of the slice's cross
// section independent of N; only its STATISTICAL UNCERTAINTY shrinks, as 1/sqrt(N). So a
// "projection" NEVER rescales a cross-section-like central value; it rescales only a
// statistics-carrying quantity (a bin error, or a raw pair count / sum-of-weight-squared).

#ifndef PT_HAT_KN45_PROJECTED_STATS_H
#define PT_HAT_KN45_PROJECTED_STATS_H

namespace PtHatKn45Projected {

inline constexpr double kNCurrent = 319999.;   // measured, both slices (see header comment)
inline constexpr double kNTarget  = 1.2e6;     // the user's requested target, both slices

// N_target / N_current: raw-count / sumw2 scale-UP factor; error/sqrt(sumw2)-equivalent
// quantities scale DOWN by this factor's sqrt (or divide by it, for a variance-like quantity).
inline constexpr double kSf = kNTarget / kNCurrent;

}  // namespace PtHatKn45Projected

#endif  // PT_HAT_KN45_PROJECTED_STATS_H
