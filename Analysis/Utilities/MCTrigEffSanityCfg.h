// MCTrigEffSanityCfg.h
//
// SINGLE SOURCE OF TRUTH for the Step-1 SANITY-CHECK truth-reco pT-match threshold
// (mc_trigger_efficiency.md §3.5).
//
// Why it is a shared header and not a literal in each file: the threshold is APPLIED by
// RDFBasedHistFilling/FillMCTrigEffHists.cxx (do_sanity mode) and PRINTED ON THE CANVAS by
// plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff.cxx. A plot must state the numerical
// value of every cut the audience needs -- and a value retyped in the plot macro would drift
// silently the first time the cut is changed, putting a WRONG number in front of the reader
// with no error anywhere.
//
// It is NOT part of the nominal MC trigger-efficiency selection: it exists only to test
// whether badly measured muons drive the forward MC/data disagreement.

#ifndef MC_TRIG_EFF_SANITY_CFG_H
#define MC_TRIG_EFF_SANITY_CFG_H

// |truth p_T - reco p_T| / truth p_T < kSanityPtMatchThr
inline constexpr double kSanityPtMatchThr = 0.10;

#endif // MC_TRIG_EFF_SANITY_CFG_H
