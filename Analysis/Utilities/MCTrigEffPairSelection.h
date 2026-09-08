#ifndef MC_TRIG_EFF_PAIR_SELECTION_H
#define MC_TRIG_EFF_PAIR_SELECTION_H

#include <string>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"

// =================================================================================================
// THE MC trigger-efficiency PAIR selection (docs/tracking/mc_trigger_efficiency.md §3.0, §3.5),
// as one builder, so a CONSUMER of the dR correction selects exactly the sample the correction was
// measured on.
//
// WHY IT EXISTS. The Step-3 dR correction eps_dR(dR) is only defined for the pairs it was measured
// from. A closure test (or any future consumer) that rebuilds "roughly the same" selection would
// apply the correction to a different population -- silently, since every histogram still fills.
//
// The three requirements below are NOT interchangeable and none of them is cosmetic:
//   (c')  TRUTH FIDUCIAL  truth pT > 4, |truth eta| < 2.4  -- removes muons that enter the reco
//         fiducial region only through mismeasurement (§3.0(c')).
//   (3.5) FORWARD LOW-pT VETO  pT > 7 || q*eta > -2  -- the r16578 production's forward-negative
//         turn-on is saturated (R8/R10/R14, a real simulation problem, not bad muons), so those
//         muons are out of the dR-correlation measurement. ASYMMETRIC: only the corner that is
//         SIMULTANEOUSLY low-pT and forward-negative.
//   (gap) FIDUCIAL GAP CUT from ParamsSet::single_mu_fiducial_gap_cuts on BOTH legs, AND the
//         PAIR-LEVEL window |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2 (user,
//         2026-09-07) -- all built from ParamsSet so no number is retyped into a JIT string.
// A PAIR is kept only if BOTH legs satisfy each of them.
//
// ⚠ MIRROR NOTICE (partly resolved 2026-09-07). RDFBasedHistFilling/FillMCTrigEffHists.cxx is the
// ORIGINAL; this header was extracted from it verbatim. `kGapPair` there now CALLS
// FiducialGapCut() below, so the gap construction can no longer drift. STILL DUPLICATED inline in
// that file: kTruthFidPair and kFwdVetoPair (and the single-muon / lg_-ot_ leg-alias variants,
// which this header does not provide at all). Migrate those the next time they are touched -- a
// drift there is exactly the class of silent error this header exists to prevent.
// =================================================================================================
namespace MCTrigEffPairSel {

// Column-name convention of the Step-3 pair node: m1_pt, m1_eta, m1_charge, m1_wp,
// m1_truth_pt, m1_truth_eta (and m2_*), all created by Alias on the pair tree.
inline std::string TruthFiducial()
{
    return "m1_truth_pt > 4 && fabs(m1_truth_eta) < 2.4 && "
           "m2_truth_pt > 4 && fabs(m2_truth_eta) < 2.4";
}

inline std::string ForwardLowPtVeto()
{
    return "(m1_pt > 7 || m1_charge * m1_eta > -2) && "
           "(m2_pt > 7 || m2_charge * m2_eta > -2)";
}

inline std::string FiducialGapCut()
{
    return ParamsSet::FiducialGapCutExpr("m1_charge * m1_eta") + " && "
         + ParamsSet::FiducialGapCutExpr("m2_charge * m2_eta") + " && "
         + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
}

// The complete Step-3 pair selection. `apply_gap_cut` mirrors the MCTRIGEFF_NO_GAPCUT escape hatch
// of the fill macro (nominal: ON).
inline std::string Step3PairSelection(bool apply_gap_cut = true)
{
    std::string sel = "m1_wp && m1_pt > 4 && fabs(m1_eta) < 2.4 && "
                      "m2_wp && m2_pt > 4 && fabs(m2_eta) < 2.4 && " + TruthFiducial()
                    + " && " + ForwardLowPtVeto();
    if (apply_gap_cut) sel += " && " + FiducialGapCut();
    return sel;
}

// The DATA-LIKE single-b RECO signal selection, byte-identical to `signal_cuts` in
// RDFBasedHistFilling/RDFBasedHistFillingPP.cxx with the m1./m2. dots replaced by the underscore
// aliases. NOT PbPb: since 2026-08-18 the PbPb signal region is still the retired one-sided
// `q*eta < 2.2` and does NOT match this (docs/signal_selection_change_impact.md §0).
// Ground truth for the values: docs/analysis_overview.md §2
// and docs/signal_selection_change_impact.md §0. There is NO dR cut (removed 2026-06-22).
// 2026-08-17: the per-muon one-sided `q*eta < 2.2` was REPLACED by the detector-gap fiducial
// cut on BOTH muons (FiducialGapCut() above) -- kept in lockstep with the crossx by
// construction, since both read ParamsSet::single_mu_fiducial_gap_cuts. 2026-09-07: the
// PAIR-LEVEL |eta^pair| < 2.2 window joined it, inside the same FiducialGapCut() helper, so
// the lockstep holds for it too.
inline std::string SingleBSignalCutsReco()
{
    return "minv > 1.08 && minv < 2.9 && pair_pt > 8 && "
         + FiducialGapCut();
}

}  // namespace MCTrigEffPairSel

#endif  // MC_TRIG_EFF_PAIR_SELECTION_H
