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
//   (gap) FIDUCIAL GAP CUT from ParamsSet::single_mu_fiducial_gap_cuts -- built from ParamsSet so
//         the windows are never retyped into a JIT string.
// A PAIR is kept only if BOTH legs satisfy each of them.
//
// ⚠ MIRROR NOTICE. RDFBasedHistFilling/FillMCTrigEffHists.cxx builds the identical strings inline
// (kTruthFidPair / kFwdVetoPair / kGapPair, ~lines 599-663 and the Step-3 block ~1109-1115). That
// file is the ORIGINAL; this header was extracted from it verbatim. It was deliberately NOT
// migrated onto this header at extraction time because a concurrent session was re-running that
// macro and an edit mid-run would have broken it. MIGRATE IT the next time it is touched, so the
// two constructions cannot drift -- a drift here is exactly the class of silent error the header
// exists to prevent.
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
         + ParamsSet::FiducialGapCutExpr("m2_charge * m2_eta");
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
// RDFBasedHistFilling/RDFBasedHistFillingPP.cxx:434 (and PbPb.cxx:977) with the m1./m2. dots
// replaced by the underscore aliases. Ground truth for the values: docs/analysis_overview.md §2
// and docs/signal_selection_change_impact.md §0. There is NO dR cut (removed 2026-06-22).
inline std::string SingleBSignalCutsReco()
{
    return "minv > 1.08 && minv < 2.9 && pair_pt > 8 && "
           "m1_charge * m1_eta < 2.2 && m2_charge * m2_eta < 2.2";
}

}  // namespace MCTrigEffPairSel

#endif  // MC_TRIG_EFF_PAIR_SELECTION_H
