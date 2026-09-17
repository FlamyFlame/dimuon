#ifndef MC_TRIG_EFF_PAIR_SELECTION_H
#define MC_TRIG_EFF_PAIR_SELECTION_H

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
// The four requirements below are NOT interchangeable and none of them is cosmetic:
//   (wp)  NOMINAL RECONSTRUCTED-MUON DEFINITION, taken at PAIR level from the ntuple-processing
//         flag `pair_pass_tight` / `pair_pass_medium` (aliased to `pair_wp` by the caller, so
//         this header stays WP-agnostic). That flag is
//             pair_pass_X = m1.pass_X && m2.pass_X && ip_pair_ok
//         (PythiaFullSimExtras.c:526-527), i.e. BOTH legs' combined/WP/IDCuts/MuonCuts bits,
//         |eta| < 2.4, pT > 4.5 GeV, one-sided dP/P < 0.12 and d0/z0 -- AND, for pp, the
//         SAME-VERTEX pair requirement (both muons within |d0| < 2 mm, |z0 sin(theta)| < 2 mm of
//         ONE common, possibly secondary vertex; Utilities/AllVertexIPSelection.h,
//         docs/tracking/pp24_all_vertex_pairs.md). The per-muon flags are deliberately loosened
//         to ANY-vertex, so selecting on `m1_wp && m2_wp` would MISS the same-vertex requirement
//         that the data applies -- the eff would be measured on a looser population than the one
//         it is applied to. `ip_pair_ok` is identically true when UseAllVertexIP() is false
//         (HIJING overlay, noovl), so the pair form is an EXACT no-op for those samples.
//         The literal `pT > 4.5 && |eta| < 2.4` terms below RESTATE two cuts already inside the
//         flag; they are kept so the selection string is self-describing.
//   (c')  TRUTH FIDUCIAL  truth pT > 4.5, |truth eta| < 2.4  -- removes muons that enter the reco
//         fiducial region only through mismeasurement (§3.0(c')). Threshold moved 4 -> 4.5 with
//         the reco one (user decision D4, docs/tracking/mu_pt45_gap125_pairpt9_adoption.md §3(d)):
//         reco and truth thresholds must never sit at different values.
//   (3.5) FORWARD LOW-pT VETO  pT > 7 || q*eta > -2  -- the r16578 production's forward-negative
//         turn-on is saturated (R8/R10/R14, a real simulation problem, not bad muons), so those
//         muons are out of the dR-correlation measurement. ASYMMETRIC: only the corner that is
//         SIMULTANEOUSLY low-pT and forward-negative.
//   (gap) FIDUCIAL GAP CUT from ParamsSet::single_mu_fiducial_gap_cuts on BOTH legs, AND the
//         PAIR-LEVEL window |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2 (user,
//         2026-09-07) -- all built from ParamsSet so no number is retyped into a JIT string.
// A PAIR is kept only if BOTH legs satisfy each of them.
//
// ⚠ MIRROR NOTICE (state as of 2026-09-08). RDFBasedHistFilling/FillMCTrigEffHists.cxx is the
// ORIGINAL; this header was extracted from it verbatim. `kGapPair` there CALLS FiducialGapCut()
// below, so the gap construction cannot drift. STILL DUPLICATED inline in that file, and every
// copy must therefore carry the SAME numbers and the SAME pair-level WP column:
//   * kTruthFidSingle / kTruthFidLeg / kTruthFidPair (FillMCTrigEffHists.cxx) -- all three now at
//     truth pT > 4.5, |truth eta| < 2.4, mirroring TruthFiducial() below.
//   * kFwdVetoSingle / kFwdVetoLeg / kFwdVetoPair -- mirror ForwardLowPtVeto() below (pT > 7 is
//     NOT a muon threshold and does NOT move with the 4 -> 4.5 change).
//   * the RECO term: `pair_wp && <leg> pt > 4.5 && fabs(<leg> eta) < 2.4` for BOTH legs, in
//     sel_pair_legs (Steps 2/4), the inline Step-3 string, and the kn-split string --
//     mirroring Step3PairSelection() below. `sel_single` (Step 1) is the ONE deliberate
//     exception: it is a SINGLE-MUON map, there is no pair, so it keeps the per-muon
//     ANY-vertex flag `pass_tight`/`pass_medium` and only the pT threshold moves.
// This header provides no single-muon or lg_/ot_ leg-alias variant, which is why those copies
// still exist. A drift there is exactly the class of silent error this header exists to prevent.
// =================================================================================================
namespace MCTrigEffPairSel {

// Column-name convention of the Step-3 pair node: m1_pt, m1_eta, m1_charge, m1_truth_pt,
// m1_truth_eta (and m2_*), plus the PAIR-level `pair_wp`, all created by Alias on the pair tree.
// `pair_wp` must alias PairWpBranch(use_tight_wp) -- see the (wp) bullet in the header comment.

// The ntuple-processing pair-level WP flag for the requested working point. It is a plain leaf of
// the split MuonPairObj branch (like pair_pt / pair_eta), so RDF can read it by this bare name.
inline std::string PairWpBranch(bool use_tight_wp)
{
    return use_tight_wp ? "pair_pass_tight" : "pair_pass_medium";
}

// Fail LOUDLY, and BEFORE the event loop, if the pair-level WP column is missing on a tree.
// ROOT swallows exceptions thrown INSIDE an RDF event loop and the job STILL exits 0, leaving a
// fresh near-empty output file (memory: reference_root_swallows_rdf_exceptions) -- so a sample
// whose pair tree predates `pair_pass_*` must be caught at booking time, never by silently
// dropping the same-vertex requirement. `available` is RDataFrame/RNode::GetColumnNames(), which
// is evaluated when it is called, i.e. before any event is read.
inline void RequirePairWpColumn(const std::vector<std::string>& available,
                                const std::string& pair_wp_col,
                                const std::string& context)
{
    if (std::find(available.begin(), available.end(), pair_wp_col) != available.end()) return;
    const std::string msg =
        "MCTrigEffPairSel: the pair-level WP column \"" + pair_wp_col + "\" is MISSING on "
        + context + ". The MC trigger-efficiency pair selection REQUIRES it: it carries the pp "
        "SAME-VERTEX pair requirement, which the per-muon pass_tight/pass_medium flags "
        "deliberately do not. Re-run the ntuple processing for this sample; do NOT fall back to "
        "the per-muon flags -- that would measure the efficiency on a looser population than the "
        "one it is applied to.";
    std::cerr << "\n##### FATAL " << msg << " #####\n" << std::endl;
    throw std::runtime_error(msg);
}

inline std::string TruthFiducial()
{
    return "m1_truth_pt > 4.5 && fabs(m1_truth_eta) < 2.4 && "
           "m2_truth_pt > 4.5 && fabs(m2_truth_eta) < 2.4";
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
    // `pair_wp` (NOT m1_wp && m2_wp): the pair flag additionally carries the pp SAME-VERTEX
    // requirement -- see the (wp) bullet in the header comment. It implies both per-muon flags by
    // construction, so this is equivalent-or-tighter, never looser.
    std::string sel = "pair_wp && "
                      "m1_pt > 4.5 && fabs(m1_eta) < 2.4 && "
                      "m2_pt > 4.5 && fabs(m2_eta) < 2.4 && " + TruthFiducial()
                    + " && " + ForwardLowPtVeto();
    if (apply_gap_cut) sel += " && " + FiducialGapCut();
    return sel;
}

// The DATA-LIKE single-b RECO signal selection, byte-identical to `signal_cuts` in
// RDFBasedHistFilling/RDFBasedHistFillingPP.cxx with the m1./m2. dots replaced by the underscore
// aliases. As of 2026-09-08 this ALSO matches Pb+Pb: its signal region was migrated onto the
// same fiducial + pair-level windows and the same ParamsSet pair-pT threshold (D1), closing
// the pp-vs-PbPb divergence that stood from 2026-08-18. It previously read the retired one-sided
// `q*eta < 2.2` (docs/signal_selection_change_impact.md §0).
// Ground truth for the values: docs/analysis_overview.md §2
// and docs/signal_selection_change_impact.md §0. There is NO dR cut (removed 2026-06-22).
// 2026-08-17: the per-muon one-sided `q*eta < 2.2` was REPLACED by the detector-gap fiducial
// cut on BOTH muons (FiducialGapCut() above) -- kept in lockstep with the crossx by
// construction, since both read ParamsSet::single_mu_fiducial_gap_cuts. 2026-09-07: the
// PAIR-LEVEL |eta^pair| < 2.2 window joined it, inside the same FiducialGapCut() helper, so
// the lockstep holds for it too.
// 2026-09-08: the signal-region pair-pT cut moved 8 -> 9 GeV (user decision D3,
// docs/tracking/mu_pt45_gap125_pairpt9_adoption.md). It is no longer retyped here: it is read
// from ParamsSet::SignalPairPtCutExpr, the same single source `signal_cuts` in
// RDFBasedHistFilling/RDFBasedHistFillingPP.cxx reads, so the two cannot drift. That matters
// because CheckSignalWindowMirror (Utilities/PairTrigEffEvaluator.h) compares only the minv
// half of this string -- a pair-pT mismatch here would have been SILENT.
inline std::string SingleBSignalCutsReco()
{
    return ParamsSet::SignalMinvCutExpr("minv") + " && "
         + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
         + FiducialGapCut();
}

}  // namespace MCTrigEffPairSel

#endif  // MC_TRIG_EFF_PAIR_SELECTION_H
