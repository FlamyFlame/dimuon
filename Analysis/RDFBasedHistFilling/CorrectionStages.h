#pragma once

#include <array>

// =============================================================================
// CorrectionStages.h
//
// Sequential efficiency/unfolding correction stages applied to the
// cross-section weight, so that each correction's impact can be visualized at
// plotting time (raw -> unfolded -> +reco-eff -> +trig-eff).
//
// The cumulative weight columns named below are Define()'d identically in the
// PP and PbPb crossx fillers (RDFBasedHistFillingPP/PbPb::FillHistogramsCrossx),
// each starting from that sample's base crossx weight:
//   cw_raw                 = base crossx weight (lumi[/T_AA]-scaled, no eff)
//   cw_unfolded            = cw_raw * w_unfold        (w_unfold == 1: PLACEHOLDER
//                            identity until det-response unfolding exists,
//                            roadmap Q4)
//   cw_unfolded_reco       = cw_unfolded * w_reco     (w_reco = 1/(eps1*eps2),
//                            Run 2 single-muon reco-eff PLACEHOLDER, see
//                            docs/tracking/reco_eff_placeholder_run2.md; to be
//                            replaced by the proper 3D pair efficiency, task_05)
//   cw_unfolded_reco_trig  = cw_unfolded_reco * w_trig (per-pair trigger eff,
//                            already established)
//
// Invariant: with w_reco == 1 and w_unfold == 1, cw_unfolded_reco_trig equals
// the existing trigger-only crossx weight, so legacy outputs are reproduced.
// =============================================================================

// Correction building blocks (the three correction TYPES), in application order.
enum class CorrectionType { Unfolding, RecoEff, TrigEff };

// Cumulative correction stages.
enum class CorrectionStage {
    Raw = 0,
    Unfolded,
    UnfoldedReco,
    UnfoldedRecoTrig
};

struct CorrectionStageInfo {
    CorrectionStage stage;
    const char*     suffix;      // appended to the crossx histogram name
    const char*     weight_col;  // RDF column carrying the cumulative weight
};

// Stage table shared by the PP and PbPb crossx fillers.
inline const std::array<CorrectionStageInfo, 4>& CrossxCorrectionStages() {
    static const std::array<CorrectionStageInfo, 4> stages = {{
        {CorrectionStage::Raw,              "_corr_raw",               "cw_raw"},
        {CorrectionStage::Unfolded,         "_corr_unfolded",          "cw_unfolded"},
        {CorrectionStage::UnfoldedReco,     "_corr_unfolded_reco",     "cw_unfolded_reco"},
        {CorrectionStage::UnfoldedRecoTrig, "_corr_unfolded_reco_trig","cw_unfolded_reco_trig"},
    }};
    return stages;
}
