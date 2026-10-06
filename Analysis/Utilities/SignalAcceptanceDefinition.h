#ifndef SIGNAL_ACCEPTANCE_DEFINITION_H
#define SIGNAL_ACCEPTANCE_DEFINITION_H

#include <string>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"

// =================================================================================================
// SignalAcceptanceDefinition -- the ONE statement of the truth signal acceptance A_sig, shared by
// the writers (RDFBasedHistFillingPythiaTruth / RDFBasedHistFillingPowhegTruth,
// FillHistogramsSignalAcceptance), the cutflow diagnostic and the plotter.
// Physics: docs/tracking/signal_acceptance_mass_window.md (Physics Procedure §2).
//
// The measured range is truth pair pT > 9 GeV and |truth eta^pair| < 2.2: it DEFINES the
// measurement and is applied to numerator AND denominator (never acceptance-corrected). The muon
// truth fiducial pT > 4.5, |eta| < 2.4 is applied at the NTuple stage, so THIS implementation is the
// INTERIM mass-window-only factor; the decided target (doc PD4/PD9, 2026-10-06) also puts the muon
// cuts in A_sig and needs a private no-muon-pT-filter sample (not yet built). The detector-gap cuts are NOT here -- they are in eps_reco
// (docs/tracking/pair_reco_eff_gap_acceptance.md). So A_sig = the truth MASS-WINDOW acceptance.
//
// Before 2026-10-05 the same histogram names held num = mass && pair pT && gap windows &&
// |eta^pair| < 2.2 over den = from_same_b (no range cut). The writers stamp this marker; the
// plotter REFUSES a file without it. Change the VALUE whenever the definition changes.
// =================================================================================================
namespace SignalAcceptanceDefinition {

inline const char* Key() { return "signal_acceptance_definition"; }

inline const char* Value() {
    return "A_sig = N[den && truth mass window] / N[den]; den = from_same_b && truth pair pT > 9 && "
           "|truth eta^pair| < 2.2 (muon pT>4.5,|eta|<2.4 at NTuple stage); NO gap cuts; "
           "eta^pair axis 44 bins over [-2.2, 2.2] "
           "(2026-10-05, docs/tracking/signal_acceptance_mass_window.md)";
}

// Truth measured range (denominator). Expressions from ParamsSet, never retyped.
inline std::string DenomExpr() {
    return std::string("from_same_b && ") + ParamsSet::SignalPairPtCutExpr("truth_pair_pt")
           + " && " + ParamsSet::PairFiducialEtaCutExpr("truth_pair_eta");
}

// Numerator = denominator && truth mass window.
inline std::string NumExpr() {
    return DenomExpr() + " && (" + ParamsSet::SignalMinvCutExpr("truth_minv") + ")";
}

}  // namespace SignalAcceptanceDefinition

#endif  // SIGNAL_ACCEPTANCE_DEFINITION_H
