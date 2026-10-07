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
// INTERIM mass-window-only factor (Value / DenomExpr / NumExpr, nominal central-sample inputs); the
// decided TARGET (doc PD4/PD9, 2026-10-06) also puts the muon cuts in A_sig and is built from the private
// no-muon-pT-filter sample (TargetValue / Target*Expr below). The detector-gap cuts are NOT here -- they are in eps_reco
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

// ---- TARGET definition (doc PD4/PD9/D5-D7, 2026-10-07) -------------------------------------
// Needs the private no-muon-pT-filter sample processed with truth_muon_fiducial_target (NTuple keeps
// both truth muons |eta| <= 2.5, no pT cut). Columns sa_mu{1,2}_{pt,eta} are aliases of m{1,2}.truth_{pt,eta}
// defined by the writer.
//   T_target = den  = from_same_b && pair pT > 9 && |eta^pair| < 2.2 && both |eta_mu| <= 2.5
//   num_mupt  = den && both truth pT >= 4.5          (single-muon pT-cut factor)
//   num_mueta = den && both |truth eta| <= 2.4       (single-muon eta-cut factor)
//   num_mu    = den && both cuts                     (muon-cut factor)
//   num_mass  = den && mass window                   (mass-window factor on T_target)
//   num       = den && both cuts && mass window      (TOTAL A_sig)
inline const char* TargetValue() {
    return "TARGET A_sig = N[den && muon pT>=4.5 && |eta_mu|<=2.4 && truth mass window] / N[den]; "
           "den = T_target = from_same_b && truth pair pT > 9 && |truth eta^pair| < 2.2 && both "
           "|truth eta_mu|<=2.5, NO muon pT cut (private no-muon-pT-filter sample, NTuple "
           "truth_muon_fiducial_target); NO gap cuts; eta^pair axis 44 bins over [-2.2, 2.2] "
           "(2026-10-07, docs/tracking/signal_acceptance_mass_window.md PD4/D7)";
}
inline std::string TargetDenomExpr() {
    return DenomExpr() + " && " + ParamsSet::TargetMuEtaCutExpr("sa_mu1_eta", "sa_mu2_eta");
}
inline std::string TargetNumMuPtExpr()  { return TargetDenomExpr() + " && " + ParamsSet::SingleMuTruthPtCutExpr("sa_mu1_pt", "sa_mu2_pt"); }
inline std::string TargetNumMuEtaExpr() { return TargetDenomExpr() + " && " + ParamsSet::SingleMuTruthEtaCutExpr("sa_mu1_eta", "sa_mu2_eta"); }
inline std::string TargetNumMuExpr() {
    return TargetDenomExpr() + " && " + ParamsSet::SingleMuTruthPtCutExpr("sa_mu1_pt", "sa_mu2_pt")
           + " && " + ParamsSet::SingleMuTruthEtaCutExpr("sa_mu1_eta", "sa_mu2_eta");
}
inline std::string TargetNumMassExpr() { return TargetDenomExpr() + " && (" + ParamsSet::SignalMinvCutExpr("truth_minv") + ")"; }
inline std::string TargetNumExpr()     { return TargetNumMuExpr() + " && (" + ParamsSet::SignalMinvCutExpr("truth_minv") + ")"; }

}  // namespace SignalAcceptanceDefinition

#endif  // SIGNAL_ACCEPTANCE_DEFINITION_H
