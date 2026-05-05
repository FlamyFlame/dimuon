#pragma once
// Shared config for PbPb event-level selection.
// Included by both the plotting code (cut derivation) and PbPbExtras (cut application)
// so that key names and cut order never drift between the two.

#include <string>
#include <unordered_set>
#include "TGraph.h"

// ---- Cut order (1-indexed; 0 unused) ----------------------------------------
enum class PbPbEvSelCut {
    kZDCFCalBanana = 1,
    kZDCTime       = 2,
    kZDCPreamp     = 3,
    kNTrkFrac      = 4,
    kNTrkFCalBand  = 5
};

// Human-readable labels indexed [0..5]; [0] unused.
// Single definition; included files must not define their own.
inline const char* kPbPbEvSelCutLabel(int idx) {
    static const char* kLabels[6] = {
        "",
        "ZDC_FCal_banana",
        "ZDC_time",
        "ZDC_preamp",
        "nTrk_frac",
        "nTrk_FCal_band"
    };
    return (idx >= 0 && idx <= 5) ? kLabels[idx] : "";
}

// ---- ROOT key names in event_sel_cuts_pbpb_20YY.root ------------------------
namespace PbPbEvSelKey {
    constexpr const char* kZDCFCalCut    = "g_ZDC_FCal_cut";
    constexpr const char* kZDCTimeCutNs  = "ZDC_time_cut_ns";
    constexpr const char* kPreampACutADC = "ZDC_preamp_A_cut_ADC";
    constexpr const char* kPreampCCutADC = "ZDC_preamp_C_cut_ADC";
    constexpr const char* kNTrkFracCutLo   = "g_ntrk_frac_cut_lo";
    constexpr const char* kNTrkFCalCutLo   = "g_ntrk_fcal_cut_lo";
    constexpr const char* kNTrkFCalCutHi   = "g_ntrk_fcal_cut_hi";
    // PbPb25 only: TTree with per-run mu+7sigma preamp cuts
    // Branches: run_number (Int_t), cut_A_ADC (Double_t), cut_C_ADC (Double_t)
    constexpr const char* kPreampPerRunTree = "t_preamp_per_run";
}

// ---- Path to per-year cuts ROOT file (nominal, fit-derived) -----------------
inline std::string PbPbEvSelCutsPath(int run_year) {
    int yr = run_year % 2000;
    return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" +
           std::to_string(yr) + "/event_sel_cuts_pbpb_20" +
           std::to_string(yr) + ".root";
}

// ---- Run-quality exclusion list (keyed by 2-digit or 4-digit year) ---------
// Returns runs excluded from all event selection derivations and analysis.
inline std::unordered_set<int> PbPbBadRuns(int run_year) {
    int yr = run_year % 2000;
    if (yr == 23) return {461674, 462964};
    return {};
}

// ---- TGraph cut evaluator with clamped extrapolation -----------------------
inline double PbPbEvSelEvalCut(const TGraph* g, double x) {
    if (!g || g->GetN() == 0) return -1.;
    if (x <= g->GetX()[0])           return g->GetY()[0];
    if (x >= g->GetX()[g->GetN()-1]) return g->GetY()[g->GetN()-1];
    return g->Eval(x);
}
