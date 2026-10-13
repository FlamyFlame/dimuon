#pragma once
// Shared config for PbPb event-level selection.
// Included by both the plotting code (cut derivation) and PbPbExtras (cut application)
// so that key names and cut order never drift between the two.

#include <string>
#include <unordered_set>
#include <vector>
#include <limits>
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
    // OPTIONAL: TTree with per-run mu+7sigma preamp cuts.  Present in the PbPb25 cuts
    // file only (as of 2026-09-10); PbPbExtras::InitEventSel loads it whenever it is
    // present, for any year, and falls back to the scalar cut when it is absent.
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

// ---- Run / LB data-quality exclusion (keyed by 2-digit or 4-digit year) ---
// Applied BEFORE ANY analysis step: the event-selection cut DERIVATION macros
// (plot_pbpb_event_sel_event_level / _cuts / plot_pbpb_fcal_comparison) and the
// NTuple processing (PbPbExtras::PassEventSel) drop these events first.
// NOTE: anything excluded here must ALSO have its luminosity removed from
// Utilities/PbPbSampledLumi.h and PbPbBaseClass.h::make_crossx_factors_pbpb_<yr>(),
// or the cross-section numerator and denominator cover different runs / LBs.
//
// 2023: two b-hadron runs.
// 2026 (USER DECISION 2026-10-12; LB boundaries = the ATLAS twiki's final selection,
// slightly different from our own N_trk/FCal LB scan): runs 522384 + 522408 entirely;
// run 522336 LB 344-521; run 522355 LB >= 448.  Same removals as the GRL
// SkimCode/xmls/physics_HI2026_50ns_noIBL_updated.xml (the 2026 skim itself was made
// with the original physics_HI2026_50ns_noIBL.xml, so they are applied here).
// See Analysis/docs/tracking/pbpb2026_analysis_support.md (D7).
inline std::unordered_set<int> PbPbBadRuns(int run_year) {
    int yr = run_year % 2000;
    if (yr == 23) return {461674, 462964};
    if (yr == 26) return {522384, 522408};
    return {};
}

// Partially-bad runs: inclusive LB ranges [lb_lo, lb_hi] to drop.
struct PbPbBadLBRange { int run; unsigned lb_lo; unsigned lb_hi; };
inline std::vector<PbPbBadLBRange> PbPbBadLBRanges(int run_year) {
    int yr = run_year % 2000;
    if (yr == 26) return {
        {522336, 344u, 521u},
        {522355, 448u, std::numeric_limits<unsigned>::max()},   // LB >= 448 (to end of run)
    };
    return {};
}

// Event-level helper for code that does not cache the lists (per-event cost is a hash
// lookup + a loop over <= a few ranges).
inline bool PbPbIsBadRunLB(const std::unordered_set<int>& bad_runs,
                           const std::vector<PbPbBadLBRange>& bad_lbs,
                           int run, unsigned lb) {
    if (bad_runs.count(run)) return true;
    for (const auto& r : bad_lbs)
        if (r.run == run && lb >= r.lb_lo && lb <= r.lb_hi) return true;
    return false;
}

// ---- TGraph cut evaluator with clamped extrapolation -----------------------
inline double PbPbEvSelEvalCut(const TGraph* g, double x) {
    if (!g || g->GetN() == 0) return -1.;
    if (x <= g->GetX()[0])           return g->GetY()[0];
    if (x >= g->GetX()[g->GetN()-1]) return g->GetY()[g->GetN()-1];
    return g->Eval(x);
}
