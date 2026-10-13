#pragma once
// Shared config for PbPb event-level selection.
// Included by both the plotting code (cut derivation) and PbPbExtras (cut application)
// so that key names and cut order never drift between the two.

#include <string>
#include <unordered_set>
#include <vector>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <cstdlib>
#include <unistd.h>
#include "TGraph.h"
#include "TFile.h"
#include "TParameter.h"

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
    // Cut 2 = per-side ZDC time window (see PbPbZdcTimeCut below), stored as 4 edges [ns].
    constexpr const char* kZDCTimeALoNs  = "ZDC_time_A_lo_ns";
    constexpr const char* kZDCTimeAHiNs  = "ZDC_time_A_hi_ns";
    constexpr const char* kZDCTimeCLoNs  = "ZDC_time_C_lo_ns";
    constexpr const char* kZDCTimeCHiNs  = "ZDC_time_C_hi_ns";
    // LEGACY symmetric |t^A|,|t^C| < x, written by derivations before 2026-10-12 (still the only cut-2 key of the
    // 2023/2024 files).  Never written any more; read ONLY through PbPbReadZdcTimeWindow.
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
// 2023: two ZDC-problem runs (data_per_run_signal_yield_qa.md; historically mislabelled "b-hadron").
// 2026 (USER DECISION 2026-10-12; LB boundaries = the ATLAS twiki's final selection,
// slightly different from our own N_trk/FCal LB scan): runs 522384 + 522408 entirely;
// run 522336 LB 344-521; run 522355 LB >= 448.  Same removals as the GRL
// SkimCode/xmls/physics_HI2026_50ns_noIBL_updated.xml (the 2026 skim itself was made
// with the original physics_HI2026_50ns_noIBL.xml, so they are applied here).
// See Analysis/docs/tracking/pbpb2026_analysis_support.md (D7).
// 2026 D8 (USER DECISION 2026-10-13, TEMPORARY pending collaborators): runs 522200 and 522041
// also dropped entirely (shifted ZDC timing; data_per_run_signal_yield_qa.md). NOT in the GRL.
inline std::unordered_set<int> PbPbBadRuns(int run_year) {
    int yr = run_year % 2000;
    if (yr == 23) return {461674, 462964};
    if (yr == 26) return {522384, 522408, 522200, 522041};   // D7 + D8
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

// ---- Cut 2: ZDC time window, per year and per side (single source) ---------
// Pass iff a_lo < t^A < a_hi AND c_lo < t^C < c_hi  (t^A = zdc_ZdcTime[1], t^C = zdc_ZdcTime[0]; strict, like the
// old |t| < 1.5).  A fixed constant, not fitted.  Default +-1.5 ns on both sides.
// USER DECISION 2026-10-12 (PROVISIONAL, being finalized from the event-selection plots):
//   2025: t^A lower edge -0.6 ns (removes the t^A ~ -1.4 ns shoulder);
//   2026: t^A upper edge +0.8 ns (removes the t^A in [1,2] ns cluster).
// Both clusters sit at nominal vertex z (not main-satellite collisions).
// See Analysis/docs/tracking/pbpb_zdc_time_cut2_window.md.
// The derivation writes these edges into the cuts file; readers MUST use PbPbReadZdcTimeWindow, which throws if the
// file and this header disagree (= stale cuts file, re-derive).
struct PbPbZdcTimeWindow { double a_lo, a_hi, c_lo, c_hi; };

inline PbPbZdcTimeWindow PbPbZdcTimeCut(int run_year) {
    const int yr = run_year % 2000;
    PbPbZdcTimeWindow w{-1.5, 1.5, -1.5, 1.5};
    if (yr == 23 || yr == 24) return w;
    if (yr == 25) { w.a_lo = -0.6; return w; }
    if (yr == 26) { w.a_hi =  0.8; return w; }
    throw std::runtime_error("PbPbZdcTimeCut: no ZDC time window defined for Pb+Pb year " + std::to_string(yr));
}

inline bool PbPbPassZdcTime(const PbPbZdcTimeWindow& w, float tA, float tC) {
    return tA > (float)w.a_lo && tA < (float)w.a_hi && tC > (float)w.c_lo && tC < (float)w.c_hi;
}

// Write the 4 edges and remove the legacy symmetric key (files opened in UPDATE mode would otherwise keep it).
inline void PbPbWriteZdcTimeWindow(TFile* f, const PbPbZdcTimeWindow& w) {
    f->cd();
    const std::pair<const char*, double> kv[4] = {{PbPbEvSelKey::kZDCTimeALoNs, w.a_lo}, {PbPbEvSelKey::kZDCTimeAHiNs, w.a_hi},
                                                  {PbPbEvSelKey::kZDCTimeCLoNs, w.c_lo}, {PbPbEvSelKey::kZDCTimeCHiNs, w.c_hi}};
    for (auto& p : kv) { TParameter<double> t(p.first, p.second); t.Write(p.first, TObject::kOverwrite); }
    f->Delete((std::string(PbPbEvSelKey::kZDCTimeCutNs) + ";*").c_str());
}

// Read the window from a cuts file and check it against PbPbZdcTimeCut(run_year).  A file written before the per-side
// window (only the legacy symmetric key) is accepted only if the header window is that same symmetric box.
inline PbPbZdcTimeWindow PbPbReadZdcTimeWindow(TFile* f, int run_year) {
    const PbPbZdcTimeWindow h = PbPbZdcTimeCut(run_year);
    auto get = [&](const char* k) { return dynamic_cast<TParameter<double>*>(f->Get(k)); };
    auto* aL = get(PbPbEvSelKey::kZDCTimeALoNs); auto* aH = get(PbPbEvSelKey::kZDCTimeAHiNs);
    auto* cL = get(PbPbEvSelKey::kZDCTimeCLoNs); auto* cH = get(PbPbEvSelKey::kZDCTimeCHiNs);
    const std::string where = std::string(" in ") + f->GetName() + " (year " + std::to_string(run_year % 2000) + ")";
    if (aL && aH && cL && cH) {
        const PbPbZdcTimeWindow r{aL->GetVal(), aH->GetVal(), cL->GetVal(), cH->GetVal()};
        if (std::abs(r.a_lo - h.a_lo) > 1e-9 || std::abs(r.a_hi - h.a_hi) > 1e-9 ||
            std::abs(r.c_lo - h.c_lo) > 1e-9 || std::abs(r.c_hi - h.c_hi) > 1e-9)
            throw std::runtime_error("PbPbReadZdcTimeWindow: ZDC time window" + where +
                                     " differs from PbPbZdcTimeCut() -- stale cuts file, re-derive the event selection");
        return r;
    }
    if (auto* old = get(PbPbEvSelKey::kZDCTimeCutNs)) {
        const double t = old->GetVal();
        if (std::abs(h.a_lo + t) < 1e-9 && std::abs(h.a_hi - t) < 1e-9 && std::abs(h.c_lo + t) < 1e-9 && std::abs(h.c_hi - t) < 1e-9)
            return h;
        throw std::runtime_error("PbPbReadZdcTimeWindow: only the legacy symmetric ZDC time cut" + where +
                                 ", but PbPbZdcTimeCut() is a different window -- re-derive the event selection");
    }
    throw std::runtime_error("PbPbReadZdcTimeWindow: no ZDC time cut" + where);
}

// ---- Raw-NTUP input path: read LOCALGROUPDISK (/pnfs) files over xrootd ----
// Many data_pbpbYY_partN.root are symlinks into /pnfs.  NFS reads of some of them hang in uninterruptible D state
// (data_pbpb23_part2.root, 2026-10-12: a Stage-0 FCal comparison stuck 2 h 45 min), so the event-selection chain opens
// /pnfs targets through the dCache xrootd door instead (needs a valid grid proxy, X509_USER_PROXY).  The door is
// overridable with PNFS_XROOTD_DOOR (it moved 1094 -> 1096 once).  Non-/pnfs paths are returned unchanged.
inline std::string PbPbRawInputPath(const std::string& path) {
    char buf[4096];
    const ssize_t n = readlink(path.c_str(), buf, sizeof(buf) - 1);
    if (n <= 0) return path;
    const std::string tgt(buf, n);
    if (tgt.rfind("/pnfs/", 0) != 0) return path;
    const char* door = std::getenv("PNFS_XROOTD_DOOR");
    return std::string(door ? door : "root://dcgftp.usatlas.bnl.gov:1096/") + tgt;
}
