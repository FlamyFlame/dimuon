// =================================================================================================
// mc_pthat_slice_mass_statistics.cxx
//
// HOW MUCH TRIGGER-EFFICIENCY SAMPLE would a NEW MC request in the two highest pT-hat slices buy,
// and where must the loose dimuon MASS FILTER sit so that it removes the back-to-back pairs
// without biasing the measurement?
//
// This macro answers both, per pT-hat slice, on the pp24-conditions Pythia fullsim FULL sample:
//
//   (1) the dimuon MASS DISTRIBUTION of the selected pairs, per slice and per pair sign,
//       inclusively and restricted to the three highest coarse pair-pT cells -- the figure that
//       shows WHERE the back-to-back peak sits relative to the proposed m < 10 GeV filter;
//   (2) the RAW PAIR COUNT and the PER-EVENT RATE in m < 10 GeV and in the signal window
//       1.08 < m < 2.9 GeV, per slice and per sign (2 x 2 table, one CSV per slice);
//   (3) the same two mass ranges, resolved in the 3 highest coarse pair-pT cells x the 3
//       |eta^pair| groups (3 x 3 table, one CSV per slice x mass range).
//
// WHY THE COUNTS ARE RAW AND UNWEIGHTED. This is a SAMPLE-SIZE question, not a cross-section
// question. Everything here is measured INSIDE ONE pT-hat slice, where the per-pair MC weight is
// one constant (sigma * eps_filt / N_slice), so the weighted and unweighted numbers carry the same
// information and only the raw one answers "how many pairs will a request of N events give me".
// The per-event rate quoted beside each count is exactly the number that multiplies a requested
// N_events. It is normalized to N_slice = the entries of the NTUP chain the ntuple processing
// loops over; see the N_slice block below for what that check can and cannot detect.
//
// THE SAMPLE IS THE TRIGGER-EFFICIENCY PROCEDURES' OWN. The selection is
// MCTrigEffPairSel::Step3PairSelection(true) -- byte-identical to the population the DeltaR
// correction (mc_trigger_efficiency.md 3.3) and the single-value pair efficiency
// (mc_trigeff_single_value_pair_eff.md) are measured on, with NO trigger requirement and no
// signal-region pair-pT cut. That is what makes these counts a statement about the statistics
// those two procedures would gain, rather than about some other population. Read from the per-kn
// pair trees `muon_pair_tree_kin{N}_sign{M}` that the ntuple processing already writes
// (PythiaAlgCoreT::fill_kn_trees_fullsim) -- ntuple-processing OUTPUT, nothing re-derived from raw
// NTUPs (.claude/CLAUDE.md NTuple-Processing Provenance).
//
// BINNINGS ARE NEVER RETYPED. Pair pT = ParamsSet::pair_pt_coarse_bins via PairTrigEff::
// PairPtEdges(); |eta^pair| = the 3 sign-independent groups via PairTrigEff::AbsEtaGroups(); the
// signal mass window = PairTrigEff::Window("sig"), guarded against the signal region by
// PairTrigEff::CheckSignalWindowMirror(). The ONE number this file introduces is the proposed
// loose filter itself, kLooseMassMax (the subject of the request).
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'mc_pthat_slice_mass_statistics.cxx+("pp_full", true)'    // Tight  (nominal)
//   root -l -b -q 'mc_pthat_slice_mass_statistics.cxx+("pp_full", false)'   // Medium (WP syst.)
// =================================================================================================

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>

#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../Utilities/MCTrigEffPairSelection.h"
#include "../../../Utilities/PairTrigEffEvaluator.h"
#include "dr_correction_sample_cfg.h"

namespace {

// ------------------------------------------------------------------ the pT-hat slices of interest
// The bounds MIRROR PythiaAlgCoreT.c:26 (`kinRanges = {8,14,24,40,70,125,300}`), which is a private
// member of the algorithm template with no shared header to read it from. The mirror is NOT taken
// on trust: `lo`/`hi` build the AMI file NAME, and the DSID inside that file is checked against
// `dsid` below -- so a wrong bound either fails to open a file or fails the DSID guard, and can
// never silently mislabel a slice. `ikin` is the index of the per-kn tree written by
// PythiaAlgCoreT::InitOutputTreesExtra_PythiaCore (`muon_pair_tree_kin{ikin}_sign{1,2}`).
// DSIDs: the pp24 fullsim FULL `_pdf` production, docs/ami_weights.md table B.
struct PtHatSlice {
    int         ikin;    // per-kn tree index
    int         lo, hi;  // pT-hat bounds [GeV], as they appear in the file names
    int         dsid;    // the DSID this slice MUST have in the FULL production
    std::string token;   // file-name token
    std::string text;    // canvas/legend text
};
const std::vector<PtHatSlice>& Slices()
{
    static const std::vector<PtHatSlice> s = {
        {4,  70, 125, 803019, "pTH70_125",  "70 < #hat{p}_{T} < 125 GeV"},
        {5, 125, 300, 803015, "pTH125_300", "125 < #hat{p}_{T} < 300 GeV"},
    };
    return s;
}

// ------------------------------------------------------------------ the mass ranges of the tables
// `sig` is NOT typed here: it is the signal window, taken from PairTrigEff::Window("sig") and
// guarded against the signal region itself by CheckSignalWindowMirror(). `loose` IS the subject of
// the request -- the proposed generator-level filter, wide enough not to bias the measurement and
// tight enough to remove the back-to-back pairs, whose different physics would otherwise average
// into the same efficiency cells.
constexpr double kLooseMassMax = 10.0;   // GeV -- the PROPOSED loose dimuon mass filter

struct MassRange {
    std::string token, csv_text, tex_text;
    double      lo, hi;
};
std::vector<MassRange> MassRanges()
{
    const PairTrigEff::MassWindow& sig = PairTrigEff::Window("sig");
    char sig_csv[128], sig_tex[160];
    snprintf(sig_csv, sizeof sig_csv, "%g < m < %g GeV", sig.lo, sig.hi);
    snprintf(sig_tex, sizeof sig_tex, "%g < m_{#mu#mu} < %g GeV", sig.lo, sig.hi);
    char loose_csv[128], loose_tex[160];
    snprintf(loose_csv, sizeof loose_csv, "m < %g GeV", kLooseMassMax);
    snprintf(loose_tex, sizeof loose_tex, "m_{#mu#mu} < %g GeV", kLooseMassMax);
    return {
        {"m_lt_10",    loose_csv, loose_tex, 0.0,    kLooseMassMax},
        {"m_sig",      sig_csv,   sig_tex,   sig.lo, sig.hi},
    };
}

// ------------------------------------------------------------------ the mass COMPOSITION bands
// FOUR EXCLUSIVE, EXHAUSTIVE bands of the top-3 pair-pT population. This is the evidence for the
// question the request turns on -- "does the proposed filter separate the single-b population
// from the back-to-back one?" -- and it exists as exact `Count()`s because the earlier answer was
// read off the display histogram. NO log display axis has an edge at 2.9, 4 or 10 -- on the axis
// in use when that error was made they fell at 2.8602 / 4.1822 / 9.8358, and they move again every
// time the axis does -- so integrating one splits the bands at the wrong masses. It moved the
// headline same-sign number by 6.2 percentage points. D4 says every counted number comes from an
// exact filter; so does this.
//
// EVERY boundary is read from a canonical source, none is typed:
//   2.9  = PairTrigEff::Window("sig").hi   -- the signal window's top
//   4.0  = PairTrigEff::Window("wide").hi  -- the template-fit window's top
//   10   = kLooseMassMax                   -- the proposed filter, the subject of the request
struct MassBand { std::string token, csv_text; double lo, hi; };
std::vector<MassBand> MassBands()
{
    const double sig_hi  = PairTrigEff::Window("sig").hi;
    const double wide_hi = PairTrigEff::Window("wide").hi;
    if (!(sig_hi < wide_hi && wide_hi < kLooseMassMax))
        throw std::runtime_error("mc_pthat_slice_mass_statistics: the mass bands are not ordered "
                                 "(signal top, template-fit top, proposed filter) -- one of the "
                                 "canonical windows moved; re-derive the bands deliberately.");
    auto lbl = [](const char* f, double a, double b) { char s[96]; snprintf(s, sizeof s, f, a, b);
                                                       return std::string(s); };
    return {
        {"lt_sig",   lbl("m < %g GeV",        sig_hi, 0.),            0.,      sig_hi},
        {"sig_wide", lbl("%g <= m < %g GeV",  sig_hi, wide_hi),       sig_hi,  wide_hi},
        {"wide_cut", lbl("%g <= m < %g GeV",  wide_hi, kLooseMassMax), wide_hi, kLooseMassMax},
        {"above",    lbl("m >= %g GeV",       kLooseMassMax, 0.),     kLooseMassMax, 1e9},
    };
}

// A DIAGNOSTIC sub-count, not a band: how much of the [2.9, 4) band is the unvetoed J/psi. Without
// it "74 % of the kept opposite-sign sample is below 4 GeV" reads as continuum when much of the
// 2.9-4 part is one resonance (this pair file carries no resonance veto, PP-3b). The window is
// ParamsSet::minv_cuts[1] = {2.9, 3.3}, the analysis's own J/psi veto -- read, not retyped.
inline std::pair<double, double> JpsiWindow()
{
    static const ParamsSet pms;
    for (const auto& w : pms.minv_cuts)
        if (std::fabs(w[0] - PairTrigEff::Window("sig").hi) < 1e-6) return {w[0], w[1]};
    throw std::runtime_error("mc_pthat_slice_mass_statistics: no ParamsSet::minv_cuts window "
                             "starts at the signal window's upper edge -- the J/psi veto window "
                             "moved; do not guess it.");
}

// ------------------------------------------------------------------ the mass DISPLAY axis
// A display axis for the mass FIGURE only. Every counted number in this macro comes from an exact
// `minv > lo && minv < hi` filter, never from integrating this histogram, so no table depends on
// these edges.
//
// It is defined here rather than read from a canonical vector because NO existing minv binning
// spans the range this figure has to show. `hist_binning_map["minv_log_bins_{ss,op}"]`
// (RDFBasedHistFillingBaseClass.cxx) stops at 60 GeV -- below the back-to-back region we are
// looking for -- and is sign-DEPENDENT, so SS and OS could not be drawn on one axis;
// `ParamsSet::hq_minvBins` starts at 2.41 GeV, above the signal window the same figure must show.
// Uniform in log(m), so the bin content IS dN/dln(m) up to one constant and the position of a
// peak on the drawn curve is not a binning artefact.
// The low edge is below the 2*m_mu = 0.2113 GeV kinematic threshold, so the first populated bin is
// physics and not a cut. The high edge is above the HIGHEST selected pair mass measured on this
// sample (358.96 GeV, pTH125_300 opposite sign), so nothing is pushed into overflow and silently
// absent from the drawn curve -- 220 GeV, the first value used here, left 30 pairs off scale
// (17 OS + 9 SS in pTH125_300, 2 + 2 in pTH70_125).
constexpr double kMassAxisLo   = 0.211;  // GeV -- just below 2*m_mu = 0.21133, the kinematic
                                         // threshold, so bin 1 is not half-empty by construction
constexpr double kMassAxisHi   = 400.0;  // GeV
constexpr int    kMassAxisBins = 80;
std::vector<double> MassAxis()
{
    std::vector<double> e(kMassAxisBins + 1);
    const double step = std::log(kMassAxisHi / kMassAxisLo) / kMassAxisBins;
    for (int i = 0; i <= kMassAxisBins; ++i) e[i] = kMassAxisLo * std::exp(i * step);
    return e;
}

// ------------------------------------------------------------------ spectrum landmarks
// WHERE IS THE BACK-TO-BACK PEAK, and is there a gap below it? The answer is the whole reason the
// figure exists, and it is EMITTED here rather than read off the rendered plot by a human.
//
// WHY IT IS CODE AND NOT A NOTE IN THE DOC. Three review iterations in a row caught a version of
// the same failure: a peak/minimum/density number transcribed into the tracking doc from one
// display axis, then left behind when the axis changed (70 bins 0.2-220 -> 80 bins 0.2-400 ->
// 80 bins 0.211-400). Every bin edge moves each time and every one of those numbers goes stale
// silently, because nothing recomputes them. Emitting them beside the figure they describe is the
// only version of this that cannot drift.
//
// The search regions are DERIVED, not tuned:
//   peak    = the tallest bin ABOVE the proposed filter -- that is what "the back-to-back peak"
//             means operationally, and it needs no hand-picked window;
//   minimum = the shallowest bin between the template-fit window's top (above which the single-b
//             population has ended) and that peak;
//   density = the last bin ENTIRELY below the proposed filter. On a log axis no bin edge lands on
//             the filter, so the straddling bin is reported separately rather than silently used.
//
// THE COMPARISON WINDOW IS PRE-REGISTERED, AND DERIVED. `drop_sigma` compares the mean density
// just BELOW the filter with the mean density just ABOVE it, over a window that is symmetric in
// log(m) about the cut:
//     below = [filter / R, filter)      above = [filter, filter * R)      R = filter / wide_hi
// With the canonical values (filter = 10, template-fit top = 4) R = 2.5, so the windows are
// [4, 10) and [10, 25) GeV. R is DERIVED from the two canonical mass scales, not tuned, and
// neither edge depends on the data.
//
// WHY THIS MATTERS AND WHAT IT REPLACES. The first version ended the "above" window at the DATA's
// own tallest bin, while claiming the regions were "fixed before looking". That excluded the
// largest bin from the mean by construction, biasing the step POSITIVE, and it made `nbins_above`
// vary from row to row so the rows compared different mass ranges as if they were alike. Caught by
// /review-analysis-code (CRITICAL) and independently by the session running the pT-4.5 adoption.
// A window chosen by the data is not a pre-registered window, and a statistic computed on one may
// not be reported as though it were.
//
// The window dependence is REPORTED, not hidden: `drop_sigma_wide`, over [filter, axis top), is
// emitted beside it. If the two disagree in sign or size, the "step" is a property of the window
// and not of the spectrum, and neither may be quoted as a result.
//
// RESOLVABILITY, AND THE LOOK-ELSEWHERE TRAP. A "peak" and a "minimum" of 25 and 7 pairs on a flat
// ~15/bin continuum are Poisson noise. Comparing the tallest and shallowest bin of a ~20-bin scan
// does NOT test that: the extremes of a flat distribution are guaranteed to be far apart, so such a
// test says "resolved" almost always. The single-bin columns are therefore reported WITH sqrt(N)
// and explicitly marked look-elsewhere biased, and the verdict this file actually stands behind is
// `drop_sigma` -- the significance of the step in MEAN DENSITY across the proposed filter,
// (below - above) / sqrt(err_below^2 + err_above^2), which is a comparison of two fixed regions
// chosen before looking and is what the physics question ("is the filter at a minimum?") means.
struct Landmark {
    int    peak_bin,  min_bin,  below_bin,  straddle_bin;
    double peak_lo, peak_hi, peak_n;
    double min_lo,  min_hi,  min_n;
    double below_lo, below_hi, below_n;
    bool   extremes_differ;   // look-elsewhere BIASED; see the header note
    double drop_sigma;        // pre-registered log-symmetric window [filter/R, filter) vs [filter, filter*R)
    double drop_sigma_wide;   // same, but [filter, axis top) -- the window-dependence check
    double above_mean_w, above_err_w; int above_nbins_w;
    // REGION-level quantities. A single bin is a Poisson draw: the top-3 pTH125_300 "minimum"
    // is a 2.5 sigma dip on a flat plateau, and quoting it as a measurement is the same error as
    // quoting the pTH70_125 one. These are the statements that survive a rebinning AND carry
    // their own error.
    double below_mean, below_err;   // mean pairs/bin over [template-fit top, filter)
    int    below_nbins;
    double above_mean, above_err;   // mean pairs/bin over [filter, peak)
    int    above_nbins;
    double gmean_above;             // geometric-mean mass of the pairs above the filter [GeV]
};
Landmark FindLandmarks(const TH1D* h, double filter_m, double search_lo)
{
    Landmark L{};
    L.peak_bin = 0; L.peak_n = -1.;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if (h->GetBinLowEdge(i) >= filter_m && h->GetBinContent(i) > L.peak_n) {
            L.peak_n = h->GetBinContent(i); L.peak_bin = i;
        }
    L.min_bin = 0; L.min_n = 1e18;
    for (int i = 1; i < L.peak_bin; ++i)
        if (h->GetBinLowEdge(i) >= search_lo && h->GetBinContent(i) < L.min_n) {
            L.min_n = h->GetBinContent(i); L.min_bin = i;
        }
    if (L.min_bin == 0) { L.min_n = 0.; }
    // the last bin entirely below the filter, and the bin that straddles it
    L.below_bin = 0; L.straddle_bin = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        if (h->GetBinLowEdge(i + 1) <= filter_m) L.below_bin = i;
        if (h->GetBinLowEdge(i) < filter_m && h->GetBinLowEdge(i + 1) > filter_m) L.straddle_bin = i;
    }
    auto edges = [&](int b, double& lo, double& hi, double& n) {
        lo = b ? h->GetBinLowEdge(b) : 0.; hi = b ? h->GetBinLowEdge(b + 1) : 0.;
        n  = b ? h->GetBinContent(b) : 0.;
    };
    edges(L.peak_bin,  L.peak_lo,  L.peak_hi,  L.peak_n);
    edges(L.min_bin,   L.min_lo,   L.min_hi,   L.min_n);
    edges(L.below_bin, L.below_lo, L.below_hi, L.below_n);
    L.extremes_differ = L.peak_bin && L.min_bin &&
                        (L.peak_n - L.min_n) > (std::sqrt(std::max(L.peak_n, 1.))
                                              + std::sqrt(std::max(L.min_n, 1.)));

    // Region means, and the binning-robust location of the high-mass population.
    auto region = [&](double lo, double hi, double& mean, double& err, int& nb) {
        double sum = 0.; nb = 0;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            if (h->GetBinLowEdge(i) >= lo && h->GetBinLowEdge(i + 1) <= hi) {
                sum += h->GetBinContent(i); ++nb;
            }
        mean = nb ? sum / nb : 0.;
        err  = nb ? std::sqrt(sum) / nb : 0.;      // Poisson on the summed count
    };
    // Pre-registered, log-symmetric about the filter. `search_lo` IS filter/R by construction
    // (it is the template-fit window's top), so the two windows are mirror images in log(m).
    const double R = filter_m / search_lo;
    region(search_lo, filter_m, L.below_mean, L.below_err, L.below_nbins);
    region(filter_m, filter_m * R, L.above_mean, L.above_err, L.above_nbins);
    // The window-dependence check: everything above the filter, to the axis top.
    region(filter_m, h->GetBinLowEdge(h->GetNbinsX() + 1),
           L.above_mean_w, L.above_err_w, L.above_nbins_w);

    // The geometric mean mass ABOVE the filter locates the back-to-back population without
    // depending on which bin happens to be tallest -- the tallest-bin comparison between the two
    // slices is a coin flip when their maxima are within a Poisson error of each other.
    double sw = 0., swl = 0.;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if (h->GetBinLowEdge(i) >= filter_m && h->GetBinContent(i) > 0.) {
            sw  += h->GetBinContent(i);
            swl += h->GetBinContent(i) * std::log(h->GetBinCenter(i));
        }
    L.gmean_above = sw > 0. ? std::exp(swl / sw) : 0.;

    // The step in mean density across the filter. POSITIVE = the density falls as the filter is
    // crossed, i.e. the filter sits ABOVE the single-b population's edge rather than in a gap.
    const double de = std::sqrt(L.below_err * L.below_err + L.above_err * L.above_err);
    L.drop_sigma = de > 0. ? (L.below_mean - L.above_mean) / de : 0.;
    const double dw = std::sqrt(L.below_err * L.below_err + L.above_err_w * L.above_err_w);
    L.drop_sigma_wide = dw > 0. ? (L.below_mean - L.above_mean_w) / dw : 0.;
    return L;
}

// ------------------------------------------------------------------ presentation order of the signs
// PairTrigEff::Signs() is the canonical list and is used unchanged for every tree name and token;
// this is ONLY a display order, because the request asks for "same sign & opposite sign" columns
// in that order. Reordering a copy, rather than typing a second list, keeps the tokens and the
// tree names single-sourced.
std::vector<PairTrigEff::PairSign> SignsSsFirst()
{
    std::vector<PairTrigEff::PairSign> v;
    for (const auto& t : {"ss", "os"}) v.push_back(PairTrigEff::Sign(t));
    return v;
}

// ------------------------------------------------------------------ small helpers
std::string Fmt(const char* f, double a)            { char b[256]; snprintf(b, sizeof b, f, a);     return b; }
std::string Fmt(const char* f, double a, double b_) { char b[256]; snprintf(b, sizeof b, f, a, b_); return b; }

// A canvas headline is ROOT latex; a CSV comment is plain text. Deliberately targeted, so an
// unexpected latex construct stays visible rather than being silently mangled.
std::string PlainText(std::string s)
{
    const std::vector<std::pair<std::string, std::string>> subs = {
        {"#sqrt{s_{NN}}", "sqrt(s_NN)"}, {"#sqrt{s}", "sqrt(s)"}, {"#hat{p}_{T}", "pT-hat"}};
    for (const auto& kv : subs)
        // Resume PAST the inserted text, not at it: a replacement that contained its own key
        // would otherwise loop forever. None does today; the cost of not relying on that is zero.
        for (size_t p = s.find(kv.first); p != std::string::npos;
             p = s.find(kv.first, p + kv.second.size()))
            s.replace(p, kv.first.size(), kv.second);
    return s;
}

// The repo convention: ONE writer, called with std::cout and with the file, so the log and the
// file cannot say different things.
void Emit(const std::string& path, const std::function<void(std::ostream&)>& write)
{
    std::cout << "\n----- " << path << " -----\n";
    write(std::cout);
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("mc_pthat_slice_mass_statistics: cannot open '" + path
                                       + "' for writing");
    write(ofs);
}

// The cell as the request asks for it: "number (percentage)".
std::string Cell(double n, double n_events)
{
    return Fmt("%.0f", n) + " (" + Fmt("%.4f", n_events > 0. ? 100. * n / n_events : 0.) + "%)";
}

// ------------------------------------------------------------------ AMI: sigma * eps_filt + DSID
// Same parse as PythiaAlgCoreT.c InitInputFullsim, and the same hard DSID guard: the AMI files are
// named by beam + pT-hat slice ONLY, so a file from ANOTHER production opens perfectly and every
// number that follows is wrong with no warning (docs/ami_weights.md, THE RULE). Here the AMI
// numbers are used only to RECONSTRUCT N_slice from the per-pair weight, but a wrong file would
// make that reconstruction wrong -- i.e. it would corrupt the denominator of every percentage in
// this macro -- so the guard is not optional.
struct AmiInfo { double sigma_nb, gen_filt_eff, sigma_eff_nb, total_events; int dsid; std::string path; };

double AmiField(const std::string& path, const std::string& key)
{
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind(key, 0) != 0) continue;                 // key must START the line
        const size_t c = line.find(':');
        if (c == std::string::npos) continue;
        double v = 0.;
        std::istringstream(line.substr(c + 1)) >> v;
        return v;
    }
    throw std::runtime_error("mc_pthat_slice_mass_statistics: field '" + key
                             + "' not found in AMI file " + path);
}

AmiInfo ReadAmi(const std::string& sample_dir, const PtHatSlice& sl)
{
    AmiInfo a;
    a.path = sample_dir + "ami_info/ami_info_mc23_5p36TeV_Py8EG_A14_pp_hQCD_DiMu_pTH"
           + std::to_string(sl.lo) + "_" + std::to_string(sl.hi) + ".txt";
    std::ifstream probe(a.path);
    if (!probe.good())
        throw std::runtime_error("mc_pthat_slice_mass_statistics: missing AMI file " + a.path);
    probe.close();
    a.sigma_nb     = AmiField(a.path, "crossSection ");
    a.gen_filt_eff = AmiField(a.path, "genFiltEff ");
    a.dsid         = static_cast<int>(AmiField(a.path, "datasetNumber"));
    // The PRODUCTION's event count, independent of anything on local disk -- the handle on an
    // incomplete NTUP farm (see NtupChainEntries).
    a.total_events = AmiField(a.path, "totalEvents ");
    a.sigma_eff_nb = a.sigma_nb * a.gen_filt_eff;
    if (a.dsid != sl.dsid)
        throw std::runtime_error(
            "mc_pthat_slice_mass_statistics: AMI PROVENANCE MISMATCH for " + sl.token + ": "
            + a.path + " has datasetNumber=" + std::to_string(a.dsid) + ", expected "
            + std::to_string(sl.dsid) + " for the pp24 fullsim FULL production. The AMI files are "
            "named by beam + pT-hat slice ONLY, so another production's file opens silently -- "
            "see docs/ami_weights.md.");
    if (!(a.sigma_eff_nb > 0.))
        throw std::runtime_error("mc_pthat_slice_mass_statistics: non-positive sigma*eps_filt in "
                                 + a.path);
    return a;
}

// ------------------------------------------------------------------ N_slice: the events actually processed
// N_slice is the denominator of EVERY percentage in this macro.
//
// IT IS READ FROM THE FILE, not inferred. `meta_tree_out` in the pair file carries
// `nproc_kin<K>_beam0` -- the number of events the ntuple processing ACTUALLY LOOPED OVER for this
// pT-hat slice -- alongside `nbeam_kin<K>_beam0` (the entries available). That is the
// authoritative number and this macro requires it.
//
// NOTE, verified against the file rather than assumed: upstream ALSO has an in-memory
// `meta_fullsim_truncated` flag (PythiaAlgCoreT.h:141) but **never Branches it**, so it does NOT
// reach the pair file -- `meta_tree_out` carries only the `nentries_`, `nproc_` and `nbeam_`
// families. The only file-level truncation detector is therefore `nproc != nbeam`, which is what
// this macro uses.
//
// HISTORY, because it is the reason the requirement is written this way. Until 2026-09-09 the
// fullsim path filled no meta tree, the per-pair weight was built from N_beam while the loop ran
// over N_proc = min(N_beam, nevents_max), and NOTHING in the pair file recorded N_proc. A run
// truncated by `nevents_max` -- which the pipeline's smoke test does, onto this same `_full` path
// -- therefore produced a file in which every available cross-check agreed at N_beam while every
// per-event rate was understated by N_proc/N_beam. `/review-analysis-code` found that this macro's
// then "denominator proved twice" was not two independent routes at all (both reduced to N_beam),
// it was raised with the session that owns PythiaAlgCoreT, and it was fixed upstream (D12): the
// weight now divides by N_proc and the meta tree records both counts. A pair file predating that
// fix has no meta tree and is REFUSED here rather than silently mis-normalised.
//
// TWO CROSS-CHECKS on the number read from the file, both required to pass:
//   (a) `nbeam_kin<K>_beam0` must equal the entries of the NTUP chain on disk -- catches a pair
//       file and an NTUP farm that are out of step, or a partially symlinked farm;
//   (b) sigma*eps_filt / w must equal N_proc, where w is the constant per-pair weight of the slice
//       and sigma*eps_filt comes from that slice's own DSID-guarded AMI file -- catches a weight
//       built from a different count than the loop used.
//   plus AMI `totalEvents`, the production's own record, which lives off this machine entirely.

Long64_t NtupChainEntries(const std::string& sample_dir, const PtHatSlice& sl)
{
    const std::string base = sample_dir + "Pythia_5p36TeV_pp_hQCD_DiMu_pTH"
                           + std::to_string(sl.lo) + "_" + std::to_string(sl.hi) + "."
                           + FullSimSampleFileTag(FullSimSampleType::pp) + ".NTUP";
    std::ifstream single(base + ".root");
    const bool have_single = single.good();
    single.close();

    TChain ch("HeavyIonD3PD");
    const int nadd = have_single ? ch.Add((base + ".root").c_str())
                                 : ch.Add((base + ".part*.root").c_str());
    if (nadd <= 0)
        throw std::runtime_error("mc_pthat_slice_mass_statistics: no NTUP input for " + sl.token
                                 + " -- looked for '" + base + ".root' and '" + base
                                 + ".part*.root'");
    // The same mutual exclusion PythiaAlgCoreT enforces: a stale hadded file shadowing the farm
    // would give a different N with an unchanged sigma.
    if (have_single) {
        TChain probe("HeavyIonD3PD");
        if (probe.Add((base + ".part*.root").c_str()) > 0)
            throw std::runtime_error("mc_pthat_slice_mass_statistics: AMBIGUOUS NTUP input for "
                                     + sl.token + ": BOTH the hadded file and the multi-part farm "
                                     "exist. Remove one.");
    }
    return ch.GetEntries();
}

// The per-slice event bookkeeping the ntuple processing now writes. A pair file that predates it
// is refused: without N_proc there is no way to tell a full run from a truncated one, and every
// per-event rate below would be silently understated.
struct MetaCounts { Long64_t nproc, nbeam; };
MetaCounts ReadMetaCounts(const std::string& pair_file, int ikin)
{
    std::unique_ptr<TFile> f(TFile::Open(pair_file.c_str(), "READ"));
    if (!f || f->IsZombie())
        throw std::runtime_error("mc_pthat_slice_mass_statistics: cannot open " + pair_file);
    TTree* t = dynamic_cast<TTree*>(f->Get("meta_tree_out"));
    const std::string stale =
        "mc_pthat_slice_mass_statistics: " + pair_file + " has no usable `meta_tree_out`, so the "
        "number of events ACTUALLY PROCESSED for each pT-hat slice is unrecorded. This pair file "
        "predates the fullsim event-bookkeeping fix (PythiaAlgCoreT, 2026-09-09). Refusing to run: "
        "a run truncated by `nevents_max` -- which the pipeline's smoke test does, onto this same "
        "`_full` path -- would understate every per-event rate here with no other symptom. Re-run "
        "the ntuple processing for this sample.";
    if (!t || t->GetEntries() < 1) throw std::runtime_error(stale);

    const std::string bp = "nproc_kin" + std::to_string(ikin) + "_beam0";
    const std::string bb = "nbeam_kin" + std::to_string(ikin) + "_beam0";
    if (!t->GetBranch(bp.c_str()) || !t->GetBranch(bb.c_str())) throw std::runtime_error(stale);

    MetaCounts m{0, 0};
    t->SetBranchAddress(bp.c_str(), &m.nproc);
    t->SetBranchAddress(bb.c_str(), &m.nbeam);
    t->GetEntry(0);
    if (m.nproc <= 0)
        throw std::runtime_error("mc_pthat_slice_mass_statistics: " + bp + " is "
                                 + std::to_string(m.nproc) + " in " + pair_file);
    return m;
}

// The single constant weight of one pT-hat slice. Every pair of the slice must carry it; a spread
// would mean the slice is not one production and N could not be reconstructed from it at all.
double ConstantSliceWeight(ROOT::RDF::RNode d, const std::string& what)
{
    auto mn = d.Min("weight");
    auto mx = d.Max("weight");
    auto n  = d.Count();
    if (*n == 0) throw std::runtime_error("mc_pthat_slice_mass_statistics: " + what + " is EMPTY");
    if (!(*mn > 0.))
        throw std::runtime_error("mc_pthat_slice_mass_statistics: non-positive pair weight on "
                                 + what);
    // Relative tolerance on the value itself -- an absolute one would never fire for a weight of
    // order 1e-6 (memory: reference_relative_tolerance_guards).
    if (*mx - *mn > 1e-6 * *mx)
        throw std::runtime_error("mc_pthat_slice_mass_statistics: the pair weight is NOT constant "
                                 "on " + what + " (min " + Fmt("%.10g", *mn) + ", max "
                                 + Fmt("%.10g", *mx) + "). One pT-hat slice must carry one "
                                 "sigma*eps_filt/N_slice, so N_slice cannot be reconstructed.");
    return *mn;
}

}  // namespace

// =================================================================================================
void mc_pthat_slice_mass_statistics(const std::string& sample = "pp_full",
                                    bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    ROOT::EnableImplicitMT();

    // The `sig` window must still BE the signal region's mass window, or the table's second row
    // would describe a different selection from the one the analysis uses.
    PairTrigEff::CheckSignalWindowMirror();

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    if (cfg.key != "pp_full")
        throw std::invalid_argument(
            "mc_pthat_slice_mass_statistics: 'pp_full' only. The pT-hat slice indices, the DSID "
            "guard and the pp-beam-only isospin weight of 1 (which is what makes "
            "N = sigma*eps_filt/weight correct) are all properties of the pp24 fullsim FULL "
            "production; got '" + cfg.key + "'.");

    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_col  = use_tight_wp ? "pass_tight" : "pass_medium";
    const std::string wp_text = use_tight_wp ? "Tight" : "Medium";

    // Same pair file the dR correction and the single-value pair efficiency read.
    const std::string pair_file = cfg.sample_dir
        + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root";
    Long_t f_id = 0, f_sz = 0, f_fl = 0, f_mt = 0;
    if (gSystem->GetPathInfo(pair_file.c_str(), &f_id, &f_sz, &f_fl, &f_mt) != 0)
        throw std::runtime_error("mc_pthat_slice_mass_statistics: cannot stat " + pair_file);
    const std::string pair_file_stamp =
        pair_file + Form(" (mtime %ld, %ld B)", f_mt, f_sz);

    // ---------------------------------------------------------------- the cells (never retyped)
    const std::vector<double> pt_edges  = PairTrigEff::PairPtEdges();
    const DrAxisGroups        eta_grp   = PairTrigEff::AbsEtaGroups();
    const std::vector<double> eta_edges = eta_grp.edges;
    const int npt  = static_cast<int>(pt_edges.size())  - 1;
    const int neta = static_cast<int>(eta_edges.size()) - 1;
    // The three highest coarse pair-pT cells = the request's cells. Expressed as the EDGE the
    // delivered range starts at, looked up in the live axis, so it follows a binning change
    // instead of being a second place to remember (PairTrigEff::FirstDeliveredPtEdge).
    const int first_cell = PairTrigEff::FirstDeliveredPtBin();
    if (npt - first_cell + 1 != 3)
        throw std::runtime_error("mc_pthat_slice_mass_statistics: the delivered pair-pT range is "
                                 + std::to_string(npt - first_cell + 1) + " cells, not the 3 the "
                                 "request's 3 x 3 tables are shaped for. The coarse pair-pT axis "
                                 "changed -- reshape the tables deliberately.");
    const double top_pt_lo = pt_edges[first_cell - 1];

    const std::vector<MassRange> ranges = MassRanges();
    const std::vector<double>    m_axis = MassAxis();

    const std::string pair_wp_col = MCTrigEffPairSel::PairWpBranch(use_tight_wp);
    const std::string base_sel    = MCTrigEffPairSel::Step3PairSelection(true);

    std::cout << "\n============ mc_pthat_slice_mass_statistics: " << sample << " / " << wp_text
              << " muons ============\n"
              << "  pair file   : " << pair_file_stamp << "\n"
              << "  selection   : " << base_sel << "\n"
              << "  pair pT     : " << npt << " coarse cells, top 3 = [" << top_pt_lo << ", "
              << pt_edges[npt] << ") GeV\n"
              << "  |eta^pair|  : " << DrGroupsDescribe(eta_grp, "|eta^pair|", "") << std::endl;

    // ---------------------------------------------------------------- book, per slice x sign
    struct Booked {
        ROOT::RDF::RResultPtr<TH1D> mass_all, mass_top;                  // the figure
        std::map<std::string, ROOT::RDF::RResultPtr<ULong64_t>> n_tot;   // per mass range, no cells
        std::map<std::string, ROOT::RDF::RResultPtr<ULong64_t>> n_band;  // exclusive bands, top-3
        std::map<std::string, ROOT::RDF::RResultPtr<TH2D>>      n_cell;  // per mass range, cells
        ROOT::RDF::RResultPtr<ULong64_t> n_sel, n_in_grid, n_top3;
    };
    std::map<std::string, Booked> B;                       // key = slice token + "_" + sign token
    std::vector<std::unique_ptr<ROOT::RDataFrame>> rdf_store;
    std::map<std::string, Long64_t> n_events;              // per slice
    std::map<std::string, AmiInfo>  ami;                   // per slice

    for (const auto& sl : Slices()) {
        // ---- the denominator: READ from the file, then cross-checked (see the N_slice block) ----
        ami[sl.token] = ReadAmi(cfg.sample_dir, sl);
        const Long64_t   n_chain = NtupChainEntries(cfg.sample_dir, sl);
        const MetaCounts meta    = ReadMetaCounts(pair_file, sl.ikin);
        // AUTHORITATIVE: the events the ntuple processing actually looped over for this slice.
        n_events[sl.token] = meta.nproc;
        if (meta.nproc != meta.nbeam)
            std::cout << "  !! " << sl.token << ": this pair file is TRUNCATED -- "
                      << meta.nproc << " of " << meta.nbeam << " events were processed"
                      << ". The per-event rates below are still correct (they divide by N_proc), "
                         "but the ABSOLUTE counts are those of a partial run." << std::endl;
        // (a) the entries available on disk must match what the processing recorded as available.
        if (meta.nbeam != n_chain)
            throw std::runtime_error(
                "mc_pthat_slice_mass_statistics: for " + sl.token + " the pair file records "
                + std::to_string(meta.nbeam) + " available events but the NTUP farm on disk holds "
                + std::to_string(n_chain) + ". The pair file and the farm are out of step -- one "
                "was regenerated without the other, or a farm symlink is missing.");

        for (const auto& S : PairTrigEff::Signs()) {
            // The sign index comes from PairTrigEff's own tree name, not from a token
            // ternary: a ternary maps ANY unexpected token to sign2 silently.
            const char sign_digit = S.tree.back();
            if (sign_digit != '1' && sign_digit != '2')
                throw std::runtime_error("mc_pthat_slice_mass_statistics: cannot read a sign index "
                                         "from PairTrigEff tree name '" + S.tree + "'");
            const std::string ktree = "muon_pair_tree_kin" + std::to_string(sl.ikin)
                                    + "_sign" + std::string(1, sign_digit);
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(ktree, pair_file));
            // BEFORE the event loop: ROOT swallows exceptions thrown inside one and still exits 0
            // (memory: reference_root_swallows_rdf_exceptions).
            MCTrigEffPairSel::RequirePairWpColumn(rdf_store.back()->GetColumnNames(), pair_wp_col,
                                                  pair_file + ":" + ktree);
            ROOT::RDF::RNode d = *rdf_store.back();
            d = d.Alias("pair_wp", pair_wp_col)
                 .Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
                 .Alias("m1_wp", "m1." + wp_col)
                 .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                 .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
                 .Alias("m2_wp", "m2." + wp_col)
                 .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta");

            // N_slice from the weight, on the UNSELECTED node: the weight is an event-level
            // constant, so reading it before any pair cut is both cheaper and immune to a
            // selection that happened to keep zero pairs.
            // Checked on BOTH sign trees, not just one: half of every reported number comes from
            // the same-sign tree, and "one constant weight per slice" is a claim about the slice,
            // not about one of its two trees.
            {
                // (b) the weight must have been built from the SAME count the loop used. Since
                // the 2026-09-09 upstream fix the weight divides by N_proc, so this is exact.
                const double w = ConstantSliceWeight(d, sl.token + " " + ktree);
                const double n_from_w = ami.at(sl.token).sigma_eff_nb / w;
                const Long64_t n_w = static_cast<Long64_t>(std::llround(n_from_w));
                if (std::llabs(n_w - meta.nproc) > 1)
                    throw std::runtime_error(
                        "mc_pthat_slice_mass_statistics: N_slice DISAGREES for " + sl.token
                        + ": the per-pair weight implies " + Fmt("%.2f", n_from_w)
                        + " events (sigma*eps_filt = " + Fmt("%.6g", ami.at(sl.token).sigma_eff_nb)
                        + " nb / w = " + Fmt("%.10g", w) + "), but the pair file records "
                        + std::to_string(meta.nproc) + " events processed. The weight was built "
                        "from a different event count than the loop used -- every per-event rate "
                        "would be wrong by their ratio.");
                // (a) vs (c): an incomplete farm. The NTUP farm is 3 part-files of ~107 k
                // events each, so 2 % (6 400 of 320 000) sits far below a dropped part (33 %)
                // and far above the single event genuinely missing from each slice (0.0003 %).
                const double ami_ev = ami.at(sl.token).total_events;
                if (!(ami_ev > 0.) || n_chain > static_cast<Long64_t>(ami_ev)
                    || ami_ev - static_cast<double>(n_chain) > 0.02 * ami_ev)
                    throw std::runtime_error(
                        "mc_pthat_slice_mass_statistics: the NTUP farm for " + sl.token
                        + " holds " + std::to_string(n_chain) + " events but AMI reports the "
                        "production has " + Fmt("%.0f", ami_ev) + ". The farm is incomplete (a "
                        "part-file missing or a broken symlink), which would understate every "
                        "per-event rate. Restore the farm; do not scale the result.");
                if (S.token == "os")
                    std::cout << "  " << sl.token << ": DSID " << ami.at(sl.token).dsid
                              << ", sigma*eps_filt = " << Fmt("%.6g", ami.at(sl.token).sigma_eff_nb)
                              << " nb, w = " << Fmt("%.6g", w) << " nb, N_slice = N_proc = "
                              << meta.nproc << " of " << meta.nbeam << " available"
                              << " (weight implies " << Fmt("%.2f", n_from_w)
                              << ", NTUP chain " << n_chain
                              << ", AMI totalEvents " << Fmt("%.0f", ami_ev) << ")" << std::endl;
            }

            d = d.Filter(base_sel, ktree + " Step-3 MC trigger-efficiency pair selection")
                 .Define("abs_pair_eta", "fabs(pair_eta)");

            Booked b;
            b.n_sel = d.Count();
            // Inside the cell grid: the |eta^pair| axis is the sign-independent FOLD, so the cut
            // is on |pair_eta| -- applying folded bounds to the SIGNED branch would silently drop
            // every negative-eta pair (mc_trigeff_dr_binning_approaches.md D11).
            auto d_grid = d.Filter(Fmt("pair_pt >= %.10g && pair_pt < %.10g", pt_edges.front(),
                                       pt_edges.back())
                                   + " && " + Fmt("abs_pair_eta < %.10g", eta_edges.back()),
                                   "inside the cell grid");
            b.n_in_grid = d_grid.Count();

            const std::string key = sl.token + "_" + S.token;
            b.mass_all = d.Histo1D({("h_minv_" + key).c_str(),
                                    ";m_{#mu#mu} [GeV];muon pairs / bin",
                                    kMassAxisBins, m_axis.data()}, "minv");
            // BOTH bounds. Without the upper one this histogram is "p_T^pair above the third
            // -highest cell's lower edge", which is NOT the 3 highest cells: 155 opposite-sign
            // and 14 same-sign pTH125_300 pairs sit above the 150 GeV axis top, so the figure
            // would describe a slightly different population from the 3 x 3 tables it is read
            // beside. Both edges come from the canonical axis.
            auto d_top3 = d.Filter(Fmt("pair_pt >= %.10g && pair_pt < %.10g",
                                       top_pt_lo, pt_edges.back()),
                                   "the 3 highest coarse pair-pT cells");
            b.mass_top = d_top3.Histo1D({("h_minv_top3_" + key).c_str(),
                                         ";m_{#mu#mu} [GeV];muon pairs / bin",
                                         kMassAxisBins, m_axis.data()}, "minv");

            // MASS COMPOSITION of the top-3 cells, from EXACT cuts.
            // This is the evidence for "does the proposed filter separate the two populations",
            // and it must NOT be read off the display histogram: none of 2.9, 4 or 10 GeV is a
            // bin edge of that axis, so integrating it splits the bands at the wrong mass -- by
            // up to 6 percentage points, measured. D4 says every counted number comes from an
            // exact filter; these are those counts.
            for (const auto& B_ : MassBands())
                b.n_band.emplace(B_.token,
                                 d_top3.Filter(Fmt("minv >= %.10g && minv < %.10g", B_.lo, B_.hi),
                                               "top-3 mass band " + B_.token).Count());
            b.n_band.emplace("jpsi",
                             d_top3.Filter(Fmt("minv >= %.10g && minv < %.10g",
                                               JpsiWindow().first, JpsiWindow().second),
                                           "top-3 J/psi window").Count());
            b.n_top3 = d_top3.Count();

            for (const auto& R : ranges) {
                const std::string mass = Fmt("minv > %.10g && minv < %.10g", R.lo, R.hi);
                b.n_tot.emplace(R.token, d.Filter(mass, "mass " + R.token).Count());
                b.n_cell.emplace(R.token,
                    d_grid.Filter(mass, "mass " + R.token + " (cells)")
                          .Histo2D({("h_cells_" + key + "_" + R.token).c_str(),
                                    ";p_{T}^{pair} [GeV];|#eta^{pair}|;muon pairs",
                                    npt, pt_edges.data(), neta, eta_edges.data()},
                                   "pair_pt", "abs_pair_eta"));
            }
            B.emplace(key, std::move(b));
        }
    }

    // ---------------------------------------------------------------- run, then prove the input held still
    // This sample is regenerated by the ntuple processing, which other sessions run on the same
    // checkout. If the pair file were replaced WHILE this process reads it, RDF would report no
    // error at all: every histogram still fills, just with fewer pairs -- and "fewer pairs" is
    // indistinguishable from the honest answer to the very question this macro asks. So the event
    // loops are forced here, and the file is re-stat'ed before anything is written; a file that
    // moved means the numbers describe a half-written sample and nothing may be published from
    // them (memory: feedback_shared_checkout_input_drift).
    for (auto& kv : B) (void)*kv.second.n_sel;      // forces every booked result of every slice
    Long_t a_id = 0, a_sz = 0, a_fl = 0, a_mt = 0;
    const bool still_there = gSystem->GetPathInfo(pair_file.c_str(), &a_id, &a_sz, &a_fl, &a_mt) == 0;
    if (!still_there || a_sz != f_sz || a_mt != f_mt)
        throw std::runtime_error(
            "mc_pthat_slice_mass_statistics: the input pair file CHANGED while it was being read.\n"
            "  before: " + pair_file_stamp + "\n"
            "  after : " + (still_there ? pair_file + Form(" (mtime %ld, %ld B)", a_mt, a_sz)
                                        : std::string("GONE")) + "\n"
            "  The ntuple processing regenerated this sample mid-run, so the counts just produced "
            "describe a partially written file. Refusing to write them -- re-run once the "
            "regeneration has finished.");

    // ---------------------------------------------------------------- output directory
    // cfg.out_base is ".../pp_trigger_efficiency/mc_based<variant>/". This product set is a
    // SIBLING of the mc_based plot tree, not part of it, so it gets its own top-level directory
    // under the same plot root -- and, per the VARIANT LAYOUT of dr_correction_sample_cfg.h, the
    // working point lives in that top-level name.
    const std::string out_dir = DrCorrSiblingTree(cfg, std::string("mc_pthat_slice_mass_stats")
                                                       + (use_tight_wp ? "" : "_medium"));
    gSystem->mkdir(out_dir.c_str(), kTRUE);

    // ---------------------------------------------------------------- provenance header (all CSVs)
    auto header = [&](std::ostream& os, const PtHatSlice& sl) {
        os << "# pp24-conditions Pythia fullsim FULL sample: muon-pair statistics of ONE pT-hat "
              "slice,\n";
        os << "# for sizing a NEW MC request with a loose dimuon mass filter.\n";
        os << "# slice = " << sl.token << " (" << PlainText(sl.text) << "), DSID "
           << ami.at(sl.token).dsid << ", " << PlainText(cfg.sample_text) << "\n";
        os << "# muon working point = " << wp_text << " (required of BOTH legs)\n";
        os << "# source = " << pair_file_stamp << "\n";
        os << "#   tree muon_pair_tree_kin" << sl.ikin << "_sign{1=same,2=opposite}, written by the\n";
        os << "#   ntuple processing (PythiaAlgCoreT::fill_kn_trees_fullsim). Nothing is re-derived\n";
        os << "#   from raw NTUPs.\n";
        os << "#\n";
        os << "# DENOMINATOR of every percentage: N_slice = " << n_events.at(sl.token)
           << " events -- the events the ntuple\n";
        os << "#   processing ACTUALLY LOOPED OVER for this slice, read from `meta_tree_out`\n";
        os << "#   (nproc_kin" << sl.ikin << "_beam0) in the pair file, not inferred. Cross-checked\n";
        os << "#   three ways, all required to pass:\n";
        os << "#   (a) the events the processing recorded as AVAILABLE (nbeam_kin" << sl.ikin
           << "_beam0) equal the\n";
        os << "#       entries of the NTUP chain on disk;\n";
        os << "#   (b) sigma*eps_filt / w = N_proc, with sigma*eps_filt = "
           << Fmt("%.6g", ami.at(sl.token).sigma_eff_nb)
           << " nb from this slice's own\n";
        os << "#       DSID-guarded AMI file and w the constant per-pair weight -- i.e. the weight\n";
        os << "#       was built from the same count the loop used;\n";
        os << "#   (c) the AMI production record, totalEvents = "
           << Fmt("%.0f", ami.at(sl.token).total_events)
           << " (a partially downloaded or\n";
        os << "#       partially symlinked farm would show up here and nowhere else).\n";
        os << "#   A percentage here is therefore the PER-EVENT PAIR RATE: multiply it by the\n";
        os << "#   number of events you request to get the expected pair count.\n";
        os << "#\n";
        os << "#   ! WHAT `N` MEANS IN THAT PROJECTION. The rate is measured on a sample with NO\n";
        os << "#   generator-level dimuon mass filter, so `rate x N` is the yield from N events\n";
        os << "#   generated BEFORE any such filter. A request whose event count is quoted AFTER a\n";
        os << "#   mass filter (the usual convention, and what AMI totalEvents records) would give\n";
        os << "#   MORE pairs than this, by 1/eps_filter. The filter's own efficiency is a\n";
        os << "#   generator-level number and is NOT measured here.\n";
        os << "# COUNTS ARE RAW AND UNWEIGHTED. Inside ONE pT-hat slice the MC weight is a single\n";
        os << "#   constant, so weighting adds no information and only the raw count answers\n";
        os << "#   \"how many pairs will N events give me\". These are NOT yields: slices are\n";
        os << "#   mixed with very different weights.\n";
        os << "#\n";
        os << "# SELECTION = MCTrigEffPairSel::Step3PairSelection(true) -- byte-identical to the\n";
        os << "#   population the DeltaR correction (docs/tracking/mc_trigger_efficiency.md 3.3)\n";
        os << "#   and the single-value pair efficiency\n";
        os << "#   (docs/tracking/mc_trigeff_single_value_pair_eff.md) are measured on:\n";
        os << "#     " << base_sel << "\n";
        os << "#   i.e. of BOTH muons: the pair-level " << wp_text << " flag (both legs' combined\n";
        os << "#   quality + WP bit + IDCuts + MuonCuts, one-sided dp/p < 0.12, d0/z0, AND the pp\n";
        os << "#   SAME-VERTEX pair requirement), pT > 4.5 GeV, |eta| < 2.4, the truth fiducial,\n";
        os << "#   the forward low-pT veto, the per-leg fiducial gap windows and the pair-level\n";
        os << "#   |eta^pair| < " << ParamsSet::pair_eta_fiducial_max << ".\n";
        os << "#   NO TRIGGER REQUIREMENT and NO signal-region pair-pT cut.\n";
        os << "#\n";
    };

    // ---------------------------------------------------------------- SET 1: 2 mass x 2 sign
    for (const auto& sl : Slices()) {
        Emit(out_dir + "mass_window_counts_" + sl.token + ".csv", [&](std::ostream& os) {
            header(os, sl);
            os << "# TABLE: rows = dimuon mass range, columns = pair charge combination.\n";
            os << "# Cell = raw pair count (percentage of N_slice = " << n_events.at(sl.token)
               << " events).\n";
            os << "# NO pair-pT or pair-eta restriction beyond the selection above, so these are\n";
            os << "# ALL selected pairs of the mass range.\n";
            os << "# The two rows OVERLAP: " << ranges[1].csv_text << " is a subset of "
               << ranges[0].csv_text << ".\n";
            os << "mass_range,same_sign,opposite_sign\n";
            for (const auto& R : ranges) {
                os << R.csv_text;
                for (const auto& S : SignsSsFirst()) {
                    const double n = *B.at(sl.token + "_" + S.token).n_tot.at(R.token);
                    os << "," << Cell(n, n_events.at(sl.token));
                }
                os << "\n";
            }
            os << "#\n";
            os << "# The same numbers as plain columns, for machine reading:\n";
            os << "# mass_range,same_sign_pairs,same_sign_percent,opposite_sign_pairs,"
                  "opposite_sign_percent\n";
            for (const auto& R : ranges) {
                os << "# " << R.csv_text;
                for (const auto& S : SignsSsFirst()) {
                    const double n = *B.at(sl.token + "_" + S.token).n_tot.at(R.token);
                    os << "," << Fmt("%.0f", n)
                       << "," << Fmt("%.6f", 100. * n / n_events.at(sl.token));
                }
                os << "\n";
            }
            os << "#\n# Context: all selected pairs of this slice, any mass -- ";
            for (const auto& S : SignsSsFirst())
                os << S.text << " " << *B.at(sl.token + "_" + S.token).n_sel << "; ";
            os << "\n";
        });
    }

    // ---------------------------------------------------------------- the bands really are exhaustive
    // "Exclusive and exhaustive" is asserted in MassBands' comment and relied on by every
    // percentage in the composition table; here it is CHECKED. A pair with a NaN or absurd minv
    // would fall in no band, the four percentages would quietly fail to sum to 100, and the
    // "all pairs in these cells" row would exceed the band sum with no diagnostic.
    for (const auto& sl : Slices())
        for (const auto& S : SignsSsFirst()) {
            auto& b = B.at(sl.token + "_" + S.token);
            double sum = 0.;
            for (const auto& Bd : MassBands()) sum += (double)*b.n_band.at(Bd.token);
            const double tot = (double)*b.n_top3;
            if (std::fabs(sum - tot) > 0.5)
                throw std::runtime_error(
                    "mc_pthat_slice_mass_statistics: the mass bands are NOT exhaustive for "
                    + sl.token + " " + S.text + ": they hold " + Fmt("%.0f", sum)
                    + " pairs but the 3 highest pair-pT cells hold " + Fmt("%.0f", tot)
                    + ". Some pair has a mass outside every band (NaN, negative, or >= the last "
                      "band's upper bound) -- every percentage of the composition table would be "
                      "computed against a total the bands do not cover.");
        }

    // ---------------------------------------------------------------- SET 1b: mass composition of the top-3 cells
    // The table the filter DECISION is argued from: how the top-3 pair-pT population divides
    // between the single-b region, the template-fit region, the 4-10 GeV continuum the proposed
    // filter would KEEP, and the back-to-back region it would REMOVE. Exact cuts, never the
    // display histogram (see MassBands).
    for (const auto& sl : Slices()) {
        Emit(out_dir + "mass_composition_top3_" + sl.token + ".csv", [&](std::ostream& os) {
            header(os, sl);
            os << "# TABLE: how the pairs of the 3 HIGHEST coarse pair-pT cells ("
               << Fmt("%.4f", top_pt_lo) << " <= pair pT < " << Fmt("%.4f", pt_edges.back())
               << " GeV,\n";
            os << "#        all 3 |eta^pair| groups) divide between four EXCLUSIVE, EXHAUSTIVE\n";
            os << "#        mass bands. Rows = band, columns = pair charge combination.\n";
            os << "#        Cell = raw pair count (percentage of the pairs in these cells).\n";
            os << "# Band edges are read from canonical sources, never typed: the signal window's\n";
            os << "#   top and the template-fit window's top (PairTrigEff::Windows()) and the\n";
            os << "#   proposed filter. The LAST band is what the filter would remove.\n";
            os << "# NOTE the percentages here are of the top-3 pair population, NOT of N_slice --\n";
            os << "#   this table is about COMPOSITION, not about yield per event. For yield see\n";
            os << "#   mass_window_counts_*.csv and cells_*.csv.\n";
            os << "mass_band,same_sign,same_sign_pct,opposite_sign,opposite_sign_pct\n";
            for (const auto& Bd : MassBands()) {
                os << Bd.csv_text;
                for (const auto& S : SignsSsFirst()) {
                    auto& b = B.at(sl.token + "_" + S.token);
                    const double n = *b.n_band.at(Bd.token), tot = *b.n_top3;
                    os << "," << Fmt("%.0f", n) << "," << Fmt("%.2f", tot > 0. ? 100. * n / tot : 0.);
                }
                os << "\n";
            }
            os << "# all pairs in these cells";
            for (const auto& S : SignsSsFirst())
                os << "," << Fmt("%.0f", (double)*B.at(sl.token + "_" + S.token).n_top3) << ",100.00";
            os << "\n#\n";
            os << "# What the proposed filter KEEPS, as a composition (the first three bands):\n";
            os << "# sign,kept_pairs,pct_below_" << Fmt("%.4g", MassBands()[1].lo)
               << ",pct_" << Fmt("%.4g", MassBands()[2].lo) << "_to_"
               << Fmt("%.4g", kLooseMassMax) << "\n";
            for (const auto& S : SignsSsFirst()) {
                auto& b = B.at(sl.token + "_" + S.token);
                const double n0 = *b.n_band.at("lt_sig"), n1 = *b.n_band.at("sig_wide"),
                             n2 = *b.n_band.at("wide_cut");
                const double kept = n0 + n1 + n2;
                os << "# " << S.text << "," << Fmt("%.0f", kept)
                   << "," << Fmt("%.2f", kept > 0. ? 100. * n0 / kept : 0.)
                   << "," << Fmt("%.2f", kept > 0. ? 100. * n2 / kept : 0.) << "\n";
            }
            os << "#\n";
            os << "# DIAGNOSTIC: the " << MassBands()[1].csv_text << " band is NOT continuum -- this\n";
            os << "#   file carries no resonance veto, so it contains the J/psi. Pairs inside the\n";
            os << "#   analysis's own J/psi veto window (ParamsSet::minv_cuts, ["
               << Fmt("%.4g", JpsiWindow().first) << ", " << Fmt("%.4g", JpsiWindow().second)
               << ") GeV):\n";
            for (const auto& S : SignsSsFirst()) {
                auto& b = B.at(sl.token + "_" + S.token);
                const double nj = *b.n_band.at("jpsi"), nb = *b.n_band.at("sig_wide");
                os << "#   " << S.text << ": " << Fmt("%.0f", nj) << " of " << Fmt("%.0f", nb)
                   << " (" << Fmt("%.1f", nb > 0. ? 100. * nj / nb : 0.) << "% of the band)\n";
            }
        });
    }

    // ---------------------------------------------------------------- SET 1c: spectrum landmarks
    // Where the back-to-back peak is, whether there is a gap below it, and whether either is
    // statistically resolvable -- emitted, never transcribed (see FindLandmarks).
    Emit(out_dir + "mass_spectrum_landmarks.csv", [&](std::ostream& os) {
        header(os, Slices().front());
        os << "# (The provenance block above names the first slice; this file covers BOTH, one\n";
        os << "#  row per slice x sign x scope.)\n";
        os << "#\n";
        os << "# TABLE: landmarks of the dimuon mass spectrum, read from the histograms this run\n";
        os << "#   produced (mc_pthat_slice_mass_stats.root), so they cannot go stale against the\n";
        os << "#   display axis they are measured on. That axis is " << kMassAxisBins
           << " bins uniform in log(m),\n";
        os << "#   " << Fmt("%.4g", kMassAxisLo) << " - " << Fmt("%.4g", kMassAxisHi)
           << " GeV. Contents are RAW pair counts per bin; because the bins are\n";
        os << "#   uniform in log(m) the content IS dN/dln(m) up to one constant, so comparing\n";
        os << "#   bin heights compares densities and a peak position is not a binning artefact.\n";
        os << "#\n";
        os << "# peak    = the tallest bin ABOVE the proposed filter (m > "
           << Fmt("%.4g", kLooseMassMax) << " GeV) -- the back-to-back peak.\n";
        os << "# min     = the shallowest bin between the template-fit window's top ("
           << Fmt("%.4g", PairTrigEff::Window("wide").hi) << " GeV, above\n";
        os << "#           which the single-b population has ended) and that peak.\n";
        os << "# below   = the last bin ENTIRELY below the proposed filter. No bin edge lands on\n";
        os << "#           the filter, so the straddling bin is reported separately -- using it\n";
        os << "#           would compare a partly-above-filter density with a below-filter one.\n";
        os << "# extremes_differ = do peak and min differ by more than sqrt(N_peak)+sqrt(N_min)?\n";
        os << "#           ** LOOK-ELSEWHERE BIASED -- do not quote it. ** It compares the extremes\n";
        os << "#           of a ~20-bin scan, and the extremes of a FLAT distribution are far apart\n";
        os << "#           by construction, so it reads YES even on pure noise.\n";
        os << "# drop_sigma = THE verdict this file stands behind: the significance of the step in\n";
        os << "#           MEAN DENSITY across the filter, (below - above)/sqrt(err^2 + err^2), over a\n";
        os << "#           PRE-REGISTERED window symmetric in log(m) about the cut -- [filter/R,\n";
        os << "#           filter) vs [filter, filter*R) with R = filter / template-fit-top = "
           << Fmt("%.3g", kLooseMassMax / PairTrigEff::Window("wide").hi) << ", i.e.\n";
        os << "#           [" << Fmt("%.4g", PairTrigEff::Window("wide").hi) << ", "
           << Fmt("%.4g", kLooseMassMax) << ") vs [" << Fmt("%.4g", kLooseMassMax) << ", "
           << Fmt("%.4g", kLooseMassMax * kLooseMassMax / PairTrigEff::Window("wide").hi)
           << ") GeV. Neither edge depends on the data.\n";
        os << "#           POSITIVE = the density FALLS across the filter, i.e. the filter sits above\n";
        os << "#           the single-b population's edge -- NOT in a gap between the two populations.\n";
        os << "# drop_sigma_wide = the same step measured against EVERYTHING above the filter (to the\n";
        os << "#           axis top). It is the WINDOW-DEPENDENCE CHECK.\n";
        os << "# window_dependence_ok = do drop_sigma and drop_sigma_wide agree in SIGN?\n";
        os << "#           If NO, the step is a property of the window and not of the spectrum and\n";
        os << "#           **NEITHER value may be quoted as a result for that row**. The magnitudes\n";
        os << "#           are NOT expected to agree: the wide window reaches into the steeply\n";
        os << "#           falling high-mass tail, so it always gives the larger number. Only the\n";
        os << "#           SIGN is the robustness claim.\n";
        os << "#           EXPECT THIS TO FAIL ON THE `all selected pairs` ROWS, and it is not a\n";
        os << "#           defect: inclusively the back-to-back peak sits INSIDE the narrow [10, 25)\n";
        os << "#           window, so that window is dominated by the peak itself and measures\n";
        os << "#           something different from the wide one. In the 3 highest pair-pT cells the\n";
        os << "#           peak is near 45-50 GeV, OUTSIDE the narrow window, so there the narrow\n";
        os << "#           window measures the shoulder just above the filter -- which IS the\n";
        os << "#           quantity the filter question asks about.\n";
        os << "#\n";
        os << "# THE PHYSICS QUESTION THIS ANSWERS: if `min` lies ABOVE the proposed filter, then\n";
        os << "#   the filter cuts through a continuum rather than through a gap between the\n";
        os << "#   single-b and back-to-back populations.\n";
        os << "# THE ROBUST COLUMNS ARE THE REGION MEANS, NOT THE SINGLE BINS. A single bin is a\n";
        os << "#   Poisson draw; quote `mean_per_bin_below_filter` vs `mean_per_bin_above_filter`\n";
        os << "#   vs `peak_N`, which carry errors and survive a change of binning. The geometric\n";
        os << "#   mean mass above the filter locates the back-to-back population without\n";
        os << "#   depending on which bin happens to be tallest.\n";
        os << "slice,sign,scope,peak_range_GeV,peak_N,peak_sqrtN,min_range_GeV,min_N,min_sqrtN,"
              "extremes_differ_LOOK_ELSEWHERE_BIASED,min_above_filter,last_bin_below_filter_GeV,its_N,pct_of_peak,"
              "straddling_bin_GeV,mean_per_bin_below_filter,its_err,nbins_below,"
              "mean_per_bin_above_filter,its_err,nbins_above,geom_mean_mass_above_filter_GeV,drop_sigma,"
              "mean_per_bin_above_wide,its_err,nbins_above_wide,drop_sigma_wide,window_dependence_ok\n";
        for (const auto& sl : Slices())
            for (const auto& S : SignsSsFirst())
                for (int top = 0; top < 2; ++top) {
                    auto& b = B.at(sl.token + "_" + S.token);
                    TH1D* h = (top ? b.mass_top : b.mass_all).GetPtr();
                    const Landmark L = FindLandmarks(h, kLooseMassMax,
                                                     PairTrigEff::Window("wide").hi);
                    os << sl.token << "," << S.text << ","
                       << (top ? "3 highest pair-pT cells" : "all selected pairs") << ","
                       << Fmt("%.2f-%.2f", L.peak_lo, L.peak_hi) << "," << Fmt("%.0f", L.peak_n)
                       << "," << Fmt("%.1f", std::sqrt(std::max(L.peak_n, 0.))) << ","
                       << Fmt("%.2f-%.2f", L.min_lo, L.min_hi) << "," << Fmt("%.0f", L.min_n)
                       << "," << Fmt("%.1f", std::sqrt(std::max(L.min_n, 0.))) << ","
                       << (L.extremes_differ ? "YES" : "NO") << ","
                       << (L.min_lo >= kLooseMassMax ? "YES" : "no") << ","
                       << Fmt("%.2f-%.2f", L.below_lo, L.below_hi) << "," << Fmt("%.0f", L.below_n)
                       << "," << Fmt("%.1f", L.peak_n > 0. ? 100. * L.below_n / L.peak_n : 0.)
                       << "," << Fmt("%.2f-%.2f", h->GetBinLowEdge(L.straddle_bin),
                                     h->GetBinLowEdge(L.straddle_bin + 1))
                       << "," << Fmt("%.1f", L.below_mean) << "," << Fmt("%.1f", L.below_err)
                       << "," << L.below_nbins
                       << "," << Fmt("%.1f", L.above_mean) << "," << Fmt("%.1f", L.above_err)
                       << "," << L.above_nbins
                       << "," << Fmt("%.1f", L.gmean_above)
                       << "," << Fmt("%.2f", L.drop_sigma)
                       << "," << Fmt("%.1f", L.above_mean_w) << "," << Fmt("%.1f", L.above_err_w)
                       << "," << L.above_nbins_w
                       << "," << Fmt("%.2f", L.drop_sigma_wide)
                       << "," << (((L.drop_sigma >= 0.) == (L.drop_sigma_wide >= 0.))
                                  ? "YES" : "NO -- do not quote this row's step")
                       << "\n";
                }
    });

    // ---------------------------------------------------------------- SETS 2-5: 3 pT x 3 |eta|
    for (const auto& sl : Slices()) {
        for (const auto& R : ranges) {
            Emit(out_dir + "cells_" + sl.token + "_" + R.token + ".csv", [&](std::ostream& os) {
                header(os, sl);
                os << "# TABLE: rows = the 3 HIGHEST coarse pair-pT cells\n";
                os << "#        (ParamsSet::pair_pt_coarse_bins, cells " << first_cell << "-"
                   << npt << " of " << npt << "),\n";
                os << "#        columns = the 3 |eta^pair| groups\n";
                os << "#        (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, folded).\n";
                os << "# MASS RANGE: " << R.csv_text << ".\n";
                os << "# Cell = raw pair count (percentage of N_slice = " << n_events.at(sl.token)
                   << " events),\n";
                os << "#        SUMMED OVER BOTH SIGNS is NOT done: one block per sign follows.\n";
                os << "# Bin edges, exactly as used:\n";
                os << "#   pair pT [GeV] :";
                for (double e : pt_edges) os << " " << Fmt("%.4f", e);
                os << "\n#   |eta^pair|    :";
                for (double e : eta_edges) os << " " << Fmt("%.4f", e);
                os << "\n#\n";
                for (const auto& S : SignsSsFirst()) {
                    const TH2D* h = B.at(sl.token + "_" + S.token).n_cell.at(R.token).GetPtr();
                    os << "# ---- " << S.text << " pairs ----\n";
                    os << "pair_pT_GeV\\|pair_eta|";
                    for (int iy = 1; iy <= neta; ++iy)
                        os << "," << Fmt("%.1f-%.1f", eta_edges[iy - 1], eta_edges[iy]);
                    os << "\n";
                    for (int ix = first_cell; ix <= npt; ++ix) {
                        os << Fmt("%.1f-%.1f", pt_edges[ix - 1], pt_edges[ix]);
                        for (int iy = 1; iy <= neta; ++iy)
                            os << "," << Cell(h->GetBinContent(ix, iy), n_events.at(sl.token));
                        os << "\n";
                    }
                    os << "# these 3 x 3 cells total,"
                       << Cell(h->Integral(first_cell, npt, 1, neta), n_events.at(sl.token))
                       << "\n";
                    os << "# all " << npt << " x " << neta << " cells of the grid,"
                       << Cell(h->Integral(1, npt, 1, neta), n_events.at(sl.token)) << "\n";
                }
                os << "#\n# Pairs OUTSIDE the grid are counted nowhere in this table. Of all "
                      "selected pairs\n";
                os << "# (any mass), inside the grid: ";
                for (const auto& S : SignsSsFirst()) {
                    auto& b = B.at(sl.token + "_" + S.token);
                    const double in = *b.n_in_grid, all = *b.n_sel;
                    os << S.text << " " << Fmt("%.0f", in) << "/" << Fmt("%.0f", all) << " ("
                       << Fmt("%.2f", all > 0. ? 100. * in / all : 0.) << "%); ";
                }
                os << "\n";
            });
        }
    }

    // ---------------------------------------------------------------- the mass figure
    gStyle->SetOptStat(0);
    auto draw_set = [&](const std::string& png, bool top_cells_only, const std::string& scope) {
        // N = 2 panels -> nrow = 1 (subplot-layout convention).
        TCanvas c(("c_" + png).c_str(), "", 1500, 640);
        c.Divide(2, 1);
        std::vector<std::unique_ptr<TLegend>> legs;
        std::vector<std::unique_ptr<TLine>>   lines;
        std::vector<std::unique_ptr<TLatex>>  texs;

        auto hist_of = [&](const PtHatSlice& sl, const PairTrigEff::PairSign& S) {
            auto& b = B.at(sl.token + "_" + S.token);
            return (top_cells_only ? b.mass_top : b.mass_all).GetPtr();
        };

        // ONE y range for BOTH panels of a figure, computed before anything is drawn. The two
        // panels are the same quantity in two pT-hat slices and the figure's second message is
        // how much MORE the harder slice delivers -- with per-panel ranges that comparison is
        // invisible to the eye and has to be read out of the tables instead.
        double ymax = 0.;
        for (const auto& sl : Slices())
            for (const auto& S : SignsSsFirst())
                ymax = std::max(ymax, hist_of(sl, S)->GetMaximum());
        if (!(ymax > 0.))
            throw std::runtime_error("mc_pthat_slice_mass_statistics: every mass distribution ("
                                     + scope + ") is EMPTY");
        // Headroom for the legend, sized PER FIGURE. The legend sits in the upper right, so what
        // it must clear is the tallest curve under its own x-range (m > ~9 GeV), not the global
        // maximum. On the inclusive figure the 20-25 GeV back-to-back peak reaches ~23 % of ymax
        // there and needs the full x60; on the top-3 figure the same region reaches only ~10 %,
        // so x60 would leave a third of the frame blank for nothing.
        const double ylo = 0.5, yhi = ymax * (top_cells_only ? 15. : 60.);

        for (size_t i = 0; i < Slices().size(); ++i) {
            const PtHatSlice& sl = Slices()[i];
            c.cd(static_cast<int>(i) + 1);
            // LOG x: the axis is uniform in log(m) and the spectrum spans three decades in mass.
            // LOG y: permitted here because the x axis is log-binned and the distribution falls by
            // several decades (memory: feedback_log_scale_plots) -- on a linear y the low-mass
            // continuum would compress the back-to-back region this figure exists to show to a
            // flat line.
            gPad->SetLogx(); gPad->SetLogy();
            // The top margin holds BOTH header lines, so neither can ever sit on the data --
            // at low mass the spectrum reaches within a decade of the frame top.
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.13);

            // DRAW OPTION, and why it differs between the two figures. The top-pair-pT figure is
            // the one whose entire purpose is a STATISTICS argument, and its same-sign series
            // holds 1-7 entries per bin: drawn as a continuous HIST line with no uncertainty it
            // would invite the reader to believe structure that is Poisson noise, and its empty
            // bins would plunge to the axis floor and read as physical holes. So it is drawn with
            // markers and Poisson bars. The inclusive figure has 10^2-10^3 entries per bin, where
            // the bars are smaller than the line width and markers only add clutter.
            const char* opt      = top_cells_only ? "E1"      : "HIST";
            const char* opt_same = top_cells_only ? "E1 SAME" : "HIST SAME";

            bool first = true;
            for (const auto& S : SignsSsFirst()) {
                TH1D* h = hist_of(sl, S);
                h->SetLineWidth(2);
                h->SetLineColor(S.token == "os" ? kAzure + 2 : kOrange + 7);
                h->SetMarkerColor(h->GetLineColor());
                h->SetMarkerStyle(S.token == "os" ? 20 : 21);
                h->SetMarkerSize(0.6);
                h->GetXaxis()->SetTitleOffset(1.1);
                h->GetYaxis()->SetTitleOffset(1.3);
                h->SetMinimum(ylo);
                h->SetMaximum(yhi);
                h->Draw(first ? opt : opt_same);
                first = false;
            }
            // The proposed filter and the signal window, as graphical markers on the axis they
            // apply to -- so the reader can see the back-to-back peak's position relative to them.
            auto vline = [&](double x, Color_t col, Style_t sty) {
                auto l = std::make_unique<TLine>(x, ylo, x, yhi);
                l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(2);
                l->Draw();
                lines.push_back(std::move(l));
                return lines.back().get();
            };
            TLine* l_loose = vline(kLooseMassMax, kRed + 1, 2);
            TLine* l_sig_lo = vline(PairTrigEff::Window("sig").lo, kGray + 2, 3);
            vline(PairTrigEff::Window("sig").hi, kGray + 2, 3);

            // UPPER RIGHT, not upper left. The spectrum falls steeply above ~100 GeV, so the
            // top-right corner is the only region of the frame that is empty in every panel of
            // both figures; on the left the opposite-sign continuum and the J/psi spike run
            // through any box tall enough to hold four entries.
            auto leg = std::make_unique<TLegend>(0.55, 0.555, 0.965, 0.855);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
            for (const auto& S : SignsSsFirst())
                leg->AddEntry(hist_of(sl, S), S.text.c_str(), top_cells_only ? "lep" : "l");
            leg->AddEntry(l_loose, Fmt("m_{#mu#mu} = %g GeV", kLooseMassMax).c_str(), "l");
            leg->AddEntry(l_sig_lo, ranges[1].tex_text.c_str(), "l");
            leg->Draw();
            legs.push_back(std::move(leg));

            auto t = std::make_unique<TLatex>();
            t->SetNDC(); t->SetTextSize(0.040); t->SetTextFont(42);
            t->DrawLatex(0.13, 0.940, (sl.text + ",  " + scope).c_str());
            texs.push_back(std::move(t));
            // "no resonance veto" states a SELECTION, not a rationale, and it is the one fact a
            // reader needs to not mistake the omega and J/psi spikes for a defect: the analysis's
            // own opposite-sign trees ARE resonance-vetoed, this MC pair file is not.
            auto t2 = std::make_unique<TLatex>();
            t2->SetNDC(); t2->SetTextSize(0.031); t2->SetTextFont(42);
            t2->DrawLatex(0.13, 0.895,
                          Form("%s, %s muons, no resonance veto",
                               cfg.sample_text.c_str(), wp_text.c_str()));
            texs.push_back(std::move(t2));
        }
        // PNG only (memory: feedback_plot_format).
        c.SaveAs((out_dir + png).c_str());
    };
    draw_set("pair_mass_by_pthat_slice.png", false, "all selected pairs");
    draw_set("pair_mass_by_pthat_slice_top3_pairpt.png", true,
             Fmt("%.2f < p_{T}^{pair} < %.0f GeV", top_pt_lo, pt_edges.back()));

    // ---------------------------------------------------------------- the histograms, for reuse
    const std::string root_out = out_dir + "mc_pthat_slice_mass_stats.root";
    TFile fout(root_out.c_str(), "RECREATE");
    if (fout.IsZombie())
        throw std::runtime_error("mc_pthat_slice_mass_statistics: cannot open " + root_out);
    for (auto& kv : B) {
        // Drop the DISPLAY range before writing. draw_set set a min/max on these objects for the
        // figure; persisting it would hand every later reader a hard-coded axis that is not even
        // its own histogram's (the cap is the max over BOTH panels). -1111 is ROOT's "unset".
        for (TH1D* h : {kv.second.mass_all.GetPtr(), kv.second.mass_top.GetPtr()}) {
            h->SetMinimum(-1111);
            h->SetMaximum(-1111);
            h->Write();
        }
        for (auto& r : kv.second.n_cell) r.second->Write();
    }
    fout.Close();

    std::cout << "\n[mc_pthat_slice_mass_statistics] " << wp_text << " WP: "
              << 1 + 2 * (2 + static_cast<int>(ranges.size())) << " CSVs + 2 PNGs + "
              << root_out << " written to " << out_dir << std::endl;
}
